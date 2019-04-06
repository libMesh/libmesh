// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ includes
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unordered_map>

// C includes
#include <sys/types.h> // for pid_t
#include <unistd.h>    // for getpid(), unlink()

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/ghosting_functor.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h" // For n_levels
#include "libmesh/parallel.h"
#include "libmesh/remote_elem.h"
#include "libmesh/namebased_io.h"
#include "libmesh/partitioner.h"
#include "libmesh/enum_order.h"
#include "libmesh/mesh_communication.h"



namespace libMesh
{


// ------------------------------------------------------------
// UnstructuredMesh class member functions
UnstructuredMesh::UnstructuredMesh (const Parallel::Communicator & comm_in,
                                    unsigned char d) :
  MeshBase (comm_in,d)
{
  libmesh_assert (libMesh::initialized());
}



void UnstructuredMesh::copy_nodes_and_elements(const UnstructuredMesh & other_mesh,
                                               const bool skip_find_neighbors,
                                               dof_id_type element_id_offset,
                                               dof_id_type node_id_offset,
                                               unique_id_type
#ifdef LIBMESH_ENABLE_UNIQUE_ID
                                                 unique_id_offset
#endif
                                               )
{
  LOG_SCOPE("copy_nodes_and_elements()", "UnstructuredMesh");

  // If we are partitioned into fewer parts than the incoming mesh,
  // then we need to "wrap" the other Mesh's processor ids to fit
  // within our range. This can happen, for example, while stitching
  // ReplicatedMeshes with small numbers of elements in parallel...
  bool wrap_proc_ids = (_n_parts < other_mesh._n_parts);

  // We're assuming our subclass data needs no copy
  libmesh_assert_equal_to (_is_prepared, other_mesh._is_prepared);

  // We're assuming the other mesh has proper element number ordering,
  // so that we add parents before their children, and that the other
  // mesh is consistently partitioned.
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_amr_elem_ids(other_mesh);
  MeshTools::libmesh_assert_valid_procids<Node>(other_mesh);
#endif

  //Copy in Nodes
  {
    //Preallocate Memory if necessary
    this->reserve_nodes(other_mesh.n_nodes());

    for (const auto & oldn : other_mesh.node_ptr_range())
      {
        processor_id_type added_pid = cast_int<processor_id_type>
          (wrap_proc_ids ? oldn->processor_id() % _n_parts : oldn->processor_id());

        // Add new nodes in old node Point locations
#ifdef LIBMESH_ENABLE_UNIQUE_ID
        Node * newn =
#endif
          this->add_point(*oldn,
                          oldn->id() + node_id_offset,
                          added_pid);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
        newn->set_unique_id() =
          oldn->unique_id() + unique_id_offset;
#endif
      }
  }

  //Copy in Elements
  {
    //Preallocate Memory if necessary
    this->reserve_elem(other_mesh.n_elem());

    // Declare a map linking old and new elements, needed to copy the neighbor lists
    typedef std::unordered_map<const Elem *, Elem *> map_type;
    map_type old_elems_to_new_elems;

    // Loop over the elements
    for (const auto & old : other_mesh.element_ptr_range())
      {
        // Build a new element
        Elem * newparent = old->parent() ?
          this->elem_ptr(old->parent()->id() + element_id_offset) :
          nullptr;
        std::unique_ptr<Elem> ap = Elem::build(old->type(), newparent);
        Elem * el = ap.release();

        el->subdomain_id() = old->subdomain_id();

#ifdef LIBMESH_ENABLE_AMR
        if (old->has_children())
          for (unsigned int c = 0, nc = old->n_children(); c != nc; ++c)
            if (old->child_ptr(c) == remote_elem)
              el->add_child(const_cast<RemoteElem *>(remote_elem), c);

        //Create the parent's child pointers if necessary
        if (newparent)
          {
            unsigned int oldc = old->parent()->which_child_am_i(old);
            newparent->add_child(el, oldc);
          }

        // Copy the refinement flags
        el->set_refinement_flag(old->refinement_flag());

        // Use hack_p_level since we may not have sibling elements
        // added yet
        el->hack_p_level(old->p_level());

        el->set_p_refinement_flag(old->p_refinement_flag());
#endif // #ifdef LIBMESH_ENABLE_AMR

        //Assign all the nodes
        for (auto i : el->node_index_range())
          el->set_node(i) =
            this->node_ptr(old->node_id(i) + node_id_offset);

        // And start it off with the same processor id (mod _n_parts).
        el->processor_id() = cast_int<processor_id_type>
          (wrap_proc_ids ? old->processor_id() % _n_parts : old->processor_id());

        // Give it the same element and unique ids
        el->set_id(old->id() + element_id_offset);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
        el->set_unique_id() =
          old->unique_id() + unique_id_offset;
#endif

        //Hold onto it
        if (!skip_find_neighbors)
          {
            for (auto s : old->side_index_range())
              if (old->neighbor_ptr(s) == remote_elem)
                el->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));
            this->add_elem(el);
          }
        else
          {
            Elem * new_el = this->add_elem(el);
            old_elems_to_new_elems[old] = new_el;
          }
      }

    // Loop (again) over the elements to fill in the neighbors
    if (skip_find_neighbors)
      {
        old_elems_to_new_elems[remote_elem] = const_cast<RemoteElem*>(remote_elem);

        for (const auto & old_elem : other_mesh.element_ptr_range())
          {
            Elem * new_elem = old_elems_to_new_elems[old_elem];
            for (auto s : old_elem->side_index_range())
              {
                const Elem * old_neighbor = old_elem->neighbor_ptr(s);
                Elem * new_neighbor = old_elems_to_new_elems[old_neighbor];
                new_elem->set_neighbor(s, new_neighbor);
              }
          }
      }
  }

  //Finally prepare the new Mesh for use.  Keep the same numbering and
  //partitioning for now.
  this->allow_renumbering(false);
  this->allow_remote_element_removal(false);

  // We should generally be able to skip *all* partitioning here
  // because we're only adding one already-consistent mesh to another.
  this->skip_partitioning(true);

  this->prepare_for_use(false, skip_find_neighbors);

  //But in the long term, use the same renumbering and partitioning
  //policies as our source mesh.
  this->allow_renumbering(other_mesh.allow_renumbering());
  this->allow_remote_element_removal(other_mesh.allow_remote_element_removal());
  this->skip_partitioning(other_mesh._skip_all_partitioning);
  this->skip_noncritical_partitioning(other_mesh._skip_noncritical_partitioning);
}



UnstructuredMesh::~UnstructuredMesh ()
{
  //  this->clear ();  // Nothing to clear at this level

  libmesh_exceptionless_assert (!libMesh::closed());
}





void UnstructuredMesh::find_neighbors (const bool reset_remote_elements,
                                       const bool reset_current_list)
{
  // We might actually want to run this on an empty mesh
  // (e.g. the boundary mesh for a nonexistent bcid!)
  // libmesh_assert_not_equal_to (this->n_nodes(), 0);
  // libmesh_assert_not_equal_to (this->n_elem(), 0);

  // This function must be run on all processors at once
  parallel_object_only();

  LOG_SCOPE("find_neighbors()", "Mesh");

  //TODO:[BSK] This should be removed later?!
  if (reset_current_list)
    for (const auto & e : this->element_ptr_range())
      for (auto s : e->side_index_range())
        if (e->neighbor_ptr(s) != remote_elem || reset_remote_elements)
          e->set_neighbor(s, nullptr);

  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  {
    // data structures -- Use the hash_multimap if available
    typedef unsigned int                    key_type;
    typedef std::pair<Elem *, unsigned char> val_type;
    typedef std::pair<key_type, val_type>   key_val_pair;

    typedef std::unordered_multimap<key_type, val_type> map_type;

    // A map from side keys to corresponding elements & side numbers
    map_type side_to_elem_map;

    // Pull objects out of the loop to reduce heap operations
    std::unique_ptr<Elem> my_side, their_side;

    for (const auto & element : this->element_ptr_range())
      {
        for (auto ms : element->side_index_range())
          {
          next_side:
            // If we haven't yet found a neighbor on this side, try.
            // Even if we think our neighbor is remote, that
            // information may be out of date.
            if (element->neighbor_ptr(ms) == nullptr ||
                element->neighbor_ptr(ms) == remote_elem)
              {
                // Get the key for the side of this element
                const unsigned int key = element->key(ms);

                // Look for elements that have an identical side key
                auto bounds = side_to_elem_map.equal_range(key);

                // May be multiple keys, check all the possible
                // elements which _might_ be neighbors.
                if (bounds.first != bounds.second)
                  {
                    // Get the side for this element
                    element->side_ptr(my_side, ms);

                    // Look at all the entries with an equivalent key
                    while (bounds.first != bounds.second)
                      {
                        // Get the potential element
                        Elem * neighbor = bounds.first->second.first;

                        // Get the side for the neighboring element
                        const unsigned int ns = bounds.first->second.second;
                        neighbor->side_ptr(their_side, ns);
                        //libmesh_assert(my_side.get());
                        //libmesh_assert(their_side.get());

                        // If found a match with my side
                        //
                        // We need special tests here for 1D:
                        // since parents and children have an equal
                        // side (i.e. a node), we need to check
                        // ns != ms, and we also check level() to
                        // avoid setting our neighbor pointer to
                        // any of our neighbor's descendants
                        if ((*my_side == *their_side) &&
                            (element->level() == neighbor->level()) &&
                            ((element->dim() != 1) || (ns != ms)))
                          {
                            // So share a side.  Is this a mixed pair
                            // of subactive and active/ancestor
                            // elements?
                            // If not, then we're neighbors.
                            // If so, then the subactive's neighbor is

                            if (element->subactive() ==
                                neighbor->subactive())
                              {
                                // an element is only subactive if it has
                                // been coarsened but not deleted
                                element->set_neighbor (ms,neighbor);
                                neighbor->set_neighbor(ns,element);
                              }
                            else if (element->subactive())
                              {
                                element->set_neighbor(ms,neighbor);
                              }
                            else if (neighbor->subactive())
                              {
                                neighbor->set_neighbor(ns,element);
                              }
                            side_to_elem_map.erase (bounds.first);

                            // get out of this nested crap
                            goto next_side;
                          }

                        ++bounds.first;
                      }
                  }

                // didn't find a match...
                // Build the map entry for this element
                key_val_pair kvp;

                kvp.first         = key;
                kvp.second.first  = element;
                kvp.second.second = cast_int<unsigned char>(ms);
                side_to_elem_map.insert (kvp);
              }
          }
      }
  }

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Here we look at all of the child elements which
   * don't already have valid neighbors.
   *
   * If a child element has a nullptr neighbor it is
   * either because it is on the boundary or because
   * its neighbor is at a different level.  In the
   * latter case we must get the neighbor from the
   * parent.
   *
   * If a child element has a remote_elem neighbor
   * on a boundary it shares with its parent, that
   * info may have become out-dated through coarsening
   * of the neighbor's parent.  In this case, if the
   * parent's neighbor is active then the child should
   * share it.
   *
   * Furthermore, that neighbor better be active,
   * otherwise we missed a child somewhere.
   *
   *
   * We also need to look through children ordered by increasing
   * refinement level in order to add new interior_parent() links in
   * boundary elements which have just been generated by refinement,
   * and fix links in boundary elements whose previous
   * interior_parent() has just been coarsened away.
   */
  const unsigned int n_levels = MeshTools::n_levels(*this);
  for (unsigned int level = 1; level < n_levels; ++level)
    {
      for (auto & current_elem : as_range(level_elements_begin(level),
                                          level_elements_end(level)))
        {
          libmesh_assert(current_elem);
          Elem * parent = current_elem->parent();
          libmesh_assert(parent);
          const unsigned int my_child_num = parent->which_child_am_i(current_elem);

          for (auto s : current_elem->side_index_range())
            {
              if (current_elem->neighbor_ptr(s) == nullptr ||
                  (current_elem->neighbor_ptr(s) == remote_elem &&
                   parent->is_child_on_side(my_child_num, s)))
                {
                  Elem * neigh = parent->neighbor_ptr(s);

                  // If neigh was refined and had non-subactive children
                  // made remote earlier, then our current elem should
                  // actually have one of those remote children as a
                  // neighbor
                  if (neigh &&
                      (neigh->ancestor() ||
                       // If neigh has subactive children which should have
                       // matched as neighbors of the current element but
                       // did not, then those likewise must be remote
                       // children.
                       (current_elem->subactive() && neigh->has_children() &&
                        (neigh->level()+1) == current_elem->level())))
                    {
#ifdef DEBUG
                      // Let's make sure that "had children made remote"
                      // situation is actually the case
                      libmesh_assert(neigh->has_children());
                      bool neigh_has_remote_children = false;
                      for (auto & child : neigh->child_ref_range())
                        if (&child == remote_elem)
                          neigh_has_remote_children = true;
                      libmesh_assert(neigh_has_remote_children);

                      // And let's double-check that we don't have
                      // a remote_elem neighboring an active local element
                      if (current_elem->active())
                        libmesh_assert_not_equal_to (current_elem->processor_id(),
                                                     this->processor_id());
#endif // DEBUG
                      neigh = const_cast<RemoteElem *>(remote_elem);
                    }
                  // If neigh and current_elem are more than one level
                  // apart, figuring out whether we have a remote
                  // neighbor here becomes much harder.
                  else if (neigh && (current_elem->subactive() &&
                                     neigh->has_children()))
                    {
                      // Find the deepest descendant of neigh which
                      // we could consider for a neighbor.  If we run
                      // out of neigh children, then that's our
                      // neighbor.  If we find a potential neighbor
                      // with remote_children and we don't find any
                      // potential neighbors among its non-remote
                      // children, then our neighbor must be remote.
                      while (neigh != remote_elem &&
                             neigh->has_children())
                        {
                          bool found_neigh = false;
                          for (unsigned int c = 0, nc = neigh->n_children();
                               !found_neigh && c != nc; ++c)
                            {
                              Elem * child = neigh->child_ptr(c);
                              if (child == remote_elem)
                                continue;
                              for (auto ncn : child->neighbor_ptr_range())
                                {
                                  if (ncn != remote_elem &&
                                      ncn->is_ancestor_of(current_elem))
                                    {
                                      neigh = ncn;
                                      found_neigh = true;
                                      break;
                                    }
                                }
                            }
                          if (!found_neigh)
                            neigh = const_cast<RemoteElem *>(remote_elem);
                        }
                    }
                  current_elem->set_neighbor(s, neigh);
#ifdef DEBUG
                  if (neigh != nullptr && neigh != remote_elem)
                    // We ignore subactive elements here because
                    // we don't care about neighbors of subactive element.
                    if ((!neigh->active()) && (!current_elem->subactive()))
                      {
                        libMesh::err << "On processor " << this->processor_id()
                                     << std::endl;
                        libMesh::err << "Bad element ID = " << current_elem->id()
                                     << ", Side " << s << ", Bad neighbor ID = " << neigh->id() << std::endl;
                        libMesh::err << "Bad element proc_ID = " << current_elem->processor_id()
                                     << ", Bad neighbor proc_ID = " << neigh->processor_id() << std::endl;
                        libMesh::err << "Bad element size = " << current_elem->hmin()
                                     << ", Bad neighbor size = " << neigh->hmin() << std::endl;
                        libMesh::err << "Bad element center = " << current_elem->centroid()
                                     << ", Bad neighbor center = " << neigh->centroid() << std::endl;
                        libMesh::err << "ERROR: "
                                     << (current_elem->active()?"Active":"Ancestor")
                                     << " Element at level "
                                     << current_elem->level() << std::endl;
                        libMesh::err << "with "
                                     << (parent->active()?"active":
                                         (parent->subactive()?"subactive":"ancestor"))
                                     << " parent share "
                                     << (neigh->subactive()?"subactive":"ancestor")
                                     << " neighbor at level " << neigh->level()
                                     << std::endl;
                        NameBasedIO(*this).write ("bad_mesh.gmv");
                        libmesh_error_msg("Problematic mesh written to bad_mesh.gmv.");
                      }
#endif // DEBUG
                }
            }

          // We can skip to the next element if we're full-dimension
          // and therefore don't have any interior parents
          if (current_elem->dim() >= LIBMESH_DIM)
            continue;

          // We have no interior parents unless we can find one later
          current_elem->set_interior_parent(nullptr);

          Elem * pip = parent->interior_parent();

          if (!pip)
            continue;

          // If there's no interior_parent children, whether due to a
          // remote element or a non-conformity, then there's no
          // children to search.
          if (pip == remote_elem || pip->active())
            {
              current_elem->set_interior_parent(pip);
              continue;
            }

          // For node comparisons we'll need a sensible tolerance
          Real node_tolerance = current_elem->hmin() * TOLERANCE;

          // Otherwise our interior_parent should be a child of our
          // parent's interior_parent.
          for (auto & child : pip->child_ref_range())
            {
              // If we have a remote_elem, that might be our
              // interior_parent.  We'll set it provisionally now and
              // keep trying to find something better.
              if (&child == remote_elem)
                {
                  current_elem->set_interior_parent
                    (const_cast<RemoteElem *>(remote_elem));
                  continue;
                }

              bool child_contains_our_nodes = true;
              for (auto & n : current_elem->node_ref_range())
                {
                  bool child_contains_this_node = false;
                  for (auto & cn : child.node_ref_range())
                    if (cn.absolute_fuzzy_equals
                        (n, node_tolerance))
                      {
                        child_contains_this_node = true;
                        break;
                      }
                  if (!child_contains_this_node)
                    {
                      child_contains_our_nodes = false;
                      break;
                    }
                }
              if (child_contains_our_nodes)
                {
                  current_elem->set_interior_parent(&child);
                  break;
                }
            }

          // We should have found *some* interior_parent at this
          // point, whether semilocal or remote.
          libmesh_assert(current_elem->interior_parent());
        }
    }

#endif // AMR


#ifdef DEBUG
  MeshTools::libmesh_assert_valid_neighbors(*this,
                                            !reset_remote_elements);
  MeshTools::libmesh_assert_valid_amr_interior_parents(*this);
#endif
}



void UnstructuredMesh::read (const std::string & name,
                             void *,
                             bool skip_renumber_nodes_and_elements,
                             bool skip_find_neighbors)
{
  // Set the skip_renumber_nodes_and_elements flag on all processors
  // if necessary.
  // This ensures that renumber_nodes_and_elements is *not* called
  // during prepare_for_use() for certain types of mesh files.
  // This is required in cases where there is an associated solution
  // file which expects a certain ordering of the nodes.
  if (name.rfind(".gmv") + 4 == name.size())
    this->allow_renumbering(false);

  NameBasedIO(*this).read(name);

  if (skip_renumber_nodes_and_elements)
    {
      // Use MeshBase::allow_renumbering() yourself instead.
      libmesh_deprecated();
      this->allow_renumbering(false);
    }

  // Done reading the mesh.  Now prepare it for use.
  this->prepare_for_use(/*skip_renumber (deprecated)*/ false,
                        skip_find_neighbors);
}



void UnstructuredMesh::write (const std::string & name)
{
  LOG_SCOPE("write()", "Mesh");

  NameBasedIO(*this).write(name);
}



void UnstructuredMesh::write (const std::string & name,
                              const std::vector<Number> & v,
                              const std::vector<std::string> & vn)
{
  LOG_SCOPE("write()", "Mesh");

  NameBasedIO(*this).write_nodal_data(name, v, vn);
}





void UnstructuredMesh::create_pid_mesh(UnstructuredMesh & pid_mesh,
                                       const processor_id_type pid) const
{

  // Issue a warning if the number the number of processors
  // currently available is less that that requested for
  // partitioning.  This is not necessarily an error since
  // you may run on one processor and still partition the
  // mesh into several partitions.
#ifdef DEBUG
  if (this->n_processors() < pid)
    {
      libMesh::out << "WARNING:  You are creating a "
                   << "mesh for a processor id (="
                   << pid
                   << ") greater than "
                   << "the number of processors available for "
                   << "the calculation. (="
                   << this->n_processors()
                   << ")."
                   << std::endl;
    }
#endif

  this->create_submesh (pid_mesh,
                        this->active_pid_elements_begin(pid),
                        this->active_pid_elements_end(pid));
}







void UnstructuredMesh::create_submesh (UnstructuredMesh & new_mesh,
                                       const const_element_iterator & it,
                                       const const_element_iterator & it_end) const
{
  // Just in case the subdomain_mesh already has some information
  // in it, get rid of it.
  new_mesh.clear();

  // If we're not serial, our submesh isn't either.
  // There are no remote elements to delete on an empty mesh, but
  // calling the method to do so marks the mesh as parallel.
  if (!this->is_serial())
    new_mesh.delete_remote_elements();

  // Fail if (*this == new_mesh), we cannot create a submesh inside ourself!
  // This may happen if the user accidentally passes the original mesh into
  // this function!  We will check this by making sure we did not just
  // clear ourself.
  libmesh_assert_not_equal_to (this->n_nodes(), 0);
  libmesh_assert_not_equal_to (this->n_elem(), 0);

  // Container to catch boundary IDs handed back by BoundaryInfo
  std::vector<boundary_id_type> bc_ids;

  for (const auto & old_elem : as_range(it, it_end))
    {
      // Add an equivalent element type to the new_mesh.
      // Copy ids for this element.
      Elem * new_elem = Elem::build(old_elem->type()).release();
      new_elem->set_id() = old_elem->id();
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      new_elem->set_unique_id() = old_elem->unique_id();
#endif
      new_elem->subdomain_id() = old_elem->subdomain_id();
      new_elem->processor_id() = old_elem->processor_id();

      new_mesh.add_elem (new_elem);

      libmesh_assert(new_elem);

      // Loop over the nodes on this element.
      for (auto n : old_elem->node_index_range())
        {
          const dof_id_type this_node_id = old_elem->node_id(n);

          // Add this node to the new mesh if it's not there already
          if (!new_mesh.query_node_ptr(this_node_id))
            {
#ifdef LIBMESH_ENABLE_UNIQUE_ID
              Node * newn =
#endif
                new_mesh.add_point (old_elem->point(n),
                                    this_node_id,
                                    old_elem->node_ptr(n)->processor_id());

#ifdef LIBMESH_ENABLE_UNIQUE_ID
              newn->set_unique_id() = old_elem->node_ptr(n)->unique_id();
#endif
            }

          // Define this element's connectivity on the new mesh
          new_elem->set_node(n) = new_mesh.node_ptr(this_node_id);
        }

      // Maybe add boundary conditions for this element
      for (auto s : old_elem->side_index_range())
        {
          this->get_boundary_info().boundary_ids(old_elem, s, bc_ids);
          new_mesh.get_boundary_info().add_side (new_elem, s, bc_ids);
        }
    } // end loop over elements

  // Prepare the new_mesh for use
  new_mesh.prepare_for_use(/*skip_renumber =*/false);
}



#ifdef LIBMESH_ENABLE_AMR
bool UnstructuredMesh::contract ()
{
  LOG_SCOPE ("contract()", "Mesh");

  // Flag indicating if this call actually changes the mesh
  bool mesh_changed = false;

#ifdef DEBUG
  for (const auto & elem : this->element_ptr_range())
    libmesh_assert(elem->active() || elem->subactive() || elem->ancestor());
#endif

  // Loop over the elements.
  for (auto & elem : this->element_ptr_range())
    {
      // Delete all the subactive ones
      if (elem->subactive())
        {
          // No level-0 element should be subactive.
          // Note that we CAN'T test elem->level(), as that
          // touches elem->parent()->dim(), and elem->parent()
          // might have already been deleted!
          libmesh_assert(elem->parent());

          // Delete the element
          // This just sets a pointer to nullptr, and doesn't
          // invalidate any iterators
          this->delete_elem(elem);

          // the mesh has certainly changed
          mesh_changed = true;
        }
      else
        {
          // Compress all the active ones
          if (elem->active())
            elem->contract();
          else
            libmesh_assert (elem->ancestor());
        }
    }

  // Strip any newly-created nullptr voids out of the element array
  this->renumber_nodes_and_elements();

  // FIXME: Need to understand why deleting subactive children
  // invalidates the point locator.  For now we will clear it explicitly
  this->clear_point_locator();

  // Allow our GhostingFunctor objects to reinit if necessary.
  for (auto & gf : as_range(this->ghosting_functors_begin(),
                            this->ghosting_functors_end()))
    {
      libmesh_assert(gf);
      gf->mesh_reinit();
    }

  return mesh_changed;
}
#endif // #ifdef LIBMESH_ENABLE_AMR



void UnstructuredMesh::all_first_order ()
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  START_LOG("all_first_order()", "Mesh");

  /**
   * Prepare to identify (and then delete) a bunch of no-longer-used nodes.
   */
  std::vector<bool> node_touched_by_me(this->max_node_id(), false);

  // Loop over the high-ordered elements.
  // First make sure they _are_ indeed high-order, and then replace
  // them with an equivalent first-order element.
  for (auto & so_elem : element_ptr_range())
    {
      libmesh_assert(so_elem);

      /*
       * build the first-order equivalent, add to
       * the new_elements list.
       */
      Elem * lo_elem = Elem::build
        (Elem::first_order_equivalent_type
         (so_elem->type()), so_elem->parent()).release();

      const unsigned short n_sides = so_elem->n_sides();

      for (unsigned short s=0; s != n_sides; ++s)
        if (so_elem->neighbor_ptr(s) == remote_elem)
          lo_elem->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));

#ifdef LIBMESH_ENABLE_AMR
      /*
       * Reset the parent links of any child elements
       */
      if (so_elem->has_children())
        for (unsigned int c = 0, nc = so_elem->n_children(); c != nc; ++c)
          {
            Elem * child = so_elem->child_ptr(c);
            child->set_parent(lo_elem);
            lo_elem->add_child(child, c);
          }

      /*
       * Reset the child link of any parent element
       */
      if (so_elem->parent())
        {
          unsigned int c =
            so_elem->parent()->which_child_am_i(so_elem);
          lo_elem->parent()->replace_child(lo_elem, c);
        }

      /*
       * Copy as much data to the new element as makes sense
       */
      lo_elem->set_p_level(so_elem->p_level());
      lo_elem->set_refinement_flag(so_elem->refinement_flag());
      lo_elem->set_p_refinement_flag(so_elem->p_refinement_flag());
#endif

      libmesh_assert_equal_to (lo_elem->n_vertices(), so_elem->n_vertices());

      /*
       * By definition the vertices of the linear and
       * second order element are identically numbered.
       * transfer these.
       */
      for (unsigned int v=0; v < so_elem->n_vertices(); v++)
        {
          lo_elem->set_node(v) = so_elem->node_ptr(v);
          node_touched_by_me[lo_elem->node_id(v)] = true;
        }

      /*
       * find_neighbors relies on remote_elem neighbor links being
       * properly maintained.
       */
      for (unsigned short s=0; s != n_sides; s++)
        {
          if (so_elem->neighbor_ptr(s) == remote_elem)
            lo_elem->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));
        }

      /**
       * If the second order element had any boundary conditions they
       * should be transferred to the first-order element.  The old
       * boundary conditions will be removed from the BoundaryInfo
       * data structure by insert_elem.
       */
      this->get_boundary_info().copy_boundary_ids
        (this->get_boundary_info(), so_elem, lo_elem);

      /*
       * The new first-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the second-order element.
       */
      lo_elem->set_id(so_elem->id());
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      lo_elem->set_unique_id() = so_elem->unique_id();
#endif
      lo_elem->processor_id() = so_elem->processor_id();
      lo_elem->subdomain_id() = so_elem->subdomain_id();
      this->insert_elem(lo_elem);
    }

  // Deleting nodes does not invalidate iterators, so this is safe.
  for (const auto & node : this->node_ptr_range())
    if (!node_touched_by_me[node->id()])
      this->delete_node(node);

  // If crazy people applied boundary info to non-vertices and then
  // deleted those non-vertices, we should make sure their boundary id
  // caches are correct.
  this->get_boundary_info().regenerate_id_sets();

  STOP_LOG("all_first_order()", "Mesh");

  // On hanging nodes that used to also be second order nodes, we
  // might now have an invalid nodal processor_id()
  Partitioner::set_node_processor_ids(*this);

  // delete or renumber nodes if desired
  this->prepare_for_use();
}



void UnstructuredMesh::all_second_order (const bool full_ordered)
{
  // This function must be run on all processors at once
  parallel_object_only();

  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  /*
   * If the mesh is empty
   * then we have nothing to do
   */
  if (!this->n_elem())
    return;

  /*
   * If the mesh is already second order
   * then we have nothing to do.
   * We have to test for this in a round-about way to avoid
   * a bug on distributed parallel meshes with more processors
   * than elements.
   */
  bool already_second_order = false;
  if (this->elements_begin() != this->elements_end() &&
      (*(this->elements_begin()))->default_order() != FIRST)
    already_second_order = true;
  this->comm().max(already_second_order);
  if (already_second_order)
    return;

  START_LOG("all_second_order()", "Mesh");

  /*
   * this map helps in identifying second order
   * nodes.  Namely, a second-order node:
   * - edge node
   * - face node
   * - bubble node
   * is uniquely defined through a set of adjacent
   * vertices.  This set of adjacent vertices is
   * used to identify already added higher-order
   * nodes.  We are safe to use node id's since we
   * make sure that these are correctly numbered.
   */
  std::map<std::vector<dof_id_type>, Node *> adj_vertices_to_so_nodes;

  /*
   * for speed-up of the \p add_point() method, we
   * can reserve memory.  Guess the number of additional
   * nodes for different dimensions
   */
  switch (this->mesh_dimension())
    {
    case 1:
      /*
       * in 1D, there can only be order-increase from Edge2
       * to Edge3.  Something like 1/2 of n_nodes() have
       * to be added
       */
      this->reserve_nodes(static_cast<unsigned int>
                          (1.5*static_cast<double>(this->n_nodes())));
      break;

    case 2:
      /*
       * in 2D, either refine from Tri3 to Tri6 (double the nodes)
       * or from Quad4 to Quad8 (again, double) or Quad9 (2.25 that much)
       */
      this->reserve_nodes(static_cast<unsigned int>
                          (2*static_cast<double>(this->n_nodes())));
      break;


    case 3:
      /*
       * in 3D, either refine from Tet4 to Tet10 (factor = 2.5) up to
       * Hex8 to Hex27 (something  > 3).  Since in 3D there _are_ already
       * quite some nodes, and since we do not want to overburden the memory by
       * a too conservative guess, use the lower bound
       */
      this->reserve_nodes(static_cast<unsigned int>
                          (2.5*static_cast<double>(this->n_nodes())));
      break;

    default:
      // Hm?
      libmesh_error_msg("Unknown mesh dimension " << this->mesh_dimension());
    }



  /*
   * form a vector that will hold the node id's of
   * the vertices that are adjacent to the son-th
   * second-order node.  Pull this outside of the
   * loop so that silly compilers don't repeatedly
   * create and destroy the vector.
   */
  std::vector<dof_id_type> adjacent_vertices_ids;

  /**
   * Loop over the low-ordered elements in the _elements vector.
   * First make sure they _are_ indeed low-order, and then replace
   * them with an equivalent second-order element.  Don't
   * forget to delete the low-order element, or else it will leak!
   */
  element_iterator
    it = elements_begin(),
    endit = elements_end();

  for (; it != endit; ++it)
    {
      // the linear-order element
      Elem * lo_elem = *it;

      libmesh_assert(lo_elem);

      // make sure it is linear order
      if (lo_elem->default_order() != FIRST)
        libmesh_error_msg("ERROR: This is not a linear element: type=" << lo_elem->type());

      // this does _not_ work for refined elements
      libmesh_assert_equal_to (lo_elem->level (), 0);

      /*
       * build the second-order equivalent, add to
       * the new_elements list.  Note that this here
       * is the only point where \p full_ordered
       * is necessary.  The remaining code works well
       * for either type of second-order equivalent, e.g.
       * Hex20 or Hex27, as equivalents for Hex8
       */
      Elem * so_elem =
        Elem::build (Elem::second_order_equivalent_type(lo_elem->type(),
                                                        full_ordered) ).release();

      libmesh_assert_equal_to (lo_elem->n_vertices(), so_elem->n_vertices());


      /*
       * By definition the vertices of the linear and
       * second order element are identically numbered.
       * transfer these.
       */
      for (unsigned int v=0; v < lo_elem->n_vertices(); v++)
        so_elem->set_node(v) = lo_elem->node_ptr(v);

      /*
       * Now handle the additional mid-side nodes.  This
       * is simply handled through a map that remembers
       * the already-added nodes.  This map maps the global
       * ids of the vertices (that uniquely define this
       * higher-order node) to the new node.
       * Notation: son = second-order node
       */
      const unsigned int son_begin = so_elem->n_vertices();
      const unsigned int son_end   = so_elem->n_nodes();


      for (unsigned int son=son_begin; son<son_end; son++)
        {
          const unsigned int n_adjacent_vertices =
            so_elem->n_second_order_adjacent_vertices(son);

          adjacent_vertices_ids.resize(n_adjacent_vertices);

          for (unsigned int v=0; v<n_adjacent_vertices; v++)
            adjacent_vertices_ids[v] =
              so_elem->node_id( so_elem->second_order_adjacent_vertex(son,v) );

          /*
           * \p adjacent_vertices_ids is now in order of the current
           * side.  sort it, so that comparisons  with the
           * \p adjacent_vertices_ids created through other elements'
           * sides can match
           */
          std::sort(adjacent_vertices_ids.begin(),
                    adjacent_vertices_ids.end());


          // does this set of vertices already have a mid-node added?
          auto pos = adj_vertices_to_so_nodes.equal_range (adjacent_vertices_ids);

          // no, not added yet
          if (pos.first == pos.second)
            {
              /*
               * for this set of vertices, there is no
               * second_order node yet.  Add it.
               *
               * compute the location of the new node as
               * the average over the adjacent vertices.
               */
              Point new_location = this->point(adjacent_vertices_ids[0]);
              for (unsigned int v=1; v<n_adjacent_vertices; v++)
                new_location += this->point(adjacent_vertices_ids[v]);

              new_location /= static_cast<Real>(n_adjacent_vertices);

              /* Add the new point to the mesh.
               * If we are on a serialized mesh, then we're doing this
               * all in sync, and the node processor_id will be
               * consistent between processors.
               * If we are on a distributed mesh, we can fix
               * inconsistent processor ids later, but only if every
               * processor gives new nodes a *locally* consistent
               * processor id, so we'll give the new node the
               * processor id of an adjacent element for now and then
               * we'll update that later if appropriate.
               */
              Node * so_node = this->add_point
                (new_location, DofObject::invalid_id,
                 lo_elem->processor_id());

              /*
               * insert the new node with its defining vertex
               * set into the map, and relocate pos to this
               * new entry, so that the so_elem can use
               * \p pos for inserting the node
               */
              adj_vertices_to_so_nodes.insert(pos.first,
                                              std::make_pair(adjacent_vertices_ids,
                                                             so_node));

              so_elem->set_node(son) = so_node;
            }
          // yes, already added.
          else
            {
              Node * so_node = pos.first->second;
              libmesh_assert(so_node);

              so_elem->set_node(son) = so_node;

              // We need to ensure that the processor who should own a
              // node *knows* they own the node.  And because
              // Node::choose_processor_id() may depend on Node id,
              // which may not yet be authoritative, we still have to
              // use a dumb-but-id-independent partitioning heuristic.
              processor_id_type chosen_pid =
                std::min (so_node->processor_id(),
                          lo_elem->processor_id());

              // Plus, if we just discovered that we own this node,
              // then on a distributed mesh we need to make sure to
              // give it a valid id, not just a placeholder id!
              if (!this->is_replicated() &&
                  so_node->processor_id() != this->processor_id() &&
                  chosen_pid == this->processor_id())
                this->own_node(*so_node);

              so_node->processor_id() = chosen_pid;
            }
        }

      /*
       * find_neighbors relies on remote_elem neighbor links being
       * properly maintained.
       */
      for (auto s : lo_elem->side_index_range())
        {
          if (lo_elem->neighbor_ptr(s) == remote_elem)
            so_elem->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));
        }

      /**
       * If the linear element had any boundary conditions they
       * should be transferred to the second-order element.  The old
       * boundary conditions will be removed from the BoundaryInfo
       * data structure by insert_elem.
       *
       * Also, prepare_for_use() will reconstruct most of our neighbor
       * links, but if we have any remote_elem links in a distributed
       * mesh, they need to be preserved.  We do that in the same loop
       * here.
       */
      this->get_boundary_info().copy_boundary_ids
        (this->get_boundary_info(), lo_elem, so_elem);

      /*
       * The new second-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the first-order element.
       */
      so_elem->set_id(lo_elem->id());
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      so_elem->set_unique_id() = lo_elem->unique_id();
#endif
      so_elem->processor_id() = lo_elem->processor_id();
      so_elem->subdomain_id() = lo_elem->subdomain_id();
      this->insert_elem(so_elem);
    }

  // we can clear the map
  adj_vertices_to_so_nodes.clear();


  STOP_LOG("all_second_order()", "Mesh");

  // On a DistributedMesh our ghost node processor ids may be bad,
  // the ids of nodes touching remote elements may be inconsistent,
  // unique_ids of newly added non-local nodes remain unset, and our
  // partitioning of new nodes may not be well balanced.
  //
  // make_nodes_parallel_consistent() will fix all this.
  if (!this->is_serial())
    MeshCommunication().make_nodes_parallel_consistent (*this);

  // renumber nodes, elements etc
  this->prepare_for_use(/*skip_renumber =*/ false);
}

} // namespace libMesh
