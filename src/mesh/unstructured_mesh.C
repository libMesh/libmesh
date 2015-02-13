// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C includes
#include <sys/types.h> // for pid_t
#include <unistd.h>    // for getpid(), unlink()

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h" // For n_levels
#include "libmesh/parallel.h"
#include "libmesh/remote_elem.h"

// For most I/O
#include "libmesh/namebased_io.h"

// for MeshData backward compatibility
#include "libmesh/unv_io.h"
#include "libmesh/tetgen_io.h"

#include LIBMESH_INCLUDE_UNORDERED_MAP


namespace libMesh
{


// ------------------------------------------------------------
// UnstructuredMesh class member functions
UnstructuredMesh::UnstructuredMesh (const Parallel::Communicator &comm_in,
                                    unsigned char d) :
  MeshBase (comm_in,d)
{
  libmesh_assert (libMesh::initialized());
}



#ifndef LIBMESH_DISABLE_COMMWORLD
UnstructuredMesh::UnstructuredMesh (unsigned char d) :
  MeshBase (d)
{
  libmesh_assert (libMesh::initialized());
}
#endif



void UnstructuredMesh::copy_nodes_and_elements
(const UnstructuredMesh& other_mesh, const bool skip_find_neighbors)
{
  // We're assuming our subclass data needs no copy
  libmesh_assert_equal_to (_n_parts, other_mesh._n_parts);
  libmesh_assert (std::equal(_elem_dims.begin(), _elem_dims.end(), other_mesh._elem_dims.begin()));
  libmesh_assert_equal_to (_is_prepared, other_mesh._is_prepared);

  // We're assuming the other mesh has proper element number ordering,
  // so that we add parents before their children.
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_amr_elem_ids(other_mesh);
#endif

  //Copy in Nodes
  {
    //Preallocate Memory if necessary
    this->reserve_nodes(other_mesh.n_nodes());

    const_node_iterator it = other_mesh.nodes_begin();
    const_node_iterator end = other_mesh.nodes_end();

    for (; it != end; ++it)
      {
        const Node *oldn = *it;

        // Add new nodes in old node Point locations
        /*Node *newn =*/ this->add_point(*oldn, oldn->id(), oldn->processor_id());

        // And start them off in the same subdomain
        //        newn->processor_id() = oldn->processor_id();
      }
  }

  //Copy in Elements
  {
    //Preallocate Memory if necessary
    this->reserve_elem(other_mesh.n_elem());

    // Declare a map linking old and new elements, needed to copy the neighbor lists
    std::map<const Elem*, Elem*> old_elems_to_new_elems;

    // Loop over the elements
    MeshBase::const_element_iterator it = other_mesh.elements_begin();
    const MeshBase::const_element_iterator end = other_mesh.elements_end();

    // FIXME: Where do we set element IDs??
    for (; it != end; ++it)
      {
        //Look at the old element
        const Elem *old = *it;
        //Build a new element
        Elem *newparent = old->parent() ?
          this->elem(old->parent()->id()) : NULL;
        AutoPtr<Elem> ap = Elem::build(old->type(), newparent);
        Elem * el = ap.release();

        el->subdomain_id() = old->subdomain_id();

        for (unsigned int s=0; s != old->n_sides(); ++s)
          if (old->neighbor(s) == remote_elem)
            el->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));

#ifdef LIBMESH_ENABLE_AMR
        if (old->has_children())
          for (unsigned int c=0; c != old->n_children(); ++c)
            if (old->child(c) == remote_elem)
              el->add_child(const_cast<RemoteElem*>(remote_elem), c);

        //Create the parent's child pointers if necessary
        if (newparent)
          {
            unsigned int oldc = old->parent()->which_child_am_i(old);
            newparent->add_child(el, oldc);
          }

        // Copy the refinement flags
        el->set_refinement_flag(old->refinement_flag());
        el->set_p_refinement_flag(old->p_refinement_flag());
#endif // #ifdef LIBMESH_ENABLE_AMR

        //Assign all the nodes
        for(unsigned int i=0;i<el->n_nodes();i++)
          el->set_node(i) = &this->node(old->node(i));

        // And start it off in the same subdomain
        el->processor_id() = old->processor_id();

        // Give it the same id
        el->set_id(old->id());

        //Hold onto it
        if(!skip_find_neighbors)
          {
            this->add_elem(el);
          }
        else
          {
            Elem* new_el = this->add_elem(el);
            old_elems_to_new_elems[old] = new_el;
          }

        // Add the link between the original element and this copy to the map
        if(skip_find_neighbors)
          old_elems_to_new_elems[old] = el;
      }

    // Loop (again) over the elements to fill in the neighbors
    if(skip_find_neighbors)
      {
        it = other_mesh.elements_begin();
        for (; it != end; ++it)
          {
            Elem* old_elem = *it;
            Elem* new_elem = old_elems_to_new_elems[old_elem];
            for (unsigned int s=0; s != old_elem->n_neighbors(); ++s)
              {
                const Elem* old_neighbor = old_elem->neighbor(s);
                Elem* new_neighbor = old_elems_to_new_elems[old_neighbor];
                new_elem->set_neighbor(s, new_neighbor);
              }
          }
      }
  }

  //Finally prepare the new Mesh for use.  Keep the same numbering and
  //partitioning but also the same renumbering and partitioning
  //policies as our source mesh.
  this->allow_renumbering(false);
  this->skip_partitioning(true);
  this->prepare_for_use(false, skip_find_neighbors);
  this->allow_renumbering(other_mesh.allow_renumbering());
  this->skip_partitioning(other_mesh.skip_partitioning());
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
  // (e.g. the boundary mesh for a nonexistant bcid!)
  // libmesh_assert_not_equal_to (this->n_nodes(), 0);
  // libmesh_assert_not_equal_to (this->n_elem(), 0);

  // This function must be run on all processors at once
  parallel_object_only();

  START_LOG("find_neighbors()", "Mesh");

  const element_iterator el_end = this->elements_end();

  //TODO:[BSK] This should be removed later?!
  if (reset_current_list)
    for (element_iterator el = this->elements_begin(); el != el_end; ++el)
      {
        Elem* e = *el;
        for (unsigned int s=0; s<e->n_neighbors(); s++)
          if (e->neighbor(s) != remote_elem ||
              reset_remote_elements)
            e->set_neighbor(s,NULL);
      }

  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  {
    // data structures -- Use the hash_multimap if available
    typedef unsigned int                    key_type;
    typedef std::pair<Elem*, unsigned char> val_type;
    typedef std::pair<key_type, val_type>   key_val_pair;

    typedef LIBMESH_BEST_UNORDERED_MULTIMAP<key_type, val_type> map_type;

    // A map from side keys to corresponding elements & side numbers
    map_type side_to_elem_map;



    for (element_iterator el = this->elements_begin(); el != el_end; ++el)
      {
        Elem* element = *el;

        for (unsigned char ms=0; ms<element->n_neighbors(); ms++)
          {
          next_side:
            // If we haven't yet found a neighbor on this side, try.
            // Even if we think our neighbor is remote, that
            // information may be out of date.
            if (element->neighbor(ms) == NULL ||
                element->neighbor(ms) == remote_elem)
              {
                // Get the key for the side of this element
                const unsigned int key = element->key(ms);

                // Look for elements that have an identical side key
                std::pair <map_type::iterator, map_type::iterator>
                  bounds = side_to_elem_map.equal_range(key);

                // May be multiple keys, check all the possible
                // elements which _might_ be neighbors.
                if (bounds.first != bounds.second)
                  {
                    // Get the side for this element
                    const AutoPtr<Elem> my_side(element->side(ms));

                    // Look at all the entries with an equivalent key
                    while (bounds.first != bounds.second)
                      {
                        // Get the potential element
                        Elem* neighbor = bounds.first->second.first;

                        // Get the side for the neighboring element
                        const unsigned int ns = bounds.first->second.second;
                        const AutoPtr<Elem> their_side(neighbor->side(ns));
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
                        if( (*my_side == *their_side) &&
                            (element->level() == neighbor->level()) &&
                            ((element->dim() != 1) || (ns != ms)) )
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
                kvp.second.second = ms;

                // use the lower bound as a hint for
                // where to put it.
#if defined(LIBMESH_HAVE_UNORDERED_MAP) || defined(LIBMESH_HAVE_TR1_UNORDERED_MAP) || defined(LIBMESH_HAVE_HASH_MAP) || defined(LIBMESH_HAVE_EXT_HASH_MAP)
                side_to_elem_map.insert (kvp);
#else
                side_to_elem_map.insert (bounds.first,kvp);
#endif
              }
          }
      }
  }

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Here we look at all of the child elements which
   * don't already have valid neighbors.
   *
   * If a child element has a NULL neighbor it is
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
   */
  const unsigned int n_levels = MeshTools::n_levels(*this);
  for (unsigned int level = 1; level < n_levels; ++level)
    {
      element_iterator end = this->level_elements_end(level);
      for (element_iterator el = this->level_elements_begin(level);
           el != end; ++el)
        {
          Elem* current_elem = *el;
          libmesh_assert(current_elem);
          Elem* parent = current_elem->parent();
          libmesh_assert(parent);
          const unsigned int my_child_num = parent->which_child_am_i(current_elem);

          for (unsigned int s=0; s < current_elem->n_neighbors(); s++)
            {
              if (current_elem->neighbor(s) == NULL ||
                  (current_elem->neighbor(s) == remote_elem &&
                   parent->is_child_on_side(my_child_num, s)))
                {
                  Elem *neigh = parent->neighbor(s);

                  // If neigh was refined and had non-subactive children
                  // made remote earlier, then a non-subactive elem should
                  // actually have one of those remote children as a
                  // neighbor
                  if (neigh && (neigh->ancestor()) && (!current_elem->subactive()))
                    {
#ifdef DEBUG
                      // Let's make sure that "had children made remote"
                      // situation is actually the case
                      libmesh_assert(neigh->has_children());
                      bool neigh_has_remote_children = false;
                      for (unsigned int c = 0; c != neigh->n_children(); ++c)
                        {
                          if (neigh->child(c) == remote_elem)
                            neigh_has_remote_children = true;
                        }
                      libmesh_assert(neigh_has_remote_children);

                      // And let's double-check that we don't have
                      // a remote_elem neighboring a local element
                      libmesh_assert_not_equal_to (current_elem->processor_id(),
                                                   this->processor_id());
#endif // DEBUG
                      neigh = const_cast<RemoteElem*>(remote_elem);
                    }

                  current_elem->set_neighbor(s, neigh);
#ifdef DEBUG
                  if (neigh != NULL && neigh != remote_elem)
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
        }
    }

#endif // AMR


#ifdef DEBUG
  MeshTools::libmesh_assert_valid_neighbors(*this);
#endif

  STOP_LOG("find_neighbors()", "Mesh");
}



void UnstructuredMesh::read (const std::string& name,
                             MeshData* mesh_data,
                             bool skip_renumber_nodes_and_elements)
{
  // Set the skip_renumber_nodes_and_elements flag on all processors
  // if necessary.
  // This ensures that renumber_nodes_and_elements is *not* called
  // during prepare_for_use() for certain types of mesh files.
  // This is required in cases where there is an associated solution
  // file which expects a certain ordering of the nodes.
  if(name.rfind(".gmv")+4==name.size())
    {
      skip_renumber_nodes_and_elements =  true;
    }

  if (mesh_data)
    {
      libmesh_deprecated();
      if (name.rfind(".unv") < name.size())
        UNVIO(*this, mesh_data).read (name);
      else if ((name.rfind(".node")  < name.size()) ||
               (name.rfind(".ele")   < name.size()))
        TetGenIO(*this,mesh_data).read (name);
    }
  else
    NameBasedIO(*this).read(name);

  if (skip_renumber_nodes_and_elements)
    {
      // Use MeshBase::allow_renumbering() yourself instead.
      libmesh_deprecated();
      this->allow_renumbering(false);
    }

  // Done reading the mesh.  Now prepare it for use.
  this->prepare_for_use();
}



void UnstructuredMesh::write (const std::string& name,
                              MeshData* mesh_data)
{
  START_LOG("write()", "Mesh");

  if (mesh_data)
    {
      libmesh_deprecated();
      if (name.rfind(".unv") < name.size())
        UNVIO(*this, mesh_data).write (name);
      else
        libmesh_error_msg("Only UNV output supports MeshData");
    }

  NameBasedIO(*this).write(name);

  STOP_LOG("write()", "Mesh");
}



void UnstructuredMesh::write (const std::string& name,
                              const std::vector<Number>& v,
                              const std::vector<std::string>& vn)
{
  START_LOG("write()", "Mesh");

  NameBasedIO(*this).write_nodal_data(name, v, vn);

  STOP_LOG("write()", "Mesh");
}





void UnstructuredMesh::create_pid_mesh(UnstructuredMesh& pid_mesh,
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

  // Create iterators to loop over the list of elements
  //   const_active_pid_elem_iterator       it(this->elements_begin(),   pid);
  //   const const_active_pid_elem_iterator it_end(this->elements_end(), pid);

  const_element_iterator       it     = this->active_pid_elements_begin(pid);
  const const_element_iterator it_end = this->active_pid_elements_end(pid);

  this->create_submesh (pid_mesh, it, it_end);
}







void UnstructuredMesh::create_submesh (UnstructuredMesh& new_mesh,
                                       const_element_iterator& it,
                                       const const_element_iterator& it_end) const
{
  // Just in case the subdomain_mesh already has some information
  // in it, get rid of it.
  new_mesh.clear();

  // Fail if (*this == new_mesh), we cannot create a submesh inside ourself!
  // This may happen if the user accidently passes the original mesh into
  // this function!  We will check this by making sure we did not just
  // clear ourself.
  libmesh_assert_not_equal_to (this->n_nodes(), 0);
  libmesh_assert_not_equal_to (this->n_elem(), 0);

  // How the nodes on this mesh will be renumbered to nodes
  // on the new_mesh.
  std::vector<dof_id_type> new_node_numbers (this->n_nodes());

  std::fill (new_node_numbers.begin(),
             new_node_numbers.end(),
             DofObject::invalid_id);



  // the number of nodes on the new mesh, will be incremented
  dof_id_type n_new_nodes = 0;
  dof_id_type n_new_elem  = 0;

  for (; it != it_end; ++it)
    {
      // increment the new element counter
      n_new_elem++;

      const Elem* old_elem = *it;

      // Add an equivalent element type to the new_mesh
      Elem* new_elem =
        new_mesh.add_elem (Elem::build(old_elem->type()).release());

      libmesh_assert(new_elem);

      // Loop over the nodes on this element.
      for (unsigned int n=0; n<old_elem->n_nodes(); n++)
        {
          libmesh_assert_less (old_elem->node(n), new_node_numbers.size());

          if (new_node_numbers[old_elem->node(n)] == DofObject::invalid_id)
            {
              new_node_numbers[old_elem->node(n)] = n_new_nodes;

              // Add this node to the new mesh
              new_mesh.add_point (old_elem->point(n));

              // Increment the new node counter
              n_new_nodes++;
            }

          // Define this element's connectivity on the new mesh
          libmesh_assert_less (new_node_numbers[old_elem->node(n)], new_mesh.n_nodes());

          new_elem->set_node(n) = new_mesh.node_ptr (new_node_numbers[old_elem->node(n)]);
        }

      // Copy ids for this element
      new_elem->subdomain_id() = old_elem->subdomain_id();
      new_elem->processor_id() = old_elem->processor_id();

      // Maybe add boundary conditions for this element
      for (unsigned short s=0; s<old_elem->n_sides(); s++)
        // We're supporting boundary ids on internal sides now
        //if (old_elem->neighbor(s) == NULL)
        {
          const std::vector<boundary_id_type>& bc_ids =
            this->get_boundary_info().boundary_ids(old_elem, s);
          for (std::vector<boundary_id_type>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
            {
              const boundary_id_type bc_id = *id_it;
              if (bc_id != this->get_boundary_info().invalid_id)
                new_mesh.get_boundary_info().add_side (new_elem, s,
                                                       bc_id);
            }
        }
    } // end loop over elements


  // Prepare the new_mesh for use
  new_mesh.prepare_for_use(/*skip_renumber =*/false);

}



#ifdef LIBMESH_ENABLE_AMR
bool UnstructuredMesh::contract ()
{
  START_LOG ("contract()", "Mesh");

  // Flag indicating if this call actually changes the mesh
  bool mesh_changed = false;

  element_iterator in        = elements_begin();
  const element_iterator end = elements_end();

#ifdef DEBUG
  for ( ; in != end; ++in)
    if (*in != NULL)
      {
        Elem* el = *in;
        libmesh_assert(el->active() || el->subactive() || el->ancestor());
      }
  in = elements_begin();
#endif

  // Loop over the elements.
  for ( ; in != end; ++in)
    if (*in != NULL)
      {
        Elem* el = *in;

        // Delete all the subactive ones
        if (el->subactive())
          {
            // No level-0 element should be subactive.
            // Note that we CAN'T test elem->level(), as that
            // touches elem->parent()->dim(), and elem->parent()
            // might have already been deleted!
            libmesh_assert(el->parent());

            // Delete the element
            // This just sets a pointer to NULL, and doesn't
            // invalidate any iterators
            this->delete_elem(el);

            // the mesh has certainly changed
            mesh_changed = true;
          }
        else
          {
            // Compress all the active ones
            if (el->active())
              el->contract();
            else
              libmesh_assert (el->ancestor());
          }
      }

  // Strip any newly-created NULL voids out of the element array
  this->renumber_nodes_and_elements();

  // FIXME: Need to understand why deleting subactive children
  // invalidates the point locator.  For now we will clear it explicitly
  this->clear_point_locator();

  STOP_LOG ("contract()", "Mesh");

  return mesh_changed;
}
#endif // #ifdef LIBMESH_ENABLE_AMR

} // namespace libMesh
