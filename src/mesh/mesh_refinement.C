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
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for isnan(), when it's defined
#include <limits>

// Local includes
#include "libmesh/libmesh_config.h"

// only compile these functions if the user requests AMR support
#ifdef LIBMESH_ENABLE_AMR

#include "libmesh/boundary_info.h"
#include "libmesh/error_vector.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/partitioner.h"
#include "libmesh/remote_elem.h"
#include "libmesh/sync_refinement_flags.h"
#include "libmesh/int_range.h"

#ifdef DEBUG
// Some extra validation for DistributedMesh
#include "libmesh/mesh_tools.h"
#endif // DEBUG

#ifdef LIBMESH_ENABLE_PERIODIC
#include "libmesh/periodic_boundaries.h"
#endif



// Anonymous namespace for helper functions
// namespace {
//
// using namespace libMesh;
//
// struct SyncCoarsenInactive
// {
//   bool operator() (const Elem * elem) const {
//     // If we're not an ancestor, there's no chance our coarsening
//     // settings need to be changed.
//     if (!elem->ancestor())
//       return false;
//
//     // If we don't have any remote children, we already know enough to
//     // determine the correct refinement flag ourselves.
//     //
//     // If we know we have a child that isn't being coarsened, that
//     // also forces a specific flag.
//     //
//     // Either way there's nothing we need to communicate.
//     bool found_remote_child = false;
//     for (auto & child : elem->child_ref_range())
//       {
//         if (child.refinement_flag() != Elem::COARSEN)
//           return false;
//         if (&child == remote_elem)
//           found_remote_child = true;
//       }
//     return found_remote_child;
//   }
// };
// }



namespace libMesh
{

//-----------------------------------------------------------------
// Mesh refinement methods
MeshRefinement::MeshRefinement (MeshBase & m) :
  ParallelObject(m),
  _mesh(m),
  _use_member_parameters(false),
  _coarsen_by_parents(false),
  _refine_fraction(0.3),
  _coarsen_fraction(0.0),
  _max_h_level(libMesh::invalid_uint),
  _coarsen_threshold(10),
  _nelem_target(0),
  _absolute_global_tolerance(0.0),
  _face_level_mismatch_limit(1),
  _edge_level_mismatch_limit(0),
  _node_level_mismatch_limit(0),
  _overrefined_boundary_limit(0),
  _underrefined_boundary_limit(0),
  _enforce_mismatch_limit_prior_to_refinement(false)
#ifdef LIBMESH_ENABLE_PERIODIC
  , _periodic_boundaries(nullptr)
#endif
{
}



#ifdef LIBMESH_ENABLE_PERIODIC
void MeshRefinement::set_periodic_boundaries_ptr(PeriodicBoundaries * pb_ptr)
{
  _periodic_boundaries = pb_ptr;
}
#endif



MeshRefinement::~MeshRefinement ()
{
  this->clear();
}



void MeshRefinement::clear ()
{
  _new_nodes_map.clear();
}



Node * MeshRefinement::add_node(Elem & parent,
                                unsigned int child,
                                unsigned int node,
                                processor_id_type proc_id)
{
  LOG_SCOPE("add_node()", "MeshRefinement");

  unsigned int parent_n = parent.as_parent_node(child, node);

  if (parent_n != libMesh::invalid_uint)
    return parent.node_ptr(parent_n);

  const std::vector<std::pair<dof_id_type, dof_id_type>>
    bracketing_nodes = parent.bracketing_nodes(child, node);

  // If we're not a parent node, we *must* be bracketed by at least
  // one pair of parent nodes
  libmesh_assert(bracketing_nodes.size());

  const dof_id_type new_node_id =
    _new_nodes_map.find(bracketing_nodes);

  // Return the node if it already exists.
  //
  // We'll leave the processor_id untouched in this case - if we're
  // repartitioning later or if this is a new unpartitioned node,
  // we'll update it then, and if not then we don't want to update it.
  if (new_node_id != DofObject::invalid_id)
    return _mesh.node_ptr(new_node_id);

  // Otherwise we need to add a new node.
  //
  // Figure out where to add the point:

  Point p; // defaults to 0,0,0

  for (auto n : parent.node_index_range())
    {
      // The value from the embedding matrix
      const float em_val = parent.embedding_matrix(child,node,n);

      if (em_val != 0.)
        {
          p.add_scaled (parent.point(n), em_val);

          // If we'd already found the node we shouldn't be here
          libmesh_assert_not_equal_to (em_val, 1);
        }
    }

  // Although we're leaving new nodes unpartitioned at first, with a
  // DistributedMesh we would need a default id based on the numbering
  // scheme for the requested processor_id.
  Node * new_node = _mesh.add_point (p, DofObject::invalid_id, proc_id);

  libmesh_assert(new_node);

  // But then we'll make sure this node is marked as unpartitioned.
  new_node->processor_id() = DofObject::invalid_processor_id;

  // Add the node to the map.
  _new_nodes_map.add_node(*new_node, bracketing_nodes);

  // Return the address of the new node
  return new_node;
}



Elem * MeshRefinement::add_elem (Elem * elem)
{
  libmesh_assert(elem);
  _mesh.add_elem (elem);
  return elem;
}



void MeshRefinement::create_parent_error_vector(const ErrorVector & error_per_cell,
                                                ErrorVector & error_per_parent,
                                                Real & parent_error_min,
                                                Real & parent_error_max)
{
  // This function must be run on all processors at once
  parallel_object_only();

  // Make sure the input error vector is valid
#ifdef DEBUG
  for (const auto & val : error_per_cell)
    {
      libmesh_assert_greater_equal (val, 0);
      // isnan() isn't standard C++ yet
#ifdef isnan
      libmesh_assert(!isnan(val));
#endif
    }

  // Use a reference to std::vector to avoid confusing
  // this->comm().verify
  std::vector<ErrorVectorReal> & epc = error_per_parent;
  libmesh_assert(this->comm().verify(epc));
#endif // #ifdef DEBUG

  // error values on uncoarsenable elements will be left at -1
  error_per_parent.clear();
  error_per_parent.resize(error_per_cell.size(), 0.0);

  {
    // Find which elements are uncoarsenable
    for (auto & elem : _mesh.active_local_element_ptr_range())
      {
        Elem * parent = elem->parent();

        // Active elements are uncoarsenable
        error_per_parent[elem->id()] = -1.0;

        // Grandparents and up are uncoarsenable
        while (parent)
          {
            parent = parent->parent();
            if (parent)
              {
                const dof_id_type parentid  = parent->id();
                libmesh_assert_less (parentid, error_per_parent.size());
                error_per_parent[parentid] = -1.0;
              }
          }
      }

    // Sync between processors.
    // Use a reference to std::vector to avoid confusing
    // this->comm().min
    std::vector<ErrorVectorReal> & epp = error_per_parent;
    this->comm().min(epp);
  }

  // The parent's error is defined as the square root of the
  // sum of the children's errors squared, so errors that are
  // Hilbert norms remain Hilbert norms.
  //
  // Because the children may be on different processors, we
  // calculate local contributions to the parents' errors squared
  // first, then sum across processors and take the square roots
  // second.
  for (auto & elem : _mesh.active_local_element_ptr_range())
    {
      Elem * parent = elem->parent();

      // Calculate each contribution to parent cells
      if (parent)
        {
          const dof_id_type parentid  = parent->id();
          libmesh_assert_less (parentid, error_per_parent.size());

          // If the parent has grandchildren we won't be able to
          // coarsen it, so forget it.  Otherwise, add this child's
          // contribution to the sum of the squared child errors
          if (error_per_parent[parentid] != -1.0)
            error_per_parent[parentid] += (error_per_cell[elem->id()] *
                                           error_per_cell[elem->id()]);
        }
    }

  // Sum the vector across all processors
  this->comm().sum(static_cast<std::vector<ErrorVectorReal> &>(error_per_parent));

  // Calculate the min and max as we loop
  parent_error_min = std::numeric_limits<double>::max();
  parent_error_max = 0.;

  for (auto i : index_range(error_per_parent))
    {
      // If this element isn't a coarsenable parent with error, we
      // have nothing to do.  Just flag it as -1 and move on
      // Note that this->comm().sum might have left uncoarsenable
      // elements with error_per_parent=-n_proc, so reset it to
      // error_per_parent=-1
      if (error_per_parent[i] < 0.)
        {
          error_per_parent[i] = -1.;
          continue;
        }

      // The error estimator might have already given us an
      // estimate on the coarsenable parent elements; if so then
      // we want to retain that estimate
      if (error_per_cell[i])
        {
          error_per_parent[i] = error_per_cell[i];
          continue;
        }
      // if not, then e_parent = sqrt(sum(e_child^2))
      else
        error_per_parent[i] = std::sqrt(error_per_parent[i]);

      parent_error_min = std::min (parent_error_min,
                                   error_per_parent[i]);
      parent_error_max = std::max (parent_error_max,
                                   error_per_parent[i]);
    }
}



void MeshRefinement::update_nodes_map ()
{
  this->_new_nodes_map.init(_mesh);
}



bool MeshRefinement::test_level_one (bool libmesh_dbg_var(libmesh_assert_pass))
{
  // This function must be run on all processors at once
  parallel_object_only();

  // We may need a PointLocator for topological_neighbor() tests
  // later, which we need to make sure gets constructed on all
  // processors at once.
  std::unique_ptr<PointLocatorBase> point_locator;

#ifdef LIBMESH_ENABLE_PERIODIC
  bool has_periodic_boundaries =
    _periodic_boundaries && !_periodic_boundaries->empty();
  libmesh_assert(this->comm().verify(has_periodic_boundaries));

  if (has_periodic_boundaries)
    point_locator = _mesh.sub_point_locator();
#endif

  bool failure = false;

#ifndef NDEBUG
  Elem * failed_elem = nullptr;
  Elem * failed_neighbor = nullptr;
#endif // !NDEBUG

  for (auto & elem : _mesh.active_local_element_ptr_range())
    for (auto n : elem->side_index_range())
      {
        Elem * neighbor =
          topological_neighbor(elem, point_locator.get(), n);

        if (!neighbor || !neighbor->active() ||
            neighbor == remote_elem)
          continue;

        if ((neighbor->level() + 1 < elem->level()) ||
            (neighbor->p_level() + 1 < elem->p_level()) ||
            (neighbor->p_level() > elem->p_level() + 1))
          {
            failure = true;
#ifndef NDEBUG
            failed_elem = elem;
            failed_neighbor = neighbor;
#endif // !NDEBUG
            break;
          }
      }

  // If any processor failed, we failed globally
  this->comm().max(failure);

  if (failure)
    {
      // We didn't pass the level one test, so libmesh_assert that
      // we're allowed not to
#ifndef NDEBUG
      if (libmesh_assert_pass)
        {
          libMesh::out << "MeshRefinement Level one failure, element: "
                       << *failed_elem
                       << std::endl;
          libMesh::out << "MeshRefinement Level one failure, neighbor: "
                       << *failed_neighbor
                       << std::endl;
        }
#endif // !NDEBUG
      libmesh_assert(!libmesh_assert_pass);
      return false;
    }
  return true;
}



bool MeshRefinement::test_unflagged (bool libmesh_dbg_var(libmesh_assert_pass))
{
  // This function must be run on all processors at once
  parallel_object_only();

  bool found_flag = false;

#ifndef NDEBUG
  Elem * failed_elem = nullptr;
#endif

  // Search for local flags
  for (auto & elem : _mesh.active_local_element_ptr_range())
    if (elem->refinement_flag() == Elem::REFINE ||
        elem->refinement_flag() == Elem::COARSEN ||
        elem->p_refinement_flag() == Elem::REFINE ||
        elem->p_refinement_flag() == Elem::COARSEN)
      {
        found_flag = true;
#ifndef NDEBUG
        failed_elem = elem;
#endif
        break;
      }

  // If we found a flag on any processor, it counts
  this->comm().max(found_flag);

  if (found_flag)
    {
#ifndef NDEBUG
      if (libmesh_assert_pass)
        {
          libMesh::out <<
            "MeshRefinement test_unflagged failure, element: " <<
            *failed_elem << std::endl;
        }
#endif
      // We didn't pass the "elements are unflagged" test,
      // so libmesh_assert that we're allowed not to
      libmesh_assert(!libmesh_assert_pass);
      return false;
    }
  return true;
}



bool MeshRefinement::refine_and_coarsen_elements ()
{
  // This function must be run on all processors at once
  parallel_object_only();

  // We can't yet turn a non-level-one mesh into a level-one mesh
  if (_face_level_mismatch_limit)
    libmesh_assert(test_level_one(true));

  // Possibly clean up the refinement flags from
  // a previous step.  While we're at it, see if this method should be
  // a no-op.
  bool elements_flagged = false;

  for (auto & elem : _mesh.element_ptr_range())
    {
      // This might be left over from the last step
      const Elem::RefinementState flag = elem->refinement_flag();

      // Set refinement flag to INACTIVE if the
      // element isn't active
      if ( !elem->active())
        {
          elem->set_refinement_flag(Elem::INACTIVE);
          elem->set_p_refinement_flag(Elem::INACTIVE);
        }
      else if (flag == Elem::JUST_REFINED)
        elem->set_refinement_flag(Elem::DO_NOTHING);
      else if (!elements_flagged)
        {
          if (flag == Elem::REFINE || flag == Elem::COARSEN)
            elements_flagged = true;
          else
            {
              const Elem::RefinementState pflag =
                elem->p_refinement_flag();
              if (pflag == Elem::REFINE || pflag == Elem::COARSEN)
                elements_flagged = true;
            }
        }
    }

  // Did *any* processor find elements flagged for AMR/C?
  _mesh.comm().max(elements_flagged);

  // If we have nothing to do, let's not bother verifying that nothing
  // is compatible with nothing.
  if (!elements_flagged)
    return false;

  // Parallel consistency has to come first, or coarsening
  // along processor boundaries might occasionally be falsely
  // prevented
#ifdef DEBUG
  bool flags_were_consistent = this->make_flags_parallel_consistent();

  libmesh_assert (flags_were_consistent);
#endif

  // Smooth refinement and coarsening flags
  _smooth_flags(true, true);

  // First coarsen the flagged elements.
  const bool coarsening_changed_mesh =
    this->_coarsen_elements ();

  // First coarsen the flagged elements.
  // FIXME: test_level_one now tests consistency across periodic
  // boundaries, which requires a point_locator, which just got
  // invalidated by _coarsen_elements() and hasn't yet been cleared by
  // prepare_for_use().

  //  libmesh_assert(this->make_coarsening_compatible());
  //  libmesh_assert(this->make_refinement_compatible());

  // FIXME: This won't pass unless we add a redundant find_neighbors()
  // call or replace find_neighbors() with on-the-fly neighbor updating
  // libmesh_assert(!this->eliminate_unrefined_patches());

  // We can't contract the mesh ourselves anymore - a System might
  // need to restrict old coefficient vectors first
  // _mesh.contract();

  // First coarsen the flagged elements.
  // Now refine the flagged elements.  This will
  // take up some space, maybe more than what was freed.
  const bool refining_changed_mesh =
    this->_refine_elements();

  // First coarsen the flagged elements.
  // Finally, the new mesh needs to be prepared for use
  if (coarsening_changed_mesh || refining_changed_mesh)
    {
#ifdef DEBUG
      _mesh.libmesh_assert_valid_parallel_ids();
#endif

      _mesh.prepare_for_use (/*skip_renumber =*/false);

      if (_face_level_mismatch_limit)
        libmesh_assert(test_level_one(true));
      libmesh_assert(test_unflagged(true));
      libmesh_assert(this->make_coarsening_compatible());
      libmesh_assert(this->make_refinement_compatible());
      // FIXME: This won't pass unless we add a redundant find_neighbors()
      // call or replace find_neighbors() with on-the-fly neighbor updating
      // libmesh_assert(!this->eliminate_unrefined_patches());

      return true;
    }
  else
    {
      if (_face_level_mismatch_limit)
        libmesh_assert(test_level_one(true));
      libmesh_assert(test_unflagged(true));
      libmesh_assert(this->make_coarsening_compatible());
      libmesh_assert(this->make_refinement_compatible());
    }

  // Otherwise there was no change in the mesh,
  // let the user know.  Also, there is no need
  // to prepare the mesh for use since it did not change.
  return false;

}







bool MeshRefinement::coarsen_elements ()
{
  // This function must be run on all processors at once
  parallel_object_only();

  // We can't yet turn a non-level-one mesh into a level-one mesh
  if (_face_level_mismatch_limit)
    libmesh_assert(test_level_one(true));

  // Possibly clean up the refinement flags from
  // a previous step
  for (auto & elem : _mesh.element_ptr_range())
    {
      // Set refinement flag to INACTIVE if the
      // element isn't active
      if (!elem->active())
        {
          elem->set_refinement_flag(Elem::INACTIVE);
          elem->set_p_refinement_flag(Elem::INACTIVE);
        }

      // This might be left over from the last step
      if (elem->refinement_flag() == Elem::JUST_REFINED)
        elem->set_refinement_flag(Elem::DO_NOTHING);
    }

  // Parallel consistency has to come first, or coarsening
  // along processor boundaries might occasionally be falsely
  // prevented
  bool flags_were_consistent = this->make_flags_parallel_consistent();

  // In theory, we should be able to remove the above call, which can
  // be expensive and should be unnecessary.  In practice, doing
  // consistent flagging in parallel is hard, it's impossible to
  // verify at the library level if it's being done by user code, and
  // we don't want to abort large parallel runs in opt mode... but we
  // do want to warn that they should be fixed.
  libmesh_assert(flags_were_consistent);
  if (!flags_were_consistent)
    {
      libMesh::out << "Refinement flags were not consistent between processors!\n"
                   << "Correcting and continuing.";
    }

  // Smooth coarsening flags
  _smooth_flags(false, true);

  // Coarsen the flagged elements.
  const bool mesh_changed =
    this->_coarsen_elements ();

  if (_face_level_mismatch_limit)
    libmesh_assert(test_level_one(true));
  libmesh_assert(this->make_coarsening_compatible());
  // FIXME: This won't pass unless we add a redundant find_neighbors()
  // call or replace find_neighbors() with on-the-fly neighbor updating
  // libmesh_assert(!this->eliminate_unrefined_patches());

  // We can't contract the mesh ourselves anymore - a System might
  // need to restrict old coefficient vectors first
  // _mesh.contract();

  // Finally, the new mesh may need to be prepared for use
  if (mesh_changed)
    _mesh.prepare_for_use (/*skip_renumber =*/false);

  return mesh_changed;
}







bool MeshRefinement::refine_elements ()
{
  // This function must be run on all processors at once
  parallel_object_only();

  if (_face_level_mismatch_limit)
    libmesh_assert(test_level_one(true));

  // Possibly clean up the refinement flags from
  // a previous step
  for (auto & elem : _mesh.element_ptr_range())
    {
      // Set refinement flag to INACTIVE if the
      // element isn't active
      if (!elem->active())
        {
          elem->set_refinement_flag(Elem::INACTIVE);
          elem->set_p_refinement_flag(Elem::INACTIVE);
        }

      // This might be left over from the last step
      if (elem->refinement_flag() == Elem::JUST_REFINED)
        elem->set_refinement_flag(Elem::DO_NOTHING);
    }



  // Parallel consistency has to come first, or coarsening
  // along processor boundaries might occasionally be falsely
  // prevented
  bool flags_were_consistent = this->make_flags_parallel_consistent();

  // In theory, we should be able to remove the above call, which can
  // be expensive and should be unnecessary.  In practice, doing
  // consistent flagging in parallel is hard, it's impossible to
  // verify at the library level if it's being done by user code, and
  // we don't want to abort large parallel runs in opt mode... but we
  // do want to warn that they should be fixed.
  libmesh_assert(flags_were_consistent);
  if (!flags_were_consistent)
    {
      libMesh::out << "Refinement flags were not consistent between processors!\n"
                   << "Correcting and continuing.";
    }

  // Smooth refinement flags
  _smooth_flags(true, false);

  // Now refine the flagged elements.  This will
  // take up some space, maybe more than what was freed.
  const bool mesh_changed =
    this->_refine_elements();

  if (_face_level_mismatch_limit)
    libmesh_assert(test_level_one(true));
  libmesh_assert(this->make_refinement_compatible());
  // FIXME: This won't pass unless we add a redundant find_neighbors()
  // call or replace find_neighbors() with on-the-fly neighbor updating
  // libmesh_assert(!this->eliminate_unrefined_patches());

  // Finally, the new mesh needs to be prepared for use
  if (mesh_changed)
    _mesh.prepare_for_use (/*skip_renumber =*/false);

  return mesh_changed;
}



bool MeshRefinement::make_flags_parallel_consistent()
{
  // This function must be run on all processors at once
  parallel_object_only();

  LOG_SCOPE ("make_flags_parallel_consistent()", "MeshRefinement");

  SyncRefinementFlags hsync(_mesh, &Elem::refinement_flag,
                            &Elem::set_refinement_flag);
  Parallel::sync_dofobject_data_by_id
    (this->comm(), _mesh.elements_begin(), _mesh.elements_end(), hsync);

  SyncRefinementFlags psync(_mesh, &Elem::p_refinement_flag,
                            &Elem::set_p_refinement_flag);
  Parallel::sync_dofobject_data_by_id
    (this->comm(), _mesh.elements_begin(), _mesh.elements_end(), psync);

  // If we weren't consistent in both h and p on every processor then
  // we weren't globally consistent
  bool parallel_consistent = hsync.parallel_consistent &&
    psync.parallel_consistent;
  this->comm().min(parallel_consistent);

  return parallel_consistent;
}



bool MeshRefinement::make_coarsening_compatible()
{
  // This function must be run on all processors at once
  parallel_object_only();

  // We may need a PointLocator for topological_neighbor() tests
  // later, which we need to make sure gets constructed on all
  // processors at once.
  std::unique_ptr<PointLocatorBase> point_locator;

#ifdef LIBMESH_ENABLE_PERIODIC
  bool has_periodic_boundaries =
    _periodic_boundaries && !_periodic_boundaries->empty();
  libmesh_assert(this->comm().verify(has_periodic_boundaries));

  if (has_periodic_boundaries)
    point_locator = _mesh.sub_point_locator();
#endif

  LOG_SCOPE ("make_coarsening_compatible()", "MeshRefinement");

  // Unless we encounter a specific situation level-one
  // will be satisfied after executing this loop just once
  bool level_one_satisfied = true;


  // Unless we encounter a specific situation we will be compatible
  // with any selected refinement flags
  bool compatible_with_refinement = true;


  // find the maximum h and p levels in the mesh
  unsigned int max_level = 0;
  unsigned int max_p_level = 0;

  {
    // First we look at all the active level-0 elements.  Since it doesn't make
    // sense to coarsen them we must un-set their coarsen flags if
    // they are set.
    for (auto & elem : _mesh.active_element_ptr_range())
      {
        max_level = std::max(max_level, elem->level());
        max_p_level =
          std::max(max_p_level,
                   static_cast<unsigned int>(elem->p_level()));

        if ((elem->level() == 0) &&
            (elem->refinement_flag() == Elem::COARSEN))
          elem->set_refinement_flag(Elem::DO_NOTHING);

        if ((elem->p_level() == 0) &&
            (elem->p_refinement_flag() == Elem::COARSEN))
          elem->set_p_refinement_flag(Elem::DO_NOTHING);
      }
  }

  // Even if there are no refined elements on this processor then
  // there may still be work for us to do on e.g. ancestor elements.
  // At the very least we need to be in the loop if a distributed mesh
  // needs to synchronize data.
#if 0
  if (max_level == 0 && max_p_level == 0)
    {
      // But we still have to check with other processors
      this->comm().min(compatible_with_refinement);

      return compatible_with_refinement;
    }
#endif

  // Loop over all the active elements.  If an element is marked
  // for coarsening we better check its neighbors.  If ANY of these neighbors
  // are marked for refinement AND are at the same level then there is a
  // conflict.  By convention refinement wins, so we un-mark the element for
  // coarsening.  Level-one would be violated in this case so we need to re-run
  // the loop.
  if (_face_level_mismatch_limit)
    {

    repeat:
      level_one_satisfied = true;

      do
        {
          level_one_satisfied = true;

          for (auto & elem : _mesh.active_element_ptr_range())
            {
              bool my_flag_changed = false;

              if (elem->refinement_flag() == Elem::COARSEN) // If the element is active and
                // the coarsen flag is set
                {
                  const unsigned int my_level = elem->level();

                  for (auto n : elem->side_index_range())
                    {
                      const Elem * neighbor =
                        topological_neighbor(elem, point_locator.get(), n);

                      if (neighbor != nullptr &&      // I have a
                          neighbor != remote_elem) // neighbor here
                        {
                          if (neighbor->active()) // and it is active
                            {
                              if ((neighbor->level() == my_level) &&
                                  (neighbor->refinement_flag() == Elem::REFINE)) // the neighbor is at my level
                                // and wants to be refined
                                {
                                  elem->set_refinement_flag(Elem::DO_NOTHING);
                                  my_flag_changed = true;
                                  break;
                                }
                            }
                          else // I have a neighbor and it is not active. That means it has children.
                            {  // While it _may_ be possible to coarsen us if all the children of
                              // that element want to be coarsened, it is impossible to know at this
                              // stage.  Forget about it for the moment...  This can be handled in
                              // two steps.
                              elem->set_refinement_flag(Elem::DO_NOTHING);
                              my_flag_changed = true;
                              break;
                            }
                        }
                    }
                }
              if (elem->p_refinement_flag() == Elem::COARSEN) // If
                // the element is active and the order reduction flag is set
                {
                  const unsigned int my_p_level = elem->p_level();

                  for (auto n : elem->side_index_range())
                    {
                      const Elem * neighbor =
                        topological_neighbor(elem, point_locator.get(), n);

                      if (neighbor != nullptr &&      // I have a
                          neighbor != remote_elem) // neighbor here
                        {
                          if (neighbor->active()) // and it is active
                            {
                              if ((neighbor->p_level() > my_p_level &&
                                   neighbor->p_refinement_flag() != Elem::COARSEN)
                                  || (neighbor->p_level() == my_p_level &&
                                      neighbor->p_refinement_flag() == Elem::REFINE))
                                {
                                  elem->set_p_refinement_flag(Elem::DO_NOTHING);
                                  my_flag_changed = true;
                                  break;
                                }
                            }
                          else // I have a neighbor and it is not active.
                            {  // We need to find which of its children
                              // have me as a neighbor, and maintain
                              // level one p compatibility with them.
                              // Because we currently have level one h
                              // compatibility, we don't need to check
                              // grandchildren

                              libmesh_assert(neighbor->has_children());
                              for (auto & subneighbor : neighbor->child_ref_range())
                                if (&subneighbor != remote_elem &&
                                    subneighbor.active() &&
                                    has_topological_neighbor(&subneighbor, point_locator.get(), elem))
                                  if ((subneighbor.p_level() > my_p_level &&
                                       subneighbor.p_refinement_flag() != Elem::COARSEN)
                                      || (subneighbor.p_level() == my_p_level &&
                                          subneighbor.p_refinement_flag() == Elem::REFINE))
                                    {
                                      elem->set_p_refinement_flag(Elem::DO_NOTHING);
                                      my_flag_changed = true;
                                      break;
                                    }
                              if (my_flag_changed)
                                break;
                            }
                        }
                    }
                }

              // If the current element's flag changed, we hadn't
              // satisfied the level one rule.
              if (my_flag_changed)
                level_one_satisfied = false;

              // Additionally, if it has non-local neighbors, and
              // we're not in serial, then we'll eventually have to
              // return compatible_with_refinement = false, because
              // our change has to propagate to neighboring
              // processors.
              if (my_flag_changed && !_mesh.is_serial())
                for (auto n : elem->side_index_range())
                  {
                    Elem * neigh =
                      topological_neighbor(elem, point_locator.get(), n);

                    if (!neigh)
                      continue;
                    if (neigh == remote_elem ||
                        neigh->processor_id() !=
                        this->processor_id())
                      {
                        compatible_with_refinement = false;
                        break;
                      }
                    // FIXME - for non-level one meshes we should
                    // test all descendants
                    if (neigh->has_children())
                      for (auto & child : neigh->child_ref_range())
                        if (&child == remote_elem ||
                            child.processor_id() !=
                            this->processor_id())
                          {
                            compatible_with_refinement = false;
                            break;
                          }
                  }
            }
        }
      while (!level_one_satisfied);

    } // end if (_face_level_mismatch_limit)


  // Next we look at all of the ancestor cells.
  // If there is a parent cell with all of its children
  // wanting to be unrefined then the element is a candidate
  // for unrefinement.  If all the children don't
  // all want to be unrefined then ALL of them need to have their
  // unrefinement flags cleared.
  for (int level = max_level; level >= 0; level--)
    for (auto & elem : as_range(_mesh.level_elements_begin(level), _mesh.level_elements_end(level)))
      if (elem->ancestor())
        {
          // right now the element hasn't been disqualified
          // as a candidate for unrefinement
          bool is_a_candidate = true;
          bool found_remote_child = false;

          for (auto & child : elem->child_ref_range())
            {
              if (&child == remote_elem)
                found_remote_child = true;
              else if ((child.refinement_flag() != Elem::COARSEN) ||
                       !child.active() )
                is_a_candidate = false;
            }

          if (!is_a_candidate && !found_remote_child)
            {
              elem->set_refinement_flag(Elem::INACTIVE);

              for (auto & child : elem->child_ref_range())
                {
                  if (&child == remote_elem)
                    continue;
                  if (child.refinement_flag() == Elem::COARSEN)
                    {
                      level_one_satisfied = false;
                      child.set_refinement_flag(Elem::DO_NOTHING);
                    }
                }
            }
        }

  if (!level_one_satisfied && _face_level_mismatch_limit) goto repeat;


  // If all the children of a parent are set to be coarsened
  // then flag the parent so that they can kill their kids.

  // On a distributed mesh, we won't always be able to determine this
  // on parent elements with remote children, even if we own the
  // parent, without communication.
  //
  // We'll first communicate *to* parents' owners when we determine
  // they cannot be coarsened, then we'll sync the final refinement
  // flag *from* the parents.

  // uncoarsenable_parents[p] live on processor id p
  const processor_id_type n_proc     = _mesh.n_processors();
  const processor_id_type my_proc_id = _mesh.processor_id();
  const bool distributed_mesh = !_mesh.is_replicated();

  std::vector<std::vector<dof_id_type>>
    uncoarsenable_parents(n_proc);

  for (auto & elem : as_range(_mesh.ancestor_elements_begin(), _mesh.ancestor_elements_end()))
    {
      // Presume all the children are flagged for coarsening and
      // then look for a contradiction
      bool all_children_flagged_for_coarsening = true;

      for (auto & child : elem->child_ref_range())
        {
          if (&child != remote_elem &&
              child.refinement_flag() != Elem::COARSEN)
            {
              all_children_flagged_for_coarsening = false;
              if (!distributed_mesh)
                break;
              if (child.processor_id() != elem->processor_id())
                {
                  uncoarsenable_parents[elem->processor_id()].push_back(elem->id());
                  break;
                }
            }
        }

      if (all_children_flagged_for_coarsening)
        elem->set_refinement_flag(Elem::COARSEN_INACTIVE);
      else
        elem->set_refinement_flag(Elem::INACTIVE);
    }

  // If we have a distributed mesh, we might need to sync up
  // INACTIVE vs. COARSEN_INACTIVE flags.
  if (distributed_mesh)
    {
      // We'd better still be in sync here
      parallel_object_only();

      Parallel::MessageTag
        uncoarsenable_tag = this->comm().get_unique_tag();
      std::vector<Parallel::Request> uncoarsenable_push_requests(n_proc-1);

      for (processor_id_type p = 0; p != n_proc; ++p)
        {
          if (p == my_proc_id)
            continue;

          Parallel::Request &request =
            uncoarsenable_push_requests[p - (p > my_proc_id)];

          _mesh.comm().send
            (p, uncoarsenable_parents[p], request, uncoarsenable_tag);
        }

      for (processor_id_type p = 1; p != n_proc; ++p)
        {
          std::vector<dof_id_type> my_uncoarsenable_parents;
          _mesh.comm().receive
            (Parallel::any_source, my_uncoarsenable_parents,
             uncoarsenable_tag);

          for (const auto & id : my_uncoarsenable_parents)
            {
              Elem & elem = _mesh.elem_ref(id);
              libmesh_assert(elem.refinement_flag() == Elem::INACTIVE ||
                             elem.refinement_flag() == Elem::COARSEN_INACTIVE);
              elem.set_refinement_flag(Elem::INACTIVE);
            }
        }

      Parallel::wait(uncoarsenable_push_requests);

      SyncRefinementFlags hsync(_mesh, &Elem::refinement_flag,
                                &Elem::set_refinement_flag);
      sync_dofobject_data_by_id
        (this->comm(), _mesh.not_local_elements_begin(),
         _mesh.not_local_elements_end(),
         // We'd like a smaller sync, but this leads to bugs?
         // SyncCoarsenInactive(),
         hsync);
    }

  // If one processor finds an incompatibility, we're globally
  // incompatible
  this->comm().min(compatible_with_refinement);

  return compatible_with_refinement;
}








bool MeshRefinement::make_refinement_compatible()
{
  // This function must be run on all processors at once
  parallel_object_only();

  // We may need a PointLocator for topological_neighbor() tests
  // later, which we need to make sure gets constructed on all
  // processors at once.
  std::unique_ptr<PointLocatorBase> point_locator;

#ifdef LIBMESH_ENABLE_PERIODIC
  bool has_periodic_boundaries =
    _periodic_boundaries && !_periodic_boundaries->empty();
  libmesh_assert(this->comm().verify(has_periodic_boundaries));

  if (has_periodic_boundaries)
    point_locator = _mesh.sub_point_locator();
#endif

  LOG_SCOPE ("make_refinement_compatible()", "MeshRefinement");

  // Unless we encounter a specific situation we will be compatible
  // with any selected coarsening flags
  bool compatible_with_coarsening = true;

  // This loop enforces the level-1 rule.  We should only
  // execute it if the user indeed wants level-1 satisfied!
  if (_face_level_mismatch_limit)
    {
      // Unless we encounter a specific situation level-one
      // will be satisfied after executing this loop just once
      bool level_one_satisfied = true;

      do
        {
          level_one_satisfied = true;

          for (auto & elem : _mesh.active_element_ptr_range())
            {
              const unsigned short n_sides = elem->n_sides();

              if (elem->refinement_flag() == Elem::REFINE)  // If the element is active and the
                // h refinement flag is set
                {
                  const unsigned int my_level = elem->level();

                  for (unsigned short side = 0; side != n_sides;
                       ++side)
                    {
                      Elem * neighbor =
                        topological_neighbor(elem, point_locator.get(), side);

                      if (neighbor != nullptr        && // I have a
                          neighbor != remote_elem && // neighbor here
                          neighbor->active()) // and it is active
                        {
                          // Case 1:  The neighbor is at the same level I am.
                          //        1a: The neighbor will be refined       -> NO PROBLEM
                          //        1b: The neighbor won't be refined      -> NO PROBLEM
                          //        1c: The neighbor wants to be coarsened -> PROBLEM
                          if (neighbor->level() == my_level)
                            {
                              if (neighbor->refinement_flag() == Elem::COARSEN)
                                {
                                  neighbor->set_refinement_flag(Elem::DO_NOTHING);
                                  if (neighbor->parent())
                                    neighbor->parent()->set_refinement_flag(Elem::INACTIVE);
                                  compatible_with_coarsening = false;
                                  level_one_satisfied = false;
                                }
                            }


                          // Case 2: The neighbor is one level lower than I am.
                          //         The neighbor thus MUST be refined to satisfy
                          //         the level-one rule, regardless of whether it
                          //         was originally flagged for refinement. If it
                          //         wasn't flagged already we need to repeat
                          //         this process.
                          else if ((neighbor->level()+1) == my_level)
                            {
                              if (neighbor->refinement_flag() != Elem::REFINE)
                                {
                                  neighbor->set_refinement_flag(Elem::REFINE);
                                  if (neighbor->parent())
                                    neighbor->parent()->set_refinement_flag(Elem::INACTIVE);
                                  compatible_with_coarsening = false;
                                  level_one_satisfied = false;
                                }
                            }
#ifdef DEBUG
                          // Note that the only other possibility is that the
                          // neighbor is already refined, in which case it isn't
                          // active and we should never get here.
                          else
                            libmesh_error_msg("ERROR: Neighbor level must be equal or 1 higher than mine.");
#endif
                        }
                    }
                }
              if (elem->p_refinement_flag() == Elem::REFINE)  // If the element is active and the
                // p refinement flag is set
                {
                  const unsigned int my_p_level = elem->p_level();

                  for (unsigned int side=0; side != n_sides; side++)
                    {
                      Elem * neighbor =
                        topological_neighbor(elem, point_locator.get(), side);

                      if (neighbor != nullptr &&      // I have a
                          neighbor != remote_elem) // neighbor here
                        {
                          if (neighbor->active()) // and it is active
                            {
                              if (neighbor->p_level() < my_p_level &&
                                  neighbor->p_refinement_flag() != Elem::REFINE)
                                {
                                  neighbor->set_p_refinement_flag(Elem::REFINE);
                                  level_one_satisfied = false;
                                  compatible_with_coarsening = false;
                                }
                              if (neighbor->p_level() == my_p_level &&
                                  neighbor->p_refinement_flag() == Elem::COARSEN)
                                {
                                  neighbor->set_p_refinement_flag(Elem::DO_NOTHING);
                                  level_one_satisfied = false;
                                  compatible_with_coarsening = false;
                                }
                            }
                          else // I have an inactive neighbor
                            {
                              libmesh_assert(neighbor->has_children());
                              for (auto & subneighbor : neighbor->child_ref_range())
                                if (&subneighbor != remote_elem && subneighbor.active() &&
                                    has_topological_neighbor(&subneighbor, point_locator.get(), elem))
                                  {
                                    if (subneighbor.p_level() < my_p_level &&
                                        subneighbor.p_refinement_flag() != Elem::REFINE)
                                      {
                                        // We should already be level one
                                        // compatible
                                        libmesh_assert_greater (subneighbor.p_level() + 2u,
                                                                my_p_level);
                                        subneighbor.set_p_refinement_flag(Elem::REFINE);
                                        level_one_satisfied = false;
                                        compatible_with_coarsening = false;
                                      }
                                    if (subneighbor.p_level() == my_p_level &&
                                        subneighbor.p_refinement_flag() == Elem::COARSEN)
                                      {
                                        subneighbor.set_p_refinement_flag(Elem::DO_NOTHING);
                                        level_one_satisfied = false;
                                        compatible_with_coarsening = false;
                                      }
                                  }
                            }
                        }
                    }
                }
            }
        }

      while (!level_one_satisfied);
    } // end if (_face_level_mismatch_limit)

  // If we're not compatible on one processor, we're globally not
  // compatible
  this->comm().min(compatible_with_coarsening);

  return compatible_with_coarsening;
}




bool MeshRefinement::_coarsen_elements ()
{
  // This function must be run on all processors at once
  parallel_object_only();

  LOG_SCOPE ("_coarsen_elements()", "MeshRefinement");

  // Flags indicating if this call actually changes the mesh
  bool mesh_changed = false;
  bool mesh_p_changed = false;

  // Clear the unused_elements data structure.
  // The elements have been packed since it was built,
  // so there are _no_ unused elements.  We cannot trust
  // any iterators currently in this data structure.
  // _unused_elements.clear();

  // Loop over the elements first to determine if the mesh will
  // undergo h-coarsening.  If it will, then we'll need to communicate
  // more ghosted elements.  We need to communicate them *before* we
  // do the coarsening; otherwise it is possible to coarsen away a
  // one-element-thick layer partition and leave the partitions on
  // either side unable to figure out how to talk to each other.
  for (auto & elem : _mesh.element_ptr_range())
    if (elem->refinement_flag() == Elem::COARSEN)
      {
        mesh_changed = true;
        break;
      }

  // If the mesh changed on any processor, it changed globally
  this->comm().max(mesh_changed);

  // And then we may need to widen the ghosting layers.
  if (mesh_changed)
    MeshCommunication().send_coarse_ghosts(_mesh);

  for (auto & elem : _mesh.element_ptr_range())
    {
      // active elements flagged for coarsening will
      // no longer be deleted until MeshRefinement::contract()
      if (elem->refinement_flag() == Elem::COARSEN)
        {
          // Huh?  no level-0 element should be active
          // and flagged for coarsening.
          libmesh_assert_not_equal_to (elem->level(), 0);

          // Remove this element from any neighbor
          // lists that point to it.
          elem->nullify_neighbors();

          // Remove any boundary information associated
          // with this element
          _mesh.get_boundary_info().remove (elem);

          // Add this iterator to the _unused_elements
          // data structure so we might fill it.
          // The _unused_elements optimization is currently off.
          // _unused_elements.push_back (it);

          // Don't delete the element until
          // MeshRefinement::contract()
          // _mesh.delete_elem(elem);
        }

      // inactive elements flagged for coarsening
      // will become active
      else if (elem->refinement_flag() == Elem::COARSEN_INACTIVE)
        {
          elem->coarsen();
          libmesh_assert (elem->active());

          // the mesh has certainly changed
          mesh_changed = true;
        }
      if (elem->p_refinement_flag() == Elem::COARSEN)
        {
          if (elem->p_level() > 0)
            {
              elem->set_p_refinement_flag(Elem::JUST_COARSENED);
              elem->set_p_level(elem->p_level() - 1);
              mesh_p_changed = true;
            }
          else
            {
              elem->set_p_refinement_flag(Elem::DO_NOTHING);
            }
        }
    }

  this->comm().max(mesh_p_changed);

  // And we may need to update DistributedMesh values reflecting the changes
  if (mesh_changed)
    _mesh.update_parallel_id_counts();

  // Node processor ids may need to change if an element of that id
  // was coarsened away
  if (mesh_changed && !_mesh.is_serial())
    {
      // Update the _new_nodes_map so that processors can
      // find requested nodes
      this->update_nodes_map ();

      MeshCommunication().make_nodes_parallel_consistent (_mesh);

      // Clear the _new_nodes_map
      this->clear();

#ifdef DEBUG
      MeshTools::libmesh_assert_valid_procids<Node>(_mesh);
#endif
    }

  // If p levels changed all we need to do is make sure that parent p
  // levels changed in sync
  if (mesh_p_changed && !_mesh.is_serial())
    {
      MeshCommunication().make_p_levels_parallel_consistent (_mesh);
    }

  return (mesh_changed || mesh_p_changed);
}



bool MeshRefinement::_refine_elements ()
{
  libmesh_assert(_mesh.is_prepared() || _mesh.is_replicated());

  // This function must be run on all processors at once
  parallel_object_only();

  // Update the _new_nodes_map so that elements can
  // find nodes to connect to.
  this->update_nodes_map ();

  LOG_SCOPE ("_refine_elements()", "MeshRefinement");

  // Iterate over the elements, counting the elements
  // flagged for h refinement.
  dof_id_type n_elems_flagged = 0;

  for (auto & elem : _mesh.element_ptr_range())
    if (elem->refinement_flag() == Elem::REFINE)
      n_elems_flagged++;

  // Construct a local vector of Elem * which have been
  // previously marked for refinement.  We reserve enough
  // space to allow for every element to be refined.
  std::vector<Elem *> local_copy_of_elements;
  local_copy_of_elements.reserve(n_elems_flagged);

  // If mesh p levels changed, we might need to synchronize parent p
  // levels on a distributed mesh.
  bool mesh_p_changed = false;

  // Iterate over the elements, looking for elements flagged for
  // refinement.

  // If we are on a ReplicatedMesh, then we just do the refinement in
  // the same order on every processor and everything stays in sync.

  // If we are on a DistributedMesh, that's impossible.
  //
  // If the mesh is distributed, we need to make sure that if we end
  // up as the owner of a new node, which might happen if that node is
  // attached to one of our own elements, then we have given it a
  // legitimate node id and our own processor id.  We generate
  // legitimate node ids and use our own processor id when we are
  // refining our own elements but not when we refine others'
  // elements.  Therefore we want to refine our own elements *first*,
  // thereby generating all nodes which might belong to us, and then
  // refine others' elements *after*, thereby generating nodes with
  // temporary ids which we know we will discard.
  //
  // Even if the DistributedMesh is serialized, we can't just treat it
  // like a ReplicatedMesh, because DistributedMesh doesn't *trust*
  // users to refine partitioned elements in a serialized way, so it
  // assigns temporary ids, so we need to synchronize ids afterward to
  // be safe anyway, so we might as well use the distributed mesh code
  // path.
  for (auto & elem : _mesh.is_replicated() ? _mesh.active_element_ptr_range() : _mesh.active_local_element_ptr_range())
    {
      if (elem->refinement_flag() == Elem::REFINE)
        local_copy_of_elements.push_back(elem);
      if (elem->p_refinement_flag() == Elem::REFINE &&
          elem->active())
        {
          elem->set_p_level(elem->p_level()+1);
          elem->set_p_refinement_flag(Elem::JUST_REFINED);
          mesh_p_changed = true;
        }
    }

  if (!_mesh.is_replicated())
    {
      for (auto & elem : as_range(_mesh.active_not_local_elements_begin(),
                                  _mesh.active_not_local_elements_end()))
        {
          if (elem->refinement_flag() == Elem::REFINE)
            local_copy_of_elements.push_back(elem);
          if (elem->p_refinement_flag() == Elem::REFINE &&
              elem->active())
            {
              elem->set_p_level(elem->p_level()+1);
              elem->set_p_refinement_flag(Elem::JUST_REFINED);
              mesh_p_changed = true;
            }
        }
    }

  // Now iterate over the local copies and refine each one.
  // This may resize the mesh's internal container and invalidate
  // any existing iterators.
  for (auto & elem : local_copy_of_elements)
    elem->refine(*this);

  // The mesh changed if there were elements h refined
  bool mesh_changed = !local_copy_of_elements.empty();

  // If the mesh changed on any processor, it changed globally
  this->comm().max(mesh_changed);
  this->comm().max(mesh_p_changed);

  // And we may need to update DistributedMesh values reflecting the changes
  if (mesh_changed)
    _mesh.update_parallel_id_counts();

  if (mesh_changed && !_mesh.is_replicated())
    {
      MeshCommunication().make_elems_parallel_consistent (_mesh);
      MeshCommunication().make_new_nodes_parallel_consistent (_mesh);
#ifdef DEBUG
      _mesh.libmesh_assert_valid_parallel_ids();
#endif
    }

  // If we're refining a ReplicatedMesh, then we haven't yet assigned
  // node processor ids.  But if we're refining a partitioned
  // ReplicatedMesh, then we *need* to assign node processor ids.
  if (mesh_changed && _mesh.is_replicated() &&
      (_mesh.unpartitioned_elements_begin() ==
       _mesh.unpartitioned_elements_end()))
    Partitioner::set_node_processor_ids(_mesh);

  if (mesh_p_changed && !_mesh.is_replicated())
    {
      MeshCommunication().make_p_levels_parallel_consistent (_mesh);
    }

  // Clear the _new_nodes_map and _unused_elements data structures.
  this->clear();

  return (mesh_changed || mesh_p_changed);
}


void MeshRefinement::_smooth_flags(bool refining, bool coarsening)
{
  // Smoothing can break in weird ways on a mesh with broken topology
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_neighbors(_mesh);
#endif

  // Repeat until flag changes match on every processor
  do
    {
      // Repeat until coarsening & refinement flags jive
      bool satisfied = false;
      do
        {
          // If we're refining or coarsening, hit the corresponding
          // face level test code.  Short-circuiting || is our friend
          const bool coarsening_satisfied =
            !coarsening ||
            this->make_coarsening_compatible();

          const bool refinement_satisfied =
            !refining ||
            this->make_refinement_compatible();

          bool smoothing_satisfied =
            !this->eliminate_unrefined_patches();// &&

          if (_edge_level_mismatch_limit)
            smoothing_satisfied = smoothing_satisfied &&
              !this->limit_level_mismatch_at_edge (_edge_level_mismatch_limit);

          if (_node_level_mismatch_limit)
            smoothing_satisfied = smoothing_satisfied &&
              !this->limit_level_mismatch_at_node (_node_level_mismatch_limit);

          if (_overrefined_boundary_limit>=0)
            smoothing_satisfied = smoothing_satisfied &&
              !this->limit_overrefined_boundary(_overrefined_boundary_limit);

          if (_underrefined_boundary_limit>=0)
            smoothing_satisfied = smoothing_satisfied &&
              !this->limit_underrefined_boundary(_underrefined_boundary_limit);

          satisfied = (coarsening_satisfied &&
                       refinement_satisfied &&
                       smoothing_satisfied);

          libmesh_assert(this->comm().verify(satisfied));
        }
      while (!satisfied);
    }
  while (!_mesh.is_serial() && !this->make_flags_parallel_consistent());
}


void MeshRefinement::uniformly_p_refine (unsigned int n)
{
  // Refine n times
  for (unsigned int rstep=0; rstep<n; rstep++)
    for (auto & elem : _mesh.active_element_ptr_range())
      {
        // P refine all the active elements
        elem->set_p_level(elem->p_level()+1);
        elem->set_p_refinement_flag(Elem::JUST_REFINED);
      }
}



void MeshRefinement::uniformly_p_coarsen (unsigned int n)
{
  // Coarsen p times
  for (unsigned int rstep=0; rstep<n; rstep++)
    for (auto & elem : _mesh.active_element_ptr_range())
      if (elem->p_level() > 0)
        {
          // P coarsen all the active elements
          elem->set_p_level(elem->p_level()-1);
          elem->set_p_refinement_flag(Elem::JUST_COARSENED);
        }
}



void MeshRefinement::uniformly_refine (unsigned int n)
{
  // Refine n times
  // FIXME - this won't work if n>1 and the mesh
  // has already been attached to an equation system
  for (unsigned int rstep=0; rstep<n; rstep++)
    {
      // Clean up the refinement flags
      this->clean_refinement_flags();

      // Flag all the active elements for refinement.
      for (auto & elem : _mesh.active_element_ptr_range())
        elem->set_refinement_flag(Elem::REFINE);

      // Refine all the elements we just flagged.
      this->_refine_elements();
    }

  // Finally, the new mesh probably needs to be prepared for use
  if (n > 0)
    _mesh.prepare_for_use (/*skip_renumber =*/false);
}



void MeshRefinement::uniformly_coarsen (unsigned int n)
{
  // Coarsen n times
  for (unsigned int rstep=0; rstep<n; rstep++)
    {
      // Clean up the refinement flags
      this->clean_refinement_flags();

      // Flag all the active elements for coarsening.
      for (auto & elem : _mesh.active_element_ptr_range())
        {
          elem->set_refinement_flag(Elem::COARSEN);
          if (elem->parent())
            elem->parent()->set_refinement_flag(Elem::COARSEN_INACTIVE);
        }

      // On a distributed mesh, we may have parent elements with
      // remote active children.  To keep flags consistent, we'll need
      // a communication step.
      if (!_mesh.is_replicated())
        {
          const processor_id_type n_proc = _mesh.n_processors();
          const processor_id_type my_proc_id = _mesh.processor_id();

          std::vector<std::vector<dof_id_type>>
            parents_to_coarsen(n_proc);

          for (const auto & elem : as_range(_mesh.ancestor_elements_begin(), _mesh.ancestor_elements_end()))
            if (elem->processor_id() != my_proc_id &&
                elem->refinement_flag() == Elem::COARSEN_INACTIVE)
              parents_to_coarsen[elem->processor_id()].push_back(elem->id());

          Parallel::MessageTag
            coarsen_tag = this->comm().get_unique_tag();
          std::vector<Parallel::Request> coarsen_push_requests(n_proc-1);

          for (processor_id_type p = 0; p != n_proc; ++p)
            {
              if (p == my_proc_id)
                continue;

              Parallel::Request &request =
                coarsen_push_requests[p - (p > my_proc_id)];

              _mesh.comm().send
                (p, parents_to_coarsen[p], request, coarsen_tag);
            }

          for (processor_id_type p = 1; p != n_proc; ++p)
            {
              std::vector<dof_id_type> my_parents_to_coarsen;
              _mesh.comm().receive
                (Parallel::any_source, my_parents_to_coarsen,
                 coarsen_tag);

              for (const auto & id : my_parents_to_coarsen)
                {
                  Elem & elem = _mesh.elem_ref(id);
                  libmesh_assert(elem.refinement_flag() == Elem::INACTIVE ||
                                 elem.refinement_flag() == Elem::COARSEN_INACTIVE);
                  elem.set_refinement_flag(Elem::COARSEN_INACTIVE);
                }
            }

          Parallel::wait(coarsen_push_requests);

          SyncRefinementFlags hsync(_mesh, &Elem::refinement_flag,
                                    &Elem::set_refinement_flag);
          sync_dofobject_data_by_id
            (this->comm(), _mesh.not_local_elements_begin(),
             _mesh.not_local_elements_end(),
             // We'd like a smaller sync, but this leads to bugs?
             // SyncCoarsenInactive(),
             hsync);
        }

      // Coarsen all the elements we just flagged.
      this->_coarsen_elements();
    }


  // Finally, the new mesh probably needs to be prepared for use
  if (n > 0)
    _mesh.prepare_for_use (/*skip_renumber =*/false);
}



Elem * MeshRefinement::topological_neighbor(Elem * elem,
                                            const PointLocatorBase * point_locator,
                                            const unsigned int side)
{
#ifdef LIBMESH_ENABLE_PERIODIC
  if (_periodic_boundaries && !_periodic_boundaries->empty())
    {
      libmesh_assert(point_locator);
      return elem->topological_neighbor(side, _mesh, *point_locator, _periodic_boundaries);
    }
#endif
  return elem->neighbor_ptr(side);
}



bool MeshRefinement::has_topological_neighbor(const Elem * elem,
                                              const PointLocatorBase * point_locator,
                                              const Elem * neighbor)
{
#ifdef LIBMESH_ENABLE_PERIODIC
  if (_periodic_boundaries && !_periodic_boundaries->empty())
    {
      libmesh_assert(point_locator);
      return elem->has_topological_neighbor(neighbor, _mesh, *point_locator, _periodic_boundaries);
    }
#endif
  return elem->has_neighbor(neighbor);
}



} // namespace libMesh


#endif
