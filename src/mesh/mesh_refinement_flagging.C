
// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/libmesh_config.h"

// only compile these functions if the user requests AMR support
#ifdef LIBMESH_ENABLE_AMR

// C++ includes
#include <algorithm> // for std::sort

// Local includes
#include "libmesh/elem.h"
#include "libmesh/error_vector.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"
#include "libmesh/remote_elem.h"

namespace libMesh
{



//-----------------------------------------------------------------
// Mesh refinement methods
void MeshRefinement::flag_elements_by_error_fraction (const ErrorVector & error_per_cell,
                                                      const Real refine_frac,
                                                      const Real coarsen_frac,
                                                      const unsigned int max_l)
{
  parallel_object_only();

  // The function arguments are currently just there for
  // backwards_compatibility
  if (!_use_member_parameters)
    {
      // If the user used non-default parameters, lets warn
      // that they're deprecated
      if (refine_frac != 0.3 ||
          coarsen_frac != 0.0 ||
          max_l != libMesh::invalid_uint)
        libmesh_deprecated();

      _refine_fraction = refine_frac;
      _coarsen_fraction = coarsen_frac;
      _max_h_level = max_l;
    }

  // Check for valid fractions..
  // The fraction values must be in [0,1]
  libmesh_assert_greater_equal (_refine_fraction, 0);
  libmesh_assert_less_equal (_refine_fraction, 1);
  libmesh_assert_greater_equal (_coarsen_fraction, 0);
  libmesh_assert_less_equal (_coarsen_fraction, 1);

  // Clean up the refinement flags.  These could be left
  // over from previous refinement steps.
  this->clean_refinement_flags();

  // We're getting the minimum and maximum error values
  // for the ACTIVE elements
  Real error_min = 1.e30;
  Real error_max = 0.;

  // And, if necessary, for their parents
  Real parent_error_min = 1.e30;
  Real parent_error_max = 0.;

  // Prepare another error vector if we need to sum parent errors
  ErrorVector error_per_parent;
  if (_coarsen_by_parents)
    {
      create_parent_error_vector(error_per_cell,
                                 error_per_parent,
                                 parent_error_min,
                                 parent_error_max);
    }

  // We need to loop over all active elements to find the minimum
  MeshBase::element_iterator       el_it  =
    _mesh.active_local_elements_begin();
  const MeshBase::element_iterator el_end =
    _mesh.active_local_elements_end();

  for (; el_it != el_end; ++el_it)
    {
      const dof_id_type id  = (*el_it)->id();
      libmesh_assert_less (id, error_per_cell.size());

      error_max = std::max (error_max, error_per_cell[id]);
      error_min = std::min (error_min, error_per_cell[id]);
    }
  this->comm().max(error_max);
  this->comm().min(error_min);

  // Compute the cutoff values for coarsening and refinement
  const Real error_delta = (error_max - error_min);
  const Real parent_error_delta = parent_error_max - parent_error_min;

  const Real refine_cutoff  = (1.- _refine_fraction)*error_max;
  const Real coarsen_cutoff = _coarsen_fraction*error_delta + error_min;
  const Real parent_cutoff = _coarsen_fraction*parent_error_delta + error_min;

  //   // Print information about the error
  //   libMesh::out << " Error Information:"                     << std::endl
  //     << " ------------------"                     << std::endl
  //     << "   min:              " << error_min      << std::endl
  //     << "   max:              " << error_max      << std::endl
  //     << "   delta:            " << error_delta    << std::endl
  //     << "     refine_cutoff:  " << refine_cutoff  << std::endl
  //     << "     coarsen_cutoff: " << coarsen_cutoff << std::endl;



  // Loop over the elements and flag them for coarsening or
  // refinement based on the element error

  MeshBase::element_iterator       e_it  =
    _mesh.active_elements_begin();
  const MeshBase::element_iterator e_end =
    _mesh.active_elements_end();
  for (; e_it != e_end; ++e_it)
    {
      Elem * elem           = *e_it;
      const dof_id_type id = elem->id();

      libmesh_assert_less (id, error_per_cell.size());

      const ErrorVectorReal elem_error = error_per_cell[id];

      if (_coarsen_by_parents)
        {
          Elem * parent           = elem->parent();
          if (parent)
            {
              const dof_id_type parentid  = parent->id();
              if (error_per_parent[parentid] >= 0. &&
                  error_per_parent[parentid] <= parent_cutoff)
                elem->set_refinement_flag(Elem::COARSEN);
            }
        }
      // Flag the element for coarsening if its error
      // is <= coarsen_fraction*delta + error_min
      else if (elem_error <= coarsen_cutoff)
        {
          elem->set_refinement_flag(Elem::COARSEN);
        }

      // Flag the element for refinement if its error
      // is >= refinement_cutoff.
      if (elem_error >= refine_cutoff)
        if (elem->level() < _max_h_level)
          elem->set_refinement_flag(Elem::REFINE);
    }
}



void MeshRefinement::flag_elements_by_error_tolerance (const ErrorVector & error_per_cell_in)
{
  parallel_object_only();

  libmesh_assert_greater (_coarsen_threshold, 0);

  // Check for valid fractions..
  // The fraction values must be in [0,1]
  libmesh_assert_greater_equal (_refine_fraction, 0);
  libmesh_assert_less_equal (_refine_fraction, 1);
  libmesh_assert_greater_equal (_coarsen_fraction, 0);
  libmesh_assert_less_equal (_coarsen_fraction, 1);

  // How much error per cell will we tolerate?
  const Real local_refinement_tolerance =
    _absolute_global_tolerance / std::sqrt(static_cast<Real>(_mesh.n_active_elem()));
  const Real local_coarsening_tolerance =
    local_refinement_tolerance * _coarsen_threshold;

  // Prepare another error vector if we need to sum parent errors
  ErrorVector error_per_parent;
  if (_coarsen_by_parents)
    {
      Real parent_error_min, parent_error_max;

      create_parent_error_vector(error_per_cell_in,
                                 error_per_parent,
                                 parent_error_min,
                                 parent_error_max);
    }

  MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem * elem = *elem_it;
      Elem * parent = elem->parent();
      const dof_id_type elem_number    = elem->id();
      const ErrorVectorReal elem_error = error_per_cell_in[elem_number];

      if (elem_error > local_refinement_tolerance &&
          elem->level() < _max_h_level)
        elem->set_refinement_flag(Elem::REFINE);

      if (!_coarsen_by_parents && elem_error <
          local_coarsening_tolerance)
        elem->set_refinement_flag(Elem::COARSEN);

      if (_coarsen_by_parents && parent)
        {
          ErrorVectorReal parent_error = error_per_parent[parent->id()];
          if (parent_error >= 0.)
            {
              const Real parent_coarsening_tolerance =
                std::sqrt(parent->n_children() *
                          local_coarsening_tolerance *
                          local_coarsening_tolerance);
              if (parent_error < parent_coarsening_tolerance)
                elem->set_refinement_flag(Elem::COARSEN);
            }
        }
    }
}



bool MeshRefinement::flag_elements_by_nelem_target (const ErrorVector & error_per_cell)
{
  parallel_object_only();

  // Check for valid fractions..
  // The fraction values must be in [0,1]
  libmesh_assert_greater_equal (_refine_fraction, 0);
  libmesh_assert_less_equal (_refine_fraction, 1);
  libmesh_assert_greater_equal (_coarsen_fraction, 0);
  libmesh_assert_less_equal (_coarsen_fraction, 1);

  // This function is currently only coded to work when coarsening by
  // parents - it's too hard to guess how many coarsenings will be
  // performed otherwise.
  libmesh_assert (_coarsen_by_parents);

  // The number of active elements in the mesh - hopefully less than
  // 2 billion on 32 bit machines
  const dof_id_type n_active_elem  = _mesh.n_active_elem();

  // The maximum number of active elements to flag for coarsening
  const dof_id_type max_elem_coarsen =
    static_cast<dof_id_type>(_coarsen_fraction * n_active_elem) + 1;

  // The maximum number of elements to flag for refinement
  const dof_id_type max_elem_refine  =
    static_cast<dof_id_type>(_refine_fraction  * n_active_elem) + 1;

  // Clean up the refinement flags.  These could be left
  // over from previous refinement steps.
  this->clean_refinement_flags();

  // The target number of elements to add or remove
  const std::ptrdiff_t n_elem_new =
    std::ptrdiff_t(_nelem_target) - std::ptrdiff_t(n_active_elem);

  // Create an vector with active element errors and ids,
  // sorted by highest errors first
  const dof_id_type max_elem_id = _mesh.max_elem_id();
  std::vector<std::pair<ErrorVectorReal, dof_id_type> > sorted_error;

  sorted_error.reserve (n_active_elem);

  // On a DistributedMesh, we need to communicate to know which remote ids
  // correspond to active elements.
  {
    std::vector<bool> is_active(max_elem_id, false);

    MeshBase::element_iterator       elem_it  = _mesh.active_local_elements_begin();
    const MeshBase::element_iterator elem_end = _mesh.active_local_elements_end();
    for (; elem_it != elem_end; ++elem_it)
      {
        const dof_id_type eid = (*elem_it)->id();
        is_active[eid] = true;
        libmesh_assert_less (eid, error_per_cell.size());
        sorted_error.push_back
          (std::make_pair(error_per_cell[eid], eid));
      }

    this->comm().max(is_active);

    this->comm().allgather(sorted_error);
  }

  // Default sort works since pairs are sorted lexicographically
  std::sort (sorted_error.begin(), sorted_error.end());
  std::reverse (sorted_error.begin(), sorted_error.end());

  // Create a sorted error vector with coarsenable parent elements
  // only, sorted by lowest errors first
  ErrorVector error_per_parent;
  std::vector<std::pair<ErrorVectorReal, dof_id_type> > sorted_parent_error;
  Real parent_error_min, parent_error_max;

  create_parent_error_vector(error_per_cell,
                             error_per_parent,
                             parent_error_min,
                             parent_error_max);

  // create_parent_error_vector sets values for non-parents and
  // non-coarsenable parents to -1.  Get rid of them.
  for (std::size_t i=0; i != error_per_parent.size(); ++i)
    if (error_per_parent[i] != -1)
      sorted_parent_error.push_back(std::make_pair(error_per_parent[i], i));

  std::sort (sorted_parent_error.begin(), sorted_parent_error.end());

  // Keep track of how many elements we plan to coarsen & refine
  dof_id_type coarsen_count = 0;
  dof_id_type refine_count = 0;

  const unsigned int dim = _mesh.mesh_dimension();
  unsigned int twotodim = 1;
  for (unsigned int i=0; i!=dim; ++i)
    twotodim *= 2;

  // First, let's try to get our element count to target_nelem
  if (n_elem_new >= 0)
    {
      // Every element refinement creates at least
      // 2^dim-1 new elements
      refine_count =
        std::min(cast_int<dof_id_type>(n_elem_new / (twotodim-1)),
                 max_elem_refine);
    }
  else
    {
      // Every successful element coarsening is likely to destroy
      // 2^dim-1 net elements.
      coarsen_count =
        std::min(cast_int<dof_id_type>(-n_elem_new / (twotodim-1)),
                 max_elem_coarsen);
    }

  // Next, let's see if we can trade any refinement for coarsening
  while (coarsen_count < max_elem_coarsen &&
         refine_count < max_elem_refine &&
         coarsen_count < sorted_parent_error.size() &&
         refine_count < sorted_error.size() &&
         sorted_error[refine_count].first >
         sorted_parent_error[coarsen_count].first * _coarsen_threshold)
    {
      coarsen_count++;
      refine_count++;
    }

  // On a DistributedMesh, we need to communicate to know which remote ids
  // correspond to refinable elements
  dof_id_type successful_refine_count = 0;
  {
    std::vector<bool> is_refinable(max_elem_id, false);

    for (std::size_t i=0; i != sorted_error.size(); ++i)
      {
        dof_id_type eid = sorted_error[i].second;
        Elem * elem = _mesh.query_elem_ptr(eid);
        if (elem && elem->level() < _max_h_level)
          is_refinable[eid] = true;
      }
    this->comm().max(is_refinable);

    if (refine_count > max_elem_refine)
      refine_count = max_elem_refine;
    for (std::size_t i=0; i != sorted_error.size(); ++i)
      {
        if (successful_refine_count >= refine_count)
          break;

        dof_id_type eid = sorted_error[i].second;
        Elem * elem = _mesh.query_elem_ptr(eid);
        if (is_refinable[eid])
          {
            if (elem)
              elem->set_refinement_flag(Elem::REFINE);
            successful_refine_count++;
          }
      }
  }

  // If we couldn't refine enough elements, don't coarsen too many
  // either
  if (coarsen_count < (refine_count - successful_refine_count))
    coarsen_count = 0;
  else
    coarsen_count -= (refine_count - successful_refine_count);

  if (coarsen_count > max_elem_coarsen)
    coarsen_count = max_elem_coarsen;

  dof_id_type successful_coarsen_count = 0;
  if (coarsen_count)
    {
      for (std::size_t i=0; i != sorted_parent_error.size(); ++i)
        {
          if (successful_coarsen_count >= coarsen_count * twotodim)
            break;

          dof_id_type parent_id = sorted_parent_error[i].second;
          Elem * parent = _mesh.query_elem_ptr(parent_id);

          // On a DistributedMesh we skip remote elements
          if (!parent)
            continue;

          libmesh_assert(parent->has_children());
          for (unsigned int c=0; c != parent->n_children(); ++c)
            {
              Elem * elem = parent->child_ptr(c);
              if (elem && elem != remote_elem)
                {
                  libmesh_assert(elem->active());
                  elem->set_refinement_flag(Elem::COARSEN);
                  successful_coarsen_count++;
                }
            }
        }
    }

  // Return true if we've done all the AMR/C we can
  if (!successful_coarsen_count &&
      !successful_refine_count)
    return true;
  // And false if there may still be more to do.
  return false;
}



void MeshRefinement::flag_elements_by_elem_fraction (const ErrorVector & error_per_cell,
                                                     const Real refine_frac,
                                                     const Real coarsen_frac,
                                                     const unsigned int max_l)
{
  parallel_object_only();

  // The function arguments are currently just there for
  // backwards_compatibility
  if (!_use_member_parameters)
    {
      // If the user used non-default parameters, lets warn
      // that they're deprecated
      if (refine_frac != 0.3 ||
          coarsen_frac != 0.0 ||
          max_l != libMesh::invalid_uint)
        libmesh_deprecated();

      _refine_fraction = refine_frac;
      _coarsen_fraction = coarsen_frac;
      _max_h_level = max_l;
    }

  // Check for valid fractions..
  // The fraction values must be in [0,1]
  libmesh_assert_greater_equal (_refine_fraction, 0);
  libmesh_assert_less_equal (_refine_fraction, 1);
  libmesh_assert_greater_equal (_coarsen_fraction, 0);
  libmesh_assert_less_equal (_coarsen_fraction, 1);

  // The number of active elements in the mesh
  const dof_id_type n_active_elem  = _mesh.n_active_elem();

  // The number of elements to flag for coarsening
  const dof_id_type n_elem_coarsen =
    static_cast<dof_id_type>(_coarsen_fraction * n_active_elem);

  // The number of elements to flag for refinement
  const dof_id_type n_elem_refine =
    static_cast<dof_id_type>(_refine_fraction  * n_active_elem);



  // Clean up the refinement flags.  These could be left
  // over from previous refinement steps.
  this->clean_refinement_flags();


  // This vector stores the error and element number for all the
  // active elements.  It will be sorted and the top & bottom
  // elements will then be flagged for coarsening & refinement
  std::vector<ErrorVectorReal> sorted_error;

  sorted_error.reserve (n_active_elem);

  // Loop over the active elements and create the entry
  // in the sorted_error vector
  MeshBase::element_iterator       elem_it  = _mesh.active_local_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_local_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    sorted_error.push_back (error_per_cell[(*elem_it)->id()]);

  this->comm().allgather(sorted_error);

  // Now sort the sorted_error vector
  std::sort (sorted_error.begin(), sorted_error.end());

  // If we're coarsening by parents:
  // Create a sorted error vector with coarsenable parent elements
  // only, sorted by lowest errors first
  ErrorVector error_per_parent, sorted_parent_error;
  if (_coarsen_by_parents)
    {
      Real parent_error_min, parent_error_max;

      create_parent_error_vector(error_per_cell,
                                 error_per_parent,
                                 parent_error_min,
                                 parent_error_max);

      sorted_parent_error = error_per_parent;
      std::sort (sorted_parent_error.begin(), sorted_parent_error.end());

      // All the other error values will be 0., so get rid of them.
      sorted_parent_error.erase (std::remove(sorted_parent_error.begin(),
                                             sorted_parent_error.end(), 0.),
                                 sorted_parent_error.end());
    }


  ErrorVectorReal top_error= 0., bottom_error = 0.;

  // Get the maximum error value corresponding to the
  // bottom n_elem_coarsen elements
  if (_coarsen_by_parents && n_elem_coarsen)
    {
      const unsigned int dim = _mesh.mesh_dimension();
      unsigned int twotodim = 1;
      for (unsigned int i=0; i!=dim; ++i)
        twotodim *= 2;

      dof_id_type n_parent_coarsen = n_elem_coarsen / (twotodim - 1);

      if (n_parent_coarsen)
        bottom_error = sorted_parent_error[n_parent_coarsen - 1];
    }
  else if (n_elem_coarsen)
    {
      bottom_error = sorted_error[n_elem_coarsen - 1];
    }

  if (n_elem_refine)
    top_error = sorted_error[sorted_error.size() - n_elem_refine];

  // Finally, let's do the element flagging
  elem_it  = _mesh.active_elements_begin();
  for (; elem_it != elem_end; ++elem_it)
    {
      Elem * elem = *elem_it;
      Elem * parent = elem->parent();

      if (_coarsen_by_parents && parent && n_elem_coarsen &&
          error_per_parent[parent->id()] <= bottom_error)
        elem->set_refinement_flag(Elem::COARSEN);

      if (!_coarsen_by_parents && n_elem_coarsen &&
          error_per_cell[elem->id()] <= bottom_error)
        elem->set_refinement_flag(Elem::COARSEN);

      if (n_elem_refine &&
          elem->level() < _max_h_level &&
          error_per_cell[elem->id()] >= top_error)
        elem->set_refinement_flag(Elem::REFINE);
    }
}



void MeshRefinement::flag_elements_by_mean_stddev (const ErrorVector & error_per_cell,
                                                   const Real refine_frac,
                                                   const Real coarsen_frac,
                                                   const unsigned int max_l)
{
  // The function arguments are currently just there for
  // backwards_compatibility
  if (!_use_member_parameters)
    {
      // If the user used non-default parameters, lets warn
      // that they're deprecated
      if (refine_frac != 0.3 ||
          coarsen_frac != 0.0 ||
          max_l != libMesh::invalid_uint)
        libmesh_deprecated();

      _refine_fraction = refine_frac;
      _coarsen_fraction = coarsen_frac;
      _max_h_level = max_l;
    }

  // Get the mean value from the error vector
  const Real mean = error_per_cell.mean();

  // Get the standard deviation.  This equals the
  // square-root of the variance
  const Real stddev = std::sqrt (error_per_cell.variance());

  // Check for valid fractions
  libmesh_assert_greater_equal (_refine_fraction, 0);
  libmesh_assert_less_equal (_refine_fraction, 1);
  libmesh_assert_greater_equal (_coarsen_fraction, 0);
  libmesh_assert_less_equal (_coarsen_fraction, 1);

  // The refine and coarsen cutoff
  const Real refine_cutoff  =  mean + _refine_fraction  * stddev;
  const Real coarsen_cutoff =  std::max(mean - _coarsen_fraction * stddev, 0.);

  // Loop over the elements and flag them for coarsening or
  // refinement based on the element error
  MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem * elem             = *elem_it;
      const dof_id_type id  = elem->id();

      libmesh_assert_less (id, error_per_cell.size());

      const ErrorVectorReal elem_error = error_per_cell[id];

      // Possibly flag the element for coarsening ...
      if (elem_error <= coarsen_cutoff)
        elem->set_refinement_flag(Elem::COARSEN);

      // ... or refinement
      if ((elem_error >= refine_cutoff) && (elem->level() < _max_h_level))
        elem->set_refinement_flag(Elem::REFINE);
    }
}



void MeshRefinement::flag_elements_by (ElementFlagging & element_flagging)
{
  element_flagging.flag_elements();
}



void MeshRefinement::switch_h_to_p_refinement ()
{
  MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      if ((*elem_it)->active())
        {
          (*elem_it)->set_p_refinement_flag((*elem_it)->refinement_flag());
          (*elem_it)->set_refinement_flag(Elem::DO_NOTHING);
        }
      else
        {
          (*elem_it)->set_p_refinement_flag((*elem_it)->refinement_flag());
          (*elem_it)->set_refinement_flag(Elem::INACTIVE);
        }
    }
}



void MeshRefinement::add_p_to_h_refinement ()
{
  MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    (*elem_it)->set_p_refinement_flag((*elem_it)->refinement_flag());
}



void MeshRefinement::clean_refinement_flags ()
{
  // Possibly clean up the refinement flags from
  // a previous step
  //   elem_iterator       elem_it (_mesh.elements_begin());
  //   const elem_iterator elem_end(_mesh.elements_end());

  MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      if ((*elem_it)->active())
        {
          (*elem_it)->set_refinement_flag(Elem::DO_NOTHING);
          (*elem_it)->set_p_refinement_flag(Elem::DO_NOTHING);
        }
      else
        {
          (*elem_it)->set_refinement_flag(Elem::INACTIVE);
          (*elem_it)->set_p_refinement_flag(Elem::INACTIVE);
        }
    }
}

} // namespace libMesh

#endif
