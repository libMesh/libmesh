// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public  License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ Includes -------------------------------------
#include <set>
#include <algorithm> // for std::fill, std::equal_range, std::max, std::lower_bound, etc.
#include <sstream>

// Local Includes -----------------------------------
#include "libmesh/coupling_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector_base.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/fe_base.h" // FEBase::build() for continuity test
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/sparsity_pattern.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/threads.h"



namespace libMesh
{

// ------------------------------------------------------------
// DofMap member functions
AutoPtr<SparsityPattern::Build> DofMap::build_sparsity
  (const MeshBase& mesh) const
{
  libmesh_assert (mesh.is_prepared());
  libmesh_assert (this->n_variables());

  START_LOG("build_sparsity()", "DofMap");

  // Compute the sparsity structure of the global matrix.  This can be
  // fed into a PetscMatrix to allocate exacly the number of nonzeros
  // necessary to store the matrix.  This algorithm should be linear
  // in the (# of elements)*(# nodes per element)

  // We can be more efficient in the threaded sparsity pattern assembly
  // if we don't need the exact pattern.  For some sparse matrix formats
  // a good upper bound will suffice.

  // See if we need to include sparsity pattern entries for coupling
  // between neighbor dofs
  bool implicit_neighbor_dofs = this->use_coupled_neighbor_dofs(mesh);

  // We can compute the sparsity pattern in parallel on multiple
  // threads.  The goal is for each thread to compute the full sparsity
  // pattern for a subset of elements.  These sparsity patterns can
  // be efficiently merged in the SparsityPattern::Build::join()
  // method, especially if there is not too much overlap between them.
  // Even better, if the full sparsity pattern is not needed then
  // the number of nonzeros per row can be estimated from the
  // sparsity patterns created on each thread.
  AutoPtr<SparsityPattern::Build> sp
    (new SparsityPattern::Build (mesh,
			         *this,
			         this->_dof_coupling,
			         implicit_neighbor_dofs,
			         need_full_sparsity_pattern));

  Threads::parallel_reduce (ConstElemRange (mesh.active_local_elements_begin(),
					    mesh.active_local_elements_end()), *sp);

  sp->parallel_sync();

#ifndef NDEBUG
  // Avoid declaring these variables unless asserts are enabled.
  const processor_id_type proc_id        = mesh.processor_id();
  const dof_id_type n_dofs_on_proc = this->n_dofs_on_processor(proc_id);
#endif
  libmesh_assert_equal_to (sp->sparsity_pattern.size(), n_dofs_on_proc);

  STOP_LOG("build_sparsity()", "DofMap");

  // Check to see if we have any extra stuff to add to the sparsity_pattern
  if (_extra_sparsity_function)
    {
      if (_augment_sparsity_pattern)
	{
	  libmesh_here();
	  libMesh::out << "WARNING:  You have specified both an extra sparsity function and object.\n"
		       << "          Are you sure this is what you meant to do??"
		       << std::endl;
	}

      _extra_sparsity_function
        (sp->sparsity_pattern, sp->n_nz,
         sp->n_oz, _extra_sparsity_context);
    }

  if (_augment_sparsity_pattern)
    _augment_sparsity_pattern->augment_sparsity_pattern
      (sp->sparsity_pattern, sp->n_nz, sp->n_oz);

  return sp;
}



DofMap::DofMap(const unsigned int number,
	       const ParallelObject &parent_decomp) :
  ParallelObject (parent_decomp),
  _dof_coupling(NULL),
  _variables(),
  _variable_groups(),
  _sys_number(number),
  _matrices(),
  _first_df(),
  _end_df(),
  _send_list(),
  _augment_sparsity_pattern(NULL),
  _extra_sparsity_function(NULL),
  _extra_sparsity_context(NULL),
  _augment_send_list(NULL),
  _extra_send_list_function(NULL),
  _extra_send_list_context(NULL),
  need_full_sparsity_pattern(false),
  _n_nz(NULL),
  _n_oz(NULL),
  _n_dfs(0),
  _n_SCALAR_dofs(0)
#ifdef LIBMESH_ENABLE_AMR
  , _n_old_dfs(0),
  _first_old_df(),
  _end_old_df()
#endif
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  , _dof_constraints()
  , _primal_constraint_values()
  , _adjoint_constraint_values()
#endif
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  , _node_constraints()
#endif
#ifdef LIBMESH_ENABLE_PERIODIC
  , _periodic_boundaries(new PeriodicBoundaries)
#endif
#ifdef LIBMESH_ENABLE_DIRICHLET
  , _dirichlet_boundaries(new DirichletBoundaries)
  , _adjoint_dirichlet_boundaries()
#endif
{
  _matrices.clear();
}



// Destructor
DofMap::~DofMap()
{
  this->clear();
#ifdef LIBMESH_ENABLE_PERIODIC
  delete _periodic_boundaries;
#endif
#ifdef LIBMESH_ENABLE_DIRICHLET
  delete _dirichlet_boundaries;
  for (unsigned int q = 0; q != _adjoint_dirichlet_boundaries.size(); ++q)
    delete _adjoint_dirichlet_boundaries[q];
#endif
}


#ifdef LIBMESH_ENABLE_PERIODIC

bool DofMap::is_periodic_boundary (const boundary_id_type boundaryid) const
{
  if (_periodic_boundaries->count(boundaryid) != 0)
    return true;

  return false;
}

#endif



// void DofMap::add_variable (const Variable &var)
// {
//   libmesh_error();
//   _variables.push_back (var);
// }



void DofMap::add_variable_group (const VariableGroup &var_group)
{
  _variable_groups.push_back(var_group);

  VariableGroup &new_var_group = _variable_groups.back();

  for (unsigned int var=0; var<new_var_group.n_variables(); var++)
    _variables.push_back (new_var_group(var));
}



void DofMap::attach_matrix (SparseMatrix<Number>& matrix)
{
  parallel_object_only();

  // We shouldn't be trying to re-attach the same matrices repeatedly
  libmesh_assert (std::find(_matrices.begin(), _matrices.end(),
                            &matrix) == _matrices.end());

  _matrices.push_back(&matrix);

  matrix.attach_dof_map (*this);

  // If we've already computed sparsity, then it's too late
  // to wait for "compute_sparsity" to help with sparse matrix
  // initialization, and we need to handle this matrix individually
  bool computed_sparsity_already =
    ((_n_nz && !_n_nz->empty()) ||
     (_n_oz && !_n_oz->empty()));
  this->comm().max(computed_sparsity_already);
  if (computed_sparsity_already &&
      matrix.need_full_sparsity_pattern())
    {
      // We'd better have already computed the full sparsity pattern
      // if we need it here
      libmesh_assert(need_full_sparsity_pattern);
      libmesh_assert(_sp.get());

      matrix.update_sparsity_pattern (_sp->sparsity_pattern);
    }

  if (matrix.need_full_sparsity_pattern())
    need_full_sparsity_pattern = true;
}



bool DofMap::is_attached (SparseMatrix<Number>& matrix)
{
  return (std::find(_matrices.begin(), _matrices.end(),
                            &matrix) != _matrices.end());
}



DofObject* DofMap::node_ptr(MeshBase& mesh, dof_id_type i) const
{
  return mesh.node_ptr(i);
}



DofObject* DofMap::elem_ptr(MeshBase& mesh, dof_id_type i) const
{
  return mesh.elem(i);
}



template <typename iterator_type>
void DofMap::set_nonlocal_dof_objects(iterator_type objects_begin,
                                      iterator_type objects_end,
                                      MeshBase &mesh,
                                      dofobject_accessor objects)
{
  // This function must be run on all processors at once
  parallel_object_only();

  // First, iterate over local objects to find out how many
  // are on each processor
  std::vector<dof_id_type>
    ghost_objects_from_proc(this->n_processors(), 0);

  iterator_type it  = objects_begin;

  for (; it != objects_end; ++it)
    {
      DofObject *obj = *it;

      if (obj)
        {
          processor_id_type obj_procid = obj->processor_id();
          // We'd better be completely partitioned by now
          libmesh_assert_not_equal_to (obj_procid, DofObject::invalid_processor_id);
          ghost_objects_from_proc[obj_procid]++;
        }
    }

  std::vector<dof_id_type> objects_on_proc(this->n_processors(), 0);
  this->comm().allgather(ghost_objects_from_proc[this->processor_id()],
				 objects_on_proc);

#ifdef DEBUG
  for (processor_id_type p=0; p != this->n_processors(); ++p)
    libmesh_assert_less_equal (ghost_objects_from_proc[p], objects_on_proc[p]);
#endif

  // Request sets to send to each processor
  std::vector<std::vector<dof_id_type> >
    requested_ids(this->n_processors());

  // We know how many of our objects live on each processor, so
  // reserve() space for requests from each.
  for (processor_id_type p=0; p != this->n_processors(); ++p)
    if (p != this->processor_id())
      requested_ids[p].reserve(ghost_objects_from_proc[p]);

  for (it = objects_begin; it != objects_end; ++it)
    {
      DofObject *obj = *it;
      if (obj->processor_id() != DofObject::invalid_processor_id)
        requested_ids[obj->processor_id()].push_back(obj->id());
    }
#ifdef DEBUG
  for (processor_id_type p=0; p != this->n_processors(); ++p)
    libmesh_assert_equal_to (requested_ids[p].size(), ghost_objects_from_proc[p]);
#endif

  // Next set ghost object n_comps from other processors
  for (processor_id_type p=1; p != this->n_processors(); ++p)
    {
      // Trade my requests with processor procup and procdown
      processor_id_type procup = (this->processor_id() + p) %
                                  this->n_processors();
      processor_id_type procdown = (this->n_processors() +
                                    this->processor_id() - p) %
                                    this->n_processors();
      std::vector<dof_id_type> request_to_fill;
      this->comm().send_receive(procup, requested_ids[procup],
					procdown, request_to_fill);

      // Fill those requests
      const unsigned int
	sys_num      = this->sys_number(),
	n_var_groups = this->n_variable_groups();

      std::vector<dof_id_type> ghost_data
        (request_to_fill.size() * 2 * n_var_groups);

      for (std::size_t i=0; i != request_to_fill.size(); ++i)
        {
          DofObject *requested = (this->*objects)(mesh, request_to_fill[i]);
          libmesh_assert(requested);
          libmesh_assert_equal_to (requested->processor_id(), this->processor_id());
          libmesh_assert_equal_to (requested->n_var_groups(sys_num), n_var_groups);
          for (unsigned int vg=0; vg != n_var_groups; ++vg)
            {
              unsigned int n_comp_g =
                requested->n_comp_group(sys_num, vg);
              ghost_data[i*2*n_var_groups+vg] = n_comp_g;
              dof_id_type my_first_dof = n_comp_g ?
                requested->vg_dof_base(sys_num, vg) : 0;
              libmesh_assert_not_equal_to (my_first_dof, DofObject::invalid_id);
              ghost_data[i*2*n_var_groups+n_var_groups+vg] = my_first_dof;
            }
        }

      // Trade back the results
      std::vector<dof_id_type> filled_request;
      this->comm().send_receive(procdown, ghost_data,
					procup, filled_request);

      // And copy the id changes we've now been informed of
      libmesh_assert_equal_to (filled_request.size(),
                              requested_ids[procup].size() * 2 * n_var_groups);
      for (std::size_t i=0; i != requested_ids[procup].size(); ++i)
        {
          DofObject *requested = (this->*objects)(mesh, requested_ids[procup][i]);
          libmesh_assert(requested);
          libmesh_assert_equal_to (requested->processor_id(), procup);
          for (unsigned int vg=0; vg != n_var_groups; ++vg)
            {
              unsigned int n_comp_g =
	        libmesh_cast_int<unsigned int>(filled_request[i*2*n_var_groups+vg]);
              requested->set_n_comp_group(sys_num, vg, n_comp_g);
              if (n_comp_g)
                {
                  dof_id_type my_first_dof =
                    filled_request[i*2*n_var_groups+n_var_groups+vg];
                  libmesh_assert_not_equal_to (my_first_dof, DofObject::invalid_id);
                  requested->set_vg_dof_base
                    (sys_num, vg, my_first_dof);
                }
            }
        }
    }

#ifdef DEBUG
  // Double check for invalid dofs
  for (it = objects_begin; it != objects_end; ++it)
    {
      DofObject *obj = *it;
      libmesh_assert (obj);
      unsigned int num_variables = obj->n_vars(this->sys_number());
      for (unsigned int v=0; v != num_variables; ++v)
        {
          unsigned int n_comp =
            obj->n_comp(this->sys_number(), v);
          dof_id_type my_first_dof = n_comp ?
            obj->dof_number(this->sys_number(), v, 0) : 0;
          libmesh_assert_not_equal_to (my_first_dof, DofObject::invalid_id);
        }
    }
#endif
}



void DofMap::reinit(MeshBase& mesh)
{
  libmesh_assert (mesh.is_prepared());

  START_LOG("reinit()", "DofMap");

  const unsigned int
    sys_num      = this->sys_number(),
    n_var_groups = this->n_variable_groups();

  // The DofObjects need to know how many variable groups we have, and
  // how many variables there are in each group.
  std::vector<unsigned int> n_vars_per_group; /**/ n_vars_per_group.reserve (n_var_groups);

  for (unsigned int vg=0; vg<n_var_groups; vg++)
    n_vars_per_group.push_back (this->variable_group(vg).n_variables());

#ifdef LIBMESH_ENABLE_AMR

  //------------------------------------------------------------
  // Clear the old_dof_objects for all the nodes
  // and elements so that we can overwrite them
  {
    MeshBase::node_iterator       node_it  = mesh.nodes_begin();
    const MeshBase::node_iterator node_end = mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      {
	(*node_it)->clear_old_dof_object();
	libmesh_assert (!(*node_it)->old_dof_object);
      }

    MeshBase::element_iterator       elem_it  = mesh.elements_begin();
    const MeshBase::element_iterator elem_end = mesh.elements_end();

    for ( ; elem_it != elem_end; ++elem_it)
      {
	(*elem_it)->clear_old_dof_object();
	libmesh_assert (!(*elem_it)->old_dof_object);
      }
  }


  //------------------------------------------------------------
  // Set the old_dof_objects for the elements that
  // weren't just created, if these old dof objects
  // had variables
  {
    MeshBase::element_iterator       elem_it  = mesh.elements_begin();
    const MeshBase::element_iterator elem_end = mesh.elements_end();

    for ( ; elem_it != elem_end; ++elem_it)
      {
	Elem* elem = *elem_it;

        // Skip the elements that were just refined
	if (elem->refinement_flag() == Elem::JUST_REFINED) continue;

	for (unsigned int n=0; n<elem->n_nodes(); n++)
	  {
	    Node* node = elem->get_node(n);

	    if (node->old_dof_object == NULL)
	      if (node->has_dofs(sys_num))
		node->set_old_dof_object();
	  }

	libmesh_assert (!elem->old_dof_object);

	if (elem->has_dofs(sys_num))
	  elem->set_old_dof_object();
      }
  }

#endif // #ifdef LIBMESH_ENABLE_AMR


  //------------------------------------------------------------
  // Then set the number of variables for each \p DofObject
  // equal to n_variables() for this system.  This will
  // handle new \p DofObjects that may have just been created
  {
    // All the nodes
    MeshBase::node_iterator       node_it  = mesh.nodes_begin();
    const MeshBase::node_iterator node_end = mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_vars_per_group(sys_num, n_vars_per_group);

    // All the elements
    MeshBase::element_iterator       elem_it  = mesh.elements_begin();
    const MeshBase::element_iterator elem_end = mesh.elements_end();

    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_vars_per_group(sys_num, n_vars_per_group);
  }


  // Zero _n_SCALAR_dofs, it will be updated below.
  this->_n_SCALAR_dofs = 0;

  //------------------------------------------------------------
  // Next allocate space for the DOF indices
  for (unsigned int vg=0; vg<n_var_groups; vg++)
    {
      const VariableGroup &vg_description = this->variable_group(vg);

      const unsigned int n_var_in_group = vg_description.n_variables();
      const FEType& base_fe_type        = vg_description.type();

      // Don't need to loop over elements for a SCALAR variable
      // Just increment _n_SCALAR_dofs
      if(base_fe_type.family == SCALAR)
	{
	  this->_n_SCALAR_dofs += base_fe_type.order*n_var_in_group;
	  continue;
	}

      // This should be constant even on p-refined elements
      const bool extra_hanging_dofs =
	FEInterface::extra_hanging_dofs(base_fe_type);

      // For all the active elements
      MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
      const MeshBase::element_iterator elem_end = mesh.active_elements_end();

      // Count vertex degrees of freedom first
      for ( ; elem_it != elem_end; ++elem_it)
	{
	  Elem* elem  = *elem_it;
	  libmesh_assert(elem);

	  // Skip the numbering if this variable is
	  // not active on this element's subdomain
	  if (!vg_description.active_on_subdomain(elem->subdomain_id()))
	    continue;

	  const ElemType type = elem->type();
          const unsigned int dim = elem->dim();

          FEType fe_type = base_fe_type;

#ifdef LIBMESH_ENABLE_AMR
          // Make sure we haven't done more p refinement than we can
          // handle
          if (elem->p_level() + base_fe_type.order >
              FEInterface::max_order(base_fe_type, type))
            {
#  ifdef DEBUG
              if (FEInterface::max_order(base_fe_type,type) <
                  static_cast<unsigned int>(base_fe_type.order))
                {
                  libMesh::err
                    << "ERROR: Finite element "
                    << Utility::enum_to_string(base_fe_type.family)
                    << " on geometric element "
                    << Utility::enum_to_string(type) << std::endl
                    << "only supports FEInterface::max_order = "
                    << FEInterface::max_order(base_fe_type,type)
                    << ", not fe_type.order = " << base_fe_type.order
                    << std::endl;

                  libmesh_error();
                }

              libMesh::err
                << "WARNING: Finite element "
                << Utility::enum_to_string(base_fe_type.family)
                << " on geometric element "
                << Utility::enum_to_string(type) << std::endl
                << "could not be p refined past FEInterface::max_order = "
                << FEInterface::max_order(base_fe_type,type)
                << std::endl;
#  endif
              elem->set_p_level(FEInterface::max_order(base_fe_type,type)
                                - base_fe_type.order);
            }
#endif

	  fe_type.order = static_cast<Order>(fe_type.order +
                                             elem->p_level());

	  // Allocate the vertex DOFs
	  for (unsigned int n=0; n<elem->n_nodes(); n++)
	    {
	      Node* node = elem->get_node(n);

	      if (elem->is_vertex(n))
	        {
	          const unsigned int old_node_dofs =
	            node->n_comp_group(sys_num, vg);

		  const unsigned int vertex_dofs =
		    std::max(FEInterface::n_dofs_at_node(dim, fe_type,
                                                         type, n),
                             old_node_dofs);

		  // Some discontinuous FEs have no vertex dofs
		  if (vertex_dofs > old_node_dofs)
		    {
		      node->set_n_comp_group(sys_num, vg,
					     vertex_dofs);

		      // Abusing dof_number to set a "this is a
		      // vertex" flag
		      node->set_vg_dof_base(sys_num, vg,
					    vertex_dofs);

		      // libMesh::out << "sys_num,vg,old_node_dofs,vertex_dofs="
		      // 	      << sys_num << ","
		      // 	      << vg << ","
		      // 	      << old_node_dofs << ","
		      // 	      << vertex_dofs << '\n',
		      // 	node->debug_buffer();

		      // libmesh_assert_equal_to (vertex_dofs, node->n_comp(sys_num, vg));
		      // libmesh_assert_equal_to (vertex_dofs, node->vg_dof_base(sys_num, vg));
		    }
	        }
	    }
	} // done counting vertex dofs

      // count edge & face dofs next
      elem_it = mesh.active_elements_begin();

      for ( ; elem_it != elem_end; ++elem_it)
	{
	  Elem* elem = *elem_it;
	  libmesh_assert(elem);

	  // Skip the numbering if this variable is
	  // not active on this element's subdomain
	  if (!vg_description.active_on_subdomain(elem->subdomain_id()))
	    continue;

	  const ElemType type = elem->type();
          const unsigned int dim = elem->dim();

          FEType fe_type = base_fe_type;
          fe_type.order = static_cast<Order>(fe_type.order +
                                             elem->p_level());

	  // Allocate the edge and face DOFs
	  for (unsigned int n=0; n<elem->n_nodes(); n++)
	    {
	      Node* node = elem->get_node(n);

	      const unsigned int old_node_dofs =
	        node->n_comp_group(sys_num, vg);

              const unsigned int vertex_dofs = old_node_dofs?
                libmesh_cast_int<unsigned int>(node->vg_dof_base (sys_num,vg)):0;

	      const unsigned int new_node_dofs =
		FEInterface::n_dofs_at_node(dim, fe_type, type, n);

	      // We've already allocated vertex DOFs
	      if (elem->is_vertex(n))
	        {
		  libmesh_assert_greater_equal (old_node_dofs, vertex_dofs);
		  // //if (vertex_dofs < new_node_dofs)
		  //   libMesh::out << "sys_num,vg,old_node_dofs,vertex_dofs,new_node_dofs="
		  //                << sys_num << ","
		  //                << vg << ","
		  //                << old_node_dofs << ","
		  //                << vertex_dofs << ","
		  //                << new_node_dofs << '\n',
		  //     node->debug_buffer();

		  libmesh_assert_greater_equal (vertex_dofs,   new_node_dofs);
		}
	      // We need to allocate the rest
	      else
	        {
		  // If this has no dofs yet, it needs no vertex
		  // dofs, so we just give it edge or face dofs
		  if (!old_node_dofs)
                    {
		      node->set_n_comp_group(sys_num, vg,
					     new_node_dofs);
		      // Abusing dof_number to set a "this has no
		      // vertex dofs" flag
                      if (new_node_dofs)
		        node->set_vg_dof_base(sys_num, vg,
					      0);
                    }

		  // If this has dofs, but has no vertex dofs,
		  // it may still need more edge or face dofs if
                  // we're p-refined.
		  else if (vertex_dofs == 0)
		    {
                      if (new_node_dofs > old_node_dofs)
                        {
		          node->set_n_comp_group(sys_num, vg,
						 new_node_dofs);

		          node->set_vg_dof_base(sys_num, vg,
						vertex_dofs);
                        }
		    }
		  // If this is another element's vertex,
		  // add more (non-overlapping) edge/face dofs if
		  // necessary
		  else if (extra_hanging_dofs)
		    {
                      if (new_node_dofs > old_node_dofs - vertex_dofs)
                        {
		          node->set_n_comp_group(sys_num, vg,
						 vertex_dofs + new_node_dofs);

		          node->set_vg_dof_base(sys_num, vg,
						vertex_dofs);
                        }
		    }
		  // If this is another element's vertex, add any
		  // (overlapping) edge/face dofs if necessary
		  else
		    {
		      libmesh_assert_greater_equal (old_node_dofs, vertex_dofs);
                      if (new_node_dofs > old_node_dofs)
                        {
		          node->set_n_comp_group(sys_num, vg,
						 new_node_dofs);

		          node->set_vg_dof_base (sys_num, vg,
						 vertex_dofs);
                        }
		    }
		}
	    }
          // Allocate the element DOFs
	  const unsigned int dofs_per_elem =
			  FEInterface::n_dofs_per_elem(dim, fe_type,
						       type);

	  elem->set_n_comp_group(sys_num, vg, dofs_per_elem);

	}
    } // end loop over variable groups

  // Calling DofMap::reinit() by itself makes little sense,
  // so we won't bother with nonlocal DofObjects.
  // Those will be fixed by distribute_dofs

  //------------------------------------------------------------
  // Finally, clear all the current DOF indices
  // (distribute_dofs expects them cleared!)
  this->invalidate_dofs(mesh);

  STOP_LOG("reinit()", "DofMap");
}



void DofMap::invalidate_dofs(MeshBase& mesh) const
{
  const unsigned int sys_num = this->sys_number();

  // All the nodes
  MeshBase::node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::node_iterator node_end = mesh.nodes_end();

  for ( ; node_it != node_end; ++node_it)
    (*node_it)->invalidate_dofs(sys_num);

  // All the elements
  MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = mesh.active_elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    (*elem_it)->invalidate_dofs(sys_num);
}



void DofMap::clear()
{
  // we don't want to clear
  // the coupling matrix!
  // It should not change...
  //_dof_coupling->clear();

  _variables.clear();
  _variable_groups.clear();
  _first_df.clear();
  _end_df.clear();
  _send_list.clear();
  this->clear_sparsity();
  need_full_sparsity_pattern = false;

#ifdef LIBMESH_ENABLE_AMR

  _dof_constraints.clear();
  _primal_constraint_values.clear();
  _adjoint_constraint_values.clear();
  _n_old_dfs = 0;
  _first_old_df.clear();
  _end_old_df.clear();

#endif

  _matrices.clear();

  _n_dfs = 0;
}



void DofMap::distribute_dofs (MeshBase& mesh)
{
  // This function must be run on all processors at once
  parallel_object_only();

  // Log how long it takes to distribute the degrees of freedom
  START_LOG("distribute_dofs()", "DofMap");

  libmesh_assert (mesh.is_prepared());

  const processor_id_type proc_id = this->processor_id();
  const processor_id_type n_proc  = this->n_processors();

//  libmesh_assert_greater (this->n_variables(), 0);
  libmesh_assert_less (proc_id, n_proc);

  // re-init in case the mesh has changed
  this->reinit(mesh);

  // By default distribute variables in a
  // var-major fashion, but allow run-time
  // specification
  bool node_major_dofs = libMesh::on_command_line ("--node_major_dofs");

  // The DOF counter, will be incremented as we encounter
  // new degrees of freedom
  dof_id_type next_free_dof = 0;

  // Clear the send list before we rebuild it
  _send_list.clear();

  // Set temporary DOF indices on this processor
  if (node_major_dofs)
    this->distribute_local_dofs_node_major (next_free_dof, mesh);
  else
    this->distribute_local_dofs_var_major (next_free_dof, mesh);

  // Get DOF counts on all processors
  std::vector<dof_id_type> dofs_on_proc(n_proc, 0);
  this->comm().allgather(next_free_dof, dofs_on_proc);

  // Resize and fill the _first_df and _end_df arrays
#ifdef LIBMESH_ENABLE_AMR
  _first_old_df = _first_df;
  _end_old_df = _end_df;
#endif

  _first_df.resize(n_proc);
  _end_df.resize (n_proc);

  // Get DOF offsets
  _first_df[0] = 0;
  for (processor_id_type i=1; i < n_proc; ++i)
    _first_df[i] = _end_df[i-1] = _first_df[i-1] + dofs_on_proc[i-1];
  _end_df[n_proc-1] = _first_df[n_proc-1] + dofs_on_proc[n_proc-1];

  // Clear all the current DOF indices
  // (distribute_dofs expects them cleared!)
  this->invalidate_dofs(mesh);

  next_free_dof = _first_df[proc_id];

  // Set permanent DOF indices on this processor
  if (node_major_dofs)
    this->distribute_local_dofs_node_major (next_free_dof, mesh);
  else
    this->distribute_local_dofs_var_major (next_free_dof, mesh);

  libmesh_assert_equal_to (next_free_dof, _end_df[proc_id]);

  //------------------------------------------------------------
  // At this point, all n_comp and dof_number values on local
  // DofObjects should be correct, but a ParallelMesh might have
  // incorrect values on non-local DofObjects.  Let's request the
  // correct values from each other processor.

  if (this->n_processors() > 1)
    {
      this->set_nonlocal_dof_objects(mesh.nodes_begin(),
                                     mesh.nodes_end(),
                                     mesh, &DofMap::node_ptr);

      this->set_nonlocal_dof_objects(mesh.elements_begin(),
                                     mesh.elements_end(),
                                     mesh, &DofMap::elem_ptr);
    }

#ifdef DEBUG
  {
    const unsigned int
      sys_num = this->sys_number();

    // Processors should all agree on DoF ids
    MeshTools::libmesh_assert_valid_dof_ids(mesh);

    // DoF processor ids should match DofObject processor ids
    MeshBase::const_node_iterator       node_it  = mesh.nodes_begin();
    const MeshBase::const_node_iterator node_end = mesh.nodes_end();
    for ( ; node_it != node_end; ++node_it)
      {
	DofObject const * const dofobj = *node_it;
	const processor_id_type proc_id = dofobj->processor_id();

	for (unsigned int v=0; v != dofobj->n_vars(sys_num); ++v)
          for (unsigned int c=0; c != dofobj->n_comp(sys_num,v); ++c)
	    {
	      const dof_id_type dofid = dofobj->dof_number(sys_num,v,c);
	      libmesh_assert_greater_equal (dofid, this->first_dof(proc_id));
	      libmesh_assert_less (dofid, this->end_dof(proc_id));
	    }
      }

    MeshBase::const_element_iterator       elem_it  = mesh.elements_begin();
    const MeshBase::const_element_iterator elem_end = mesh.elements_end();
    for ( ; elem_it != elem_end; ++elem_it)
      {
	DofObject const * const dofobj = *elem_it;
	const processor_id_type proc_id = dofobj->processor_id();

	for (unsigned int v=0; v != dofobj->n_vars(sys_num); ++v)
          for (unsigned int c=0; c != dofobj->n_comp(sys_num,v); ++c)
	    {
	      const dof_id_type dofid = dofobj->dof_number(sys_num,v,c);
	      libmesh_assert_greater_equal (dofid, this->first_dof(proc_id));
	      libmesh_assert_less (dofid, this->end_dof(proc_id));
	    }
      }
  }
#endif

  // Set the total number of degrees of freedom
#ifdef LIBMESH_ENABLE_AMR
  _n_old_dfs = _n_dfs;
#endif
  _n_dfs = _end_df[n_proc-1];

  STOP_LOG("distribute_dofs()", "DofMap");

  // Note that in the add_neighbors_to_send_list nodes on processor
  // boundaries that are shared by multiple elements are added for
  // each element.
  this->add_neighbors_to_send_list(mesh);

  // Here we used to clean up that data structure; now System and
  // EquationSystems call that for us, after we've added constraint
  // dependencies to the send_list too.
  // this->sort_send_list ();
}


void DofMap::distribute_local_dofs_node_major(dof_id_type &next_free_dof,
                                              MeshBase& mesh)
{
  const unsigned int sys_num       = this->sys_number();
  const unsigned int n_var_groups  = this->n_variable_groups();

  //-------------------------------------------------------------------------
  // First count and assign temporary numbers to local dofs
  MeshBase::element_iterator       elem_it  = mesh.active_local_elements_begin();
  const MeshBase::element_iterator elem_end = mesh.active_local_elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      // Only number dofs connected to active
      // elements on this processor.
      Elem* elem                 = *elem_it;
      const unsigned int n_nodes = elem->n_nodes();

      // First number the nodal DOFS
      for (unsigned int n=0; n<n_nodes; n++)
        {
          Node* node = elem->get_node(n);

          for (unsigned vg=0; vg<n_var_groups; vg++)
	    {
	      const VariableGroup &vg_description(this->variable_group(vg));

	      if( (vg_description.type().family != SCALAR) &&
                  (vg_description.active_on_subdomain(elem->subdomain_id())) )
		{
		  // assign dof numbers (all at once) if this is
		  // our node and if they aren't already there
		  if ((node->n_comp_group(sys_num,vg) > 0) &&
		      (node->processor_id() == this->processor_id()) &&
		      (node->vg_dof_base(sys_num,vg) ==
		       DofObject::invalid_id))
		    {
		      node->set_vg_dof_base(sys_num,
					    vg,
					    next_free_dof);
		      next_free_dof += (vg_description.n_variables()*
					node->n_comp_group(sys_num,vg));
		      //node->debug_buffer();
		    }
		}
	    }
        }

      // Now number the element DOFS
      for (unsigned vg=0; vg<n_var_groups; vg++)
	{
	  const VariableGroup &vg_description(this->variable_group(vg));

	  if ( (vg_description.type().family != SCALAR) &&
	       (vg_description.active_on_subdomain(elem->subdomain_id())) )
	    if (elem->n_comp_group(sys_num,vg) > 0)
	      {
		libmesh_assert_equal_to (elem->vg_dof_base(sys_num,vg),
					 DofObject::invalid_id);

		elem->set_vg_dof_base(sys_num,
				      vg,
				      next_free_dof);

		next_free_dof += (vg_description.n_variables()*
				  elem->n_comp(sys_num,vg));
	      }
	}
    } // done looping over elements


  // we may have missed assigning DOFs to nodes that we own
  // but to which we have no connected elements matching our
  // variable restriction criterion.  this will happen, for example,
  // if variable V is restricted to subdomain S.  We may not own
  // any elements which live in S, but we may own nodes which are
  // *connected* to elements which do.  in this scenario these nodes
  // will presently have unnumbered DOFs. we need to take care of
  // them here since we own them and no other processor will touch them.
  {
    MeshBase::node_iterator       node_it  = mesh.local_nodes_begin();
    const MeshBase::node_iterator node_end = mesh.local_nodes_end();

    for (; node_it != node_end; ++node_it)
      {
	Node *node = *node_it;
	libmesh_assert(node);

	for (unsigned vg=0; vg<n_var_groups; vg++)
	  {
	    const VariableGroup &vg_description(this->variable_group(vg));

	    if (node->n_comp_group(sys_num,vg))
	      if (node->vg_dof_base(sys_num,vg) == DofObject::invalid_id)
		{
		  node->set_vg_dof_base (sys_num,
					 vg,
					 next_free_dof);

		  next_free_dof += (vg_description.n_variables()*
				    node->n_comp(sys_num,vg));
		}
	  }
      }
  }

  // Finally, count up the SCALAR dofs
  this->_n_SCALAR_dofs = 0;
  for (unsigned vg=0; vg<n_var_groups; vg++)
    {
      const VariableGroup &vg_description(this->variable_group(vg));

      if( vg_description.type().family == SCALAR )
        {
          this->_n_SCALAR_dofs += (vg_description.n_variables()*
				   vg_description.type().order);
          continue;
        }
    }

  // Only increment next_free_dof if we're on the processor
  // that holds this SCALAR variable
  if ( this->processor_id() == (this->n_processors()-1) )
    next_free_dof += _n_SCALAR_dofs;

#ifdef DEBUG
  {
    // libMesh::out << "next_free_dof=" << next_free_dof << std::endl
    // 	      << "_n_SCALAR_dofs=" << _n_SCALAR_dofs << std::endl;

    // Make sure we didn't miss any nodes
    MeshTools::libmesh_assert_valid_procids<Node>(mesh);

    MeshBase::node_iterator       node_it  = mesh.local_nodes_begin();
    const MeshBase::node_iterator node_end = mesh.local_nodes_end();
    for (; node_it != node_end; ++node_it)
      {
	Node *obj = *node_it;
	libmesh_assert(obj);
	unsigned int n_var_g = obj->n_var_groups(this->sys_number());
	for (unsigned int vg=0; vg != n_var_g; ++vg)
	  {
	    unsigned int n_comp_g =
	      obj->n_comp_group(this->sys_number(), vg);
	    dof_id_type my_first_dof = n_comp_g ?
	      obj->vg_dof_base(this->sys_number(), vg) : 0;
	    libmesh_assert_not_equal_to (my_first_dof, DofObject::invalid_id);
	  }
      }
  }
#endif // DEBUG
}



void DofMap::distribute_local_dofs_var_major(dof_id_type &next_free_dof,
                                             MeshBase& mesh)
{
  const unsigned int sys_num      = this->sys_number();
  const unsigned int n_var_groups = this->n_variable_groups();

  //-------------------------------------------------------------------------
  // First count and assign temporary numbers to local dofs
  for (unsigned vg=0; vg<n_var_groups; vg++)
    {
      const VariableGroup &vg_description(this->variable_group(vg));

      const unsigned int n_vars_in_group = vg_description.n_variables();

      // Skip the SCALAR dofs
      if (vg_description.type().family == SCALAR)
        continue;

      MeshBase::element_iterator       elem_it  = mesh.active_local_elements_begin();
      const MeshBase::element_iterator elem_end = mesh.active_local_elements_end();

      for ( ; elem_it != elem_end; ++elem_it)
        {
          // Only number dofs connected to active
          // elements on this processor.
          Elem* elem  = *elem_it;

	  // ... and only variables which are active on
	  // on this element's subdomain
	  if (!vg_description.active_on_subdomain(elem->subdomain_id()))
	    continue;

          const unsigned int n_nodes = elem->n_nodes();

          // First number the nodal DOFS
          for (unsigned int n=0; n<n_nodes; n++)
            {
              Node* node = elem->get_node(n);

              // assign dof numbers (all at once) if this is
              // our node and if they aren't already there
              if ((node->n_comp_group(sys_num,vg) > 0) &&
                  (node->processor_id() == this->processor_id()) &&
                  (node->vg_dof_base(sys_num,vg) ==
                   DofObject::invalid_id))
                {
                  node->set_vg_dof_base(sys_num,
					vg,
					next_free_dof);

                  next_free_dof += (n_vars_in_group*
				    node->n_comp_group(sys_num,vg));
                }
            }

          // Now number the element DOFS
          if (elem->n_comp_group(sys_num,vg) > 0)
            {
              libmesh_assert_equal_to (elem->vg_dof_base(sys_num,vg),
                                       DofObject::invalid_id);

              elem->set_vg_dof_base(sys_num,
				    vg,
				    next_free_dof);

              next_free_dof += (n_vars_in_group*
				elem->n_comp_group(sys_num,vg));
            }
        } // end loop on elements

      // we may have missed assigning DOFs to nodes that we own
      // but to which we have no connected elements matching our
      // variable restriction criterion.  this will happen, for example,
      // if variable V is restricted to subdomain S.  We may not own
      // any elements which live in S, but we may own nodes which are
      // *connected* to elements which do.  in this scenario these nodes
      // will presently have unnumbered DOFs. we need to take care of
      // them here since we own them and no other processor will touch them.
      {
	MeshBase::node_iterator       node_it  = mesh.local_nodes_begin();
	const MeshBase::node_iterator node_end = mesh.local_nodes_end();

	for (; node_it != node_end; ++node_it)
	  {
	    Node *node = *node_it;
	    libmesh_assert(node);

	    if (node->n_comp_group(sys_num,vg))
	      if (node->vg_dof_base(sys_num,vg) == DofObject::invalid_id)
		{
		  node->set_vg_dof_base (sys_num,
					 vg,
					 next_free_dof);

		  next_free_dof += (n_vars_in_group*
				    node->n_comp_group(sys_num,vg));
		}
	  }
      }
    } // end loop on variable groups

  // Finally, count up the SCALAR dofs
  this->_n_SCALAR_dofs = 0;
  for (unsigned vg=0; vg<n_var_groups; vg++)
    {
      const VariableGroup &vg_description(this->variable_group(vg));

      if( vg_description.type().family == SCALAR )
        {
          this->_n_SCALAR_dofs += (vg_description.n_variables()*
				   vg_description.type().order);
          continue;
        }
    }

  // Only increment next_free_dof if we're on the processor
  // that holds this SCALAR variable
  if ( this->processor_id() == (this->n_processors()-1) )
    next_free_dof += _n_SCALAR_dofs;

#ifdef DEBUG
  {
    // Make sure we didn't miss any nodes
    MeshTools::libmesh_assert_valid_procids<Node>(mesh);

    MeshBase::node_iterator       node_it  = mesh.local_nodes_begin();
    const MeshBase::node_iterator node_end = mesh.local_nodes_end();
    for (; node_it != node_end; ++node_it)
      {
	Node *obj = *node_it;
	libmesh_assert(obj);
	unsigned int n_var_g = obj->n_var_groups(this->sys_number());
	for (unsigned int vg=0; vg != n_var_g; ++vg)
	  {
	    unsigned int n_comp_g =
	      obj->n_comp_group(this->sys_number(), vg);
	    dof_id_type my_first_dof = n_comp_g ?
	      obj->vg_dof_base(this->sys_number(), vg) : 0;
	    libmesh_assert_not_equal_to (my_first_dof, DofObject::invalid_id);
	  }
      }
  }
#endif // DEBUG
}



void DofMap::add_neighbors_to_send_list(MeshBase& mesh)
{
  START_LOG("add_neighbors_to_send_list()", "DofMap");

  const unsigned int sys_num = this->sys_number();

  //-------------------------------------------------------------------------
  // We need to add the DOFs from elements that live on neighboring processors
  // that are neighbors of the elements on the local processor
  //-------------------------------------------------------------------------

  MeshBase::const_element_iterator       local_elem_it
    = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator local_elem_end
    = mesh.active_local_elements_end();

  std::vector<bool> node_on_processor(mesh.max_node_id(), false);
  std::vector<dof_id_type> di;
  std::vector<const Elem *> family;

  // Loop over the active local elements, adding all active elements
  // that neighbor an active local element to the send list.
  for ( ; local_elem_it != local_elem_end; ++local_elem_it)
    {
      const Elem* elem = *local_elem_it;

      for (unsigned int n=0; n!=elem->n_nodes(); n++)
        {
	  // Flag all the nodes of active local elements as seen, so
	  // we can add nodal neighbor dofs to the send_list later.
          node_on_processor[elem->node(n)] = true;

          // Add all remote dofs on these nodes to the send_list.
          // This is necessary in case those dofs are *not* also dofs
          // on neighbors; e.g. in the case of a HIERARCHIC's local
          // side which is only a vertex on the neighbor that owns it.
          const Node* node = elem->get_node(n);
          const unsigned n_vars = node->n_vars(sys_num);
          for (unsigned int v=0; v != n_vars; ++v)
            {
              const unsigned int n_comp = node->n_comp(sys_num, v);
              for (unsigned int c=0; c != n_comp; ++c)
                {
                  const dof_id_type dof_index = node->dof_number(sys_num, v, c);
		  if (dof_index < this->first_dof() || dof_index >= this->end_dof())
		    {
		      _send_list.push_back(dof_index);
		      // libmesh_here();
		      // libMesh::out << "sys_num,v,c,dof_index="
		      // 		<< sys_num << ", "
		      // 		<< v << ", "
		      // 		<< c << ", "
		      // 		<< dof_index << '\n';
		      // node->debug_buffer();
		    }
                }
            }
        }

      // Loop over the neighbors of those elements
      for (unsigned int s=0; s<elem->n_neighbors(); s++)
	if (elem->neighbor(s) != NULL)
	  {
	    family.clear();

            // Find all the active elements that neighbor elem
#ifdef LIBMESH_ENABLE_AMR
            if (!elem->neighbor(s)->active())
              elem->neighbor(s)->active_family_tree_by_neighbor(family, elem);
            else
#endif
              family.push_back(elem->neighbor(s));

            for (dof_id_type i=0; i!=family.size(); ++i)
	      // If the neighbor lives on a different processor
	      if (family[i]->processor_id() != this->processor_id())
		{
		  // Get the DOF indices for this neighboring element
		  this->dof_indices (family[i], di);

	          // Insert the remote DOF indices into the send list
                  for (std::size_t j=0; j != di.size(); ++j)
                    if (di[j] < this->first_dof() ||
                        di[j] >= this->end_dof())
	              _send_list.push_back(di[j]);
	        }
	  }
    }

  // Now loop over all non_local active elements and add any missing
  // nodal-only neighbors.  This will also get any dofs from nonlocal
  // nodes on local elements, because every nonlocal node exists on a
  // nonlocal nodal neighbor element.
  MeshBase::const_element_iterator       elem_it
    = mesh.active_elements_begin();
  const MeshBase::const_element_iterator elem_end
    = mesh.active_elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;

      // If this is one of our elements, we've already added it
      if (elem->processor_id() == this->processor_id())
        continue;

      // Do we need to add the element DOFs?
      bool add_elem_dofs = false;

      // Check all the nodes of the element to see if it
      // shares a node with us
      for (unsigned int n=0; n!=elem->n_nodes(); n++)
        if (node_on_processor[elem->node(n)])
	  add_elem_dofs = true;

      // Add the element degrees of freedom if it shares at
      // least one node.
      if (add_elem_dofs)
	{
	  // Get the DOF indices for this neighboring element
	  this->dof_indices (elem, di);

	  // Insert the remote DOF indices into the send list
            for (std::size_t j=0; j != di.size(); ++j)
              if (di[j] < this->first_dof() ||
                  di[j] >= this->end_dof())
	        _send_list.push_back(di[j]);
	}
    }

  STOP_LOG("add_neighbors_to_send_list()", "DofMap");
}



void DofMap::prepare_send_list ()
{
  START_LOG("prepare_send_list()", "DofMap");

  // Check to see if we have any extra stuff to add to the send_list
  if (_extra_send_list_function)
    {
      if (_augment_send_list)
	{
	  libmesh_here();
	  libMesh::out << "WARNING:  You have specified both an extra send list function and object.\n"
		       << "          Are you sure this is what you meant to do??"
		       << std::endl;
	}

      _extra_send_list_function(_send_list, _extra_send_list_context);
    }

  if (_augment_send_list)
    _augment_send_list->augment_send_list (_send_list);

  // First sort the send list.  After this
  // duplicated elements will be adjacent in the
  // vector
  std::sort(_send_list.begin(), _send_list.end());

  // Now use std::unique to remove duplicate entries
  std::vector<dof_id_type>::iterator new_end =
    std::unique (_send_list.begin(), _send_list.end());

  // Remove the end of the send_list.  Use the "swap trick"
  // from Effective STL
  std::vector<dof_id_type> (_send_list.begin(), new_end).swap (_send_list);

  STOP_LOG("prepare_send_list()", "DofMap");
}



bool DofMap::use_coupled_neighbor_dofs(const MeshBase& mesh) const
{
  // If we were asked on the command line, then we need to
  // include sensitivities between neighbor degrees of freedom
  bool implicit_neighbor_dofs =
    libMesh::on_command_line ("--implicit_neighbor_dofs");

  // look at all the variables in this system.  If every one is
  // discontinuous then the user must be doing DG/FVM, so be nice
  // and  force implicit_neighbor_dofs=true
  {
    bool all_discontinuous_dofs = true;

    for (unsigned int var=0; var<this->n_variables(); var++)
      if (FEAbstract::build (mesh.mesh_dimension(),
			     this->variable_type(var))->get_continuity() !=  DISCONTINUOUS)
	all_discontinuous_dofs = false;

    if (all_discontinuous_dofs)
      implicit_neighbor_dofs = true;
  }

  return implicit_neighbor_dofs;
}



void DofMap::compute_sparsity(const MeshBase& mesh)
{
  _sp = this->build_sparsity(mesh);

  // It is possible that some \p SparseMatrix implementations want to
  // see it.  Let them see it before we throw it away.
  std::vector<SparseMatrix<Number>* >::const_iterator
    pos = _matrices.begin(),
    end = _matrices.end();

  // If we need the full sparsity pattern, then we share a view of its
  // arrays, and we pass it in to the matrices.
  if (need_full_sparsity_pattern)
    {
      _n_nz = &_sp->n_nz;
      _n_oz = &_sp->n_oz;

      for (; pos != end; ++pos)
        (*pos)->update_sparsity_pattern (_sp->sparsity_pattern);
    }
  // If we don't need the full sparsity pattern anymore, steal the
  // arrays we do need and free the rest of the memory
  else
    {
      if (!_n_nz)
        _n_nz = new std::vector<dof_id_type>();
      _n_nz->swap(_sp->n_nz);
      if (!_n_oz)
        _n_oz = new std::vector<dof_id_type>();
      _n_oz->swap(_sp->n_oz);

      _sp.reset();
    }
}



void DofMap::clear_sparsity()
{
  if (need_full_sparsity_pattern)
    {
      libmesh_assert(_sp.get());
      libmesh_assert(!_n_nz || _n_nz == &_sp->n_nz);
      libmesh_assert(!_n_oz || _n_oz == &_sp->n_oz);
      _sp.reset();
    }
  else
    {
      libmesh_assert(!_sp.get());
      delete _n_nz;
      delete _n_oz;
    }
  _n_nz = NULL;
  _n_oz = NULL;
}



void DofMap::extract_local_vector (const NumericVector<Number>& Ug,
				   const std::vector<dof_id_type>& dof_indices_in,
				   DenseVectorBase<Number>& Ue) const
{
#ifdef LIBMESH_ENABLE_AMR

  // Trivial mapping
  libmesh_assert_equal_to (dof_indices_in.size(), Ue.size());
  bool has_constrained_dofs = false;

  for (unsigned int il=0;
       il != libmesh_cast_int<unsigned int>(dof_indices_in.size()); il++)
    {
      const dof_id_type ig = dof_indices_in[il];

      if (this->is_constrained_dof (ig)) has_constrained_dofs = true;

      libmesh_assert_less (ig, Ug.size());

      Ue.el(il) = Ug(ig);
    }

  // If the element has any constrained DOFs then we need
  // to account for them in the mapping.  This will handle
  // the case that the input vector is not constrained.
  if (has_constrained_dofs)
    {
      // Copy the input DOF indices.
      std::vector<dof_id_type> constrained_dof_indices(dof_indices_in);

      DenseMatrix<Number> C;
      DenseVector<Number> H;

      this->build_constraint_matrix_and_vector (C, H, constrained_dof_indices);

      libmesh_assert_equal_to (dof_indices_in.size(), C.m());
      libmesh_assert_equal_to (constrained_dof_indices.size(), C.n());

      // zero-out Ue
      Ue.zero();

      // compute Ue = C Ug, with proper mapping.
      const unsigned int n_original_dofs =
        libmesh_cast_int<unsigned int>(dof_indices_in.size());
      for (unsigned int i=0; i != n_original_dofs; i++)
	{
	  Ue.el(i) = H(i);

          const unsigned int n_constrained =
            libmesh_cast_int<unsigned int>(constrained_dof_indices.size());
	  for (unsigned int j=0; j<n_constrained; j++)
	    {
	      const dof_id_type jg = constrained_dof_indices[j];

//          If Ug is a serial or ghosted vector, then this assert is
//          overzealous.  If Ug is a parallel vector, then this assert
//          is redundant.
//	    libmesh_assert ((jg >= Ug.first_local_index()) &&
//		    (jg <  Ug.last_local_index()));

	      Ue.el(i) += C(i,j)*Ug(jg);
	    }
	}
    }

#else

  // Trivial mapping

  const unsigned int n_original_dofs =
    libmesh_cast_int<unsigned int>(dof_indices_in.size());

  libmesh_assert_equal_to (n_original_dofs, Ue.size());

  for (unsigned int il=0; il<n_original_dofs; il++)
    {
      const dof_id_type ig = dof_indices_in[il];

      libmesh_assert ((ig >= Ug.first_local_index()) && (ig <  Ug.last_local_index()));

      Ue.el(il) = Ug(ig);
    }

#endif
}

void DofMap::dof_indices (const Elem* const elem,
			  std::vector<dof_id_type>& di,
			  const unsigned int vn) const
{
  START_LOG("dof_indices()", "DofMap");

  libmesh_assert(elem);

  const unsigned int n_nodes = elem->n_nodes();
  const ElemType type        = elem->type();
  const unsigned int sys_num = this->sys_number();
  const unsigned int n_vars  = this->n_variables();
  const unsigned int dim     = elem->dim();

  // Clear the DOF indices vector
  di.clear();

#ifdef DEBUG
  // Check that sizes match in DEBUG mode
  unsigned int tot_size = 0;
#endif

  // Create a vector to indicate which
  // SCALAR variables have been requested
  std::vector<unsigned int> SCALAR_var_numbers;
  SCALAR_var_numbers.clear();

  // Get the dof numbers
  for (unsigned int v=0; v<n_vars; v++)
    if ((v == vn) || (vn == libMesh::invalid_uint))
    {
      if(this->variable(v).type().family == SCALAR)
      {
        // We asked for this variable, so add it to the vector.
        SCALAR_var_numbers.push_back(v);

#ifdef DEBUG
        FEType fe_type = this->variable_type(v);
        tot_size += FEInterface::n_dofs(dim,
                                        fe_type,
					type);
#endif
      }
      else
      if (this->variable(v).active_on_subdomain(elem->subdomain_id()))
	{ // Do this for all the variables if one was not specified
	  // or just for the specified variable

	  // Increase the polynomial order on p refined elements
	  FEType fe_type = this->variable_type(v);
	  fe_type.order = static_cast<Order>(fe_type.order +
					     elem->p_level());

	  const bool extra_hanging_dofs =
	    FEInterface::extra_hanging_dofs(fe_type);

#ifdef DEBUG
	  tot_size += FEInterface::n_dofs(dim,
					  fe_type,
					  type);
#endif

	  // Get the node-based DOF numbers
	  for (unsigned int n=0; n<n_nodes; n++)
	    {
	      const Node* node      = elem->get_node(n);

	      // There is a potential problem with h refinement.  Imagine a
	      // quad9 that has a linear FE on it.  Then, on the hanging side,
	      // it can falsely identify a DOF at the mid-edge node. This is why
	      // we call FEInterface instead of node->n_comp() directly.
	      const unsigned int nc = FEInterface::n_dofs_at_node (dim,
								   fe_type,
								   type,
								   n);

	      // If this is a non-vertex on a hanging node with extra
	      // degrees of freedom, we use the non-vertex dofs (which
	      // come in reverse order starting from the end, to
	      // simplify p refinement)
	      if (extra_hanging_dofs && !elem->is_vertex(n))
		{
		  const int dof_offset = node->n_comp(sys_num,v) - nc;

		  // We should never have fewer dofs than necessary on a
		  // node unless we're getting indices on a parent element,
		  // and we should never need the indices on such a node
		  if (dof_offset < 0)
		    {
		      libmesh_assert(!elem->active());
		      di.resize(di.size() + nc, DofObject::invalid_id);
		    }
		  else
		    for (int i=node->n_comp(sys_num,v)-1; i>=dof_offset; i--)
		      {
			libmesh_assert_not_equal_to (node->dof_number(sys_num,v,i),
					             DofObject::invalid_id);
			di.push_back(node->dof_number(sys_num,v,i));
		      }
		}
	      // If this is a vertex or an element without extra hanging
	      // dofs, our dofs come in forward order coming from the
	      // beginning
	      else
		for (unsigned int i=0; i<nc; i++)
		  {
		    libmesh_assert_not_equal_to (node->dof_number(sys_num,v,i),
				                 DofObject::invalid_id);
		    di.push_back(node->dof_number(sys_num,v,i));
		  }
	    }

	  // If there are any element-based DOF numbers, get them
	  const unsigned int nc = FEInterface::n_dofs_per_elem(dim,
							       fe_type,
							       type);
	  // We should never have fewer dofs than necessary on an
	  // element unless we're getting indices on a parent element,
	  // and we should never need those indices
	  if (nc != 0)
	    {
	      if (elem->n_systems() > sys_num &&
		  nc <= elem->n_comp(sys_num,v))
		{
		  for (unsigned int i=0; i<nc; i++)
		    {
		      libmesh_assert_not_equal_to (elem->dof_number(sys_num,v,i),
				                   DofObject::invalid_id);

		      di.push_back(elem->dof_number(sys_num,v,i));
		    }
		}
	      else
		{
		  libmesh_assert(!elem->active() || fe_type.family == LAGRANGE);
		  di.resize(di.size() + nc, DofObject::invalid_id);
		}
	    }
	}
      } // end loop over variables

  // Finally append any SCALAR dofs that we asked for.
  std::vector<dof_id_type> di_new;
  std::vector<unsigned int>::iterator it           = SCALAR_var_numbers.begin();
  std::vector<unsigned int>::const_iterator it_end = SCALAR_var_numbers.end();
  for( ; it != it_end; ++it)
  {
    this->SCALAR_dof_indices(di_new,*it);
    di.insert( di.end(), di_new.begin(), di_new.end());
  }

#ifdef DEBUG
    libmesh_assert_equal_to (tot_size, di.size());
#endif

  STOP_LOG("dof_indices()", "DofMap");
}



void DofMap::SCALAR_dof_indices (std::vector<dof_id_type>& di,
			         const unsigned int vn,
#ifdef LIBMESH_ENABLE_AMR
			         const bool old_dofs
#else
			         const bool
#endif
				 ) const
{
  START_LOG("SCALAR_dof_indices()", "DofMap");

  if(this->variable(vn).type().family != SCALAR)
  {
    libMesh::err << "ERROR: SCALAR_dof_indices called for a non-SCALAR variable."
                  << std::endl;
  }

  // Clear the DOF indices vector
  di.clear();

  // First we need to find out the first dof
  // index for each SCALAR.
#ifdef LIBMESH_ENABLE_AMR
  dof_id_type first_SCALAR_dof_index = (old_dofs ? n_old_dofs() : n_dofs()) - n_SCALAR_dofs();
#else
  dof_id_type first_SCALAR_dof_index = n_dofs() - n_SCALAR_dofs();
#endif
  std::map<unsigned int, dof_id_type> SCALAR_first_dof_index;
  SCALAR_first_dof_index.clear();

  // Iterate over _all_ of the SCALARs and store each one's first dof index
  // We need to do this since the SCALAR dofs are packed contiguously
  for (unsigned int v=0; v<this->n_variables(); v++)
    if(this->variable(v).type().family == SCALAR)
    {
      unsigned int current_n_SCALAR_dofs = this->variable(v).type().order;
      SCALAR_first_dof_index.insert(
        std::pair<unsigned int, dof_id_type>(v,first_SCALAR_dof_index) );
      first_SCALAR_dof_index += current_n_SCALAR_dofs;
    }

  // Now use vn to index into SCALAR_first_dof_index
  std::map<unsigned int, dof_id_type>::const_iterator iter =
    SCALAR_first_dof_index.find(vn);

#ifdef DEBUG
  libmesh_assert (iter != SCALAR_first_dof_index.end());
#endif

  dof_id_type current_first_SCALAR_dof_index = iter->second;

  // Also, get the number of SCALAR dofs from the variable order
  unsigned int current_n_SCALAR_dofs = this->variable(vn).type().order;

  for(unsigned int j=0; j<current_n_SCALAR_dofs; j++)
  {
    dof_id_type index = current_first_SCALAR_dof_index+j;
    di.push_back(index);
  }

  STOP_LOG("SCALAR_dof_indices()", "DofMap");
}



bool DofMap::all_semilocal_indices (const std::vector<dof_id_type>& dof_indices_in) const
{
  // We're all semilocal unless we find a counterexample
  for (std::size_t i=0; i != dof_indices_in.size(); ++i)
    {
      const dof_id_type di = dof_indices_in[i];
      // If it's not in the local indices
      if (di < this->first_dof() ||
          di >= this->end_dof())
        {
          // and if it's not in the ghost indices, then we're not
          // semilocal
          if (!std::binary_search(_send_list.begin(), _send_list.end(), di))
            return false;
        }
    }

  return true;
}


#ifdef LIBMESH_ENABLE_AMR

void DofMap::old_dof_indices (const Elem* const elem,
			      std::vector<dof_id_type>& di,
			      const unsigned int vn) const
{
  START_LOG("old_dof_indices()", "DofMap");

  libmesh_assert(elem);
  libmesh_assert(elem->old_dof_object);


  const unsigned int n_nodes = elem->n_nodes();
  const ElemType type        = elem->type();
  const unsigned int sys_num = this->sys_number();
  const unsigned int n_vars  = this->n_variables();
  const unsigned int dim     = elem->dim();

  // Clear the DOF indices vector.
  di.clear();

#ifdef DEBUG
  // Check that sizes match
  unsigned int tot_size = 0;
#endif

  // Create a vector to indicate which
  // SCALAR variables have been requested
  std::vector<unsigned int> SCALAR_var_numbers;
  SCALAR_var_numbers.clear();

  // Get the dof numbers
  for (unsigned int v=0; v<n_vars; v++)
    if ((v == vn) || (vn == libMesh::invalid_uint))
    {
      if(this->variable(v).type().family == SCALAR)
        {
          // We asked for this variable, so add it to the vector.
          SCALAR_var_numbers.push_back(v);

#ifdef DEBUG
	  FEType fe_type = this->variable_type(v);
	  tot_size += FEInterface::n_dofs(dim,
				          fe_type,
					  type);
#endif
        }
      else
      if (this->variable(v).active_on_subdomain(elem->subdomain_id()))
	{ // Do this for all the variables if one was not specified
	  // or just for the specified variable

	  // Increase the polynomial order on p refined elements,
	  // but make sure you get the right polynomial order for
	  // the OLD degrees of freedom
	  int p_adjustment = 0;
	  if (elem->p_refinement_flag() == Elem::JUST_REFINED)
	    {
	      libmesh_assert_greater (elem->p_level(), 0);
	      p_adjustment = -1;
	    }
	  else if (elem->p_refinement_flag() == Elem::JUST_COARSENED)
	    {
	      p_adjustment = 1;
	    }
	  FEType fe_type = this->variable_type(v);
	  fe_type.order = static_cast<Order>(fe_type.order +
					     elem->p_level() +
					     p_adjustment);

	  const bool extra_hanging_dofs =
	    FEInterface::extra_hanging_dofs(fe_type);

#ifdef DEBUG
	  tot_size += FEInterface::n_dofs(dim,
					  fe_type,
					  type);
#endif

	  // Get the node-based DOF numbers
	  for (unsigned int n=0; n<n_nodes; n++)
	    {
	      const Node* node      = elem->get_node(n);

	      // There is a potential problem with h refinement.  Imagine a
	      // quad9 that has a linear FE on it.  Then, on the hanging side,
	      // it can falsely identify a DOF at the mid-edge node. This is why
	      // we call FEInterface instead of node->n_comp() directly.
	      const unsigned int nc = FEInterface::n_dofs_at_node (dim,
								   fe_type,
								   type,
								   n);
	      libmesh_assert(node->old_dof_object);

	      // If this is a non-vertex on a hanging node with extra
	      // degrees of freedom, we use the non-vertex dofs (which
	      // come in reverse order starting from the end, to
	      // simplify p refinement)
	      if (extra_hanging_dofs && !elem->is_vertex(n))
		{
		  const int dof_offset =
		    node->old_dof_object->n_comp(sys_num,v) - nc;

		  // We should never have fewer dofs than necessary on a
		  // node unless we're getting indices on a parent element
		  // or a just-coarsened element
		  if (dof_offset < 0)
		    {
		      libmesh_assert(!elem->active() || elem->refinement_flag() ==
				     Elem::JUST_COARSENED);
		      di.resize(di.size() + nc, DofObject::invalid_id);
		    }
		  else
		    for (int i=node->old_dof_object->n_comp(sys_num,v)-1;
			 i>=dof_offset; i--)
		      {
			libmesh_assert_not_equal_to (node->old_dof_object->dof_number(sys_num,v,i),
					             DofObject::invalid_id);
			di.push_back(node->old_dof_object->dof_number(sys_num,v,i));
		      }
		}
	      // If this is a vertex or an element without extra hanging
	      // dofs, our dofs come in forward order coming from the
	      // beginning
	      else
		for (unsigned int i=0; i<nc; i++)
		  {
		    libmesh_assert_not_equal_to (node->old_dof_object->dof_number(sys_num,v,i),
				                 DofObject::invalid_id);
		    di.push_back(node->old_dof_object->dof_number(sys_num,v,i));
		  }
	    }

	  // If there are any element-based DOF numbers, get them
	  const unsigned int nc = FEInterface::n_dofs_per_elem(dim,
							       fe_type,
							       type);

	  // We should never have fewer dofs than necessary on an
	  // element unless we're getting indices on a parent element
	  // or a just-coarsened element
	  if (nc != 0)
	    {
	      if (elem->old_dof_object->n_systems() > sys_num &&
		  nc <= elem->old_dof_object->n_comp(sys_num,v))
		{
		  libmesh_assert(elem->old_dof_object);

		  for (unsigned int i=0; i<nc; i++)
		    {
		      libmesh_assert_not_equal_to (elem->old_dof_object->dof_number(sys_num,v,i),
				                   DofObject::invalid_id);

		      di.push_back(elem->old_dof_object->dof_number(sys_num,v,i));
		    }
		}
	      else
		{
		  libmesh_assert(!elem->active() || fe_type.family == LAGRANGE ||
				 elem->refinement_flag() == Elem::JUST_COARSENED);
		  di.resize(di.size() + nc, DofObject::invalid_id);
		}
	    }
	}
    } // end loop over variables

  // Finally append any SCALAR dofs that we asked for.
  std::vector<dof_id_type> di_new;
  std::vector<unsigned int>::iterator it           = SCALAR_var_numbers.begin();
  std::vector<unsigned int>::const_iterator it_end = SCALAR_var_numbers.end();
  for( ; it != it_end; ++it)
  {
    this->SCALAR_dof_indices(di_new,*it,true);
    di.insert( di.end(), di_new.begin(), di_new.end());
  }

#ifdef DEBUG
  libmesh_assert_equal_to (tot_size, di.size());
#endif

  STOP_LOG("old_dof_indices()", "DofMap");
}



/*
void DofMap::augment_send_list_for_projection(const MeshBase& mesh)
{
  // Loop over the active local elements in the mesh.
  // If the element was just refined and its parent lives
  // on a different processor then we need to augment the
  // _send_list with the parent's DOF indices so that we
  // can access the parent's data for computing solution
  // projections, etc...

  // The DOF indices for the parent
  std::vector<dof_id_type> di;

  // Flag telling us if we need to re-sort the send_list.
  // Note we won't need to re-sort unless a child with a
  // parent belonging to a different processor is found
  bool needs_sorting = false;


  MeshBase::const_element_iterator       elem_it  = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      const Elem* const elem   = *elem_it;

      // We only need to consider the children that
      // were just refined
      if (elem->refinement_flag() != Elem::JUST_REFINED) continue;

      const Elem* const parent = elem->parent();

      // If the parent lives on another processor
      // than the child
      if (parent != NULL)
	if (parent->processor_id() != elem->processor_id())
	  {
	    // Get the DOF indices for the parent
	    this->dof_indices (parent, di);

	    // Insert the DOF indices into the send list
	    _send_list.insert (_send_list.end(),
			       di.begin(), di.end());

	    // We will need to re-sort the send list
	    needs_sorting = true;
	  }
    }

  // The send-list might need to be sorted again.
  if (needs_sorting) this->sort_send_list ();
}
*/

#endif // LIBMESH_ENABLE_AMR


#ifdef LIBMESH_ENABLE_CONSTRAINTS

void DofMap::find_connected_dofs (std::vector<dof_id_type>& elem_dofs) const
{
  typedef std::set<dof_id_type> RCSet;

  // First insert the DOFS we already depend on into the set.
  RCSet dof_set (elem_dofs.begin(), elem_dofs.end());

  bool done = true;

  // Next insert any dofs those might be constrained in terms
  // of.  Note that in this case we may not be done:  Those may
  // in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  for (unsigned int i=0; i<elem_dofs.size(); i++)
    if (this->is_constrained_dof(elem_dofs[i]))
      {
	// If the DOF is constrained
	DofConstraints::const_iterator
	  pos = _dof_constraints.find(elem_dofs[i]);

	libmesh_assert (pos != _dof_constraints.end());

	const DofConstraintRow& constraint_row = pos->second;

// adaptive p refinement currently gives us lots of empty constraint
// rows - we should optimize those DoFs away in the future.  [RHS]
//	libmesh_assert (!constraint_row.empty());

	DofConstraintRow::const_iterator it     = constraint_row.begin();
	DofConstraintRow::const_iterator it_end = constraint_row.end();


	// Add the DOFs this dof is constrained in terms of.
	// note that these dofs might also be constrained, so
	// we will need to call this function recursively.
	for ( ; it != it_end; ++it)
	  if (!dof_set.count (it->first))
	    {
	      dof_set.insert (it->first);
	      done = false;
	    }
	}


  // If not done then we need to do more work
  // (obviously :-) )!
  if (!done)
    {
      // Fill the vector with the contents of the set
      elem_dofs.clear();
      elem_dofs.insert (elem_dofs.end(),
			dof_set.begin(), dof_set.end());


      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      this->find_connected_dofs (elem_dofs);

    } // end if (!done)
}

#endif // LIBMESH_ENABLE_CONSTRAINTS



#if defined(__GNUC__) && (__GNUC__ < 4) && !defined(__INTEL_COMPILER)

void SparsityPattern::_dummy_function(void)
{
}

#endif



void SparsityPattern::Build::operator()(const ConstElemRange &range)
{
  // Compute the sparsity structure of the global matrix.  This can be
  // fed into a PetscMatrix to allocate exacly the number of nonzeros
  // necessary to store the matrix.  This algorithm should be linear
  // in the (# of elements)*(# nodes per element)
  const processor_id_type proc_id           = mesh.processor_id();
  const dof_id_type n_dofs_on_proc    = dof_map.n_dofs_on_processor(proc_id);
  const dof_id_type first_dof_on_proc = dof_map.first_dof(proc_id);
  const dof_id_type end_dof_on_proc   = dof_map.end_dof(proc_id);

  sparsity_pattern.resize(n_dofs_on_proc);

  // If the user did not explicitly specify the DOF coupling
  // then all the DOFS are coupled to each other.  Furthermore,
  // we can take a shortcut and do this more quickly here.  So
  // we use an if-test.
  if ((dof_coupling == NULL) || (dof_coupling->empty()))
    {
      std::vector<dof_id_type>
	element_dofs,
	neighbor_dofs,
	dofs_to_add;

      std::vector<const Elem*> active_neighbors;

      for (ConstElemRange::const_iterator elem_it = range.begin() ; elem_it != range.end(); ++elem_it)
	{
	  const Elem* const elem = *elem_it;

	  // Get the global indices of the DOFs with support on this element
	  dof_map.dof_indices (elem, element_dofs);
#ifdef LIBMESH_ENABLE_CONSTRAINTS
	  dof_map.find_connected_dofs (element_dofs);
#endif

	  // We can be more efficient if we sort the element DOFs
	  // into increasing order
	  std::sort(element_dofs.begin(), element_dofs.end());

	  const unsigned int n_dofs_on_element =
	    libmesh_cast_int<unsigned int>(element_dofs.size());

	  for (unsigned int i=0; i<n_dofs_on_element; i++)
	    {
	      const dof_id_type ig = element_dofs[i];

	      SparsityPattern::Row *row;

	      // We save non-local row components for now so we can
	      // communicate them to other processors later.

	      if ((ig >= first_dof_on_proc) &&
	          (ig <  end_dof_on_proc))
		{
		  // This is what I mean
		  // libmesh_assert_greater_equal ((ig - first_dof_on_proc), 0);
		  // but do the test like this because ig and
		  // first_dof_on_proc are unsigned ints
		  libmesh_assert_greater_equal (ig, first_dof_on_proc);
		  libmesh_assert_less ((ig - first_dof_on_proc), sparsity_pattern.size());

		  row = &sparsity_pattern[ig - first_dof_on_proc];
                }
              else
                {
                  row = &nonlocal_pattern[ig];
                }

              // If the row is empty we will add *all* the element DOFs,
              // so just do that.
              if (row->empty())
                {
                   row->insert(row->end(),
                               element_dofs.begin(),
                               element_dofs.end());
                }
              else
                {
                  // Build a list of the DOF indices not found in the
                  // sparsity pattern
                  dofs_to_add.clear();

                  // Cache iterators.  Low will move forward, subsequent
                  // searches will be on smaller ranges
                  SparsityPattern::Row::iterator
                    low  = std::lower_bound (row->begin(), row->end(), element_dofs.front()),
                    high = std::upper_bound (low,          row->end(), element_dofs.back());

                  for (unsigned int j=0; j<n_dofs_on_element; j++)
                    {
                      const dof_id_type jg = element_dofs[j];

                      // See if jg is in the sorted range
                      std::pair<SparsityPattern::Row::iterator,
                                SparsityPattern::Row::iterator>
                        pos = std::equal_range (low, high, jg);

                      // Must add jg if it wasn't found
                      if (pos.first == pos.second)
                        dofs_to_add.push_back(jg);

                      // pos.first is now a valid lower bound for any
                      // remaining element DOFs. (That's why we sorted them.)
                      // Use it for the next search
                      low = pos.first;
                    }

                  // Add to the sparsity pattern
                  if (!dofs_to_add.empty())
                    {
                      const std::size_t old_size = row->size();

                      row->insert (row->end(),
                                   dofs_to_add.begin(),
                                   dofs_to_add.end());

		      SparsityPattern::sort_row
                        (row->begin(), row->begin()+old_size, row->end());
                    }
                }

              // Now (possibly) add dofs from neighboring elements
              // TODO:[BSK] optimize this like above!
              if (implicit_neighbor_dofs)
                for (unsigned int s=0; s<elem->n_sides(); s++)
                  if (elem->neighbor(s) != NULL)
                    {
                      const Elem* const neighbor_0 = elem->neighbor(s);
#ifdef LIBMESH_ENABLE_AMR
                      neighbor_0->active_family_tree_by_neighbor(active_neighbors,elem);
#else
                      active_neighbors.clear();
                      active_neighbors.push_back(neighbor_0);
#endif

                      for (std::size_t a=0; a != active_neighbors.size(); ++a)
                        {
                          const Elem *neighbor = active_neighbors[a];

                          dof_map.dof_indices (neighbor, neighbor_dofs);
#ifdef LIBMESH_ENABLE_CONSTRAINTS
                          dof_map.find_connected_dofs (neighbor_dofs);
#endif
                          const std::size_t n_dofs_on_neighbor = neighbor_dofs.size();

                          for (std::size_t j=0; j<n_dofs_on_neighbor; j++)
                            {
                              const dof_id_type jg = neighbor_dofs[j];

                              // See if jg is in the sorted range
                              std::pair<SparsityPattern::Row::iterator,
                                        SparsityPattern::Row::iterator>
                                pos = std::equal_range (row->begin(), row->end(), jg);

                              // Insert jg if it wasn't found
                              if (pos.first == pos.second)
                                row->insert (pos.first, jg);
                            }
                        }
                    }
	    }
	}
    }

  // This is what we do in the case that the user has specified
  // explicit DOF coupling.
  else
    {
      libmesh_assert(dof_coupling);
      libmesh_assert_equal_to (dof_coupling->size(),
		               dof_map.n_variables());

      const unsigned int n_var = dof_map.n_variables();

      std::vector<dof_id_type>
	element_dofs_i,
	element_dofs_j,
	neighbor_dofs,
	dofs_to_add;


      std::vector<const Elem*> active_neighbors;
      for (ConstElemRange::const_iterator elem_it = range.begin() ; elem_it != range.end(); ++elem_it)
	for (unsigned int vi=0; vi<n_var; vi++)
	  {
	    const Elem* const elem = *elem_it;

	    // Find element dofs for variable vi
	    dof_map.dof_indices (elem, element_dofs_i, vi);
#ifdef LIBMESH_ENABLE_CONSTRAINTS
	    dof_map.find_connected_dofs (element_dofs_i);
#endif

	    // We can be more efficient if we sort the element DOFs
	    // into increasing order
	    std::sort(element_dofs_i.begin(), element_dofs_i.end());
	    const unsigned int n_dofs_on_element_i =
	      libmesh_cast_int<unsigned int>(element_dofs_i.size());

	    for (unsigned int vj=0; vj<n_var; vj++)
	      if ((*dof_coupling)(vi,vj)) // If vi couples to vj
		{
		  // Find element dofs for variable vj, note that
		  // if vi==vj we already have the dofs.
		  if (vi != vj)
		    {
		      dof_map.dof_indices (elem, element_dofs_j, vj);
#ifdef LIBMESH_ENABLE_CONSTRAINTS
		      dof_map.find_connected_dofs (element_dofs_j);
#endif

		      // We can be more efficient if we sort the element DOFs
		      // into increasing order
		      std::sort (element_dofs_j.begin(), element_dofs_j.end());
		    }
		  else
		    element_dofs_j = element_dofs_i;

		  const unsigned int n_dofs_on_element_j =
		    libmesh_cast_int<unsigned int>(element_dofs_j.size());

		  // there might be 0 dofs for the other variable on the same element (when subdomain variables do not overlap) and that's when we do not do anything
		  if (n_dofs_on_element_j > 0)
		  {
		    for (unsigned int i=0; i<n_dofs_on_element_i; i++)
		    {
		      const dof_id_type ig = element_dofs_i[i];

	              SparsityPattern::Row *row;

	              // We save non-local row components for now so we can
	              // communicate them to other processors later.

		      if ((ig >= first_dof_on_proc) &&
			  (ig <  end_dof_on_proc))
			{
			  // This is what I mean
			  // libmesh_assert_greater_equal ((ig - first_dof_on_proc), 0);
			  // but do the test like this because ig and
			  // first_dof_on_proc are unsigned ints
			  libmesh_assert_greater_equal (ig, first_dof_on_proc);
			  libmesh_assert_less (ig, (sparsity_pattern.size() +
					first_dof_on_proc));

			  row = &sparsity_pattern[ig - first_dof_on_proc];
                        }
                      else
                        {
                          row = &nonlocal_pattern[ig];
                        }

                      // If the row is empty we will add *all* the element j DOFs,
                      // so just do that.
                      if (row->empty())
                        {
                          row->insert(row->end(),
                                      element_dofs_j.begin(),
                                      element_dofs_j.end());
                        }
                      else
                        {
                          // Build a list of the DOF indices not found in the
                          // sparsity pattern
                          dofs_to_add.clear();

                          // Cache iterators.  Low will move forward, subsequent
                          // searches will be on smaller ranges
                          SparsityPattern::Row::iterator
                            low  = std::lower_bound
                              (row->begin(), row->end(), element_dofs_j.front()),
                            high = std::upper_bound
                              (low,          row->end(), element_dofs_j.back());

                          for (unsigned int j=0; j<n_dofs_on_element_j; j++)
                            {
                              const dof_id_type jg = element_dofs_j[j];

                              // See if jg is in the sorted range
                              std::pair<SparsityPattern::Row::iterator,
                                        SparsityPattern::Row::iterator>
                                pos = std::equal_range (low, high, jg);

                              // Must add jg if it wasn't found
                              if (pos.first == pos.second)
                                dofs_to_add.push_back(jg);

                              // pos.first is now a valid lower bound for any
                              // remaining element j DOFs. (That's why we sorted them.)
                              // Use it for the next search
                              low = pos.first;
                            }

                          // Add to the sparsity pattern
                          if (!dofs_to_add.empty())
                            {
                              const std::size_t old_size = row->size();

                              row->insert (row->end(),
                                           dofs_to_add.begin(),
                                           dofs_to_add.end());

                              SparsityPattern::sort_row
				(row->begin(), row->begin()+old_size,
                                 row->end());
                            }
                        }
                       // Now (possibly) add dofs from neighboring elements
                       // TODO:[BSK] optimize this like above!
                       if (implicit_neighbor_dofs)
                         for (unsigned int s=0; s<elem->n_sides(); s++)
                           if (elem->neighbor(s) != NULL)
                             {
                               const Elem* const neighbor_0 = elem->neighbor(s);
#ifdef LIBMESH_ENABLE_AMR
                               neighbor_0->active_family_tree_by_neighbor(active_neighbors,elem);
#else
                               active_neighbors.clear();
                               active_neighbors.push_back(neighbor_0);
#endif

                               for (std::size_t a=0; a != active_neighbors.size(); ++a)
                                 {
                                   const Elem *neighbor = active_neighbors[a];

                                   dof_map.dof_indices (neighbor, neighbor_dofs);
#ifdef LIBMESH_ENABLE_CONSTRAINTS
                                   dof_map.find_connected_dofs (neighbor_dofs);
#endif
                                   const unsigned int n_dofs_on_neighbor =
				     libmesh_cast_int<unsigned int>(neighbor_dofs.size());

                                   for (unsigned int j=0; j<n_dofs_on_neighbor; j++)
                                     {
                                       const dof_id_type jg = neighbor_dofs[j];

                                       // See if jg is in the sorted range
                                       std::pair<SparsityPattern::Row::iterator,
                                                 SparsityPattern::Row::iterator>
                                         pos = std::equal_range (row->begin(), row->end(), jg);

                                       // Insert jg if it wasn't found
                                       if (pos.first == pos.second)
                                         row->insert (pos.first, jg);
                                     }
                                 }
                             }
		    }
		  }
		}
	  }
    }

  // Now a new chunk of sparsity structure is built for all of the
  // DOFs connected to our rows of the matrix.

  // If we're building a full sparsity pattern, then we've got
  // complete rows to work with, so we can just count them from
  // scratch.
  if (need_full_sparsity_pattern)
    {
      n_nz.clear();
      n_oz.clear();
    }

  n_nz.resize (n_dofs_on_proc, 0);
  n_oz.resize (n_dofs_on_proc, 0);

  for (dof_id_type i=0; i<n_dofs_on_proc; i++)
    {
      // Get the row of the sparsity pattern
      SparsityPattern::Row &row = sparsity_pattern[i];

      for (dof_id_type j=0; j<row.size(); j++)
        if ((row[j] < first_dof_on_proc) || (row[j] >= end_dof_on_proc))
          n_oz[i]++;
        else
          n_nz[i]++;

      // If we're not building a full sparsity pattern, then we want
      // to avoid overcounting these entries as much as possible.
      if (!need_full_sparsity_pattern)
        row.clear();
    }
}



void SparsityPattern::Build::join (const SparsityPattern::Build &other)
{
  const processor_id_type proc_id           = mesh.processor_id();
  const dof_id_type       n_global_dofs     = dof_map.n_dofs();
  const dof_id_type       n_dofs_on_proc    = dof_map.n_dofs_on_processor(proc_id);
  const dof_id_type       first_dof_on_proc = dof_map.first_dof(proc_id);
  const dof_id_type       end_dof_on_proc   = dof_map.end_dof(proc_id);

  libmesh_assert_equal_to (sparsity_pattern.size(), other.sparsity_pattern.size());
  libmesh_assert_equal_to (n_nz.size(), sparsity_pattern.size());
  libmesh_assert_equal_to (n_oz.size(), sparsity_pattern.size());

  for (dof_id_type r=0; r<n_dofs_on_proc; r++)
    {
      // increment the number of on and off-processor nonzeros in this row
      // (note this will be an upper bound unless we need the full sparsity pattern)
      if (need_full_sparsity_pattern)
	{
	  SparsityPattern::Row       &my_row    = sparsity_pattern[r];
	  const SparsityPattern::Row &their_row = other.sparsity_pattern[r];

	  // simple copy if I have no dofs
	  if (my_row.empty())
	    my_row = their_row;

	  // otherwise add their DOFs to mine, resort, and re-unique the row
	  else if (!their_row.empty()) // do nothing for the trivial case where
	    {                          // their row is empty
	      my_row.insert (my_row.end(),
			     their_row.begin(),
			     their_row.end());

	      // We cannot use SparsityPattern::sort_row() here because it expects
	      // the [begin,middle) [middle,end) to be non-overlapping.  This is not
	      // necessarily the case here, so use std::sort()
	      std::sort (my_row.begin(), my_row.end());

	      my_row.erase(std::unique (my_row.begin(), my_row.end()), my_row.end());
	    }

	  // fix the number of on and off-processor nonzeros in this row
	  n_nz[r] = n_oz[r] = 0;

	  for (std::size_t j=0; j<my_row.size(); j++)
	    if ((my_row[j] < first_dof_on_proc) || (my_row[j] >= end_dof_on_proc))
	      n_oz[r]++;
	    else
	      n_nz[r]++;
	}
      else
        {
          n_nz[r] += other.n_nz[r];
	  n_nz[r] = std::min(n_nz[r], n_dofs_on_proc);
          n_oz[r] += other.n_oz[r];
	  n_oz[r] =std::min(n_oz[r], static_cast<dof_id_type>(n_global_dofs-n_nz[r]));
	}
    }

  // Move nonlocal row information to ourselves; the other thread
  // won't need it in the map after that.
  NonlocalGraph::const_iterator it = other.nonlocal_pattern.begin();
  for (; it != other.nonlocal_pattern.end(); ++it)
    {
      const dof_id_type dof_id = it->first;

#ifndef NDEBUG
      processor_id_type dbg_proc_id = 0;
      while (dof_id >= dof_map.end_dof(dbg_proc_id))
        dbg_proc_id++;
      libmesh_assert (dbg_proc_id != this->processor_id());
#endif

      const SparsityPattern::Row &their_row = it->second;

      // We should have no empty values in a map
      libmesh_assert (!their_row.empty());

      NonlocalGraph::iterator my_it = nonlocal_pattern.find(it->first);
      if (my_it == nonlocal_pattern.end())
        {
//          nonlocal_pattern[it->first].swap(their_row);
          nonlocal_pattern[it->first] = their_row;
        }
      else
        {
          SparsityPattern::Row &my_row = my_it->second;

	  my_row.insert (my_row.end(),
			 their_row.begin(),
			 their_row.end());

	  // We cannot use SparsityPattern::sort_row() here because it expects
	  // the [begin,middle) [middle,end) to be non-overlapping.  This is not
	  // necessarily the case here, so use std::sort()
	  std::sort (my_row.begin(), my_row.end());

	  my_row.erase(std::unique (my_row.begin(), my_row.end()), my_row.end());
        }
    }
}



void SparsityPattern::Build::parallel_sync ()
{
  parallel_object_only();
  this->comm().verify(need_full_sparsity_pattern);

  const dof_id_type n_global_dofs   = dof_map.n_dofs();
  const dof_id_type n_dofs_on_proc  = dof_map.n_dofs_on_processor(this->processor_id());
  const dof_id_type local_first_dof = dof_map.first_dof();
  const dof_id_type local_end_dof   = dof_map.end_dof();

  // Trade sparsity rows with other processors
  for (processor_id_type p=1; p != this->n_processors(); ++p)
    {
      // Push to processor procup while receiving from procdown
      processor_id_type procup = (this->processor_id() + p) %
                                  this->n_processors();
      processor_id_type procdown = (this->n_processors() +
                                    this->processor_id() - p) %
                                    this->n_processors();

      // Pack the sparsity pattern rows to push to procup
      std::vector<dof_id_type> pushed_row_ids,
                               pushed_row_ids_to_me;
      std::vector<std::vector<dof_id_type> > pushed_rows,
                                             pushed_rows_to_me;

      // Move nonlocal row information to a structure to send it from;
      // we don't need it in the map after that.
      NonlocalGraph::iterator it = nonlocal_pattern.begin();
      while (it != nonlocal_pattern.end())
        {
          const dof_id_type dof_id = it->first;
          processor_id_type proc_id = 0;
          while (dof_id >= dof_map.end_dof(proc_id))
            proc_id++;

          libmesh_assert (proc_id != this->processor_id());

          if (proc_id == procup)
            {
              pushed_row_ids.push_back(dof_id);

	      // We can't just do the swap trick here, thanks to the
	      // differing vector allocators?
              pushed_rows.push_back(std::vector<dof_id_type>());
	      pushed_rows.back().assign
                (it->second.begin(), it->second.end());

              nonlocal_pattern.erase(it++);
            }
          else
            ++it;
        }

      this->comm().send_receive(procup, pushed_row_ids,
					procdown, pushed_row_ids_to_me);
      this->comm().send_receive(procup, pushed_rows,
					procdown, pushed_rows_to_me);
      pushed_row_ids.clear();
      pushed_rows.clear();

      const std::size_t n_rows = pushed_row_ids_to_me.size();
      for (std::size_t i=0; i != n_rows; ++i)
        {
          const dof_id_type r = pushed_row_ids_to_me[i];
          const dof_id_type my_r = r - local_first_dof;

          std::vector<dof_id_type> &their_row = pushed_rows_to_me[i];

          if (need_full_sparsity_pattern)
            {
              SparsityPattern::Row &my_row =
                sparsity_pattern[my_r];

              // They wouldn't have sent an empty row
              libmesh_assert(!their_row.empty());

              // We can end up with an empty row on a dof that touches our
              // inactive elements but not our active ones
              if (my_row.empty())
                {
                  my_row.assign (their_row.begin(),
                                 their_row.end());
                }
              else
                {
                  my_row.insert (my_row.end(),
                                 their_row.begin(),
                                 their_row.end());

                  // We cannot use SparsityPattern::sort_row() here because it expects
                  // the [begin,middle) [middle,end) to be non-overlapping.  This is not
                  // necessarily the case here, so use std::sort()
                  std::sort (my_row.begin(), my_row.end());

                  my_row.erase(std::unique (my_row.begin(), my_row.end()), my_row.end());
                }

              // fix the number of on and off-processor nonzeros in this row
              n_nz[my_r] = n_oz[my_r] = 0;

              for (std::size_t j=0; j<my_row.size(); j++)
                if ((my_row[j] < local_first_dof) || (my_row[j] >= local_end_dof))
                  n_oz[my_r]++;
                else
                  n_nz[my_r]++;
            }
          else
            {
              for (std::size_t j=0; j<their_row.size(); j++)
	        if ((their_row[j] < local_first_dof) || (their_row[j] >= local_end_dof))
                  n_oz[my_r]++;
                else
                  n_nz[my_r]++;

              n_nz[my_r] = std::min(n_nz[my_r], n_dofs_on_proc);
              n_oz[my_r] = std::min(n_oz[my_r],
				    static_cast<dof_id_type>(n_global_dofs-n_nz[my_r]));
            }
        }
    }

    // We should have sent everything at this point.
    libmesh_assert (nonlocal_pattern.empty());
}



void DofMap::print_info(std::ostream& os) const
{
  os << this->get_info();
}



std::string DofMap::get_info() const
{
  std::ostringstream os;

  // If we didn't calculate the exact sparsity pattern, the threaded
  // sparsity pattern assembly may have just given us an upper bound
  // on sparsity.
  const char* may_equal = " <= ";

  // If we calculated the exact sparsity pattern, then we can report
  // exact bandwidth figures:
  std::vector<SparseMatrix<Number>* >::const_iterator
    pos = _matrices.begin(),
    end = _matrices.end();

  for (; pos != end; ++pos)
    if ((*pos)->need_full_sparsity_pattern())
      may_equal = " = ";

  dof_id_type max_n_nz = 0, max_n_oz = 0;
  long double avg_n_nz = 0, avg_n_oz = 0;

  if (_n_nz)
    {
      for (std::size_t i = 0; i != _n_nz->size(); ++i)
        {
          max_n_nz = std::max(max_n_nz, (*_n_nz)[i]);
          avg_n_nz += (*_n_nz)[i];
        }

      std::size_t n_nz_size = _n_nz->size();

      this->comm().max(max_n_nz);
      this->comm().sum(avg_n_nz);
      this->comm().sum(n_nz_size);

      avg_n_nz /= std::max(n_nz_size,std::size_t(1));

      libmesh_assert(_n_oz);

      for (std::size_t i = 0; i != (*_n_oz).size(); ++i)
        {
          max_n_oz = std::max(max_n_oz, (*_n_oz)[i]);
          avg_n_oz += (*_n_oz)[i];
        }

      std::size_t n_oz_size = _n_oz->size();

      this->comm().max(max_n_oz);
      this->comm().sum(avg_n_oz);
      this->comm().sum(n_oz_size);

      avg_n_oz /= std::max(n_oz_size,std::size_t(1));
    }

  os << "    DofMap Sparsity\n      Average  On-Processor Bandwidth"
     << may_equal << avg_n_nz << '\n';

  os << "      Average Off-Processor Bandwidth"
     << may_equal << avg_n_oz << '\n';

  os << "      Maximum  On-Processor Bandwidth"
     << may_equal << max_n_nz << '\n';

  os << "      Maximum Off-Processor Bandwidth"
     << may_equal << max_n_oz << std::endl;

#ifdef LIBMESH_ENABLE_CONSTRAINTS

  std::size_t n_constraints = 0, max_constraint_length = 0,
              n_rhss = 0;
  long double avg_constraint_length = 0.;

  for (DofConstraints::const_iterator it=_dof_constraints.begin();
       it != _dof_constraints.end(); ++it)
    {
      // Only count local constraints, then sum later
      const dof_id_type constrained_dof = it->first;
      if (constrained_dof < this->first_dof() ||
          constrained_dof >= this->end_dof())
        continue;

      const DofConstraintRow& row = it->second;
      std::size_t rowsize = row.size();

      max_constraint_length = std::max(max_constraint_length,
                                       rowsize);
      avg_constraint_length += rowsize;
      n_constraints++;

      if (_primal_constraint_values.count(constrained_dof))
        n_rhss++;
    }

  this->comm().sum(n_constraints);
  this->comm().sum(n_rhss);
  this->comm().sum(avg_constraint_length);
  this->comm().max(max_constraint_length);

  os << "    DofMap Constraints\n      Number of DoF Constraints = "
     << n_constraints;
  if (n_rhss)
    os << '\n'
       << "      Number of Heterogenous Constraints= " << n_rhss;
  if (n_constraints)
    {
      avg_constraint_length /= n_constraints;

      os << '\n'
         << "      Average DoF Constraint Length= " << avg_constraint_length;
    }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  std::size_t n_node_constraints = 0, max_node_constraint_length = 0,
              n_node_rhss = 0;
  long double avg_node_constraint_length = 0.;

  for (NodeConstraints::const_iterator it=_node_constraints.begin();
       it != _node_constraints.end(); ++it)
    {
      // Only count local constraints, then sum later
      const Node *node = it->first;
      if (node->processor_id() != this->processor_id())
        continue;

      const NodeConstraintRow& row = it->second.first;
      std::size_t rowsize = row.size();

      max_node_constraint_length = std::max(max_node_constraint_length,
                                            rowsize);
      avg_node_constraint_length += rowsize;
      n_node_constraints++;

      if (it->second.second != Point(0))
        n_node_rhss++;
    }

  this->comm().sum(n_node_constraints);
  this->comm().sum(n_node_rhss);
  this->comm().sum(avg_node_constraint_length);
  this->comm().max(max_node_constraint_length);

  os << "\n      Number of Node Constraints = " << n_node_constraints;
  if (n_node_rhss)
    os << '\n'
       << "      Number of Heterogenous Node Constraints= " << n_node_rhss;
  if (n_node_constraints)
    {
      avg_node_constraint_length /= n_node_constraints;
      os << "\n      Maximum Node Constraint Length= " << max_node_constraint_length
         << '\n'
         << "      Average Node Constraint Length= " << avg_node_constraint_length;
    }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  os << std::endl;

#endif // LIBMESH_ENABLE_CONSTRAINTS

  return os.str();
}

} // namespace libMesh
