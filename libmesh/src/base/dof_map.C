// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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

// Local Includes -----------------------------------
#include "coupling_matrix.h"
#include "dense_matrix.h"
#include "dense_vector_base.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_interface.h"
#include "fe_type.h"
#include "fe_base.h" // FEBase::build() for continuity test
#include "libmesh_logging.h"
#include "mesh_base.h"
#include "mesh_tools.h"
#include "numeric_vector.h"
#include "parallel.h"
#include "sparse_matrix.h"
#include "string_to_enum.h"
#include "threads_allocators.h"


// ------------------------------------------------------------
// DofMap member functions

// Destructor
DofMap::~DofMap()
{
  this->clear();
}



void DofMap::add_variable (const System::Variable &var)
{
  _variables.push_back (var);
}



const System::Variable & DofMap::variable (const unsigned int c) const
{
  libmesh_assert (c < _variables.size());

  return _variables[c];
}



Order DofMap::variable_order (const unsigned int c) const
{
  libmesh_assert (c < _variables.size());

  return _variables[c].type().order;
}



const FEType& DofMap::variable_type (const unsigned int c) const
{
  libmesh_assert (c < _variables.size());

  return _variables[c].type();
}



void DofMap::attach_matrix (SparseMatrix<Number>& matrix)
{
  _matrices.push_back(&matrix);
  
  matrix.attach_dof_map (*this);
}



DofObject* DofMap::node_ptr(MeshBase& mesh, unsigned int i) const
{
  return mesh.node_ptr(i);
}



DofObject* DofMap::elem_ptr(MeshBase& mesh, unsigned int i) const
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
  parallel_only();

  // First, iterate over local objects to find out how many
  // are on each processor
  std::vector<unsigned int>
    ghost_objects_from_proc(libMesh::n_processors(), 0);

  iterator_type it  = objects_begin;

  for (; it != objects_end; ++it)
    {
      DofObject *obj = *it;

      if (obj)
        {
          unsigned int obj_procid = obj->processor_id();
          // We'd better be completely partitioned by now
          libmesh_assert(obj_procid != DofObject::invalid_processor_id);
          ghost_objects_from_proc[obj_procid]++;
        }
    }

  std::vector<unsigned int> objects_on_proc(libMesh::n_processors(), 0);
  Parallel::allgather(ghost_objects_from_proc[libMesh::processor_id()],
                      objects_on_proc);

#ifdef DEBUG
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    libmesh_assert(ghost_objects_from_proc[p] <= objects_on_proc[p]);
#endif

  // Request sets to send to each processor
  std::vector<std::vector<unsigned int> >
    requested_ids(libMesh::n_processors());

  // We know how many of our objects live on each processor, so
  // reserve() space for requests from each.
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      requested_ids[p].reserve(ghost_objects_from_proc[p]);

  for (it = objects_begin; it != objects_end; ++it)
    {
      DofObject *obj = *it;
      if (obj->processor_id() != DofObject::invalid_processor_id)
        requested_ids[obj->processor_id()].push_back(obj->id());
    }
#ifdef DEBUG
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    libmesh_assert(requested_ids[p].size() == ghost_objects_from_proc[p]);
#endif

  // Next set ghost object n_comps from other processors
  for (unsigned int p=1; p != libMesh::n_processors(); ++p)
    {
      // Trade my requests with processor procup and procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();
      std::vector<unsigned int> request_to_fill;
      Parallel::send_receive(procup, requested_ids[procup],
                             procdown, request_to_fill);

      // Fill those requests
      const unsigned int n_variables = this->n_variables();

      std::vector<unsigned int> ghost_data
        (request_to_fill.size() * 2 * n_variables);

      for (unsigned int i=0; i != request_to_fill.size(); ++i)
        {
          DofObject *requested = (this->*objects)(mesh, request_to_fill[i]);
          libmesh_assert(requested);
          libmesh_assert(requested->processor_id() == libMesh::processor_id());
          libmesh_assert(requested->n_vars(this->sys_number()) == n_variables);
          for (unsigned int v=0; v != n_variables; ++v)
            {
              unsigned int n_comp =
                requested->n_comp(this->sys_number(), v);
              ghost_data[i*2*n_variables+v] = n_comp;
              unsigned int first_dof = n_comp ?
                requested->dof_number(this->sys_number(), v, 0) : 0;
              libmesh_assert(first_dof != DofObject::invalid_id);
              ghost_data[i*2*n_variables+n_variables+v] = first_dof;
            }
        }

      // Trade back the results
      std::vector<unsigned int> filled_request;
      Parallel::send_receive(procdown, ghost_data,
                             procup, filled_request);

      // And copy the id changes we've now been informed of
      libmesh_assert(filled_request.size() ==
             requested_ids[procup].size() * 2 * n_variables);
      for (unsigned int i=0; i != requested_ids[procup].size(); ++i)
        {
          DofObject *requested = (this->*objects)(mesh, requested_ids[procup][i]);
          libmesh_assert(requested);
          libmesh_assert (requested->processor_id() == procup);
          for (unsigned int v=0; v != n_variables; ++v)
            {
              unsigned int n_comp = filled_request[i*2*n_variables+v];
              requested->set_n_comp(this->sys_number(), v, n_comp);
              if (n_comp)
                {
                  unsigned int first_dof = 
                    filled_request[i*2*n_variables+n_variables+v];
                  libmesh_assert(first_dof != DofObject::invalid_id);
                  requested->set_dof_number
                    (this->sys_number(), v, 0, first_dof);

		  // don't forget to add these remote dofs to the _send_list
		  for (unsigned int comp=0; comp!=n_comp; ++comp)
		    _send_list.push_back(first_dof+comp);
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
      unsigned int n_variables = obj->n_vars(this->sys_number());
      for (unsigned int v=0; v != n_variables; ++v)
        {
          unsigned int n_comp =
            obj->n_comp(this->sys_number(), v);
          unsigned int first_dof = n_comp ?
            obj->dof_number(this->sys_number(), v, 0) : 0;
          libmesh_assert(first_dof != DofObject::invalid_id);
        }
    }
#endif
}



void DofMap::reinit(MeshBase& mesh)
{
  libmesh_assert (mesh.is_prepared());
  
  START_LOG("reinit()", "DofMap");
  
  //this->clear();

  const unsigned int n_var = this->n_variables();

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
	libmesh_assert ((*node_it)->old_dof_object == NULL);
      }
    
    MeshBase::element_iterator       elem_it  = mesh.elements_begin();
    const MeshBase::element_iterator elem_end = mesh.elements_end();
    
    for ( ; elem_it != elem_end; ++elem_it)
      {
	(*elem_it)->clear_old_dof_object();
	libmesh_assert ((*elem_it)->old_dof_object == NULL);
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
	      if (node->has_dofs())
		node->set_old_dof_object();
	  }

	libmesh_assert (elem->old_dof_object == NULL);

	if (elem->has_dofs())
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
      (*node_it)->set_n_vars(this->sys_number(),n_var);

    // All the elements
    MeshBase::element_iterator       elem_it  = mesh.elements_begin();
    const MeshBase::element_iterator elem_end = mesh.elements_end();

    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_vars(this->sys_number(),n_var);
  }


  // Zero _n_SCALAR_dofs, it will be updated below.
  this->_n_SCALAR_dofs = 0;

  //------------------------------------------------------------
  // Next allocate space for the DOF indices
  for (unsigned int var=0; var<this->n_variables(); var++)
    {
      const System::Variable &var_description =	this->variable(var);
      const FEType& base_fe_type              = this->variable_type(var);

      // Skip dof counting for a SCALAR variable
      if(base_fe_type.family == SCALAR)
      {
        // TODO: This only works in serial at the moment
        if(libMesh::n_processors() > 1)
          libmesh_not_implemented();

        this->_n_SCALAR_dofs += base_fe_type.order;
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
	  libmesh_assert (elem != NULL);
	  
	  // Skip the numbering if this variable is 
	  // not active on this element's subdoman
	  if (!var_description.active_on_subdomain(elem->subdomain_id()))
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
                  std::cerr << "ERROR: Finite element "
                    << Utility::enum_to_string(base_fe_type.family)
                    << " on geometric element "
                    << Utility::enum_to_string(type) << std::endl
                    << "only supports FEInterface::max_order = " 
                    << FEInterface::max_order(base_fe_type,type)
                    << ", not fe_type.order = " << base_fe_type.order
                    << std::endl;

                  libmesh_error();
                }

              std::cerr << "WARNING: Finite element "
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
	            node->n_comp(this->sys_number(), var);

		  const unsigned int vertex_dofs =
		    std::max(FEInterface::n_dofs_at_node(dim, fe_type,
                                                         type, n),
                             old_node_dofs);
		  
		  // Some discontinuous FEs have no vertex dofs
		  if (vertex_dofs > old_node_dofs)
		    {
		      node->set_n_comp(this->sys_number(), var,
				       vertex_dofs);
		      // Abusing dof_number to set a "this is a
		      // vertex" flag
		      node->set_dof_number(this->sys_number(),
					   var, 0, vertex_dofs);
		    }
	        }
	    }
	} // done counting vertex dofs

      // count edge & face dofs next
      elem_it = mesh.active_elements_begin();

      for ( ; elem_it != elem_end; ++elem_it)
	{
	  Elem* elem = *elem_it;
	  libmesh_assert (elem != NULL);
	  
	  // Skip the numbering if this variable is 
	  // not active on this element's subdoman
	  if (!var_description.active_on_subdomain(elem->subdomain_id()))
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
	        node->n_comp(this->sys_number(), var);

              const unsigned int vertex_dofs = old_node_dofs?
                node->dof_number (this->sys_number(),var,0):0;

	      const unsigned int new_node_dofs =
		FEInterface::n_dofs_at_node(dim, fe_type, type, n);

	      // We've already allocated vertex DOFs
	      if (elem->is_vertex(n))
	        {
		  libmesh_assert(old_node_dofs >= vertex_dofs);
		  libmesh_assert(vertex_dofs >= new_node_dofs);
		}
	      // We need to allocate the rest
	      else
	        {
		  // If this has no dofs yet, it needs no vertex
		  // dofs, so we just give it edge or face dofs
		  if (!old_node_dofs)
                    {
		      node->set_n_comp(this->sys_number(), var,
				       new_node_dofs);
		      // Abusing dof_number to set a "this has no
		      // vertex dofs" flag
                      if (new_node_dofs)
		        node->set_dof_number(this->sys_number(),
					     var, 0, 0);
                    }

		  // If this has dofs, but has no vertex dofs,
		  // it may still need more edge or face dofs if
                  // we're p-refined.
		  else if (vertex_dofs == 0)
		    {
                      if (new_node_dofs > old_node_dofs)
                        {
		          node->set_n_comp(this->sys_number(), var,
				           new_node_dofs);
		          node->set_dof_number(this->sys_number(),
					       var, 0, vertex_dofs);
                        }
		    }
		  // If this is another element's vertex,
		  // add more (non-overlapping) edge/face dofs if
		  // necessary
		  else if (extra_hanging_dofs)
		    {
                      if (new_node_dofs > old_node_dofs - vertex_dofs)
                        {
		          node->set_n_comp(this->sys_number(), var,
				           vertex_dofs + new_node_dofs);
		          node->set_dof_number(this->sys_number(),
					       var, 0, vertex_dofs);
                        }
		    }
		  // If this is another element's vertex, add any
		  // (overlapping) edge/face dofs if necessary
		  else
		    {
		      libmesh_assert(old_node_dofs >= vertex_dofs);
                      if (new_node_dofs > old_node_dofs)
                        {
		          node->set_n_comp(this->sys_number(), var,
				           new_node_dofs);
		          node->set_dof_number(this->sys_number(),
					       var, 0, vertex_dofs);
                        }
		    }
		}
	    }
          // Allocate the element DOFs
	  const unsigned int dofs_per_elem =
			  FEInterface::n_dofs_per_elem(dim, fe_type,
						       type);

	  elem->set_n_comp(this->sys_number(), var, dofs_per_elem);

	}
    }

  // Calling DofMap::reinit() by itself makes little sense,
  // so we won't bother with nonlocal DofObjects.
  // Those will be fixed by distribute_dofs
/*
  //------------------------------------------------------------
  // At this point, all n_comp and dof_number values on local
  // DofObjects should be correct, but a ParallelMesh might have
  // incorrect values on non-local DofObjects.  Let's request the
  // correct values from each other processor.

  this->set_nonlocal_n_comps(mesh.nodes_begin(),
                             mesh.nodes_end(),
                             mesh,
                             &DofMap::node_ptr);

  this->set_nonlocal_n_comps(mesh.elements_begin(),
                             mesh.elements_end(),
                             mesh,
                             &DofMap::elem_ptr);
*/

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
  _first_df.clear();
  _end_df.clear();
  _var_first_local_df.clear();
  _send_list.clear();
  _n_nz.clear();
  _n_oz.clear();


#ifdef LIBMESH_ENABLE_AMR

  _dof_constraints.clear();
  _n_old_dfs = 0;

#endif

  _matrices.clear();

  _n_dfs = 0;
}



void DofMap::distribute_dofs (MeshBase& mesh)
{
  // This function must be run on all processors at once
  parallel_only();

  // Log how long it takes to distribute the degrees of freedom
  START_LOG("distribute_dofs()", "DofMap");

  libmesh_assert (mesh.is_prepared());

  const unsigned int proc_id = libMesh::processor_id();
  const unsigned int n_proc  = libMesh::n_processors();
  
//  libmesh_assert (this->n_variables() > 0);
  libmesh_assert (proc_id < n_proc);
  
  // re-init in case the mesh has changed
  this->reinit(mesh);
  
  // By default distribute variables in a
  // var-major fashion, but allow run-time
  // specification
  bool node_major_dofs = libMesh::on_command_line ("--node_major_dofs");

  // The DOF counter, will be incremented as we encounter
  // new degrees of freedom
  unsigned int next_free_dof = 0;

  // Set temporary DOF indices on this processor
  if (node_major_dofs)
    this->distribute_local_dofs_node_major (next_free_dof, mesh);
  else
    this->distribute_local_dofs_var_major (next_free_dof, mesh);

  // Get DOF counts on all processors
  std::vector<unsigned int> dofs_on_proc(n_proc, 0);
  Parallel::allgather(next_free_dof, dofs_on_proc);

  // Resize the _first_df and _end_df arrays
  _first_df.resize(n_proc);
  _end_df.resize (n_proc);

  // Get DOF offsets
  _first_df[0] = 0;
  for (unsigned int i=1; i < n_proc; ++i)
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

  libmesh_assert(next_free_dof == _end_df[proc_id]);

  //------------------------------------------------------------
  // At this point, all n_comp and dof_number values on local
  // DofObjects should be correct, but a ParallelMesh might have
  // incorrect values on non-local DofObjects.  Let's request the
  // correct values from each other processor.

  if (libMesh::n_processors() > 1)
    {
      this->set_nonlocal_dof_objects(mesh.nodes_begin(),
                                     mesh.nodes_end(),
                                     mesh, &DofMap::node_ptr);

      this->set_nonlocal_dof_objects(mesh.elements_begin(),
                                     mesh.elements_end(),
                                     mesh, &DofMap::elem_ptr);
    }
  
  // Set the total number of degrees of freedom
#ifdef LIBMESH_ENABLE_AMR
  _n_old_dfs = _n_dfs;
#endif
  _n_dfs = _end_df[n_proc-1];

  STOP_LOG("distribute_dofs()", "DofMap");

  _send_list.clear();

  // Note that in the add_neighbors_to_send_list nodes on processor
  // boundaries that are shared by multiple elements are added for
  // each element.
  this->add_neighbors_to_send_list(mesh);
  
  // Here we used to clean up that data structure; now System and
  // EquationSystems call that for us, after we've added constraint
  // dependencies to the send_list too.
  // this->sort_send_list ();
}


void DofMap::distribute_local_dofs_node_major(unsigned int &next_free_dof,
                                              MeshBase& mesh)
{
  const unsigned int sys_num = this->sys_number();
  const unsigned int n_vars  = this->n_variables();

  // We now only add remote dofs to the _send_list
  // unsigned int send_list_size = 0;

  // _var_first_local_df does not work with node_major dofs
  _var_first_local_df.resize(n_vars+1);
  std::fill (_var_first_local_df.begin(),
	     _var_first_local_df.end(),
	     DofObject::invalid_id);

  // Set _n_SCALAR_dofs
  this->_n_SCALAR_dofs = 0;
  for (unsigned var=0; var<n_vars; var++)
  {
    if( this->variable(var).type().family == SCALAR )
    {
      // TODO: This only works in serial at the moment
      if(libMesh::n_processors() > 1)
        libmesh_not_implemented();

      this->_n_SCALAR_dofs += this->variable(var).type().order;
      continue;
    }
  }
  
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

          for (unsigned var=0; var<n_vars; var++)
          {
	      if( (this->variable(var).type().family != SCALAR) &&
                  (this->variable(var).active_on_subdomain(elem->subdomain_id())) )
	      {
		// assign dof numbers (all at once) if this is
		// our node and if they aren't already there
		if ((node->n_comp(sys_num,var) > 0) &&
		    (node->processor_id() == libMesh::processor_id()) &&
		    (node->dof_number(sys_num,var,0) ==
		     DofObject::invalid_id))
		  {
		    node->set_dof_number(sys_num,
					 var,
					 0,
					 next_free_dof);
		    next_free_dof += node->n_comp(sys_num,var);
		  }
	      }
          }
        }

      // Now number the element DOFS
      for (unsigned var=0; var<n_vars; var++)
	if ( (this->variable(var).type().family != SCALAR) &&
             (this->variable(var).active_on_subdomain(elem->subdomain_id())) )
	  if (elem->n_comp(sys_num,var) > 0)
	    {
	      libmesh_assert (elem->dof_number(sys_num,var,0) ==
			      DofObject::invalid_id);

	      elem->set_dof_number(sys_num,
				   var,
				   0,
				   next_free_dof);

	      next_free_dof += elem->n_comp(sys_num,var);
	    }
    } // done looping over elements

#ifdef DEBUG
// Make sure we didn't miss any nodes
  MeshTools::libmesh_assert_valid_node_procids(mesh);

  MeshBase::node_iterator       node_it  = mesh.local_nodes_begin();
  const MeshBase::node_iterator node_end = mesh.local_nodes_end();
  for (; node_it != node_end; ++node_it)
    {
      Node *obj = *node_it;
      libmesh_assert(obj);
      unsigned int n_variables = obj->n_vars(this->sys_number());
      for (unsigned int v=0; v != n_variables; ++v)
        {
          unsigned int n_comp =
            obj->n_comp(this->sys_number(), v);
          unsigned int first_dof = n_comp ?
            obj->dof_number(this->sys_number(), v, 0) : 0;
          libmesh_assert(first_dof != DofObject::invalid_id);
        }
    }
#endif // DEBUG
}



void DofMap::distribute_local_dofs_var_major(unsigned int &next_free_dof,
                                             MeshBase& mesh)
{
  const unsigned int sys_num = this->sys_number();
  const unsigned int n_vars  = this->n_variables();

  // We now only add remote dofs to the _send_list
  // unsigned int send_list_size = 0;

  // We will cache the first local index for each variable
  _var_first_local_df.clear();

  // Zero _n_SCALAR_dofs, we will set it below
  this->_n_SCALAR_dofs = 0;

  //-------------------------------------------------------------------------
  // First count and assign temporary numbers to local dofs
  for (unsigned var=0; var<n_vars; var++)
    {
      _var_first_local_df.push_back(next_free_dof);

      const System::Variable var_description = this->variable(var);

      // Skip dof counting for a SCALAR variable
      if(var_description.type().family == SCALAR)
      {
        // TODO: This only works in serial at the moment
        if(libMesh::n_processors() > 1)
          libmesh_not_implemented();

        this->_n_SCALAR_dofs += var_description.type().order;
        continue;
      }

      MeshBase::element_iterator       elem_it  = mesh.active_local_elements_begin();
      const MeshBase::element_iterator elem_end = mesh.active_local_elements_end();

      for ( ; elem_it != elem_end; ++elem_it)
        {
          // Only number dofs connected to active
          // elements on this processor.
          Elem* elem  = *elem_it;

	  // ... and only variables which are active on
	  // on this element's subdomain
	  if (!var_description.active_on_subdomain(elem->subdomain_id()))
	    continue;

          const unsigned int n_nodes = elem->n_nodes();
          
          // First number the nodal DOFS
          for (unsigned int n=0; n<n_nodes; n++)
            {
              Node* node = elem->get_node(n);
              
              // assign dof numbers (all at once) if this is
              // our node and if they aren't already there
              if ((node->n_comp(sys_num,var) > 0) &&
                  (node->processor_id() == libMesh::processor_id()) &&
                  (node->dof_number(sys_num,var,0) ==
                   DofObject::invalid_id))
                {
                  node->set_dof_number(sys_num,
                                       var,
                                       0,
                                       next_free_dof);
                  next_free_dof += node->n_comp(sys_num,var);
                }
            }
                  
          // Now number the element DOFS
          if (elem->n_comp(sys_num,var) > 0)
            {
              libmesh_assert (elem->dof_number(sys_num,var,0) ==
                      DofObject::invalid_id);

              elem->set_dof_number(sys_num,
                                   var,
                                   0,
                                   next_free_dof);

              next_free_dof += elem->n_comp(sys_num,var);
            }
        } // end loop on elements
    } // end loop on variables

  // Cache the last local dof number too
  _var_first_local_df.push_back(next_free_dof);

#ifdef DEBUG
  // Make sure we didn't miss any nodes
  MeshTools::libmesh_assert_valid_node_procids(mesh);

  MeshBase::node_iterator       node_it  = mesh.local_nodes_begin();
  const MeshBase::node_iterator node_end = mesh.local_nodes_end();
  for (; node_it != node_end; ++node_it)
    {
      Node *obj = *node_it;
      libmesh_assert(obj);
      unsigned int n_variables = obj->n_vars(this->sys_number());
      for (unsigned int v=0; v != n_variables; ++v)
        {
          unsigned int n_comp =
            obj->n_comp(this->sys_number(), v);
          unsigned int first_dof = n_comp ?
            obj->dof_number(this->sys_number(), v, 0) : 0;
          libmesh_assert(first_dof != DofObject::invalid_id);
        }
    }
#endif // DEBUG
}



void DofMap::add_neighbors_to_send_list(MeshBase& mesh)
{
  START_LOG("add_neighbors_to_send_list()", "DofMap");
  
  //-------------------------------------------------------------------------
  // We need to add the DOFs from elements that live on neighboring processors
  // that are neighbors of the elements on the local processor
  //-------------------------------------------------------------------------

  MeshBase::const_element_iterator       local_elem_it
    = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator local_elem_end
    = mesh.active_local_elements_end(); 

  std::vector<bool> node_on_processor(mesh.max_node_id(), false);
  std::vector<unsigned int> di;
  std::vector<const Elem *> family;

  // Loop over the active local elements, adding all active elements
  // that neighbor an active local element to the send list.
  for ( ; local_elem_it != local_elem_end; ++local_elem_it)
    {
      const Elem* elem = *local_elem_it;

      // Flag all the nodes of active local elements as seen
      for (unsigned int n=0; n!=elem->n_nodes(); n++)
        node_on_processor[elem->node(n)] = true;

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

            for (unsigned int i=0; i!=family.size(); ++i)
	      // If the neighbor lives on a different processor
	      if (family[i]->processor_id() != libMesh::processor_id())
		{
		  // Get the DOF indices for this neighboring element
		  this->dof_indices (family[i], di);

	          // Insert the remote DOF indices into the send list
                  for (unsigned int j=0; j != di.size(); ++j)
                    if (di[j] < this->first_dof() ||
                        di[j] >= this->end_dof())
	              _send_list.push_back(di[j]);
	        }
	  }
    }

  // Now loop over all non_local active elements and add any missing
  // nodal-only neighbors
  MeshBase::const_element_iterator       elem_it
    = mesh.active_elements_begin();
  const MeshBase::const_element_iterator elem_end
    = mesh.active_elements_end();
  
  for ( ; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;

      // If this is one of our elements, we've already added it
      if (elem->processor_id() == libMesh::processor_id())
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
            for (unsigned int j=0; j != di.size(); ++j)
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
  
  // First sort the send list.  After this
  // duplicated elements will be adjacent in the
  // vector
  std::sort(_send_list.begin(), _send_list.end());

  // Now use std::unique to remove duplicate entries  
  std::vector<unsigned int>::iterator new_end =
    std::unique (_send_list.begin(), _send_list.end());

  // Remove the end of the send_list.  Use the "swap trick"
  // from Effective STL
  std::vector<unsigned int> (_send_list.begin(), new_end).swap (_send_list);
  
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
      if (FEBase::build (mesh.mesh_dimension(),
			 this->variable_type(var))->get_continuity() !=  DISCONTINUOUS)
	all_discontinuous_dofs = false;

    if (all_discontinuous_dofs)
      implicit_neighbor_dofs = true;
  }
  
  return implicit_neighbor_dofs;
}



void DofMap::compute_sparsity(const MeshBase& mesh)
{
  libmesh_assert (mesh.is_prepared());
  libmesh_assert (this->n_variables());

  START_LOG("compute_sparsity()", "DofMap");

  // Compute the sparsity structure of the global matrix.  This can be
  // fed into a PetscMatrix to allocate exacly the number of nonzeros
  // necessary to store the matrix.  This algorithm should be linear
  // in the (# of elements)*(# nodes per element)

  // We can be more efficient in the threaded sparsity pattern assembly
  // if we don't need the exact pattern.  For some sparse matrix formats
  // a good upper bound will suffice.
  bool need_full_sparsity_pattern=false;
  std::vector<SparseMatrix<Number>* >::iterator
    pos = _matrices.begin(),
    end = _matrices.end();
      
  for (; pos != end; ++pos)
    if ((*pos)->need_full_sparsity_pattern())
      need_full_sparsity_pattern = true;

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
  SparsityPattern::Build sp (mesh,
			     *this,
			     _dof_coupling,
			     implicit_neighbor_dofs,
			     need_full_sparsity_pattern);
  
  Threads::parallel_reduce (ConstElemRange (mesh.active_elements_begin(),
					    mesh.active_elements_end()), sp);

#ifndef NDEBUG
  // Avoid declaring these variables unless asserts are enabled.
  const unsigned int proc_id        = mesh.processor_id();
  const unsigned int n_dofs_on_proc = this->n_dofs_on_processor(proc_id);
#endif
  libmesh_assert (sp.sparsity_pattern.size() == n_dofs_on_proc);

  // steal the n_nz and n_oz arrays from sp -- it won't need them any more,
  // and this is more efficient than copying them.
  _n_nz.swap(sp.n_nz);
  _n_oz.swap(sp.n_oz);
  
  STOP_LOG("compute_sparsity()", "DofMap");
  
  // We are done with the sparsity_pattern.  However, quite a
  // lot has gone into computing it.  It is possible that some
  // \p SparseMatrix implementations want to see it.  Let them
  // see it before we throw it away.
  pos = _matrices.begin();
  end = _matrices.end();
      
  for (; pos != end; ++pos)
    (*pos)->update_sparsity_pattern (sp.sparsity_pattern);     
}



void DofMap::extract_local_vector (const NumericVector<Number>& Ug,
				   const std::vector<unsigned int>& dof_indices,
				   DenseVectorBase<Number>& Ue) const
{
#ifdef LIBMESH_ENABLE_AMR

  // Trivial mapping
  libmesh_assert (dof_indices.size() == Ue.size());
  bool has_constrained_dofs = false;

  for (unsigned int il=0; il<dof_indices.size(); il++)
    {
      const unsigned int ig = dof_indices[il];

      if (this->is_constrained_dof (ig)) has_constrained_dofs = true;
      
      libmesh_assert ((il >= Ug.first_local_index()) &&
	      (il <  Ug.last_local_index()));

      Ue.el(il) = Ug(ig);
    }

  // If the element has any constrained DOFs then we need
  // to account for them in the mapping.  This will handle
  // the case that the input vector is not constrained.
  if (has_constrained_dofs)
    {
      // Copy the input DOF indices.
      std::vector<unsigned int> constrained_dof_indices(dof_indices);

      DenseMatrix<Number> C;

      this->build_constraint_matrix (C, constrained_dof_indices);

      libmesh_assert (dof_indices.size()             == C.m());
      libmesh_assert (constrained_dof_indices.size() == C.n());

      // zero-out Ue
      Ue.zero();

      // compute Ue = C Ug, with proper mapping.
      for (unsigned int i=0; i<dof_indices.size(); i++)
	for (unsigned int j=0; j<constrained_dof_indices.size(); j++)
	  {
	    const unsigned int jg = constrained_dof_indices[j];

	    libmesh_assert ((jg >= Ug.first_local_index()) &&
		    (jg <  Ug.last_local_index()));
	    
	    Ue.el(i) += C(i,j)*Ug(jg); 	    
	  }
    }  
   
#else
  
  // Trivial mapping
  libmesh_assert (dof_indices.size() == Ue.size());
  
  for (unsigned int il=0; il<dof_indices.size(); il++)
    {
      const unsigned int ig = dof_indices[il];
      
      libmesh_assert ((ig >= Ug.first_local_index()) && (ig <  Ug.last_local_index()));

      Ue.el(il) = Ug(ig);
    }
  
#endif
}



void DofMap::dof_indices (const Elem* const elem,
			  std::vector<unsigned int>& di,
			  const unsigned int vn) const
{
  START_LOG("dof_indices()", "DofMap");
  
  libmesh_assert (elem != NULL);
  
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

  // Create a vector to store the variable numbers
  // of SCALAR variables that we've requested
  std::vector<unsigned int> SCALAR_var_numbers;
  SCALAR_var_numbers.clear();

  // Get the dof numbers
  for (unsigned int v=0; v<n_vars; v++)
    if ((v == vn) || (vn == libMesh::invalid_uint))
    {
      if(this->variable(v).type().family == SCALAR)
      {
        // TODO: This only works in serial at the moment
        if(libMesh::n_processors() > 1)
          libmesh_not_implemented();

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
			libmesh_assert (node->dof_number(sys_num,v,i) !=
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
		    libmesh_assert (node->dof_number(sys_num,v,i) !=
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
		      libmesh_assert (elem->dof_number(sys_num,v,i) !=
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
  // First we need to find out the first dof index for each
  // SCALAR.
  unsigned int first_SCALAR_dof_index = n_dofs() - n_SCALAR_dofs();
  std::map<unsigned int, unsigned int> SCALAR_first_dof_index;
  SCALAR_first_dof_index.clear();

  // iterate over _all_ of the SCALARs and store each one's first dof index
  for (unsigned int v=0; v<n_vars; v++)
    if(this->variable(v).type().family == SCALAR)
    {
      unsigned int current_n_SCALAR_dofs = this->variable(v).type().order;
      SCALAR_first_dof_index.insert(
        std::pair<unsigned int, unsigned int>(v,first_SCALAR_dof_index) );
      first_SCALAR_dof_index += current_n_SCALAR_dofs;
    }

  for(unsigned int i=0; i<SCALAR_var_numbers.size(); i++)
  {
    unsigned int var_number = SCALAR_var_numbers[i];

    // Now use current_SCALAR_var_number to index into
    // SCALAR_first_dof_index
    std::map<unsigned int, unsigned int>::const_iterator iter =
      SCALAR_first_dof_index.find(var_number);

#ifdef DEBUG
    libmesh_assert(iter != SCALAR_first_dof_index.end());
#endif

    unsigned int current_first_SCALAR_dof_index = iter->second;

    // Also, get the number of SCALAR dofs from the variable order
    unsigned int current_n_SCALAR_dofs = this->variable(var_number).type().order;

    for(unsigned int j=0; j<current_n_SCALAR_dofs; j++)
    {
      unsigned int index = current_first_SCALAR_dof_index+j;
      di.push_back(index);
    }
  }

#ifdef DEBUG
  libmesh_assert (tot_size == di.size());
#endif
  
  STOP_LOG("dof_indices()", "DofMap");  
}
 


#ifdef LIBMESH_ENABLE_AMR

void DofMap::old_dof_indices (const Elem* const elem,
			      std::vector<unsigned int>& di,
			      const unsigned int vn) const
{
  START_LOG("old_dof_indices()", "DofMap");
  
  libmesh_assert (elem != NULL);
  libmesh_assert (elem->old_dof_object != NULL);
	    

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

  // Create a vector to store the variable numbers
  // of SCALAR variables that we've requested
  std::vector<unsigned int> SCALAR_var_numbers;
  SCALAR_var_numbers.clear();
  
  // Get the dof numbers
  for (unsigned int v=0; v<n_vars; v++)
    if ((v == vn) || (vn == libMesh::invalid_uint))
    {
      if(this->variable(v).type().family == SCALAR)
        {
          // TODO: This only works in serial at the moment
          if(libMesh::n_processors() > 1)
            libmesh_not_implemented();

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
	      libmesh_assert (elem->p_level() > 0);
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
	      libmesh_assert (node->old_dof_object != NULL);
	    
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
			libmesh_assert (node->old_dof_object->dof_number(sys_num,v,i) !=
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
		    libmesh_assert (node->old_dof_object->dof_number(sys_num,v,i) !=
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
		  libmesh_assert (elem->old_dof_object != NULL);
		
		  for (unsigned int i=0; i<nc; i++)
		    {
		      libmesh_assert (elem->old_dof_object->dof_number(sys_num,v,i) !=
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
  // First we need to find out the first dof index for each
  // SCALAR.
  unsigned int first_SCALAR_dof_index = n_old_dofs() - n_SCALAR_dofs();
  std::map<unsigned int, unsigned int> SCALAR_first_dof_index;
  SCALAR_first_dof_index.clear();

  // iterate over _all_ of the SCALARs and store each one's first dof index
  for (unsigned int v=0; v<n_vars; v++)
    if(this->variable(v).type().family == SCALAR)
    {
      unsigned int current_n_SCALAR_dofs = this->variable(v).type().order;
      SCALAR_first_dof_index.insert(
        std::pair<unsigned int, unsigned int>(v,first_SCALAR_dof_index) );
      first_SCALAR_dof_index += current_n_SCALAR_dofs;
    }

  for(unsigned int i=0; i<SCALAR_var_numbers.size(); i++)
  {
    unsigned int var_number = SCALAR_var_numbers[i];

    // Now use current_SCALAR_var_number to index into
    // SCALAR_first_dof_index
    std::map<unsigned int, unsigned int>::const_iterator iter =
      SCALAR_first_dof_index.find(var_number);

#ifdef DEBUG
    libmesh_assert(iter != SCALAR_first_dof_index.end());
#endif

    unsigned int current_first_SCALAR_dof_index = iter->second;

    // Also, get the number of SCALAR dofs from the variable order
    unsigned int current_n_SCALAR_dofs = this->variable(var_number).type().order;

    for(unsigned int j=0; j<current_n_SCALAR_dofs; j++)
    {
      unsigned int index = current_first_SCALAR_dof_index+j;
      di.push_back(index);
    }
  }
  
#ifdef DEBUG
  libmesh_assert (tot_size == di.size());
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
  std::vector<unsigned int> di;

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


#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)

void DofMap::find_connected_dofs (std::vector<unsigned int>& elem_dofs) const
{
  typedef std::set<unsigned int> RCSet;

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

#endif // LIBMESH_ENABLE_AMR || LIBMESH_ENABLE_PERIODIC



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
  const unsigned int proc_id           = mesh.processor_id();
  const unsigned int n_dofs_on_proc    = dof_map.n_dofs_on_processor(proc_id);
  const unsigned int first_dof_on_proc = dof_map.first_dof(proc_id);
  const unsigned int end_dof_on_proc   = dof_map.end_dof(proc_id);

  sparsity_pattern.resize(n_dofs_on_proc);
  
  // If the user did not explicitly specify the DOF coupling
  // then all the DOFS are coupled to each other.  Furthermore,
  // we can take a shortcut and do this more quickly here.  So
  // we use an if-test.
  if ((dof_coupling == NULL) || (dof_coupling->empty()))
    {
      std::vector<unsigned int>
	element_dofs,
	neighbor_dofs,
	dofs_to_add;

      std::vector<const Elem*> active_neighbors;

      for (ConstElemRange::const_iterator elem_it = range.begin() ; elem_it != range.end(); ++elem_it)
	{
	  const Elem* const elem = *elem_it;

	  // Get the global indices of the DOFs with support on this element
	  dof_map.dof_indices (elem, element_dofs);
#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)
	  dof_map.find_connected_dofs (element_dofs);
#endif

	  // We can be more efficient if we sort the element DOFs
	  // into increasing order
	  std::sort(element_dofs.begin(), element_dofs.end());
	  
	  const unsigned int n_dofs_on_element = element_dofs.size();
	    
	  for (unsigned int i=0; i<n_dofs_on_element; i++)
	    {
	      const unsigned int ig = element_dofs[i];
	      
	      // Only bother if this matrix row will be stored
	      // on this processor.
	      if ((ig >= first_dof_on_proc) &&
		  (ig <  end_dof_on_proc))
		{
		  // This is what I mean
		  // libmesh_assert ((ig - first_dof_on_proc) >= 0);
		  // but do the test like this because ig and
		  // first_dof_on_proc are unsigned ints
		  libmesh_assert (ig >= first_dof_on_proc);
		  libmesh_assert ((ig - first_dof_on_proc) < sparsity_pattern.size());
		  
		  SparsityPattern::Row &row = sparsity_pattern[ig - first_dof_on_proc];

		  // If the row is empty we will add *all* the element DOFs,
		  // so just do that.
		  if (row.empty())
		    {
 		      row.insert(row.end(),
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
			low  = std::lower_bound (row.begin(), row.end(), element_dofs.front()),
			high = std::upper_bound (low,         row.end(), element_dofs.back());
		  
		      for (unsigned int j=0; j<n_dofs_on_element; j++)
			{
			  const unsigned int jg = element_dofs[j];
			  
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
			  const unsigned int old_size = row.size();
			  
			  row.insert (row.end(),
				      dofs_to_add.begin(),
				      dofs_to_add.end());
		      
			  SparsityPattern::sort_row (row.begin(), row.begin()+old_size, row.end());
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
			  
			  for (unsigned int a=0; a != active_neighbors.size(); ++a)
			    {
			      const Elem *neighbor = active_neighbors[a];
			      
			      dof_map.dof_indices (neighbor, neighbor_dofs);
#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)
			      dof_map.find_connected_dofs (neighbor_dofs);
#endif			      
			      const unsigned int n_dofs_on_neighbor = neighbor_dofs.size();
			      
			      for (unsigned int j=0; j<n_dofs_on_neighbor; j++)
				{
				  const unsigned int jg = neighbor_dofs[j];
				  
				  // See if jg is in the sorted range
				  std::pair<SparsityPattern::Row::iterator,
				            SparsityPattern::Row::iterator>
				    pos = std::equal_range (row.begin(), row.end(), jg);
			        
				  // Insert jg if it wasn't found
				  if (pos.first == pos.second)
				    row.insert (pos.first, jg);
				}    
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
      libmesh_assert (dof_coupling != NULL);
      libmesh_assert (dof_coupling->size() ==
		      dof_map.n_variables());
      
      const unsigned int n_var = dof_map.n_variables();
      
      std::vector<unsigned int>
	element_dofs_i,
	element_dofs_j,
	dofs_to_add;


      for (ConstElemRange::const_iterator elem_it = range.begin() ; elem_it != range.end(); ++elem_it)
	for (unsigned int vi=0; vi<n_var; vi++)
	  {
	    const Elem* const elem = *elem_it;
	    
	    // Find element dofs for variable vi
	    dof_map.dof_indices (elem, element_dofs_i, vi);
#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)
	    dof_map.find_connected_dofs (element_dofs_i);
#endif

	    // We can be more efficient if we sort the element DOFs
	    // into increasing order
	    std::sort(element_dofs_i.begin(), element_dofs_i.end());
	    const unsigned int n_dofs_on_element_i = element_dofs_i.size();
	    
	    for (unsigned int vj=0; vj<n_var; vj++)
	      if ((*dof_coupling)(vi,vj)) // If vi couples to vj
		{
		  // Find element dofs for variable vj, note that
		  // if vi==vj we already have the dofs.
		  if (vi != vj)
		    {
		      dof_map.dof_indices (elem, element_dofs_j, vj);
#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)
		      dof_map.find_connected_dofs (element_dofs_j);
#endif

		      // We can be more efficient if we sort the element DOFs
		      // into increasing order
		      std::sort (element_dofs_j.begin(), element_dofs_j.end());
		    }
		  else
		    element_dofs_j = element_dofs_i;
		  
		  const unsigned int n_dofs_on_element_j =
		    element_dofs_j.size();
		  
		  for (unsigned int i=0; i<n_dofs_on_element_i; i++)
		    {
		      const unsigned int ig = element_dofs_i[i];
		      
		      // We are only interested in computing the sparsity pattern
		      // for the rows in [range.begin(),range.end()).  So, if this
		      // element's DOFs are not in that range just go ahead to the
		      // next one.
		      //	      
		      // Also, only bother if this matrix row will be stored
		      // on this processor.
		      if ((ig >= first_dof_on_proc) &&
			  (ig <  end_dof_on_proc))			
			{
			  // This is what I mean
			  // libmesh_assert ((ig - first_dof_on_proc) >= 0);
			  // but do the test like this because ig and
			  // first_dof_on_proc are unsigned ints
			  libmesh_assert (ig >= first_dof_on_proc);
			  libmesh_assert (ig < (sparsity_pattern.size() +
					first_dof_on_proc));
			  
			  SparsityPattern::Row &row = sparsity_pattern[ig - first_dof_on_proc];

			  // If the row is empty we will add *all* the element j DOFs,
			  // so just do that.
			  if (row.empty())
			    {
			      row.insert(row.end(),
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
				low  = std::lower_bound (row.begin(), row.end(), element_dofs_j.front()),
				high = std::upper_bound (low,         row.end(), element_dofs_j.back());

			      for (unsigned int j=0; j<n_dofs_on_element_j; j++)
				{
				  const unsigned int jg = element_dofs_j[j];
				  
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
				  const unsigned int old_size = row.size();
				  
				  row.insert (row.end(),
					      dofs_to_add.begin(),
					      dofs_to_add.end());
		      
				  SparsityPattern::sort_row (row.begin(), row.begin()+old_size, row.end());
				}
			    }
			}
		    }
		}
	  }
    }
  
  // Now the full sparsity structure is built for all of the
  // DOFs connected to our rows of the matrix.
  n_nz.resize (n_dofs_on_proc);
  n_oz.resize (n_dofs_on_proc);

  // First zero the counters.
  std::fill(n_nz.begin(), n_nz.end(), 0);
  std::fill(n_oz.begin(), n_oz.end(), 0);
  
  for (unsigned int i=0; i<n_dofs_on_proc; i++)
    {
      // Get the row of the sparsity pattern
      const SparsityPattern::Row &row = sparsity_pattern[i];

      for (unsigned int j=0; j<row.size(); j++)
	if ((row[j] < first_dof_on_proc) || (row[j] >= end_dof_on_proc))
	  n_oz[i]++;
	else
	  n_nz[i]++;
    }
}



void SparsityPattern::Build::join (const SparsityPattern::Build &other)
{
  const unsigned int proc_id           = mesh.processor_id();
  const unsigned int n_global_dofs     = dof_map.n_dofs();
  const unsigned int n_dofs_on_proc    = dof_map.n_dofs_on_processor(proc_id);
  const unsigned int first_dof_on_proc = dof_map.first_dof(proc_id);
  const unsigned int end_dof_on_proc   = dof_map.end_dof(proc_id);

  libmesh_assert (sparsity_pattern.size() ==  other.sparsity_pattern.size());
  libmesh_assert (n_nz.size() == sparsity_pattern.size());
  libmesh_assert (n_oz.size() == sparsity_pattern.size());
  
  for (unsigned int r=0; r<n_dofs_on_proc; r++)
    {
      // incriment the number of on and off-processor nonzeros in this row
      // (note this will be an upper bound unless we need the full sparsity pattern)
      n_nz[r] += other.n_nz[r];  n_nz[r] = std::min(n_nz[r], n_dofs_on_proc);
      n_oz[r] += other.n_oz[r];  n_oz[r] = std::min(n_oz[r], n_global_dofs-n_nz[r]);

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
	  
	  for (unsigned int j=0; j<my_row.size(); j++)
	    if ((my_row[j] < first_dof_on_proc) || (my_row[j] >= end_dof_on_proc))
	      n_oz[r]++;
	    else
	      n_nz[r]++;
	}
    }
}
