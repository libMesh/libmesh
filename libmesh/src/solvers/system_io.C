// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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


#include "libmesh_common.h"
#include "parallel.h"

// C++ Includes
#include <cstdio> // for std::sprintf
#include <set>

// Local Includes
#include "system.h"
#include "mesh.h"
#include "mesh_tools.h"
#include "elem.h"
#include "xdr_cxx.h"
#include "numeric_vector.h"



// Anonymous namespace for implementation details.
namespace {
  
  // Comments:
  // ---------
  // - The io_blksize governs how many nodes or elements will be treated as a single block
  //   when performing parallel IO.
  // - This parameter only loosely affects the size of the actual IO buffer as this depends
  //   on the number of components a given variable has for the nodes/elements in the block.
  // - When reading/writing each processor uses an ID map which is 3*io_blksize*sizeof(unsigned int) bytes
  //   long, so if io_blksize=256000 we would expect that buffer alone to be ~3Mb.
  // - In general, an increase in io_blksize should increase the efficiency of the parallel
  //   read/writes by reducing the number of MPI messages at the expense of memory.
  // - If the library exhausts memory during IO you might reduce this parameter.
  
  const unsigned int io_blksize = 256000;
}



// ------------------------------------------------------------
// System class implementation
void System::read_header (Xdr& io,
			  const bool read_header,
			  const bool read_additional_data)
{
  // This method implements the input of a
  // System object, embedded in the output of
  // an EquationSystems<T_sys>.  This warrants some 
  // documentation.  The output file essentially
  // consists of 5 sections:
  //
  // for this system
  //
  //   5.) The number of variables in the system (unsigned int)
  //
  //   for each variable in the system
  //     
  //     6.) The name of the variable (string)
  //     
  //     7.) Combined in an FEType:
  //         - The approximation order(s) of the variable 
  //           (Order Enum, cast to int/s)
  //         - The finite element family/ies of the variable 
  //           (FEFamily Enum, cast to int/s)
  // 
  //   end variable loop
  //
  //   8.) The number of additional vectors (unsigned int),      
  //
  //     for each additional vector in the system object
  // 
  //     9.) the name of the additional vector  (string)
  //
  // end system
  assert (io.reading());
  
  // Possibly clear data structures and start from scratch.
  if (read_header)
    this->clear ();
  
  {
    // 5.) 
    // Read the number of variables in the system
    unsigned int n_vars=0;
    if (libMesh::processor_id() == 0) io.data (n_vars);
    Parallel::broadcast(n_vars);
      
    for (unsigned int var=0; var<n_vars; var++)
      {	            
	// 6.)
	// Read the name of the var-th variable
	std::string var_name;	  
	if (libMesh::processor_id() == 0) io.data (var_name);
	Parallel::broadcast(var_name);
	      	
	// 7.)
	// Read the approximation order(s) of the var-th variable 
	int order=0;	  
	if (libMesh::processor_id() == 0) io.data (order);
	Parallel::broadcast(order);
	
#ifdef ENABLE_INFINITE_ELEMENTS
	
	// do the same for radial_order
	int rad_order=0;	  
	if (libMesh::processor_id() == 0) io.data(rad_order);
	Parallel::broadcast(rad_order);

#endif

	// Read the finite element type of the var-th variable 
	int fam=0;
	if (libMesh::processor_id() == 0) io.data (fam);
	Parallel::broadcast(fam);
	FEType type;
	type.order  = static_cast<Order>(order);
	type.family = static_cast<FEFamily>(fam);
	
#ifdef ENABLE_INFINITE_ELEMENTS
	
	// Read additional information for infinite elements	
	int radial_fam=0;
	int i_map=0;
	
	if (libMesh::processor_id() == 0) io.data (radial_fam);
	Parallel::broadcast(radial_fam);
	if (libMesh::processor_id() == 0) io.data (i_map);
	Parallel::broadcast(i_map);
	
	type.radial_order  = static_cast<Order>(rad_order);
	type.radial_family = static_cast<FEFamily>(radial_fam);
	type.inf_map       = static_cast<InfMapType>(i_map);	  

#endif

	if (read_header) 
	  this->add_variable (var_name, type);
      }
  }

  // 8.)  
  // Read the number of additional vectors.  
  unsigned int n_vectors=0;  
  if (libMesh::processor_id() == 0) io.data (n_vectors);
  Parallel::broadcast(n_vectors);
  
  // If n_vectors > 0, this means that write_additional_data
  // was true when this file was written.  We will need to
  // make use of this fact later.
  if (n_vectors > 0)
    this->_additional_data_written = true;  
  
  for (unsigned int vec=0; vec<n_vectors; vec++)
    {
      // 9.)
      // Read the name of the vec-th additional vector
      std::string vec_name;      
      if (libMesh::processor_id() == 0) io.data (vec_name);
      Parallel::broadcast(vec_name);
      
      if (read_additional_data)
	{
	  // sanity checks
	  assert(this->_can_add_vectors);
	  assert(this->_vectors.count(vec_name) == 0);

	  this->add_vector(vec_name);
	}
    }
}



void System::read_legacy_data (Xdr& io,
			       const bool read_additional_data)
{
  deprecated();
  
  // This method implements the output of the vectors
  // contained in this System object, embedded in the 
  // output of an EquationSystems<T_sys>. 
  //
  //   10.) The global solution vector, re-ordered to be node-major 
  //       (More on this later.)                                    
  //                                                                
  //      for each additional vector in the object          
  //                                                                
  //      11.) The global additional vector, re-ordered to be       
  //           node-major (More on this later.)
  assert (io.reading());

  // read and reordering buffers
  std::vector<Number> global_vector;
  std::vector<Number> reordered_vector;

  // 10.)
  // Read and set the solution vector
  {	
    if (libMesh::processor_id() == 0) io.data (global_vector);	  
    Parallel::broadcast(global_vector);
    
    // Remember that the stored vector is node-major.
    // We need to put it into whatever application-specific
    // ordering we may have using the dof_map.
    reordered_vector.resize(global_vector.size());

    //std::cout << "global_vector.size()=" << global_vector.size() << std::endl;
    //std::cout << "this->n_dofs()=" << this->n_dofs() << std::endl;
    
    assert (global_vector.size() == this->n_dofs());
	
    unsigned int cnt=0;

    const unsigned int sys     = this->number();
    const unsigned int n_vars  = this->n_vars();

    for (unsigned int var=0; var<n_vars; var++)
      {	
	// First reorder the nodal DOF values
	{
	  MeshBase::node_iterator
	    it  = this->get_mesh().nodes_begin(),
	    end = this->get_mesh().nodes_end();
	
  	  for (; it != end; ++it)
	    for (unsigned int index=0; index<(*it)->n_comp(sys,var); index++)
	      {
	        assert ((*it)->dof_number(sys, var, index) !=
			DofObject::invalid_id);
		
                assert (cnt < global_vector.size());
                
                reordered_vector[(*it)->dof_number(sys, var, index)] =
	  	  global_vector[cnt++]; 
	      }
	}
	
	// Then reorder the element DOF values
	{
	  MeshBase::element_iterator
	    it  = this->get_mesh().active_elements_begin(),
	    end = this->get_mesh().active_elements_end();

	  for (; it != end; ++it)
	    for (unsigned int index=0; index<(*it)->n_comp(sys,var); index++)
	      {  
		assert ((*it)->dof_number(sys, var, index) !=
			DofObject::invalid_id);
		
		assert (cnt < global_vector.size());
		
		reordered_vector[(*it)->dof_number(sys, var, index)] =
		  global_vector[cnt++]; 
	      }
	}
      }
	    
    *(this->solution) = reordered_vector;
  }

  // For each additional vector, simply go through the list.
  // ONLY attempt to do this IF additional data was actually
  // written to the file for this system (controlled by the
  // _additional_data_written flag).  
  if (this->_additional_data_written)
    {
      std::map<std::string, NumericVector<Number>* >::iterator
	pos = this->_vectors.begin();
  
      for (; pos != this->_vectors.end(); ++pos)
	{
	  // 11.) 	   
	  // Read the values of the vec-th additional vector.
	  // Prior do _not_ clear, but fill with zero, since the
	  // additional vectors _have_ to have the same size
	  // as the solution vector
	  std::fill (global_vector.begin(), global_vector.end(), libMesh::zero);

	  if (libMesh::processor_id() == 0) io.data (global_vector);
	  Parallel::broadcast(global_vector);

	  // If read_additional_data==true, then we will keep this vector, otherwise
	  // we are going to throw it away.
	  if (read_additional_data)
	    {
	      // Remember that the stored vector is node-major.
	      // We need to put it into whatever application-specific
	      // ordering we may have using the dof_map.
	      std::fill (reordered_vector.begin(),
			 reordered_vector.end(),
			 libMesh::zero);
	
	      reordered_vector.resize(global_vector.size());	

	      assert (global_vector.size() == this->n_dofs());
	
	      unsigned int cnt=0;

	      const unsigned int sys     = this->number();
	      const unsigned int n_vars  = this->n_vars();
	
	      for (unsigned int var=0; var<n_vars; var++)
		{
		  // First reorder the nodal DOF values
		  {
		    MeshBase::node_iterator
		      it  = this->get_mesh().nodes_begin(),
		      end = this->get_mesh().nodes_end();

		    for (; it!=end; ++it)
		      for (unsigned int index=0; index<(*it)->n_comp(sys,var); index++)
			{
			  assert ((*it)->dof_number(sys, var, index) !=
				  DofObject::invalid_id);
			  
			  assert (cnt < global_vector.size());
			  
			  reordered_vector[(*it)->dof_number(sys, var, index)] =
			    global_vector[cnt++]; 
			}
		  }

		  // Then reorder the element DOF values
		  {
		    MeshBase::element_iterator
		      it  = this->get_mesh().active_elements_begin(),
		      end = this->get_mesh().active_elements_end();

		    for (; it!=end; ++it)
		      for (unsigned int index=0; index<(*it)->n_comp(sys,var); index++)
			{  
			  assert ((*it)->dof_number(sys, var, index) !=
				  DofObject::invalid_id);
			  
			  assert (cnt < global_vector.size());
			  
			  reordered_vector[(*it)->dof_number(sys, var, index)] =
			    global_vector[cnt++]; 
			}
		  }
		}
	    
	      // use the overloaded operator=(std::vector) to assign the values
	      *(pos->second) = reordered_vector;
	    }
	}
    } // end if (_additional_data_written)    
}



void System::read_parallel_data (Xdr &io,
				 const bool read_additional_data)
{
  /**
   * This method implements the output of the vectors
   * contained in this System object, embedded in the 
   * output of an EquationSystems<T_sys>. 
   *
   *   9.) The global solution vector, re-ordered to be node-major 
   *       (More on this later.)                                    
   *                                                                
   *      for each additional vector in the object          
   *                                                                
   *      10.) The global additional vector, re-ordered to be       
   *           node-major (More on this later.)
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */ 
  assert (io.reading());
  assert (io.is_open());
  
  std::vector<Number> io_buffer;
  
  // 9.)
  //
  // Actually read the solution components
  // for the ith system to disk
  io.data(io_buffer);
	  
  const unsigned int sys_num = this->number();
  const unsigned int n_vars  = this->n_vars();

  unsigned int cnt=0;
  
  // Loop over each variable and each node, and read out the value.
  for (unsigned int var=0; var<n_vars; var++)
    {
      // First read the node DOF values
      {	
	MeshBase::const_node_iterator
	  it  = this->get_mesh().local_nodes_begin(),
	  end = this->get_mesh().local_nodes_end();
	
	for (; it != end; ++it)
	  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
	    {
	      assert ((*it)->dof_number(sys_num, var, comp) !=
		      DofObject::invalid_id);
	      assert (cnt < io_buffer.size());
	      this->solution->set((*it)->dof_number(sys_num, var, comp), io_buffer[cnt++]);
	    }
      }

      // Then read the element DOF values
      {	
	MeshBase::const_element_iterator
	  it  = this->get_mesh().local_elements_begin(),
	  end = this->get_mesh().local_elements_end();
	
	for (; it != end; ++it)
	  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
	    {
	      assert ((*it)->dof_number(sys_num, var, comp) !=
		      DofObject::invalid_id);
	      assert (cnt < io_buffer.size());
	      this->solution->set((*it)->dof_number(sys_num, var, comp), io_buffer[cnt++]);
	    }
      }
    }
  
  // Only read additional vectors if wanted  
  if (read_additional_data)
    {	  
      std::map<std::string, NumericVector<Number>* >::const_iterator
	pos = _vectors.begin();
  
      for(; pos != this->_vectors.end(); ++pos)
        {
	  cnt=0;
	  io_buffer.clear();
	  
	  // 10.)
	  //
	  // Actually read the additional vector components
	  // for the ith system to disk
	  io.data(io_buffer);
	  
	  // Loop over each variable and each node, and read out the value.
	  for (unsigned int var=0; var<n_vars; var++)
	    {
	      // First read the node DOF values
	      {	
		MeshBase::const_node_iterator
		  it  = this->get_mesh().local_nodes_begin(),
		  end = this->get_mesh().local_nodes_end();
		
		for (; it != end; ++it)
		  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
		    {
		      assert ((*it)->dof_number(sys_num, var, comp) !=
			      DofObject::invalid_id);
		      assert (cnt < io_buffer.size());
		      this->solution->set((*it)->dof_number(sys_num, var, comp), io_buffer[cnt++]);
		    }
	      }
	      
	      // Then read the element DOF values
	      {	
		MeshBase::const_element_iterator
		  it  = this->get_mesh().local_elements_begin(),
		  end = this->get_mesh().local_elements_end();
		
		for (; it != end; ++it)
		  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
		    {
		      assert ((*it)->dof_number(sys_num, var, comp) !=
			      DofObject::invalid_id);
		      assert (cnt < io_buffer.size());
		      this->solution->set((*it)->dof_number(sys_num, var, comp), io_buffer[cnt++]);
		    }
	      }
	    }
	}
    }
}



void System::read_serialized_data (Xdr& io,
				   const bool read_additional_data)
{
  // This method implements the input of the vectors
  // contained in this System object, embedded in the 
  // output of an EquationSystems<T_sys>. 
  //
  //   10.) The global solution vector, re-ordered to be node-major 
  //       (More on this later.)                                    
  //                                                                
  //      for each additional vector in the object          
  //                                                                
  //      11.) The global additional vector, re-ordered to be       
  //          node-major (More on this later.)
  parallel_only();
  std::string comment;

  // 10.)
  // Read the global solution vector
  {
    this->read_serialized_vector(io, *this->solution); 

    // get the comment
    if (libMesh::processor_id() == 0)
      io.comment (comment);  
  }
  
  // 11.)
  // Only read additional vectors if wanted  
  if (read_additional_data)
    {	  
      std::map<std::string, NumericVector<Number>* >::const_iterator
	pos = _vectors.begin();
  
      for(; pos != this->_vectors.end(); ++pos)
        {
	  this->read_serialized_vector(io, *pos->second);

	  // get the comment
	  if (libMesh::processor_id() == 0)
	    io.comment (comment);	  
	    
	}
    }
}



template <typename iterator_type>
unsigned int System::read_serialized_blocked_dof_objects (const unsigned int var,
							  const unsigned int n_objects,
							  const iterator_type begin,
							  const iterator_type end,
							  Xdr &io,
							  NumericVector<Number> &vec) const
{
  const unsigned int sys_num = this->number();
  
  std::vector<Number> input_buffer;        // buffer to hold the input block read from io.               
  std::vector<Number> local_values;
  std::vector<std::vector<unsigned int> >  // The IDs from each processor which map to the objects
    recv_ids(libMesh::n_processors());     //  read in the current block
  std::vector<unsigned int> idx_map;       // Reordering map to traverse entry-wise rather than processor-wise
  unsigned int n_assigned_vals = 0;        // the number of values assigned, this will be returned.
  
  //-----------------------------------
  // Collect the values for all objects
  unsigned int first_object=0, last_object=0;

  for (unsigned int blk=0; last_object<n_objects; blk++)
    {
      //std::cout << "Reading object block " << blk << std::endl;
      
      // Each processor should build up its transfer buffers for its
      // local objects in [first_object,last_object).
      first_object = blk*io_blksize;
      last_object  = std::min((blk+1)*io_blksize,n_objects);
      
      // Clear the transfer buffers for this block.
      recv_ids[libMesh::processor_id()].clear();
      unsigned int n_local_dofs=0;
      for (iterator_type it=begin; it!=end; ++it)
	if (((*it)->id() >= first_object) && // object in [first_object,last_object)
	    ((*it)->id() <   last_object) &&
	    (*it)->n_comp(sys_num,var))      // var has a nonzero # of components on this object
	  {
	    recv_ids[libMesh::processor_id()].push_back((*it)->id());
	    recv_ids[libMesh::processor_id()].push_back((*it)->n_comp(sys_num, var));
	    recv_ids[libMesh::processor_id()].push_back(n_local_dofs);
	    n_local_dofs += (*it)->n_comp(sys_num, var);	    
	  }
      
      // Get the recv_ids for all other processors.
      {
	const unsigned int curr_vec_size = recv_ids[libMesh::processor_id()].size();
	std::vector<unsigned int> recv_id_sizes(libMesh::n_processors());
	Parallel::allgather(curr_vec_size, recv_id_sizes);
	for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	  {
	    recv_ids[pid].resize(recv_id_sizes[pid]);
	    Parallel::broadcast(recv_ids[pid], pid);
	  }
      }

      // create the idx map for all processors -- this will match the ordering
      // in the input buffer chunk which we are about to read.
      idx_map.resize(3*io_blksize); std::fill (idx_map.begin(), idx_map.end(), libMesh::invalid_uint);
      unsigned int tot_n_comp=0;
      for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	for (unsigned int idx=0; idx<recv_ids[pid].size(); idx+=3)
	  {
	    const unsigned int local_idx = recv_ids[pid][idx+0]-first_object;
	    assert (local_idx < std::min(io_blksize,n_objects));
	    const unsigned int n_comp    = recv_ids[pid][idx+1];
	    const unsigned int start_pos = recv_ids[pid][idx+2];
	    
	    idx_map[3*local_idx+0] = pid;
	    idx_map[3*local_idx+1] = n_comp;
	    idx_map[3*local_idx+2] = start_pos; // this tells us where the values
	                                        // for this object will live in the local_values buffer
	    tot_n_comp += n_comp;
	  }
      
      // Processor 0 will read the block from the buffer stream and send it to all other processors
      input_buffer.resize (tot_n_comp);
      if (libMesh::processor_id() == 0) io.data_stream(input_buffer.empty() ? NULL : &input_buffer[0], tot_n_comp);
      Parallel::broadcast(input_buffer);

      // Extract the values corresponding to the local objects from the input_buffer
      // and place them into the local_values temporary buffer.
      {
	unsigned int processed_size=0;
	std::vector<Number>::const_iterator next_value = input_buffer.begin();
	local_values.clear(); local_values.resize(n_local_dofs);

	for (unsigned int idx=0; idx<idx_map.size(); idx+=3)
	  if (idx_map[idx] != libMesh::invalid_uint) // this could happen when an object
	    {                                        // has no components for the current variable
	      const unsigned int pid       = idx_map[idx+0];
	      const unsigned int n_comp    = idx_map[idx+1];
	      const unsigned int start_pos = idx_map[idx+2];
	      
	      for (unsigned int comp=0; comp<n_comp; comp++)
		{
		  assert (next_value != input_buffer.end());
		  if (pid == libMesh::processor_id())
		    {
		      assert ((start_pos+comp) < local_values.size());
		      local_values[start_pos+comp] = *next_value;
		    }
		  ++next_value;
		  ++processed_size;				
		}
	    }
	assert (processed_size == input_buffer.size());
      }
      
      // A subset of the components (potentially null set) will match our objects in
      // [first_object,last_object), and we will assign the corresponding values from
      // the local_values buffer.
      for (iterator_type it=begin; it!=end; ++it)
	if (((*it)->id() >= first_object) && // object in [first_object,last_object)
	    ((*it)->id() <   last_object) &&
	    (*it)->n_comp(sys_num,var))      // var has a nonzero # of components on this object
	  {
	    const unsigned int local_idx = (*it)->id()-first_object;
	    assert (local_idx < std::min(io_blksize,n_objects));

	    const unsigned int pid       = idx_map[3*local_idx+0];
	    const unsigned int n_comp    = idx_map[3*local_idx+1];
	    const unsigned int start_pos = idx_map[3*local_idx+2];
	    
	    assert (pid == libMesh::processor_id());
	    assert (n_comp == (*it)->n_comp(sys_num, var));
	    
	    for (unsigned int comp=0; comp<n_comp; comp++)
	      {
		assert ((start_pos+comp) < local_values.size());
		const Number &value = local_values[start_pos+comp];				
		const unsigned int dof_index = (*it)->dof_number (sys_num, var, comp);
		assert (dof_index >= vec.first_local_index());
		assert (dof_index <  vec.last_local_index());
		//std::cout << "di=" << dof_index << ", val=" << value << std::endl;
		vec.set (dof_index, value);
		++n_assigned_vals;
	      }
	  }
    }
      
  return n_assigned_vals;
}



void System::read_serialized_vector (Xdr& io, NumericVector<Number>& vec)
{
  parallel_only();

#ifndef NDEBUG
  // In parallel we better be reading a parallel vector -- if not
  // we will not set all of its components below!!
  if (libMesh::n_processors() > 1)
    assert (vec.size() != vec.local_size());
   
  // If this is not the same on all processors we're in trouble!
  Parallel::verify(io_blksize);
#endif
  
  assert (io.reading());

  // vector length
  unsigned int vector_length=0, n_assigned_vals=0;

  // Get the buffer size
  if (libMesh::processor_id() == 0) io.data(vector_length, "# vector length");
  Parallel::broadcast(vector_length);
  
  // Loop over each variable in the system, and then each node/element in the mesh.
  for (unsigned int var=0; var<this->n_vars(); var++)
    {      
      //---------------------------------
      // Collect the values for all nodes
      n_assigned_vals +=
	this->read_serialized_blocked_dof_objects (var,
						   this->get_mesh().n_nodes(),
						   this->get_mesh().local_nodes_begin(),
						   this->get_mesh().local_nodes_end(),
						   io,
						   vec);
      
      
      //------------------------------------
      // Collect the values for all elements
      n_assigned_vals +=
	this->read_serialized_blocked_dof_objects (var,
						   this->get_mesh().n_elem(),
						   this->get_mesh().local_elements_begin(),
						   this->get_mesh().local_elements_end(),
						   io,
						   vec);
    } // end variable loop
  Parallel::sum (n_assigned_vals);
  assert (n_assigned_vals == vector_length);
}



void System::write_header (Xdr& io,
			   const bool write_additional_data) const
{
  /**
   * This method implements the output of a
   * System object, embedded in the output of
   * an EquationSystems<T_sys>.  This warrants some 
   * documentation.  The output of this part 
   * consists of 5 sections:
   *
   * for this system
   *
   *   5.) The number of variables in the system (unsigned int)
   *
   *   for each variable in the system
   *     
   *     6.) The name of the variable (string)
   *     
   *     7.) Combined in an FEType:
   *         - The approximation order(s) of the variable 
   *           (Order Enum, cast to int/s)
   *         - The finite element family/ies of the variable 
   *           (FEFamily Enum, cast to int/s)
   * 
   *   end variable loop
   *
   *   8.) The number of additional vectors (unsigned int),      
   *
   *     for each additional vector in the system object
   * 
   *     9.) the name of the additional vector  (string)
   *
   * end system
   */ 
  assert (io.writing());


  // Only write the header information
  // if we are processor 0.
  if (this->get_mesh().processor_id() != 0)
    return;
  
  std::string comment;
  char buf[80];

  // 5.) 
  // Write the number of variables in the system

  {
    // set up the comment
    comment = "# No. of Variables in System \"";
    comment += this->name();
    comment += "\"";    
	  
    unsigned int n_vars = this->n_vars();   
    io.data (n_vars, comment.c_str());
  }


  for (unsigned int var=0; var<this->n_vars(); var++)
    {
      // 6.)
      // Write the name of the var-th variable
      {
	// set up the comment
	comment  = "#   Name, Variable No. ";
	std::sprintf(buf, "%d", var);
	comment += buf;
	comment += ", System \"";
	comment += this->name();
	comment += "\"";
	      
	std::string var_name = this->variable_name(var);	     
	io.data (var_name, comment.c_str());
      }
      
      // 7.)
      // Write the approximation order of the var-th variable 
      // in this system
      {
	// set up the comment
	comment = "#     Approximation Order, Variable \"";
	std::sprintf(buf, "%s", this->variable_name(var).c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
	      
	int order = static_cast<int>(this->variable_type(var).order);	      
	io.data (order, comment.c_str());
      }


#ifdef ENABLE_INFINITE_ELEMENTS
      
      // do the same for radial_order
      {
	comment = "#     Radial Approximation Order, Variable \"";
	std::sprintf(buf, "%s", this->variable_name(var).c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
      
	int rad_order = static_cast<int>(this->variable_type(var).radial_order);	      
	io.data (rad_order, comment.c_str());
      }

#endif
      
      // Write the Finite Element type of the var-th variable 
      // in this System
      {
	// set up the comment
	comment = "#     FE Family, Variable \"";
	std::sprintf(buf, "%s", this->variable_name(var).c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
      
	const FEType& type = this->variable_type(var);	      
	int fam = static_cast<int>(type.family);	      
	io.data (fam, comment.c_str());

#ifdef ENABLE_INFINITE_ELEMENTS
	
	comment = "#     Radial FE Family, Variable \"";
	std::sprintf(buf, "%s", this->variable_name(var).c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
      
	int radial_fam = static_cast<int>(type.radial_family);
	io.data (radial_fam, comment.c_str());	
	
	comment = "#     Infinite Mapping Type, Variable \"";
	std::sprintf(buf, "%s", this->variable_name(var).c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";

	int i_map = static_cast<int>(type.inf_map);	      
	io.data (i_map, comment.c_str());
#endif
      }
    } // end of the variable loop
 
  // 8.) 
  // Write the number of additional vectors in the System.
  // If write_additional_data==false, then write zero for
  // the number of additional vectors.
  {	  
    {
      // set up the comment
      comment = "# No. of Additional Vectors, System \"";
      comment += this->name();
      comment += "\"";
    
      unsigned int n_vectors = write_additional_data ? this->n_vectors () : 0;
      io.data (n_vectors, comment.c_str());
    }  
      
    if (write_additional_data)
      {
	std::map<std::string, NumericVector<Number>* >::const_iterator
	  vec_pos = this->_vectors.begin();
	unsigned int cnt=0;
	
	for (; vec_pos != this->_vectors.end(); ++vec_pos)
	  {
	    // 9.)
	    // write the name of the cnt-th additional vector
	    comment =  "# Name of ";
	    std::sprintf(buf, "%d", cnt++);
	    comment += buf;
	    comment += "th vector";
	    std::string vec_name = vec_pos->first;
	    
	    io.data (vec_name, comment.c_str());
	  }
      }
  }
}



// void System::write_data (Xdr& io,
// 			 const bool write_additional_data) const
// {
//   // This is deprecated -- use write_serialized_data() instead.  This will be kept for reference
//   // for a little while, then dropped.  There is no need call this method any more.
//   deprecated();
  
//   /**
//    * This method implements the output of the vectors
//    * contained in this System object, embedded in the 
//    * output of an EquationSystems<T_sys>. 
//    *
//    *   9.) The global solution vector, re-ordered to be node-major 
//    *       (More on this later.)                                    
//    *                                                                
//    *      for each additional vector in the object          
//    *                                                                
//    *      10.) The global additional vector, re-ordered to be       
//    *           node-major (More on this later.)
//    *
//    * Note that the actual IO is handled through the Xdr class 
//    * (to be renamed later?) which provides a uniform interface to 
//    * both the XDR (eXternal Data Representation) interface and standard
//    * ASCII output.  Thus this one section of code will read XDR or ASCII
//    * files with no changes.
//    */ 
//   assert (io.writing());

//   const unsigned int proc_id = this->get_mesh().processor_id();
 
//   std::string comment;
  
//   // All processors contribute numeric vector values
//   std::vector<Number> global_vector;

//   // Collect the global solution on one processor
//   this->solution->localize_to_one (global_vector, 0);       

//   // Only processor 0 actually writes out the soltuion vector.  
//   if (proc_id == 0)
//     {	       
//       // First we need to re-order the solution so that it
//       // is dof_map agnostic.  This is necessary so that the 
//       // vector might be re-read with a different partitioning
//       // or DOF distribution.  
//       //
//       // Currently the vector is written in node-major order.
//       std::vector<Number> reordered_soln(global_vector.size());
	  
//       unsigned int cnt=0;

//       const unsigned int sys_num = this->number();
//       const unsigned int n_vars  = this->n_vars();

//       // Build a set of non subactive node indices. 
//       std::set<unsigned int> not_subactive_node_ids;
//       MeshTools::get_not_subactive_node_ids(this->get_mesh(), not_subactive_node_ids);

//       // Loop over each variable and each node, and write out the value.
//       for (unsigned int var=0; var<n_vars; var++)
//         {		
// 	  // First write the nodal DOF values
//           std::set<unsigned int>::iterator it = not_subactive_node_ids.begin();
//           const std::set<unsigned int>::iterator end = not_subactive_node_ids.end();
          
//           for (; it != end; ++it)
//           {
//             // Get the global index of this node
//             const unsigned int node = *it;
// 	    for (unsigned int index=0; index<this->get_mesh().node(node).n_comp(sys_num, var); index++)
// 	      {
//                 assert (this->get_mesh().node(node).id() == node);
// 		assert (this->get_mesh().node(node).dof_number(sys_num, var, index) !=
// 			DofObject::invalid_id);
// 		assert (cnt < reordered_soln.size());
		
// 		reordered_soln[cnt++] = 
// 		  global_vector[this->get_mesh().node(node).dof_number(sys_num, var, index)];
// 	      }
//           }

// 	  // Then write the element DOF values
// 	  {
// 	    MeshBase::const_element_iterator
// 	      it  = this->get_mesh().active_elements_begin(),
// 	      end = this->get_mesh().active_elements_end();

// 	    for (; it!=end; ++it)
// 	      for (unsigned int index=0; index<(*it)->n_comp(sys_num, var); index++)
// 	        {
// 		  assert ((*it)->dof_number(sys_num, var, index) !=
// 			  DofObject::invalid_id);
			
// 		  assert (cnt < reordered_soln.size());
			
// 		  reordered_soln[cnt++] = 
// 		      global_vector[(*it)->dof_number(sys_num, var, index)];
// 		}
// 	  }
// 	}

//       // 9.)
//       //
//       // Actually write the reordered solution vector 
//       // for the ith system to disk

//       // set up the comment
//       {
// 	comment = "# System \"";
// 	comment += this->name();
// 	comment += "\" Solution Vector";
//       }

//       io.data (reordered_soln, comment.c_str());	  
//     }
  
//   // Only write additional vectors if wanted  
//   if (write_additional_data)
//     {	  
//       std::map<std::string, NumericVector<Number>* >::const_iterator
// 	pos = _vectors.begin();
  
//       for(; pos != this->_vectors.end(); ++pos)
//         {
// 	  // fill with zero.  In general, a resize is not necessary	       
// 	  std::fill (global_vector.begin(), global_vector.end(), libMesh::zero);

// 	  // Collect the global solution on one processor
// 	  pos->second->localize_to_one (global_vector, 0);       
    
// 	  // Only processor 0 actually writes out the  vector.  
// 	  if (proc_id == 0)
// 	    {	  
// 	      // First we need to re-order the solution so that it
// 	      // is dof_map agnostic.  This is necessary so that the 
// 	      // vector might be re-read with a different partitioning
// 	      // or DOF distribution.  
// 	      //
// 	      // Currently the vector is written in node-major order.
// 	      std::vector<Number> reordered_soln(global_vector.size());
	  
// 	      unsigned int cnt=0;
	      
// 	      const unsigned int sys_num = this->number();
// 	      const unsigned int n_vars  = this->n_vars();
	      
// 	      for (unsigned int var=0; var<n_vars; var++)
// 		{		
// 		  // First write the nodal DOF values
// 		  {
// 		    MeshBase::const_node_iterator
// 		      it  = this->get_mesh().nodes_begin(),
// 		      end = this->get_mesh().nodes_end();
		    
// 		    for (; it !=end; ++it)
// 		      for (unsigned int index=0; index<(*it)->n_comp(sys_num, var); index++)
// 			{
// 			  assert ((*it)->dof_number(sys_num, var, index) !=
// 				  DofObject::invalid_id);			  
// 			  assert (cnt < reordered_soln.size());
			  
// 			  reordered_soln[cnt++] = 
// 			    global_vector[(*it)->dof_number(sys_num, var, index)];
// 			}
// 		  }
		  
// 		  // Then write the element DOF values
// 		  {
// 		    MeshBase::const_element_iterator
// 		      it  = this->get_mesh().active_elements_begin(),
// 		      end = this->get_mesh().active_elements_end();
		    
// 		    for (; it!=end; ++it)
// 		      for (unsigned int index=0; index<(*it)->n_comp(sys_num, var); index++)
// 			{
// 			  assert ((*it)->dof_number(sys_num, var, index) !=
// 				  DofObject::invalid_id);			    
// 			  assert (cnt < reordered_soln.size());
			  
// 			  reordered_soln[cnt++] = 
// 			    global_vector[(*it)->dof_number(sys_num, var, index)];
// 			}
// 		  }
// 		}

	    
// 	      // 10.)
// 	      //
// 	      // Actually write the reordered additional vector 
// 	      // for this system to disk

// 	      // set up the comment
// 	      {
// 		comment = "# System \"";
// 		comment += this->name(); 
// 		comment += "\" Additional Vector \"";
// 		comment += pos->first;
// 		comment += "\"";
// 	      }
	      
// 	      io.data (reordered_soln, comment.c_str());	  
// 	    }
// 	}
//     }
// }



void System::write_parallel_data (Xdr &io,
				  const bool write_additional_data) const
{
  /**
   * This method implements the output of the vectors
   * contained in this System object, embedded in the 
   * output of an EquationSystems<T_sys>. 
   *
   *   9.) The global solution vector, re-ordered to be node-major 
   *       (More on this later.)                                    
   *                                                                
   *      for each additional vector in the object          
   *                                                                
   *      10.) The global additional vector, re-ordered to be       
   *           node-major (More on this later.)
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */ 
  std::string comment;
  
  assert (io.writing());
  
  std::vector<Number> io_buffer; io_buffer.reserve(this->solution->local_size());
	  
  const unsigned int sys_num = this->number();
  const unsigned int n_vars  = this->n_vars();
  
  // Loop over each variable and each node, and write out the value.
  for (unsigned int var=0; var<n_vars; var++)
    {
      // First write the node DOF values
      {	
	MeshBase::const_node_iterator
	  it  = this->get_mesh().local_nodes_begin(),
	  end = this->get_mesh().local_nodes_end();
	
	for (; it != end; ++it)
	  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
	    {
	      assert ((*it)->dof_number(sys_num, var, comp) !=
		      DofObject::invalid_id);
	      
	      io_buffer.push_back((*this->solution)((*it)->dof_number(sys_num, var, comp)));
	    }
      }

      // Then write the element DOF values
      {	
	MeshBase::const_element_iterator
	  it  = this->get_mesh().local_elements_begin(),
	  end = this->get_mesh().local_elements_end();
	
	for (; it != end; ++it)
	  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
	    {
	      assert ((*it)->dof_number(sys_num, var, comp) !=
		      DofObject::invalid_id);
	      
	      io_buffer.push_back((*this->solution)((*it)->dof_number(sys_num, var, comp)));    
	    }
      }
    }
  // 9.)
  //
  // Actually write the reordered solution vector 
  // for the ith system to disk
  
  // set up the comment
  {
    comment = "# System \"";
    comment += this->name();
    comment += "\" Solution Vector";
  }
  
  io.data (io_buffer, comment.c_str());	  
  
  // Only write additional vectors if wanted  
  if (write_additional_data)
    {	  
      std::map<std::string, NumericVector<Number>* >::const_iterator
	pos = _vectors.begin();
  
      for(; pos != this->_vectors.end(); ++pos)
        {
	  io_buffer.clear(); io_buffer.reserve( pos->second->local_size());
	  
	  // Loop over each variable and each node, and write out the value.
	  for (unsigned int var=0; var<n_vars; var++)
	    {
	      // First write the node DOF values
	      {	
		MeshBase::const_node_iterator
		  it  = this->get_mesh().local_nodes_begin(),
		  end = this->get_mesh().local_nodes_end();
	
		for (; it != end; ++it)
		  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
		    {
		      assert ((*it)->dof_number(sys_num, var, comp) !=
			      DofObject::invalid_id);
		      
		      io_buffer.push_back((*pos->second)((*it)->dof_number(sys_num, var, comp)));   
		    }
	      }

	      // Then write the element DOF values
	      {	
		MeshBase::const_element_iterator
		  it  = this->get_mesh().local_elements_begin(),
		  end = this->get_mesh().local_elements_end();
	
		for (; it != end; ++it)
		  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
		    {
		      assert ((*it)->dof_number(sys_num, var, comp) !=
			      DofObject::invalid_id);
	      
		      io_buffer.push_back((*pos->second)((*it)->dof_number(sys_num, var, comp)));
		    }
	      }
	    }
	  // 10.)
	  //
	  // Actually write the reordered additional vector 
	  // for this system to disk
	  
	  // set up the comment
	  {
	    comment = "# System \"";
	    comment += this->name(); 
	    comment += "\" Additional Vector \"";
	    comment += pos->first;
	    comment += "\"";
	  }
	      
	  io.data (io_buffer, comment.c_str());	  	
	}
    }
}



void System::write_serialized_data (Xdr& io,
				    const bool write_additional_data) const
{
  /**
   * This method implements the output of the vectors
   * contained in this System object, embedded in the 
   * output of an EquationSystems<T_sys>. 
   *
   *   9.) The global solution vector, re-ordered to be node-major 
   *       (More on this later.)                                    
   *                                                                
   *      for each additional vector in the object          
   *                                                                
   *      10.) The global additional vector, re-ordered to be       
   *          node-major (More on this later.)
   */ 
  parallel_only();
  std::string comment;

  this->write_serialized_vector(io, *this->solution); 

  // set up the comment
  if (libMesh::processor_id() == 0)
    {
      comment = "# System \"";
      comment += this->name();
      comment += "\" Solution Vector";

      io.comment (comment);
    }

  // Only write additional vectors if wanted  
  if (write_additional_data)
    {	  
      std::map<std::string, NumericVector<Number>* >::const_iterator
	pos = _vectors.begin();
  
      for(; pos != this->_vectors.end(); ++pos)
        {
	  this->write_serialized_vector(io, *pos->second);

	  // set up the comment
	  if (libMesh::processor_id() == 0)
	    {
	      comment = "# System \"";
	      comment += this->name(); 
	      comment += "\" Additional Vector \"";
	      comment += pos->first;
	      comment += "\"";
	      io.comment (comment);	  
	    }
	}
    }
}



template <typename iterator_type>
unsigned int System::write_serialized_blocked_dof_objects (const NumericVector<Number> &vec,
							   const unsigned int var,
							   const unsigned int n_objects,
							   const iterator_type begin,
							   const iterator_type end,
							   Xdr &io) const
{
  
  const unsigned int sys_num = this->number();
  
  unsigned int written_length=0;                       // The numer of values written.  This will be returned 
  std::vector<unsigned int> xfer_ids;                  // The global IDs and # of components for the local objects in the current block
  std::vector<Number>       xfer_vals;                 // The raw values for the local objects in the current block
  std::vector<std::vector<unsigned int> >              // The global ID and # of components received from each processor
    recv_ids (libMesh::n_processors());                //  for the current block
  std::vector<std::vector<Number> >                    // The raw values received from each processor
    recv_vals(libMesh::n_processors());                //  for the current block
  std::vector<std::vector<Number>::iterator>           // The next value on each processor for the current block
    val_iters;
  val_iters.reserve(libMesh::n_processors());
  std::vector<unsigned int> &idx_map     = xfer_ids;   // map to traverse entry-wise rather than processor-wise (renamed for notational convenience)
  std::vector<Number>       &output_vals = xfer_vals;  // The output buffer for the current block (renamed for notational convenience)
  
  //---------------------------------
  // Collect the values for all objects
  unsigned int first_object=0, last_object=0;
  
  for (unsigned int blk=0; last_object<n_objects; blk++)
    {
      //std::cout << "Writing object block " << blk << " for var " << var << std::endl;
      
      // Each processor should build up its transfer buffers for its
      // local objects in [first_object,last_object).
      first_object = blk*io_blksize;
      last_object  = std::min((blk+1)*io_blksize,n_objects);
	  
      // Clear the transfer buffers for this block.
      xfer_ids.clear(); xfer_vals.clear();
      
      for (iterator_type it=begin; it!=end; ++it)
	if (((*it)->id() >= first_object) && // object in [first_object,last_object)
	    ((*it)->id() <   last_object) &&
	    (*it)->n_comp(sys_num,var)  )    // var has a nonzero # of components on this object
	  {
	    xfer_ids.push_back((*it)->id());
	    xfer_ids.push_back((*it)->n_comp(sys_num, var));
	    
	    for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
	      {
		assert ((*it)->dof_number(sys_num, var, comp) >= vec.first_local_index());
		assert ((*it)->dof_number(sys_num, var, comp) <  vec.last_local_index());
		xfer_vals.push_back(vec((*it)->dof_number(sys_num, var, comp)));
	      }
	  }

      //-----------------------------------------
      // Send the transfer buffers to processor 0.
      
      // Get the size of the incoming buffers -- optionally
      // we could over-size the recv buffers based on
      // some maximum size to avoid these communications
      std::vector<unsigned int> ids_size, vals_size;
      const unsigned int my_ids_size  = xfer_ids.size();
      const unsigned int my_vals_size = xfer_vals.size();
      
      Parallel::gather (0, my_ids_size,  ids_size);
      Parallel::gather (0, my_vals_size, vals_size);
      
      // Note that we will actually send/receive to ourself if we are
      // processor 0, so let's use nonblocking receives.
      std::vector<Parallel::request>
	id_request_handles(libMesh::n_processors()),
	val_request_handles(libMesh::n_processors());
      
#ifdef HAVE_MPI
      const unsigned int id_tag=0, val_tag=1;
      
      // Post the receives -- do this on processor 0 only.
      if (libMesh::processor_id() == 0)
	for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	  {
	    recv_ids[pid].resize(ids_size[pid]);
	    recv_vals[pid].resize(vals_size[pid]);
	    
	    Parallel::irecv (pid, recv_ids[pid],  id_request_handles[pid],  id_tag);
	    Parallel::irecv (pid, recv_vals[pid], val_request_handles[pid], val_tag);
	  }
      
      // Send -- do this on all processors.
      Parallel::send(0, xfer_ids,  id_tag);
      Parallel::send(0, xfer_vals, val_tag);
#else
      // On one processor there's nothing to send
      recv_ids[0] = xfer_ids;
      recv_vals[0] = xfer_vals;
#endif
      
      // -------------------------------------------------------
      // Receive the messages and write the output on processor 0.
      if (libMesh::processor_id() == 0)
	{
#ifdef HAVE_MPI
	  // Wait for all the receives to complete. We have no
	  // need for the statuses since we already know the
	  // buffer sizes.
	  MPI_Waitall (libMesh::n_processors(),
		       &id_request_handles[0],
		       MPI_STATUSES_IGNORE);
	  MPI_Waitall (libMesh::n_processors(),
		       &val_request_handles[0],
		       MPI_STATUSES_IGNORE);
#endif
	  
	  // Write the values in this block.
	  unsigned int tot_id_size=0, tot_val_size=0;
	  val_iters.clear();
	  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	    {
	      tot_id_size  += recv_ids[pid].size();		    
	      tot_val_size += recv_vals[pid].size();
	      val_iters.push_back(recv_vals[pid].begin());
	    }
	  
	  assert (tot_id_size <= 2*std::min(io_blksize,n_objects));
	  
	  // Create a map to avoid searching.  This will allow us to
	  // traverse the received values in [first_object,last_object) order.
	  idx_map.resize(3*io_blksize); std::fill (idx_map.begin(), idx_map.end(), libMesh::invalid_uint);
	  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	    for (unsigned int idx=0; idx<recv_ids[pid].size(); idx+=2)
	      {
		const unsigned int local_idx = recv_ids[pid][idx+0]-first_object;
		assert (local_idx < std::min(io_blksize,n_objects));
		const unsigned int n_comp    = recv_ids[pid][idx+1];
		
		idx_map[3*local_idx+0] = pid;
		idx_map[3*local_idx+1] = n_comp;
		idx_map[3*local_idx+2] = std::distance(recv_vals[pid].begin(), val_iters[pid]);
		val_iters[pid] += n_comp;
	      }

	  output_vals.clear(); output_vals.reserve (tot_val_size);
	  for (unsigned int idx=0; idx<idx_map.size(); idx+=3)
	    if (idx_map[idx] != libMesh::invalid_uint) // this could happen when a local object
	      {                                        // has no components for the current variable
		const unsigned int pid       = idx_map[idx+0];
		const unsigned int n_comp    = idx_map[idx+1];
		const unsigned int first_pos = idx_map[idx+2];
		
		for (unsigned int comp=0; comp<n_comp; comp++)
		  {
		    assert (first_pos + comp < recv_vals[pid].size());
		    output_vals.push_back(recv_vals[pid][first_pos + comp]);
		  }
	      }
	  assert (output_vals.size() == tot_val_size);
	  
	  // write the stream
	  io.data_stream (output_vals.empty() ? NULL : &output_vals[0], output_vals.size());
	  written_length += output_vals.size();
	} // end processor 0 conditional block	  
    } // end object block loop

  return written_length;
}



void System::write_serialized_vector (Xdr& io, const NumericVector<Number>& vec) const
{
  parallel_only();
  
  // If this is not the same on all processors we're in trouble!
  assert (Parallel::verify(io_blksize));
  assert (io.writing());
  
  unsigned int vec_length = vec.size();
  if (libMesh::processor_id() == 0) io.data (vec_length, "# vector length");
  
  unsigned int written_length = 0;
  
  // Loop over each variable in the system, and then each node/element in the mesh.
  for (unsigned int var=0; var<this->n_vars(); var++)
    {
      //---------------------------------
      // Collect the values for all nodes
      written_length +=
	this->write_serialized_blocked_dof_objects (vec,
						    var,
						    this->get_mesh().n_nodes(),
						    this->get_mesh().local_nodes_begin(),
						    this->get_mesh().local_nodes_end(),
						    io);
      
      //------------------------------------
      // Collect the values for all elements
      written_length +=
	this->write_serialized_blocked_dof_objects (vec,
						    var,
						    this->get_mesh().n_elem(),
						    this->get_mesh().local_elements_begin(),
						    this->get_mesh().local_elements_end(),
						    io);
    } // end variable loop

  if (libMesh::processor_id() == 0)
    assert(written_length == vec_length);
}


