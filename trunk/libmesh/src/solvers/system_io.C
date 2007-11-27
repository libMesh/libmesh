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
#warning "Increase to something usable for final checkin"
  static const unsigned int io_blksize = 100;
}


// ------------------------------------------------------------
// System class implementation
void System::read (Xdr& io,
		   const bool read_header,
		   const bool read_additional_data)
{
  /**
   * This program implements the output of a
   * System object, embedded in the output of
   * an EquationSystems<T_sys>.  This warrants some 
   * documentation.  The output file essentially
   * consists of 4 sections:
   *
   * for this system
   *
   *   1.) The number of variables in the system (unsigned int)
   *
   *   for each variable in the system
   *     
   *     2.) The name of the variable (string)
   *     
   *     3.) Combined in an FEType:
   *         - The approximation order(s) of the variable 
   *           (Order Enum, cast to int/s)
   *         - The finite element family/ies of the variable 
   *           (FEFamily Enum, cast to int/s)
   * 
   *   end variable loop
   *
   *   4.) The number of additional vectors (unsigned int),      
   *
   *     for each additional vector in the system object
   * 
   *     5.) the name of the additional vector  (string)
   *
   * end system
   *
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */
  assert (io.reading());

  /**
   * Possibly clear data structures and start from scratch.
   */
  if (read_header)
    this->clear ();


  /**
   * 1.) 
   *
   * Read the number of variables in the system
   */
  {
    unsigned int n_vars=0;
      
    io.data (n_vars);
      
    for (unsigned int var=0; var<n_vars; var++)
      {	            
	/**
	 * 2.)
	 *
	 * Read the name of the var-th variable
	 */
	std::string var_name;
	  
	io.data (var_name);
	      
	/**
	 * 3.)
	 *
	 * Read the approximation order(s) of the var-th variable 
	 */
	int order=0;
	  
	io.data (order);

#ifdef ENABLE_INFINITE_ELEMENTS

	/**
	 * do the same for radial_order
	 */
	int rad_order=0;
	  
	io.data(rad_order);

#endif

	/**
	 * Read the finite element type of the var-th variable 
	 */
	int fam=0;
	  
	io.data (fam);

	FEType type;

	type.order  = static_cast<Order>(order);
	type.family = static_cast<FEFamily>(fam);
	
#ifdef ENABLE_INFINITE_ELEMENTS

	/**
	 * Read additional information for infinite elements
	 */
	int radial_fam=0;
	int i_map=0;
	  
	io.data (radial_fam);
	io.data (i_map);

	type.radial_order  = static_cast<Order>(rad_order);
	type.radial_family = static_cast<FEFamily>(radial_fam);
	type.inf_map       = static_cast<InfMapType>(i_map);	  

#endif

	if (read_header) 
	  this->add_variable (var_name,
			      type);
      }
  }



  /**
   * 4.)  
   *
   * Read the number of additional vectors.  
   */
  unsigned int n_vectors=0;
  
  io.data (n_vectors);

  /**
   * If n_vectors > 0, this means that write_additional_data
   * was true when this file was written.  We will need to
   * make use of this fact later.
   */
  if (n_vectors > 0)
    this->_additional_data_written = true;  
  
  for (unsigned int vec=0; vec<n_vectors; vec++)
    {
      /**
       * 5.)
       *
       * Read the name of the vec-th additional vector
       */
      std::string vec_name;
      
      io.data (vec_name);

	
      if (read_additional_data)
	{
	  // sanity checks
	  assert(this->_can_add_vectors);
	  assert(this->_vectors.count(vec_name) == 0);

	  this->add_vector(vec_name);
	}
    }
}



void System::read_data (Xdr& io,
			const bool read_additional_data)
{
  /**
   * This program implements the output of the vectors
   * contained in this System object, embedded in the 
   * output of an EquationSystems<T_sys>. 
   *
   *   14.) The global solution vector, re-ordered to be node-major 
   *       (More on this later.)                                    
   *                                                                
   *      for each additional vector in the object          
   *                                                                
   *      15.) The global additional vector, re-ordered to be       
   *           node-major (More on this later.)
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */

  assert (io.reading());


  /**
   * read and reordering buffers
   */
  std::vector<Number> global_vector;
  std::vector<Number> reordered_vector;



  /**
   * 14.)
   *
   * Read and set the solution vector
   */
  {
	
    io.data (global_vector);	  
	
    /**
     * Remember that the stored vector is node-major.
     * We need to put it into whatever application-specific
     * ordering we may have using the dof_map.
     */
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





  /**
   * For each additional vector, simply go through the list.
   * ONLY attempt to do this IF additional data was actually
   * written to the file for this system (controlled by the
   * _additional_data_written flag).  
   */
  if (this->_additional_data_written)
    {
      std::map<std::string, NumericVector<Number>* >::iterator
	pos = this->_vectors.begin();
  
      for (; pos != this->_vectors.end(); ++pos)
	{
	  /**
	   * 15.)
	   *
	   * Read the values of the vec-th additional vector.
	   * Prior do _not_ clear, but fill with zero, since the
	   * additional vectors _have_ to have the same size
	   * as the solution vector
	   */
	  std::fill (global_vector.begin(), global_vector.end(), libMesh::zero);


	  io.data (global_vector);	  


	  // If read_additional_data==true, then we will keep this vector, otherwise
	  // we are going to throw it away.
	  if (read_additional_data)
	    {
	      /**
	       * Remember that the stored vector is node-major.
	       * We need to put it into whatever application-specific
	       * ordering we may have using the dof_map.
	       */
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



void System::read_parallel_data (Xdr& io,
				 const bool read_additional_data)
{
  parallel_only();
  error();
}



void System::write(Xdr& io,
		   const bool write_additional_data) const
{
  /**
   * This program implements the output of a
   * System object, embedded in the output of
   * an EquationSystems<T_sys>.  This warrants some 
   * documentation.  The output of this part 
   * consists of 5 sections:
   *
   * for this system
   *
   *   1.) The number of variables in the system (unsigned int)
   *
   *   for each variable in the system
   *     
   *     2.) The name of the variable (string)
   *     
   *     3.) Combined in an FEType:
   *         - The approximation order(s) of the variable 
   *           (Order Enum, cast to int/s)
   *         - The finite element family/ies of the variable 
   *           (FEFamily Enum, cast to int/s)
   * 
   *   end variable loop
   *
   *   4.) The number of additional vectors (unsigned int),      
   *
   *     for each additional vector in the system object
   * 
   *     5.) the name of the additional vector  (string)
   *
   * end system
   *
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will write XDR or ASCII
   * files with no changes.
   */
 
  assert (io.writing());


  /**
   * Only write the header information
   * if we are processor 0.
   */
  assert (this->get_mesh().processor_id() == 0);
  
  std::string comment;
  char buf[80];



  /**
   * 1.) 
   *
   * Write the number of variables in the system
   */
	  
  // set up the comment
  {
    comment = "# No. of Variables in System \"";
    comment += this->name();
    comment += "\"";
  }
	  
  unsigned int n_vars = this->n_vars();
    
  io.data (n_vars, comment.c_str());



  for (unsigned int var=0; var<n_vars; var++)
    {
      /**
       * 2.)
       *
       * Write the name of the var-th variable
       */

	// set up the comment
      {
	comment  = "# Name, Variable No. ";
	std::sprintf(buf, "%d", var);
	comment += buf;
	comment += ", System \"";
	comment += this->name();
	comment += "\"";
      }
	      
      std::string var_name = this->variable_name(var);
	     
      io.data (var_name, comment.c_str());
	      
	      
      /**
       * 3.)
       *
       * Write the approximation order of the var-th variable 
       * in this system
       */

      // set up the comment
      {
	comment = "# Approximation Order, Variable \"";
	std::sprintf(buf, "%s", var_name.c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
      }
	      
      int order = static_cast<int>(this->variable_type(var).order);
	      
      io.data (order, comment.c_str());


#ifdef ENABLE_INFINITE_ELEMENTS

      /**
       * do the same for radial_order
       */
      {
	comment = "# Radial Approximation Order, Variable \"";
	std::sprintf(buf, "%s", var_name.c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
      }
      int rad_order = static_cast<int>(this->variable_type(var).radial_order);
	      
      io.data (rad_order, comment.c_str());

#endif
   


      /**
       *
       * Write the Finite Element type of the var-th variable 
       * in this System
       */

      // set up the comment
      {
	comment = "# FE Family, Variable \"";
	std::sprintf(buf, "%s", var_name.c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
      }
      
      const FEType& type = this->variable_type(var);
	      
      int fam = static_cast<int>(type.family);
	      
      io.data (fam, comment.c_str());


#ifdef ENABLE_INFINITE_ELEMENTS

      {
	comment = "# Radial FE Family, Variable \"";
	std::sprintf(buf, "%s", var_name.c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
      }

      int radial_fam = static_cast<int>(type.radial_family);
      int i_map = static_cast<int>(type.inf_map);
	      
      io.data (radial_fam, comment.c_str());

      {
	comment = "# Infinite Mapping Type, Variable \"";
	std::sprintf(buf, "%s", var_name.c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
      }

      io.data (i_map, comment.c_str());

#endif


    }



  /**
   * 4.) 
   *
   * Write the number of additional vectors in the System.
   * If write_additional_data==false, then write zero for
   * the number of additional vectors.
   */
  {	  
    // set up the comment
    {
      comment = "# No. of Additional Vectors, System \"";
      comment += this->name();
      comment += "\"";
    }
    
    unsigned int n_vectors = this->n_vectors ();

    if (write_additional_data == false)
      n_vectors = 0;
    
    io.data (n_vectors, comment.c_str());
  }
  
      
  if (write_additional_data)
    {
      {
	std::map<std::string, NumericVector<Number>* >::const_iterator
	  vec_pos = this->_vectors.begin();
	
	unsigned int cnt=0;
	
	for(;  vec_pos != this->_vectors.end(); ++vec_pos)
	  {
	    /**
	     * 5.)
	     *
	     * write the name of the cnt-th additional vector
	     */
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



void System::write_data (Xdr& io,
			 const bool write_additional_data) const
{
  /**
   * This program implements the output of the vectors
   * contained in this System object, embedded in the 
   * output of an EquationSystems<T_sys>. 
   *
   *   14.) The global solution vector, re-ordered to be node-major 
   *       (More on this later.)                                    
   *                                                                
   *      for each additional vector in the object          
   *                                                                
   *      15.) The global additional vector, re-ordered to be       
   *           node-major (More on this later.)
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */
 
  assert (io.writing());

  const unsigned int proc_id = this->get_mesh().processor_id();
 
  std::string comment;


  /**
   * All processors contribute numeric vector values
   */
  std::vector<Number> global_vector;

	
  /**
   * Collect the global solution on one processor
   */
  this->solution->localize_to_one (global_vector, 0);       

    

  /**
   * Only processor 0 actually writes out the soltuion
   * vector.  
   */
  if (proc_id == 0)
    {	  
      /**
       * First we need to re-order the solution so that it
       * is dof_map agnostic.  This is necessary so that the 
       * vector might be re-read with a different partitioning
       * or DOF distribution.  
       *
       * Currently the vector is written in node-major order.
       * Obviously, a value is only written out if it corresponds
       * to a global DOF.  The code should make this clear.
       */
      std::vector<Number> reordered_soln(global_vector.size());
	  
      unsigned int cnt=0;

      const unsigned int sys_num = this->number();
      const unsigned int n_vars  = this->n_vars();

      // Build a set of non subactive node indices. 
      std::set<unsigned int> not_subactive_node_ids;
      MeshTools::get_not_subactive_node_ids(this->get_mesh(), not_subactive_node_ids);

      // Loop over each variable and each node, and write out the value.
      for (unsigned int var=0; var<n_vars; var++)
        {		
	  // First write the nodal DOF values
          std::set<unsigned int>::iterator it = not_subactive_node_ids.begin();
          const std::set<unsigned int>::iterator end = not_subactive_node_ids.end();
          
          for (; it != end; ++it)
          {
            // Get the global index of this node
            const unsigned int node = *it;
            // std::cout << "Setting value for node " << node << std::endl;
	    for (unsigned int index=0; index<this->get_mesh().node(node).n_comp(sys_num, var); index++)
            {
                assert (this->get_mesh().node(node).id() == node);

		assert (this->get_mesh().node(node).dof_number(sys_num, var, index) !=
			DofObject::invalid_id);
		      
		assert (cnt < reordered_soln.size());
		    
		reordered_soln[cnt++] = 
		    global_vector[this->get_mesh().node(node).dof_number(sys_num, var, index)];
            }
          }

	  // Then write the element DOF values
	  {
	    MeshBase::const_element_iterator
	      it  = this->get_mesh().active_elements_begin(),
	      end = this->get_mesh().active_elements_end();

	    for (; it!=end; ++it)
	      for (unsigned int index=0; index<(*it)->n_comp(sys_num, var); index++)
	        {
		  assert ((*it)->dof_number(sys_num, var, index) !=
			  DofObject::invalid_id);
			
		  assert (cnt < reordered_soln.size());
			
		  reordered_soln[cnt++] = 
		      global_vector[(*it)->dof_number(sys_num, var, index)];
		}
	  }
	}

      /**
       * 14.)
       *
       * Actually write the reordered solution vector 
       * for the ith system to disk
       */

      // set up the comment
      {
	comment = "# System \"";
	comment += this->name();
	comment += "\" Solution Vector";
      }

      io.data (reordered_soln, comment.c_str());	  
    }



  /**
   * Only write additional vectors if wanted  
   */
  if (write_additional_data)
    {	  

      std::map<std::string, NumericVector<Number>* >::const_iterator
	pos = _vectors.begin();
  
      for(; pos != this->_vectors.end(); ++pos)
        {
	  /**
	   * fill with zero.  In general, a resize is not necessary
	   */	    
	  std::fill (global_vector.begin(), global_vector.end(), libMesh::zero);

	    
	  /**
	   * Collect the global solution on one processor
	   */
	  pos->second->localize_to_one (global_vector, 0);       
    

	  /**
	   * Only processor 0 actually writes out the soltuion
	   * vector.  
	   */
	  if (proc_id == 0)
	    {	  
	      /**
	       * First we need to re-order the solution so that it
	       * is dof_map agnostic.  This is necessary so that the 
	       * vector might be re-read with a different partitioning
	       * or DOF distribution.  
	       *
	       * Currently the vector is written in node-major order.
	       * Obviously, a value is only written out if it corresponds
	       * to a global DOF.  The code should make this clear.
	       */
		std::vector<Number> reordered_soln(global_vector.size());
	  
		unsigned int cnt=0;

		const unsigned int sys_num = this->number();
		const unsigned int n_vars  = this->n_vars();

		for (unsigned int var=0; var<n_vars; var++)
		  {		
		    // First write the nodal DOF values
		    {
		      MeshBase::const_node_iterator
			it  = this->get_mesh().nodes_begin(),
			end = this->get_mesh().nodes_end();

		      for (; it !=end; ++it)
			for (unsigned int index=0; index<(*it)->n_comp(sys_num, var); index++)
			  {
			    assert ((*it)->dof_number(sys_num, var, index) !=
				    DofObject::invalid_id);
			    
			    assert (cnt < reordered_soln.size());
			    
			    reordered_soln[cnt++] = 
			      global_vector[(*it)->dof_number(sys_num, var, index)];
			  }
		    }
		    
		    // Then write the element DOF values
		    {
		      MeshBase::const_element_iterator
			it  = this->get_mesh().active_elements_begin(),
			end = this->get_mesh().active_elements_end();
		      
		      for (; it!=end; ++it)
			for (unsigned int index=0; index<(*it)->n_comp(sys_num, var); index++)
			  {
			    assert ((*it)->dof_number(sys_num, var, index) !=
				    DofObject::invalid_id);
			    
			    assert (cnt < reordered_soln.size());
			    
			    reordered_soln[cnt++] = 
			      global_vector[(*it)->dof_number(sys_num, var, index)];
			  }
		    }
		  }

	    
		/**
		 * 15.)
		 *
		 * Actually write the reordered additional vector 
		 * for this system to disk
		 */

		// set up the comment
		{
		  comment = "# System \"";
		  comment += this->name(); 
		  comment += "\" Additional Vector \"";
		  comment += pos->first;
		  comment += "\"";
		}

		io.data (reordered_soln, comment.c_str());	  
	    }
	}
    }
}



void System::write_parallel_data (Xdr& io,
				  const bool write_additional_data) const
{
  parallel_only();
  
#ifndef NDEBUG
  // If this is not the same on all processors we're in trouble!
  Parallel::verify(io_blksize);
#endif
  
  /**
   * This program implements the output of the vectors
   * contained in this System object, embedded in the 
   * output of an EquationSystems<T_sys>. 
   *
   *   1.) The global solution vector, re-ordered to be node-major 
   *       (More on this later.)                                    
   *                                                                
   *      for each additional vector in the object          
   *                                                                
   *      2.) The global additional vector, re-ordered to be       
   *          node-major (More on this later.)
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will write XDR or ASCII
   * files with no changes.
   */ 
  assert (io.writing());

  std::string comment;

  const unsigned int sys_num = this->number();
  const unsigned int n_vars  = this->n_vars();

  std::vector<unsigned int> xfer_ids;
  std::vector<Number>       xfer_vals;

  std::vector<std::vector<unsigned int> > recv_ids (libMesh::n_processors());
  std::vector<std::vector<Number> >       recv_vals(libMesh::n_processors());
  std::vector<std::vector<Number>::const_iterator>  val_iters; val_iters.reserve(libMesh::n_processors());

  // Loop over each variable in the system, and then each node/element in the mesh.
  for (unsigned int var=0; var<n_vars; var++)
    {
      const unsigned int n_nodes = this->get_mesh().n_nodes();      
      unsigned int first_node=0, last_node=0;

      //---------------------------------
      // Collect the values for all nodes
      for (unsigned int blk=0; last_node<n_nodes; blk++)
	{
	  std::cout << "Writing node block " << blk << std::endl;
	  
	  // Each processor should build up its transfer buffers for its
	  // local nodes in [first_node,last_node).
	  first_node = blk*io_blksize;
	  last_node  = std::min((blk+1)*io_blksize,n_nodes);
	  
	  // Clear the transfer buffers for this block.
	  xfer_ids.clear(); xfer_vals.clear();

	  MeshBase::const_node_iterator
	    it  = this->get_mesh().local_nodes_begin(),
	    end = this->get_mesh().local_nodes_end();

	  for (; it!=end; ++it)
	    if (((*it)->id() >= first_node) && /* node in [first_node,last_node) */
		((*it)->id() <   last_node) &&
		(*it)->n_comp(sys_num,var))    /* var has n_comp components on this node */
	      {
		xfer_ids.push_back((*it)->id());
		xfer_ids.push_back((*it)->n_comp(sys_num, var));

		for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
		  {
		    assert ((*it)->dof_number(sys_num, var, comp) >= this->solution->first_local_index());
		    assert ((*it)->dof_number(sys_num, var, comp) <  this->solution->last_local_index());
		    xfer_vals.push_back((*this->solution)((*it)->dof_number(sys_num, var, comp)));
		  }
	      }

	  //-----------------------------------------
	  // Send the transfer buffers to processor 0.
	  
	  // Get the size of the incoming buffers -- optionally
	  // we could over-size the recv buffers based on
	  // some maximum size to avoid these communications
	  std::vector<unsigned int> ids_size, vals_size;
	  Parallel::gather (0, static_cast<unsigned int>(xfer_ids.size()),  ids_size);
	  Parallel::gather (0, static_cast<unsigned int>(xfer_vals.size()), vals_size);
	  
	  // Note that we will actually send/receive to ourself if we are
	  // processor 0, so let's use nonblocking receives.
	  std::vector<MPI_Request>
	    id_request_handles(libMesh::n_processors()),
	    val_request_handles(libMesh::n_processors());
	    
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
	    
	  // Wait for all the receives to complete -- do this on processor 0 only.
	  // We have no need for the statuses since we already know the buffer sizes.
	  if (libMesh::processor_id() == 0)
	    {
	      MPI_Waitall (libMesh::n_processors(), &id_request_handles[0],  MPI_STATUSES_IGNORE);
	      MPI_Waitall (libMesh::n_processors(), &val_request_handles[0], MPI_STATUSES_IGNORE);
	      
	      // Write the values in this block.
	      unsigned int tot_id_size = 0;
	      val_iters.clear();
	      for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
		{
		  tot_id_size += recv_ids[pid].size();
		  
		  val_iters.push_back(recv_vals[pid].begin());
		}

	      assert (tot_id_size <= 2*io_blksize);

	      // Create a useful map to avoid searching 
	      std::vector<unsigned int> idx_map(tot_id_size);
	      for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
		for (unsigned int idx=0; idx<recv_ids[pid].size(); idx+=2)
		  {
		    const unsigned int local_idx = recv_ids[pid][idx]-first_node;
		    assert (local_idx < std::min(io_blksize,n_nodes));
		    const unsigned int n_comp    = recv_ids[pid][idx+1];

		    idx_map[2*local_idx+0] = pid;
		    idx_map[2*local_idx+1] = n_comp;
		  }

	      for (unsigned int idx=0; idx<idx_map.size(); idx+=2)
		{
		  const unsigned int pid    = idx_map[idx+0];
		  const unsigned int n_comp = idx_map[idx+1];

		  for (unsigned int comp=0; comp<n_comp; comp++)
		    {
		      assert (val_iters[pid] != recv_vals[pid].end());
		      Number value = *val_iters[pid];
		      io.data_stream (&value, 1);
		      ++val_iters[pid];
		    }
		}
	    }
	  
	}            
    }
  
  
}


