// $Id: system_io.C,v 1.3 2004-08-05 15:58:44 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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


// C++ Includes
#include <stdio.h> // for sprintf

// Local Includes
#include "system.h"
#include "libmesh.h"
#include "mesh.h"
#include "xdr_cxx.h"

// Forward Declarations




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
   * Read the number of additional vectors
   */
  {
    unsigned int n_vectors=0;
  
    io.data (n_vectors);

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
	
    assert (global_vector.size() == this->n_dofs());
	
    unsigned int cnt=0;

    const unsigned int sys     = this->number();
    const unsigned int n_vars  = this->n_vars();
    const unsigned int n_nodes = _mesh.n_nodes();
    const unsigned int n_elem  = _mesh.n_elem();
	
    for (unsigned int var=0; var<n_vars; var++)
      {
	// First reorder the nodal DOF values
	for (unsigned int node=0; node<n_nodes; node++)
	  for (unsigned int index=0; index<_mesh.node(node).n_comp(sys,var); index++)
	    {
	      assert (_mesh.node(node).dof_number(sys, var, index) !=
		      DofObject::invalid_id);

	      assert (cnt < global_vector.size());
		  
	      reordered_vector[_mesh.node(node).dof_number(sys, var, index)] =
		  global_vector[cnt++]; 
	    }

	// Then reorder the element DOF values
	for (unsigned int elem=0; elem<n_elem; elem++)
	  for (unsigned int index=0; index<_mesh.elem(elem)->n_comp(sys,var); index++)
	    {  
	      assert (_mesh.elem(elem)->dof_number(sys, var, index) !=
		      DofObject::invalid_id);
		  
	      assert (cnt < global_vector.size());
		  
	      reordered_vector[_mesh.elem(elem)->dof_number(sys, var, index)] =
		global_vector[cnt++]; 
	    }
      }
	    
    *(this->solution) = reordered_vector;
  }





  /**
   * for each additional vector,
   * simply go through the list
   */

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
	  const unsigned int n_nodes = _mesh.n_nodes();
	  const unsigned int n_elem  = _mesh.n_elem();
	
	  for (unsigned int var=0; var<n_vars; var++)
	    {
	      // First reorder the nodal DOF values
	      for (unsigned int node=0; node<n_nodes; node++)
		for (unsigned int index=0; index<_mesh.node(node).n_comp(sys,var); index++)
		  {
		    assert (_mesh.node(node).dof_number(sys, var, index) !=
			    DofObject::invalid_id);

		    assert (cnt < global_vector.size());
		  
		    reordered_vector[_mesh.node(node).dof_number(sys, var, index)] =
			global_vector[cnt++]; 
		  }

	      // Then reorder the element DOF values
	      for (unsigned int elem=0; elem<n_elem; elem++)
		for (unsigned int index=0; index<_mesh.elem(elem)->n_comp(sys,var); index++)
		  {  
		    assert (_mesh.elem(elem)->dof_number(sys, var, index) !=
			    DofObject::invalid_id);
		  
		    assert (cnt < global_vector.size());
		  
		    reordered_vector[_mesh.elem(elem)->dof_number(sys, var, index)] =
			global_vector[cnt++]; 
		  }
	    }
	    
	  // use the overloaded operator=(std::vector) to assign the values
	  *(pos->second) = reordered_vector;
	}
    }
}









void System::write(Xdr& io) const
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
  assert (_mesh.processor_id() == 0);
  
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
	sprintf(buf, "%d", var);
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
	sprintf(buf, "%s", var_name.c_str());
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
	sprintf(buf, "%s", var_name.c_str());
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
	sprintf(buf, "%s", var_name.c_str());
	comment += buf;
	comment += "\", System \"";
	comment += this->name();
	comment += "\"";
      }
      
      FEType type = this->variable_type(var);
	      
      int fam = static_cast<int>(type.family);
	      
      io.data (fam, comment.c_str());


#ifdef ENABLE_INFINITE_ELEMENTS

      {
	comment = "# Radial FE Family, Variable \"";
	sprintf(buf, "%s", var_name.c_str());
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
	sprintf(buf, "%s", var_name.c_str());
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
   * Write the number of additional vectors in the System
   */
  {	  
    // set up the comment
    {
      comment = "# No. of Additional Vectors, System \"";
      comment += this->name();
      comment += "\"";
    }
	  
    unsigned int n_vectors = this->n_vectors ();
    io.data (n_vectors, comment.c_str());
  }



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
	sprintf(buf, "%d", cnt++);
	comment += buf;
	comment += "th vector";
	std::string vec_name = vec_pos->first;

	io.data (vec_name, comment.c_str());
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

  const unsigned int proc_id = _mesh.processor_id();
 
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
      const unsigned int n_nodes = _mesh.n_nodes();
      const unsigned int n_elem  = _mesh.n_elem();

      for (unsigned int var=0; var<n_vars; var++)
        {		
	  // First write the nodal DOF values
	  for (unsigned int node=0; node<n_nodes; node++)
	    for (unsigned int index=0; index<_mesh.node(node).n_comp(sys_num, var); index++)
	      {
		assert (_mesh.node(node).dof_number(sys_num, var, index) !=
			DofObject::invalid_id);
		      
		assert (cnt < reordered_soln.size());
		      
		reordered_soln[cnt++] = 
		    global_vector[_mesh.node(node).dof_number(sys_num, var, index)];
	      }

	  // Then write the element DOF values
	  for (unsigned int elem=0; elem<n_elem; elem++)
	    if (_mesh.elem(elem)->active())
	      for (unsigned int index=0; index<_mesh.elem(elem)->n_comp(sys_num, var); index++)
	        {
		  assert (_mesh.elem(elem)->dof_number(sys_num, var, index) !=
			  DofObject::invalid_id);
			
		  assert (cnt < reordered_soln.size());
			
		  reordered_soln[cnt++] = 
		      global_vector[_mesh.elem(elem)->dof_number(sys_num, var, index)];
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
		const unsigned int n_nodes = _mesh.n_nodes();
		const unsigned int n_elem  = _mesh.n_elem();

		for (unsigned int var=0; var<n_vars; var++)
		  {		
		    // First write the nodal DOF values
		    for (unsigned int node=0; node<n_nodes; node++)
		      for (unsigned int index=0; index<_mesh.node(node).n_comp(sys_num, var); index++)
		        {
			  assert (_mesh.node(node).dof_number(sys_num, var, index) !=
				  DofObject::invalid_id);
		      
			  assert (cnt < reordered_soln.size());
		      
			  reordered_soln[cnt++] = 
			      global_vector[_mesh.node(node).dof_number(sys_num, var, index)];
			}

		    // Then write the element DOF values
		    for (unsigned int elem=0; elem<n_elem; elem++)
		      if (_mesh.elem(elem)->active())
			for (unsigned int index=0; index<_mesh.elem(elem)->n_comp(sys_num, var); index++)
			  {
			    assert (_mesh.elem(elem)->dof_number(sys_num, var, index) !=
				    DofObject::invalid_id);
			
			    assert (cnt < reordered_soln.size());
			
			    reordered_soln[cnt++] = 
				global_vector[_mesh.elem(elem)->dof_number(sys_num, var, index)];
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


