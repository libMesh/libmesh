// $Id: equation_systems_io.C,v 1.11 2003-02-12 02:03:49 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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


#include "mesh_common.h"


// System Includes
#ifdef HAVE_RPC_RPC_H
# include <rpc/rpc.h>
#endif



// Local Includes
#include "fe_type.h"
#include "petsc_interface.h"
#include "equation_systems.h"
#include "general_system.h"

// Forward Declarations




// ------------------------------------------------------------
// EquationSystem class implementation
void EquationSystems::read(const std::string& name,
			   const Xdr::XdrMODE mode,
			   const bool read_header,
			   const bool read_data)
{
  /**
   * This program implements the output of an 
   * EquationSystems object.  This warrants some 
   * documentation.  The output file essentially
   * consists of 6 sections:
   *
   * 1.) The number of individual equation systems (unsigned int)
   * 
   * for each system
   *
   *   2.) The name of the system (string) 
   *   3.) The number of variables in the system (unsigned int)
   *
   *   for each variable in the system
   *     
   *     4.) The name of the variable (string)
   *     5. & 6.) Combined in an FEType:
   *              - The approximation order of the variable (Order Enum, cast to int)
   *              - The finite element family/ies of the variable (FEType struct, cast to int/s)
   *
   *   end variable loop
   * end system loop
   *
   * for each system
   *   
   *   7.) The global solution vector, re-ordered to be node-major
   *       (More on this later.)
   *
   * end system loop
   *
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */

  
  Xdr io (name, mode);

  assert (io.reading());


  /**
   * Possibly clear data structures and start from scratch.
   */
  if (read_header)
    clear ();
      
  /**
   * 1.)  
   *
   * Read the number of equation systms
   */
  unsigned int n_sys=0;
  
  io.data (n_sys);
  
  for (unsigned int sys=0; sys<n_sys; sys++)
    {
      /**
       * 2.)
       *
       * Read the name of the ith system
       */
      std::string sys_name;
      
      io.data (sys_name);
      
      if (read_header) add_system (sys_name);
	  
      GeneralSystem& new_system = (*this)(sys_name);
	  
      /**
       * 3.) 
       *
       * Read the number of variables in the ith system
       */
      unsigned int n_vars=0;
      
      io.data (n_vars);
      
      for (unsigned int var=0; var<n_vars; var++)
	{	            
	  /**
	   * 4.)
	   *
	   * Read the name of the jth variable in the ith system
	   */
	  std::string var_name;
	  
	  io.data (var_name);
	      
	  /**
	   * 5.)
	   *
	   * Read the approximation order(s) of the jth variable 
	   * in the ith system
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
	   * 6.)
	   *
	   * Read the finite element type of the jth variable 
	   * in the ith system
	   */
	  int fam=0;
	  
	  io.data (fam);

	  FEType type;

#ifndef ENABLE_INFINITE_ELEMENTS

	  type.order  = static_cast<Order>(order);
	  type.family = static_cast<FEFamily>(fam);

#else

//TODO:[DD] flag the use of infinite elements, somewhere in the outfile?
	  int radial_fam=0;
	  int i_map=0;
	  
	  io.data (radial_fam);
	  io.data (i_map);

	  type.order         = static_cast<Order>(order);
	  type.radial_order  = static_cast<Order>(rad_order);
	  type.family        = static_cast<FEFamily>(fam);
	  type.radial_family = static_cast<FEFamily>(radial_fam);
	  type.inf_map       = static_cast<InfMapType>(i_map);	  
#endif


	  if (read_header) new_system.add_variable (var_name,
						    type);



	};
    };
      
  /**
   * Now we are ready to initialize the underlying data
   * structures. This will initialize the vectors for 
   * storage, the dof_map, etc...
   */ 
  if (read_header) init();

  

  /**
   * 7.)
   *
   * Read and set the numeric vector values
   */
  if (read_data)
    for (unsigned int sys=0; sys<n_systems(); sys++)
      {
	GeneralSystem& system = (*this)(sys);
	DofMap&     dof_map   = system.get_dof_map(); 
	std::vector<Complex> global_soln;
	std::vector<Complex> reordered_soln;
	
	io.data (global_soln);	  
	
	/**
	 * Remember that the stored vector is node-major.
	 * We need to put it into whatever application-specific
	 * ordering we may have using the dof_map.
	 */
	reordered_soln.resize(global_soln.size());
	
	assert (global_soln.size() == system.n_dofs());
	
	unsigned int cnt=0;

	const unsigned int n_vars  = system.n_vars();
	const unsigned int n_nodes = _mesh.n_nodes();
	const unsigned int n_elem  = _mesh.n_elem();
	
	for (unsigned int var=0; var<n_vars; var++)
	  {
	    // First reorder the nodal DOF values
	    for (unsigned int node=0; node<n_nodes; node++)
	      for (unsigned int index=0;
		   index<dof_map.n_dofs_at_node(node, var); index++)
		{
		  
		  assert (dof_map.node_dof_number(node, var, index) !=
			  dof_map.invalid_number);

		  assert (cnt < global_soln.size());
		  
		  reordered_soln[dof_map.node_dof_number(node, var, index)] =
		    global_soln[cnt++]; 
		};

	    // Then reorder the element DOF values
	    for (unsigned int elem=0; elem<n_elem; elem++)
	      for (unsigned int index=0;
		   index<dof_map.n_dofs_on_elem(elem, var); index++)
		{
		  
		  assert (dof_map.elem_dof_number(elem, var, index) !=
			  dof_map.invalid_number);

		  assert (cnt < global_soln.size());
		  
		  reordered_soln[dof_map.elem_dof_number(elem, var, index)] =
		    global_soln[cnt++]; 
		};
	  };
	    
	*(system.solution) = reordered_soln;
      };
};



void EquationSystems::write(const std::string& name,
			    const Xdr::XdrMODE mode,
			    const bool write_data)
{
  /**
   * This program implements the output of an 
   * EquationSystems object.  This warrants some 
   * documentation.  The output file essentially
   * consists of 6 sections:
   *
   * 1.) The number of individual equation systems (unsigned int)
   * 
   * for each system
   *
   *   2.) The name of the system (string) 
   *   3.) The number of variables in the system (unsigned int)
   *
   *   for each variable in the system
   *     
   *     4.) The name of the variable (string)
   *     5. & 6.) Combined in an FEType:
   *              - The approximation order(s) of the variable (Order Enum, cast to int/s)
   *              - The finite element family/ies of the variable (FEType struct, cast to int/s)
   *
   *   end variable loop
   * end system loop
   *
   * for each system
   *   
   *   7.) The global solution vector, re-ordered to be node-major
   *       (More on this later.)
   *
   * end system loop
   *
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */

  Xdr io(name, mode);

  assert (io.writing());

  const unsigned int proc_id = _mesh.processor_id();
  unsigned int n_sys         = n_systems();

  std::map<std::string, GeneralSystem*>::iterator
    pos = _systems.begin();
  
  unsigned int sys_num=0;
  std::string comment;
  char buf[80];



  /**
   * Only write the header information
   * if we are processor 0.
   */
  if (proc_id == 0) 
    {
      /**
       * 1.)  
       *
       * Write the number of equation systms
       */
      io.data (n_sys, "# The number of equation systems");
        
      while (pos != _systems.end())
	{
	  std::string sys_name  = pos->first;
	  GeneralSystem& system = *pos->second;
	  


	  /**
	   * 2.)
	   *
	   * Write the name of the ith system
	   */

	  // set up the comment
	  {
	    comment =  "# System ";
	    sprintf(buf, "%d", sys_num);
	    comment += buf;
	    comment += " name";
	  }
	  
	  io.data (sys_name, comment.c_str());


	  
	  /**
	   * 3.) 
	   *
	   * Write the number of variables in the ith system
	   */
	  
	  // set up the comment
	  {
	    comment = "# No. of variables in system \"";
	    sprintf(buf, "%s\"", sys_name.c_str());
	    comment += buf;
	  }
	  
	  unsigned int n_vars = system.n_vars();

	  io.data (n_vars, comment.c_str());


	  
	  for (unsigned int var=0; var<n_vars; var++)
	    {
	      /**
	       * 4.)
	       *
	       * Write the name of the jth variable in the ith system
	       */

	      // set up the comment
	      {
		comment  = "# Variable No. ";
		sprintf(buf, "%d", var);
		comment += buf;
		comment += " name, system \"";
		sprintf(buf, "%s\"", sys_name.c_str());
		comment += buf;
	      }
	      
	      std::string var_name = system.variable_name(var);
	     
	      io.data (var_name, comment.c_str());
	      
	      


	      /**
	       * 5.)
	       *
	       * Write the approximation order of the jth variable 
	       * in the ith system
	       */

	      // set up the comment
	      {
		comment = "# Variable \"";
		sprintf(buf, "%s", var_name.c_str());
		comment += buf;
		comment += "\", system \"";
		sprintf(buf, "%s\"", sys_name.c_str());
		comment += buf;
		comment += ", approximation order";
	      }
	      
	      int order = static_cast<int>(system.variable_type(var).order);
	      
	      io.data (order, comment.c_str());
	   

#ifdef ENABLE_INFINITE_ELEMENTS
	      /**
	       * do the same for radial_order
	       */
	      int rad_order = static_cast<int>(system.variable_type(var).radial_order);
	      
	      io.data (rad_order);

#endif
   


	      /**
	       * 6.)
	       *
	       * Write the Finite Element type of the jth variable 
	       * in the ith system
	       */

	      // set up the comment
	      {
		comment = "# Variable \"";
		sprintf(buf, "%s", var_name.c_str());
		comment += buf;
		comment += "\", system \"";
		sprintf(buf, "%s\"", sys_name.c_str());
		comment += buf;
		comment += ", finite element type";
	      }

	      FEType type = system.variable_type(var);
	      
	      int fam = static_cast<int>(type.family);
	      
	      io.data (fam, comment.c_str());


#ifdef ENABLE_INFINITE_ELEMENTS

	      int radial_fam = static_cast<int>(type.radial_family);
	      int i_map = static_cast<int>(type.inf_map);
	      
	      io.data (radial_fam);
	      io.data (i_map);

#endif


	    };

	  ++pos;
	  ++sys_num;
	};      
    };

  pos     = _systems.begin();
  sys_num = 0;

  /**
   * All processors contribute numeric vector values
   */
  if (write_data)
    while (pos != _systems.end())
      {
	// Convenient references
	std::string sys_name  = pos->first;
	GeneralSystem& system = *pos->second;
	DofMap&     dof_map   = system.get_dof_map();

	std::vector<Complex> global_soln;
	
	/**
	 * Collect the global solution on one processor
	 */
	system.solution->localize_to_one (global_soln, 0);       
      

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
	    std::vector<Complex> reordered_soln(global_soln.size());
	  
	    unsigned int cnt=0;

	    const unsigned int n_vars  = system.n_vars();
	    const unsigned int n_nodes = _mesh.n_nodes();
	    const unsigned int n_elem  = _mesh.n_elem();

	    for (unsigned int var=0; var<n_vars; var++)
	      {		
		// First write the nodal DOF values
		for (unsigned int node=0; node<n_nodes; node++)
		  for (unsigned int index=0;
		       index<dof_map.n_dofs_at_node(node, var); index++)
		    {
		      assert (dof_map.node_dof_number(node, var, index) !=
			      dof_map.invalid_number);
		      
		      assert (cnt < reordered_soln.size());
		      
		      reordered_soln[cnt++] = 
			global_soln[dof_map.node_dof_number(node, var, index)];
		    };

		// Then write the element DOF values
		for (unsigned int elem=0; elem<n_elem; elem++)
		  if (_mesh.elem(elem)->active())
		    for (unsigned int index=0;
			 index<dof_map.n_dofs_on_elem(elem, var); index++)
		      {
			assert (dof_map.elem_dof_number(elem, var, index) !=
				dof_map.invalid_number);
			
			assert (cnt < reordered_soln.size());
			
			reordered_soln[cnt++] = 
			  global_soln[dof_map.elem_dof_number(elem, var, index)];
		      };
	      };
	    
	    /**
	     * 7.)
	     *
	     * Actually write the reordered solution vector 
	     * for the ith system to disk
	     */

	    // set up the comment
	    {
	      comment = "# System \"";
	      sprintf(buf, "%s\"", sys_name.c_str());
	      comment += buf;
	      comment += " solution vector";
	    }

	    io.data (reordered_soln, comment.c_str());	  
	  };

	++pos;
	++sys_num;
      };
};
