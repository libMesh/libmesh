// $Id: equation_systems_io.C,v 1.23 2003-04-09 01:20:21 benkirk Exp $

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
#include "frequency_system.h"
#include "thin_system.h"

// Forward Declarations




// ------------------------------------------------------------
// EquationSystem class implementation
template <typename T_sys>
void EquationSystems<T_sys>::read (const std::string& name,
				   const Xdr::XdrMODE mode,
				   const bool read_header,
				   const bool read_data,
				   const bool read_additional_data)
{
  /**
   * This program implements the output of an 
   * EquationSystems object.  This warrants some 
   * documentation.  The output file essentially
   * consists of 12 sections:
   *
   * 1.) The type of system handled (string),
   * 2.) The number of flags that are set (unsigned int),
   *
   * for each flag in the equation system object
   *
   *   3.) the name (string)
   *
   * end flag loop
   *
   * 4.) The number of parameters that are set (unsigned int),
   *
   * for each parameter in the equation system object
   *
   *   5.) the name of the parameter  (string)
   *   6.) the value of the parameter (real)
   *
   * end parameter loop
   * 
   * 7.) The number of individual equation systems (unsigned int)
   * 
   *   for each system
   *                                                      
   *    8.)  The name of the system (string)            
   *
   *    handled through SystemBase::read():
   *
   * +-------------------------------------------------------------+
   * |  9.) The number of variables in the system (unsigned int)   |
   * |                                                             |
   * |   for each variable in the system                           |
   * |                                                             |
   * |    10.) The name of the variable (string)                   |
   * |                                                             |
   * |    11.) Combined in an FEType:                              |
   * |         - The approximation order(s) of the variable (Order |
   * |           Enum, cast to int/s)                              |
   * |         - The finite element family/ies of the variable     |
   * |           (FEFamily Enum, cast to int/s)                    |
   * |                                                             |
   * |   end variable loop                                         |
   * |                                                             |
   * | 12.) The number of additional vectors (unsigned int),       |
   * |                                                             |
   * |    for each additional vector in the equation system object |
   * |                                                             |
   * |    13.) the name of the additional vector  (string)         |
   * +-------------------------------------------------------------+
   *
   * end system loop
   *
   *
   *   for each system, handled through SystemBase::read_data():
   *   
   * +-------------------------------------------------------------+
   * | 14.) The global solution vector, re-ordered to be node-major|
   * |     (More on this later.)                                   |
   * |                                                             |
   * |    for each additional vector in the equation system object |
   * |                                                             |
   * |    15.) The global additional vector, re-ordered to be      |
   * |         node-major (More on this later.)                    |
   * +-------------------------------------------------------------+
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
   * Read the type of system handled
   */
  std::string sys_type;
      
  io.data (sys_type);
      
  if (sys_type != T_sys::system_type())
    {
      // wrong T_sys for this file
      std::cerr << "ERROR: System mismatch: This EquationSystems object handles" 
		<< std::endl
		<< " systems of type " << T_sys::system_type() 
		<< ", while the file" << std::endl
		<< " contains systems of type " << sys_type << std::endl;
      error();
    }


  /**
   * 2.)  
   *
   * Read the number of flags that are set
   */
  {
    unsigned int n_flags=0;
  
    io.data (n_flags);
  
    for (unsigned int flags=0; flags<n_flags; flags++)
      {
	/**
	 * 3.)
	 *
	 * Read the name of the ith flag
	 */
	std::string flag_name;
     
	io.data (flag_name);
       
	if (read_header)
	  this->set_flag (flag_name);
      }
  }


  /**
   * 4.)  
   *
   * Read the number of params that are set
   */
  {
    unsigned int n_params=0;
  
    io.data (n_params);
  
    for (unsigned int params=0; params<n_params; params++)
      {
        /**
	 * 5.)
	 *
	 * Read the name of the ith param
	 */
	std::string param_name;
     
	io.data (param_name);
 
	/**
	 * 6.)
	 *
	 * Read the value of the ith param
	 */
	Real param_value;
     
	io.data (param_value);

	if (read_header)
	  this->set_parameter (param_name) = param_value;
      }
  }


	  
  /**
   * 7.)  
   *
   * Read the number of equation systems
   */
  {
    unsigned int n_sys=0;
  
    io.data (n_sys);
  
    for (unsigned int sys=0; sys<n_sys; sys++)
      {
	/**
	 * 8.)
	 *
	 * Read the name of the sys-th equation system
	 */
	std::string sys_name;
      
	io.data (sys_name);
      
	if (read_header)
	  this->add_system (sys_name);

	
	/**
	 * 9.) - 11.)
	 *
	 * Let SystemBase::read() do the job
	 */
	T_sys& new_system = (*this)(sys_name);
	  
	new_system.read (io,
			 read_header,
			 read_additional_data);
      }
  }
      


  /**
   * Now we are ready to initialize the underlying data
   * structures. This will initialize the vectors for 
   * storage, the dof_map, etc...
   */ 
  if (read_header) 
    init();



  /**
   * 12.)
   *
   * Read and set the numeric vector values
   */
  if (read_data)
    for (unsigned int sys=0; sys<this->n_systems(); sys++)
      {
	T_sys& system = (*this)(sys);

	system.read_data (io,
			  read_additional_data);

      }
}















template <typename T_sys>
void EquationSystems<T_sys>::write(const std::string& name,
				   const Xdr::XdrMODE mode,
				   const bool write_data,
				   const bool write_additional_data)
{
  /**
   * This program implements the output of an 
   * EquationSystems object.  This warrants some 
   * documentation.  The output file essentially
   * consists of 12 sections:
   *
   * 1.) The type of system handled (string),
   * 2.) The number of flags that are set (unsigned int),
   *
   * for each flag in the equation system object
   *
   *   3.) the name (string)
   *
   * end flag loop
   *
   * 4.) The number of parameters that are set (unsigned int),
   *
   * for each parameter in the equation system object
   *
   *   5.) the name of the parameter  (string)
   *   6.) the value of the parameter (real)
   *
   * end parameter loop
   * 
   * 7.) The number of individual equation systems (unsigned int)
   * 
   *   for each system
   *                                                      
   *    8.)  The name of the system (string)            
   *
   *    handled through SystemBase::read():
   *
   * +-------------------------------------------------------------+
   * |  9.) The number of variables in the system (unsigned int)   |
   * |                                                             |
   * |   for each variable in the system                           |
   * |                                                             |
   * |    10.) The name of the variable (string)                   |
   * |                                                             |
   * |    11.) Combined in an FEType:                              |
   * |         - The approximation order(s) of the variable (Order |
   * |           Enum, cast to int/s)                              |
   * |         - The finite element family/ies of the variable     |
   * |           (FEFamily Enum, cast to int/s)                    |
   * |                                                             |
   * |   end variable loop                                         |
   * |                                                             |
   * | 12.) The number of additional vectors (unsigned int),       |
   * |                                                             |
   * |    for each additional vector in the equation system object |
   * |                                                             |
   * |    13.) the name of the additional vector  (string)         |
   * +-------------------------------------------------------------+
   *
   * end system loop
   *
   *
   *   for each system, handled through SystemBase::read_data():
   *   
   * +-------------------------------------------------------------+
   * | 14.) The global solution vector, re-ordered to be node-major|
   * |     (More on this later.)                                   |
   * |                                                             |
   * |    for each additional vector in the equation system object |
   * |                                                             |
   * |    15.) The global additional vector, re-ordered to be      |
   * |         node-major (More on this later.)                    |
   * +-------------------------------------------------------------+
   *
   * end system loop
   *
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will write XDR or ASCII
   * files with no changes.
   */

  Xdr io(name, mode);

  assert (io.writing());

  const unsigned int proc_id = _mesh.processor_id();
  unsigned int n_sys         = this->n_systems();

  typename std::map<std::string, T_sys*>::iterator pos = _systems.begin();
  
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
       * Write the type of system handled
       */
      {
        // set up the comment
	comment =  "# System Type";
	std::string sys_type = T_sys::system_type();
	io.data (sys_type, comment.c_str());
      }



      /**
       * 2.)  
       *
       * Write the number of flags
       */
      {
        unsigned int n_flags = this->_flags.size();
	io.data (n_flags, "# No. of Flags");
      }



      /**
       * 3.)  
       *
       * Write the flags
       */
      {
        std::set<std::string>::const_iterator flag_pos = _flags.begin();
	std::set<std::string>::const_iterator flag_end = _flags.end();
	unsigned int cnt=0;
	for (; flag_pos!= flag_end; ++flag_pos)
          {
	    comment =  "# Name, Flag ";
	    sprintf(buf, "%d", cnt++);
	    comment += buf;
	    std::string flag_name = *flag_pos;
	    io.data (flag_name, comment.c_str());
	  }
      }



      /**
       * 4.)  
       *
       * Write the number of parameters
       */
      {
        unsigned int n_params = this->_parameters.size();
	io.data (n_params, "# No. of Parameters");
      }



      /**
       * 5.) + 6.)
       *
       * Write the parameter names and values
       */
      {
        std::map<std::string, Real>::const_iterator param_pos = _parameters.begin();
	std::map<std::string, Real>::const_iterator param_end = _parameters.end();
	unsigned int cnt=0;
	for (; param_pos!= param_end; ++param_pos)
          {
	    comment =  "# Name,  Parameter No. ";
	    sprintf(buf, "%d", cnt);
	    comment += buf;
	    std::string param_name = param_pos->first;
	    io.data (param_name, comment.c_str());

	    comment = "# Value, Parameter No. ";
	    sprintf(buf, "%d", cnt++);
	    comment += buf;
	    Real param_value = param_pos->second;
	    io.data (param_value, comment.c_str());
	  }
      }



      /**
       * 7.)  
       *
       * Write the number of equation systems
       */
      io.data (n_sys, "# No. of Equation Systems");
        

      while (pos != _systems.end())
	{
	  /**
	   * 8.)
	   *
	   * Write the name of the sys_num-th system
	   */
	  {
	    const unsigned int sys_num = pos->second->number();
	    std::string sys_name       = pos->first;

	    comment =  "# Name, System No. ";
	    sprintf(buf, "%d", sys_num);
	    comment += buf;
	  
	    io.data (sys_name, comment.c_str());
	  }

	
	  /**
	   * 9.) - 13.)
	   *
	   * Let SystemBase::write() do the job
	   */
	  pos->second->write (io);

	  ++pos;
	}
    }




  /**
   * Start from the first system, again,
   * to write vectors to disk, if wanted
   */
  pos = _systems.begin();

  if (write_data)
    while (pos != _systems.end())
      {
	/**
	 * 14.) + 15.)
	 *
	 * Let SystemBase::write_data() do the job
	 */
	pos->second->write_data (io,
				 write_additional_data);
	
	++pos;
      }
}




//--------------------------------------------------------------
// Explicit instantiations
template class EquationSystems<GeneralSystem>;
template class EquationSystems<ThinSystem>;

#if defined(USE_COMPLEX_NUMBERS)
template class EquationSystems<FrequencySystem>;
#endif


