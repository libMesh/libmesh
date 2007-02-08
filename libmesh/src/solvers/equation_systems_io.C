// $Id: equation_systems_io.C,v 1.12 2007-02-08 22:43:18 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include <cstdio> // for std::sprintf

// Local Includes
#include "equation_systems.h"
#include "xdr_cxx.h"

// Forward Declarations




// ------------------------------------------------------------
// EquationSystem class implementation
void EquationSystems::read (const std::string& name,
			    const libMeshEnums::XdrMODE mode,
                            const unsigned int read_flags)
{
  /**
   * This program implements the output of an 
   * EquationSystems object.  This warrants some 
   * documentation.  The output file essentially
   * consists of 13 sections:
   *
   * 1.) The number of flags that are set (unsigned int),
   *
   * for each flag in the equation system object
   *
   *   2.) the name (string)
   *
   * end flag loop
   *
   * 3.) The number of parameters that are set (unsigned int),
   *
   * for each parameter in the equation system object
   *
   *   4.) the name of the parameter  (string)
   *   5.) the value of the parameter (real)
   *
   * end parameter loop
   * 
   * 6.) The number of individual equation systems (unsigned int)
   * 
   *   for each system
   *                                                      
   *    7.)  The name of the system (string)
   *    8.)  The type of the system (string)
   *
   *    handled through System::read():
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
   *   for each system, handled through System::read_data():
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

   // Set booleans from the read_flags argument
   const bool read_header = read_flags & EquationSystems::READ_HEADER;
   const bool read_data   = read_flags & EquationSystems::READ_DATA;
   const bool read_additional_data 
                          = read_flags & EquationSystems::READ_ADDITIONAL_DATA;

  
  Xdr io (name, mode);

  assert (io.reading());


  /**
   * Possibly clear data structures and start from scratch.
   */
  if (read_header)
    this->clear ();


//   /**
//    * 1.)  
//    *
//    * Read the number of flags that are set
//    */
//   {
//     unsigned int n_flags=0;
  
//     io.data (n_flags);
  
//     for (unsigned int flags=0; flags<n_flags; flags++)
//       {
// 	/**
// 	 * 2.)
// 	 *
// 	 * Read the name of the ith flag
// 	 */
// 	std::string flag_name;
     
// 	io.data (flag_name);
       
// 	if (read_header)
// 	  this->parameters.set<bool> (flag_name) = true;
//       }
//   }


//   /**
//    * 3.)  
//    *
//    * Read the number of params that are set
//    */
//   {
//     unsigned int n_params=0;
  
//     io.data (n_params);
  
//     for (unsigned int params=0; params<n_params; params++)
//       {
//         /**
// 	 * 4.)
// 	 *
// 	 * Read the name of the ith param
// 	 */
// 	std::string param_name;
     
// 	io.data (param_name);
 
// 	/**
// 	 * 5.)
// 	 *
// 	 * Read the value of the ith param
// 	 */
// 	Real param_value;
     
// 	io.data (param_value);

// 	if (read_header)
// 	  this->parameters.set<Real> (param_name) = param_value;
//       }
//   }


	  
  /**
   * 6.)  
   *
   * Read the number of equation systems
   */
  {
    unsigned int n_sys=0;
  
    io.data (n_sys);

    for (unsigned int sys=0; sys<n_sys; sys++)
      {
	/**
	 * 7.)
	 *
	 * Read the name of the sys-th equation system
	 */
	std::string sys_name;
      
	io.data (sys_name);
      
	/**
	 * 8.)
	 *
	 * Read the type of the sys-th equation system
	 */
	std::string sys_type;
	
	io.data (sys_type);
	
	if (read_header)
	  this->add_system (sys_type, sys_name);

	/**
	 * 9.) - 11.)
	 *
	 * Let System::read() do the job
	 */
	System& new_system = this->get_system(sys_name);
	  
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
    this->init();



  /**
   * 12.)
   *
   * Read and set the numeric vector values
   */
  if (read_data)
    {
      std::map<std::string, System*>::iterator
	pos = _systems.begin();
      
      for (; pos != _systems.end(); ++pos)
	pos->second->read_data (io,
				read_additional_data);       
    }  

  // Localize each system's data
  this->update();
}















void EquationSystems::write(const std::string& name,
			    const libMeshEnums::XdrMODE mode,
                            const unsigned int write_flags) const
{
  /**
   * This program implements the output of an 
   * EquationSystems object.  This warrants some 
   * documentation.  The output file essentially
   * consists of 13 sections:
   *
   * 1.) The number of flags that are set (unsigned int),
   *
   * for each flag in the equation system object
   *
   *   2.) the name (string)
   *
   * end flag loop
   *
   * 3.) The number of parameters that are set (unsigned int),
   *
   * for each parameter in the equation system object
   *
   *   4.) the name of the parameter  (string)
   *   5.) the value of the parameter (real)
   *
   * end parameter loop
   * 
   * 6.) The number of individual equation systems (unsigned int)
   * 
   *   for each system
   *                                                      
   *    7.)  The name of the system (string)            
   *    8.)  The type of the system (string)            
   *
   *    handled through System::read():
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
   *   for each system, handled through System::read_data():
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

   // set booleans from write_flags argument
   const bool write_data = write_flags & EquationSystems::WRITE_DATA;
   const bool write_additional_data 
                         = write_flags & EquationSystems::WRITE_ADDITIONAL_DATA;

  Xdr io(name, mode);

  assert (io.writing());

  const unsigned int proc_id = libMesh::processor_id();
  unsigned int n_sys         = this->n_systems();

  std::map<std::string, System*>::const_iterator
    pos = _systems.begin();
  
  std::string comment;
  char buf[80];



  /**
   * Only write the header information
   * if we are processor 0.
   */
  if (proc_id == 0) 
    {
//       /**
//        * 1.)  
//        *
//        * Write the number of flags
//        */
//       {
//         unsigned int n_flags = this->_flags.size();
// 	io.data (n_flags, "# No. of Flags");
//       }



//       /**
//        * 2.)  
//        *
//        * Write the flags
//        */
//       {
//         std::set<std::string>::const_iterator flag_pos = _flags.begin();
// 	std::set<std::string>::const_iterator flag_end = _flags.end();
// 	unsigned int cnt=0;
// 	for (; flag_pos!= flag_end; ++flag_pos)
//           {
// 	    comment =  "# Name, Flag ";
// 	    std::sprintf(buf, "%d", cnt++);
// 	    comment += buf;
// 	    std::string flag_name = *flag_pos;
// 	    io.data (flag_name, comment.c_str());
// 	  }
//       }



//       /**
//        * 3.)  
//        *
//        * Write the number of parameters
//        */
//       {
//         unsigned int n_params = this->_parameters.size();
// 	io.data (n_params, "# No. of Parameters");
//       }



//       /**
//        * 4.) + 5.)
//        *
//        * Write the parameter names and values
//        */
//       {
//         std::map<std::string, Real>::const_iterator param_pos = _parameters.begin();
// 	std::map<std::string, Real>::const_iterator param_end = _parameters.end();
// 	unsigned int cnt=0;
// 	for (; param_pos!= param_end; ++param_pos)
//           {
// 	    comment =  "# Name,  Parameter No. ";
// 	    std::sprintf(buf, "%d", cnt);
// 	    comment += buf;
// 	    std::string param_name = param_pos->first;
// 	    io.data (param_name, comment.c_str());

// 	    comment = "# Value, Parameter No. ";
// 	    std::sprintf(buf, "%d", cnt++);
// 	    comment += buf;
// 	    Real param_value = param_pos->second;
// 	    io.data (param_value, comment.c_str());
// 	  }
//       }



      /**
       * 6.)  
       *
       * Write the number of equation systems
       */
      io.data (n_sys, "# No. of Equation Systems");
        

      while (pos != _systems.end())
	{
	  /**
	   * 7.)
	   *
	   * Write the name of the sys_num-th system
	   */
	  {
	    const unsigned int sys_num = pos->second->number();
	    std::string sys_name       = pos->first;

	    comment =  "# Name, System No. ";
	    std::sprintf(buf, "%d", sys_num);
	    comment += buf;
	  
	    io.data (sys_name, comment.c_str());
	  }



	  /**
	   * 8.)
	   *
	   * Write the type of system handled
	   */
	  {
	    const unsigned int sys_num = pos->second->number();
	    std::string sys_type       = pos->second->system_type();

	    comment =  "# Type, System No. ";
	    std::sprintf(buf, "%d", sys_num);
	    comment += buf;
	  
	    io.data (sys_type, comment.c_str());
	  }


	
	  /**
	   * 9.) - 13.)
	   *
	   * Let System::write() do the job
	   */
	  pos->second->write (io, write_additional_data);

	  ++pos;
	}
    }




  /**
   * Start from the first system, again,
   * to write vectors to disk, if wanted
   */
  pos = _systems.begin();

  if (write_data)
    for (; pos != _systems.end(); ++pos) 
      {
	/**
	 * 14.) + 15.)
	 *
	 * Let System::write_data() do the job
	 */
	pos->second->write_data (io,
				 write_additional_data);
      }
}
