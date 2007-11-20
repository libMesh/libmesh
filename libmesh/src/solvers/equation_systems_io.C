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
#include "libmesh_logging.h"


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
   * consists of 10 sections:
     \verbatim
     1.) The number of individual equation systems (unsigned int)
     
       for each system
                                                          
        2.)  The name of the system (string)
        3.)  The type of the system (string)
    
        handled through System::read():
    
     +-------------------------------------------------------------+
     |  4.) The number of variables in the system (unsigned int)   |
     |                                                             |
     |   for each variable in the system                           |
     |                                                             |
     |    5.) The name of the variable (string)                   |
     |                                                             |
     |    6.) Combined in an FEType:                              |
     |         - The approximation order(s) of the variable (Order |
     |           Enum, cast to int/s)                              |
     |         - The finite element family/ies of the variable     |
     |           (FEFamily Enum, cast to int/s)                    |
     |                                                             |
     |   end variable loop                                         |
     |                                                             |
     | 7.) The number of additional vectors (unsigned int),       |
     |                                                             |
     |    for each additional vector in the equation system object |
     |                                                             |
     |    8.) the name of the additional vector  (string)         |
     +-------------------------------------------------------------+
    
     end system loop
    
    
       for each system, handled through System::read_data():
       
     +-------------------------------------------------------------+
     | 9.) The global solution vector, re-ordered to be node-major|
     |     (More on this later.)                                   |
     |                                                             |
     |    for each additional vector in the equation system object |
     |                                                             |
     |    10.) The global additional vector, re-ordered to be      |
     |         node-major (More on this later.)                    |
     +-------------------------------------------------------------+
    
     end system loop
   \endverbatim
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

  // Nasty hack for reading/writing zipped files
  std::string new_name = name;
  if (name.size() - name.rfind(".bz2") == 4)
    {
      new_name.erase(new_name.end() - 4, new_name.end());
      START_LOG("system(bunzip2)", "EquationSystems");
      if (libMesh::processor_id() == 0)
        {
          std::string system_string = "bunzip2 -f -k ";
          system_string += name;
          std::system(system_string.c_str());
        }
#ifdef HAVE_MPI
      MPI_Barrier(libMesh::COMM_WORLD);
#endif // HAVE_MPI
      STOP_LOG("system(bunzip2)", "EquationSystems");
    }

  
  Xdr io (new_name, mode);

  assert (io.reading());
	  
  // 1.)  
  //
  // Read the number of equation systems
  {
    unsigned int n_sys=0;
  
    io.data (n_sys);

    for (unsigned int sys=0; sys<n_sys; sys++)
      {
	// 2.)
	// Read the name of the sys-th equation system
	std::string sys_name;
      
	io.data (sys_name);
      
	// 3.)
	// Read the type of the sys-th equation system
	std::string sys_type;
	
	io.data (sys_type);
	
	if (read_header)
	  this->add_system (sys_type, sys_name);

	// 4.) - 8.)
	// Let System::read() do the job
	System& new_system = this->get_system(sys_name);
	  
	new_system.read (io,
			 read_header,
			 read_additional_data);
      }
  }
      


  // Now we are ready to initialize the underlying data
  // structures. This will initialize the vectors for 
  // storage, the dof_map, etc...
  if (read_header) 
    this->init();



  // 9.) & 10.)
  // Read and set the numeric vector values
  if (read_data)
    {
      std::map<std::string, System*>::iterator
	pos = _systems.begin();
      
      for (; pos != _systems.end(); ++pos)
	pos->second->read_data (io,
				read_additional_data);       
    }  

  // If we temporarily decompressed a .bz2 file, remove the
  // uncompressed version
  if (name.size() - name.rfind(".bz2") == 4)
    std::remove(new_name.c_str());

  // Localize each system's data
  this->update();
}















void EquationSystems::write(const std::string& name,
			    const libMeshEnums::XdrMODE mode,
                            const unsigned int write_flags) const
{
  // Currently we only support I/O in serial
  // So we're going to horribly abuse const_cast and change
  // that
  const_cast<EquationSystems*>(this)->allgather();

  /**
   * This program implements the output of an 
   * EquationSystems object.  This warrants some 
   * documentation.  The output file essentially
   * consists of 10 sections:
   \verbatim
     1.) The number of individual equation systems (unsigned int)
     
       for each system
                                                          
        2.)  The name of the system (string)            
        3.)  The type of the system (string)            
   
        handled through System::read():
   
     +-------------------------------------------------------------+
     |  4.) The number of variables in the system (unsigned int)   |
     |                                                             |
     |   for each variable in the system                           |
     |                                                             |
     |    5.) The name of the variable (string)                   |
     |                                                             |
     |    6.) Combined in an FEType:                              |
     |         - The approximation order(s) of the variable (Order |
     |           Enum, cast to int/s)                              |
     |         - The finite element family/ies of the variable     |
     |           (FEFamily Enum, cast to int/s)                    |
     |                                                             |
     |   end variable loop                                         |
     |                                                             |
     | 7.) The number of additional vectors (unsigned int),       |
     |                                                             |
     |    for each additional vector in the equation system object |
     |                                                             |
     |    8.) the name of the additional vector  (string)         |
     +-------------------------------------------------------------+
   
    end system loop
   
   
    for each system, handled through System::read_data():
       
     +-------------------------------------------------------------+
     | 9.) The global solution vector, re-ordered to be node-major|
     |     (More on this later.)                                   |
     |                                                             |
     |    for each additional vector in the equation system object |
     |                                                             |
     |    10.) The global additional vector, re-ordered to be      |
     |         node-major (More on this later.)                    |
     +-------------------------------------------------------------+

    end system loop
   \endverbatim
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

  // Nasty hack for reading/writing zipped files
  std::string new_name = name;
  if (name.size() - name.rfind(".bz2") == 4)
    new_name.erase(new_name.end() - 4, new_name.end());

  // New scope so that io will close before we try to zip the file
  {
  Xdr io(new_name, mode);

  assert (io.writing());

  const unsigned int proc_id = libMesh::processor_id();
  unsigned int n_sys         = this->n_systems();

  std::map<std::string, System*>::const_iterator
    pos = _systems.begin();
  
  std::string comment;
  char buf[80];

  // Only write the header information
  // if we are processor 0.
  if (proc_id == 0) 
    {
      // 1.)  
      // Write the number of equation systems
      io.data (n_sys, "# No. of Equation Systems");
        

      while (pos != _systems.end())
	{
	  // 2.)
	  // Write the name of the sys_num-th system
	  {
	    const unsigned int sys_num = pos->second->number();
	    std::string sys_name       = pos->first;

	    comment =  "# Name, System No. ";
	    std::sprintf(buf, "%d", sys_num);
	    comment += buf;
	  
	    io.data (sys_name, comment.c_str());
	  }


	  
	  // 3.)
	  // Write the type of system handled
	  {
	    const unsigned int sys_num = pos->second->number();
	    std::string sys_type       = pos->second->system_type();

	    comment =  "# Type, System No. ";
	    std::sprintf(buf, "%d", sys_num);
	    comment += buf;
	  
	    io.data (sys_type, comment.c_str());
	  }


	
	  // 4.) - 8.)
	  //
	  // Let System::write() do the job
	  pos->second->write (io, write_additional_data);

	  ++pos;
	}
    }




  /**
   // Start from the first system, again,
   // to write vectors to disk, if wanted
   */
  pos = _systems.begin();

  if (write_data)
    for (; pos != _systems.end(); ++pos) 
      {
	// 9.) + 10.)
	// Let System::write_data() do the job
	pos->second->write_data (io,
				 write_additional_data);
      }
  }

  // Nasty hack for reading/writing zipped files
  if (name.size() - name.rfind(".bz2") == 4)
    {
      START_LOG("system(bzip2)", "EquationSystems");
#ifdef HAVE_MPI
      MPI_Barrier(libMesh::COMM_WORLD);
#endif
      if (libMesh::processor_id() == 0)
        {
          std::string system_string = "bzip2 -f ";
          system_string += new_name;
          std::system(system_string.c_str());
        }
#ifdef HAVE_MPI
      MPI_Barrier(libMesh::COMM_WORLD);
#endif
      STOP_LOG("system(bzip2)", "EquationSystems");
    }
}
