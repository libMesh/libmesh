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
#include "mesh_base.h"
#include "mesh_tools.h"
#include "parallel_mesh.h"
#include "parallel.h"
#include "serial_mesh.h"
#include "xdr_cxx.h"

// Forward Declarations

// Anonymous namespace for implementation details.
namespace {
  std::string local_file_name (const std::string &name)
  {
    std::string basename(name);
    char buf[256];

    if (basename.size() - basename.rfind(".bz2") == 4)
      {
	basename.erase(basename.end()-4, basename.end());
	std::sprintf(buf, "%s.%04d.bz2", basename.c_str(), libMesh::processor_id());
      }
    else if (basename.size() - basename.rfind(".gz") == 3)
      {
	basename.erase(basename.end()-3, basename.end());
	std::sprintf(buf, "%s.%04d.gz", basename.c_str(), libMesh::processor_id());
      }
    else
      std::sprintf(buf, "%s.%04d", basename.c_str(), libMesh::processor_id());
      
    return std::string(buf);
  }
}




// ------------------------------------------------------------
// EquationSystem class implementation
void EquationSystems::read (const std::string& name,
			    const libMeshEnums::XdrMODE mode,
                            const unsigned int read_flags)
{
#ifdef LIBMESH_ENABLE_EXCEPTIONS

  // If we have exceptions enabled we can be considerate and try
  // to read old restart files which contain infinite element
  // information but do not have the " with infinite elements"
  // string in the version information.
  
  // First try the read the user requested
  try
    {
      this->_read_impl (name, mode, read_flags);
    }

  // If that fails, try it again but explicitly request we look for infinite element info      
  catch (...)
    {
      std::cout << "\n*********************************************************************\n"
		<< "READING THE FILE \"" << name << "\" FAILED.\n"
		<< "It is possible this file contains infinite element information,\n"
		<< "but the version string does not contain \" with infinite elements\"\n"
		<< "Let's try this again, but looking for infinite element information...\n"
		<< "*********************************************************************\n"
		<< std::endl;

      try
	{
	  this->_read_impl (name, mode, read_flags | EquationSystems::TRY_READ_IFEMS);
	}

      // If all that failed, we are out of ideas here...
      catch (...)
	{
	  std::cout << "\n*********************************************************************\n"
		    << "Well, at least we tried!\n"
		    << "Good Luck!!\n"
		    << "*********************************************************************\n"
		    << std::endl;
	  throw;
	}
    }

#else

  // no exceptions - cross your fingers...
  this->_read_impl (name, mode, read_flags);
  
#endif // #ifdef LIBMESH_ENABLE_EXCEPTIONS
}



void EquationSystems::_read_impl (const std::string& name,
				  const libMeshEnums::XdrMODE mode,
				  const unsigned int read_flags)
{
  /**
   * This program implements the output of an 
   * EquationSystems object.  This warrants some 
   * documentation.  The output file essentially
   * consists of 11 sections:
     \verbatim
     1.) A version header (for non-'legacy' formats, libMesh-0.7.0 and greater). 
     2.) The number of individual equation systems (unsigned int)
     
       for each system
                                                          
        3.)  The name of the system (string)
        4.)  The type of the system (string)
    
        handled through System::read():
    
     +-------------------------------------------------------------+
     |  5.) The number of variables in the system (unsigned int)   |
     |                                                             |
     |   for each variable in the system                           |
     |                                                             |
     |    6.) The name of the variable (string)                    |
     |                                                             |
     |    7.) Combined in an FEType:                               |
     |         - The approximation order(s) of the variable (Order |
     |           Enum, cast to int/s)                              |
     |         - The finite element family/ies of the variable     |
     |           (FEFamily Enum, cast to int/s)                    |
     |                                                             |
     |   end variable loop                                         |
     |                                                             |
     | 8.) The number of additional vectors (unsigned int),        |
     |                                                             |
     |    for each additional vector in the equation system object |
     |                                                             |
     |    9.) the name of the additional vector  (string)          |
     +-------------------------------------------------------------+
    
     end system loop
    
    
       for each system, handled through System::read_{serialized,parallel}_data():
       
     +--------------------------------------------------------------+
     | 10.) The global solution vector, re-ordered to be node-major |
     |     (More on this later.)                                    |
     |                                                              |
     |    for each additional vector in the equation system object  |
     |                                                              |
     |    11.) The global additional vector, re-ordered to be       |
     |         node-major (More on this later.)                     |
     +--------------------------------------------------------------+
    
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
   const bool read_header          = read_flags & EquationSystems::READ_HEADER;
   const bool read_data            = read_flags & EquationSystems::READ_DATA;
   const bool read_additional_data = read_flags & EquationSystems::READ_ADDITIONAL_DATA;
   const bool read_legacy_format   = read_flags & EquationSystems::READ_LEGACY_FORMAT;
   const bool try_read_ifems       = read_flags & EquationSystems::TRY_READ_IFEMS;
         bool read_parallel_files  = false;
	 
  // This will unzip a file with .bz2 as the extension, otherwise it
  // simply returns the name if the file need not be unzipped.
  Xdr io ((libMesh::processor_id() == 0) ? name : "", mode);
  libmesh_assert (io.reading());
	  
  {
    // 1.)
    // Read the version header.
    std::string version = "legacy";
    if (!read_legacy_format)
      {
	if (libMesh::processor_id() == 0) io.data(version);	
	Parallel::broadcast(version);
	
	// All processors have the version header, if it does not contain
	// "libMesh" something then it is a legacy file.
	if (!(version.find("libMesh") < version.size()))
	  {
	    io.close();
	    
	    // Recursively call this read() function but with the 
	    // EquationSystems::READ_LEGACY_FORMAT bit set.
	    this->read (name, mode, (read_flags | EquationSystems::READ_LEGACY_FORMAT));
	    return;
	  }

	read_parallel_files = (version.rfind(" parallel") < version.size());

	// If requested that we try to read infinite element information,
	// and the string " with infinite elements" is not in the version,
	// then tack it on.  This is for compatibility reading ifem
	// files written prior to 11/10/2008 - BSK
	if (try_read_ifems)
	  if (!(version.rfind(" with infinite elements") < version.size()))
	    version += " with infinite elements";
	
      }
    else
      deprecated();

  START_LOG("read()","EquationSystems");
      
    // 2.)  
    // Read the number of equation systems
    unsigned int n_sys=0;
    if (libMesh::processor_id() == 0) io.data (n_sys);
    Parallel::broadcast(n_sys);

    for (unsigned int sys=0; sys<n_sys; sys++)
      {
	// 3.)
	// Read the name of the sys-th equation system
	std::string sys_name;      
	if (libMesh::processor_id() == 0) io.data (sys_name);
	Parallel::broadcast(sys_name);
	
	// 4.)
	// Read the type of the sys-th equation system
	std::string sys_type;	
	if (libMesh::processor_id() == 0) io.data (sys_type);
	Parallel::broadcast(sys_type);
	
	if (read_header)
	  this->add_system (sys_type, sys_name);

	// 5.) - 9.)
	// Let System::read_header() do the job
	System& new_system = this->get_system(sys_name);	  
	new_system.read_header (io,
				version,
				read_header,
				read_additional_data,
				read_legacy_format);
      }
  }
      


  // Now we are ready to initialize the underlying data
  // structures. This will initialize the vectors for 
  // storage, the dof_map, etc...
  if (read_header) 
    this->init();

  // 10.) & 11.)
  // Read and set the numeric vector values
  if (read_data)
    {
      // the EquationSystems::read() method should look constant from the mesh
      // perspective, but we need to assign a temporary numbering to the nodes
      // and elements in the mesh, which requires that we abuse const_cast
      if (!read_legacy_format)
	{
	  MeshBase &mesh = const_cast<MeshBase&>(this->get_mesh());
	  MeshTools::Private::globally_renumber_nodes_and_elements(mesh);
	}

      Xdr local_io (read_parallel_files ? local_file_name(name) : "", mode);
   
      std::map<std::string, System*>::iterator
	pos = _systems.begin();
      
      for (; pos != _systems.end(); ++pos)
	if (read_legacy_format)
	  {
	    deprecated();
	    pos->second->read_legacy_data (io, read_additional_data);
	  }
	else
	  if (read_parallel_files)
	    pos->second->read_parallel_data   (local_io, read_additional_data);
	  else
	    pos->second->read_serialized_data (io, read_additional_data);

      
      // Undo the temporary numbering.
      if (!read_legacy_format)
	{
	  if (dynamic_cast<ParallelMesh*>(const_cast<MeshBase*>(&_mesh)))
	    {
	      ParallelMesh *mesh = dynamic_cast<ParallelMesh*>(const_cast<MeshBase*>(&_mesh));    
	      MeshTools::Private::fix_broken_node_and_element_numbering(*mesh);
	    }
	  else if (dynamic_cast<SerialMesh*>(const_cast<MeshBase*>(&_mesh)))
	    {
	      SerialMesh *mesh = dynamic_cast<SerialMesh*>(const_cast<MeshBase*>(&_mesh));    
	      MeshTools::Private::fix_broken_node_and_element_numbering(*mesh);
	    }
	  else
	    {
	      std::cerr << "ERROR:  dynamic_cast<> to ParallelMesh and SerialMesh failed!"
			<< std::endl;
	      libmesh_error();
	    }	  
	}
    }  

  STOP_LOG("read()","EquationSystems");
  
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
   * consists of 11 sections:
   \verbatim
     1.) The version header.
     2.) The number of individual equation systems (unsigned int)
     
       for each system
                                                          
        3.)  The name of the system (string)            
        4.)  The type of the system (string)            
   
        handled through System::read():
   
     +-------------------------------------------------------------+
     |  5.) The number of variables in the system (unsigned int)   |
     |                                                             |
     |   for each variable in the system                           |
     |                                                             |
     |    6.) The name of the variable (string)                    |
     |                                                             |
     |    7.) Combined in an FEType:                               |
     |         - The approximation order(s) of the variable (Order |
     |           Enum, cast to int/s)                              |
     |         - The finite element family/ies of the variable     |
     |           (FEFamily Enum, cast to int/s)                    |
     |                                                             |
     |   end variable loop                                         |
     |                                                             |
     | 8.) The number of additional vectors (unsigned int),        |
     |                                                             |
     |    for each additional vector in the equation system object |
     |                                                             |
     |    9.) the name of the additional vector  (string)          |
     +-------------------------------------------------------------+
   
    end system loop
   
   
    for each system, handled through System::write_{serialized,parallel}_data():
       
     +--------------------------------------------------------------+
     | 10.) The global solution vector, re-ordered to be node-major |
     |     (More on this later.)                                    |
     |                                                              |
     |    for each additional vector in the equation system object  |
     |                                                              |
     |    11.) The global additional vector, re-ordered to be       |
     |         node-major (More on this later.)                     |
     +--------------------------------------------------------------+

    end system loop
   \endverbatim
   *
   * Note that the actual IO is handled through the Xdr class 
   * (to be renamed later?) which provides a uniform interface to 
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will write XDR or ASCII
   * files with no changes.
   */
  
  // the EquationSystems::write() method should look constant,
  // but we need to assign a temporary numbering to the nodes
  // and elements in the mesh, which requires that we abuse const_cast
  {
    MeshBase &mesh = const_cast<MeshBase&>(this->get_mesh());
    MeshTools::Private::globally_renumber_nodes_and_elements(mesh);
  }
  
   // set booleans from write_flags argument
   const bool write_data            = write_flags & EquationSystems::WRITE_DATA;
   const bool write_parallel_files  = write_flags & EquationSystems::WRITE_PARALLEL_FILES;
   const bool write_additional_data = write_flags & EquationSystems::WRITE_ADDITIONAL_DATA;

  // New scope so that io will close before we try to zip the file
  {
    Xdr io((libMesh::processor_id()==0) ? name : "", mode);    
    libmesh_assert (io.writing());
    
    START_LOG("write()","EquationSystems");

    const unsigned int proc_id = libMesh::processor_id();
    unsigned int n_sys         = this->n_systems();
    
    std::map<std::string, System*>::const_iterator
      pos = _systems.begin();
    
    std::string comment;
    char buf[256];
    
    // Only write the header information
    // if we are processor 0.
    if (proc_id == 0) 
      {
	// 1.)
	// Write the version header
	std::string version = "libMesh-0.7.0+";
	if (write_parallel_files) version += " parallel";
	
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
	version += " with infinite elements";
#endif
	io.data (version, "# File Format Identifier");
	
	// 2.)  
	// Write the number of equation systems
	io.data (n_sys, "# No. of Equation Systems");
        
	while (pos != _systems.end())
	  {
	    // 3.)
	    // Write the name of the sys_num-th system
	    {
	      const unsigned int sys_num = pos->second->number();
	      std::string sys_name       = pos->first;
	      
	      comment =  "# Name, System No. ";
	      std::sprintf(buf, "%d", sys_num);
	      comment += buf;
	  
	      io.data (sys_name, comment.c_str());
	    }
	  
	    // 4.)
	    // Write the type of system handled
	    {
	      const unsigned int sys_num = pos->second->number();
	      std::string sys_type       = pos->second->system_type();

	      comment =  "# Type, System No. ";
	      std::sprintf(buf, "%d", sys_num);
	      comment += buf;
	  
	      io.data (sys_type, comment.c_str());
	    }
	
	    // 5.) - 9.)
	    // Let System::write_header() do the job
	    pos->second->write_header (io, version, write_additional_data);
	    
	    ++pos;
	  }
      }

    // Start from the first system, again,
    // to write vectors to disk, if wanted
    if (write_data)
      {
	// open a parallel buffer if warranted.
	Xdr local_io (write_parallel_files ? local_file_name(name) : "", mode);
	
	for (pos = _systems.begin(); pos != _systems.end(); ++pos) 
	  {
	    // 10.) + 11.)
	    if (write_parallel_files)
	      pos->second->write_parallel_data (local_io,write_additional_data);
	    else
	      pos->second->write_serialized_data (io,write_additional_data);
	  }
      }
	
    STOP_LOG("write()","EquationSystems");
  }
 
  // the EquationSystems::write() method should look constant,
  // but we need to undo the temporary numbering of the nodes
  // and elements in the mesh, which requires that we abuse const_cast
  if (dynamic_cast<ParallelMesh*>(const_cast<MeshBase*>(&_mesh)))
    {
      ParallelMesh *mesh = dynamic_cast<ParallelMesh*>(const_cast<MeshBase*>(&_mesh));    
      MeshTools::Private::fix_broken_node_and_element_numbering(*mesh);
    }
  else if (dynamic_cast<SerialMesh*>(const_cast<MeshBase*>(&_mesh)))
    {
      SerialMesh *mesh = dynamic_cast<SerialMesh*>(const_cast<MeshBase*>(&_mesh));    
      MeshTools::Private::fix_broken_node_and_element_numbering(*mesh);
    }
  else
    {
      std::cerr << "ERROR:  dynamic_cast<> to ParallelMesh and SerialMesh failed!"
		<< std::endl;
      libmesh_error();
    }
}
