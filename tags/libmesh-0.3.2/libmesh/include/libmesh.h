// $Id: libmesh.h,v 1.3 2003-02-24 14:35:50 benkirk Exp $

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



#ifndef __libmesh_h__
#define __libmesh_h__


// C++ includes

// Local includes
#include "mesh_common.h"
#include "libmesh_base.h"
#include "perf_log.h"





/**
 * The \p libMesh class provides an interface to certain functionality
 * in the library.  It provides a uniform \p init() method that
 * initializes any other dependent libraries (e.g. MPI or PETSC),
 * and a \p close() method for closing those libraries.  It also
 * provides a centralized place for performance logging and other
 * functionality.
 */
class libMesh : public libMeshBase
{

private:

  /**
   * This class only contains static members and should never be
   * instantiated nor derived from, so it has a private default constructor.
   */
  libMesh () {}
  
public:
  
  /**
   * Initialize the library for use.  This will call
   * PetscInitialize if PETSC is available.  You must call
   * this method before using any of the library functionality.
   */
  static void init (int & argc, char** & argv);

  /**
   * Checks that the \p init() member has been called.  If it
   * hasn't an error message is printed and the code aborts.
   * It is useful to \p assert(libMesh::initialized()) in library
   * object constructors.
   */
  static bool initialized ();
  
  /**
   * Stop using the mesh library.  This will call PetscFinalize()
   * if PETSC is available.  This method should be called after
   * all other library objects have gone out of scope, as it
   * interrogates the \p ReferenceCounter object to look for memory
   * leaks.
   */
  static int close ();

  /**
   * Checks that the \p close() member has been called. This should
   * always return false when called from a library object.
   * It is useful to \p assert(!libMesh::closed()) in library
   * object destructors.
   */
  static bool closed ();
  
  /**
   * A \p PerfLog object to log performance.  If the library is configured
   * with \p --enable-perflog then it will log key functions.
   */
  static PerfLog log;

private:

  /**
   * Flag that tells if \p init() has been called.
   */
  static bool _is_initialized;
};



// ------------------------------------------------------------
// libMesh inline member functions
inline
bool libMesh::initialized()
{
  if (!_is_initialized)
    {
      std::cerr << "ERROR: You must call libMesh::init() before using any of the library!"
		<< std::endl
		<< std::endl;
      error();      
    }

  // Note we can only get here if the library is initialized
  return true;
}



inline
bool libMesh::closed()
{

  if (!_is_initialized)
    {
      std::cerr << "ERROR: libMesh::close() called prematurely.  Library objects still in scope!"
		<< std::endl
		<< std::endl;
      
      error();
    }

  // Note we can only get here if the library is initialized
  return false;
}



#endif // #define __libmesh_h__
