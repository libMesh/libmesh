// $Id: libmesh.h,v 1.2 2004-01-03 15:37:41 benkirk Exp $

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



#ifndef __libmesh_h__
#define __libmesh_h__


// C++ includes
#include <iostream>
#include <string>

// Local includes
#include "libmesh_common.h"
#include "libmesh_base.h"
#include "perf_log.h"
#include "enum_solver_package.h"


// Forward declaration, required when included
// in perf_log.{C,h} because the preceeding
// #include "perf_log.h" is ineffective.
// Multiple inclusion avoidance problem...
// __perf_log_h__ already #define'd, but the
// class has not been declared yet!.
class PerfLog;



/**
 * The \p libMesh namespace provides an interface to certain functionality
 * in the library.  It provides a uniform \p init() method that
 * initializes any other dependent libraries (e.g. MPI or PETSC),
 * and a \p close() method for closing those libraries.  It also
 * provides a centralized place for performance logging and other
 * functionality.
 */
namespace libMesh {
  /**
   * Initialize the library for use.  This will call
   * PetscInitialize if PETSC is available.  You must call
   * this method before using any of the library functionality.
   */
  void init (int & argc, char** & argv);

  /**
   * Checks that the \p init() member has been called.  If it
   * hasn't an error message is printed and the code aborts.
   * It is useful to \p assert(libMesh::initialized()) in library
   * object constructors.
   */
  bool initialized ();
  
  /**
   * Stop using the mesh library.  This will call PetscFinalize()
   * if PETSC is available.  This method should be called after
   * all other library objects have gone out of scope, as it
   * interrogates the \p ReferenceCounter object to look for memory
   * leaks.
   */
  int close ();

  /**
   * Checks that the \p close() member has been called. This should
   * always return false when called from a library object.
   * It is useful to \p assert(!libMesh::closed()) in library
   * object destructors.
   */
  bool closed ();
  
  /**
   * A \p PerfLog object to log performance.  If the library is configured
   * with \p --enable-perflog then it will log key functions.
   */
  extern PerfLog log;

  /**
   * @returns true if the argument \p arg was specified on the command line,
   * \p false otherwise.
   */
  bool on_command_line (const std::string& arg);  

  /**
   * The imaginary unit, \f$ \sqrt{-1} \f$.
   */
#ifdef USE_COMPLEX_NUMBERS
  extern const Number imaginary;
#endif

  /**
   * @returns the default solver interface to use.  The value depends on
   * which solver packages  were available when the library was configured.
   * The command-line is also checked, allowing the user to override the
   * compiled default.  For example, \p --use-petsc will force the use of
   * PETSc solvers, and \p --use-laspack will force the use of LASPACK
   * solvers.   
   */
  SolverPackage default_solver_package ();

  /**
   * \f$ \pi=3.14159... \f$.
   */
  extern const Real pi;

  /**
   * \f$ zero=0. \f$.
   */
  extern const Number zero;

  /**
   * A number which is used quite often to represent
   * an invalid or uninitialized value.
   */
  extern const unsigned int invalid_uint;
  
  /**
   * Namespaces don't provide private data,
   * so let's take the data we would like
   * private and put it in an obnoxious
   * namespace.  At least that way it is a
   * pain to use, thus discouraging errors.
   */
  namespace libMeshPrivateData {
    
    /**
     * Flag that tells if \p init() has been called.
     */
    extern bool _is_initialized;
    
    /**
     * The default solver package to use.
     */
    extern SolverPackage _solver_package;
  }
}



// ------------------------------------------------------------
// libMesh inline member functions
inline
bool libMesh::initialized()
{
  return libMeshPrivateData::_is_initialized;
}



inline
bool libMesh::closed()
{
  return !libMesh::initialized();
}


#endif // #define __libmesh_h__
