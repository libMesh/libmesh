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



#ifndef __libmesh_h__
#define __libmesh_h__


// C++ includes
#include <string>

// Local includes
#include "libmesh_common.h"
#include "libmesh_base.h"
#include "perf_log.h"
#include "enum_solver_package.h"

/**
 * The \p libMesh namespace provides an interface to certain functionality
 * in the library.  It provides a uniform \p init() method that
 * initializes any other dependent libraries (e.g. MPI or PETSC),
 * and a \p close() method for closing those libraries.  It also
 * provides a centralized place for performance logging and other
 * functionality.
 */
namespace libMesh
{


// Forward declaration, required when included
// in perf_log.{C,h} because the preceeding
// #include "perf_log.h" is ineffective.
// Multiple inclusion avoidance problem...
// __perf_log_h__ already #define'd, but the
// class has not been declared yet!.
class PerfLog;



/**
 * The \p LibMeshInit class, when constructed, initializes
 * the dependent libraries (e.g. MPI or PETSC) and does the 
 * command line parsing needed by libMesh.  The LibMeshInit
 * destructor closes those libraries properly.
 *
 * For most users, a single LibMeshInit object should be created at
 * the start of your main() function.  This object replaces the
 * previous libMesh::init()/libMesh::close() methods, which are
 * now deprecated.
 */
class LibMeshInit
{
public:
#ifdef LIBMESH_HAVE_MPI
  /**
   * Initialize the library for use, with the command line options
   * provided.  This will e.g. call PetscInitialize if PETSC is
   * available.  You must create a LibMeshInit object before using any
   * of the library functionality.  This method may take an optional
   * parameter to use a user-specified MPI communicator.
   */
  LibMeshInit(int & argc, char** & argv,
	      MPI_Comm COMM_WORLD_IN=MPI_COMM_WORLD);
#else
  LibMeshInit(int & argc, char** & argv);
#endif
	
  virtual ~LibMeshInit();
};


#ifndef LIBMESH_HAVE_MPI
  
  /**
   * Initialize the library for use.  This will call
   * PetscInitialize if PETSC is available.
   *
   * You must perform an initialization before using any of the
   * library functionality, but libMesh::init() is a deprecated
   * way to do so.  Create a LibMeshInit object instead.
   */
  void init (int & argc, char** & argv);

#else
  
  /**
   * Initialize the library for use.  This will call
   * PetscInitialize if PETSC is available.
   * This method takes an optional parameter 
   *
   * You must perform an initialization before using any of the
   * library functionality, but libMesh::init() is a deprecated
   * way to do so.  Create a LibMeshInit object instead.
   */
  void init (int & argc, char** & argv,
	     MPI_Comm COMM_WORLD_IN=MPI_COMM_WORLD);

#endif

  /**
   * Checks that library initialization has been done.  If it
   * hasn't an error message is printed and the code aborts.
   * It is useful to \p libmesh_assert(libMesh::initialized()) in library
   * object constructors.
   */
  bool initialized ();
  
  /**
   * Stop using the mesh library.  This will call PetscFinalize()
   * if PETSC is available.  This method should be called after
   * all other library objects have gone out of scope, as it
   * interrogates the \p ReferenceCounter object to look for memory
   * leaks.
   *
   * libMesh::init() and libMesh::close() are a deprecated method
   * of library initialization.  Create a LibMeshInit object to
   * begin using the library; when the LibMeshInit object is
   * destroyed the library will be closed.
   */
  int close ();

  /**
   * Checks that the library has been closed.  This should
   * always return false when called from a library object.
   * It is useful to \p libmesh_assert(!libMesh::closed()) in library
   * object destructors.
   */
  bool closed ();
  
  /**
   * A \p PerfLog object to log performance.  If the library is configured
   * with \p --enable-perflog then it will log key functions.
   */
  extern PerfLog perflog;

  /**
   * @returns true if the argument \p arg was specified on the command line,
   * \p false otherwise.
   */
  bool on_command_line (const std::string& arg);  

  /**
   * \returns the value associated with name on the command line if it is specified,
   * otherwise return the default, provided value.
   */
  template <typename T>
  T command_line_value (const std::string &, T);  

  /**
   * The imaginary unit, \f$ \sqrt{-1} \f$.
   */
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
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


} // namespace libMesh



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
