// $Id: libmesh.h,v 1.8 2003-04-08 03:21:35 benkirk Exp $

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
#include <iostream>

// Local includes
#include "mesh_common.h"
#include "libmesh_base.h"
#include "perf_log.h"
#include "auto_ptr.h"



// Forward declaration, required when included
// in perf_log.{C,h} because the preceeding
// #include "perf_log.h" is ineffective.
// Multiple inclusion avoidance problem...
// __perf_log_h__ already #define'd, but the
// class has not been declared yet!.
class PerfLog;

// Forward declaration required for the GetPot
// command line parser.
class GetPot;



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

  /**
   * @returns true if the argument \p arg was specified on the command line,
   * \p false otherwise.
   */
  static bool on_command_line (const std::string& arg);
  
  /**
   * @returns an \p std::ostream to print messages with.
   */
  static std::ostream & msg ();
  
  /**
   * @returns an \p std::ostream to print messages with.
   */
  static std::ostream & err ();
  

#if defined(USE_COMPLEX_NUMBERS)
  /**
   * The imaginary unit, \f$ \sqrt{-1} \f$.
   */
  static const Number imaginary;
  
#endif

  /**
   * \f$ \pi=3.14159... \f$.
   */
  static const Real pi;


  /**
   * \f$ \zero=0. \f$.
   */
  static const Number zero;

  
private:

  
  /**
   * A \p GetPot object that contains the command-line arguments
   * as passed to the \p init member.  This object can be used to
   * parse the command line for user-specified run-time options.
   * This is stored as an \p AutoPtr so this header file does
   * not need to actually include the heavyweight \p GetPot header. 
   */
  static AutoPtr<GetPot> _command_line;

  /**
   * A \p std::ostream to print messages into.  This is useful
   * since it allows messages to be printed just from processor 0.
   */
  static std::ostream *_out;
  
  /**
   * A \p std::ostream to print messages into.  This is useful
   * since it allows messages to be printed just from processor 0.
   */
  static std::ostream *_err;

  
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



inline
std::ostream & libMesh::msg ()
{
  assert (_out != NULL);

  return *_out;
}



inline
std::ostream & libMesh::err ()
{
  assert (_err != NULL);

  return *_err;
}


#endif // #define __libmesh_h__
