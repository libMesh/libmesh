// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_LIBMESH_H
#define LIBMESH_LIBMESH_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_base.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/parallel.h"

// C++ includes
#include <string>
#include <vector>

// For dealing with MPI stuff in VTK.
#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
class vtkMPIController;
#endif

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


/**
 * The \p LibMeshInit class, when constructed, initializes
 * the dependent libraries (e.g. MPI or PETSC) and does the
 * command line parsing needed by libMesh.  The LibMeshInit
 * destructor closes those libraries properly.
 *
 * For most users, a single LibMeshInit object should be created at
 * the start of your main() function.  This object replaces the
 * previous \code libMesh::init()/libMesh::close() \endcode methods,
 * which are now deprecated.
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
  LibMeshInit(int argc, const char * const * argv,
              MPI_Comm COMM_WORLD_IN=MPI_COMM_WORLD);
#else
  LibMeshInit(int argc, const char * const * argv);
#endif

  virtual ~LibMeshInit();

  const Parallel::Communicator & comm() const { return _comm; }

  Parallel::Communicator & comm() { return _comm; }

private:
  Parallel::Communicator _comm;

#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
  // VTK object for dealing with MPI stuff in VTK.
  // This can't be a UniquePtr because VTK makes the destructor
  // protected and forces us to use a named destructor manually
  vtkMPIController * _vtk_mpi_controller;
#endif
};

/**
 * Checks that library initialization has been done.  If it
 * hasn't an error message is printed and the code aborts.
 * It is useful to \p libmesh_assert(libMesh::initialized()) in library
 * object constructors.
 */
bool initialized ();

/**
 * Checks that the library has been closed.  This should
 * always return false when called from a library object.
 * It is useful to \p libmesh_assert(!libMesh::closed()) in library
 * object destructors.
 */
bool closed ();

/**
 * Toggle hardware trap floating point exceptions
 */
void enableFPE(bool on);

/**
 * Toggle libMesh reporting of segmentation faults
 */
void enableSEGV(bool on);

/**
 * @returns true if the argument \p arg was specified on the command line,
 * \p false otherwise.
 */
bool on_command_line (const std::string & arg);

/**
 * \returns the value associated with name on the command line if it is specified,
 * otherwise return the default, provided value.  A second template function is provided
 * to support recognizing multiple variations of a given option
 */
template <typename T>
T command_line_value (const std::string &, T);
template <typename T>
T command_line_value (const std::vector<std::string> &, T);

/**
 * Use GetPot's search()/next() functions to get following arguments
 * from the command line.
 */
template <typename T>
T command_line_next (const std::string &, T);

/**
 * \returns the array of values associated with name on the command line if it is specified,
 * otherwise return the default, provided array.
 */
template <typename T>
void command_line_vector (const std::string &, std::vector<T> &);

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
const Real pi =
  static_cast<Real>(3.1415926535897932384626433832795029L);

/**
 * \f$ zero=0. \f$.
 */
const Number zero = 0.;

/**
 * A number which is used quite often to represent
 * an invalid or uninitialized value.
 */
const unsigned int invalid_uint = static_cast<unsigned int>(-1);

} // namespace libMesh

#endif // LIBMESH_LIBMESH_H
