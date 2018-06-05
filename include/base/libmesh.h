// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/parallel.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum SolverPackage : int;
}
#else
#include "libmesh/enum_solver_package.h"
#endif

// C++ includes
#include <string>
#include <vector>

// For dealing with MPI stuff in VTK.
#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
class vtkMPIController;
#endif

/**
 * The \p libMesh namespace provides an interface to certain functionality
 * in the library.  Here, it provides a LibMeshInit class which uses
 * the RAII (Resource Acquisition Is Initialization) idiom to ensure
 * initialization of any other dependent libraries (e.g. MPI or PETSC),
 * and to close those libraries when it goes out of scope.  It also
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
 * the start of your main() function.
 *
 * All libMesh functionality should be used only when a LibMeshInit
 * object exists.  Dependent library functionality, likewise, except
 * in codes which manually initialize those libraries before
 * LibMeshInit creation and finalize them after LibMeshInit
 * destruction.
 *
 * Since "it is best not to perform much more than a return rc after
 * calling MPI_Finalize", applications which want to do anything after
 * LibMeshInit destruction should manage MPI initialization and
 * finalization manually.
 */
class LibMeshInit
{
public:
#ifdef LIBMESH_HAVE_MPI
  /**
   * Initialize the library for use, with the command line options
   * provided.  This will e.g. call MPI_Init if MPI is available and
   * enabled and has not already been initialized; similar
   * initialization may take place for Petsc, Slepc, multithreading
   * support, libMesh Singleton objects, the libMesh::out/err IO
   * streams, and any libMesh handlers for floating-point exceptions,
   * signals, and/or C++ aborts.
   *
   * You must create a LibMeshInit object before using any
   * of the library functionality.  This method may take an optional
   * parameter to use a user-specified MPI communicator.
   */
  LibMeshInit(int argc, const char * const * argv,
              MPI_Comm COMM_WORLD_IN=MPI_COMM_WORLD);
#else
  LibMeshInit(int argc, const char * const * argv);
#endif

  /**
   * Destructor.  Cleans up libMesh Singleton objects, and thread
   * manager if threading is in use.  Prints reference count and
   * performance logging information if enabled.  Restores
   * pre-LibMeshInit terminate handler and floating-point-exception
   * handling.  Finalizes any of Slepc, Petsc, and MPI which were
   * initialized by LibMeshInit.
   */
  virtual ~LibMeshInit();

  /**
   * Returns the Communicator created by this libMesh object, which
   * will be a compatibility shim if MPI is not enabled, or a wrapper
   * for the user-input MPI_Comm if we were constructed with one, or a
   * wrapper for MPI_COMM_WORLD by default.
   */
  const Parallel::Communicator & comm() const { return _comm; }

  Parallel::Communicator & comm() { return _comm; }

private:
  Parallel::Communicator _comm;

#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
  // VTK object for dealing with MPI stuff in VTK.
  // This can't be a std::unique_ptr because VTK makes the destructor
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
 * \returns \p true if the argument \p arg was specified on the command line,
 * \p false otherwise.
 *
 * For backwards compatibility with past option naming conventions,
 * libMesh searches for the given argument first in its original form,
 * then with all underscores changed to dashes, then with all dashes
 * (except any leading dashes) changed to underscores, and returns
 * true if any of the above finds a match.
 *
 * This routine manipulates the command_line cursor and should not be
 * called concurrently with similar utilities in multiple threads.
 */
bool on_command_line (std::string arg);

/**
 * \returns The value associated with name on the command line if it is specified,
 * otherwise return the default, provided value.  A second template function is provided
 * to support recognizing multiple variations of a given option
 *
 * This routine manipulates the command_line cursor and should not be
 * called concurrently with similar utilities in multiple threads.
 */
template <typename T>
T command_line_value (const std::string &, T);
template <typename T>
T command_line_value (const std::vector<std::string> &, T);

/**
 * Use GetPot's search()/next() functions to get following arguments
 * from the command line.
 *
 * For backwards compatibility with past option naming conventions,
 * libMesh searches for the given argument first in its original form,
 * then with all underscores changed to dashes, then with all dashes
 * (except any leading dashes) changed to underscores, and returns
 * true if any of the above finds a match.
 *
 * This routine manipulates the command_line cursor and should not be
 * called concurrently with similar utilities in multiple threads.
 */
template <typename T>
T command_line_next (std::string name, T default_value);

/**
 * \returns The array of values associated with name on the command line if it is specified,
 * otherwise return the default, provided array.
 *
 * This routine manipulates the command_line cursor and should not be
 * called concurrently with similar utilities in multiple threads.
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
 * \returns The default solver interface to use.  The value depends on
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
