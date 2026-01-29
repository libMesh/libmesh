// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_LIBMESH_EXCEPTIONS_H
#define LIBMESH_LIBMESH_EXCEPTIONS_H

#include "libmesh/libmesh_config.h"

#include "libmesh/libmesh_abort.h"

#include <stdexcept>
#include <string>
#include <sstream>

namespace libMesh {

/**
 * A terminate handler.  libMesh sets this to handle uncaught
 * exceptions; it can also be called manually to cleanup, print
 * any diagnostics, do cleanup, and abort.
 *
 * If an uncaught exception is a TerminationException, as thrown by
 * libmesh_terminate(), the handler avoids any diagnostic output.
 *
 * If an uncaught exception is a std::exception, its message is
 * printed, followed by stack trace and performance log output.
 */
void libmesh_terminate_handler();

/**
 * Toggle hardware trap floating point exceptions
 */
void enableFPE(bool on);

/**
 * Toggle libMesh reporting of segmentation faults
 */
void enableSEGV(bool on);

/**
 * A class to represent the internal "this should never happen"
 * errors, to be thrown by "libmesh_error();"
 */
class LogicError : public std::logic_error
{
public:
  LogicError() : std::logic_error( "Error in libMesh internal logic" ) {}
  LogicError(const std::string & msg) : std::logic_error( msg ) {}
};


/**
 * A class to stub for features that should be in libMesh, but
 * haven't been written yet, to be thrown by
 * "libmesh_not_implemented();"
 */
class NotImplemented : public std::logic_error
{
public:
  NotImplemented(std::string msg="") : std::logic_error( "Error: feature not implemented!\n" + msg ) {}
};


/**
 * A class representing a failed attempt by the library to open a
 * file (or construct an fstream, etc), to be thrown by
 * "libmesh_file_error(filename);" For ease of debugging, "filename"
 * should include any (absolute or relative or implicit) pathname
 * that was part of the failed open.
 */
class FileError : public std::runtime_error
{
public:
  FileError(const std::string & filename, const std::string msg="") :
    std::runtime_error("Error with file `" + filename + "'\n" + std::move(msg)) {}
};


/**
 * A class representing the detection of an unexpected degeneracy,
 * e.g. a negative-determinant Jacobian in a map expected to be
 * positive, or a non-trivial kernel in a map expected to be a
 * bijection (such as a singular matrix).
 *
 * libMesh::FEMap throws this if it encounters a point xi in an
 * element's "master space" at which the mapping to physical space has
 * a too-small (negative, or zero, or nearly zero) Jacobian
 * determinant, where "too-small" is determined by a particular
 * library method's assigned tolerance.
 *
 * libMesh::DenseMatrix throws this if it is asked to solve a system
 * with a singular matrix and a method (such as lu_solve()) that
 * cannot handle singularities.
 */
class DegenerateMap : public std::runtime_error
{
public:
  DegenerateMap(std::string msg="") :
    std::runtime_error( "Degenerate map, e.g. negative Jacobian or singular matrix.\n" + msg ) {}
};


/**
 * A class representing a solver's failure to converge, to be thrown
 * by "libmesh_convergence_failure();"  This should be a last
 * resort; more often, a solve which has failed should be
 * reattempted after switching to a smaller timestep, adding
 * underrelaxation, taking a smaller continuation step, etc.
 */
class ConvergenceFailure : public std::runtime_error
{
public:
  ConvergenceFailure(const std::string & err_msg="Unrecoverable failure to converge") : std::runtime_error( err_msg ) {}
};


/**
 * A class representing that a dynamic cast failed to produce expected output.
 */
class DynamicCastFailure:  public std::runtime_error
{
public:
  DynamicCastFailure() : std::runtime_error( "Failed dynamic cast!" ) {}
};

/**
 * A class representing a floating point exception.
 */
class FloatingPointException: public std::runtime_error
{
public:
  FloatingPointException() : std::runtime_error( "libmesh FPE!" ) {}
};

/**
 * A class representing an exception during a solve.
 */
class SolverException: public std::exception
{
public:
  SolverException(int error_code_in) :
    std::exception(),
    error_code(error_code_in)
  {
    std::ostringstream oss;
    oss << "Error code " << error_code << " during solve." << std::endl;
    what_message = oss.str();
  }

  /**
   * Virtual destructor, gotta have one of those.
   */
  virtual ~SolverException() = default;

  /**
   * Override the what() function to provide a generic error message.
   */
  virtual const char * what() const noexcept override
  {
    // std::string::c_str() is noexcept in C++11, so it's safe to call
    // in what() because it can't throw.
    return what_message.c_str();
  }

  /**
   * The error code generated by the solver.
   */
  int error_code;

  /**
   * string which holds the message built in the constructor.
   */
  std::string what_message;
};

/**
 * A class representing an exception used only to send a program to
 * the terminate handler for abort after cleanup, while bypassing the
 * usual debugging output (performance logs, stack traces,
 * "terminating" messages) that the handler does to ease debugging of
 * uncaught error exceptions.
 *
 * We don't even inherit from std::exception here, to avoid being
 * caught as that type.
 */
class TerminationException
{
public:
  TerminationException() {}

  const char * what() const noexcept { return "libMesh termination requested"; }
};

}

#ifdef LIBMESH_ENABLE_EXCEPTIONS
#define libmesh_noexcept noexcept

#define LIBMESH_THROW(e) do { throw e; } while (0)
#define libmesh_rethrow throw
#define libmesh_try try
#define libmesh_catch(e) catch(e)

#else

#define LIBMESH_THROW(e) do { libMesh::err << e.what(); libMesh::libmesh_abort(); } while (0)
#define libmesh_rethrow
#define libmesh_try
#define libmesh_catch(e) if (0)

#endif // LIBMESH_ENABLE_EXCEPTIONS

#endif // LIBMESH_LIBMESH_EXCEPTIONS_H
