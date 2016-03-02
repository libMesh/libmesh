// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SOLVER_CONFIGURATION_H
#define LIBMESH_SOLVER_CONFIGURATION_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"

namespace libMesh
{

/**
 * This class stores solver configuration data, e.g. tolerances
 * and iteration limits. Override the pure virtual function
 * configure_solver() in order to provide specific behavior for
 * particular types of solvers (e.g. LinearSolver, EigenSolver, etc).
 *
 * \author David Knezevic
 * \date 2015
 */
class SolverConfiguration : public ReferenceCountedObject<SolverConfiguration>
{
public:

  /**
   *  Constructor.
   */
  SolverConfiguration () {}

  /**
   * Destructor.
   */
  virtual ~SolverConfiguration () {}

  /**
   * Apply options during initialization of a solver.
   * Default is a no-op. Override in subclasses to provide
   * specific behavior.
   */
  virtual void set_options_during_init() {}

  /**
   * Apply solver options to a particular solver.
   * Override in subclasses to provide specific behavior.
   */
  virtual void configure_solver() = 0;

  /**
   * This method can be called after the solver has failed (e.g. on
   * catching an exception thrown by a solver). It allows an appropriate
   * response to the failure, e.g. terminate the program, or change
   * configuation options and try again.
   * \p solve_failure_count specifies the number of times we've tried to
   * recover from a failure.
   */
  virtual void respond_to_solve_failure(unsigned int /*solve_failure_count*/) {}

  /**
   * Store real-valued solver parameters in this map, e.g. solver tolerances.
   */
  std::map<std::string, Real> real_valued_data;

  /**
   * Store integer solver parameters in this map, e.g. iteration limits.
   */
  std::map<std::string, int> int_valued_data;

  /**
   * Store string data in this map, e.g. "solver_type" : "gmres".
   */
  std::map<std::string, std::string> string_data;

};

} // namespace libMesh

#endif // LIBMESH_SOLVER_CONFIGURATION_H
