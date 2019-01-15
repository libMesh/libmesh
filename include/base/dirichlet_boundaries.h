// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DIRICHLET_BOUNDARIES_H
#define LIBMESH_DIRICHLET_BOUNDARIES_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_DIRICHLET

// Local Includes
#include "libmesh/auto_ptr.h" // deprecated
#include "libmesh/id_types.h"
#include "libmesh/vector_value.h"

// C++ Includes
#include <cstddef>
#include <set>
#include <vector>
#include <memory>

namespace libMesh
{

// Forward declarations
class System;

template <typename Output> class FEMFunctionBase;
template <typename Output> class FunctionBase;

/**
 * Dirichlet functions may be indexed either by "system variable
 * order" or "local variable order", depending on how the
 * DirichletBoundary object is constructed.  For example, suppose a
 * system has variables {a, b, c, d}, and a DirichletBoundary is set
 * for variables {b, d} (i.e. variables_in is {1, 3}).  If the
 * boundary is constructed to use "system variable order", input
 * function(s) will be queried for components 1 and 3; this is useful
 * for reusing input functions as both exact solutions and Dirichlet
 * boundaries in benchmark problems.  If the boundary is constructed
 * to use "local variable order", input function(s) will be queried
 * for components 0 and 1; this is useful for flexibly constructing
 * Dirichlet boundaries in multiphysics codes or from user input
 * files.
 */
enum VariableIndexing { SYSTEM_VARIABLE_ORDER = 0,
                        LOCAL_VARIABLE_ORDER };

/**
 * This class allows one to associate Dirichlet boundary values with
 * a given set of mesh boundary ids and system variable ids.
 *
 * Dirichlet values must be supplied as the input function "f"; when
 * using some specialized elements, gradient values must be supplied
 * via the input function "g".
 *
 * Dirichlet functions may be subclasses of FunctionBase or
 * FEMFunctionBase; in the latter case the user must also supply a
 * reference to the System on which the FEMFunctionBase will be
 * evaluated.
 *
 * Dirichlet functions are allowed to return NaN; if this is
 * encountered, then the degree of freedom values in a patch around
 * the location of the returned NaN will be left unconstrained.  E.g.
 * a NaN on a boundary edge in 3D would leave that edge and the two
 * adjoining face interiors unconstrained, but would still permit the
 * other edge and node DoFs around those faces to be constrained.
 *
 * \author Roy Stogner
 * \date 2012
 * \brief Class for specifying Dirichlet boundary conditions as constraints.
 */
class DirichletBoundary
{
public:

  /**
   * Constructor for a system-variable-order boundary using
   * pointers-to-functors.
   */
  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const FunctionBase<Number> * f_in,
                    const FunctionBase<Gradient> * g_in = nullptr);

  /**
   * Constructor for a boundary from reference-to-functor.
   *
   * Defaults to system variable indexing for backwards compatibility,
   * but most users will prefer local indexing.
   */
  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const FunctionBase<Number> & f_in,
                    VariableIndexing type = SYSTEM_VARIABLE_ORDER);

  /**
   * Constructor for a system-variable-order boundary from
   * references-to-functors.
   *
   * Defaults to system variable indexing for backwards compatibility,
   * but most users will prefer local indexing.
   */
  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const FunctionBase<Number> & f_in,
                    const FunctionBase<Gradient> & g_in,
                    VariableIndexing type = SYSTEM_VARIABLE_ORDER);

  /**
   * Constructor for a system-variable-order boundary from
   * pointers-to-fem-functors.
   */
  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const System & f_sys_in,
                    const FEMFunctionBase<Number> * f_in,
                    const FEMFunctionBase<Gradient> * g_in = nullptr);

  /**
   * Constructor for a system-variable-order boundary from
   * reference-to-fem-functor.
   *
   * Defaults to system variable indexing for backwards compatibility,
   * but most users will prefer local indexing.
   */
  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const System & f_sys_in,
                    const FEMFunctionBase<Number> & f_in,
                    VariableIndexing type = SYSTEM_VARIABLE_ORDER);

  /**
   * Constructor for a system-variable-order boundary from
   * references-to-fem-functors.
   *
   * Defaults to system variable indexing for backwards compatibility,
   * but most users will prefer local indexing.
   */
  DirichletBoundary(const std::set<boundary_id_type> & b_in,
                    const std::vector<unsigned int> & variables_in,
                    const System & f_sys_in,
                    const FEMFunctionBase<Number> & f_in,
                    const FEMFunctionBase<Gradient> & g_in,
                    VariableIndexing type = SYSTEM_VARIABLE_ORDER);

  /**
   * Copy constructor.  Deep copies (clones) functors; shallow copies
   * any System reference
   */
  DirichletBoundary (const DirichletBoundary & dirichlet_in);

  /**
   * Standard destructor
   */
  ~DirichletBoundary ();

  std::set<boundary_id_type> b;
  std::vector<unsigned int> variables;

  std::unique_ptr<FunctionBase<Number>> f;
  std::unique_ptr<FunctionBase<Gradient>> g;

  std::unique_ptr<FEMFunctionBase<Number>> f_fem;
  std::unique_ptr<FEMFunctionBase<Gradient>> g_fem;

  const System * f_system;
};


/**
 * We're using a class instead of a typedef to allow forward
 * declarations and future flexibility.
 *
 * \note \p std::map has no virtual destructor, so downcasting here
 * would be dangerous.
 */
class DirichletBoundaries : public std::vector<DirichletBoundary *>
{
public:
  DirichletBoundaries() {}

  ~DirichletBoundaries();
};

} // namespace libMesh

#endif // LIBMESH_ENABLE_DIRICHLET

#endif // LIBMESH_DIRICHLET_BOUNDARIES_H
