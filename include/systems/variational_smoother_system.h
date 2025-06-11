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

#ifndef LIBMESH_VARIATIONAL_SMOOTHER_SYSTEM_H
#define LIBMESH_VARIATIONAL_SMOOTHER_SYSTEM_H

// libMesh includes
#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/fem_system.h"
#include "libmesh/libmesh_common.h"

// C++ includes
#include <map>
#include <memory>

namespace libMesh
{

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class VariationalSmootherSystem : public libMesh::FEMSystem
{
/**
 * This is an FEMSystem to solve the optimization probelem posed by the
 * VariationalMeshSmoother class.
 *
 * The residual is coded as the gradient of the distortion-dilation metric, and
 * the jacobian as analytically coded as the Hessian of the metric.
 *
 * The nodes of the system mesh are updated during the solve.
 */
public:
  VariationalSmootherSystem(libMesh::EquationSystems & es,
                const std::string & name,
                const unsigned int number)
  : libMesh::FEMSystem(es, name, number),
    input_system(nullptr),
    _fe_family("LAGRANGE"),
    _fe_order(1),
    _epsilon_squared(1e-10),
    _ref_vol(1.),
    _dilation_weight(0.5)
  {}

  // Default destructor
  ~VariationalSmootherSystem() override;

  /**
   * Assembly method to update the mesh based on the smoother solve.
   */
  virtual void assembly (bool get_residual,
                         bool get_jacobian,
                         bool apply_heterogeneous_constraints = false,
                         bool apply_no_constraints = false) override;

  Real & get_ref_vol() { return _ref_vol; }
  Real & get_dilation_weight() { return _dilation_weight; }

  std::string & fe_family() { return _fe_family; }
  unsigned int & fe_order() { return _fe_order; }

  // We want to be able to project functions based on *other* systems'
  // values.  For that we need not only a FEMFunction but also a
  // reference to the system where it applies and a separate context
  // object (or multiple separate context objects, in the threaded
  // case) for that system.
  libMesh::System * input_system;

protected:
  std::map<libMesh::FEMContext *, std::unique_ptr<libMesh::FEMContext>>
    input_contexts;

  // System initialization
  virtual void init_data () override;

  // Context initialization
  virtual void init_context (libMesh::DiffContext & context) override;

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        libMesh::DiffContext & context) override;

  /// The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  /// The small nonzero constant to prevent zero denominators (degenerate elements only)
  const Real _epsilon_squared;

  /// The reference volume on an ideal element
  Real _ref_vol;

  /// The relative weight to give the dilation metric. The distortion metric is given weight 1 - _dilation_weight.
  Real _dilation_weight;
};

} // namespace libMesh

#endif // LIBMESH_VARIATIONAL_SMOOTHER_SYSTEM_H
