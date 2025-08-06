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

/**
 * Struct to hold smoother-relevant information about the mesh quality.
 */
struct MeshQualityInfo
{
  // dof_id_types will hold the id of the Elem where the metric occurs

  std::pair<dof_id_type, Real> max_elem_distortion{DofObject::invalid_id,
                                                   std::numeric_limits<Real>::lowest()};
  std::pair<dof_id_type, Real> min_elem_distortion{DofObject::invalid_id,
                                                   std::numeric_limits<Real>::max()};
  Real total_distortion = 0.;

  std::pair<dof_id_type, Real> max_elem_dilation{DofObject::invalid_id,
                                                 std::numeric_limits<Real>::lowest()};
  std::pair<dof_id_type, Real> min_elem_dilation{DofObject::invalid_id,
                                                 std::numeric_limits<Real>::max()};
  Real total_dilation = 0.;

  std::pair<dof_id_type, Real> max_elem_combined{DofObject::invalid_id,
                                                 std::numeric_limits<Real>::lowest()};
  std::pair<dof_id_type, Real> min_elem_combined{DofObject::invalid_id,
                                                 std::numeric_limits<Real>::max()};
  Real total_combined = 0.;

  std::pair<dof_id_type, Real> max_elem_det_S{DofObject::invalid_id,
                                              std::numeric_limits<Real>::lowest()};
  std::pair<dof_id_type, Real> min_elem_det_S{DofObject::invalid_id,
                                              std::numeric_limits<Real>::max()};
  Real total_det_S = 0.;
  Real max_qp_det_S = std::numeric_limits<Real>::lowest();
  Real min_qp_det_S = std::numeric_limits<Real>::max();

  bool mesh_is_tangled = false;
  bool initialized = false;
};

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
      _epsilon_squared(TOLERANCE),
      _ref_vol(0.),
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

  Real & get_dilation_weight() { return _dilation_weight; }

  /**
   * Get the target element for a given element type.
   * @param type Element type
   * @return a std::pair containing the target element for type and the
   * corresponding nodes that must be kept in scope while the target element is
   * used.
   */
  static std::pair<std::unique_ptr<Elem>, std::vector<std::unique_ptr<Node>>>
  get_target_elem(const ElemType & type);

  /**
   * Get the jacobians (and determinants) of the target-to-reference element mapping.
   * @param target_elem Target element.
   * @param femcontext Context used to build mapping.
   * @param jacobian Vector in which to store the jacobians for each quadrature point.
   * @param jacobian_dets Vector in which to store the determinant of the jacobians
   * for each quadrature point.
   */
  static void get_target_to_reference_jacobian(const Elem * const target_elem,
                                               const FEMContext & femcontext,
                                               std::vector<RealTensor> & jacobians,
                                               std::vector<Real> & jacobian_dets);

  /**
   * Getter for the _mesh_info attribute. If this attribute has not yet been
   * initialized, compute_mesh_quality_info is called to initialize it.
   */
  const MeshQualityInfo & get_mesh_info();

  /*
   * Computes information about the mesh quality and sets the _mesh_info attribute.
   */
  void compute_mesh_quality_info();

protected:

  // System initialization
  virtual void init_data () override;

  // Context initialization
  virtual void init_context (libMesh::DiffContext & context) override;

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        libMesh::DiffContext & context) override;

  /* Computes the element reference volume used in the dilation metric
   * The reference value is set to the averaged value of all elements' average
   * |J|. Also computes any applicable target element inverse Jacobians. Target
   * elements are relavant when the reference element does not minimize the
   * distortion metric.
   */
  void prepare_for_smoothing();

  /**
  * The small nonzero constant to prevent zero denominators (degenerate meshes only)
  */
  const Real _epsilon_squared;

  /**
   * Epsilon squared value determined at runtime during each assembly. The value
   * depends on whether the mesh is tangled.
   */
  Real _epsilon_squared_assembly;

  /**
  * The reference volume for each element
  */
  Real _ref_vol;

  /**
  * The relative weight to give the dilation metric. The distortion metric is given weight 1 - _dilation_weight.
  */
  Real _dilation_weight;

  /* Map to hold target qp-dependent element target-to-reference mapping
   * Jacobians, if any
   */
  std::map<ElemType, std::vector<RealTensor>> _target_jacobians;

  /*
   * Map to hold the determinants of _target_jacobians.
   */
  std::map<ElemType, std::vector<Real>> _target_jacobian_dets;

  /**
   * Information about the mesh quality.
   */
  MeshQualityInfo _mesh_info;
};

} // namespace libMesh

#endif // LIBMESH_VARIATIONAL_SMOOTHER_SYSTEM_H
