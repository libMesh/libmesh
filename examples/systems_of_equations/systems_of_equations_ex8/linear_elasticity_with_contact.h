#ifndef LINEAR_ELASTICITY_WITH_CONTACT_H
#define LINEAR_ELASTICITY_WITH_CONTACT_H

#include "libmesh/dof_map.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "augment_sparsity_on_contact.h"

// define the subdomain IDs
#define TOP_SUBDOMAIN 2
#define BOTTOM_SUBDOMAIN 1

// define the boundary IDs in the mesh
#define MIN_Z_BOUNDARY 1
#define MAX_Z_BOUNDARY 2
#define CONTACT_BOUNDARY_LOWER 3
#define CONTACT_BOUNDARY_UPPER 4

using libMesh::DofMap;
using libMesh::NonlinearImplicitSystem;
using libMesh::dof_id_type;
using libMesh::Point;
using libMesh::Real;
using libMesh::Number;
using libMesh::MeshBase;
using libMesh::NumericVector;
using libMesh::SparseMatrix;

/**
 * Convenient class to encapsulate (element,side,qp) data.
 */
class QuadraturePointOnSideId
{
public:

  /**
   * Need a default constructor so that we can use this class in STL containers.
   */
  QuadraturePointOnSideId()
  :
    _element_id(0),
    _side_index(0),
    _qp(0)
  {}

  /**
   * Constructor to set the members.
   */
  QuadraturePointOnSideId(
    dof_id_type element_id,
    unsigned char side_index,
    unsigned int qp)
  :
    _element_id(element_id),
    _side_index(side_index),
    _qp(qp)
  {
  }

  /**
   * The element ID for this quadrature point.
   */
  dof_id_type _element_id;

  /**
   * The side of the element that this quadrature point belongs to.
   */
  unsigned char _side_index;

  /**
   * The quadrature point index on the side.
   */
  unsigned int _qp;
};

/**
 * Need operator< for QuadraturePointOnSideId so that we can use it as
 * the key in a map.
 */
bool operator< (QuadraturePointOnSideId const& a, QuadraturePointOnSideId const& b);

/**
 * Convenient class to encapsulate data related to the
 * line that defines the distance between elements in contact.
 */
class IntersectionPointData
{
public:

  /**
   * Need a default constructor so that we can use this class in STL containers.
   */
  IntersectionPointData()
  :
    _neighbor_element_id(0),
    _neighbor_side_index(0)
  {}

  /**
   * Constructor to set the members.
   */
  IntersectionPointData(
    dof_id_type neighbor_element_id,
    unsigned char neighbor_side_index,
    Point intersection_point,
    Point inverse_mapped_intersection_point,
    Point line_point,
    Point line_direction)
  :
    _neighbor_element_id(neighbor_element_id),
    _neighbor_side_index(neighbor_side_index),
    _intersection_point(intersection_point),
    _inverse_mapped_intersection_point(inverse_mapped_intersection_point),
    _line_point(line_point),
    _line_direction(line_direction)
  {
  }

  /**
   * The element ID for this intersection point.
   */
  dof_id_type _neighbor_element_id;

  /**
   * The side of the element that this intersection point belongs to.
   */
  unsigned char _neighbor_side_index;

  /**
   * The intersection point.
   */
  Point _intersection_point;

  /**
   * The intersection point mapped back to the reference element.
   */
  Point _inverse_mapped_intersection_point;

  /**
   * The point on the "master" element.
   */
  Point _line_point;

  /**
   * The normal on the "master" element, which
   * defines the direction of the line
   */
  Point _line_direction;

};

/**
 * This class encapsulate all functionality required for assembling
 * and solving a linear elastic model with contact.
 */
class LinearElasticityWithContact :
  public NonlinearImplicitSystem::ComputeResidualandJacobian
{
private:

  /**
   * Keep a reference to the NonlinearImplicitSystem.
   */
  NonlinearImplicitSystem& _sys;

  /**
   * Penalize overlapping elements.
   */
  Real _contact_penalty;

  /**
   * Tolerance relative to element size for detecting nearby elements
   * which might be in contact.
   */
  Real _contact_proximity_tol;

  /**
   * Tolerance relative to element size for detecting if a point
   * is inside an element.
   */
  Real _contains_point_tol;

  /**
   * This maps from (element ID,side,qp) to the corresponding IntersectionPointData.
   */
  std::map<QuadraturePointOnSideId, IntersectionPointData> _contact_intersection_data;

  /**
   * The object that handles augmenting the sparsity pattern.
   */
  AugmentSparsityOnContact _augment_sparsity;

public:

  /**
   * Constructor.
   */
  LinearElasticityWithContact(
    NonlinearImplicitSystem &sys_in,
    Real contact_penalty_in,
    Real contact_proximity_tol);

  /**
   * @return a reference to the object for augmenting the sparsity pattern.
   */
  AugmentSparsityOnContact& get_augment_sparsity();

  /**
   * Update the penalty parameter.
   */
  void set_contact_penalty(Real contact_penalty_in);

  /**
   * Get the penalty parameter.
   */
  Real get_contact_penalty() const;

  /**
   * Clear the non-zero contact forces from a previous iteration.
   */
  void clear_contact_data();

  /**
   * Set non-zero contact force for the specified element, side, quadrature point.
   * Also store the "nearest point" on the other element.
   */
  void set_contact_data(
    dof_id_type element_id,
    unsigned char side_index,
    unsigned int qp,
    IntersectionPointData intersection_pt_data);

  /**
   * @return true is we have detected contact for the specified element,side,qp.
   */
  bool is_contact_detected(
    dof_id_type element_id,
    unsigned char side_index,
    unsigned int qp);

  /**
   * Get an entry from the contact intersection point data structure.
   * Throw an error if the entry doesn't exist.
   */
  IntersectionPointData get_contact_data(
    dof_id_type element_id,
    unsigned char side_index,
    unsigned int qp);

  /**
   * Kronecker delta function.
   */
  Real kronecker_delta(
    unsigned int i,
    unsigned int j);

  /**
   * Evaluate the fourth order tensor (C_ijkl) that relates stress to strain.
   */
  Real elasticity_tensor(
    Real young_modulus,
    Real poisson_ratio,
    unsigned int i,
    unsigned int j,
    unsigned int k,
    unsigned int l);

  /**
   * Move the mesh nodes of \p input_mesh based on the displacement field
   * in \p input_solution.
   */
  void move_mesh(
    MeshBase& input_mesh,
    const NumericVector<Number>& input_solution);

  /**
   * Evaluate the Jacobian of the nonlinear system.
   */
  virtual void residual_and_jacobian (
    const NumericVector<Number>& soln,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian,
    NonlinearImplicitSystem& /*sys*/);

  /**
   * Compute the Cauchy stress for the current solution.
   */
  void compute_stresses();

};

#endif
