#ifndef LINEAR_ELASTICITY_WITH_CONTACT_H
#define LINEAR_ELASTICITY_WITH_CONTACT_H

#include "libmesh/dof_map.h"
#include "libmesh/nonlinear_implicit_system.h"

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
  NonlinearImplicitSystem & _sys;

  /**
   * Penalize overlapping elements.
   */
  Real _contact_penalty;

  /**
   * Store the intermediate values of lambda plus penalty. The dof IDs refer
   * to the nodes on the upper contact surface.
   */
  std::map<dof_id_type, Real> _lambda_plus_penalty_values;

  /**
   * Augmented Lagrangian values at each contact node. The dof IDs refer
   * to the nodes on the upper contact surface.
   */
  std::map<dof_id_type, Real> _lambdas;

  /**
   * This provides a map between contact nodes.
   */
  std::map<dof_id_type, dof_id_type> _contact_node_map;

public:

  /**
   * Constructor.
   */
  LinearElasticityWithContact(NonlinearImplicitSystem & sys_in,
                              Real contact_penalty_in);

  /**
   * Update the penalty parameter.
   */
  void set_contact_penalty(Real contact_penalty_in);

  /**
   * Get the penalty parameter.
   */
  Real get_contact_penalty() const;

  /**
   * Kronecker delta function.
   */
  Real kronecker_delta(unsigned int i,
                       unsigned int j);

  /**
   * Evaluate the fourth order tensor (C_ijkl) that relates stress to strain.
   */
  Real elasticity_tensor(Real young_modulus,
                         Real poisson_ratio,
                         unsigned int i,
                         unsigned int j,
                         unsigned int k,
                         unsigned int l);

  /**
   * Move the mesh nodes of input_mesh based on the displacement field
   * in input_solution.
   */
  void move_mesh(MeshBase & input_mesh,
                 const NumericVector<Number> & input_solution);

  /**
   * Set up the load paths on the contact surfaces.
   */
  void initialize_contact_load_paths();

  /**
   * Add edge elements into the mesh based on the contact load paths.
   * This ensure proper parallel communication of data, e.g. if a node
   * on one side of the contact surface has a constraint on it, then
   * adding contact elements into the mesh ensures that the constraint
   * will be enforced properly.
   */
  void add_contact_edge_elements();

  /**
   * Evaluate the Jacobian of the nonlinear system.
   */
  virtual void residual_and_jacobian (const NumericVector<Number> & soln,
                                      NumericVector<Number> * residual,
                                      SparseMatrix<Number> * jacobian,
                                      NonlinearImplicitSystem & /*sys*/);

  /**
   * Compute the Cauchy stress for the current solution.
   */
  void compute_stresses();

  /**
   * Update the lambda parameters in the augmented Lagrangian
   * method.
   * @return the largest change in the lambdas, and the largest
   * lambda value.
   */
  std::pair<Real, Real> update_lambdas();

  /**
   * @return the least and max gap function values for the current solution.
   */
  std::pair<Real, Real> get_least_and_max_gap_function();
};

#endif
