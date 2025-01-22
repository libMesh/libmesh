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

// <h1>Vector Finite Elements Example 8 - Single Face Hybridizable Discontinuous Galerkin</h1>
// \author Alexander Lindsay
// \date 2023
//
// This example is identical to Example 7 with the exception that the local
// solvers are built using a Single Face Discontinuous Galerkin
// formulation as opposed to a Raviart-Thomas Galerkin formulation. Equal order
// polynomials are used to discretize the scalar and vector fields. "Single
// Face" (SF) refers to the stabilization term added to the trace of the
// flux. Whereas in example 7, the trace of the flux \hat{u} is simply set to
// the flux u, in this SF-DG formulation \hat{u} = u + \tau * (p - \lambda)
// where \tau is only nonzero on a single element face

// Basic utilities.
#include "libmesh/string_to_enum.h"

// The solver packages supported by libMesh.
#include "libmesh/enum_solver_package.h"

// The mesh object and mesh generation and modification utilities.
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"

// Matrix and vector types.
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/numeric_vector.h"

// The finite element object and the geometric element type.
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"

// Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// The dof map, which handles degree of freedom indexing.
#include "libmesh/dof_map.h"

// The system of equations.
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"

// The exact solution and error computation.
#include "libmesh/exact_solution.h"
#include "libmesh/enum_norm_type.h"
#include "solution_function.h"

// I/O utilities.
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"

#ifdef LIBMESH_HAVE_EIGEN_DENSE

#include <Eigen/Dense>

using namespace libMesh;
using namespace Eigen;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
typedef MatrixXcd EigenMatrix;
typedef VectorXcd EigenVector;
#else
typedef MatrixXd EigenMatrix;
typedef VectorXd EigenVector;
#endif

void fe_assembly(EquationSystems & es, bool global_solve);
void assemble_hdg(EquationSystems & es, const std::string & system_name);
void alternative_fe_assembly(EquationSystems & es, bool global_solve);
void alternative_assemble_hdg(EquationSystems & es, const std::string & system_name);

int
main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init(argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Parse the input file.
  GetPot infile("vector_fe_ex8.in");

  // But allow the command line to override it.
  infile.parse_command_line(argc, argv);

  // Read in parameters from the command line and the input file.
  const unsigned int dimension = infile("dim", 2);
  const unsigned int grid_size = infile("grid_size", 15);
  const unsigned int order = infile("order", 1);
  const std::string family_str = infile("family", "L2_LAGRANGE");
  const FEFamily family = Utility::string_to_enum<FEFamily>(family_str);
  const std::string vec_family_str = infile("vecfamily", "L2_LAGRANGE_VEC");
  const FEFamily vec_family = Utility::string_to_enum<FEFamily>(vec_family_str);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build.
  libmesh_example_requires(dimension <= LIBMESH_DIM, dimension << "D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the cube [-1,1]^D. To accomodate first order side hierarchics, we must
  // use TRI6/7 or QUAD8/9 elements in 2d, or TET14 or HEX27 in 3d.
  const std::string elem_str = infile("element_type", std::string("TRI6"));

  libmesh_error_msg_if((dimension == 2 && elem_str != "TRI6" && elem_str != "TRI7" &&
                        elem_str != "QUAD8" && elem_str != "QUAD9") ||
                           (dimension == 3 && elem_str != "TET14" && elem_str != "HEX27"),
                       "You selected "
                           << elem_str
                           << " but this example must be run with TRI6, TRI7, QUAD8, or QUAD9 in 2d"
                           << " or with TET14, or HEX27 in 3d.");

  if (dimension == 2)
    MeshTools::Generation::build_square(
        mesh, grid_size, grid_size, -1., 1., -1., 1., Utility::string_to_enum<ElemType>(elem_str));
  else if (dimension == 3)
    MeshTools::Generation::build_cube(mesh,
                                      grid_size,
                                      grid_size,
                                      grid_size,
                                      -1.,
                                      1.,
                                      -1.,
                                      1.,
                                      -1.,
                                      1.,
                                      Utility::string_to_enum<ElemType>(elem_str));

  // Make sure the code is robust against nodal reorderings.
  MeshTools::Modification::permute_elements(mesh);

  // Create an equation systems object.
  EquationSystems equation_systems(mesh);

  // Declare the system "Mixed" and its variables.
  auto & system = equation_systems.add_system<System>("Mixed");

  // Add the LM system
  auto & lm_system = equation_systems.add_system<LinearImplicitSystem>("Lambda");

  // Adds the variable "q" and "u" to "Mixed". "q" will be our vector field
  // whereas "u" will be the scalar field.
  system.add_variable("q", static_cast<Order>(order), vec_family);
  system.add_variable("u", static_cast<Order>(order), family);

  // We also add a higher order version of our 'p' variable whose solution we
  // will postprocess using the Lagrange multiplier, 'u', and the low order 'p'
  // solution
  system.add_variable("u_enriched", static_cast<Order>(order+1), family);

  // Add our Lagrange multiplier to the implicit system
  lm_system.add_variable("lambda", static_cast<Order>(order), SIDE_HIERARCHIC);

  // Give the system a pointer to the matrix assembly
  // function. This will be called when needed by the library.
  lm_system.attach_assemble_function(assemble_hdg);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Solve the implicit system for the Lagrange multiplier
  lm_system.solve();

  // Armed with our Lagrange multiplier solution, we can now compute the vector and scalar solutions
  fe_assembly(equation_systems, /*global_solve=*/false);

  //
  // Now we will compute our solution approximation errors
  //

  ExactSolution exact_sol(equation_systems);

  if (dimension == 2)
  {
    SolutionFunction<2> soln_func;
    SolutionGradient<2> soln_grad;

    // Build FunctionBase* containers to attach to the ExactSolution object.
    std::vector<FunctionBase<Number> *> sols(1, &soln_func);
    std::vector<FunctionBase<Gradient> *> grads(1, &soln_grad);

    exact_sol.attach_exact_values(sols);
    exact_sol.attach_exact_derivs(grads);
  }
  else if (dimension == 3)
  {
    SolutionFunction<3> soln_func;
    SolutionGradient<3> soln_grad;

    // Build FunctionBase* containers to attach to the ExactSolution object.
    std::vector<FunctionBase<Number> *> sols(1, &soln_func);
    std::vector<FunctionBase<Gradient> *> grads(1, &soln_grad);

    exact_sol.attach_exact_values(sols);
    exact_sol.attach_exact_derivs(grads);
  }

  // Use higher quadrature order for more accurate error results.
  int extra_error_quadrature = infile("extra_error_quadrature", 2);
  if (extra_error_quadrature)
    exact_sol.extra_quadrature_order(extra_error_quadrature);

  // Compute the error.
  exact_sol.compute_error("Mixed", "q");
  exact_sol.compute_error("Mixed", "u");
  exact_sol.compute_error("Mixed", "u_enriched");

  // Print out the error values.
  libMesh::out << "L2 error is: " << exact_sol.l2_error("Mixed", "q") << std::endl;
  libMesh::out << "L2 error for u is: " << exact_sol.l2_error("Mixed", "u") << std::endl;
  libMesh::out << "L2 error u_enriched is: " << exact_sol.l2_error("Mixed", "u_enriched")
               << std::endl;

#ifdef LIBMESH_HAVE_EXODUS_API

  // We write the file in the ExodusII format.
  ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // Allright let's dry a different implementation and then ensure we get the same results
  // auto mixed_soln = system.solution->clone();
  // auto lm_soln = lm_system.solution->clone();
  lm_system.attach_assemble_function(alternative_assemble_hdg);
  lm_system.solve();
  // Armed with our Lagrange multiplier solution, we can now compute the vector and scalar solutions
  alternative_fe_assembly(equation_systems, /*global_solve=*/false);

  // Compute the error.
  exact_sol.compute_error("Mixed", "q");
  exact_sol.compute_error("Mixed", "u");
  exact_sol.compute_error("Mixed", "u_enriched");
  // Print out the error values.
  libMesh::out << "L2 error is: " << exact_sol.l2_error("Mixed", "q") << std::endl;
  libMesh::out << "L2 error for u is: " << exact_sol.l2_error("Mixed", "u") << std::endl;
  libMesh::out << "L2 error u_enriched is: " << exact_sol.l2_error("Mixed", "u_enriched")
               << std::endl;

  // All done.
  return 0;
}

// compute a solution indexable at quadrature points composed from the local degree of freedom
// solution vector and associated basis functions
template <typename SolnType, typename PhiType, typename LocalSolutionVector>
void
compute_qp_soln(std::vector<SolnType> & qp_vec,
                const unsigned int n_qps,
                const std::vector<std::vector<PhiType>> & phi,
                const LocalSolutionVector & soln)
{
  libmesh_assert(cast_int<std::size_t>(soln.size()) == phi.size());
  qp_vec.resize(n_qps);
  for (auto & val : qp_vec)
    val = {};
  for (const auto i : index_range(phi))
  {
    const auto & qp_phis = phi[i];
    libmesh_assert(qp_phis.size() == n_qps);
    const auto sol = soln(i);
    for (const auto qp : make_range(n_qps))
      qp_vec[qp] += qp_phis[qp] * sol;
  }
}

void
compute_enriched_soln(const MeshBase & mesh,
                      const DofMap & dof_map,
                      System & system,
                      const Elem * const elem,
                      const EigenVector & vector_soln,
                      const EigenVector & scalar_soln,
                      const EigenVector & Lambda,
                      FEVectorBase & vector_fe,
                      FEVectorBase & vector_fe_face,
                      FEBase & scalar_fe,
                      FEBase & scalar_fe_face,
                      FEBase & lambda_fe_face,
                      QBase & qrule,
                      QBase & qface)
{
  const unsigned int dim = mesh.mesh_dimension();

  // Create FE objects
  const FEType enriched_scalar_fe_type =
      dof_map.variable_type(system.variable_number("u_enriched"));
  std::unique_ptr<FEBase> enriched_scalar_fe(FEBase::build(dim, enriched_scalar_fe_type));
  std::unique_ptr<FEBase> enriched_scalar_fe_face(FEBase::build(dim, enriched_scalar_fe_type));

  // Attach quadrature rules
  enriched_scalar_fe->attach_quadrature_rule(&qrule);
  enriched_scalar_fe_face->attach_quadrature_rule(&qface);

  // pre-request our required volumetric data
  const auto & JxW = vector_fe.get_JxW();
  const auto & q_point = vector_fe.get_xyz();
  const auto & scalar_phi = scalar_fe.get_phi();
  const auto & enriched_scalar_phi = enriched_scalar_fe->get_phi();
  const auto & enriched_scalar_dphi = enriched_scalar_fe->get_dphi();

  // pre-request our required element face data
  const auto & vector_phi_face = vector_fe_face.get_phi();
  const auto & scalar_phi_face = scalar_fe_face.get_phi();
  const auto & enriched_scalar_phi_face = enriched_scalar_fe_face->get_phi();
  const auto & lambda_phi_face = lambda_fe_face.get_phi();
  const auto & JxW_face = scalar_fe_face.get_JxW();
  const auto & normals = vector_fe_face.get_normals();

  // Data structures for computing the enriched scalar solution
  DenseMatrix<Number> K_enriched_scalar;
  DenseVector<Number> F_enriched_scalar, U_enriched_scalar;
  std::vector<dof_id_type> enriched_scalar_dof_indices;

  // The lambda solution at the quadrature points, used for computing the enriched scalar solution
  std::vector<Number> lambda_qps;
  // The scalar solution at the quadrature points, used for computing the enriched scalar solution
  std::vector<Number> scalar_qps;
  // The vector solution at the quadrature points, used for computing the enriched scalar solution
  std::vector<Gradient> vector_qps;

  dof_map.dof_indices(elem, enriched_scalar_dof_indices, system.variable_number("u_enriched"));
  const auto enriched_scalar_n_dofs = enriched_scalar_dof_indices.size();
  // We have to add one for the mean value constraint
  const auto m = enriched_scalar_n_dofs + 1;
  const auto n = enriched_scalar_n_dofs + 1;

  K_enriched_scalar.resize(m, n);
  F_enriched_scalar.resize(m);
  U_enriched_scalar.resize(n);

  enriched_scalar_fe->reinit(elem);

  // We need the u solution for getting the correct mean value
  compute_qp_soln(scalar_qps, qrule.n_points(), scalar_phi, scalar_soln);

  //
  // We solve a modified diffusion problem
  //

  for (const auto qp : make_range(qrule.n_points()))
  {
    // Diffusion kernel
    for (const auto i : make_range(enriched_scalar_n_dofs))
      for (const auto j : make_range(enriched_scalar_n_dofs))
        K_enriched_scalar(i, j) +=
            JxW[qp] * (enriched_scalar_dphi[i][qp] * enriched_scalar_dphi[j][qp]);

    // Forcing function kernel
    {
      const Real x = q_point[qp](0);
      const Real y = q_point[qp](1);
      const Real z = q_point[qp](2);

      // "f" is the forcing function for the Poisson equation, which is
      // just the divergence of the exact solution for the vector field.
      // This is the well-known "method of manufactured solutions".
      Real f = 0;
      if (dim == 2)
        f = MixedExactSolution().forcing(x, y);
      else if (dim == 3)
        f = MixedExactSolution().forcing(x, y, z);

      // Scalar equation RHS
      for (const auto i : make_range(enriched_scalar_n_dofs))
        F_enriched_scalar(i) += JxW[qp] * enriched_scalar_phi[i][qp] * f;
    }

    // Mean value part
    {
      // u dependence on LM
      for (const auto i : make_range(enriched_scalar_n_dofs))
        K_enriched_scalar(i, enriched_scalar_n_dofs) += JxW[qp] * enriched_scalar_phi[i][qp];

      // LM dependence on u
      for (const auto j : make_range(enriched_scalar_n_dofs))
        K_enriched_scalar(enriched_scalar_n_dofs, j) += JxW[qp] * enriched_scalar_phi[j][qp];

      // And RHS of LM equation
      F_enriched_scalar(enriched_scalar_n_dofs) += JxW[qp] * scalar_qps[qp];
    }
  }

  bool tau_found = false;
  for (auto side : elem->side_index_range())
  {
    // Reinit our face FE objects
    vector_fe_face.reinit(elem, side);
    compute_qp_soln(vector_qps, qface.n_points(), vector_phi_face, vector_soln);

    enriched_scalar_fe_face->reinit(elem, side);

    if (elem->neighbor_ptr(side) == nullptr)
      for (const auto qp : make_range(qface.n_points()))
        // Now do the external boundary term for <\hat{q} \cdot \vec{n}, \omega> ->
        // <q + \tau (u - g), \omega> ->
        // <q, \omega> + <\tau u, \omega> - <\tau g, \omega>
        // BUT u = g on the boundary so we can drop those terms and just end up with
        // <q, \omega>
        for (const auto i : make_range(enriched_scalar_n_dofs))
          F_enriched_scalar(i) -=
              JxW_face[qp] * enriched_scalar_phi_face[i][qp] * vector_qps[qp] * normals[qp];
    else
    {
      scalar_fe_face.reinit(elem, side);
      compute_qp_soln(scalar_qps, qface.n_points(), scalar_phi_face, scalar_soln);
      lambda_fe_face.reinit(elem, side);
      compute_qp_soln(lambda_qps, qface.n_points(), lambda_phi_face, Lambda);

      // Stabilization parameter. In the single face discretization, only a single face has a
      // non-zero value of tau
      const Real tau = tau_found ? 0 : 1 / elem->hmin();
      tau_found = true;

      for (const auto qp : make_range(qface.n_points()))
      {
        const auto normal = normals[qp];
        const auto qhat = vector_qps[qp] + tau * (scalar_qps[qp] - lambda_qps[qp]) * normal;

        // Now do the internal boundary term for <\hat{q} \cdot \vec{n}, \omega> ->
        for (const auto i : make_range(enriched_scalar_n_dofs))
          F_enriched_scalar(i) -= JxW_face[qp] * enriched_scalar_phi_face[i][qp] * qhat * normal;
      }
    }
  }

  K_enriched_scalar.lu_solve(F_enriched_scalar, U_enriched_scalar);

  // Our solution for the local enriched scalar dofs is complete. Insert into the global vector
  // after eliminating the mean value constraint dof
  DenseVector<Number> U_insertion(enriched_scalar_n_dofs);
  for (const auto i : make_range(enriched_scalar_n_dofs))
    U_insertion(i) = U_enriched_scalar(i);
  system.solution->insert(U_insertion, enriched_scalar_dof_indices);
}

// We will perform finite element assembly twice. The first time to setup the global implicit system
// for the Lagrange multiplier degrees of freedom. And then the second time to compute the vector
// and scalar field solutions using the already-compute Lagrange multiplier solution
void
fe_assembly(EquationSystems & es, const bool global_solve)
{
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // The mixed, e.g. vector-scalar system
  auto & system = es.get_system<System>("Mixed");
  // Our implicit Lagrange multiplier system
  auto & lambda_system = es.get_system<LinearImplicitSystem>("Lambda");

  const auto & dof_map = system.get_dof_map();
  const auto & lambda_dof_map = lambda_system.get_dof_map();

  const FEType vector_fe_type = dof_map.variable_type(system.variable_number("q"));
  const FEType scalar_fe_type = dof_map.variable_type(system.variable_number("u"));
  const FEType lambda_fe_type =
      lambda_dof_map.variable_type(lambda_system.variable_number("lambda"));

  // Volumetric FE objects
  std::unique_ptr<FEVectorBase> vector_fe(FEVectorBase::build(dim, vector_fe_type));
  std::unique_ptr<FEBase> scalar_fe(FEBase::build(dim, scalar_fe_type));

  // Volumetric quadrature rule
  QGauss qrule(dim, scalar_fe_type.default_quadrature_order());

  // Attach quadrature rules for the FE objects that we will reinit within the element "volume"
  vector_fe->attach_quadrature_rule(&qrule);
  scalar_fe->attach_quadrature_rule(&qrule);

  // Declare finite element objects for boundary integration
  std::unique_ptr<FEVectorBase> vector_fe_face(FEVectorBase::build(dim, vector_fe_type));
  std::unique_ptr<FEBase> scalar_fe_face(FEBase::build(dim, scalar_fe_type));
  std::unique_ptr<FEBase> lambda_fe_face(FEBase::build(dim, lambda_fe_type));

  // Boundary integration requires one quadrature rule with dimensionality one
  // less than the dimensionality of the element.
  QGauss qface(dim - 1, scalar_fe_type.default_quadrature_order());

  // Attach quadrature rules for the FE objects that we will reinit on the element faces
  vector_fe_face->attach_quadrature_rule(&qface);
  scalar_fe_face->attach_quadrature_rule(&qface);
  lambda_fe_face->attach_quadrature_rule(&qface);

  // pre-request our required volumetric data
  const auto & JxW = vector_fe->get_JxW();
  const auto & q_point = vector_fe->get_xyz();
  const auto & vector_phi = vector_fe->get_phi();
  const auto & scalar_phi = scalar_fe->get_phi();
  const auto & grad_scalar_phi = scalar_fe->get_dphi();
  const auto & div_vector_phi = vector_fe->get_div_phi();

  // pre-request our required element face data
  const auto & vector_phi_face = vector_fe_face->get_phi();
  const auto & scalar_phi_face = scalar_fe_face->get_phi();
  const auto & lambda_phi_face = lambda_fe_face->get_phi();
  const auto & JxW_face = scalar_fe_face->get_JxW();
  const auto & qface_point = vector_fe_face->get_xyz();
  const auto & normals = vector_fe_face->get_normals();

  //
  // We will need "Eigen" versions of many of the matrices/vectors because the
  // libMesh DenseMatrix doesn't have an inverse API
  //

  // LM matrix and RHS after eliminating vector and scalar dofs
  DenseMatrix<Number> K_libmesh;
  DenseVector<Number> P_libmesh;
  EigenMatrix K;
  EigenVector P;
  // (A  B  C)(q) = (G)
  // (Bt D  E)(u) = (F)
  // (Ct Et H)(l) = (0)
  // Sinv is the inverse of a Schur complement S which is given by
  // S = Bt * A^{-1} * B.
  EigenMatrix A, Ainv, B, Bt, C, Ct, Sinv, D, E, Et, H;
  EigenVector G, F;

  // Lambda eigen vector for constructing vector and scalar solutions
  EigenVector Lambda;

  // Containers for dof indices
  std::vector<dof_id_type> vector_dof_indices;
  std::vector<dof_id_type> scalar_dof_indices;
  std::vector<dof_id_type> lambda_dof_indices;
  std::vector<Number> lambda_solution_std_vec;

  // The global system matrix
  SparseMatrix<Number> & matrix = lambda_system.get_system_matrix();

  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Retrive our dof indices for all fields
    dof_map.dof_indices(elem, vector_dof_indices, system.variable_number("q"));
    dof_map.dof_indices(elem, scalar_dof_indices, system.variable_number("u"));
    lambda_dof_map.dof_indices(elem, lambda_dof_indices, lambda_system.variable_number("lambda"));

    const auto vector_n_dofs = vector_dof_indices.size();
    const auto scalar_n_dofs = scalar_dof_indices.size();
    const auto lambda_n_dofs = lambda_dof_indices.size();

    // Reinit our volume FE objects
    vector_fe->reinit(elem);
    scalar_fe->reinit(elem);

    libmesh_assert_equal_to(vector_n_dofs, vector_phi.size());
    libmesh_assert_equal_to(scalar_n_dofs, scalar_phi.size());

    // prepare our matrix/vector data structures for the vector equation
    A.setZero(vector_n_dofs, vector_n_dofs);
    B.setZero(vector_n_dofs, scalar_n_dofs);
    C.setZero(vector_n_dofs, lambda_n_dofs);
    G.setZero(vector_n_dofs);
    // and for the scalar equation
    Bt.setZero(scalar_n_dofs, vector_n_dofs);
    D.setZero(scalar_n_dofs, scalar_n_dofs);
    E.setZero(scalar_n_dofs, lambda_n_dofs);
    F.setZero(scalar_n_dofs);
    // and for the LM equation
    Ct.setZero(lambda_n_dofs, vector_n_dofs);
    Et.setZero(lambda_n_dofs, scalar_n_dofs);
    H.setZero(lambda_n_dofs, lambda_n_dofs);

    for (const auto qp : make_range(qrule.n_points()))
    {
      // Vector equation dependence on vector dofs
      for (const auto i : make_range(vector_n_dofs))
        for (const auto j : make_range(vector_n_dofs))
          A(i, j) -= JxW[qp] * (vector_phi[i][qp] * vector_phi[j][qp]);

      // Vector equation dependence on scalar dofs
      for (const auto i : make_range(vector_n_dofs))
        for (const auto j : make_range(scalar_n_dofs))
          B(i, j) += JxW[qp] * (div_vector_phi[i][qp] * scalar_phi[j][qp]);

      // Scalar equation dependence on vector dofs
      for (const auto i : make_range(scalar_n_dofs))
        for (const auto j : make_range(vector_n_dofs))
          Bt(i, j) += JxW[qp] * (grad_scalar_phi[i][qp] * vector_phi[j][qp]);

      // This is the end of the matrix summation loop
      // Now we build the element right-hand-side contribution.
      // This involves a single loop in which we integrate the "forcing
      // function" in the PDE against the scalar test functions (k).
      {
        const Real x = q_point[qp](0);
        const Real y = q_point[qp](1);
        const Real z = q_point[qp](2);

        // "f" is the forcing function for the Poisson equation, which is
        // just the divergence of the exact solution for the vector field.
        // This is the well-known "method of manufactured solutions".
        Real f = 0;
        if (dim == 2)
          f = MixedExactSolution().forcing(x, y);
        else if (dim == 3)
          f = MixedExactSolution().forcing(x, y, z);

        // Scalar equation RHS
        for (const auto i : make_range(scalar_n_dofs))
          F(i) -= JxW[qp] * scalar_phi[i][qp] * f;
      }
    }

    // At the beginning of the loop, we mark that we haven't found our "Single-Face" yet
    bool tau_found = false;
    for (auto side : elem->side_index_range())
    {
      // Reinit our face FE objects
      vector_fe_face->reinit(elem, side);
      scalar_fe_face->reinit(elem, side);
      lambda_fe_face->reinit(elem, side);

      if (elem->neighbor_ptr(side) == nullptr)
        for (const auto qp : make_range(qface.n_points()))
        {
          const Real xf = qface_point[qp](0);
          const Real yf = qface_point[qp](1);
          const Real zf = qface_point[qp](2);

          // The boundary value for scalar field.
          Real scalar_value = 0;
          if (dim == 2)
            scalar_value = MixedExactSolution().scalar(xf, yf);
          else if (dim == 3)
            scalar_value = MixedExactSolution().scalar(xf, yf, zf);

          // External boundary -> Dirichlet faces -> Vector equation RHS
          for (const auto i : make_range(vector_n_dofs))
            G(i) += JxW_face[qp] * (vector_phi_face[i][qp] * normals[qp]) * scalar_value;

          // Now do the external boundary term for <\hat{q} \cdot \vec{n}, \omega> ->
          // <q + \tau (u - g), \omega> ->
          // <q, \omega> + <\tau u, \omega> - <\tau g, \omega>
          // BUT u = g on the boundary so we can drop those terms and just end up with
          // <q, \omega>
          for (const auto i : make_range(scalar_n_dofs))
            for (const auto j : make_range(vector_n_dofs))
              Bt(i, j) -=
                  JxW_face[qp] * scalar_phi_face[i][qp] * (vector_phi_face[j][qp] * normals[qp]);

          // Need to do something with the external boundary LM dofs to prevent
          // the matrix from being singular. Use a minus sign because the
          // formulation by Cockburn results in negative diagonals for the LM
          for (const auto i : make_range(lambda_n_dofs))
            for (const auto j : make_range(lambda_n_dofs))
              H(i, j) += JxW_face[qp] * lambda_phi_face[i][qp] * lambda_phi_face[j][qp];
        }
      else
      {
        // If we haven't found our "Single-Face" yet, then we assign tau to
        // something non-zero. Else we have previously designated our nonzero
        // "Single-Face" and so we mark tau as zero
        const Real tau = tau_found ? 0 : 1 / elem->hmin();
        tau_found = true;

        for (const auto qp : make_range(qface.n_points()))
        {
          const auto normal = normals[qp];
          const auto normal_sq = normal * normal;

          // Vector equation dependence on LM dofs
          for (const auto i : make_range(vector_n_dofs))
            for (const auto j : make_range(lambda_n_dofs))
              C(i, j) -=
                  JxW_face[qp] * (vector_phi_face[i][qp] * normals[qp]) * lambda_phi_face[j][qp];

          // Now do the internal boundary term for <\hat{q} \cdot \vec{n}, \omega> ->
          // <q + \tau (u - \lambda), \omega> ->
          // <q, \omega> + <\tau u, \omega> - <\tau \lambda, \omega>
          for (const auto i : make_range(scalar_n_dofs))
          {
            for (const auto j : make_range(vector_n_dofs))
              Bt(i, j) -= JxW_face[qp] * scalar_phi_face[i][qp] * (vector_phi_face[j][qp] * normal);

            for (const auto j : make_range(scalar_n_dofs))
              D(i, j) -=
                  JxW_face[qp] * scalar_phi_face[i][qp] * tau * scalar_phi_face[j][qp] * normal_sq;

            for (const auto j : make_range(lambda_n_dofs))
              E(i, j) +=
                  JxW_face[qp] * scalar_phi_face[i][qp] * tau * lambda_phi_face[j][qp] * normal_sq;
          }

          // Now do the internal boundary term for <\hat{q} \cdot \vec{n}, \mu> ->
          // <q + \tau (u - \lambda), \mu> ->
          // <q, \mu> + <\tau u, \mu> - <\tau \lambda, \mu>
          for (const auto i : make_range(lambda_n_dofs))
          {
            for (const auto j : make_range(vector_n_dofs))
              Ct(i, j) -= JxW_face[qp] * lambda_phi_face[i][qp] * (vector_phi_face[j][qp] * normal);

            for (const auto j : make_range(scalar_n_dofs))
              Et(i, j) -=
                  JxW_face[qp] * lambda_phi_face[i][qp] * tau * scalar_phi_face[j][qp] * normal_sq;

            for (const auto j : make_range(lambda_n_dofs))
              H(i, j) +=
                  JxW_face[qp] * lambda_phi_face[i][qp] * tau * lambda_phi_face[j][qp] * normal_sq;
          }
        }
      }
    }

    Ainv = A.inverse();
    // Compute the Schur complement inverse
    Sinv = (D - Bt * Ainv * B).inverse();
    const auto BtAinvCMinusE = Bt * Ainv * C - E;
    const auto BtAinvGMinusF = Bt * Ainv * G - F;
    const auto CtAinv = Ct * Ainv;
    const auto BSinv = B * Sinv;
    const auto EtSinv = Et * Sinv;
    K = -CtAinv * (BSinv * BtAinvCMinusE + C) + EtSinv * BtAinvCMinusE + H;
    P = -CtAinv * (BSinv * BtAinvGMinusF + G) + EtSinv * BtAinvGMinusF;

    // Build our libMesh data structures from the Eigen ones
    K_libmesh.resize(K.rows(), K.cols());
    P_libmesh.resize(P.size());
    libmesh_assert((K.rows() == K.cols()) && (K.rows() == P.size()));
    for (const auto i : make_range(K.rows()))
    {
      P_libmesh(i) = P(i);
      for (const auto j : make_range(K.cols()))
        K_libmesh(i, j) = K(i, j);
    }

    if (global_solve)
    {
      // We were performing our finite element assembly for the implicit solve step of our
      // example. Add our local element vectors/matrices into the global system
      dof_map.constrain_element_matrix_and_vector(K_libmesh, P_libmesh, lambda_dof_indices);
      matrix.add_matrix(K_libmesh, lambda_dof_indices);
      lambda_system.rhs->add_vector(P_libmesh, lambda_dof_indices);
    }
    else
    {
      //
      // We are doing our finite element assembly for the second time. We now know the Lagrange
      // multiplier solution. With that and the local element matrices and vectors we can compute
      // the vector and scalar solutions
      //

      Lambda.resize(lambda_n_dofs);
      lambda_system.current_local_solution->get(lambda_dof_indices, lambda_solution_std_vec);
      for (const auto i : make_range(lambda_n_dofs))
        Lambda(i) = lambda_solution_std_vec[i];
      const auto scalar_soln = Sinv * (BtAinvCMinusE * Lambda - BtAinvGMinusF);
      const auto vector_soln = Ainv * (G - B * scalar_soln - C * Lambda);
      for (const auto i : make_range(vector_n_dofs))
        system.solution->set(vector_dof_indices[i], vector_soln(i));
      for (const auto i : make_range(scalar_n_dofs))
        system.solution->set(scalar_dof_indices[i], scalar_soln(i));

      // Now solve for the enriched scalar solution using our Lagrange
      // multiplier solution, q, and our low-order u
      compute_enriched_soln(mesh,
                            dof_map,
                            system,
                            elem,
                            vector_soln,
                            scalar_soln,
                            Lambda,
                            *vector_fe,
                            *vector_fe_face,
                            *scalar_fe,
                            *scalar_fe_face,
                            *lambda_fe_face,
                            qrule,
                            qface);
    }
  }

  if (!global_solve)
  {
    system.solution->close();
    // Scatter solution into the current_solution which is used in error computation
    system.update();
  }
}

// Call this assembly function when assembling the global implicit system
void
assemble_hdg(EquationSystems & es, const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "Lambda");
  fe_assembly(es, /*global_solve=*/true);
}

// Compute the stabilization parameter
Real
compute_tau(const bool internal_face, bool & tau_found, const Elem * const elem)
{
  if (!internal_face)
    // if we're an external face then tau is 0
    return 0;
  else if (tau_found)
    // if we've already applied a non-zero tau on another face, then we also return 0
    return 0;
  else
  {
    // This is our first internal face. We will apply our non-zero tau (and mark that we have)
    tau_found = true;
    return 1 / elem->hmin();
  }
}

// We will perform finite element assembly twice. The first time to setup the global implicit system
// for the Lagrange multiplier degrees of freedom. And then the second time to compute the vector
// and scalar field solutions using the already-compute Lagrange multiplier solution
void
alternative_fe_assembly(EquationSystems & es, const bool global_solve)
{
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // The mixed, e.g. vector-scalar system
  auto & system = es.get_system<System>("Mixed");
  // Our implicit Lagrange multiplier system
  auto & lambda_system = es.get_system<LinearImplicitSystem>("Lambda");

  const auto & dof_map = system.get_dof_map();
  const auto & lambda_dof_map = lambda_system.get_dof_map();

  const FEType vector_fe_type = dof_map.variable_type(system.variable_number("q"));
  const FEType scalar_fe_type = dof_map.variable_type(system.variable_number("u"));
  const FEType lambda_fe_type =
      lambda_dof_map.variable_type(lambda_system.variable_number("lambda"));

  // Volumetric FE objects
  std::unique_ptr<FEVectorBase> vector_fe(FEVectorBase::build(dim, vector_fe_type));
  std::unique_ptr<FEBase> scalar_fe(FEBase::build(dim, scalar_fe_type));

  // Volumetric quadrature rule
  QGauss qrule(dim, scalar_fe_type.default_quadrature_order());

  // Attach quadrature rules for the FE objects that we will reinit within the element "volume"
  vector_fe->attach_quadrature_rule(&qrule);
  scalar_fe->attach_quadrature_rule(&qrule);

  // Declare finite element objects for boundary integration
  std::unique_ptr<FEVectorBase> vector_fe_face(FEVectorBase::build(dim, vector_fe_type));
  std::unique_ptr<FEBase> scalar_fe_face(FEBase::build(dim, scalar_fe_type));
  std::unique_ptr<FEBase> lambda_fe_face(FEBase::build(dim, lambda_fe_type));

  // Boundary integration requires one quadrature rule with dimensionality one
  // less than the dimensionality of the element.
  QGauss qface(dim - 1, scalar_fe_type.default_quadrature_order());

  // Attach quadrature rules for the FE objects that we will reinit on the element faces
  vector_fe_face->attach_quadrature_rule(&qface);
  scalar_fe_face->attach_quadrature_rule(&qface);
  lambda_fe_face->attach_quadrature_rule(&qface);

  // pre-request our required volumetric data
  const auto & JxW = vector_fe->get_JxW();
  const auto & q_point = vector_fe->get_xyz();
  const auto & vector_phi = vector_fe->get_phi();
  const auto & scalar_phi = scalar_fe->get_phi();
  const auto & grad_scalar_phi = scalar_fe->get_dphi();
  const auto & div_vector_phi = vector_fe->get_div_phi();

  // pre-request our required element face data
  const auto & vector_phi_face = vector_fe_face->get_phi();
  const auto & scalar_phi_face = scalar_fe_face->get_phi();
  const auto & lambda_phi_face = lambda_fe_face->get_phi();
  const auto & JxW_face = scalar_fe_face->get_JxW();
  const auto & qface_point = vector_fe_face->get_xyz();
  const auto & normals = vector_fe_face->get_normals();

  //
  // We will need "Eigen" versions of many of the matrices/vectors because the
  // libMesh DenseMatrix doesn't have an inverse API
  //

  // LM matrix and RHS after eliminating vector and scalar dofs
  DenseMatrix<Number> K_lm_libmesh;
  DenseVector<Number> F_lm_libmesh;
  EigenMatrix K_mixed;
  EigenMatrix Kinv_mixed;
  EigenVector F_mixed;

  // Lambda eigen vector for constructing vector and scalar solutions
  EigenVector Lambda;

  // Containers for dof indices
  std::vector<dof_id_type> vector_dof_indices;
  std::vector<dof_id_type> scalar_dof_indices;
  std::vector<dof_id_type> lambda_dof_indices;
  std::vector<Number> lambda_solution_std_vec;

  // Helper container for storing a given "mu" (function in the Langrange multiplier space) with
  // quadrature point evaluations
  std::vector<Number> mu;

  // The global system matrix
  auto & matrix = lambda_system.get_system_matrix();

  auto compute_and_invert_K =
      [&](const auto vector_n_dofs_in, const auto scalar_n_dofs_in, const Elem * const elem_in)
  {
    const auto mixed_size = vector_n_dofs_in + scalar_n_dofs_in;

    K_mixed.setZero(mixed_size, mixed_size);

    for (const auto qp : make_range(qrule.n_points()))
    {
      // Vector equation dependence on vector dofs
      for (const auto i : make_range(vector_n_dofs_in))
        for (const auto j : make_range(vector_n_dofs_in))
          K_mixed(i, j) += JxW[qp] * (vector_phi[i][qp] * vector_phi[j][qp]);

      // Vector equation dependence on scalar dofs
      for (const auto i : make_range(vector_n_dofs_in))
        for (const auto j : make_range(scalar_n_dofs_in))
          K_mixed(i, j + vector_n_dofs_in) -= JxW[qp] * (div_vector_phi[i][qp] * scalar_phi[j][qp]);

      // Scalar equation dependence on vector dofs
      for (const auto i : make_range(scalar_n_dofs_in))
        for (const auto j : make_range(vector_n_dofs_in))
          K_mixed(i + vector_n_dofs_in, j) -= JxW[qp] * (grad_scalar_phi[i][qp] * vector_phi[j][qp]);
    }

    // At the beginning of the loop, we mark that we haven't found our "Single-Face" yet
    bool tau_found = false;
    for (auto side : elem_in->side_index_range())
    {
      // Reinit our face FE objects
      vector_fe_face->reinit(elem_in, side);
      scalar_fe_face->reinit(elem_in, side);
      const bool internal_face = elem_in->neighbor_ptr(side);
      const auto tau = compute_tau(internal_face, tau_found, elem_in);

      for (const auto qp : make_range(qface.n_points()))
      {
        const auto normal = normals[qp];
        const auto normal_sq = normal * normal;

        // Now do the internal boundary term for <\hat{q} \cdot \vec{n}, \omega> ->
        // <q + \tau (u - \lambda), \omega> ->
        // <q, \omega> + <\tau u, \omega> - <\tau \lambda, \omega>
        for (const auto i : make_range(scalar_n_dofs_in))
        {
          for (const auto j : make_range(vector_n_dofs_in))
            K_mixed(i + vector_n_dofs_in, j) +=
                JxW_face[qp] * scalar_phi_face[i][qp] * (vector_phi_face[j][qp] * normal);

          if (tau) // Don't do unnecessary math ops if tau is 0
            for (const auto j : make_range(scalar_n_dofs_in))
              K_mixed(i + vector_n_dofs_in, j + vector_n_dofs_in) +=
                  JxW_face[qp] * scalar_phi_face[i][qp] * tau * scalar_phi_face[j][qp] * normal_sq;
        }
      }
    }

    Kinv_mixed = K_mixed.inverse();
  };

  auto compute_rhs = [&](const auto vector_n_dofs_in,
                         const auto scalar_n_dofs_in,
                         const Elem * const elem_in,
                         const unsigned int shape_function)
  {
    const auto mixed_size = vector_n_dofs_in + scalar_n_dofs_in;
    F_mixed.setZero(mixed_size);

    // If the approximate LM solution was passed in, then we are solving for the elemental solution
    // of q and u, not the mappings of individual LM shape functions. Consequently, we include the
    // contribution of the forcing function
    if (shape_function == libMesh::invalid_uint)
      for (const auto qp : make_range(qrule.n_points()))
      {
        const Real x = q_point[qp](0);
        const Real y = q_point[qp](1);
        const Real z = q_point[qp](2);

        Real f = 0;
        if (dim == 2)
          f = MixedExactSolution().forcing(x, y);
        else if (dim == 3)
          f = MixedExactSolution().forcing(x, y, z);
        for (const auto ii : make_range(scalar_n_dofs_in))
          F_mixed(ii + vector_n_dofs_in) += JxW[qp] * f * scalar_phi[ii][qp];
      }

    // At the beginning of the loop, we mark that we haven't found our "Single-Face" yet
    bool tau_found = false;
    std::vector<Number> g;

    for (auto side : elem_in->side_index_range())
    {
      // Reinit our face FE objects
      vector_fe_face->reinit(elem_in, side);
      scalar_fe_face->reinit(elem_in, side);
      lambda_fe_face->reinit(elem_in, side);

      const auto & qp_mu = [&]()
      {
        if (shape_function == libMesh::invalid_uint)
        {
          if (elem_in->neighbor_ptr(side))
          {
            compute_qp_soln(lambda_solution_std_vec, qface.n_points(), lambda_phi_face, Lambda);
            return lambda_solution_std_vec;
          }
          else
          {
            g.resize(qface.n_points());
            for (const auto qp : make_range(qface.n_points()))
            {
              auto & g_qp = g[qp];
              const Real xf = qface_point[qp](0);
              const Real yf = qface_point[qp](1);
              const Real zf = qface_point[qp](2);

              // The boundary value for scalar field.
              if (dim == 2)
                g_qp = MixedExactSolution().scalar(xf, yf);
              else if (dim == 3)
                g_qp = MixedExactSolution().scalar(xf, yf, zf);
            }
            return g;
          }
        }
        else
        {
          auto & real_mu = lambda_phi_face[shape_function];
          mu.resize(qface.n_points());
          for (const auto qp : make_range(qface.n_points()))
            mu[qp] = real_mu[qp];
          return mu;
        }
      }();

      const bool internal_face = elem_in->neighbor_ptr(side);
      const auto tau = compute_tau(internal_face, tau_found, elem_in);

      for (const auto qp : make_range(qface.n_points()))
      {
        const auto normal = normals[qp];
        const auto normal_sq = normal * normal;

        // Vector equation dependence on LM/mu
        for (const auto ii : make_range(vector_n_dofs_in))
          F_mixed(ii) -= JxW_face[qp] * (vector_phi_face[ii][qp] * normal) * qp_mu[qp];

        // Now do the boundary term for <\hat{q} \cdot \vec{n}, \omega> ->
        // <q + \tau (u - \lambda), \omega> ->
        // <q, \omega> + <\tau u, \omega> - <\tau \lambda, \omega>
        for (const auto ii : make_range(scalar_n_dofs_in))
          F_mixed(ii + vector_n_dofs_in) +=
              JxW_face[qp] * scalar_phi_face[ii][qp] * tau * qp_mu[qp] * normal_sq;
      }
    }
  };

  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    std::vector<EigenVector> local_solns;
    std::vector<unsigned int> dofs_on_side;
    std::unordered_set<unsigned int> external_boundary_indices;
    std::vector<std::vector<Gradient>> volumetric_q;
    std::vector<std::vector<Number>> volumetric_u;
    std::vector<std::vector<Gradient>> face_q;

    // Retrive our dof indices for all fields
    dof_map.dof_indices(elem, vector_dof_indices, system.variable_number("q"));
    dof_map.dof_indices(elem, scalar_dof_indices, system.variable_number("u"));
    lambda_dof_map.dof_indices(elem, lambda_dof_indices, lambda_system.variable_number("lambda"));

    const auto vector_n_dofs = vector_dof_indices.size();
    const auto scalar_n_dofs = scalar_dof_indices.size();
    const auto lambda_n_dofs = lambda_dof_indices.size();
    const std::size_t n_mu_funcs = global_solve ? lambda_n_dofs : 1;

    if (global_solve)
    {
      K_lm_libmesh.resize(lambda_n_dofs, lambda_n_dofs);
      F_lm_libmesh.resize(lambda_n_dofs);

      for (const auto s : elem->side_index_range())
        if (!elem->neighbor_ptr(s))
        {
          FEInterface::dofs_on_side(elem, dim, lambda_fe_type, s, dofs_on_side);
          external_boundary_indices.insert(dofs_on_side.begin(), dofs_on_side.end());
        }
    }

    // Reinit our volume FE objects
    vector_fe->reinit(elem);
    scalar_fe->reinit(elem);

    libmesh_assert_equal_to(vector_n_dofs, vector_phi.size());
    libmesh_assert_equal_to(scalar_n_dofs, scalar_phi.size());

    compute_and_invert_K(vector_n_dofs, scalar_n_dofs, elem);
    local_solns.resize(n_mu_funcs);

    if (global_solve)
    {
      volumetric_q.resize(lambda_n_dofs);
      volumetric_u.resize(lambda_n_dofs);
      face_q.resize(lambda_n_dofs);

      for (const auto i : make_range(lambda_n_dofs))
      {
        if (external_boundary_indices.count(i))
          continue;

        compute_rhs(vector_n_dofs, scalar_n_dofs, elem, i);
        auto & local_soln = local_solns[i];
        local_soln = Kinv_mixed * F_mixed;
        const auto local_q_soln = local_soln.head(vector_n_dofs);
        const auto local_u_soln = local_soln.tail(scalar_n_dofs);
        compute_qp_soln(volumetric_q[i], qrule.n_points(), vector_phi, local_q_soln);
        compute_qp_soln(volumetric_u[i], qrule.n_points(), scalar_phi, local_u_soln);
      }

      // Create the bilinear form for lambda
      for (const auto i : make_range(lambda_n_dofs))
        if (!external_boundary_indices.count(i))
          for (const auto j : make_range(lambda_n_dofs))
            if (!external_boundary_indices.count(j))
              for (const auto qp : make_range(qrule.n_points()))
                K_lm_libmesh(i, j) += JxW[qp] * volumetric_q[i][qp] * volumetric_q[j][qp];

      // Now for the volumetric portion of the RHS
      for (const auto qp : make_range(qrule.n_points()))
      {
        const Real x = q_point[qp](0);
        const Real y = q_point[qp](1);
        const Real z = q_point[qp](2);

        // "f" is the forcing function for the Poisson equation, which is
        // just the divergence of the exact solution for the vector field.
        // This is the well-known "method of manufactured solutions".
        Real f = 0;
        if (dim == 2)
          f = MixedExactSolution().forcing(x, y);
        else if (dim == 3)
          f = MixedExactSolution().forcing(x, y, z);
        for (const auto i : make_range(lambda_n_dofs))
          if (!external_boundary_indices.count(i))
            F_lm_libmesh(i) += JxW[qp] * f * volumetric_u[i][qp];
      }

      // Now for the Dirichlet boundary portion of our RHS
      for (const auto s : elem->side_index_range())
      {
        lambda_fe_face->reinit(elem, s);
        for (const auto i : make_range(lambda_n_dofs))
          if (external_boundary_indices.count(i))
            for (const auto j : make_range(lambda_n_dofs))
              if (external_boundary_indices.count(j))
                for (const auto qp : make_range(qface.n_points()))
                  K_lm_libmesh(i, j) +=
                      JxW_face[qp] * lambda_phi_face[i][qp] * lambda_phi_face[j][qp];

        if (!elem->neighbor_ptr(s))
        {
          vector_fe_face->reinit(elem, s);
          for (const auto i : make_range(lambda_n_dofs))
          {
            if (external_boundary_indices.count(i))
              continue;

            const auto local_q_soln = local_solns[i].head(vector_n_dofs);
            compute_qp_soln(face_q[i], qface.n_points(), vector_phi_face, local_q_soln);
          }

          for (const auto qp : make_range(qface.n_points()))
          {
            const Real xf = qface_point[qp](0);
            const Real yf = qface_point[qp](1);
            const Real zf = qface_point[qp](2);

            // The boundary value for scalar field.
            Real scalar_value = 0;
            if (dim == 2)
              scalar_value = MixedExactSolution().scalar(xf, yf);
            else if (dim == 3)
              scalar_value = MixedExactSolution().scalar(xf, yf, zf);
            for (const auto i : make_range(lambda_n_dofs))
              if (!external_boundary_indices.count(i))
                F_lm_libmesh(i) += JxW_face[qp] * scalar_value * face_q[i][qp] * normals[qp];
          }
        }
      }

      // We were performing our finite element assembly for the implicit solve step of our
      // example. Add our local element vectors/matrices into the global system
      dof_map.constrain_element_matrix_and_vector(K_lm_libmesh, F_lm_libmesh, lambda_dof_indices);
      matrix.add_matrix(K_lm_libmesh, lambda_dof_indices);
      lambda_system.rhs->add_vector(F_lm_libmesh, lambda_dof_indices);
    }
    else
    {
      // Must populate Lambda before calling compute_rhs
      Lambda.resize(lambda_n_dofs);
      lambda_system.current_local_solution->get(lambda_dof_indices, lambda_solution_std_vec);
      for (const auto i : make_range(lambda_n_dofs))
        Lambda(i) = lambda_solution_std_vec[i];

      compute_rhs(vector_n_dofs, scalar_n_dofs, elem, libMesh::invalid_uint);
      auto & local_soln = local_solns[0];
      local_soln = Kinv_mixed * F_mixed;
      const auto vector_soln = local_soln.head(vector_n_dofs);
      const auto scalar_soln = local_soln.tail(scalar_n_dofs);
      for (const auto i : make_range(vector_n_dofs))
        system.solution->set(vector_dof_indices[i], vector_soln(i));
      for (const auto i : make_range(scalar_n_dofs))
        system.solution->set(scalar_dof_indices[i], scalar_soln(i));

      // Now solve for the enriched scalar solution using our Lagrange
      // multiplier solution, q, and our low-order u
      compute_enriched_soln(mesh,
                            dof_map,
                            system,
                            elem,
                            vector_soln,
                            scalar_soln,
                            Lambda,
                            *vector_fe,
                            *vector_fe_face,
                            *scalar_fe,
                            *scalar_fe_face,
                            *lambda_fe_face,
                            qrule,
                            qface);
    }
  }

  if (!global_solve)
  {
    system.solution->close();
    // Scatter solution into the current_solution which is used in error computation
    system.update();
  }
}

// Call this assembly function when assembling the global implicit system
void
alternative_assemble_hdg(EquationSystems & es, const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "Lambda");
  alternative_fe_assembly(es, /*global_solve=*/true);
}

#else

int
main()
{
  return 0;
}

#endif
