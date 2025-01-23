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


// <h1>Vector Finite Elements Example 7 - Hybridized Raviart-Thomas elements</h1>
// \author Alexander Lindsay
// \date 2023
//
// This example hybridizes Raviart-Thomas elements to solve a model div-grad
// problem in both 2d and 3d. Before hybridization, the mixed problem is simply
// a div-grad formulation, \vec{u} = -\nabla p, and \nabla \cdot \vec{u} = f, of
// the Poisson problem in Introduction Example 3, \nabla^2 p = -f. A standard
// (non-hybridized) Raviart-Thomas discretization of the same problem can be
// found in Vector Example 6.
//
// One of the first references for the hybridized Raviart-Thomas (RT)
// formulation is Arnold's and Brezzi's "Mixed and nonconforming finite element
// methods: implementation, postprocessing and error estimates". The key piece
// to the hybridized formulation is breaking the continuity of the RT space,
// e.g. we no longer require that the normal component of the RT space be
// continuous across interelement boundaries. In the notation of Arnold and
// Brezzi this means a change from RT_0^k -> RT_{-1}^k where k is the polynomial
// order of the basis (note that here k denotes the polynomial order of the
// *divergence* of the vector shape functions, so for k = 0, the RT basis in
// libMesh is FIRST). This breakage in continuity is accomplished by changing
// from an FEFamily of RAVIART_THOMAS to L2_RAVIART_THOMAS. Instead of being
// enforced by the finite element space itself, continuity of the normal
// component of the vector field is enforced through Lagrange multipliers that
// live on the sides of the elements. Let's introduce a subspace of L2, M_{-1}^k,
// where the functions in M_{-1}^k are polynomials of degree of k or less. These
// polynomials can live on our elements (I_h) or on our internal element faces
// (E0_h); we will denote these polynomial subspaces respectively as M_{-1}^k (I_h)
//  and M_{-1}^k (E0_h). We will denote boundary faces by E. The hybridized
// problem can be summarized as follows:
//
// find the approximate solutions (u_h, p_h, lambda_h) in
// RT_{-1}^k x M_{-1}^k (I_h) x M^{-1}^k (E0_h)
// such that
//
// (u, tau) - (p, div(tau)) + <lambda, tau*n>_E0 = -<g, tau*n>_E for all tau in RT_{-1}^k (I_h)
// -(v, div(u)) = -(f, v) for all v in M_{-1}^k (I_h)
// <mu, u*n> = 0 for all mu in M_{-1}^k (E0_h)
//
// with p = g on the boundary (Dirichlet problem). In the above, n denotes the
// normal vector facing outward from the element for which we are doing local
// assembly
//
// We can write this in a matrix form:
// (A  B  C)(u)        (G)
// (Bt 0  0)(p)      = (F)
// (Ct 0  0)(lambda)   (0)
//
// The matrix A is block diagonal, e.g. it is entirely localized within an
// element due to the breakage of the continuity requirements of the RT
// space. This allows us to write:
//
// u = A^{-1} (G - Bp - C lambda)
//
// We can further eliminate p to end up with a global system that depends only on lambda:
//
// E lambda = H
//
// where
//
// E = Ct * A^{-1} * (A - B * S^{-1} * Bt) * A^{-1} * C
// S = Bt * A^{-1} * B
// H = Hg + Hf
// Hg = Ct * A^{-1} * (A - B * S^{-1} * Bt) * A^{-1} * G
// Hf = Ct * A^{-1} * B * S^{-1} * F
//
// Here in our example we compose local element matrices A, B, C, Bt, Ct, G, and
// F using finite element assembly. Then due to the small size of the local
// element matrices, we actually compute the required inverses and build E and H
// which is then fed into the global matrix and vector respectively. The
// resulting global matrix is symmetric positive definite which is a nice change
// over the saddle point standard RT discretization. Moreover, the global system
// size of the hybridized problem is less than standard RT.
//
// Once the global system is solved. We go through finite element assembly a
// second time to locally construct the u and p solutions from lambda. Finally,
// lambda (and p for non-simplices) are used to reconstruct a higher-order
// approximation of p, e.g. a solution using a basis of polynomials of degree k + 1.
// Lambda is used in this reconstruction because it represents a
// projection of the true solution for p onto M_{-1}^k (E0_h) ; e.g. as is so
// often the case, the Lagrange multipliers have a physical meaning

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
void assemble_divgrad(EquationSystems & es, const std::string & system_name);

int
main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init(argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Parse the input file.
  GetPot infile("vector_fe_ex7.in");

  // But allow the command line to override it.
  infile.parse_command_line(argc, argv);

  // Read in parameters from the command line and the input file.
  const unsigned int dimension = infile("dim", 2);
  const unsigned int grid_size = infile("grid_size", 15);

  // Skip higher-dimensional examples on a lower-dimensional libMesh build.
  libmesh_example_requires(dimension <= LIBMESH_DIM, dimension << "D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the cube [-1,1]^D. To accomodate Raviart-Thomas elements, we must
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
  const bool simplicial = (elem_str == "TRI6") || (elem_str == "TRI7") || (elem_str == "TET14");
  equation_systems.parameters.set<bool>("simplicial") = simplicial;


  // Declare the system  "DivGrad" and its variables.
  auto & system = equation_systems.add_system<System>("DivGrad");

  // Add the LM system
  auto & lm_system = equation_systems.add_system<LinearImplicitSystem>("Lambda");

  // Adds the variable "u" and "p" to "DivGrad". "u" will be our vector field
  // whereas "p" will be the scalar field.
  system.add_variable("u", FIRST, L2_RAVIART_THOMAS);
  system.add_variable("p", CONSTANT, MONOMIAL);
  // We also add a higher order version of our 'p' variable whose solution we
  // will compute using the Lagrange multiplier field and, for non-simplexes,
  // the low order 'p' solution
  system.add_variable("p_enriched", FIRST, MONOMIAL);

  // Add our Lagrange multiplier to the implicit system
  lm_system.add_variable("lambda", CONSTANT, SIDE_HIERARCHIC);

  // Give the system a pointer to the matrix assembly
  // function. This will be called when needed by the library.
  lm_system.attach_assemble_function(assemble_divgrad);

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
  exact_sol.compute_error("DivGrad", "u");
  exact_sol.compute_error("DivGrad", "p");
#if !defined(LIBMESH_HAVE_PETSC) || !defined(LIBMESH_USE_REAL_NUMBERS)
  if (simplicial)
#endif
    exact_sol.compute_error("DivGrad", "p_enriched");

  // Print out the error values.
  libMesh::out << "L2 error is: " << exact_sol.l2_error("DivGrad", "u") << std::endl;
  libMesh::out << "HDiv semi-norm error is: " << exact_sol.error_norm("DivGrad", "u", HDIV_SEMINORM)
               << std::endl;
  libMesh::out << "HDiv error is: " << exact_sol.hdiv_error("DivGrad", "u") << std::endl;
  libMesh::out << "L2 error for p is: " << exact_sol.l2_error("DivGrad", "p") << std::endl;
#if !defined(LIBMESH_HAVE_PETSC) || !defined(LIBMESH_USE_REAL_NUMBERS)
  if (simplicial)
#endif
    libMesh::out << "L2 error p_enriched is: " << exact_sol.l2_error("DivGrad", "p_enriched")
                 << std::endl;

#ifdef LIBMESH_HAVE_EXODUS_API

  // We write the file in the ExodusII format.
  ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}

// We will perform finite element assembly twice. The first time to setup the global implicit system
// for the Lagrange multiplier degrees of freedom. And then the second time to compute the vector
// and scalar field solutions using the already-compute Lagrange multiplier solution
void
fe_assembly(EquationSystems & es, const bool global_solve)
{
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  // Are our elements simplicial?
  const bool simplicial = es.parameters.get<bool>("simplicial");

  // The div-grad, e.g. vector-scalar system
  auto & system = es.get_system<System>("DivGrad");
  // Our implicit Lagrange multiplier system
  auto & lambda_system = es.get_system<LinearImplicitSystem>("Lambda");

  const auto & dof_map = system.get_dof_map();
  const auto & lambda_dof_map = lambda_system.get_dof_map();

  const FEType vector_fe_type = dof_map.variable_type(system.variable_number("u"));
  const FEType scalar_fe_type = dof_map.variable_type(system.variable_number("p"));
  const FEType enriched_scalar_fe_type =
      dof_map.variable_type(system.variable_number("p_enriched"));
  const FEType lambda_fe_type =
      lambda_dof_map.variable_type(lambda_system.variable_number("lambda"));

  // Volumetric FE objects
  std::unique_ptr<FEVectorBase> vector_fe(FEVectorBase::build(dim, vector_fe_type));
  std::unique_ptr<FEBase> scalar_fe(FEBase::build(dim, scalar_fe_type));
  std::unique_ptr<FEBase> enriched_scalar_fe(FEBase::build(dim, enriched_scalar_fe_type));

  // Volumetric quadrature rule
  QGauss qrule(dim, FIFTH);

  // Attach quadrature rules for the FE objects that we will reinit within the element "volume"
  vector_fe->attach_quadrature_rule(&qrule);
  scalar_fe->attach_quadrature_rule(&qrule);
  enriched_scalar_fe->attach_quadrature_rule(&qrule);

  // Declare finite element objects for boundary integration
  std::unique_ptr<FEVectorBase> vector_fe_face(FEVectorBase::build(dim, vector_fe_type));
  std::unique_ptr<FEBase> enriched_scalar_fe_face(FEBase::build(dim, enriched_scalar_fe_type));
  std::unique_ptr<FEBase> lambda_fe_face(FEBase::build(dim, lambda_fe_type));

  // Boundary integration requires one quadrature rule with dimensionality one
  // less than the dimensionality of the element.
  QGauss qface(dim - 1, FIFTH);

  // Attach quadrature rules for the FE objects that we will reinit on the element faces
  vector_fe_face->attach_quadrature_rule(&qface);
  enriched_scalar_fe_face->attach_quadrature_rule(&qface);
  lambda_fe_face->attach_quadrature_rule(&qface);

  // pre-request our required volumetric data
  const auto & JxW = vector_fe->get_JxW();
  const auto & q_point = vector_fe->get_xyz();
  const auto & vector_phi = vector_fe->get_phi();
  const auto & scalar_phi = scalar_fe->get_phi();
  const auto & enriched_scalar_phi = enriched_scalar_fe->get_phi();
  const auto & div_vector_phi = vector_fe->get_div_phi();

  // pre-request our required element face data
  const auto & vector_phi_face = vector_fe_face->get_phi();
  const auto & enriched_scalar_phi_face = enriched_scalar_fe_face->get_phi();
  const auto & lambda_phi_face = lambda_fe_face->get_phi();
  const auto & JxW_face = vector_fe_face->get_JxW();
  const auto & qface_point = vector_fe_face->get_xyz();
  const auto & normals = vector_fe_face->get_normals();

  //
  // We follow the notation of Cockburn in "A Charactrization of Hybridized
  // Mixed methods for Second Order Elliptic problems"for the element
  // matrices/vectors. We will need "Eigen" versions of many of the
  // matrices/vectors because the libMesh DenseMatrix doesn't have an inverse
  // API
  //

  // LM matrix and RHS after eliminating vector and scalar dofs
  DenseMatrix<Number> E_libmesh;
  DenseVector<Number> H_libmesh;
  EigenMatrix E;
  EigenMatrix H;
  // Auxiliary matrices and RHS. A is vector-vector. B is vector-scalar. C is
  // vector-LM. Sinv is the inverse of the Schur complement S which is given by
  // S = Bt * A^{-1} * B. G is the RHS of the vector equation (e.g. the
  // Dirichlet condition). F is the RHS of the scalar equation (resulting from
  // the body force or MMS force in this case). Hg is the post-elimination
  // version of G, and Hf is the post-elimination version of F
  EigenMatrix A, Ainv, B, Bt, C, Ct, Sinv;
  EigenVector G, F, Hg, Hf;
  // element matrix for boundary LM dofs. We have to have this because our
  // SIDE_HIERARCHIC variables live on *all* element faces, not just interior
  // ones.
  EigenMatrix L;

  // Lambda eigen vector for constructing vector and scalar solutions
  EigenVector Lambda;
  // The lambda solution at the quadrature points
  std::vector<Number> lambda_qps;
  /// The scalar solution at the quadrature points
  std::vector<Number> scalar_qps;

  /// Data structures for computing the enriched scalar solution
  DenseMatrix<Number> K_enriched_scalar;
  DenseVector<Number> F_enriched_scalar, U_enriched_scalar;

  // Containers for dof indices
  std::vector<dof_id_type> vector_dof_indices;
  std::vector<dof_id_type> scalar_dof_indices;
  std::vector<dof_id_type> enriched_scalar_dof_indices;
  std::vector<dof_id_type> lambda_dof_indices;
  std::vector<Number> lambda_solution_std_vec;

  // The global system matrix
  SparseMatrix<Number> & matrix = lambda_system.get_system_matrix();

  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Retrive our dof indices for all fields
    dof_map.dof_indices(elem, vector_dof_indices, system.variable_number("u"));
    dof_map.dof_indices(elem, scalar_dof_indices, system.variable_number("p"));
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
    F.setZero(scalar_n_dofs);
    // and for the LM equation
    Ct.setZero(lambda_n_dofs, vector_n_dofs);
    L.setZero(lambda_n_dofs, lambda_n_dofs);

    for (const auto qp : make_range(qrule.n_points()))
    {
      // Vector equation dependence on vector dofs
      for (const auto i : make_range(vector_n_dofs))
        for (const auto j : make_range(vector_n_dofs))
          A(i, j) += JxW[qp] * (vector_phi[i][qp] * vector_phi[j][qp]);

      // Vector equation dependence on scalar dofs
      for (const auto i : make_range(vector_n_dofs))
        for (const auto j : make_range(scalar_n_dofs))
          B(i, j) -= JxW[qp] * (div_vector_phi[i][qp] * scalar_phi[j][qp]);

      // Scalar equation dependence on vector dofs
      for (const auto i : make_range(scalar_n_dofs))
        for (const auto j : make_range(vector_n_dofs))
          Bt(i, j) -= JxW[qp] * (scalar_phi[i][qp] * div_vector_phi[j][qp]);

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
          f = DivGradExactSolution().forcing(x, y);
        else if (dim == 3)
          f = DivGradExactSolution().forcing(x, y, z);

        // Scalar equation RHS
        for (const auto i : make_range(scalar_n_dofs))
          F(i) -= JxW[qp] * scalar_phi[i][qp] * f;
      }
    }

    {
      for (auto side : elem->side_index_range())
      {
        // Reinit our face FE objects
        vector_fe_face->reinit(elem, side);
        libmesh_assert_equal_to(vector_n_dofs, vector_phi_face.size());
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
              scalar_value = DivGradExactSolution().scalar(xf, yf);
            else if (dim == 3)
              scalar_value = DivGradExactSolution().scalar(xf, yf, zf);

            // External boundary -> Dirichlet faces -> Vector equaton RHS
            for (const auto i : make_range(vector_n_dofs))
              G(i) -= JxW_face[qp] * (vector_phi_face[i][qp] * normals[qp]) * scalar_value;

            // Need to do something with the external boundary LM dofs to prevent the matrix from
            // being singular
            for (const auto i : make_range(lambda_n_dofs))
              for (const auto j : make_range(lambda_n_dofs))
                L(i, j) += JxW_face[qp] * lambda_phi_face[i][qp] * lambda_phi_face[j][qp];
          }
        else
          for (const auto qp : make_range(qface.n_points()))
          {
            // Vector equation dependence on LM dofs
            for (const auto i : make_range(vector_n_dofs))
              for (const auto j : make_range(lambda_n_dofs))
                C(i, j) +=
                    JxW_face[qp] * (vector_phi_face[i][qp] * normals[qp]) * lambda_phi_face[j][qp];

            // LM equation dependence on vector dofs
            for (const auto i : make_range(lambda_n_dofs))
              for (const auto j : make_range(vector_n_dofs))
                Ct(i, j) +=
                    JxW_face[qp] * lambda_phi_face[i][qp] * (vector_phi_face[j][qp] * normals[qp]);
          }
      }
    }

    Ainv = A.inverse();
    // Compute the Schur complement inverse
    Sinv = (Bt * Ainv * B).inverse();
    // These equations can be derived at by eliminating the u and p dofs from the system
    // (A  B  C)(u)        (G)
    // (Bt 0  0)(p)      = (F)
    // (Ct 0  0)(lambda)   (0)
    E = Ct * Ainv * (A - B * Sinv * Bt) * Ainv * C;
    E = E + L;
    Hg = Ct * Ainv * (A - B * Sinv * Bt) * Ainv * G;
    Hf = Ct * Ainv * B * Sinv * F;
    H = Hg + Hf;

    // Build our libMesh data structures from the Eigen ones
    E_libmesh.resize(E.rows(), E.cols());
    H_libmesh.resize(H.size());
    libmesh_assert((E.rows() == E.cols()) && (E.rows() == H.size()));
    for (const auto i : make_range(E.rows()))
    {
      H_libmesh(i) = H(i);
      for (const auto j : make_range(E.cols()))
        E_libmesh(i, j) = E(i, j);
    }

    if (global_solve)
    {
      // We were performing our finite element assembly for the implicit solve step of our
      // example. Add our local element vectors/matrices into the global system
      dof_map.constrain_element_matrix_and_vector(E_libmesh, H_libmesh, lambda_dof_indices);
      matrix.add_matrix(E_libmesh, lambda_dof_indices);
      lambda_system.rhs->add_vector(H_libmesh, lambda_dof_indices);
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
      const auto scalar_soln = Sinv * Bt * Ainv * G - Sinv * F - Sinv * Bt * Ainv * C * Lambda;
      const auto vector_soln = Ainv * (G - B * scalar_soln - C * Lambda);
      for (const auto i : make_range(vector_n_dofs))
        system.solution->set(vector_dof_indices[i], vector_soln(i));
      for (const auto i : make_range(scalar_n_dofs))
        system.solution->set(scalar_dof_indices[i], scalar_soln(i));

#if !defined(LIBMESH_HAVE_PETSC) || !defined(LIBMESH_USE_REAL_NUMBERS)
      if (!simplicial)
        // We don't support SVD solves in this configuration
        continue;
#endif

      //
      // Now solve for the enriched scalar solution using our Lagrange multiplier solution and, for
      // non-simplexes, the lower-order scalar solution. Note that the Lagrange multiplier
      // represents the trace of p so it is a logical choice to leverage in this postprocessing
      // stage!
      //

      dof_map.dof_indices(elem, enriched_scalar_dof_indices, system.variable_number("p_enriched"));
      const auto enriched_scalar_n_dofs = enriched_scalar_dof_indices.size();
      const auto m = simplicial ? lambda_n_dofs : lambda_n_dofs + 1;
      const auto n = enriched_scalar_n_dofs;

      K_enriched_scalar.resize(m, n);
      F_enriched_scalar.resize(m);
      U_enriched_scalar.resize(n);

      // L2 projection of the enriched scalar into the LM space
      for (const auto side : elem->side_index_range())
      {
        vector_fe_face->reinit(elem, side); // for JxW_face and qface_point
        enriched_scalar_fe_face->reinit(elem, side);
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
              scalar_value = DivGradExactSolution().scalar(xf, yf);
            else if (dim == 3)
              scalar_value = DivGradExactSolution().scalar(xf, yf, zf);

            for (const auto i : make_range(lambda_n_dofs))
            {
              F_enriched_scalar(i) += JxW_face[qp] * lambda_phi_face[i][qp] * scalar_value;
              for (const auto j : make_range(enriched_scalar_n_dofs))
                K_enriched_scalar(i, j) +=
                    JxW_face[qp] * lambda_phi_face[i][qp] * enriched_scalar_phi_face[j][qp];
            }
          }
        else
        {
          // compute local face lambda solution
          lambda_qps.resize(qface.n_points());
          for (auto & lambda_qp : lambda_qps)
            lambda_qp = 0;
          for (const auto qp : make_range(qface.n_points()))
            for (const auto i : index_range(lambda_solution_std_vec))
              lambda_qps[qp] += lambda_solution_std_vec[i] * lambda_phi_face[i][qp];

          for (const auto qp : make_range(qface.n_points()))
            for (const auto i : make_range(lambda_n_dofs))
            {
              F_enriched_scalar(i) += JxW_face[qp] * lambda_phi_face[i][qp] * lambda_qps[qp];
              for (const auto j : make_range(enriched_scalar_n_dofs))
                K_enriched_scalar(i, j) +=
                    JxW_face[qp] * lambda_phi_face[i][qp] * enriched_scalar_phi_face[j][qp];
            }
        }
      }

      if (simplicial)
        K_enriched_scalar.lu_solve(F_enriched_scalar, U_enriched_scalar);
      else
      {
        //
        // For tensor product elements, the system of equations coming from the L2 projection into
        // the Lagrange multiplier space is singular. To make the system nonsingular, we add an L2
        // projection into the low-order scalar solution space.
        //

        enriched_scalar_fe->reinit(elem);
        libmesh_assert(scalar_n_dofs == 1);
        libmesh_assert(scalar_soln.size() == 1);
        libmesh_assert(scalar_phi.size() == 1);

        scalar_qps.resize(qrule.n_points());
        for (auto & scalar_qp : scalar_qps)
          scalar_qp = 0;
        for (const auto qp : make_range(qrule.n_points()))
          scalar_qps[qp] += scalar_soln(0) * scalar_phi[0][qp];
        for (const auto qp : make_range(qrule.n_points()))
        {
          F_enriched_scalar(n) += JxW[qp] * scalar_phi[0][qp] * scalar_qps[qp];
          for (const auto j : make_range(n))
            K_enriched_scalar(n, j) += JxW[qp] * scalar_phi[0][qp] * enriched_scalar_phi[j][qp];
        }
        // We have more rows than columns so an LU solve isn't an option
        K_enriched_scalar.svd_solve(F_enriched_scalar, U_enriched_scalar);
      }

      // Our solution for the local enriched scalar dofs is complete. Insert into the global vector
      system.solution->insert(U_enriched_scalar, enriched_scalar_dof_indices);
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
assemble_divgrad(EquationSystems & es, const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "Lambda");
  fe_assembly(es, /*global_solve=*/true);
}

#else

int
main()
{
  return 0;
}

#endif
