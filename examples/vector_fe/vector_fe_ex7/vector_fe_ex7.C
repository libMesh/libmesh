// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include <Eigen/Dense>

using namespace libMesh;
using namespace Eigen;

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

  // Declare the system  "DivGrad" and its variables.
  auto & system = equation_systems.add_system<System>("DivGrad");

  // Add the LM system
  auto & lm_system = equation_systems.add_system<LinearImplicitSystem>("Lambda");

  // Adds the variable "u" and "p" to "DivGrad". "u" will be our vector field
  // whereas "p" will be the scalar field.
  system.add_variable("u", FIRST, L2_RAVIART_THOMAS);
  system.add_variable("p", CONSTANT, MONOMIAL);

  lm_system.add_variable("lambda", CONSTANT, SIDE_HIERARCHIC);

  // Give the system a pointer to the matrix assembly
  // function. This will be called when needed by the library.
  lm_system.attach_assemble_function(assemble_divgrad);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  lm_system.solve();
  // Now populate the vector and scalar solutions
  fe_assembly(equation_systems, /*global_solve=*/false);

  ExactSolution exact_sol(equation_systems);

  if (dimension == 2)
  {
    SolutionFunction<2> soln_func(system.variable_number("u"));
    SolutionGradient<2> soln_grad(system.variable_number("u"));

    // Build FunctionBase* containers to attach to the ExactSolution object.
    std::vector<FunctionBase<Number> *> sols(1, &soln_func);
    std::vector<FunctionBase<Gradient> *> grads(1, &soln_grad);

    exact_sol.attach_exact_values(sols);
    exact_sol.attach_exact_derivs(grads);
  }
  else if (dimension == 3)
  {
    SolutionFunction<3> soln_func(system.variable_number("u"));
    SolutionGradient<3> soln_grad(system.variable_number("u"));

    // Build FunctionBase* containers to attach to the ExactSolution object.
    std::vector<FunctionBase<Number> *> sols(1, &soln_func);
    std::vector<FunctionBase<Gradient> *> grads(1, &soln_grad);

    exact_sol.attach_exact_values(sols);
    exact_sol.attach_exact_derivs(grads);
  }

  // Use higher quadrature order for more accurate error results.
  int extra_error_quadrature = infile("extra_error_quadrature", 2);
  exact_sol.extra_quadrature_order(extra_error_quadrature);

  // Compute the error.
  exact_sol.compute_error("DivGrad", "u");

  // Print out the error values.
  libMesh::out << "L2 error is: " << exact_sol.l2_error("DivGrad", "u") << std::endl;
  libMesh::out << "HDiv semi-norm error is: " << exact_sol.error_norm("DivGrad", "u", HDIV_SEMINORM)
               << std::endl;
  libMesh::out << "HDiv error is: " << exact_sol.hdiv_error("DivGrad", "u") << std::endl;

#ifdef LIBMESH_HAVE_EXODUS_API

  // We write the file in the ExodusII format.
  ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}

void
fe_assembly(EquationSystems & es, const bool global_solve)
{
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  auto & system = es.get_system<System>("DivGrad");
  auto & lambda_system = es.get_system<LinearImplicitSystem>("Lambda");

  const auto & dof_map = system.get_dof_map();
  const auto & lambda_dof_map = lambda_system.get_dof_map();

  const FEType vector_fe_type = dof_map.variable_type(system.variable_number("u"));
  const FEType scalar_fe_type = dof_map.variable_type(system.variable_number("p"));
  const FEType lambda_fe_type =
      lambda_dof_map.variable_type(lambda_system.variable_number("lambda"));

  // Volumetric FE objects
  std::unique_ptr<FEVectorBase> vector_fe(FEVectorBase::build(dim, vector_fe_type));
  std::unique_ptr<FEBase> scalar_fe(FEBase::build(dim, scalar_fe_type));

  // Volumetric quadrature rule
  QGauss qrule(dim, FIFTH);

  vector_fe->attach_quadrature_rule(&qrule);
  scalar_fe->attach_quadrature_rule(&qrule);

  // Declare finite element objects for boundary integration
  std::unique_ptr<FEVectorBase> vector_fe_face(FEVectorBase::build(dim, vector_fe_type));
  std::unique_ptr<FEBase> lambda_fe_face(FEBase::build(dim, lambda_fe_type));

  // Boundary integration requires one quadrature rule with dimensionality one
  // less than the dimensionality of the element.
  QGauss qface(dim - 1, FIFTH);

  vector_fe_face->attach_quadrature_rule(&qface);
  lambda_fe_face->attach_quadrature_rule(&qface);

  const auto & JxW = vector_fe->get_JxW();
  const auto & q_point = vector_fe->get_xyz();
  const auto & vector_phi = vector_fe->get_phi();
  const auto & scalar_phi = scalar_fe->get_phi();
  const auto & div_vector_phi = vector_fe->get_div_phi();
  const auto & vector_phi_face = vector_fe_face->get_phi();
  const auto & lambda_phi_face = lambda_fe_face->get_phi();
  const auto & JxW_face = vector_fe_face->get_JxW();
  const auto & qface_point = vector_fe_face->get_xyz();
  const auto & normals = vector_fe_face->get_normals();

  // We follow the notation of Cockburn
  // LM matrix and RHS
  DenseMatrix<Number> E_libmesh;
  DenseVector<Number> H_libmesh;
  MatrixXd E;
  MatrixXd H;
  // Auxiliary matrices and RHS
  MatrixXd A, Ainv, B, Bt, C, Ct, Sinv;
  VectorXd G, F, Hg, Hf;
  // element matrix for boundary LM dofs
  MatrixXd L;

  // Lambda eigen vector for constructing vector and scalar solutions
  VectorXd Lambda;

  auto zero_mat = [](auto & mat)
  {
    for (const auto i : make_range(mat.rows()))
      for (const auto j : make_range(mat.cols()))
        mat(i, j) = 0;
  };
  auto zero_vec = [](auto & vec)
  {
    for (const auto i : make_range(vec.size()))
      vec(i) = 0;
  };

  std::vector<dof_id_type> vector_dof_indices;
  std::vector<dof_id_type> scalar_dof_indices;
  std::vector<dof_id_type> lambda_dof_indices;
  std::vector<Number> lambda_solution_std_vec;

  // The global system matrix
  SparseMatrix<Number> & matrix = lambda_system.get_system_matrix();

  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    dof_map.dof_indices(elem, vector_dof_indices, system.variable_number("u"));
    dof_map.dof_indices(elem, scalar_dof_indices, system.variable_number("p"));
    lambda_dof_map.dof_indices(elem, lambda_dof_indices, lambda_system.variable_number("lambda"));

    const auto vector_n_dofs = vector_dof_indices.size();
    const auto scalar_n_dofs = scalar_dof_indices.size();
    const auto lambda_n_dofs = lambda_dof_indices.size();

    vector_fe->reinit(elem);
    scalar_fe->reinit(elem);

    libmesh_assert_equal_to(vector_n_dofs, vector_phi.size());
    libmesh_assert_equal_to(scalar_n_dofs, scalar_phi.size());

    // vector equation
    A.resize(vector_n_dofs, vector_n_dofs);
    B.resize(vector_n_dofs, scalar_n_dofs);
    C.resize(vector_n_dofs, lambda_n_dofs);
    G.resize(vector_n_dofs);
    zero_mat(A);
    zero_mat(B);
    zero_mat(C);
    zero_vec(G);
    // scalar equation
    Bt.resize(scalar_n_dofs, vector_n_dofs);
    F.resize(scalar_n_dofs);
    zero_mat(Bt);
    zero_vec(F);
    // lm equation
    Ct.resize(lambda_n_dofs, vector_n_dofs);
    L.resize(lambda_n_dofs, lambda_n_dofs);
    zero_mat(Ct);
    zero_mat(L);

    for (const auto qp : make_range(qrule.n_points()))
    {
      for (const auto i : make_range(vector_n_dofs))
        for (const auto j : make_range(vector_n_dofs))
          A(i, j) += JxW[qp] * (vector_phi[i][qp] * vector_phi[j][qp]);

      for (const auto i : make_range(vector_n_dofs))
        for (const auto j : make_range(scalar_n_dofs))
          B(i, j) -= JxW[qp] * (div_vector_phi[i][qp] * scalar_phi[j][qp]);

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

        for (const auto i : make_range(scalar_n_dofs))
          F(i) -= JxW[qp] * scalar_phi[i][qp] * f;
      }
    }

    {
      for (auto side : elem->side_index_range())
      {
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
            for (const auto i : make_range(vector_n_dofs))
              for (const auto j : make_range(lambda_n_dofs))
                C(i, j) +=
                    JxW_face[qp] * (vector_phi_face[i][qp] * normals[qp]) * lambda_phi_face[j][qp];

            for (const auto i : make_range(lambda_n_dofs))
              for (const auto j : make_range(vector_n_dofs))
                Ct(i, j) +=
                    JxW_face[qp] * lambda_phi_face[i][qp] * (vector_phi_face[j][qp] * normals[qp]);
          }
      }
    }

    Ainv = A.inverse();
    Sinv = (Bt * Ainv * B).inverse();
    E = Ct * Ainv * (A - B * Sinv * Bt) * Ainv * C;
    E = E + L;
    Hg = Ct * Ainv * (A - B * Sinv * Bt) * Ainv * G;
    Hf = Ct * Ainv * B * Sinv * F;
    H = Hg + Hf;

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
      dof_map.constrain_element_matrix_and_vector(E_libmesh, H_libmesh, lambda_dof_indices);
      matrix.add_matrix(E_libmesh, lambda_dof_indices);
      lambda_system.rhs->add_vector(H_libmesh, lambda_dof_indices);
    }
    else
    {
      Lambda.resize(lambda_n_dofs);
      lambda_system.solution->get(lambda_dof_indices, lambda_solution_std_vec);
      for (const auto i : make_range(lambda_n_dofs))
        Lambda(i) = lambda_solution_std_vec[i];
      const auto scalar_soln = Sinv * Bt * Ainv * G - Sinv * F - Sinv * Bt * Ainv * C * Lambda;
      const auto vector_soln = Ainv * (G - B * scalar_soln - C * Lambda);
      for (const auto i : make_range(vector_n_dofs))
        system.solution->set(vector_dof_indices[i], vector_soln(i));
      for (const auto i : make_range(scalar_n_dofs))
        system.solution->set(scalar_dof_indices[i], scalar_soln(i));
    }
  }

  if (!global_solve)
  {
    system.solution->close();
    // Scatter solution into the current_solution which is used in error computation
    system.update();
  }
}

// We now define the matrix assembly function for the
// div-grad system. We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void
assemble_divgrad(EquationSystems & es, const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "Lambda");
  fe_assembly(es, /*global_solve=*/true);
}
