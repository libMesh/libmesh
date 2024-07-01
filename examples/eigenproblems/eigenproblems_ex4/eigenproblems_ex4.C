// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// <h1>Eigenproblems Example 4 - Constrained Eigen Problem Solved with Newton Method</h1>
// \author Alexander Lindsay
// \date 2024
//
// This example illustrates how to solve an eigen problem with constraints from hanging
// nodes using Newton's method, which is implemented in SLEPc as a variant of the power method.
// Unlike for the standard power method, the Bx function is computed using the current solution.
// This example will only work if the library is compiled with SLEPc support enabled.
//
// The Ax function is composed of diffusion with vacuum boundary conditions (-n * grad u
// = u). The Bx function corresponds to a mass term

// libMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/enum_eigen_solver_type.h"
#include "libmesh/petsc_shell_matrix.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/slepc_macro.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifdef LIBMESH_HAVE_SLEPC
#include <petscsnes.h>

PetscErrorCode form_functionA(SNES snes, Vec x, Vec Ax, void * ctx);
PetscErrorCode form_matrixA(SNES snes, Vec x, Mat jac, Mat pc, void * ctx);
PetscErrorCode form_functionB(SNES snes, Vec x, Vec Bx, void * ctx);
#endif

int
main(int argc, char ** argv)
{
  // Initialize libMesh and the dependent libraries.
  LibMeshInit init(argc, argv);

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc currently gives us a nasty crash with Real==float
  libmesh_example_requires(false, "--disable-singleprecision");
#endif

#ifndef LIBMESH_HAVE_SLEPC
  libmesh_example_requires(false, "--enable-slepc with a slepc version >=3.13");
#else
#if SLEPC_VERSION_LESS_THAN(3, 13, 0)
  libmesh_example_requires(false, "--enable-slepc with a slepc version >=3.13");
#else
  // Tell the user what we are doing.
  libMesh::out << "Running " << argv[0];

  for (int i = 1; i < argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

  // We only solve for a single eigenvalue with the nonlinear power method
  constexpr int nev = 1;

  // Possibly get the mesh size from -nx and -ny
  const int nx = libMesh::command_line_next("-nx", 20), ny = libMesh::command_line_next("-ny", 20);

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Create a mesh, with dimension to be overridden later, on the
  // default MPI communicator.
  Mesh mesh(init.comm());

  // Use the internal mesh generator to create a uniform
  // 2D grid on a square.
  MeshTools::Generation::build_square(mesh, nx, ny, -1., 1., -1., 1., QUAD4);

#ifdef LIBMESH_ENABLE_AMR
  for (auto * elem : mesh.active_element_ptr_range())
    if (elem->vertex_average()(0) < 0)
      elem->set_refinement_flag(Elem::REFINE);

  MeshRefinement refinedmesh(mesh);
  refinedmesh.refine_elements();
#endif

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems(mesh);

  // Create a CondensedEigenSystem named "Eigensystem" and (for convenience)
  // use a reference to the system we create.
  auto & eigen_system = equation_systems.add_system<CondensedEigenSystem>("Condensed Eigensystem");

  // Declare the system variables.
  // Adds the variable "p" to "Eigensystem".   "p"
  // will be approximated using second-order approximation.
  eigen_system.add_variable("p", FIRST);

  // Set necessary parameters used in EigenSystem::solve(),
  // i.e. the number of requested eigenpairs nev and the number
  // of basis vectors ncv used in the solution algorithm. Note that
  // ncv >= nev must hold and ncv >= 2*nev is recommended.
  equation_systems.parameters.set<unsigned int>("eigenpairs") = nev;
  equation_systems.parameters.set<unsigned int>("basis vectors") = nev * 3;

  eigen_system.set_eigenproblem_type(GNHEP);
  eigen_system.use_shell_matrices(true);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  auto ierr = PetscOptionsSetValue(LIBMESH_PETSC_NULLPTR, "-eps_type", "power");
  CHKERRQ(ierr);
  ierr = PetscOptionsSetValue(LIBMESH_PETSC_NULLPTR, "-eps_power_update", "1");
  CHKERRQ(ierr);
  ierr = PetscOptionsSetValue(LIBMESH_PETSC_NULLPTR, "-eps_power_nonlinear", "1");
  CHKERRQ(ierr);
  ierr = PetscOptionsSetValue(LIBMESH_PETSC_NULLPTR, "-eps_max_it", "1");
  CHKERRQ(ierr);
  ierr = PetscOptionsSetValue(LIBMESH_PETSC_NULLPTR, "-eps_power_snes_mf_operator", "1");
  CHKERRQ(ierr);
  ierr = PetscOptionsSetValue(LIBMESH_PETSC_NULLPTR, "-eps_power_pc_type", "lu");
  CHKERRQ(ierr);

  //
  // Set function/operator callback functions
  //

  auto & A = static_cast<PetscShellMatrix<Number> &>(eigen_system.get_shell_matrix_A());
  auto Amat = A.mat();
  ierr = PetscObjectComposeFunction((PetscObject)Amat, "formFunction", form_functionA);
  CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)Amat, "formJacobian", form_matrixA);
  CHKERRQ(ierr);

  auto & B = static_cast<PetscShellMatrix<Number> &>(eigen_system.get_shell_matrix_B());
  auto Bmat = B.mat();
  ierr = PetscObjectComposeFunction((PetscObject)Bmat, "formFunction", form_functionB);
  CHKERRQ(ierr);

  //
  // Set function/operator callback function contexts
  //

  PetscContainer container;
  ierr = PetscContainerCreate(equation_systems.comm().get(), &container);
  CHKERRQ(ierr);
  ierr = PetscContainerSetPointer(container, &equation_systems);
  CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)Amat, "formFunctionCtx", (PetscObject)container);
  CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)Amat, "formJacobianCtx", (PetscObject)container);
  CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)Bmat, "formFunctionCtx", (PetscObject)container);
  CHKERRQ(ierr);
  ierr = PetscContainerDestroy(&container);
  CHKERRQ(ierr);

  // Set the initial space
  Vec initial_space;
  ierr = MatCreateVecs(Amat, &initial_space, nullptr);
  CHKERRQ(ierr);
  PetscVector<Number> wrapped_initial_space(initial_space, equation_systems.comm());
  wrapped_initial_space.add(1);
  wrapped_initial_space.close();
  eigen_system.get_eigen_solver().set_initial_space(wrapped_initial_space);

  // Solve the system "Condensed Eigensystem"
  eigen_system.initialize_condensed_dofs();
  eigen_system.dont_create_submatrices_in_solve();
  eigen_system.assemble_before_solve = false;
  eigen_system.get_eigen_solver().set_close_matrix_before_solve(false);
  eigen_system.solve();

  // Get the number of converged eigen pairs.
  unsigned int nconv = eigen_system.get_n_converged();

  // Get the last converged eigenpair
  if (nconv != 0)
  {
    auto [re, im] = eigen_system.get_eigenpair(0);
    libMesh::out << "The converged eigenvalue is " << re << " + " << im << "i" << std::endl;

#ifdef LIBMESH_HAVE_EXODUS_API
    // Write the eigen vector to file.
    ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  }
  else
    libMesh::out << "WARNING: Solver did not converge!\n" << nconv << std::endl;

  ierr = VecDestroy(&initial_space);
  CHKERRQ(ierr);

  // All done.
  return 0;
#endif // SLEPC version >= 3.13.0
#endif // LIBMESH_HAVE_SLEPC
}

#ifdef LIBMESH_HAVE_SLEPC
void
update_current_local_solution(CondensedEigenSystem & sys, Vec x)
{
  auto & dof_map = sys.get_dof_map();

  PetscVector<Number> X_global(x, sys.comm());

  if (dof_map.n_constrained_dofs())
  {
    sys.copy_sub_to_super(X_global, *sys.solution);
    // Set the constrained dof values
    dof_map.enforce_constraints_exactly(sys);
    sys.update();
  }
  else
  {
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());

    // Use the system's update() to get a good local version of the
    // parallel solution.  This operation does not modify the incoming
    // "x" vector, it only localizes information from "x" into
    // sys.current_local_solution.
    X_global.swap(X_sys);
    sys.update();
    X_global.swap(X_sys);
  }
}

std::unique_ptr<NumericVector<Number>>
create_wrapped_function(CondensedEigenSystem & sys, Vec f)
{
  auto & dof_map = sys.get_dof_map();

  if (dof_map.n_constrained_dofs())
    return sys.solution->zero_clone();
  else
  {
    auto F = std::make_unique<PetscVector<Number>>(f, sys.comm());
    F->zero();
    return F;
  }
}

PetscErrorCode
form_functionA(SNES /*snes*/, Vec x, Vec Ax, void * ctx)
{
  PetscFunctionBegin;

  auto & es = *static_cast<EquationSystems *>(ctx);

  // Get a reference to our system.
  auto & eigen_system = es.get_system<CondensedEigenSystem>("Condensed Eigensystem");

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  auto & dof_map = eigen_system.get_dof_map();

  // Create our data structures correponding to both constrained and unconstrained degrees of
  // freedom
  update_current_local_solution(eigen_system, x);
  auto AX = create_wrapped_function(eigen_system, Ax);

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  QGauss qrule(dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule(&qrule);

  // Do the same for faces
  std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));
  QGauss qrule_face(dim - 1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule(&qrule_face);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const auto & dphi = fe->get_dphi();

  // The element degree of freedom values
  std::vector<Number> Ue;

  // The diffusion residual
  DenseVector<Number> Ke;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users
  // later modify this program to include refinement, we will
  // be safe and will only consider the active elements;
  // hence we use a variant of the active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit(elem);

    // Get the element degree of freedom values
    eigen_system.current_local_solution->get(dof_indices, Ue);

    // Zero the element matrices and rhs before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).
    const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());
    Ke.resize(n_dofs);

    // Now loop over the quadrature points.  This handles
    // the numeric integration.
    //
    // We will build the element vector
    for (const auto qp : make_range(qrule.n_points()))
    {
      // Build the solution at the quadrature point
      Gradient grad_u_qp = 0;
      for (const auto i : make_range(n_dofs))
        grad_u_qp += dphi[i][qp] * Ue[i];
      // Now compute the action of the mass operator on the solution
      for (const auto i : make_range(n_dofs))
        Ke(i) += JxW[qp] * dphi[i][qp] * grad_u_qp;
    }

    for (const auto s : elem->side_index_range())
      if (!elem->neighbor_ptr(s))
      {
        const auto & phi_face = fe_face->get_phi();
        const auto & JxW_face = fe_face->get_JxW();
        // vacuum boundary condition
        fe_face->reinit(elem, s);
        for (const auto qp : make_range(qrule_face.n_points()))
        {
          Number u_qp = 0;
          for (const auto i : make_range(n_dofs))
            u_qp += phi_face[i][qp] * Ue[i];
          for (const auto i : make_range(n_dofs))
            Ke(i) += JxW_face[qp] * phi_face[i][qp] * u_qp;
        }
      }

    // On an unrefined mesh, constrain_element_vector does
    // nothing.  If the assembly function is
    // run on a refined mesh, getting the hanging node constraints
    // right will be important
    dof_map.constrain_element_vector(Ke, dof_indices);

    // Finally, simply add the element contribution to the
    // overall matrix.
    AX->add_vector(Ke, dof_indices);
  } // end of element loop

  AX->close();

  if (dof_map.n_constrained_dofs())
  {
    PetscVector<Number> wrapped_Ax(Ax, eigen_system.comm());
    eigen_system.copy_super_to_sub(*AX, wrapped_Ax);
  }

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

PetscErrorCode
form_functionB(SNES /*snes*/, Vec x, Vec Bx, void * ctx)
{
  PetscFunctionBegin;

  auto & es = *static_cast<EquationSystems *>(ctx);

  // Get a reference to our system.
  auto & eigen_system = es.get_system<CondensedEigenSystem>("Condensed Eigensystem");

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  auto & dof_map = eigen_system.get_dof_map();

  // Create our data structures correponding to both constrained and unconstrained degrees of
  // freedom
  update_current_local_solution(eigen_system, x);
  auto BX = create_wrapped_function(eigen_system, Bx);

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  QGauss qrule(dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule(&qrule);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // The element degree of freedom values
  std::vector<Number> Ue;

  // The element mass matrix residual
  DenseVector<Number> Me;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users
  // later modify this program to include refinement, we will
  // be safe and will only consider the active elements;
  // hence we use a variant of the active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit(elem);

    // Get the element degree of freedom values
    eigen_system.current_local_solution->get(dof_indices, Ue);

    // Zero the element matrices and rhs before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).
    const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());
    Me.resize(n_dofs);

    // Now loop over the quadrature points.  This handles
    // the numeric integration.
    //
    // We will build the element vector
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      // Build the solution at the quadrature point
      Number Uqp = 0;
      for (const auto i : make_range(n_dofs))
        Uqp += phi[i][qp] * Ue[i];
      // Now compute the action of the mass operator on the solution
      for (unsigned int i = 0; i != n_dofs; i++)
        Me(i) += JxW[qp] * phi[i][qp] * Uqp;
    }

    // On an unrefined mesh, constrain_element_vector does
    // nothing.  If the assembly function is
    // run on a refined mesh, getting the hanging node constraints
    // right will be important
    dof_map.constrain_element_vector(Me, dof_indices);

    // Finally, simply add the element contribution to the
    // overall matrix.
    BX->add_vector(Me, dof_indices);
  } // end of element loop

  BX->close();

  if (dof_map.n_constrained_dofs())
  {
    PetscVector<Number> wrapped_Bx(Bx, eigen_system.comm());
    eigen_system.copy_super_to_sub(*BX, wrapped_Bx);
  }

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

PetscErrorCode
form_matrixA(SNES /*snes*/, Vec x, Mat jac, Mat pc, void * ctx)
{
  PetscFunctionBegin;

  auto & es = *static_cast<EquationSystems *>(ctx);

  // Get a reference to our system.
  auto & eigen_system = es.get_system<CondensedEigenSystem>("Condensed Eigensystem");

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  auto & dof_map = eigen_system.get_dof_map();

  //
  // Create our data structures correponding to both constrained and unconstrained degrees of
  // freedom
  //

  update_current_local_solution(eigen_system, x);

  PetscBool pisshell, jismffd;

  auto ierr = PetscObjectTypeCompare((PetscObject)pc, MATSHELL, &pisshell);
  CHKERRQ(ierr);
  if (pisshell)
    libmesh_error_msg("Generic preconditioning requires that an explicit matrix representation of "
                      "the preconditioner be formed");
  ierr = PetscObjectTypeCompare((PetscObject)jac, MATMFFD, &jismffd);
  CHKERRQ(ierr);
  if (!jismffd)
    libmesh_error_msg("The operator should be formed matrix free");

  auto & pc_super = eigen_system.get_precond_matrix();

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  QGauss qrule(dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule(&qrule);

  // Do the same for faces
  std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));
  QGauss qrule_face(dim - 1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule(&qrule_face);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape function gradients evaluated at the quadrature points.
  const auto & dphi = fe->get_dphi();

  // The element diffusion matrix
  DenseMatrix<Number> Ke;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users
  // later modify this program to include refinement, we will
  // be safe and will only consider the active elements;
  // hence we use a variant of the active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit(elem);

    // Zero the element matrices and rhs before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).
    const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());
    Ke.resize(n_dofs, n_dofs);

    // Now loop over the quadrature points.  This handles
    // the numeric integration.
    //
    // We will build the element matrix
    for (const auto qp : make_range(qrule.n_points()))
      // Now compute the mass operator
      for (const auto i : make_range(n_dofs))
        for (const auto j : make_range(n_dofs))
          Ke(i, j) += JxW[qp] * dphi[i][qp] * dphi[j][qp];

    for (const auto s : elem->side_index_range())
      if (!elem->neighbor_ptr(s))
      {
        const auto & phi_face = fe_face->get_phi();
        const auto & JxW_face = fe_face->get_JxW();
        // vacuum boundary condition
        fe_face->reinit(elem, s);
        for (const auto qp : make_range(qrule_face.n_points()))
        {
          for (const auto i : make_range(n_dofs))
            for (const auto j : make_range(n_dofs))
              Ke(i, j) += JxW_face[qp] * phi_face[i][qp] * phi_face[j][qp];
        }
      }

    // On an unrefined mesh, constrain_element_vector does
    // nothing.  If the assembly function is
    // run on a refined mesh, getting the hanging node constraints
    // right will be important
    dof_map.constrain_element_matrix(Ke, dof_indices);

    // Finally, simply add the element contribution to the
    // overall matrix.
    pc_super.add_matrix(Ke, dof_indices, dof_indices);
  } // end of element loop

  pc_super.close();

  if (dof_map.n_constrained_dofs())
  {
    PetscMatrix<Number> sub(pc, eigen_system.comm());
    eigen_system.copy_super_to_sub(pc_super, sub);
  }

  // The MFFD Jac still must have assemble called on it
  ierr = MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#endif
