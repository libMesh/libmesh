// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// <h1>Optimization Example 2 - Optimization with constraints</h1>
// \author David Knezevic
// \date 2015
//
// In this example we extend example 1 to demonstrate how to use
// OptimizationSystem's interface for imposing equality and inequality
// constraints.

// C++ include files that we need
#include <iostream>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/optimization_system.h"
#include "libmesh/optimization_solver.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

/**
 * This class encapsulate all functionality required for assembling
 * the objective function, gradient, hessian, and constraints.
 */
class AssembleOptimization :
  public OptimizationSystem::ComputeObjective,
  public OptimizationSystem::ComputeGradient,
  public OptimizationSystem::ComputeHessian,
  public OptimizationSystem::ComputeEqualityConstraints,
  public OptimizationSystem::ComputeEqualityConstraintsJacobian,
  public OptimizationSystem::ComputeInequalityConstraints,
  public OptimizationSystem::ComputeInequalityConstraintsJacobian,
  public OptimizationSystem::ComputeLowerAndUpperBounds
{
private:

  /**
   * Keep a reference to the OptimizationSystem.
   */
  OptimizationSystem & _sys;

public:

  /**
   * Constructor.
   */
  AssembleOptimization(OptimizationSystem & sys_in);

  /**
   * The optimization problem we consider here is:
   *   min_U 0.5*U^T A U - U^T F.
   * In this method, we assemble A and F.
   */
  void assemble_A_and_F();

  /**
   * Evaluate the objective function.
   */
  virtual Number objective (const NumericVector<Number> & soln,
                            OptimizationSystem & /*sys*/);

  /**
   * Evaluate the gradient.
   */
  virtual void gradient (const NumericVector<Number> & soln,
                         NumericVector<Number> & grad_f,
                         OptimizationSystem & /*sys*/);

  /**
   * Evaluate the Hessian.
   */
  virtual void hessian (const NumericVector<Number> & soln,
                        SparseMatrix<Number> & H_f,
                        OptimizationSystem & /*sys*/);

  /**
   * Evaluate the equality constraints.
   */
  virtual void equality_constraints (const NumericVector<Number> & X,
                                     NumericVector<Number> & C_eq,
                                     OptimizationSystem & /*sys*/);

  /**
   * Evaluate the equality constraints Jacobian.
   */
  virtual void equality_constraints_jacobian (const NumericVector<Number> & X,
                                              SparseMatrix<Number> & C_eq_jac,
                                              OptimizationSystem & /*sys*/);

  /**
   * Evaluate the inequality constraints.
   */
  virtual void inequality_constraints (const NumericVector<Number> & X,
                                       NumericVector<Number> & C_ineq,
                                       OptimizationSystem & /*sys*/);

  /**
   * Evaluate the inequality constraints Jacobian.
   */
  virtual void inequality_constraints_jacobian (const NumericVector<Number> & X,
                                                SparseMatrix<Number> & C_ineq_jac,
                                                OptimizationSystem & /*sys*/);

  /**
   * Evaluate the lower and upper bounds vectors.
   */
  virtual void lower_and_upper_bounds (OptimizationSystem & sys);

  /**
   * Sparse matrix for storing the matrix A. We use
   * this to facilitate computation of objective, gradient
   * and hessian.
   */
  SparseMatrix<Number> * A_matrix;

  /**
   * Vector for storing F. We use this to facilitate
   * computation of objective, gradient and hessian.
   */
  NumericVector<Number> * F_vector;
};

AssembleOptimization::AssembleOptimization(OptimizationSystem & sys_in) :
  _sys(sys_in)
{}

void AssembleOptimization::assemble_A_and_F()
{
  A_matrix->zero();
  F_vector->zero();

  const MeshBase & mesh = _sys.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int u_var = _sys.variable_number ("u");

  const DofMap & dof_map = _sys.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  std::vector<dof_id_type> dof_indices;

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      dof_map.dof_indices (elem, dof_indices);

      const unsigned int n_dofs = dof_indices.size();

      fe->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          for (unsigned int dof_i=0; dof_i<n_dofs; dof_i++)
            {
              for (unsigned int dof_j=0; dof_j<n_dofs; dof_j++)
                {
                  Ke(dof_i, dof_j) += JxW[qp] * (dphi[dof_j][qp]* dphi[dof_i][qp]);
                }
              Fe(dof_i) += JxW[qp] * phi[dof_i][qp];
            }
        }

      A_matrix->add_matrix (Ke, dof_indices);
      F_vector->add_vector (Fe, dof_indices);
    }

  A_matrix->close();
  F_vector->close();
}

Number AssembleOptimization::objective (const NumericVector<Number> & soln,
                                        OptimizationSystem & /*sys*/)
{
  UniquePtr<NumericVector<Number> > AxU = soln.zero_clone();

  A_matrix->vector_mult(*AxU, soln);
  Number UTxAxU = AxU->dot(soln);

  Number UTxF = F_vector->dot(soln);

  return 0.5 * UTxAxU - UTxF;
}

void AssembleOptimization::gradient (const NumericVector<Number> & soln,
                                     NumericVector<Number> & grad_f,
                                     OptimizationSystem & /*sys*/)
{
  grad_f.zero();

  A_matrix->vector_mult(grad_f, soln);
  grad_f.add(-1, *F_vector);
}


void AssembleOptimization::hessian (const NumericVector<Number> & /*soln*/,
                                    SparseMatrix<Number> & H_f,
                                    OptimizationSystem & sys)
{
  H_f.zero();
  H_f.add(1., *A_matrix);

  // We also need to add the Hessian of the inequality and equality constraints,
  //  \sum_i^n_eq lambda_eq_i H_eq_i + \sum_i^n_ineq lambda_ineq_i H_ineq_i
  // where lambda_eq and lambda_ineq denote Lagrange multipliers associated
  // with the equality and inequality constraints, respectively.
  // In this example, our equality constraints are linear, hence H_eq_i = 0.
  // However, our inequality constraint is nonlinear, so it will contribute
  // to the Hessian matrix.
  sys.optimization_solver->get_dual_variables();
  dof_id_type ineq_index = 0;
  Number lambda_ineq_0 = 0.;
  unsigned int lambda_rank = 0;
  if ((sys.lambda_ineq->first_local_index() <= ineq_index) &&
      (ineq_index < sys.lambda_ineq->last_local_index()))
    {
      lambda_ineq_0 = (*sys.lambda_ineq)(0);
      lambda_rank = sys.comm().rank();
    }

  // Sync lambda_rank across all processors.
  sys.comm().sum(lambda_rank);
  sys.comm().broadcast(lambda_rank, lambda_rank);

  if ((sys.get_dof_map().first_dof() <= 200) && (200 < sys.get_dof_map().end_dof()))
    H_f.add(200, 200, 2. * lambda_ineq_0);
}

void AssembleOptimization::equality_constraints (const NumericVector<Number> & X,
                                                 NumericVector<Number> & C_eq,
                                                 OptimizationSystem & /*sys*/)
{
  C_eq.zero();

  UniquePtr<NumericVector<Number> > X_localized =
    NumericVector<Number>::build(X.comm());
  X_localized->init(X.size(), false, SERIAL);
  X.localize(*X_localized);

  std::vector<Number> constraint_values(3);
  constraint_values[0] = (*X_localized)(17);
  constraint_values[1] = (*X_localized)(23);
  constraint_values[2] = (*X_localized)(98) + (*X_localized)(185);

  for (std::size_t i=0; i<constraint_values.size(); i++)
    if ((C_eq.first_local_index() <= i) &&
        (i < C_eq.last_local_index()))
      C_eq.set(i, constraint_values[i]);
}

void AssembleOptimization::equality_constraints_jacobian (const NumericVector<Number> & /*X*/,
                                                          SparseMatrix<Number> & C_eq_jac,
                                                          OptimizationSystem & sys)
{
  C_eq_jac.zero();

  std::vector< std::vector<Number> > constraint_jac_values(3);
  std::vector< std::vector<dof_id_type> > constraint_jac_indices(3);

  constraint_jac_values[0].resize(1);
  constraint_jac_indices[0].resize(1);
  constraint_jac_values[0][0] = 1.;
  constraint_jac_indices[0][0] = 17;

  constraint_jac_values[1].resize(1);
  constraint_jac_indices[1].resize(1);
  constraint_jac_values[1][0] = 1.;
  constraint_jac_indices[1][0] = 23;

  constraint_jac_values[2].resize(2);
  constraint_jac_indices[2].resize(2);
  constraint_jac_values[2][0] = 1.;
  constraint_jac_values[2][1] = 1.;
  constraint_jac_indices[2][0] = 98;
  constraint_jac_indices[2][1] = 185;

  for (std::size_t i=0; i<constraint_jac_values.size(); i++)
    for (std::size_t j=0; j<constraint_jac_values[i].size(); j++)
      if ((sys.C_eq->first_local_index() <= i) &&
          (i < sys.C_eq->last_local_index()))
        {
          dof_id_type col_index = constraint_jac_indices[i][j];
          Number value = constraint_jac_values[i][j];
          C_eq_jac.set(i, col_index, value);
        }
}

void AssembleOptimization::inequality_constraints (const NumericVector<Number> & X,
                                                   NumericVector<Number> & C_ineq,
                                                   OptimizationSystem & /*sys*/)
{
  C_ineq.zero();

  UniquePtr<NumericVector<Number> > X_localized =
    NumericVector<Number>::build(X.comm());
  X_localized->init(X.size(), false, SERIAL);
  X.localize(*X_localized);

  std::vector<Number> constraint_values(1);
  constraint_values[0] = (*X_localized)(200)*(*X_localized)(200) + (*X_localized)(201) - 5.;

  for (std::size_t i=0; i<constraint_values.size(); i++)
    if ((C_ineq.first_local_index() <= i) && (i < C_ineq.last_local_index()))
      C_ineq.set(i, constraint_values[i]);
}

void AssembleOptimization::inequality_constraints_jacobian (const NumericVector<Number> & X,
                                                            SparseMatrix<Number> & C_ineq_jac,
                                                            OptimizationSystem & sys)
{
  C_ineq_jac.zero();

  UniquePtr<NumericVector<Number> > X_localized =
    NumericVector<Number>::build(X.comm());
  X_localized->init(X.size(), false, SERIAL);
  X.localize(*X_localized);

  std::vector< std::vector<Number> > constraint_jac_values(1);
  std::vector< std::vector<dof_id_type> > constraint_jac_indices(1);

  constraint_jac_values[0].resize(2);
  constraint_jac_indices[0].resize(2);
  constraint_jac_values[0][0] = 2.* (*X_localized)(200);
  constraint_jac_values[0][1] = 1.;
  constraint_jac_indices[0][0] = 200;
  constraint_jac_indices[0][1] = 201;

  for (std::size_t i=0; i<constraint_jac_values.size(); i++)
    for (std::size_t j=0; j<constraint_jac_values[i].size(); j++)
      if ((sys.C_ineq->first_local_index() <= i) &&
          (i < sys.C_ineq->last_local_index()))
        {
          dof_id_type col_index = constraint_jac_indices[i][j];
          Number value = constraint_jac_values[i][j];
          C_ineq_jac.set(i, col_index, value);
        }

}

void AssembleOptimization::lower_and_upper_bounds (OptimizationSystem & sys)
{
  for (unsigned int i=sys.get_dof_map().first_dof(); i<sys.get_dof_map().end_dof(); i++)
    {
      sys.get_vector("lower_bounds").set(i, -2.);
      sys.get_vector("upper_bounds").set(i, 2.);
    }
}

int main (int argc, char ** argv)
{
  LibMeshInit init (argc, argv);

#ifndef LIBMESH_HAVE_PETSC_TAO

  libmesh_example_requires(false, "PETSc >= 3.5.0 with built-in TAO support");

#elif LIBMESH_USE_COMPLEX_NUMBERS

  // According to
  // http://www.mcs.anl.gov/research/projects/tao/documentation/installation.html
  // TAO & PETSc-complex are currently mutually exclusive
  libmesh_example_requires(false, "PETSc >= 3.5.0 with built-in TAO support & real-numbers only");

#endif

  // We use a 2D domain.
  libmesh_example_requires(LIBMESH_DIM > 1, "--disable-1D-only");

  // TAO is giving us problems in parallel?
  if (init.comm().size() != 1)
    {
      libMesh::out << "This example can currently only be run in serial." << std::endl;
      return 77;
    }

  GetPot infile("optimization_ex2.in");
  const std::string approx_order = infile("approx_order", "FIRST");
  const std::string fe_family = infile("fe_family", "LAGRANGE");
  const unsigned int n_elem = infile("n_elem", 10);

  Mesh mesh(init.comm());
  MeshTools::Generation::build_square (mesh,
                                       n_elem,
                                       n_elem,
                                       -1., 1.,
                                       -1., 1.,
                                       QUAD9);

  mesh.print_info();

  EquationSystems equation_systems (mesh);

  OptimizationSystem & system =
    equation_systems.add_system<OptimizationSystem> ("Optimization");

  // The default is to use PETSc/Tao solvers, but let the user change
  // the optimization solver package on the fly.
  {
    const std::string optimization_solver_type = infile("optimization_solver_type",
                                                        "PETSC_SOLVERS");
    SolverPackage sp = Utility::string_to_enum<SolverPackage>(optimization_solver_type);
    UniquePtr<OptimizationSolver<Number> > new_solver =
      OptimizationSolver<Number>::build(system, sp);
    system.optimization_solver.reset(new_solver.release());
  }

  // Set tolerances and maximum iteration counts directly on the optimization solver.
  system.optimization_solver->max_objective_function_evaluations = 128;
  system.optimization_solver->objective_function_relative_tolerance = 1.e-4;
  system.optimization_solver->verbose = true;

  AssembleOptimization assemble_opt(system);

  system.add_variable("u",
                      Utility::string_to_enum<Order>   (approx_order),
                      Utility::string_to_enum<FEFamily>(fe_family));

  system.optimization_solver->objective_object = &assemble_opt;
  system.optimization_solver->gradient_object = &assemble_opt;
  system.optimization_solver->hessian_object = &assemble_opt;
  system.optimization_solver->equality_constraints_object = &assemble_opt;
  system.optimization_solver->equality_constraints_jacobian_object = &assemble_opt;
  system.optimization_solver->inequality_constraints_object = &assemble_opt;
  system.optimization_solver->inequality_constraints_jacobian_object = &assemble_opt;
  system.optimization_solver->lower_and_upper_bounds_object = &assemble_opt;

  // system.matrix and system.rhs are used for the gradient and Hessian,
  // so in this case we add an extra matrix and vector to store A and F.
  // This makes it easy to write the code for evaluating the objective,
  // gradient, and hessian.
  system.add_matrix("A_matrix");
  system.add_vector("F_vector");
  assemble_opt.A_matrix = &system.get_matrix("A_matrix");
  assemble_opt.F_vector = &system.get_vector("F_vector");


  equation_systems.init();
  equation_systems.print_info();

  assemble_opt.assemble_A_and_F();

  {
    std::vector< std::set<numeric_index_type> > constraint_jac_sparsity;
    std::set<numeric_index_type> sparsity_row;
    sparsity_row.insert(17);
    constraint_jac_sparsity.push_back(sparsity_row);
    sparsity_row.clear();

    sparsity_row.insert(23);
    constraint_jac_sparsity.push_back(sparsity_row);
    sparsity_row.clear();

    sparsity_row.insert(98);
    sparsity_row.insert(185);
    constraint_jac_sparsity.push_back(sparsity_row);

    system.initialize_equality_constraints_storage(constraint_jac_sparsity);
  }

  {
    std::vector< std::set<numeric_index_type> > constraint_jac_sparsity;
    std::set<numeric_index_type> sparsity_row;
    sparsity_row.insert(200);
    sparsity_row.insert(201);
    constraint_jac_sparsity.push_back(sparsity_row);

    system.initialize_inequality_constraints_storage(constraint_jac_sparsity);
  }

  // We need to close the matrix so that we can use it to store the
  // Hessian during the solve.
  system.matrix->close();
  system.solve();

  // Print convergence information
  system.optimization_solver->print_converged_reason();

  // Write the solution to file if the optimization solver converged
  if (system.optimization_solver->get_converged_reason() > 0)
    {
      std::stringstream filename;
      ExodusII_IO (mesh).write_equation_systems("optimization_soln.exo",
                                                equation_systems);
    }

  return 0;
}
