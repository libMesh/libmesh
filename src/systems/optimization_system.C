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



// C++ includes

// Local includes
#include "libmesh/equation_systems.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/optimization_solver.h"
#include "libmesh/optimization_system.h"

namespace libMesh
{

// ------------------------------------------------------------
// OptimizationSystem implementation
OptimizationSystem::OptimizationSystem (EquationSystems & es,
                                        const std::string & name_in,
                                        const unsigned int number_in) :

  Parent(es, name_in, number_in),
  optimization_solver(OptimizationSolver<Number>::build(*this)),
  C_eq(NumericVector<Number>::build(this->comm())),
  C_eq_jac(SparseMatrix<Number>::build(this->comm())),
  C_ineq(NumericVector<Number>::build(this->comm())),
  C_ineq_jac(SparseMatrix<Number>::build(this->comm())),
  lambda_eq(NumericVector<Number>::build(this->comm())),
  lambda_ineq(NumericVector<Number>::build(this->comm()))
{
}



OptimizationSystem::~OptimizationSystem ()
{
  // Clear data
  this->clear();
}



void OptimizationSystem::clear ()
{
  // clear the optimization solver
  optimization_solver->clear();

  // clear the parent data
  Parent::clear();
}


void OptimizationSystem::init_data ()
{
  this->add_vector("lower_bounds");
  this->add_vector("upper_bounds");

  Parent::init_data();

  optimization_solver->clear();
}


void OptimizationSystem::reinit ()
{
  optimization_solver->clear();

  Parent::reinit();
}


void OptimizationSystem::
initialize_equality_constraints_storage(const std::vector<std::set<numeric_index_type>> & constraint_jac_sparsity)
{
  unsigned int n_eq_constraints = constraint_jac_sparsity.size();

  // Assign rows to each processor as evenly as possible
  unsigned int n_procs = comm().size();
  unsigned int n_local_rows = n_eq_constraints / n_procs;
  if (comm().rank() < (n_eq_constraints % n_procs))
    n_local_rows++;

  C_eq->init(n_eq_constraints, n_local_rows, false, PARALLEL);
  lambda_eq->init(n_eq_constraints, n_local_rows, false, PARALLEL);

  // Get the maximum number of non-zeros per row
  unsigned int max_nnz = 0;
  for (unsigned int i=0; i<n_eq_constraints; i++)
    {
      unsigned int nnz = constraint_jac_sparsity[i].size();
      if (nnz > max_nnz)
        max_nnz = nnz;
    }

  C_eq_jac->init(n_eq_constraints,
                 get_dof_map().n_dofs(),
                 n_local_rows,
                 get_dof_map().n_local_dofs(),
                 max_nnz,
                 max_nnz);

  eq_constraint_jac_sparsity = constraint_jac_sparsity;
}


void OptimizationSystem::
initialize_inequality_constraints_storage(const std::vector<std::set<numeric_index_type>> & constraint_jac_sparsity)
{
  unsigned int n_ineq_constraints = constraint_jac_sparsity.size();

  // Assign rows to each processor as evenly as possible
  unsigned int n_procs = comm().size();
  unsigned int n_local_rows = n_ineq_constraints / n_procs;
  if (comm().rank() < (n_ineq_constraints % n_procs))
    n_local_rows++;

  C_ineq->init(n_ineq_constraints, n_local_rows, false, PARALLEL);
  lambda_ineq->init(n_ineq_constraints, n_local_rows, false, PARALLEL);

  // Get the maximum number of non-zeros per row
  unsigned int max_nnz = 0;
  for (unsigned int i=0; i<n_ineq_constraints; i++)
    {
      unsigned int nnz = constraint_jac_sparsity[i].size();
      if (nnz > max_nnz)
        max_nnz = nnz;
    }

  C_ineq_jac->init(n_ineq_constraints,
                   get_dof_map().n_dofs(),
                   n_local_rows,
                   get_dof_map().n_local_dofs(),
                   max_nnz,
                   max_nnz);

  ineq_constraint_jac_sparsity = constraint_jac_sparsity;
}


void OptimizationSystem::solve ()
{
  START_LOG("solve()", "OptimizationSystem");

  optimization_solver->init();
  optimization_solver->solve ();

  STOP_LOG("solve()", "OptimizationSystem");

  this->update();
}


} // namespace libMesh
