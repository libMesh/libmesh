// $Id: eigen_system.C 4015 2010-10-01 21:03:53Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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

#include "libmesh_config.h"

// Currently, the EigenSystem should only be available
// if SLEPc support is enabled. 
#if defined(LIBMESH_HAVE_SLEPC)

#include "condensed_eigen_system.h"
#include "libmesh_logging.h"
#include "numeric_vector.h"
#include "equation_systems.h"
#include "getpot.h"
#include "rb_system.h"
#include "parallel.h"

namespace libMesh
{

CondensedEigenSystem::CondensedEigenSystem (EquationSystems& es,
                                            const std::string& name,
                                            const unsigned int number)
  : Parent(es, name, number),
    condensed_matrix_A(SparseMatrix<Number>::build()),
    condensed_matrix_B(SparseMatrix<Number>::build())
{
}

void CondensedEigenSystem::set_local_non_condensed_dofs(std::vector<unsigned int>& non_condensed_dofs_in)
{
  local_non_condensed_dofs_vector = non_condensed_dofs_in;
}


void CondensedEigenSystem::solve()
{
  START_LOG("solve()", "CondensedEigenSystem");

  // Make sure we have initialized the non-condensed DOFs
  libmesh_assert(!local_non_condensed_dofs_vector.empty());

  // Need to get a global copy of local_non_condensed_dofs_vector on all processors
  std::vector<unsigned int> global_non_condensed_dofs_vector = local_non_condensed_dofs_vector;
  Parallel::allgather(global_non_condensed_dofs_vector);

  // Now condense the matrices
  matrix_A->create_submatrix(*condensed_matrix_A,
                             local_non_condensed_dofs_vector,
                             global_non_condensed_dofs_vector);

  matrix_B->create_submatrix(*condensed_matrix_B,
                             local_non_condensed_dofs_vector,
                             global_non_condensed_dofs_vector);

  // A reference to the EquationSystems
  EquationSystems& es = this->get_equation_systems();

  // check that necessary parameters have been set
  libmesh_assert (es.parameters.have_parameter<unsigned int>("eigenpairs"));
  libmesh_assert (es.parameters.have_parameter<unsigned int>("basis vectors"));

  // Get the tolerance for the solver and the maximum
  // number of iterations. Here, we simply adopt the linear solver
  // specific parameters.
  const Real tol            =
    es.parameters.get<Real>("linear solver tolerance");

  const unsigned int maxits =
    es.parameters.get<unsigned int>("linear solver maximum iterations");

  const unsigned int nev    =
    es.parameters.get<unsigned int>("eigenpairs");

  const unsigned int ncv    =
    es.parameters.get<unsigned int>("basis vectors");

  std::pair<unsigned int, unsigned int> solve_data;

  // call the solver depending on the type of eigenproblem
  if ( generalized() )
    {
      //in case of a generalized eigenproblem
      solve_data = eigen_solver->solve_generalized
        (*condensed_matrix_A,*condensed_matrix_B, nev, ncv, tol, maxits);
    }

  else
    {
      libmesh_assert (matrix_B == NULL);

      //in case of a standard eigenproblem
      solve_data = eigen_solver->solve_standard (*condensed_matrix_A, nev, ncv, tol, maxits);
    }

  set_n_converged(solve_data.first);
  set_n_iterations(solve_data.second);

  STOP_LOG("solve()", "CondensedEigenSystem");
}



std::pair<Real, Real> CondensedEigenSystem::get_eigenpair(unsigned int i)
{
  START_LOG("get_eigenpair()", "CondensedEigenSystem");

  // Make sure we have initialized the non-condensed DOFs
  libmesh_assert(!local_non_condensed_dofs_vector.empty());

  // This function assumes that condensed_solve has just been called.
  // If this is not the case, then we will trip an asset in get_eigenpair
  AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build();
  unsigned int n_local = local_non_condensed_dofs_vector.size();
  unsigned int n       = n_local;
  Parallel::sum(n);

  temp->init (n, n_local, false, libMeshEnums::PARALLEL);

  std::pair<Real, Real> eval = eigen_solver->get_eigenpair (i, *temp);

  // Now map temp to solution. Loop over local entries of local_non_condensed_dofs_vector
  this->solution->zero();
  for(unsigned int i=0; i<local_non_condensed_dofs_vector.size(); i++)
  {
    unsigned int index = local_non_condensed_dofs_vector[i];
    solution->set(index,(*temp)(temp->first_local_index()+i));
  }

  solution->close();
  this->update();

  STOP_LOG("get_eigenpair()", "CondensedEigenSystem");

  return eval;
}




} // namespace libMesh


#endif // LIBMESH_HAVE_SLEPC