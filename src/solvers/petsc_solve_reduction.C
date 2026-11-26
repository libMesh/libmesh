// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_solve_reduction.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/sparse_matrix.h"

namespace libMesh
{

//------------------------------------------------------------------
// PetscSolveReduction

template <typename T>
PetscSolveReduction<T>::PetscSolveReduction(NonlinearImplicitSystem & sys)
  : _sys(sys), _solution(nullptr), _snes_solution(nullptr), _snes_residual(nullptr)
{
}

template <typename T>
PetscSolveReduction<T>::~PetscSolveReduction() = default;

template <typename T>
void PetscSolveReduction<T>::jacobian(Vec y,
                                      SparseMatrix<T> & jac,
                                      SparseMatrix<T> & pc,
                                      const AssembleJacobian & assemble)
{
  this->prolong(y);
  assemble(pc);
  pc.close();
  this->finish_jacobian(pc);

  if (&jac != &pc)
    {
      this->check_jacobian_operator(jac);
      jac.close();
    }
}

//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscSolveReduction<Number>;

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
