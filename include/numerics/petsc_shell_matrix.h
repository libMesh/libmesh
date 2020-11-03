// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PETSC_SHELL_MATRIX_H
#define LIBMESH_PETSC_SHELL_MATRIX_H


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_solver_exception.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/wrapped_petsc.h"

// Petsc include files.
#ifdef I
# define LIBMESH_SAW_I
#endif
#include <petscmat.h>
#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif

namespace libMesh
{

/**
 * This class allows to use a PETSc shell matrix.
 * All overridden virtual functions are documented in
 * shell_matrix.h.
 *
 * \author Fande Kong (fdkong.jd@gmail.com)
 * \date 2019
 */
template <typename T>
class PetscShellMatrix : public ShellMatrix<T>
{
public:

  PetscShellMatrix (const Parallel::Communicator & comm_in);

  virtual ~PetscShellMatrix () = default;

  virtual numeric_index_type m () const override;

  virtual numeric_index_type n () const override;

  virtual numeric_index_type local_m () const;

  virtual numeric_index_type local_n () const;

  virtual void vector_mult (NumericVector<T> & dest,
                            const NumericVector<T> & arg) const override;

  virtual void vector_mult_add (NumericVector<T> & dest,
                                const NumericVector<T> & arg) const override;

  virtual void get_diagonal (NumericVector<T> & dest) const override;

  virtual void clear () override;

  virtual void init () override;

  /**
   * \returns \p true if the matrix has been initialized,
   * \p false otherwise.
   */
  virtual bool initialized() const;

  /**
   * Returns a pointer to the underlying PETSc Mat object. Must call
   * init() before this.
   */
  Mat mat();

protected:

  /**
   * Petsc Shell Matrix
   */
  WrappedPetsc<Mat> _mat;

  bool _is_initialized;
};


//-----------------------------------------------------------------------
// PetscShellMatrix inline members
template <typename T>
inline
PetscShellMatrix<T>::PetscShellMatrix (const Parallel::Communicator & comm_in):
  ShellMatrix<T>(comm_in),
  _is_initialized(false)
{}



template <typename T>
inline
numeric_index_type PetscShellMatrix<T>::m () const
{
  PetscErrorCode ierr;
  PetscInt m;

  ierr = MatGetSize(_mat, &m, nullptr);
  LIBMESH_CHKERR(ierr);

  return m;
}



template <typename T>
inline
numeric_index_type PetscShellMatrix<T>::n () const
{
  PetscErrorCode ierr;
  PetscInt n;

  ierr = MatGetSize(_mat, nullptr, &n);
  LIBMESH_CHKERR(ierr);

  return n;
}


template <typename T>
inline
numeric_index_type PetscShellMatrix<T>::local_m () const
{
  PetscErrorCode ierr;
  PetscInt m;

  ierr = MatGetLocalSize(_mat, &m, nullptr);
  LIBMESH_CHKERR(ierr);

  return m;
}



template <typename T>
inline
numeric_index_type PetscShellMatrix<T>::local_n () const
{
  PetscErrorCode ierr;
  PetscInt n;

  ierr = MatGetLocalSize(_mat, nullptr, &n);
  LIBMESH_CHKERR(ierr);

  return n;
}


template <typename T>
inline
void PetscShellMatrix<T>::get_diagonal (NumericVector<T> & dest) const
{
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T> & petsc_dest = cast_ref<PetscVector<T> &>(dest);

  PetscErrorCode ierr = MatGetDiagonal(_mat, petsc_dest.vec());
  LIBMESH_CHKERR(ierr);
}



} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_SPARSE_SHELL_MATRIX_H
