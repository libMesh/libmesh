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

#include "libmesh/dense_matrix_base_impl.h"
#include "libmesh/dense_matrix_impl.h"

namespace libMesh {
template <>
void DenseMatrix<GeomReal>::_multiply_blas(const DenseMatrixBase<GeomReal> &,
                                           _BLAS_Multiply_Flag) {
  libmesh_not_implemented();
}

template <> void DenseMatrix<GeomReal>::_lu_decompose_lapack() {
  libmesh_not_implemented();
}

template <> void DenseMatrix<GeomReal>::_svd_lapack(DenseVector<Real> &) {
  libmesh_not_implemented();
}

template <>
void DenseMatrix<GeomReal>::_svd_lapack(DenseVector<Real> &,
                                        DenseMatrix<Number> &,
                                        DenseMatrix<Number> &) {
  libmesh_not_implemented();
}

template <>
void DenseMatrix<GeomReal>::_svd_helper(char, char, std::vector<Real> &,
                                        std::vector<Number> &,
                                        std::vector<Number> &) {
  libmesh_not_implemented();
}

template <>
void DenseMatrix<GeomReal>::_svd_solve_lapack(
    const DenseVector<GeomReal> & /*rhs*/, DenseVector<GeomReal> & /*x*/,
    Real /*rcond*/) const {
  libmesh_not_implemented();
}

template <>
void DenseMatrix<GeomReal>::_evd_lapack(DenseVector<GeomReal> &,
                                        DenseVector<GeomReal> &,
                                        DenseMatrix<GeomReal> *,
                                        DenseMatrix<GeomReal> *) {
  libmesh_not_implemented();
}

template <>
void DenseMatrix<GeomReal>::_lu_back_substitute_lapack(
    const DenseVector<GeomReal> &, DenseVector<GeomReal> &) {
  libmesh_not_implemented();
}

template <>
void DenseMatrix<GeomReal>::_matvec_blas(GeomReal, GeomReal,
                                         DenseVector<GeomReal> &,
                                         const DenseVector<GeomReal> &,
                                         bool) const {
  libmesh_not_implemented();
}

template <>
DenseMatrix<GeomReal>::DenseMatrix(const unsigned int new_m,
                                   const unsigned int new_n)
    : DenseMatrixBase<GeomReal>(new_m, new_n), use_blas_lapack(false), _val(),
      _decomposition_type(NONE) {
  this->resize(new_m, new_n);
}

template class DenseMatrixBase<GeomReal>;
template class DenseMatrix<GeomReal>;

template void
DenseMatrix<GeomReal>::vector_mult_add(DenseVector<GeomReal> &, const int,
                                       const DenseVector<GeomReal> &) const;
template void
DenseMatrix<GeomReal>::vector_mult_add(DenseVector<GeomReal> &, const double,
                                       const DenseVector<GeomReal> &) const;

template void
DenseMatrix<Real>::vector_mult(DenseVector<GeomReal> &,
                               const DenseVector<GeomReal> &) const;
template void
DenseMatrix<GeomReal>::vector_mult(DenseVector<GeomReal> &,
                                   const DenseVector<Real> &) const;

template void
DenseMatrix<GeomReal>::cholesky_solve(const DenseVector<GeomReal> &b,
                                      DenseVector<GeomReal> &x);
} // namespace libMesh
