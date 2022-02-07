// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local Includes
#include "libmesh/dense_matrix_impl.h"

namespace libMesh
{

//--------------------------------------------------------------
// Explicit instantiations
#define LIBMESH_VMA_INSTANTIATE(T1,T2,T3)               \
  template void DenseMatrix<T1>::vector_mult_add        \
  (DenseVector<                                         \
   CompareTypes<T1,                                     \
   CompareTypes<T2,T3>::supertype>::supertype> & dest,  \
   const T2 factor,                                     \
   const DenseVector<T3> & arg) const

template class DenseMatrix<Real>;
template void DenseMatrix<Real>::cholesky_solve(const DenseVector<Real> &, DenseVector<Real> &);
template void DenseMatrix<Real>::_cholesky_back_substitute(const DenseVector<Real> &, DenseVector<Real> &) const;
template void DenseMatrix<Real>::cholesky_solve(const DenseVector<Complex> &, DenseVector<Complex> &);
template void DenseMatrix<Real>::_cholesky_back_substitute(const DenseVector<Complex> &, DenseVector<Complex> &) const;
LIBMESH_VMA_INSTANTIATE(Real,int,Real);
#ifndef LIBMESH_DEFAULT_SINGLE_PRECISION
LIBMESH_VMA_INSTANTIATE(Real,float,Real);
#endif
#ifndef LIBMESH_DEFAULT_DOUBLE_PRECISION
LIBMESH_VMA_INSTANTIATE(Real,double,Real);
#endif

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template class DenseMatrix<Complex>;
template void DenseMatrix<Complex>::cholesky_solve(const DenseVector<Complex> &,DenseVector<Complex> &);
template void DenseMatrix<Complex>::_cholesky_back_substitute(const DenseVector<Complex> &, DenseVector<Complex> &) const;
template void DenseMatrix<Real>::vector_mult (DenseVector<CompareTypes<Real,Complex>::supertype> & dest,
                                              const DenseVector<Complex> & arg) const;
template void DenseMatrix<Real>::vector_mult_transpose (DenseVector<CompareTypes<Real,Complex>::supertype> & dest,
                                                        const DenseVector<Complex> & arg) const;
LIBMESH_VMA_INSTANTIATE(Real,int,Complex);
LIBMESH_VMA_INSTANTIATE(Complex,int,Complex);
LIBMESH_VMA_INSTANTIATE(Complex,int,Real);

// complex<int> and complex<float_foo> don't interact well
//LIBMESH_VMA_INSTANTIATE(Real,std::complex<int>,Complex);
//LIBMESH_VMA_INSTANTIATE(Complex,std::complex<int>,Complex);
//LIBMESH_VMA_INSTANTIATE(Complex,std::complex<int>,Real);

LIBMESH_VMA_INSTANTIATE(Real,float,Complex);
LIBMESH_VMA_INSTANTIATE(Complex,float,Complex);
LIBMESH_VMA_INSTANTIATE(Complex,float,Real);
LIBMESH_VMA_INSTANTIATE(Real,std::complex<float>,Complex);
#ifndef LIBMESH_DEFAULT_SINGLE_PRECISION
LIBMESH_VMA_INSTANTIATE(Complex,std::complex<float>,Complex);
#endif
LIBMESH_VMA_INSTANTIATE(Complex,std::complex<float>,Real);

LIBMESH_VMA_INSTANTIATE(Real,double,Complex);
LIBMESH_VMA_INSTANTIATE(Complex,double,Complex);
LIBMESH_VMA_INSTANTIATE(Complex,double,Real);
LIBMESH_VMA_INSTANTIATE(Real,std::complex<double>,Complex);
#ifndef LIBMESH_DEFAULT_DOUBLE_PRECISION
LIBMESH_VMA_INSTANTIATE(Complex,std::complex<double>,Complex);
#endif
LIBMESH_VMA_INSTANTIATE(Complex,std::complex<double>,Real);
#endif

} // namespace libMesh
