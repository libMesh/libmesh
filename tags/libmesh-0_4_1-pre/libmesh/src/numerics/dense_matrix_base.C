// $Id: dense_matrix_base.C,v 1.5 2003-09-02 18:02:43 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
#include "dense_matrix_base.h"



template<typename T>
void DenseMatrixBase<T>::multiply (DenseMatrixBase<T>& M1,
				   const DenseMatrixBase<T>& M2,
				   const DenseMatrixBase<T>& M3)
{
  // Assertions to make sure we have been
  // passed matrices of the correct dimension.
  assert (M1.m() == M2.m());
  assert (M1.n() == M3.n());
  assert (M2.n() == M3.m());
	  
  const unsigned int m_s = M2.m();
  const unsigned int p_s = M2.n(); 
  const unsigned int n_s = M1.n();
  
  // Do it this way because there is a
  // decent chance (at least for constraint matrices)
  // that M3(k,j) = 0. when right-multiplying.
  for (unsigned int k=0; k<p_s; k++)
    for (unsigned int j=0; j<n_s; j++)
      if (M3.el(k,j) != 0.)
	for (unsigned int i=0; i<m_s; i++)
	  M1.el(i,j) += M2.el(i,k) * M3.el(k,j);	          
}











//--------------------------------------------------------------
// Explicit instantiations
template class DenseMatrixBase<Real>;

#ifdef USE_COMPLEX_NUMBERS
template class DenseMatrixBase<Complex>;
#endif
