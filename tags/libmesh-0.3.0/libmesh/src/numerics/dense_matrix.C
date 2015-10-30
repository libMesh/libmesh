// $Id: dense_matrix.C,v 1.10 2003-02-21 21:03:56 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "dense_matrix.h"




// ------------------------------------------------------------
// Dense Matrix member functions
template<typename T>
void DenseMatrix<T>::left_multiply (const DenseMatrix<T>& A,
				    const bool transpose)
{
  // C = A*B  C (mxn), A (mxp), B (pxn)
  
  DenseMatrix<T> B(*this);

  DenseMatrix<T>& C = *this;

  C.zero();
  
  if (!transpose)
    {
      resize (A.m(), B.n());
      
      assert (A.n() == B.m());
      assert (m() == A.m());
      assert (n() == B.n());
      
      const unsigned int m_s = A.m();
      const unsigned int n_s = C.n();
      const unsigned int p_s = A.n(); 

      // Do it this way because there is a
      // decent chance (at least for constraint matrices)
      // that A(i,k) = 0.
      for (unsigned int i=0; i<m_s; i++)
	for (unsigned int k=0; k<p_s; k++)
	  if (A(i,k) != 0.)
	    for (unsigned int j=0; j<n_s; j++)
	      C(i,j) += A(i,k)*B(k,j);	           
    }
  else
    {
      /*
	std::cout << "A=" << std::endl;
	A.print();
	std::cout << "B=" << std::endl;
	B.print();
      */
      
      resize (A.n(), B.n());
      
      assert (A.m() == B.m());
      assert (m() == A.n());
      assert (n() == B.n());
      
      const unsigned int m_s = A.n();
      const unsigned int n_s = C.n();
      const unsigned int p_s = A.m(); 

      // Do it this way because there is a
      // decent chance (at least for constraint matrices)
      // that A.transpose(i,k) = 0.
      for (unsigned int i=0; i<m_s; i++)
	for (unsigned int k=0; k<p_s; k++)
	  if (A.transpose(i,k) != 0.)
	    for (unsigned int j=0; j<n_s; j++)
	      C(i,j) += A.transpose(i,k)*B(k,j);

      /*
	std::cout << "C=" << std::endl;
	C.print();
      */
    }    	     
}



template<typename T>
void DenseMatrix<T>::right_multiply (const DenseMatrix<T>& B,
				     const bool transpose)
{
  // C = A*B  C (mxn), A (mxp), B (pxn)
  
  DenseMatrix<T> A(*this);
  
  DenseMatrix<T>& C = *this;
  
  C.zero();
  
  if (!transpose)
    {
      /*
	std::cout << "A=" << std::endl;
	A.print();
	std::cout << "B=" << std::endl;
	B.print();
      */
      resize (A.m(), B.n());
      
      assert (A.n() == B.m());
      assert (m() == A.m());
      assert (n() == B.n());
      
      const unsigned int m_s = A.m();
      const unsigned int n_s = C.n();
      const unsigned int p_s = A.n(); 

      // Do it this way because there is a
      // decent chance (at least for constraint matrices)
      // that B(k,j) = 0.
      for (unsigned int j=0; j<n_s; j++)
	for (unsigned int k=0; k<p_s; k++)
	  if (B(k,j) != 0.)
	    for (unsigned int i=0; i<m_s; i++)
	      C(i,j) += A(i,k)*B(k,j);

      /*
	std::cout << "C=" << std::endl;
	C.print();
      */
    }
  else
    {
      resize (A.m(), B.m());
      
      assert (A.n() == B.n());
      assert (m() == A.m());
      assert (n() == B.m());
      
      const unsigned int m_s = A.m();
      const unsigned int n_s = C.n();
      const unsigned int p_s = A.n(); 

      // Do it this way because there is a
      // decent chance (at least for constraint matrices)
      // that B.transpose(k,j) = 0.
      for (unsigned int j=0; j<n_s; j++)
	for (unsigned int k=0; k<p_s; k++)
	  if (B.transpose(k,j) != 0.)
	    for (unsigned int i=0; i<m_s; i++)
	      C(i,j) += A(i,k)*B.transpose(k,j);	           
    }    	     
}


//--------------------------------------------------------------
// Explicit instantiations
template class DenseMatrix<Real>;
template class DenseMatrix<Complex>;
