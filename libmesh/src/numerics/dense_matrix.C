// $Id: dense_matrix.C,v 1.6 2003-02-03 03:51:50 ddreyer Exp $

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
template<typename Tp>
void DenseMatrix<Tp>::left_multiply (const DenseMatrix<Tp>& A,
				     const bool transpose)
{
  // C = A*B  C (mxn), A (mxp), B (pxn)
  
  DenseMatrix<Tp> B(*this);

  DenseMatrix<Tp>& C = *this;

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

      for (unsigned int i=0; i<m_s; i++)
	for (unsigned int j=0; j<n_s; j++)
	  for (unsigned int k=0; k<p_s; k++)
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

      for (unsigned int i=0; i<m_s; i++)
	for (unsigned int j=0; j<n_s; j++)
	  for (unsigned int k=0; k<p_s; k++)
	    C(i,j) += A.transpose(i,k)*B(k,j);

      /*
      std::cout << "C=" << std::endl;
      C.print();
      */
    };    	     
};



template<typename Tp>
void DenseMatrix<Tp>::right_multiply (const DenseMatrix<Tp>& B,
				      const bool transpose)
{
  // C = A*B  C (mxn), A (mxp), B (pxn)
  
  DenseMatrix<Tp> A(*this);
  
  DenseMatrix<Tp>& C = *this;
  
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

      for (unsigned int i=0; i<m_s; i++)
	for (unsigned int j=0; j<n_s; j++)
	  for (unsigned int k=0; k<p_s; k++)
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

      for (unsigned int i=0; i<m_s; i++)
	for (unsigned int j=0; j<n_s; j++)
	  for (unsigned int k=0; k<p_s; k++)
	    C(i,j) += A(i,k)*B.transpose(k,j);	           
    };    	     
};


//--------------------------------------------------------------
// Explicit instantiations
template class DenseMatrix<Complex>;
#ifdef USE_COMPLEX_NUMBERS
 /* Avoid double instantiation in case of real-only */
 template class DenseMatrix<Real>;
#endif
