// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __tensor_value_h__
#define __tensor_value_h__

// C++ includes

// Local includes
#include "type_tensor.h"




/**
 * This class defines a tensor in DIM dimensional Real or Complex
 * space.  The typedef RealTensorValue always defines a real-valued tensor,
 * and NumberTensorValue defines a real or complex-valued tensor depending
 * on how the library was configured.
 * 
 * \author Roy H. Stogner 2004
 */
template <typename T>
class TensorValue : public TypeTensor<T>
{
public:

  /**
   * Constructor.  By default sets all entries to 0.
   */
  TensorValue  (const T xx=0.,
		const T xy=0.,
		const T xz=0.,
		const T yx=0.,
		const T yy=0.,
		const T yz=0.,
		const T zx=0.,
		const T zy=0.,
		const T zz=0.);
  
  /**
   * Copy-constructor.
   */
  TensorValue (const TensorValue<T>& p);

  /**
   * Copy-constructor.
   */
  TensorValue (const TypeTensor<T>& p);


#ifdef USE_COMPLEX_NUMBERS
  /**
   * Constructor that takes two \p TypeTensor<Real>
   * representing the real and imaginary part as
   * arguments.
   */
  TensorValue (const TypeTensor<Real>& p_re,
	       const TypeTensor<Real>& p_im);
#endif

  
private:

  
};



/**
 * Useful typedefs to allow transparent switching
 * between Real and Complex data types.
 */
typedef TensorValue<Real>   RealTensorValue;
typedef TensorValue<Number> NumberTensorValue;
typedef RealTensorValue     RealTensor;
typedef NumberTensorValue   Tensor;



//------------------------------------------------------
// Inline functions
template <typename T>
inline
TensorValue<T>::TensorValue (const T xx,
			     const T xy,
			     const T xz,
			     const T yx,
			     const T yy,
			     const T yz,
			     const T zx,
			     const T zy,
			     const T zz) :
  TypeTensor<T> (xx,xy,xz,yx,yy,yz,zx,zy,zz)
{
}



template <typename T>
inline
TensorValue<T>::TensorValue (const TensorValue<T>& p) :
  TypeTensor<T> (p)
{
}



template <typename T>
inline
TensorValue<T>::TensorValue (const TypeTensor<T>& p) :
  TypeTensor<T> (p)
{
}


#ifdef USE_COMPLEX_NUMBERS
template <typename T>
inline
TensorValue<T>::TensorValue (const TypeTensor<Real>& p_re,
			     const TypeTensor<Real>& p_im) :
  TypeTensor<T> (Complex (p_re(0,0), p_im(0,0)),
		 Complex (p_re(0,1), p_im(0,1)),
		 Complex (p_re(0,2), p_im(0,2)),
		 Complex (p_re(1,0), p_im(1,0)),
		 Complex (p_re(1,1), p_im(1,1)),
		 Complex (p_re(1,2), p_im(1,2)),
		 Complex (p_re(2,0), p_im(2,0)),
		 Complex (p_re(2,1), p_im(2,1)),
		 Complex (p_re(2,2), p_im(2,2)))
{
}
#endif


#endif // #define __tensor_value_h__
