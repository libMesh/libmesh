// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __raw_accessor_h__
#define __raw_accessor_h__

namespace libMesh
{
  // Local includes
#include "libmesh_common.h"

  /**
   * This class provides single index access to FieldType (i.e. Number, Gradient, Tensor, etc.).
   */
  template <typename FieldType>
  class RawAccessor
  {
  public:
    
    RawAccessor( FieldType& data, const unsigned int dim )
      : _data(data),
	_dim(dim)
    {}

    ~RawAccessor(){};

    Number& operator()( unsigned int i );
    const Number& operator()( unsigned int i ) const;

  private:
    RawAccessor();

    FieldType& _data;
    const unsigned int _dim;
  };

  // Specialize for specific cases
  template<>
  Number& RawAccessor<Number>::operator()( unsigned int i )
  {
    libmesh_assert(i == 0);
    return this->_data;
  }

  template<>
  Number& RawAccessor<Gradient>::operator()( unsigned int i )
  {
    libmesh_assert(i < this->_dim);
    return this->_data(i);
  }

  template<>
  Number& RawAccessor<Tensor>::operator()( unsigned int k )
  {
    libmesh_assert(k < this->_dim*this->_dim);

    // For tensors, each row is filled first, i.e. for 2-D
    // [ 0 1; 2 3]
    // Thus, k(i,j) = i + j*dim
    unsigned int jj = k/_dim;
    unsigned int ii = k - jj*_dim;

    return this->_data(ii,jj);
  }

}
#endif //__raw_accessor_h__
