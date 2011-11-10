// $Id$

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


#ifndef __parallel_algebra_h__
#define __parallel_algebra_h__

// This class contains all the functionality for bin sorting
// Templated on the type of keys you will be sorting and the
// type of iterator you will be using.


#include "libmesh_config.h"

#include "parallel.h"
#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"

namespace libMesh {
namespace Parallel {
  // StandardType<> specializations to return a derived MPI datatype
  // to handle communication of LIBMESH_DIM-vectors.  We use a
  // singleton pattern here because a global variable would have tried
  // to call MPI functions before MPI got initialized.
  template <typename T>
  class StandardType<TypeVector<T> > : public DataType
  {
  public:
    inline StandardType(const TypeVector<T> *example=NULL) :
      DataType(StandardType<T>(example ?  &((*example)(0)) : NULL), LIBMESH_DIM) {}

    inline ~StandardType() { this->free(); }
  };

  template <typename T>
  class StandardType<VectorValue<T> > : public DataType
  {
  public:
    inline StandardType(const VectorValue<T> *example=NULL) :
      DataType(StandardType<T>(example ? &((*example)(0)) : NULL), LIBMESH_DIM) {}

    inline ~StandardType() { this->free(); }
  };

  template <>
  class StandardType<Point> : public DataType
  {
  public:
    inline StandardType(const Point *example=NULL) :
      DataType(StandardType<Real>(example ? &((*example)(0)) : NULL), LIBMESH_DIM) {}

    inline ~StandardType() { this->free(); }
  };

  // StandardType<> specializations to return a derived MPI datatype
  // to handle communication of LIBMESH_DIM*LIBMESH_DIM-tensors.  We use a
  // singleton pattern here because a global variable would have tried
  // to call MPI functions before MPI got initialized.
  template <typename T>
  class StandardType<TypeTensor<T> > : public DataType
  {
  public:
    inline StandardType(const TypeTensor<T> *example=NULL) :
      DataType(StandardType<T>(example ?  &((*example)(0,0)) : NULL), LIBMESH_DIM*LIBMESH_DIM) {}

    inline ~StandardType() { this->free(); }
  };

  template <typename T>
  class StandardType<TensorValue<T> > : public DataType
  {
  public:
    inline StandardType(const TensorValue<T> *example=NULL) :
      DataType(StandardType<T>(example ?  &((*example)(0,0)) : NULL), LIBMESH_DIM*LIBMESH_DIM) {}

    inline ~StandardType() { this->free(); }
  };
} // namespace Parallel
} // namespace libMesh

#endif // __parallel_algebra_h__

