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


#ifndef LIBMESH_PARALLEL_ALGEBRA_H
#define LIBMESH_PARALLEL_ALGEBRA_H

// This class contains all the functionality for bin sorting
// Templated on the type of keys you will be sorting and the
// type of iterator you will be using.


// Local Includes
#include "libmesh/libmesh_config.h"

#include "libmesh/auto_ptr.h"
#include "libmesh/parallel.h"
#include "libmesh/point.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

// C++ includes
#include <cstddef>

namespace libMesh {
namespace Parallel {
  // StandardType<> specializations to return a derived MPI datatype
  // to handle communication of LIBMESH_DIM-vectors.
  //
  // We use static variables to minimize the number of MPI datatype
  // construction calls executed over the course of the program.
  //
  // We use a singleton pattern because a global variable would
  // have tried to call MPI functions before MPI got initialized.
  //
  // We use MPI_Create_struct here because our vector classes might
  // have vptrs, and I'd rather not have the datatype break in those
  // cases.
  template <typename T>
  class StandardType<TypeVector<T> > : public DataType
  {
  public:
    explicit
    StandardType(const TypeVector<T> *example=NULL) {
      // We need an example for MPI_Address to use
      TypeVector<T> *ex;
      AutoPtr<TypeVector<T> > temp;
      if (example)
        ex = const_cast<TypeVector<T> *>(example);
      else
        {
          temp.reset(new TypeVector<T>());
          ex = temp.get();
        }

      // _static_type never gets freed, but it only gets committed once
      // per T, so it's not a *huge* memory leak...
      static data_type _static_type;
      static bool _is_initialized = false;
      if (!_is_initialized)
        {
#ifdef LIBMESH_HAVE_MPI
          StandardType<T> T_type(&((*ex)(0)));

#if MPI_VERSION == 1

          int blocklengths[LIBMESH_DIM+2];
          MPI_Aint displs[LIBMESH_DIM+2];
          MPI_Datatype types[LIBMESH_DIM+2];
          MPI_Aint start, later;

          MPI_Address(ex, &start);
          blocklengths[0] = 1;
          displs[0] = 0;
          types[0] = MPI_LB;
          for (unsigned int i=0; i != LIBMESH_DIM; ++i)
            {
              MPI_Address(&((*ex)(i)), &later);
              blocklengths[i+1] = 1;
              displs[i+1] = later - start;
              types[i+1] = T_type;
            }
          MPI_Address((ex+1), &later);
          blocklengths[LIBMESH_DIM+1] = 1;
          displs[LIBMESH_DIM+1] = later - start;
          types[LIBMESH_DIM+1] = MPI_UB;

          MPI_Type_struct (LIBMESH_DIM+2, blocklengths, displs, types, &_static_type);

#else // MPI_VERSION >= 2

          int blocklength = LIBMESH_DIM;
          MPI_Aint displs, start;
          MPI_Datatype tmptype, type = T_type;

	  MPI_Get_address (ex,   &start);
	  MPI_Get_address (&((*ex)(0)), &displs);

	  // subtract off offset to first value from the beginning of the structure
	  displs -= start;

	  // create a prototype structure
	  MPI_Type_create_struct (1, &blocklength, &displs, &type, &tmptype);

	  // resize the structure type to account for padding, if any
	  MPI_Type_create_resized (tmptype, 0, sizeof(TypeVector<T>), &_static_type);
#endif

          MPI_Type_commit (&_static_type);
#endif // #ifdef LIBMESH_HAVE_MPI

          _is_initialized = true;
        }
      _datatype = _static_type;
    }
  };

  template <typename T>
  class StandardType<VectorValue<T> > : public DataType
  {
  public:
    explicit
    StandardType(const VectorValue<T> *example=NULL) {
      // We need an example for MPI_Address to use
      VectorValue<T> *ex;
      AutoPtr<VectorValue<T> > temp;
      if (example)
        ex = const_cast<VectorValue<T> *>(example);
      else
        {
          temp.reset(new VectorValue<T>());
          ex = temp.get();
        }

      // _static_type never gets freed, but it only gets committed once
      // per T, so it's not a *huge* memory leak...
      static data_type _static_type;
      static bool _is_initialized = false;
      if (!_is_initialized)
        {
#ifdef LIBMESH_HAVE_MPI
          StandardType<T> T_type(&((*ex)(0)));

#if MPI_VERSION == 1

          int blocklengths[LIBMESH_DIM+2];
          MPI_Aint displs[LIBMESH_DIM+2];
          MPI_Datatype types[LIBMESH_DIM+2];
          MPI_Aint start, later;

          MPI_Address(ex, &start);
          blocklengths[0] = 1;
          displs[0] = 0;
          types[0] = MPI_LB;
          for (unsigned int i=0; i != LIBMESH_DIM; ++i)
            {
              MPI_Address(&((*ex)(i)), &later);
              blocklengths[i+1] = 1;
              displs[i+1] = later - start;
              types[i+1] = T_type;
            }
          MPI_Address((ex+1), &later);
          blocklengths[LIBMESH_DIM+1] = 1;
          displs[LIBMESH_DIM+1] = later - start;
          types[LIBMESH_DIM+1] = MPI_UB;

          MPI_Type_struct (LIBMESH_DIM+2, blocklengths, displs, types, &_static_type);

#else // MPI_VERSION >= 2

          int blocklength = LIBMESH_DIM;
          MPI_Aint displs, start;
          MPI_Datatype tmptype, type = T_type;

	  MPI_Get_address (ex,   &start);
	  MPI_Get_address (&((*ex)(0)), &displs);

	  // subtract off offset to first value from the beginning of the structure
	  displs -= start;

	  // create a prototype structure
	  MPI_Type_create_struct (1, &blocklength, &displs, &type, &tmptype);

	  // resize the structure type to account for padding, if any
	  MPI_Type_create_resized (tmptype, 0, sizeof(VectorValue<T>), &_static_type);
#endif

          MPI_Type_commit (&_static_type);
#endif // #ifdef LIBMESH_HAVE_MPI

          _is_initialized = true;
        }
      _datatype = _static_type;
    }
  };

  template <>
  class StandardType<Point> : public DataType
  {
  public:
    explicit
    StandardType(const Point *example=NULL) {
      // We need an example for MPI_Address to use
      Point *ex;
      AutoPtr<Point> temp;
      if (example)
        ex = const_cast<Point *>(example);
      else
        {
          temp.reset(new Point());
          ex = temp.get();
        }

      // _static_type never gets freed, but it only gets committed once
      // per T, so it's not a *huge* memory leak...
      static data_type _static_type;
      static bool _is_initialized = false;
      if (!_is_initialized)
        {
#ifdef LIBMESH_HAVE_MPI
          StandardType<Real> T_type(&((*ex)(0)));

#if MPI_VERSION == 1

          int blocklengths[LIBMESH_DIM+2];
          MPI_Aint displs[LIBMESH_DIM+2];
          MPI_Datatype types[LIBMESH_DIM+2];
          MPI_Aint start, later;

          MPI_Address(ex, &start);
          blocklengths[0] = 1;
          displs[0] = 0;
          types[0] = MPI_LB;
          for (unsigned int i=0; i != LIBMESH_DIM; ++i)
            {
              MPI_Address(&((*ex)(i)), &later);
              blocklengths[i+1] = 1;
              displs[i+1] = later - start;
              types[i+1] = T_type;
            }
          MPI_Address((ex+1), &later);
          blocklengths[LIBMESH_DIM+1] = 1;
          displs[LIBMESH_DIM+1] = later - start;
          types[LIBMESH_DIM+1] = MPI_UB;

          MPI_Type_struct (LIBMESH_DIM+2, blocklengths, displs, types, &_static_type);

#else // MPI_VERSION >= 2

          int blocklength = LIBMESH_DIM;
          MPI_Aint displs, start;
          MPI_Datatype tmptype, type = T_type;

	  MPI_Get_address (ex,   &start);
	  MPI_Get_address (&((*ex)(0)), &displs);

	  // subtract off offset to first value from the beginning of the structure
	  displs -= start;

	  // create a prototype structure
	  MPI_Type_create_struct (1, &blocklength, &displs, &type, &tmptype);

	  // resize the structure type to account for padding, if any
	  MPI_Type_create_resized (tmptype, 0, sizeof(Point), &_static_type);
#endif

          MPI_Type_commit (&_static_type);
#endif // #ifdef LIBMESH_HAVE_MPI

          _is_initialized = true;
        }
      _datatype = _static_type;
    }
  };

  // StandardType<> specializations to return a derived MPI datatype
  // to handle communication of LIBMESH_DIM*LIBMESH_DIM-tensors.
  //
  // We use a singleton pattern here because a global variable would
  // have tried to call MPI functions before MPI got initialized.
  //
  // We assume contiguous storage here
  template <typename T>
  class StandardType<TypeTensor<T> > : public DataType
  {
  public:
    explicit
    StandardType(const TypeTensor<T> *example=NULL) :
      DataType(StandardType<T>(example ?  &((*example)(0,0)) : NULL), LIBMESH_DIM*LIBMESH_DIM) {}

    inline ~StandardType() { this->free(); }
  };

  template <typename T>
  class StandardType<TensorValue<T> > : public DataType
  {
  public:
    explicit
    StandardType(const TensorValue<T> *example=NULL) :
      DataType(StandardType<T>(example ?  &((*example)(0,0)) : NULL), LIBMESH_DIM*LIBMESH_DIM) {}

    inline ~StandardType() { this->free(); }
  };
} // namespace Parallel
} // namespace libMesh

#endif // LIBMESH_PARALLEL_ALGEBRA_H
