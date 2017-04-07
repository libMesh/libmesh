// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
  StandardType(const TypeVector<T> * example=libmesh_nullptr) {
    // We need an example for MPI_Address to use
    TypeVector<T> * ex;
    UniquePtr<TypeVector<T> > temp;
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

        int blocklengths[3] = {1, LIBMESH_DIM, 1};
        MPI_Aint displs[3];
        MPI_Datatype types[3] = {MPI_LB, T_type, MPI_UB};
        MPI_Aint start, later;

        libmesh_call_mpi
          (MPI_Address(ex, &start));
        displs[0] = 0;
        libmesh_call_mpi
          (MPI_Address(&((*ex)(0)), &later));
        displs[1] = later - start;
        libmesh_call_mpi
          (MPI_Address((ex+1), &later));
        displs[2] = later - start;

        libmesh_call_mpi
          (MPI_Type_struct (3, blocklengths, displs, types,
                            &_static_type));

#else // MPI_VERSION >= 2

        int blocklength = LIBMESH_DIM;
        MPI_Aint displs, start;
        MPI_Datatype tmptype, type = T_type;

        libmesh_call_mpi
          (MPI_Get_address (ex, &start));
        libmesh_call_mpi
          (MPI_Get_address (&((*ex)(0)), &displs));

        // subtract off offset to first value from the beginning of the structure
        displs -= start;

        // create a prototype structure
        libmesh_call_mpi
          (MPI_Type_create_struct (1, &blocklength, &displs, &type,
                                   &tmptype));
        libmesh_call_mpi
          (MPI_Type_commit (&tmptype));

        // resize the structure type to account for padding, if any
        libmesh_call_mpi
          (MPI_Type_create_resized (tmptype, 0, sizeof(TypeVector<T>),
                                    &_static_type));
#endif

        libmesh_call_mpi
          (MPI_Type_commit (&_static_type));
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
  StandardType(const VectorValue<T> * example=libmesh_nullptr) {
    // We need an example for MPI_Address to use
    VectorValue<T> * ex;
    UniquePtr<VectorValue<T> > temp;
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

        int blocklengths[3] = {1, LIBMESH_DIM, 1};
        MPI_Aint displs[3];
        MPI_Datatype types[3] = {MPI_LB, T_type, MPI_UB};
        MPI_Aint start, later;

        libmesh_call_mpi
          (MPI_Address(ex, &start));
        displs[0] = 0;
        libmesh_call_mpi
          (MPI_Address(&((*ex)(0)), &later));
        displs[1] = later - start;
        libmesh_call_mpi
          (MPI_Address((ex+1), &later));
        displs[2] = later - start;

        libmesh_call_mpi
          (MPI_Type_struct (3, blocklengths, displs, types,
                            &_static_type));

#else // MPI_VERSION >= 2

        int blocklength = LIBMESH_DIM;
        MPI_Aint displs, start;
        MPI_Datatype tmptype, type = T_type;

        libmesh_call_mpi
          (MPI_Get_address (ex, &start));
        libmesh_call_mpi
          (MPI_Get_address (&((*ex)(0)), &displs));

        // subtract off offset to first value from the beginning of the structure
        displs -= start;

        // create a prototype structure
        libmesh_call_mpi
          (MPI_Type_create_struct (1, &blocklength, &displs, &type,
                                   &tmptype));
        libmesh_call_mpi
          (MPI_Type_commit (&tmptype));

        // resize the structure type to account for padding, if any
        libmesh_call_mpi
          (MPI_Type_create_resized (tmptype, 0,
                                    sizeof(VectorValue<T>),
                                    &_static_type));
#endif

        libmesh_call_mpi
          (MPI_Type_commit (&_static_type));
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
  StandardType(const Point * example=libmesh_nullptr)
  {
    // Prevent unused variable warnings when !LIBMESH_HAVE_MPI
    libmesh_ignore(example);

    // _static_type never gets freed, but it only gets committed once
    // per T, so it's not a *huge* memory leak...
    static data_type _static_type;
    static bool _is_initialized = false;
    if (!_is_initialized)
      {
#ifdef LIBMESH_HAVE_MPI

        // We need an example for MPI_Address to use
        Point * ex;

        UniquePtr<Point> temp;
        if (example)
          ex = const_cast<Point *>(example);
        else
          {
            temp.reset(new Point());
            ex = temp.get();
          }

        StandardType<Real> T_type(&((*ex)(0)));

#if MPI_VERSION == 1

        int blocklengths[3] = {1, LIBMESH_DIM, 1};
        MPI_Aint displs[3];
        MPI_Datatype types[3] = {MPI_LB, T_type, MPI_UB};
        MPI_Aint start, later;

        libmesh_call_mpi
          (MPI_Address(ex, &start));
        displs[0] = 0;
        libmesh_call_mpi
          (MPI_Address(&((*ex)(0)), &later));
        displs[1] = later - start;
        libmesh_call_mpi
          (MPI_Address((ex+1), &later));
        displs[2] = later - start;

        libmesh_call_mpi
          (MPI_Type_struct (3, blocklengths, displs, types,
                            &_static_type));

#else // MPI_VERSION >= 2

        int blocklength = LIBMESH_DIM;
        MPI_Aint displs, start;
        MPI_Datatype tmptype, type = T_type;

        libmesh_call_mpi
          (MPI_Get_address (ex, &start));
        libmesh_call_mpi
          (MPI_Get_address (&((*ex)(0)), &displs));

        // subtract off offset to first value from the beginning of the structure
        displs -= start;

        // create a prototype structure
        libmesh_call_mpi
          (MPI_Type_create_struct (1, &blocklength, &displs, &type,
                                   &tmptype));
        libmesh_call_mpi
          (MPI_Type_commit (&tmptype));

        // resize the structure type to account for padding, if any
        libmesh_call_mpi
          (MPI_Type_create_resized (tmptype, 0, sizeof(Point),
                                    &_static_type));
#endif

        libmesh_call_mpi
          (MPI_Type_commit (&_static_type));
#endif // #ifdef LIBMESH_HAVE_MPI

        _is_initialized = true;
      }
    _datatype = _static_type;
  }
};

// OpFunction<> specializations to return an MPI_Op version of the
// reduction operations on LIBMESH_DIM-vectors.
//
// We use static variables to minimize the number of MPI datatype
// construction calls executed over the course of the program.
//
// We use a singleton pattern because a global variable would
// have tried to call MPI functions before MPI got initialized.
//
// min() and max() are applied component-wise; this makes them useful
// for bounding box reduction operations.
template <typename V>
class TypeVectorOpFunction
{
public:
#ifdef LIBMESH_HAVE_MPI
  static void vector_max (void *invec, void *inoutvec, int *len, MPI_Datatype *)
  {
    V *in = static_cast<V *>(invec);
    V *inout = static_cast<V *>(inoutvec);
    for (int i=0; i != *len; ++i)
      for (int d=0; d != LIBMESH_DIM; ++d)
        inout[i](d) = std::max(in[i](d), inout[i](d));
  }

  static void vector_min (void *invec, void *inoutvec, int *len, MPI_Datatype *)
  {
    V *in = static_cast<V *>(invec);
    V *inout = static_cast<V *>(inoutvec);
    for (int i=0; i != *len; ++i)
      for (int d=0; d != LIBMESH_DIM; ++d)
        inout[i](d) = std::min(in[i](d), inout[i](d));
  }

  static void vector_sum (void *invec, void *inoutvec, int *len, MPI_Datatype *)
  {
    V *in = static_cast<V *>(invec);
    V *inout = static_cast<V *>(inoutvec);
    for (int i=0; i != *len; ++i)
      for (int d=0; d != LIBMESH_DIM; ++d)
        inout[i](d) += in[i](d);
  }

  static MPI_Op max()
  {
    // _static_op never gets freed, but it only gets committed once
    // per T, so it's not a *huge* memory leak...
    static MPI_Op _static_op;
    static bool _is_initialized = false;
    if (!_is_initialized)
      {
        libmesh_call_mpi
          (MPI_Op_create (vector_max, /*commute=*/ true,
                          &_static_op));

        _is_initialized = true;
      }

    return _static_op;
  }
  static MPI_Op min()
  {
    // _static_op never gets freed, but it only gets committed once
    // per T, so it's not a *huge* memory leak...
    static MPI_Op _static_op;
    static bool _is_initialized = false;
    if (!_is_initialized)
      {
        libmesh_call_mpi
          (MPI_Op_create (vector_min, /*commute=*/ true,
                          &_static_op));

        _is_initialized = true;
      }

    return _static_op;
  }
  static MPI_Op sum()
  {
    // _static_op never gets freed, but it only gets committed once
    // per T, so it's not a *huge* memory leak...
    static MPI_Op _static_op;
    static bool _is_initialized = false;
    if (!_is_initialized)
      {
        libmesh_call_mpi
          (MPI_Op_create (vector_sum, /*commute=*/ true,
                          &_static_op));

        _is_initialized = true;
      }

    return _static_op;
  }

#endif // LIBMESH_HAVE_MPI
};

template <typename T>
class OpFunction<TypeVector<T> > : public TypeVectorOpFunction<TypeVector<T> > {};

template <typename T>
class OpFunction<VectorValue<T> > : public TypeVectorOpFunction<VectorValue<T> > {};

template <>
class OpFunction<Point> : public TypeVectorOpFunction<Point> {};

// StandardType<> specializations to return a derived MPI datatype
// to handle communication of LIBMESH_DIM*LIBMESH_DIM-tensors.
//
// We assume contiguous storage here
template <typename T>
class StandardType<TypeTensor<T> > : public DataType
{
public:
  explicit
  StandardType(const TypeTensor<T> * example=libmesh_nullptr) :
    DataType(StandardType<T>(example ?  &((*example)(0,0)) : libmesh_nullptr), LIBMESH_DIM*LIBMESH_DIM) {}

  inline ~StandardType() { this->free(); }
};

template <typename T>
class StandardType<TensorValue<T> > : public DataType
{
public:
  explicit
  StandardType(const TensorValue<T> * example=libmesh_nullptr) :
    DataType(StandardType<T>(example ?  &((*example)(0,0)) : libmesh_nullptr), LIBMESH_DIM*LIBMESH_DIM) {}

  inline ~StandardType() { this->free(); }
};
} // namespace Parallel
} // namespace libMesh

#endif // LIBMESH_PARALLEL_ALGEBRA_H
