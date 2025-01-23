// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_PARALLEL_FE_TYPE_H
#define LIBMESH_PARALLEL_FE_TYPE_H


// libMesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/fe_type.h"

// TIMPI includes
#include "timpi/attributes.h"
#include "timpi/op_function.h"
#include "timpi/standard_type.h"

// C++ includes
#include <cstddef>
#include <memory>
#include <type_traits>

namespace TIMPI {

TIMPI_STANDARD_TYPE(libMesh::Order,MPI_INT);
TIMPI_STANDARD_TYPE(libMesh::FEFamily,MPI_INT);
TIMPI_PARALLEL_INTEGER_OPS(libMesh::Order);

// OpFunction<> specializations to return an MPI_Op version of the
// reduction operations on FETypes
//
// We use static variables to minimize the number of MPI datatype
// construction calls executed over the course of the program.
//
// We use a singleton pattern because a global variable would
// have tried to call MPI functions before MPI got initialized.
//
// Doing this to the FEFamily part of a type doesn't make much sense,
// *except* that it's how we verify that a requested family matches
// across multiple processors.

template <>
class OpFunction<libMesh::FEType>
{
public:
#ifdef LIBMESH_HAVE_MPI
  static void fetype_max (void * intype, void * inouttype, int * len, MPI_Datatype *)
  {
    libMesh::FEType *in = static_cast<libMesh::FEType *>(intype);
    libMesh::FEType *inout = static_cast<libMesh::FEType *>(inouttype);
    for (int i=0; i != *len; ++i)
      {
        inout[i].family = std::max(in[i].family, inout[i].family);
        inout[i].order = std::max(in[i].order, inout[i].order);
      }
  }

  static void fetype_min (void * intype, void * inouttype, int * len, MPI_Datatype *)
  {
    libMesh::FEType *in = static_cast<libMesh::FEType *>(intype);
    libMesh::FEType *inout = static_cast<libMesh::FEType *>(inouttype);
    for (int i=0; i != *len; ++i)
      {
        inout[i].family = std::min(in[i].family, inout[i].family);
        inout[i].order = std::min(in[i].order, inout[i].order);
      }
  }

  static void fetype_sum (void *, void *, int *, MPI_Datatype *)
  {
    libmesh_not_implemented();
  }

  static MPI_Op max()
  {
    // _static_op never gets freed, but it only gets committed once,
    // so it's not a *huge* memory leak...
    static MPI_Op _static_op;
    static bool _is_initialized = false;
    if (!_is_initialized)
      {
        timpi_call_mpi
          (MPI_Op_create (fetype_max, /*commute=*/ true,
                          &_static_op));

        _is_initialized = true;
      }

    return _static_op;
  }
  static MPI_Op min()
  {
    // _static_op never gets freed, but it only gets committed once,
    // so it's not a *huge* memory leak...
    static MPI_Op _static_op;
    static bool _is_initialized = false;
    if (!_is_initialized)
      {
        timpi_call_mpi
          (MPI_Op_create (fetype_min, /*commute=*/ true,
                          &_static_op));

        _is_initialized = true;
      }

    return _static_op;
  }
  static MPI_Op sum()
  {
    libmesh_not_implemented();
  }

#endif // LIBMESH_HAVE_MPI
};


template <>
struct Attributes<libMesh::FEType>
{
  static const bool has_min_max = true;
  static void set_lowest(libMesh::FEType & x) {
    x.family = libMesh::FEFamily(0);
    x.order = libMesh::Order(0);
  }
  static void set_highest(libMesh::FEType & x) {
    x.family = libMesh::INVALID_FE;
    x.order = libMesh::INVALID_ORDER;
  }
};


// StandardType<> specializations to return a derived MPI datatype
// to handle communication of FETypes.
template <>
class StandardType<libMesh::FEType> : public DataType
{
public:
  explicit
  StandardType(const libMesh::FEType * example=nullptr)
  {
#ifdef LIBMESH_HAVE_MPI
    using libMesh::FEType;

    // We need an example for MPI_Address to use
    FEType * ex;
    std::unique_ptr<FEType> temp;
    if (example)
      ex = const_cast<FEType *>(example);
    else
      {
        temp = std::make_unique<FEType>();
        ex = temp.get();
      }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    constexpr std::size_t structsize = 5;
#else
    constexpr std::size_t structsize = 2;
#endif

    // Our FEType enums pack into int
    MPI_Datatype int_types[structsize];

    // We require MPI-2 here:
    int blocklengths[structsize];
    MPI_Aint start, displs[structsize];
    MPI_Datatype tmptype;

    timpi_call_mpi
      (MPI_Get_address (ex, &start));
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    timpi_call_mpi
      (MPI_Get_address (&(ex->order), &displs[0]));
    timpi_call_mpi
      (MPI_Get_address (&(ex->radial_order), &displs[1]));
    timpi_call_mpi
      (MPI_Get_address (&(ex->family), &displs[2]));
    timpi_call_mpi
      (MPI_Get_address (&(ex->radial_family), &displs[3]));
    timpi_call_mpi
      (MPI_Get_address (&(ex->inf_map), &displs[4]));
#else // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    timpi_call_mpi
      (MPI_Get_address (&(ex->order), &displs[0]));
    timpi_call_mpi
      (MPI_Get_address (&(ex->family), &displs[1]));
#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

    // subtract off beginning of the structure to get offsets.  We'll
    // be creating a struct type later to account for any weird
    // padding issues.
    for (std::size_t i = 0; i != structsize; ++i)
    {
      displs[i] -= start;
      int_types[i] = MPI_INT;
      blocklengths[i] = 1;
    }

    // create a prototype structure
    timpi_call_mpi
      (MPI_Type_create_struct (structsize, blocklengths, displs,
                               int_types, &tmptype));
    timpi_call_mpi
      (MPI_Type_commit (&tmptype));

    // resize the structure type to account for padding, if any
    timpi_call_mpi
      (MPI_Type_create_resized (tmptype, 0, sizeof(FEType),
                                &_datatype));

    timpi_call_mpi
      (MPI_Type_commit (&_datatype));

    timpi_call_mpi
      (MPI_Type_free (&tmptype));
#else
    libmesh_ignore(example);
#endif // #ifdef LIBMESH_HAVE_MPI
  }

  StandardType(const StandardType<libMesh::FEType> & timpi_mpi_var(t))
    : DataType()
  {
    timpi_call_mpi (MPI_Type_dup (t._datatype, &_datatype));
  }

  ~StandardType() { this->free(); }

  static const bool is_fixed_type = true;
};

} // namespace TIMPI

#endif // LIBMESH_PARALLEL_FE_TYPE_H
