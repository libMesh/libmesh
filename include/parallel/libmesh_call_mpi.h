// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_CALL_MPI_H
#define LIBMESH_CALL_MPI_H

// libMesh Includes
#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_MPI

/**
 * Macros to test MPI return values
 */
#ifndef NDEBUG
#define libmesh_assert_mpi_success(error_code)                          \
  do                                                                    \
    {                                                                   \
      if (error_code != MPI_SUCCESS)                                    \
        {                                                               \
          char libmesh_mpi_error_string[MPI_MAX_ERROR_STRING+1];        \
          int libmesh_mpi_error_string_len;                             \
          MPI_Error_string(error_code, libmesh_mpi_error_string,        \
                           &libmesh_mpi_error_string_len);              \
          libmesh_assert_equal_to_msg(error_code, MPI_SUCCESS,          \
                                      libmesh_mpi_error_string);        \
        }                                                               \
    }                                                                   \
  while (0)

#else

#define libmesh_assert_mpi_success(error_code)  ((void) 0)

#endif



// Only catch MPI return values when asserts are active.
#ifndef NDEBUG
#define libmesh_call_mpi(mpi_call)                              \
  do                                                            \
    {                                                           \
      unsigned int libmesh_mpi_error_code = mpi_call;           \
      libmesh_assert_mpi_success (libmesh_mpi_error_code);      \
    }                                                           \
  while (0)

#else

#define libmesh_call_mpi(mpi_call)              \
  do                                            \
    {                                           \
      mpi_call;                                 \
    }                                           \
  while (0)
#endif

#endif // LIBMESH_HAVE_MPI


#endif // LIBMESH_CALL_MPI_H
