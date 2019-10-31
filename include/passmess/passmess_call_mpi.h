// The PassMess Message-Passing Parallelism Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef PASSMESS_CALL_MPI_H
#define PASSMESS_CALL_MPI_H

// libMesh Includes
#include "libmesh/libmesh_config.h"

// PassMess Includes
#include "libmesh/passmess_assert.h"

#ifdef LIBMESH_HAVE_MPI
#  include "mpi.h"

/**
 * Macro to use arguments only when MPI is in use
 */
#define passmess_mpi_var(var) var

/**
 * Macros to test MPI return values
 */
#ifndef NDEBUG
#define passmess_assert_mpi_success(error_code)                         \
  do                                                                    \
    {                                                                   \
      if (error_code != MPI_SUCCESS)                                    \
        {                                                               \
          char passmess_mpi_error_string[MPI_MAX_ERROR_STRING+1];       \
          int passmess_mpi_error_string_len;                            \
          MPI_Error_string(error_code, passmess_mpi_error_string,       \
                           &passmess_mpi_error_string_len);             \
          passmess_assert_equal_to_msg(error_code, MPI_SUCCESS,         \
                                       passmess_mpi_error_string);      \
        }                                                               \
    }                                                                   \
  while (0)

#else

#define passmess_assert_mpi_success(error_code)  ((void) 0)

#endif



// Only catch MPI return values when asserts are active.
#ifndef NDEBUG
#define passmess_call_mpi(mpi_call)                             \
  do                                                            \
    {                                                           \
      unsigned int passmess_mpi_error_code = mpi_call;          \
      passmess_assert_mpi_success (passmess_mpi_error_code);    \
    }                                                           \
  while (0)

#else

#define passmess_call_mpi(mpi_call)             \
  do                                            \
    {                                           \
      mpi_call;                                 \
    }                                           \
  while (0)
#endif

#else // LIBMESH_HAVE_MPI

#define passmess_mpi_var(var)

#define passmess_call_mpi(mpi_call)             \
  do {}                                         \
  while (0)

#endif // LIBMESH_HAVE_MPI


#endif // PASSMESS_CALL_MPI_H
