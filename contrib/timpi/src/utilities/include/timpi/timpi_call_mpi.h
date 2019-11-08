// The TIMPI Message-Passing Parallelism Library.
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


#ifndef TIMPI_CALL_MPI_H
#define TIMPI_CALL_MPI_H

// TIMPI includes
#include "timpi/timpi_assert.h"
#include "timpi/timpi_config.h"

// C/C++ includes
#ifdef TIMPI_HAVE_MPI
#  include "timpi/ignore_warnings.h"
#  include "mpi.h"
#  include "timpi/restore_warnings.h"

/**
 * Macro to use arguments only when MPI is in use
 */
#define timpi_mpi_var(var) var

/**
 * Macros to test MPI return values
 */
#ifndef NDEBUG
#define timpi_assert_mpi_success(error_code)                         \
  do                                                                    \
    {                                                                   \
      if (error_code != MPI_SUCCESS)                                    \
        {                                                               \
          char timpi_mpi_error_string[MPI_MAX_ERROR_STRING+1];       \
          int timpi_mpi_error_string_len;                            \
          MPI_Error_string(error_code, timpi_mpi_error_string,       \
                           &timpi_mpi_error_string_len);             \
          timpi_assert_equal_to_msg(error_code, MPI_SUCCESS,         \
                                       timpi_mpi_error_string);      \
        }                                                               \
    }                                                                   \
  while (0)

#else

#define timpi_assert_mpi_success(error_code)  ((void) 0)

#endif



// Only catch MPI return values when asserts are active.
#ifndef NDEBUG
#define timpi_call_mpi(mpi_call)                             \
  do                                                            \
    {                                                           \
      unsigned int timpi_mpi_error_code = mpi_call;          \
      timpi_assert_mpi_success (timpi_mpi_error_code);    \
    }                                                           \
  while (0)

#else

#define timpi_call_mpi(mpi_call)             \
  do                                            \
    {                                           \
      mpi_call;                                 \
    }                                           \
  while (0)
#endif

#else // TIMPI_HAVE_MPI

#define timpi_mpi_var(var)

#define timpi_call_mpi(mpi_call)             \
  do {}                                         \
  while (0)

#endif // TIMPI_HAVE_MPI


#endif // TIMPI_CALL_MPI_H
