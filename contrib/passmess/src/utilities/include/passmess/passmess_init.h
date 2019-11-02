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



#ifndef PASSMESS_INIT_H
#define PASSMESS_INIT_H


// Local includes
#include "passmess/passmess_config.h"

// C/C++ includes

#ifdef PASSMESS_HAVE_MPI
# include "passmess/ignore_warnings.h"
# include <mpi.h>
# include "passmess/restore_warnings.h"
#endif // #ifdef PASSMESS_HAVE_MPI

namespace PassMess
{

// Forward declarations
  class Communicator;

/**
 * The \p PassMessInit class, when constructed, initializes
 * any dependent libraries (e.g. MPI).
 *
 * For many users, a single PassMessInit object should be created at
 * the start of your main() function.
 *
 * Since "it is best not to perform much more than a return rc after
 * calling MPI_Finalize", applications which want to do anything after
 * PassMessInit destruction should manage MPI initialization and
 * finalization manually.
 */
class PassMessInit
{
public:
#ifdef PASSMESS_HAVE_MPI
  /**
   * Initialize the library for use, with the command line options
   * provided.  This will e.g. call MPI_Init_thread if MPI is
   * available and enabled and has not already been initialized.
   *
   * using_threads should be set to true when MPI_THREAD_FUNNELED
   * initialization is desired.
   *
   * handle_mpi_errors should be set to true to use the
   * "passmess_not_implemented()" behavior (e.g. throwing an
   * exception) as an MPI error handler, which can aid in debugging.
   *
   * When building with MPI, this method may take an optional
   * parameter to use a user-specified MPI communicator.
   */
  PassMessInit(int argc, const char * const * argv,
               bool using_threads=false,
               bool handle_mpi_errors=false,
               MPI_Comm COMM_WORLD_IN=MPI_COMM_WORLD);
#else
  PassMessInit(int argc, const char * const * argv,
               bool using_threads=false,
               bool handle_mpi_errors=false);
#endif

  /**
   * Destructor.
   * Finalizes MPI if that was initialized by PassMessInit.
   */
  virtual ~PassMessInit();

  /**
   * Returns the Communicator created by this object, which will be a
   * compatibility shim if MPI is not enabled, or a wrapper for the
   * user-input MPI_Comm if we were constructed with one, or a wrapper
   * for MPI_COMM_WORLD by default.
   */
  const Communicator & comm() const { return *_comm; }

  Communicator & comm() { return *_comm; }

private:
  Communicator * _comm;

  bool i_initialized_mpi;

#ifdef PASSMESS_HAVE_MPI
  MPI_Errhandler my_errhandler;
  bool err_handler_set;
#endif
};

} // namespace PassMess

#endif // PASSMESS_INIT_H
