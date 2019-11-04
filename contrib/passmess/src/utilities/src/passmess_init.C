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


// Local includes
#include "passmess/passmess_init.h"

// PassMess includes
#include "passmess/communicator.h"
#include "passmess/passmess_assert.h"



#ifdef PASSMESS_HAVE_MPI
void PassMess_MPI_Handler (MPI_Comm *, int *, ...)
{
  passmess_not_implemented();
}
#endif


namespace PassMess
{

#ifdef PASSMESS_HAVE_MPI
PassMessInit::PassMessInit (int argc, const char * const * argv,
                            bool using_threads,
                            bool handle_mpi_errors,
                            MPI_Comm COMM_WORLD_IN) :
  i_initialized_mpi(false),
  err_handler_set(false)
{
  // Check whether the calling program has already initialized
  // MPI, and avoid duplicate Init/Finalize
  int flag;
  passmess_call_mpi(MPI_Initialized (&flag));

  if (!flag)
    {
      int mpi_thread_provided;
      const int mpi_thread_requested = using_threads ?
        MPI_THREAD_FUNNELED :
        MPI_THREAD_SINGLE;

      passmess_call_mpi
        (MPI_Init_thread (&argc, const_cast<char ***>(&argv),
                          mpi_thread_requested, &mpi_thread_provided));

      if (using_threads &&
          (mpi_thread_provided < MPI_THREAD_FUNNELED))
        {
          // Ideally, if an MPI stack tells us it's unsafe for us
          // to use threads, we should scream and die or at least
          // disable threads.
          //
          // In practice, we've encountered one MPI stack (an mvapich2
          // configuration) that returned MPI_THREAD_SINGLE as a
          // proper warning, two stacks that handle
          // MPI_THREAD_FUNNELED properly, and two current stacks plus
          // a couple old stacks that return MPI_THREAD_SINGLE but
          // support threaded runs anyway, so we just emit a warning.
          //
          passmess_warning("Warning: MPI failed to guarantee MPI_THREAD_FUNNELED\n"
                           << "for a threaded run.\n"
                           << "Be sure your library is funneled-thread-safe..."
                           << std::endl);
        }
      this->i_initialized_mpi = true;
    }

  // Duplicate the input communicator for internal use
  // And get a Communicator copy too, to use
  // as a default for that API
  this->_comm = new Communicator(COMM_WORLD_IN);

  // Set up an MPI error handler if requested.  This helps us get
  // into a debugger with a proper stack when an MPI error occurs.
  if (handle_mpi_errors)
    {
      passmess_call_mpi
        (MPI_Comm_create_errhandler(PassMess_MPI_Handler, &my_errhandler));
      passmess_call_mpi
        (MPI_Comm_set_errhandler(COMM_WORLD_IN, my_errhandler));
      passmess_call_mpi
        (MPI_Comm_set_errhandler(MPI_COMM_WORLD, my_errhandler));
      err_handler_set = true;
    }
}
#else
PassMessInit::PassMessInit (int /* argc */, const char * const * /* argv */,
                            bool /* handle_mpi_errors */,
                            bool /* using_threads */) :
  i_initialized_mpi(false)
{
  this->_comm = new Communicator(); // So comm() doesn't dereference null
}
#endif



PassMessInit::~PassMessInit()
{
  // Every processor had better be ready to exit at the same time.
  // This would be a passmess_parallel_only() function, except that
  // passmess_parallel_only() uses passmess_assert() which throws an
  // exception which causes compilers to scream about exceptions
  // inside destructors.

  // Even if we're not doing parallel_only debugging, we don't want
  // one processor to try to exit until all others are done working.
  this->comm().barrier();

#ifdef PASSMESS_HAVE_MPI
  if (err_handler_set)
    {
      unsigned int error_code =
        MPI_Errhandler_free(&my_errhandler);
      if (error_code != MPI_SUCCESS)
        {
          std::cerr <<
            "Failure when freeing MPI_Errhandler! Continuing..." <<
            std::endl;
        }
    }

  this->comm().clear();
  delete this->_comm;

  if (this->i_initialized_mpi)
    {
      // We can't just passmess_assert here because destructor,
      // but we ought to report any errors
      unsigned int error_code = MPI_Finalize();
      if (error_code != MPI_SUCCESS)
        {
          char error_string[MPI_MAX_ERROR_STRING+1];
          int error_string_len;
          MPI_Error_string(error_code, error_string,
                           &error_string_len);
          std::cerr << "Failure from MPI_Finalize():\n"
                    << error_string << std::endl;
        }
    }
#else
  delete this->_comm;
#endif
}

} // namespace PassMess
