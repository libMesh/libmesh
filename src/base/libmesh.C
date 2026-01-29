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


// Local includes
#include "libmesh/libmesh.h"

// libMesh includes
#include "libmesh/getpot.h"
#include "libmesh/reference_counter.h"
#include "libmesh/libmesh_singleton.h"
#include "libmesh/remote_elem.h"
#include "libmesh/threads.h"
#include "libmesh/parallel_only.h"
#include "libmesh/print_trace.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/perf_log.h"
#include "libmesh/thread_buffered_syncbuf.h"

// TIMPI includes
#include "timpi/communicator.h"
#include "timpi/timpi_init.h"

// C/C++ includes
#include <cstdlib>
#include <iostream>
#include <fstream>

#ifdef LIBMESH_ENABLE_EXCEPTIONS
#include <exception>
#include <optional>
#endif

#ifdef LIBMESH_HAVE_OPENMP
#include <omp.h>
#endif

#include "stdlib.h" // C, not C++ - we need setenv() from POSIX


#if defined(LIBMESH_HAVE_MPI)
# include "libmesh/ignore_warnings.h"
# include <mpi.h>
# include "libmesh/restore_warnings.h"
#endif // #if defined(LIBMESH_HAVE_MPI)

#if defined(LIBMESH_HAVE_PETSC)
# include "libmesh/petsc_solver_exception.h"
# include <petsc.h>
# include <petscerror.h>
# include "libmesh/petscdmlibmesh.h"
# if defined(LIBMESH_HAVE_SLEPC)
// Ignore unused variable warnings from SLEPc
#  include "libmesh/ignore_warnings.h"
#  include "libmesh/slepc_macro.h"
#  include <slepc.h>
#  include "libmesh/restore_warnings.h"
# endif // #if defined(LIBMESH_HAVE_SLEPC)
#endif // #if defined(LIBMESH_HAVE_PETSC)

#ifdef LIBMESH_HAVE_NETGEN
// We need the nglib namespace, because it's used everywhere in nglib
// and we don't get binary compatibility without it.
//
// We need the nglib namespace *here*, because somehow nobody ever
// figured out to just put it in nglib.h?
namespace nglib {
#include "netgen/nglib/nglib.h"
}
#endif

// If we're using MPI and VTK has been detected, we need to do some
// MPI initialize/finalize stuff for VTK.
#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
#include "libmesh/ignore_warnings.h"
# include "vtkMPIController.h"
#include "libmesh/restore_warnings.h"
#endif

#include <mutex>

// --------------------------------------------------------
// Local anonymous namespace to hold miscellaneous bits
namespace {

std::unique_ptr<GetPot> command_line;

std::set<std::string> command_line_name_set;

std::unique_ptr<std::ofstream> _ofstream;
// If std::cout and std::cerr are redirected, we need to
// be a little careful and save the original streambuf objects,
// replacing them in the destructor before program termination.
std::streambuf * out_buf (nullptr);
std::streambuf * err_buf (nullptr);

std::unique_ptr<libMesh::Threads::task_scheduler_init> task_scheduler;
#if defined(LIBMESH_HAVE_PETSC)
bool libmesh_initialized_petsc = false;
#endif
#if defined(LIBMESH_HAVE_SLEPC)
bool libmesh_initialized_slepc = false;
#endif

} // anonymous namespace



#ifdef LIBMESH_HAVE_MPI
void libMesh_MPI_Handler (MPI_Comm *, int *, ...)
{
  libmesh_not_implemented();
}
#endif


namespace libMesh
{

/**
 * Namespaces don't provide private data,
 * so let's take the data we would like
 * private and put it in an obnoxious
 * namespace.  At least that way it is a
 * pain to use, thus discouraging errors.
 */
namespace libMeshPrivateData {

/**
 * Flag that tells if \p init() has been called.
 */
extern bool _is_initialized;

/**
 * The default solver package to use.
 */
extern SolverPackage _solver_package;
}


// ------------------------------------------------------------
// libMesh data initialization
#ifdef LIBMESH_HAVE_MPI
MPI_Comm           GLOBAL_COMM_WORLD = MPI_COMM_NULL;
#else
int                GLOBAL_COMM_WORLD = 0;
#endif

#ifdef LIBMESH_ENABLE_EXCEPTIONS
std::terminate_handler LibMeshInit::_old_terminate_handler;
#endif

OStreamProxy out(std::cout);
OStreamProxy err(std::cerr);

// Our thread-safe wrappers and the “pre-wrap” sink pointers
std::unique_ptr<ThreadBufferedSyncbuf> _out_syncd_thread_buffer;
std::unique_ptr<ThreadBufferedSyncbuf> _err_syncd_thread_buffer;
std::streambuf * _out_prewrap_buf = nullptr;
std::streambuf * _err_prewrap_buf = nullptr;

// Install after redirection/drop decisions
void install_thread_buffered_sync()
{
  auto check_stored_buffer = [](const auto * const libmesh_dbg_var(stored_buffer)) {
    libmesh_assert_msg(!stored_buffer,
                       "Oops, we've already stored a prewrapped buffer. We must be calling this "
                       "function a second time in which case we're going to lose the already "
                       "stored prewrapped buffer forever");
  };

  // libMesh::out
  if (auto * ob = libMesh::out.rdbuf(); ob)
    {
      check_stored_buffer(_out_prewrap_buf);
      _out_prewrap_buf = ob;
      _out_syncd_thread_buffer =
          std::make_unique<ThreadBufferedSyncbuf>(*ob, /*flush_on_newline=*/true);
      libMesh::out.rdbuf(_out_syncd_thread_buffer.get());
    }

  // libMesh::err
  if (auto * eb = libMesh::err.rdbuf(); eb)
    {
      check_stored_buffer(_err_prewrap_buf);
      _err_prewrap_buf = eb;
      _err_syncd_thread_buffer =
          std::make_unique<ThreadBufferedSyncbuf>(*eb, /*flush_on_newline=*/true);
      libMesh::err.rdbuf(_err_syncd_thread_buffer.get());
    }
}

// Uninstall before restoring/redirection teardown
void uninstall_thread_buffered_sync()
{
  if (_out_syncd_thread_buffer)
    {
      // flush any thread-local leftovers on this thread
      libMesh::out << std::flush;
      libMesh::out.rdbuf(_out_prewrap_buf);
      _out_syncd_thread_buffer.reset();
      _out_prewrap_buf = nullptr;
    }
  if (_err_syncd_thread_buffer)
    {
      libMesh::err << std::flush;
      libMesh::err.rdbuf(_err_prewrap_buf);
      _err_syncd_thread_buffer.reset();
      _err_prewrap_buf = nullptr;
    }
}


/**
 * Helper to do cleanup from both destructor and terminate
 */
void cleanup_stream_buffers()
{
  // Before resetting the stream buffers, let's remove our thread wrappers
  if (!libMesh::on_command_line ("--disable-thread-safe-output"))
    uninstall_thread_buffered_sync();

  if (libMesh::on_command_line ("--redirect-stdout") ||
      libMesh::on_command_line ("--redirect-output"))
    {
      // If stdout/stderr were redirected to files, reset them now.
      libMesh::out.rdbuf (out_buf);
      libMesh::err.rdbuf (err_buf);
    }

  // If we built our own output streams, we want to clean them up.
  if (libMesh::on_command_line ("--separate-libmeshout"))
    {
      delete libMesh::out.get();
      delete libMesh::err.get();

      libMesh::out.reset(std::cout);
      libMesh::err.reset(std::cerr);
    }
}


bool warned_about_auto_ptr(false);

PerfLog            perflog ("libMesh",
#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
                            true
#else
                            false
#endif
                            );


#ifdef LIBMESH_USE_COMPLEX_NUMBERS
const Number       imaginary (0., 1.);
// const Number       zero      (0., 0.);
#else
// const Number       zero = 0.;
#endif

// This is now a static constant in the header; no reason not to let
// the compiler inline it.

// const unsigned int invalid_uint = static_cast<unsigned int>(-1);



// ------------------------------------------------------------
// libMesh::libMeshPrivateData data initialization
#ifdef LIBMESH_HAVE_MPI
MPI_Errhandler libmesh_errhandler;

processor_id_type libMesh::libMeshPrivateData::_n_processors = 1;
processor_id_type libMesh::libMeshPrivateData::_processor_id = 0;
#endif
int           libMesh::libMeshPrivateData::_n_threads = 1; /* Threads::task_scheduler_init::automatic; */
bool          libMesh::libMeshPrivateData::_is_initialized = false;
SolverPackage libMesh::libMeshPrivateData::_solver_package =
#if   defined(LIBMESH_HAVE_PETSC)    // PETSc is the default
  PETSC_SOLVERS;
#elif defined(LIBMESH_TRILINOS_HAVE_AZTECOO) // Use Trilinos if PETSc isn't there
TRILINOS_SOLVERS;
#elif defined(LIBMESH_HAVE_EIGEN)    // Use Eigen if neither are there
EIGEN_SOLVERS;
#elif defined(LIBMESH_HAVE_LASPACK)  // Use LASPACK as a last resort
LASPACK_SOLVERS;
#else                        // No valid linear solver package at compile time
INVALID_SOLVER_PACKAGE;
#endif



// ------------------------------------------------------------
// libMesh functions

bool initialized()
{
  return libMeshPrivateData::_is_initialized;
}



bool closed()
{
  return !libMeshPrivateData::_is_initialized;
}



void libmesh_abort()
{
  libMesh::perflog.clear();

  // Now that we're done with output we should clean up our stream
  // buffers; if we fail to uninstall_thread_buffered_sync() when
  // needed we can end up seeing a segfault in iostreams destructors
  cleanup_stream_buffers();

  // If we have MPI and it has been initialized, we need to be sure
  // and call MPI_Abort instead of std::abort, so that the parallel
  // job can die nicely.
#if defined(LIBMESH_HAVE_MPI)
  int mpi_initialized;
  MPI_Initialized (&mpi_initialized);

  if (mpi_initialized)
    MPI_Abort(libMesh::GLOBAL_COMM_WORLD, 1);
#endif

#ifdef LIBMESH_ENABLE_EXCEPTIONS
  // The system terminate_handler may do useful things, or the user
  // may have set their own terminate handler that we want to call.
  LibMeshInit::_old_terminate_handler();
#endif

  // The last attempt to die if nothing else has killed us
  std::abort();
}



LibMeshInit::LibMeshInit (int argc, const char * const * argv,
                          TIMPI::communicator COMM_WORLD_IN, int n_threads)
{
  // should _not_ be initialized already.
  libmesh_assert (!libMesh::initialized());

  // Build a command-line parser.
  command_line = std::make_unique<GetPot>(argc, argv);

  // Disable performance logging upon request
  {
    if (libMesh::on_command_line ("--disable-perflog"))
      libMesh::perflog.disable_logging();
  }

  auto check_empty_command_line_value = [](auto & cl, const std::string & option) {
    libmesh_error_msg_if(cl.search(option),
                         "Detected option " << option << " with no value.  Did you forget '='?");
  };

  // Build a task scheduler
  {
    // Get the requested number of threads, defaults to 1 to avoid MPI and
    // multithreading competition.  If you would like to use MPI and multithreading
    // at the same time then (n_mpi_processes_per_node)x(n_threads) should be the
    //  number of processing cores per node.
    std::vector<std::string> n_threads_opt(2);
    n_threads_opt[0] = "--n_threads";
    n_threads_opt[1] = "--n-threads";
    libMesh::libMeshPrivateData::_n_threads =
      libMesh::command_line_value(n_threads_opt, n_threads);

    if (libMesh::libMeshPrivateData::_n_threads == -1)
      {
        for (auto & option : n_threads_opt)
          check_empty_command_line_value(*command_line, option);

        libMesh::libMeshPrivateData::_n_threads = 1;
      }

    // If there's no threading model active, force _n_threads==1
#if !LIBMESH_USING_THREADS
    if (libMesh::libMeshPrivateData::_n_threads != 1)
      {
        libMesh::libMeshPrivateData::_n_threads = 1;
        libmesh_warning("Warning: You requested --n-threads>1 but no threading model is active!\n"
                        << "Forcing --n-threads==1 instead!");
      }
#endif

    // Set the number of OpenMP threads to the same as the number of threads libMesh is going to use
#ifdef LIBMESH_HAVE_OPENMP
    omp_set_num_threads(libMesh::libMeshPrivateData::_n_threads);
#endif

    task_scheduler = std::make_unique<Threads::task_scheduler_init>(libMesh::n_threads());
  }

  // Construct singletons who may be at risk of the
  // "static initialization order fiasco"
  Singleton::setup();

  // Make sure the construction worked
  libmesh_assert(remote_elem);

  const bool using_threads = libMesh::n_threads() > 1;
  const bool handle_mpi_errors = false; // libMesh does this

#if defined(LIBMESH_HAVE_MPI)

  // Allow the user to bypass MPI initialization
  if (!libMesh::on_command_line ("--disable-mpi"))
    {
      int mpi_thread_request = using_threads;
      const auto mpi_thread_type = libMesh::command_line_value("--mpi-thread-type", std::string(""));
      if (mpi_thread_type.empty())
        {
          check_empty_command_line_value(*command_line, "--mpi-thread-type");
#if defined(PETSC_HAVE_STRUMPACK) && defined(PETSC_HAVE_SLATE)
          // For GPU computations, the solver strumpack uses slate which always requests
          // MPI_THREAD_MULTIPLE. The solution here is not perfect because the run may never be
          // using strumpack, but we believe it's better to force the MPI library to use locks
          // whenever it accesses the message queue, that is, when processing any sends and receive,
          // than it is to require users to pre-announce/signal what solvers they are using through
          // --mpi-thread-type
          mpi_thread_request = 3;
#endif
        }
      else
        {
          if (mpi_thread_type == "single")
            {
              if (using_threads)
                libmesh_error_msg("We are using threads, so we require more mpi thread support "
                                  "than '--mpi-thread-type=single'");
              mpi_thread_request = 0;
            }
          else if (mpi_thread_type == "funneled")
            mpi_thread_request = 1;
          else if (mpi_thread_type == "serialized")
            mpi_thread_request = 2;
          else if (mpi_thread_type == "multiple")
            mpi_thread_request = 3;
          else
            libmesh_error_msg(
                "Unsupported mpi thread type '"
                << mpi_thread_type
                << "'. Allowed options are 'single', 'funneled', 'serialized', and 'multiple'");
        }

      this->_timpi_init =
        new TIMPI::TIMPIInit(argc, argv, mpi_thread_request,
                             handle_mpi_errors, COMM_WORLD_IN);
      _comm = new Parallel::Communicator(this->_timpi_init->comm().get());

      const std::string timpi_sync =
        libMesh::command_line_value("--timpi-sync", std::string("nbx"));
      _comm->sync_type(timpi_sync);

      libMesh::GLOBAL_COMM_WORLD = COMM_WORLD_IN;

      //MPI_Comm_set_name not supported in at least SGI MPT's MPI implementation
      //MPI_Comm_set_name (libMesh::COMM_WORLD, "libMesh::COMM_WORLD");

      libMeshPrivateData::_processor_id =
        cast_int<processor_id_type>(this->comm().rank());
      libMeshPrivateData::_n_processors =
        cast_int<processor_id_type>(this->comm().size());

      // Set up an MPI error handler if requested.  This helps us get
      // into a debugger with a proper stack when an MPI error occurs.
      if (libMesh::on_command_line ("--handle-mpi-errors"))
        {
          timpi_call_mpi
            (MPI_Comm_create_errhandler(libMesh_MPI_Handler, &libmesh_errhandler));
          timpi_call_mpi
            (MPI_Comm_set_errhandler(libMesh::GLOBAL_COMM_WORLD, libmesh_errhandler));
          timpi_call_mpi
            (MPI_Comm_set_errhandler(MPI_COMM_WORLD, libmesh_errhandler));
        }
    }

  // Could we have gotten bad values from the above calls?
  libmesh_assert_greater (libMeshPrivateData::_n_processors, 0);

  // The cast_int already tested _processor_id>=0
  // libmesh_assert_greater_equal (libMeshPrivateData::_processor_id, 0);

  // Let's be sure we properly initialize on every processor at once:
  libmesh_parallel_only(this->comm());

#else
  libmesh_ignore(COMM_WORLD_IN);
  this->_timpi_init =
    new TIMPI::TIMPIInit(argc, argv, using_threads,
                         handle_mpi_errors);
  _comm = new Parallel::Communicator(this->_timpi_init->comm().get());
#endif

#if defined(LIBMESH_HAVE_PETSC)

  // Allow the user to bypass PETSc initialization
  if (!libMesh::on_command_line ("--disable-petsc")

#if defined(LIBMESH_HAVE_MPI)
      // If the user bypassed MPI, we'd better be safe and assume that
      // PETSc was built to require it; otherwise PETSc initialization
      // dies.
      && !libMesh::on_command_line ("--disable-mpi")
#endif
      )
    {
#ifdef LIBMESH_HAVE_MPI
      PETSC_COMM_WORLD = libMesh::GLOBAL_COMM_WORLD;
#else
      // PETSc --with-mpi=0 doesn't like our default "communicator" 0
      this->_comm->get() = PETSC_COMM_SELF;
#endif

      // Check whether the calling program has already initialized
      // PETSc, and avoid duplicate Initialize/Finalize
      PetscBool petsc_already_initialized;
      LibmeshPetscCallA(libMesh::GLOBAL_COMM_WORLD, PetscInitialized(&petsc_already_initialized));
      if (petsc_already_initialized != PETSC_TRUE)
        libmesh_initialized_petsc = true;
# if defined(LIBMESH_HAVE_SLEPC)

      // If SLEPc allows us to check whether the calling program
      // has already initialized it, we do that, and avoid
      // duplicate Initialize/Finalize.
      // We assume that SLEPc will handle PETSc appropriately,
      // which it does in the versions we've checked.
      if (!SlepcInitializeCalled)
        {
          LibmeshPetscCallA(libMesh::GLOBAL_COMM_WORLD, SlepcInitialize  (&argc, const_cast<char ***>(&argv), nullptr, nullptr));
          libmesh_initialized_slepc = true;
        }
# else
      if (libmesh_initialized_petsc)
        LibmeshPetscCallA(libMesh::GLOBAL_COMM_WORLD, PetscInitialize (&argc, const_cast<char ***>(&argv), nullptr, nullptr));
# endif
      // Register the reference implementation of DMlibMesh
      LibmeshPetscCallA(libMesh::GLOBAL_COMM_WORLD, DMRegister(DMLIBMESH, DMCreate_libMesh));
    }
#endif

#ifdef LIBMESH_HAVE_NETGEN
  nglib::Ng_Init();
#endif

#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
  // Do MPI initialization for VTK.
  _vtk_mpi_controller = vtkMPIController::New();
  _vtk_mpi_controller->Initialize(&argc, const_cast<char ***>(&argv), /*initialized_externally=*/1);
  _vtk_mpi_controller->SetGlobalController(_vtk_mpi_controller);
#endif

  // Re-parse the command-line arguments.  Note that PETSc and MPI
  // initialization above may have removed command line arguments
  // that are not relevant to this application in the above calls.
  // We don't want a false-positive by detecting those arguments.
  //
  // Note: this seems overly paranoid/like it should be unnecessary,
  // plus we were doing it wrong for many years and not clearing the
  // existing GetPot object before re-parsing the command line, so all
  // the command line arguments appeared twice in the GetPot object...
  command_line = std::make_unique<GetPot>(argc, argv);

  // Not syncing with stdio is an optimization when simultaneous
  // C and C++ style access to output streams is not required.
  // The amount of benefit which occurs is probably implementation
  // defined, and may be nothing.  On the other hand, I have seen
  // some IO tests where IO performance improves by a factor of two.
  // However, not syncing the streams means thread safety is no longer
  // guaranteed and we have indeed seen thread sanitizer warnings
  // during concurrent writes to std::cout when toggling this boolean
  // to true. We default to preferring safety and correctness instead
  // of performance
  if (libMesh::on_command_line ("--dont-sync-with-stdio"))
    std::ios::sync_with_stdio(false);

  // Honor the --separate-libmeshout command-line option.
  // When this is specified, the library uses an independent ostream
  // for libMesh::out/libMesh::err messages, and
  // std::cout and std::cerr are untouched by any other options
  if (libMesh::on_command_line ("--separate-libmeshout"))
    {
      // Redirect.  We'll share streambufs with cout/cerr for now, but
      // presumably anyone using this option will want to replace the
      // bufs later.
      std::ostream * newout = new std::ostream(std::cout.rdbuf());
      libMesh::out = *newout;
      std::ostream * newerr = new std::ostream(std::cerr.rdbuf());
      libMesh::err = *newerr;
    }

  // Process command line arguments for redirecting stdout/stderr.
  bool
    cmdline_has_redirect_stdout = libMesh::on_command_line ("--redirect-stdout"),
    cmdline_has_redirect_output = libMesh::on_command_line ("--redirect-output");

  // The --redirect-stdout command-line option has been deprecated in
  // favor of "--redirect-output basename".
  if (cmdline_has_redirect_stdout)
    libmesh_warning("The --redirect-stdout command line option has been deprecated. "
                    "Use '--redirect-output basename' instead.");

  // Honor the "--redirect-stdout" and "--redirect-output basename"
  // command-line options.  When one of these is specified, each
  // processor sends libMesh::out/libMesh::err messages to
  // stdout.processor.#### (default) or basename.processor.####.
  if (cmdline_has_redirect_stdout || cmdline_has_redirect_output)
    {
      std::string basename = "stdout";

      // Look for following argument if using new API
      if (cmdline_has_redirect_output)
        {
          // Set the cursor to the correct location in the list of command line arguments.
          command_line->search(1, "--redirect-output");

          // Get the next option on the command line as a string.
          std::string next_string = "";
          next_string = command_line->next(next_string);

          // If the next string starts with a dash, we assume it's
          // another flag and not a file basename requested by the
          // user.
          if (next_string.size() > 0 && next_string.find_first_of("-") != 0)
            basename = next_string;
        }

      std::ostringstream filename;
      filename << basename << ".processor." << libMesh::global_processor_id();
      _ofstream = std::make_unique<std::ofstream>(filename.str().c_str());

      // Redirect, saving the original streambufs!
      out_buf = libMesh::out.rdbuf (_ofstream->rdbuf());
      err_buf = libMesh::err.rdbuf (_ofstream->rdbuf());
    }

  // redirect libMesh::out to nothing on all
  // other processors unless explicitly told
  // not to via the --keep-cout command-line argument.
  if (!libMesh::on_command_line ("--keep-cout"))
    if (libMesh::global_processor_id() != 0)
      libMesh::out.rdbuf (nullptr);

  // Similarly, the user can request to drop cerr on all non-0 ranks.
  // By default, errors are printed on all ranks, but this can lead to
  // interleaved/unpredictable outputs when doing parallel regression
  // testing, which this option is designed to support.
  if (libMesh::on_command_line ("--drop-cerr"))
    if (libMesh::global_processor_id() != 0)
      libMesh::err.rdbuf (nullptr);

  // Now that we've finished setting our stream buffers, let's wrap them (if they're non-null) in our thread-safe wrapper
  if (!libMesh::on_command_line ("--disable-thread-safe-output"))
    install_thread_buffered_sync();

  // Check command line to override printing
  // of reference count information.
  if (libMesh::on_command_line("--disable-refcount-printing"))
    ReferenceCounter::disable_print_counter_info();

#ifdef LIBMESH_ENABLE_EXCEPTIONS
  // Set our terminate handler to write stack traces in the event of a
  // crash
  _old_terminate_handler = std::set_terminate(libmesh_terminate_handler);
#endif


  if (libMesh::on_command_line("--enable-fpe"))
    libMesh::enableFPE(true);

  if (libMesh::on_command_line("--enable-segv"))
    libMesh::enableSEGV(true);

#if defined(LIBMESH_HAVE_HDF5) && !defined(_MSC_VER)
  // We may be running with ExodusII configured not to use HDF5 (in
  // which case user code which wants to lock files has to do flock()
  // itself) or with ExodusII configured to use HDF5 (which will
  // helpfully try to get an exclusive flock() itself, and then scream
  // and die if user code has already locked the file.  To get
  // consistent behavior, we need to disable file locking.  The only
  // reliable way I can see to do this is via an HDF5 environment
  // variable.
  //
  // If the user set this environment variable then we'll trust that
  // they know what they're doing.  If not then we'll set FALSE,
  // because that's a better default for us than the unset default.
  setenv("HDF5_USE_FILE_LOCKING", "FALSE", /*overwrite=false*/0);
#endif // LIBMESH_HAVE_HDF5

  // The library is now ready for use
  libMeshPrivateData::_is_initialized = true;


  // Make sure these work.  Library methods
  // depend on these being implemented properly,
  // so this is a good time to test them!
  libmesh_assert (libMesh::initialized());
  libmesh_assert (!libMesh::closed());
}



LibMeshInit::~LibMeshInit()
{
  // Every processor had better be ready to exit at the same time.
  // Even if we're not doing parallel_only debugging, we don't want
  // one processor to try to exit until all others are done working.

  // We could be destructing here because we're unwinding the stack
  // due to a thrown exception, though.  It's possible that an
  // application is catching exceptions outside of the LibMeshInit
  // scope, or that we're using a C++ compiler that does unwinding
  // for uncaught exceptions (the standard says whether to go straight
  // to terminate() or unwind first is "implementation-defined").  If
  // *that* is the case then we can't safely communicate with other
  // processors that might not all be unwinding too.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
  if (!std::uncaught_exceptions())
#endif
    this->comm().barrier();

  // We can't delete, finalize, etc. more than once without
  // reinitializing in between
  libmesh_exceptionless_assert(!libMesh::closed());

  // Delete reference counted singleton(s)
  Singleton::cleanup();

  // Clear the thread task manager we started
  task_scheduler.reset();

  // Force the \p ReferenceCounter to print
  // its reference count information.  This allows
  // us to find memory leaks.  By default the
  // \p ReferenceCounter only prints its information
  // when the last created object has been destroyed.
  // That does no good if we are leaking memory!
  ReferenceCounter::print_info ();


  // Print an informative message if we detect a memory leak
  if (ReferenceCounter::n_objects() != 0)
    {
      libMesh::err << "Memory leak detected!"
                   << std::endl;

#if !defined(LIBMESH_ENABLE_REFERENCE_COUNTING) || defined(NDEBUG)

      libMesh::err << "Compile in DEBUG mode with --enable-reference-counting"
                   << std::endl
                   << "for more information"
                   << std::endl;
#endif

    }

  //  print the perflog to individual processor's file.
  libMesh::perflog.print_log();

  // Now clear the logging object, we don't want it to print
  // a second time during the PerfLog destructor.
  libMesh::perflog.clear();

  // Reconnect the output streams
  // (don't do this, or we will get messages from objects
  //  that go out of scope after the following return)
  //std::cout.rdbuf(std::cerr.rdbuf());


  // Set the initialized() flag to false
  libMeshPrivateData::_is_initialized = false;

  cleanup_stream_buffers();

#ifdef LIBMESH_ENABLE_EXCEPTIONS
  // Reset the old terminate handler; maybe the user code wants to
  // keep doing C++ stuff after closing libMesh stuff.
  std::set_terminate(_old_terminate_handler);
#endif

#ifdef LIBMESH_HAVE_NETGEN
  nglib::Ng_Exit();
#endif

  if (libMesh::on_command_line("--enable-fpe"))
    libMesh::enableFPE(false);

#if defined(LIBMESH_HAVE_PETSC)
  // Let PETSc know about all the command line objects that we
  // consumed without asking them, so we don't get unused option
  // warnings about those.
  std::vector<std::string> cli_names = command_line_names();
  for (const auto & name : cli_names)
    if (!name.empty() && name[0] == '-')
      {
        // Newer PETSc can give me a double-free I'm having trouble
        // replicating; let's protect against trying to clear
        // already-used values
        PetscBool used = PETSC_FALSE;
        auto ierr = PetscOptionsUsed(NULL, name.c_str(), &used);
        CHKERRABORT(libMesh::GLOBAL_COMM_WORLD, ierr);
        if (used == PETSC_FALSE)
          {
            ierr = PetscOptionsClearValue(NULL, name.c_str());
            CHKERRABORT(libMesh::GLOBAL_COMM_WORLD, ierr);
          }
      }

  // Allow the user to bypass PETSc finalization
  if (!libMesh::on_command_line ("--disable-petsc")
#if defined(LIBMESH_HAVE_MPI)
      && !libMesh::on_command_line ("--disable-mpi")
#endif
      )
    {
      PetscErrorCode ierr = LIBMESH_PETSC_SUCCESS;
# if defined(LIBMESH_HAVE_SLEPC)
      if (libmesh_initialized_slepc)
        ierr = SlepcFinalize();
# else
      if (libmesh_initialized_petsc)
        ierr = PetscFinalize();
# endif
      CHKERRABORT(libMesh::GLOBAL_COMM_WORLD, ierr);
    }
#endif

#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
  _vtk_mpi_controller->Finalize(/*finalized_externally=*/1);
  _vtk_mpi_controller->Delete();
#endif

  delete this->_comm;

#if defined(LIBMESH_HAVE_MPI)
  // Allow the user to bypass MPI finalization
  if (!libMesh::on_command_line ("--disable-mpi"))
    {
      delete this->_timpi_init;
    }
#else
  delete this->_timpi_init;
#endif
}


PerfLog & LibMeshInit::perf_log()
{
  return libMesh::perflog;
}



void add_command_line_name(const std::string & name)
{
  // Users had better not be asking about an empty string
  libmesh_assert(!name.empty());

  static std::mutex command_line_names_mutex;
  std::scoped_lock lock(command_line_names_mutex);

  command_line_name_set.insert(name);
}



void add_command_line_names(const GetPot & getpot)
{
  for (auto & getter : {&GetPot::get_requested_arguments,
                        &GetPot::get_requested_variables,
                        &GetPot::get_requested_sections})
    for (const std::string & name : (getpot.*getter)())
      add_command_line_name(name);
}


std::vector<std::string> command_line_names()
{
  return std::vector<std::string>(command_line_name_set.begin(),
                                  command_line_name_set.end());
}



bool on_command_line (std::string arg)
{
  // Make sure the command line parser is ready for use.  If it's not,
  // then we'll have to treat the command line as empty, for maximum
  // compatibility with programs that don't use LibMeshInit but
  // indirectly (e.g. via error handling code) query the command line.
  if (!command_line.get())
    return false;

  // Keep track of runtime queries, for later
  add_command_line_name(arg);

  bool found_it = command_line->search(arg);

  if (!found_it)
    {
      // Try with all dashes instead of underscores
      std::replace(arg.begin(), arg.end(), '_', '-');
      found_it = command_line->search(arg);
    }

  if (!found_it)
    {
      // OK, try with all underscores instead of dashes
      auto name_begin = arg.begin();
      while (*name_begin == '-')
        ++name_begin;
      std::replace(name_begin, arg.end(), '-', '_');
      found_it = command_line->search(arg);
    }

  return found_it;
}



template <typename T>
T command_line_value (const std::string & name, T value)
{
  // Make sure the command line parser is ready for use
  libmesh_assert(command_line.get());

  // only if the variable exists in the file
  if (command_line->have_variable(name))
    {
      value = (*command_line)(name, value);

      // Keep track of runtime queries, for later.  GetPot splits
      // foo=bar into a separate name=value, so we can query for the
      // name, but as far as PETSc is concerned that's one CLI
      // argument.  We'll store it that way.
      const std::string stringvalue =
        (*command_line)(name, std::string());
      add_command_line_name(name+"="+stringvalue);
    }

  return value;
}

template <typename T>
T command_line_value (const std::vector<std::string> & names, T value)
{
  // Make sure the command line parser is ready for use
  libmesh_assert(command_line.get());

  // Keep track of runtime queries, for later
  for (const auto & entry : names)
    {
      // Keep track of runtime queries, for later.  GetPot splits
      // foo=bar into a separate name=value, so we can query for the
      // name, but as far as PETSc is concerned that's one CLI
      // argument.  We'll store it that way.
      const std::string stringvalue =
        (*command_line)(entry, std::string());
      add_command_line_name(entry+"="+stringvalue);
    }

  // Check for multiple options (return the first that matches)
  for (const auto & entry : names)
    if (command_line->have_variable(entry))
      {
        value = (*command_line)(entry, value);
        break;
      }

  return value;
}



template <typename T>
T command_line_next (std::string name, T value)
{
  // Make sure the command line parser is ready for use
  libmesh_assert(command_line.get());

  // Keep track of runtime queries, for later
  add_command_line_name(name);

  // on_command_line also puts the command_line cursor in the spot we
  // need
  if (on_command_line(name))
    return command_line->next(value);

  return value;
}



template <typename T>
void command_line_vector (const std::string & name, std::vector<T> & vec)
{
  // Make sure the command line parser is ready for use
  libmesh_assert(command_line.get());

  // Keep track of runtime queries, for later
  add_command_line_name(name);

  // only if the variable exists on the command line
  if (command_line->have_variable(name))
    {
      unsigned size = command_line->vector_variable_size(name);
      vec.resize(size);

      for (unsigned i=0; i<size; ++i)
        vec[i] = (*command_line)(name, vec[i], i);
    }
}


SolverPackage default_solver_package ()
{
  libmesh_assert (libMesh::initialized());

  static bool called = false;

  // Check the command line.  Since the command line is
  // unchanging it is sufficient to do this only once.
  if (!called)
    {
      called = true;

#ifdef LIBMESH_HAVE_PETSC
      if (libMesh::on_command_line ("--use-petsc"))
        libMeshPrivateData::_solver_package = PETSC_SOLVERS;
#endif

#ifdef LIBMESH_TRILINOS_HAVE_AZTECOO
      if (libMesh::on_command_line ("--use-trilinos") ||
          libMesh::on_command_line ("--disable-petsc"))
        libMeshPrivateData::_solver_package = TRILINOS_SOLVERS;
#endif

#ifdef LIBMESH_HAVE_EIGEN
      if (libMesh::on_command_line ("--use-eigen"  ) ||
#if defined(LIBMESH_HAVE_MPI)
          // If the user bypassed MPI, we disable PETSc and Trilinos
          // too
          libMesh::on_command_line ("--disable-mpi") ||
#endif
          libMesh::on_command_line ("--disable-petsc"))
        libMeshPrivateData::_solver_package = EIGEN_SOLVERS;
#endif

#ifdef LIBMESH_HAVE_LASPACK
      if (libMesh::on_command_line ("--use-laspack"  ) ||
#if defined(LIBMESH_HAVE_MPI)
          // If the user bypassed MPI, we disable PETSc and Trilinos
          // too
          libMesh::on_command_line ("--disable-mpi") ||
#endif
          libMesh::on_command_line ("--disable-petsc"))
        libMeshPrivateData::_solver_package = LASPACK_SOLVERS;
#endif

      if (libMesh::on_command_line ("--disable-laspack") &&
          libMesh::on_command_line ("--disable-trilinos") &&
          libMesh::on_command_line ("--disable-eigen") &&
          (
#if defined(LIBMESH_HAVE_MPI)
           // If the user bypassed MPI, we disable PETSc too
           libMesh::on_command_line ("--disable-mpi") ||
#endif
           libMesh::on_command_line ("--disable-petsc")))
        libMeshPrivateData::_solver_package = INVALID_SOLVER_PACKAGE;
    }


  return libMeshPrivateData::_solver_package;
}



//-------------------------------------------------------------------------------
template LIBMESH_EXPORT unsigned char  command_line_value<unsigned char>  (const std::string &, unsigned char);
template LIBMESH_EXPORT unsigned short command_line_value<unsigned short> (const std::string &, unsigned short);
template LIBMESH_EXPORT unsigned int   command_line_value<unsigned int>   (const std::string &, unsigned int);
template LIBMESH_EXPORT char           command_line_value<char>           (const std::string &, char);
template LIBMESH_EXPORT short          command_line_value<short>          (const std::string &, short);
template LIBMESH_EXPORT int            command_line_value<int>            (const std::string &, int);
template LIBMESH_EXPORT float          command_line_value<float>          (const std::string &, float);
template LIBMESH_EXPORT double         command_line_value<double>         (const std::string &, double);
template LIBMESH_EXPORT long double    command_line_value<long double>    (const std::string &, long double);
template LIBMESH_EXPORT std::string    command_line_value<std::string>    (const std::string &, std::string);

template LIBMESH_EXPORT unsigned char  command_line_value<unsigned char>  (const std::vector<std::string> &, unsigned char);
template LIBMESH_EXPORT unsigned short command_line_value<unsigned short> (const std::vector<std::string> &, unsigned short);
template LIBMESH_EXPORT unsigned int   command_line_value<unsigned int>   (const std::vector<std::string> &, unsigned int);
template LIBMESH_EXPORT char           command_line_value<char>           (const std::vector<std::string> &, char);
template LIBMESH_EXPORT short          command_line_value<short>          (const std::vector<std::string> &, short);
template LIBMESH_EXPORT int            command_line_value<int>            (const std::vector<std::string> &, int);
template LIBMESH_EXPORT float          command_line_value<float>          (const std::vector<std::string> &, float);
template LIBMESH_EXPORT double         command_line_value<double>         (const std::vector<std::string> &, double);
template LIBMESH_EXPORT long double    command_line_value<long double>    (const std::vector<std::string> &, long double);
template LIBMESH_EXPORT std::string    command_line_value<std::string>    (const std::vector<std::string> &, std::string);

template LIBMESH_EXPORT unsigned char  command_line_next<unsigned char>   (std::string, unsigned char);
template LIBMESH_EXPORT unsigned short command_line_next<unsigned short>  (std::string, unsigned short);
template LIBMESH_EXPORT unsigned int   command_line_next<unsigned int>    (std::string, unsigned int);
template LIBMESH_EXPORT char           command_line_next<char>            (std::string, char);
template LIBMESH_EXPORT short          command_line_next<short>           (std::string, short);
template LIBMESH_EXPORT int            command_line_next<int>             (std::string, int);
template LIBMESH_EXPORT float          command_line_next<float>           (std::string, float);
template LIBMESH_EXPORT double         command_line_next<double>          (std::string, double);
template LIBMESH_EXPORT long double    command_line_next<long double>     (std::string, long double);
template LIBMESH_EXPORT std::string    command_line_next<std::string>     (std::string, std::string);

template LIBMESH_EXPORT void           command_line_vector<unsigned char> (const std::string &, std::vector<unsigned char> &);
template LIBMESH_EXPORT void           command_line_vector<unsigned short>(const std::string &, std::vector<unsigned short> &);
template LIBMESH_EXPORT void           command_line_vector<unsigned int>  (const std::string &, std::vector<unsigned int> &);
template LIBMESH_EXPORT void           command_line_vector<char>          (const std::string &, std::vector<char> &);
template LIBMESH_EXPORT void           command_line_vector<short>         (const std::string &, std::vector<short> &);
template LIBMESH_EXPORT void           command_line_vector<int>           (const std::string &, std::vector<int> &);
template LIBMESH_EXPORT void           command_line_vector<float>         (const std::string &, std::vector<float> &);
template LIBMESH_EXPORT void           command_line_vector<double>        (const std::string &, std::vector<double> &);
template LIBMESH_EXPORT void           command_line_vector<long double>   (const std::string &, std::vector<long double> &);

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
template LIBMESH_EXPORT Real           command_line_value<Real>           (const std::string &, Real);
template LIBMESH_EXPORT Real           command_line_value<Real>           (const std::vector<std::string> &, Real);
template LIBMESH_EXPORT Real           command_line_next<Real>            (std::string, Real);
template LIBMESH_EXPORT void           command_line_vector<Real>          (const std::string &, std::vector<Real> &);
#endif

} // namespace libMesh
