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


// Local includes
#include "libmesh/libmesh.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/getpot.h"
#include "libmesh/parallel.h"
#include "libmesh/reference_counter.h"
#include "libmesh/libmesh_singleton.h"
#include "libmesh/remote_elem.h"
#include "libmesh/threads.h"
#include "libmesh/print_trace.h"


// C/C++ includes
#include <iostream>
#include <fstream>

#ifdef LIBMESH_ENABLE_EXCEPTIONS
#include <exception>
#endif

#ifdef LIBMESH_HAVE_OPENMP
#include <omp.h>
#endif

#include "signal.h"


// floating-point exceptions
#ifdef LIBMESH_HAVE_FENV_H
#  include <fenv.h>
#endif
#ifdef LIBMESH_HAVE_XMMINTRIN_H
#  include <xmmintrin.h>
#endif


#if defined(LIBMESH_HAVE_MPI)
# include "libmesh/ignore_warnings.h"
# include <mpi.h>
# include "libmesh/restore_warnings.h"
#endif // #if defined(LIBMESH_HAVE_MPI)

#if defined(LIBMESH_HAVE_PETSC)
# include "libmesh/petsc_macro.h"
# include <petsc.h>
# include <petscerror.h>
#if !PETSC_RELEASE_LESS_THAN(3,3,0)
#include "libmesh/petscdmlibmesh.h"
#endif
# if defined(LIBMESH_HAVE_SLEPC)
// Ignore unused variable warnings from SLEPc
#  include "libmesh/ignore_warnings.h"
#  include "libmesh/slepc_macro.h"
#  include <slepc.h>
#  include "libmesh/restore_warnings.h"
# endif // #if defined(LIBMESH_HAVE_SLEPC)
#endif // #if defined(LIBMESH_HAVE_PETSC)

// If we're using MPI and VTK has been detected, we need to do some
// MPI initialize/finalize stuff for VTK.
#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
#include "libmesh/ignore_warnings.h"
# include "vtkMPIController.h"
#include "libmesh/restore_warnings.h"
#endif

// --------------------------------------------------------
// Local anonymous namespace to hold miscelaneous bits
namespace {

libMesh::UniquePtr<GetPot> command_line;
libMesh::UniquePtr<std::ofstream> _ofstream;
// If std::cout and std::cerr are redirected, we need to
// be a little careful and save the original streambuf objects,
// replacing them in the destructor before program termination.
std::streambuf * out_buf (libmesh_nullptr);
std::streambuf * err_buf (libmesh_nullptr);

libMesh::UniquePtr<libMesh::Threads::task_scheduler_init> task_scheduler;
#if defined(LIBMESH_HAVE_MPI)
bool libmesh_initialized_mpi = false;
#endif
#if defined(LIBMESH_HAVE_PETSC)
bool libmesh_initialized_petsc = false;
#endif
#if defined(LIBMESH_HAVE_SLEPC)
bool libmesh_initialized_slepc = false;
#endif



/**
 * Floating point exception handler -- courtesy of Cody Permann & MOOSE team
 */
void libmesh_handleFPE(int /*signo*/, siginfo_t * info, void * /*context*/)
{
  libMesh::err << std::endl;
  libMesh::err << "Floating point exception signaled (";
  switch (info->si_code)
    {
    case FPE_INTDIV: libMesh::err << "integer divide by zero"; break;
    case FPE_INTOVF: libMesh::err << "integer overflow"; break;
    case FPE_FLTDIV: libMesh::err << "floating point divide by zero"; break;
    case FPE_FLTOVF: libMesh::err << "floating point overflow"; break;
    case FPE_FLTUND: libMesh::err << "floating point underflow"; break;
    case FPE_FLTRES: libMesh::err << "floating point inexact result"; break;
    case FPE_FLTINV: libMesh::err << "invalid floating point operation"; break;
    case FPE_FLTSUB: libMesh::err << "subscript out of range"; break;
    default:         libMesh::err << "unrecognized"; break;
    }
  libMesh::err << ")!" << std::endl;

  libmesh_error_msg("\nTo track this down, compile in debug mode, then in gdb do:\n" \
                    << "  break libmesh_handleFPE\n"                    \
                    << "  run ...\n"                                    \
                    << "  bt");
}


void libmesh_handleSEGV(int /*signo*/, siginfo_t * info, void * /*context*/)
{
  libMesh::err << std::endl;
  libMesh::err << "Segmentation fault exception signaled (";
  switch (info->si_code)
    {
    case SEGV_MAPERR: libMesh::err << "Address not mapped"; break;
    case SEGV_ACCERR: libMesh::err << "Invalid permissions"; break;
    default:         libMesh::err << "unrecognized"; break;
    }
  libMesh::err << ")!" << std::endl;

  libmesh_error_msg("\nTo track this down, compile in debug mode, then in gdb do:\n" \
                    << "  break libmesh_handleSEGV\n"                    \
                    << "  run ...\n"                                    \
                    << "  bt");
}
}



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
#ifndef LIBMESH_DISABLE_COMMWORLD
MPI_Comm           COMM_WORLD = MPI_COMM_NULL;
#endif
MPI_Comm           GLOBAL_COMM_WORLD = MPI_COMM_NULL;
#else
#ifndef LIBMESH_DISABLE_COMMWORLD
int                COMM_WORLD = 0;
#endif
int                GLOBAL_COMM_WORLD = 0;
#endif

#ifdef LIBMESH_DISABLE_COMMWORLD
Parallel::FakeCommunicator CommWorld;
Parallel::FakeCommunicator & Parallel::Communicator_World = CommWorld;
#else
Parallel::Communicator CommWorld;
Parallel::Communicator & Parallel::Communicator_World = CommWorld;
#endif


OStreamProxy out(std::cout);
OStreamProxy err(std::cerr);

bool warned_about_auto_ptr(false);

PerfLog            perflog ("libMesh",
#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
                            true
#else
                            false
#endif
                            );


// const Real         pi = 3.1415926535897932384626433832795029L;

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


#ifdef LIBMESH_ENABLE_EXCEPTIONS
std::terminate_handler old_terminate_handler;

void libmesh_terminate_handler()
{
  // If this got called then we're probably crashing; let's print a
  // stack trace.  The trace files that are ultimately written depend on:
  // 1.) Who throws the exception.
  // 2.) Whether the C++ runtime unwinds the stack before the
  //     terminate_handler is called (this is implementation defined).
  //
  // The various cases are summarized in the table below:
  //
  //                        | libmesh exception | other exception
  //                        -------------------------------------
  // stack unwinds          |        A          |       B
  // stack does not unwind  |        C          |       D
  //
  // Case A: There will be two stack traces in the file: one "useful"
  //         one, and one nearly empty one due to stack unwinding.
  // Case B: You will get one nearly empty stack trace (not great, Bob!)
  // Case C: You will get two nearly identical stack traces, ignore one of them.
  // Case D: You will get one useful stack trace.
  //
  // Cases A and B (where the stack unwinds when an exception leaves
  // main) appear to be non-existent in practice.  I don't have a
  // definitive list, but the stack does not unwind for GCC on either
  // Mac or Linux.  I think there's good reasons for this behavior too:
  // it's much easier to get a stack trace when the stack doesn't
  // unwind, for example.
  libMesh::write_traceout();

  // We may care about performance data pre-crash; it would be sad to
  // throw that away.
  libMesh::perflog.print_log();
  libMesh::perflog.clear();

  // If we have MPI and it has been initialized, we need to be sure
  // and call MPI_Abort instead of std::abort, so that the parallel
  // job can die nicely.
#if defined(LIBMESH_HAVE_MPI)
  int mpi_initialized;
  MPI_Initialized (&mpi_initialized);

  if (mpi_initialized)
    MPI_Abort(libMesh::GLOBAL_COMM_WORLD, 1);
  else
#endif
    // The system terminate_handler may do useful things like printing
    // uncaught exception information, or the user may have created
    // their own terminate handler that we want to call.
    old_terminate_handler();
}
#endif



#ifndef LIBMESH_HAVE_MPI
LibMeshInit::LibMeshInit (int argc, const char * const * argv)
#else
LibMeshInit::LibMeshInit (int argc, const char * const * argv,
                          MPI_Comm COMM_WORLD_IN)
#endif
{
  // should _not_ be initialized already.
  libmesh_assert (!libMesh::initialized());

  // Build a command-line parser.
  command_line.reset (new GetPot (argc, argv));

  // Disable performance logging upon request
  {
    if (libMesh::on_command_line ("--disable-perflog"))
      libMesh::perflog.disable_logging();
  }

  // Build a task scheduler
  {
    // Get the requested number of threads, defaults to 1 to avoid MPI and
    // multithreading competition.  If you would like to use MPI and multithreading
    // at the same time then (n_mpi_processes_per_node)x(n_threads) should be the
    //  number of processing cores per node.
    std::vector<std::string> n_threads(2);
    n_threads[0] = "--n_threads";
    n_threads[1] = "--n-threads";
    libMesh::libMeshPrivateData::_n_threads =
      libMesh::command_line_value (n_threads, 1);

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

    task_scheduler.reset (new Threads::task_scheduler_init(libMesh::n_threads()));
  }

  // Construct singletons who may be at risk of the
  // "static initialization order fiasco"
  Singleton::setup();

  // Make sure the construction worked
  libmesh_assert(remote_elem);

#if defined(LIBMESH_HAVE_MPI)

  // Allow the user to bypass MPI initialization
  if (!libMesh::on_command_line ("--disable-mpi"))
    {
      // Check whether the calling program has already initialized
      // MPI, and avoid duplicate Init/Finalize
      int flag;
      libmesh_call_mpi(MPI_Initialized (&flag));

      if (!flag)
        {
#if MPI_VERSION > 1
          int mpi_thread_provided;
          const int mpi_thread_requested = libMesh::n_threads() > 1 ?
            MPI_THREAD_FUNNELED :
            MPI_THREAD_SINGLE;

          libmesh_call_mpi
            (MPI_Init_thread (&argc, const_cast<char ***>(&argv),
                              mpi_thread_requested, &mpi_thread_provided));

          if ((libMesh::n_threads() > 1) &&
              (mpi_thread_provided < MPI_THREAD_FUNNELED))
            {
              libmesh_warning("Warning: MPI failed to guarantee MPI_THREAD_FUNNELED\n"
                              << "for a threaded run.\n"
                              << "Be sure your library is funneled-thread-safe..."
                              << std::endl);

              // Ideally, if an MPI stack tells us it's unsafe for us
              // to use threads, we shouldn't use threads.
              // In practice, we've encountered one MPI stack (an
              // mvapich2 configuration) that returned
              // MPI_THREAD_SINGLE as a proper warning, two stacks
              // that handle MPI_THREAD_FUNNELED properly, and two
              // current stacks plus a couple old stacks that return
              // MPI_THREAD_SINGLE but support libMesh threaded runs
              // anyway.

              // libMesh::libMeshPrivateData::_n_threads = 1;
              // task_scheduler.reset (new Threads::task_scheduler_init(libMesh::n_threads()));
            }
#else
          if (libMesh::libMeshPrivateData::_n_threads > 1)
            {
              libmesh_warning("Warning: using MPI1 for threaded code.\n" <<
                              "Be sure your library is funneled-thread-safe..." <<
                              std::endl);
            }

          libmesh_call_mpi
            (MPI_Init (&argc, const_cast<char ***>(&argv)));
#endif
          libmesh_initialized_mpi = true;
        }

      // Duplicate the input communicator for internal use
      // And get a Parallel::Communicator copy too, to use
      // as a default for that API
      this->_comm = COMM_WORLD_IN;

      libMesh::GLOBAL_COMM_WORLD = COMM_WORLD_IN;

#ifndef LIBMESH_DISABLE_COMMWORLD
      libMesh::COMM_WORLD = COMM_WORLD_IN;
      Parallel::Communicator_World = COMM_WORLD_IN;
#endif

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
#if MPI_VERSION > 1
          libmesh_call_mpi
            (MPI_Comm_create_errhandler(libMesh_MPI_Handler, &libmesh_errhandler));
          libmesh_call_mpi
            (MPI_Comm_set_errhandler(libMesh::GLOBAL_COMM_WORLD, libmesh_errhandler));
          libmesh_call_mpi
            (MPI_Comm_set_errhandler(MPI_COMM_WORLD, libmesh_errhandler));
#else
          libmesh_call_mpi
            (MPI_Errhandler_create(libMesh_MPI_Handler, &libmesh_errhandler));
          libmesh_call_mpi
            (MPI_Errhandler_set(libMesh::GLOBAL_COMM_WORLD, libmesh_errhandler));
          libmesh_call_mpi
            (MPI_Errhandler_set(MPI_COMM_WORLD, libmesh_errhandler));
#endif // #if MPI_VERSION > 1
        }
    }

  // Could we have gotten bad values from the above calls?
  libmesh_assert_greater (libMeshPrivateData::_n_processors, 0);

  // The cast_int already tested _processor_id>=0
  // libmesh_assert_greater_equal (libMeshPrivateData::_processor_id, 0);

  // Let's be sure we properly initialize on every processor at once:
  libmesh_parallel_only(this->comm());

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
      int ierr=0;

      PETSC_COMM_WORLD = libMesh::GLOBAL_COMM_WORLD;

      // Check whether the calling program has already initialized
      // PETSc, and avoid duplicate Initialize/Finalize
      PetscBool petsc_already_initialized;
      ierr = PetscInitialized(&petsc_already_initialized);
      CHKERRABORT(libMesh::GLOBAL_COMM_WORLD,ierr);
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
          ierr = SlepcInitialize  (&argc, const_cast<char ***>(&argv), libmesh_nullptr, libmesh_nullptr);
          CHKERRABORT(libMesh::GLOBAL_COMM_WORLD,ierr);
          libmesh_initialized_slepc = true;
        }
# else
      if (libmesh_initialized_petsc)
        {
          ierr = PetscInitialize (&argc, const_cast<char ***>(&argv), libmesh_nullptr, libmesh_nullptr);
          CHKERRABORT(libMesh::GLOBAL_COMM_WORLD,ierr);
        }
# endif
#if !PETSC_RELEASE_LESS_THAN(3,3,0)
      // Register the reference implementation of DMlibMesh
#if PETSC_RELEASE_LESS_THAN(3,4,0)
      ierr = DMRegister(DMLIBMESH, PETSC_NULL, "DMCreate_libMesh", DMCreate_libMesh); CHKERRABORT(libMesh::GLOBAL_COMM_WORLD,ierr);
#else
      ierr = DMRegister(DMLIBMESH, DMCreate_libMesh); CHKERRABORT(libMesh::GLOBAL_COMM_WORLD,ierr);
#endif

#endif
    }
#endif

#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
  // Do MPI initializtion for VTK.
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
  command_line.reset (new GetPot (argc, argv));

  // The following line is an optimization when simultaneous
  // C and C++ style access to output streams is not required.
  // The amount of benefit which occurs is probably implementation
  // defined, and may be nothing.  On the other hand, I have seen
  // some IO tests where IO peformance improves by a factor of two.
  if (!libMesh::on_command_line ("--sync-with-stdio"))
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
      _ofstream.reset (new std::ofstream (filename.str().c_str()));

      // Redirect, saving the original streambufs!
      out_buf = libMesh::out.rdbuf (_ofstream->rdbuf());
      err_buf = libMesh::err.rdbuf (_ofstream->rdbuf());
    }

  // redirect libMesh::out to nothing on all
  // other processors unless explicitly told
  // not to via the --keep-cout command-line argument.
  if (libMesh::global_processor_id() != 0)
    if (!libMesh::on_command_line ("--keep-cout"))
      libMesh::out.rdbuf (libmesh_nullptr);

  // Similarly, the user can request to drop cerr on all non-0 ranks.
  // By default, errors are printed on all ranks, but this can lead to
  // interleaved/unpredictable outputs when doing parallel regression
  // testing, which this option is designed to support.
  if (libMesh::global_processor_id() != 0)
    if (libMesh::on_command_line ("--drop-cerr"))
      libMesh::err.rdbuf (libmesh_nullptr);

  // Check command line to override printing
  // of reference count information.
  if(libMesh::on_command_line("--disable-refcount-printing") )
    ReferenceCounter::disable_print_counter_info();

#ifdef LIBMESH_ENABLE_EXCEPTIONS
  // Set our terminate handler to write stack traces in the event of a
  // crash
  old_terminate_handler = std::set_terminate(libmesh_terminate_handler);
#endif


  if (libMesh::on_command_line("--enable-fpe"))
    libMesh::enableFPE(true);

  if (libMesh::on_command_line("--enable-segv"))
    libMesh::enableSEGV(true);

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
  // This would be a libmesh_parallel_only() function, except that
  // libmesh_parallel_only() uses libmesh_assert() which throws an
  // exception() which causes compilers to scream about exceptions
  // inside destructors.

  // Even if we're not doing parallel_only debugging, we don't want
  // one processor to try to exit until all others are done working.
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

#ifdef LIBMESH_ENABLE_EXCEPTIONS
  // Reset the old terminate handler; maybe the user code wants to
  // keep doing C++ stuff after closing libMesh stuff.
  std::set_terminate(old_terminate_handler);
#endif


  if (libMesh::on_command_line("--enable-fpe"))
    libMesh::enableFPE(false);

#if defined(LIBMESH_HAVE_PETSC)
  // Allow the user to bypass PETSc finalization
  if (!libMesh::on_command_line ("--disable-petsc")
#if defined(LIBMESH_HAVE_MPI)
      && !libMesh::on_command_line ("--disable-mpi")
#endif
      )
    {
# if defined(LIBMESH_HAVE_SLEPC)
      if (libmesh_initialized_slepc)
        SlepcFinalize();
# else
      if (libmesh_initialized_petsc)
        PetscFinalize();
# endif
    }
#endif

#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
  _vtk_mpi_controller->Finalize(/*finalized_externally=*/1);
  _vtk_mpi_controller->Delete();
#endif

#if defined(LIBMESH_HAVE_MPI)
  // Allow the user to bypass MPI finalization
  if (!libMesh::on_command_line ("--disable-mpi"))
    {
      this->_comm.clear();
#ifndef LIBMESH_DISABLE_COMMWORLD
      Parallel::Communicator_World.clear();
#endif

      if (libmesh_initialized_mpi)
        {
          // We can't just libmesh_assert here because destructor,
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
    }
#endif
}



/**
 * Toggle floating point exceptions -- courtesy of Cody Permann & MOOSE team
 */
void enableFPE(bool on)
{
#if !defined(LIBMESH_HAVE_FEENABLEEXCEPT) && defined(LIBMESH_HAVE_XMMINTRIN_H) && !defined(__SUNPRO_CC)
  static int flags = 0;
#endif

  if (on)
    {
      struct sigaction new_action, old_action;

#ifdef LIBMESH_HAVE_FEENABLEEXCEPT
      feenableexcept(FE_DIVBYZERO | FE_INVALID);
#elif  LIBMESH_HAVE_XMMINTRIN_H
#  ifndef __SUNPRO_CC
      flags = _MM_GET_EXCEPTION_MASK();           // store the flags
      _MM_SET_EXCEPTION_MASK(flags & ~_MM_MASK_INVALID);
#  endif
#endif


      // Set up the structure to specify the new action.
      new_action.sa_sigaction = libmesh_handleFPE;
      sigemptyset (&new_action.sa_mask);
      new_action.sa_flags = SA_SIGINFO;

      sigaction (SIGFPE, libmesh_nullptr, &old_action);
      if (old_action.sa_handler != SIG_IGN)
        sigaction (SIGFPE, &new_action, libmesh_nullptr);
    }
  else
    {
#ifdef LIBMESH_HAVE_FEDISABLEEXCEPT
      fedisableexcept(FE_DIVBYZERO | FE_INVALID);
#elif  LIBMESH_HAVE_XMMINTRIN_H
#  ifndef __SUNPRO_CC
      _MM_SET_EXCEPTION_MASK(flags);
#  endif
#endif
      signal(SIGFPE, SIG_DFL);
    }
}


// Enable handling of SIGSEGV by libMesh
// (potentially instead of PETSc)
void enableSEGV(bool on)
{
  static struct sigaction old_action;
  static bool was_on = false;

  if (on)
    {
      struct sigaction new_action;
      was_on = true;

      // Set up the structure to specify the new action.
      new_action.sa_sigaction = libmesh_handleSEGV;
      sigemptyset (&new_action.sa_mask);
      new_action.sa_flags = SA_SIGINFO;

      sigaction (SIGSEGV, &new_action, &old_action);
    }
  else if (was_on)
    {
      was_on = false;
      sigaction (SIGSEGV, &old_action, libmesh_nullptr);
    }
}



bool on_command_line (const std::string & arg)
{
  // Make sure the command line parser is ready for use
  libmesh_assert(command_line.get());

  return command_line->search (arg);
}



template <typename T>
T command_line_value (const std::string & name, T value)
{
  // Make sure the command line parser is ready for use
  libmesh_assert(command_line.get());

  // only if the variable exists in the file
  if (command_line->have_variable(name.c_str()))
    value = (*command_line)(name.c_str(), value);

  return value;
}

template <typename T>
T command_line_value (const std::vector<std::string> & name, T value)
{
  // Make sure the command line parser is ready for use
  libmesh_assert(command_line.get());

  // Check for multiple options (return the first that matches)
  for (std::vector<std::string>::const_iterator i=name.begin(); i != name.end(); ++i)
    if (command_line->have_variable(i->c_str()))
      {
        value = (*command_line)(i->c_str(), value);
        break;
      }

  return value;
}



template <typename T>
T command_line_next (const std::string & name, T value)
{
  // Make sure the command line parser is ready for use
  libmesh_assert(command_line.get());

  if (command_line->search(1, name.c_str()))
    value = command_line->next(value);

  return value;
}



template <typename T>
void command_line_vector (const std::string & name, std::vector<T> & vec)
{
  // Make sure the command line parser is ready for use
  libmesh_assert(command_line.get());

  // only if the variable exists on the command line
  if (command_line->have_variable(name.c_str()))
    {
      unsigned size = command_line->vector_variable_size(name.c_str());
      vec.resize(size);

      for (unsigned i=0; i<size; ++i)
        vec[i] = (*command_line)(name.c_str(), vec[i], i);
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
template int          command_line_value<int>         (const std::string &, int);
template float        command_line_value<float>       (const std::string &, float);
template double       command_line_value<double>      (const std::string &, double);
template long double  command_line_value<long double> (const std::string &, long double);
template std::string  command_line_value<std::string> (const std::string &, std::string);

template int          command_line_value<int>         (const std::vector<std::string> &, int);
template float        command_line_value<float>       (const std::vector<std::string> &, float);
template double       command_line_value<double>      (const std::vector<std::string> &, double);
template long double  command_line_value<long double> (const std::vector<std::string> &, long double);
template std::string  command_line_value<std::string> (const std::vector<std::string> &, std::string);

template int          command_line_next<int>         (const std::string &, int);
template float        command_line_next<float>       (const std::string &, float);
template double       command_line_next<double>      (const std::string &, double);
template long double  command_line_next<long double> (const std::string &, long double);
template std::string  command_line_next<std::string> (const std::string &, std::string);

template void         command_line_vector<int>         (const std::string &, std::vector<int> &);
template void         command_line_vector<float>       (const std::string &, std::vector<float> &);
template void         command_line_vector<double>      (const std::string &, std::vector<double> &);
template void         command_line_vector<long double> (const std::string &, std::vector<long double> &);

} // namespace libMesh
