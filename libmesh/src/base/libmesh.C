// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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

 

// C++ includes
#include <iostream>
#include <fstream>

#ifdef LIBMESH_ENABLE_EXCEPTIONS
#include <exception>
#endif

// Local includes
#include "libmesh.h"
#include "auto_ptr.h"
#include "getpot.h"
#include "parallel.h"
#include "reference_counter.h"
#include "remote_elem.h"
#include "threads.h"



#if defined(LIBMESH_HAVE_MPI)
# include <mpi.h>
#endif // #if defined(LIBMESH_HAVE_MPI)

#if defined(LIBMESH_HAVE_PETSC)
# include "petsc_macro.h"
EXTERN_C_FOR_PETSC_BEGIN
# include <petsc.h>
# include <petscerror.h>
EXTERN_C_FOR_PETSC_END
# if defined(LIBMESH_HAVE_SLEPC)
#  include "slepc_macro.h"
EXTERN_C_FOR_PETSC_BEGIN
#  include <slepc.h>
EXTERN_C_FOR_PETSC_END
# endif // #if defined(LIBMESH_HAVE_SLEPC)
#endif // #if defined(LIBMESH_HAVE_PETSC)


// --------------------------------------------------------
// Local anonymous namespace to hold miscelaneous variables
namespace {

  using libMesh::AutoPtr;

  AutoPtr<GetPot> command_line (NULL);
  AutoPtr<std::ofstream> _ofstream (NULL);
  // If std::cout and std::cerr are redirected, we need to
  // be a little careful and save the original streambuf objects,
  // replacing them in the destructor before program termination.
  std::streambuf* out_buf (NULL);
  std::streambuf* err_buf (NULL);
  
  AutoPtr<libMesh::Threads::task_scheduler_init> task_scheduler (NULL);
#if defined(LIBMESH_HAVE_MPI)
  bool libmesh_initialized_mpi = false;
#endif
#if defined(LIBMESH_HAVE_PETSC)
  bool libmesh_initialized_petsc = false;
#endif
#if defined(LIBMESH_HAVE_SLEPC)
  bool libmesh_initialized_slepc = false;
#endif
}


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
// libMeshdata initialization
#ifdef LIBMESH_HAVE_MPI
MPI_Comm           COMM_WORLD = MPI_COMM_NULL;
#else
int                COMM_WORLD = 0;
#endif

OStreamProxy out(std::cout);
OStreamProxy err(std::cerr);

Parallel::Communicator Parallel::Communicator_World;

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
int           libMesh::libMeshPrivateData::_n_processors = 1;
int           libMesh::libMeshPrivateData::_processor_id = 0;
#endif
int           libMesh::libMeshPrivateData::_n_threads = 1; /* Threads::task_scheduler_init::automatic; */
bool          libMesh::libMeshPrivateData::_is_initialized = false;
SolverPackage libMesh::libMeshPrivateData::_solver_package =
#if   defined(LIBMESH_HAVE_PETSC)    // PETSc is the default
                                                       PETSC_SOLVERS;
#elif defined(LIBMESH_HAVE_TRILINOS) // Use Trilinos if PETSc isn't there
                                                       TRILINOS_SOLVERS;
#elif defined(LIBMESH_HAVE_LASPACK)  // Use LASPACK if neither are there
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
  // stack trace.
  libmesh_write_traceout();

  // The system terminate_handler may do useful things like printing
  // uncaught exception information, or the user may have created
  // their own terminate handler that we want to call.
  old_terminate_handler();
}
#endif



#ifndef LIBMESH_HAVE_MPI
void _init (int &argc, char** & argv)
#else
void _init (int &argc, char** & argv,
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
    libMesh::libMeshPrivateData::_n_threads = 
      libMesh::command_line_value ("--n_threads", 1);

    task_scheduler.reset (new Threads::task_scheduler_init(libMesh::n_threads()));
  }

  // Construct singletons who may be at risk of the
  // "static initialization order fiasco"
  //
  // RemoteElem depends on static reference counting data
  remote_elem = new RemoteElem();

#if defined(LIBMESH_HAVE_MPI)

  // Allow the user to bypass MPI initialization
  if (!libMesh::on_command_line ("--disable-mpi"))
    {
      // Check whether the calling program has already initialized
      // MPI, and avoid duplicate Init/Finalize
      int flag;
      MPI_Initialized (&flag);

      if (!flag)
	{
	  MPI_Init (&argc, &argv);
	  libmesh_initialized_mpi = true;
	}
      
      // Duplicate the input communicator for internal use
      MPI_Comm_dup (COMM_WORLD_IN, &libMesh::COMM_WORLD);

      // And get a Parallel::Communicator copy too, to use
      // as a default for that API
      Parallel::Communicator_World = COMM_WORLD_IN;

      //MPI_Comm_set_name not supported in at least SGI MPT's MPI implementation
      //MPI_Comm_set_name (libMesh::COMM_WORLD, "libMesh::COMM_WORLD");
      
      libMeshPrivateData::_processor_id = Parallel::Communicator_World.rank();
      libMeshPrivateData::_n_processors = Parallel::Communicator_World.size();
    }

  // Could we have gotten bad values from the above calls?
  libmesh_assert (libMeshPrivateData::_n_processors >  0);
  libmesh_assert (libMeshPrivateData::_processor_id >= 0);

  // Let's be sure we properly initialize on every processor at once:
  parallel_only();

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
	  
      PETSC_COMM_WORLD = libMesh::COMM_WORLD;

      // Check whether the calling program has already initialized
      // PETSc, and avoid duplicate Initialize/Finalize
      PetscBool petsc_already_initialized;
      ierr = PetscInitialized(&petsc_already_initialized);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
      if (petsc_already_initialized != PETSC_TRUE)
        libmesh_initialized_petsc = true;
# if defined(LIBMESH_HAVE_SLEPC)

      // If SLEPc allows us to check whether the calling program
      // has already initialized it, we do that, and avoid
      // duplicate Initialize/Finalize.
      // We assume that SLEPc will handle PETSc appropriately,
      // which it does in the versions we've checked.
#  if !SLEPC_VERSION_LESS_THAN(2,3,3)
      if (!SlepcInitializeCalled)
#  endif
        {
          ierr = SlepcInitialize  (&argc, &argv, NULL, NULL);
	         CHKERRABORT(libMesh::COMM_WORLD,ierr);
          libmesh_initialized_slepc = true;
        }
# else
      if (libmesh_initialized_petsc)
        {
          ierr = PetscInitialize (&argc, &argv, NULL, NULL);
	         CHKERRABORT(libMesh::COMM_WORLD,ierr);
        }
# endif
    }
#endif

  // Re-parse the command-line arguments.  Note that PETSc and MPI
  // initialization above may have removed command line arguments
  // that are not relevant to this application in the above calls.
  // We don't want a false-positive by detecting those arguments.
  command_line->parse_command_line (argc, argv);

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
      std::ostream* newout = new std::ostream(std::cout.rdbuf());
      libMesh::out = *newout;
      std::ostream* newerr = new std::ostream(std::cerr.rdbuf());
      libMesh::err = *newerr;
    }
  
  // Honor the --redirect-stdout command-line option.
  // When this is specified each processor sends
  // libMesh::out/libMesh::err messages to
  // stdout.processor.####
  if (libMesh::on_command_line ("--redirect-stdout"))
    {
      char filechar[80];
      sprintf (filechar, "stdout.processor.%04d",
	       libMesh::processor_id());
      _ofstream.reset (new std::ofstream (filechar));
      // Redirect, saving the original streambufs!
      out_buf = libMesh::out.rdbuf (_ofstream->rdbuf());
      err_buf = libMesh::err.rdbuf (_ofstream->rdbuf());
    }

  // redirect libMesh::out to nothing on all
  // other processors unless explicitly told
  // not to via the --keep-cout command-line argument.
  if (libMesh::processor_id() != 0)
    if (!libMesh::on_command_line ("--keep-cout"))
      libMesh::out.rdbuf (NULL);

  // Check command line to override printing
  // of reference count information.
  if(libMesh::on_command_line("--disable-refcount-printing") )
    ReferenceCounter::disable_print_counter_info();

#ifdef LIBMESH_ENABLE_EXCEPTIONS
  // Set our terminate handler to write stack traces in the event of a
  // crash
  old_terminate_handler = std::set_terminate(libmesh_terminate_handler);
#endif
  
  // The library is now ready for use
  libMeshPrivateData::_is_initialized = true;

  
  // Make sure these work.  Library methods
  // depend on these being implemented properly,
  // so this is a good time to test them!
  libmesh_assert (libMesh::initialized());
  libmesh_assert (!libMesh::closed());
}



int _close ()
{
  // We can't delete, finalize, etc. more than once without
  // reinitializing in between
  libmesh_assert(!libMesh::closed());

  // Delete reference counted singleton(s)
  delete remote_elem;

  // Clear the thread task manager we started
  task_scheduler.reset();

  // Let's be sure we properly close on every processor at once:
  parallel_only();

#if defined(LIBMESH_HAVE_MPI)
  // We may be here in only one process,
  // because an uncaught libmesh_error() exception
  // called the LibMeshInit destructor.
  //
  // If that's the case, we need to MPI_Abort(),
  // not just wait for other processes that
  // might never get to MPI_Finalize()
  if (libmesh_initialized_mpi &&
      std::uncaught_exception())
    {
      libMesh::err << "Uncaught exception - aborting" << std::endl;
      if (libMesh::on_command_line ("--disable-mpi"))
        exit(1);
      else
        MPI_Abort(libMesh::COMM_WORLD,1);
    }
#endif

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

#if defined(LIBMESH_HAVE_MPI)
  // Allow the user to bypass MPI finalization
  if (!libMesh::on_command_line ("--disable-mpi"))
    {
      Parallel::Communicator_World.clear();
      MPI_Comm_free (&libMesh::COMM_WORLD);

      if (libmesh_initialized_mpi)
	MPI_Finalize();
    }
#endif

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


  // Reconnect the output streams
  // (don't do this, or we will get messages from objects
  //  that go out of scope after the following return)
  //std::cout.rdbuf(std::cerr.rdbuf());

  
  // Set the initialized() flag to false
  libMeshPrivateData::_is_initialized = false;

  if (libMesh::on_command_line ("--redirect-stdout"))
    {
      // Before handing back the std stream buffers, print the
      // perflog to the individual processor's files.
      libMesh::perflog.print_log();

      // Now clear the logging object, we don't want it to print
      // a second time during the PerfLog destructor.
      libMesh::perflog.clear();
      
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

  // Return the number of outstanding objects.
  // This is equivalent to return 0 if all of
  // the reference counted objects have been
  // deleted.
  return static_cast<int>(ReferenceCounter::n_objects());
}



#ifndef LIBMESH_HAVE_MPI
void init (int &argc, char** & argv)
{
  libmesh_deprecated();  // Use LibMeshInit instead
  libMesh::_init(argc, argv);
}
#else
void init (int &argc, char** & argv,
		    MPI_Comm COMM_WORLD_IN)
{
  libmesh_deprecated();  // Use LibMeshInit instead
  libMesh::_init(argc, argv, COMM_WORLD_IN);
}
#endif




int close ()
{
  libmesh_deprecated();  // Use LibMeshInit instead
  return libMesh::_close();
}



#ifndef LIBMESH_HAVE_MPI
LibMeshInit::LibMeshInit (int &argc, char** & argv)
{
  libMesh::_init(argc, argv);
}
#else
LibMeshInit::LibMeshInit (int &argc, char** & argv,
		          MPI_Comm COMM_WORLD_IN)
{
  libMesh::_init(argc, argv, COMM_WORLD_IN);
}
#endif

LibMeshInit::~LibMeshInit()
{
  libMesh::_close();
}




bool on_command_line (const std::string& arg)
{
  // Make sure the command line parser is ready for use
  libmesh_assert (command_line.get() != NULL);
  
  return command_line->search (arg);
}



template <typename T>
T command_line_value (const std::string &name, T value)
{
  // Make sure the command line parser is ready for use
  libmesh_assert (command_line.get() != NULL);
  
    // only if the variable exists in the file
    if (command_line->have_variable(name.c_str()))
      value = (*command_line)(name.c_str(), value);      

    return value;
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

#ifdef LIBMESH_HAVE_TRILINOS
      if (libMesh::on_command_line ("--use-trilinos") ||
	  libMesh::on_command_line ("--disable-petsc"))
	libMeshPrivateData::_solver_package = TRILINOS_SOLVERS;
#endif
      
#ifdef LIBMESH_HAVE_LASPACK
      if (libMesh::on_command_line ("--use-laspack"  ) ||
	  libMesh::on_command_line ("--disable-petsc"))
	libMeshPrivateData::_solver_package = LASPACK_SOLVERS;
#endif

      if (libMesh::on_command_line ("--disable-laspack") &&
	  libMesh::on_command_line ("--disable-trilinos") &&
	  libMesh::on_command_line ("--disable-petsc"))
	libMeshPrivateData::_solver_package = INVALID_SOLVER_PACKAGE;
    }
  
  
  return libMeshPrivateData::_solver_package;  
}



//-------------------------------------------------------------------------------
template int          command_line_value<int>         (const std::string&, int);
template Real         command_line_value<Real>        (const std::string&, Real);
template std::string  command_line_value<std::string> (const std::string&, std::string);

} // namespace libMesh
