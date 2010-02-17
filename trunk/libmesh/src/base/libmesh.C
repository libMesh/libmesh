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
# include "petsc_macro.h"
# include "slepc_macro.h"
# if defined(LIBMESH_HAVE_PETSC)
EXTERN_C_FOR_PETSC_BEGIN
#   include <petsc.h>
#   include <petscerror.h>
EXTERN_C_FOR_PETSC_END
# endif // #if defined(LIBMESH_HAVE_PETSC)
# if defined(LIBMESH_HAVE_SLEPC)
EXTERN_C_FOR_PETSC_BEGIN
#   include <slepc.h>
EXTERN_C_FOR_PETSC_END
# endif // #if defined(LIBMESH_HAVE_SLEPC)
#endif // #if defined(LIBMESH_HAVE_MPI)


// --------------------------------------------------------
// Local anonymous namespace to hold miscelaneous variables
namespace {
  AutoPtr<GetPot> command_line (NULL);
  AutoPtr<std::ofstream> _ofstream (NULL);
  // If std::cout and std::cerr are redirected, we need to
  // be a little careful and save the original streambuf objects,
  // replacing them in the destructor before program termination.
  std::streambuf* cout_buf (NULL);
  std::streambuf* cerr_buf (NULL);
  
  AutoPtr<Threads::task_scheduler_init> task_scheduler (NULL);
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



// ------------------------------------------------------------
// libMeshdata initialization
#ifdef LIBMESH_HAVE_MPI
MPI_Comm           libMesh::COMM_WORLD = MPI_COMM_NULL;
#endif

Parallel::Communicator Parallel::Communicator_World;

PerfLog            libMesh::perflog ("libMesh",
#ifdef LIBMESH_ENABLE_PERFORMANCE_LOGGING
				     true
#else
				     false
#endif
				     );


const Real         libMesh::pi = 3.1415926535897932384626433832795029L;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
const Number       libMesh::imaginary (0., 1.);
const Number       libMesh::zero      (0., 0.);
#else
const Number       libMesh::zero = 0.;
#endif

const unsigned int libMesh::invalid_uint = static_cast<unsigned int>(-1);



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
namespace libMesh {
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

  // Allow the user to bypass PETSc initialization
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
      
# if defined(LIBMESH_HAVE_PETSC)
      
      if (!libMesh::on_command_line ("--disable-petsc"))
	{
	  int ierr=0;
	  
	  PETSC_COMM_WORLD = libMesh::COMM_WORLD;

          // Check whether the calling program has already initialized
          // PETSc, and avoid duplicate Initialize/Finalize
          PetscTruth petsc_already_initialized;
          ierr = PetscInitialized(&petsc_already_initialized);
	         CHKERRABORT(libMesh::COMM_WORLD,ierr);
          if (petsc_already_initialized != PETSC_TRUE)
            libmesh_initialized_petsc = true;
#  if defined(LIBMESH_HAVE_SLEPC)

          // If SLEPc allows us to check whether the calling program
          // has already initialized it, we do that, and avoid
          // duplicate Initialize/Finalize.
          // We assume that SLEPc will handle PETSc appropriately,
          // which it does in the versions we've checked.
#   if !SLEPC_VERSION_LESS_THAN(2,3,3)
          if (!SlepcInitializeCalled)
#   endif
            {
	      ierr = SlepcInitialize  (&argc, &argv, NULL, NULL);
	             CHKERRABORT(libMesh::COMM_WORLD,ierr);
              libmesh_initialized_slepc = true;
            }
#  else
          if (libmesh_initialized_petsc)
            {
	      ierr = PetscInitialize (&argc, &argv, NULL, NULL);
	             CHKERRABORT(libMesh::COMM_WORLD,ierr);
            }
#  endif
	}
# endif
    }

  // Could we have gotten bad values from the above calls?
  libmesh_assert (libMeshPrivateData::_n_processors >  0);
  libmesh_assert (libMeshPrivateData::_processor_id >= 0);

#else

  // No MPI, can only be uniprocessor
  // libmesh_assert (libMeshPrivateData::_n_processors == 1);
  // libmesh_assert (libMeshPrivateData::_processor_id == 0);
  
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
  
  // redirect std::cout to nothing on all
  // other processors unless explicitly told
  // not to via the --keep-cout command-line argument.
  if (libMesh::processor_id() != 0)
    if (!libMesh::on_command_line ("--keep-cout"))
      std::cout.rdbuf (NULL);
  
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

#if defined(LIBMESH_HAVE_MPI)
  // Allow the user to bypass MPI finalization
  if (!libMesh::on_command_line ("--disable-mpi"))
    {
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
          std::cerr << "Uncaught exception - aborting" << std::endl;
          MPI_Abort(libMesh::COMM_WORLD,1);
        }
# if defined(LIBMESH_HAVE_PETSC) 

      // Allow the user to bypass PETSc finalization
      if (!libMesh::on_command_line ("--disable-petsc"))
	{
#  if defined(LIBMESH_HAVE_SLEPC)
          if (libmesh_initialized_slepc)
	    SlepcFinalize();
#  else
          if (libmesh_initialized_petsc)
	    PetscFinalize();
#  endif
	}
# endif
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
      std::cerr << "Memory leak detected!"
		<< std::endl;
      
#if !defined(LIBMESH_ENABLE_REFERENCE_COUNTING) || defined(NDEBUG)

      std::cerr << "Compile in DEBUG mode with --enable-reference-counting"
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
      std::cout.rdbuf (cout_buf);
      std::cerr.rdbuf (cerr_buf);
    }
  
  // Return the number of outstanding objects.
  // This is equivalent to return 0 if all of
  // the reference counted objects have been
  // deleted.
  return static_cast<int>(ReferenceCounter::n_objects());
}
}



#ifndef LIBMESH_HAVE_MPI
void libMesh::init (int &argc, char** & argv)
{
  libmesh_deprecated();  // Use LibMeshInit instead
  libMesh::_init(argc, argv);
}
#else
void libMesh::init (int &argc, char** & argv,
		    MPI_Comm COMM_WORLD_IN)
{
  libmesh_deprecated();  // Use LibMeshInit instead
  libMesh::_init(argc, argv, COMM_WORLD_IN);
}
#endif




int libMesh::close ()
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

  // Honor the --redirect-stdout command-line option.
  // When this is specified each processor sends
  // std::cout/std::cerr messages to
  // stdout.processor.####
  if (libMesh::on_command_line ("--redirect-stdout"))
    {
      char filechar[80];
      sprintf (filechar, "stdout.processor.%04d",
	       libMesh::processor_id());
      _ofstream.reset (new std::ofstream (filechar));
      // Redirect, saving the original streambufs!
      cout_buf = std::cout.rdbuf (_ofstream->rdbuf());
      cerr_buf = std::cerr.rdbuf (_ofstream->rdbuf());
    }
}
#endif

LibMeshInit::~LibMeshInit()
{
  libMesh::_close();
}




bool libMesh::on_command_line (const std::string& arg)
{
  // Make sure the command line parser is ready for use
  libmesh_assert (command_line.get() != NULL);
  
  return command_line->search (arg);
}



template <typename T>
T libMesh::command_line_value (const std::string &name, T value)
{
  // Make sure the command line parser is ready for use
  libmesh_assert (command_line.get() != NULL);
  
    // only if the variable exists in the file
    if (command_line->have_variable(name.c_str()))
      value = (*command_line)(name.c_str(), value);      

    return value;
}


SolverPackage libMesh::default_solver_package ()
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
template int          libMesh::command_line_value<int>         (const std::string&, int);
template Real         libMesh::command_line_value<Real>        (const std::string&, Real);
template std::string  libMesh::command_line_value<std::string> (const std::string&, std::string);
