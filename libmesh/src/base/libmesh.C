// $Id: libmesh.C,v 1.24 2004-08-06 16:48:40 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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

// Local includes
#include "libmesh.h"
#include "auto_ptr.h"
#include "getpot.h"
#include "reference_counter.h"



#if defined(HAVE_PETSC)


#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petsc.h>
}
#else
# include <petsc.h>
#endif


#elif defined(HAVE_MPI)
# include <mpi.h>
#endif



// ------------------------------------------------------------
// Local anonymous namespace to hold a GetPot object
namespace {
  AutoPtr<GetPot> command_line (NULL);
}



// ------------------------------------------------------------
// libMesh data initialization
PerfLog            libMesh::perflog ("libMesh",
#ifdef ENABLE_PERFORMANCE_LOGGING
				     true
#else
				     false
#endif
				     );

#ifdef USE_COMPLEX_NUMBERS
const Number       libMesh::imaginary (0., 1.);
#endif

const Real         libMesh::pi = 3.1415926535897932384626433832795029L;

#ifdef USE_COMPLEX_NUMBERS
const Number       libMesh::zero (0., 0.);
#else
const Number       libMesh::zero = 0.;
#endif

const unsigned int libMesh::invalid_uint = static_cast<unsigned int>(-1);



// ------------------------------------------------------------
// libMesh::libMeshPrivateData data initialization
int           libMesh::libMeshPrivateData::_n_processors = 1;
int           libMesh::libMeshPrivateData::_processor_id = 0;
bool          libMesh::libMeshPrivateData::_is_initialized = false;
SolverPackage libMesh::libMeshPrivateData::_solver_package =
#if   defined(HAVE_PETSC)   // PETSc is the default
                                                       PETSC_SOLVERS;
#elif defined(HAVE_LASPACK) // Use LASPACK if PETSc isn't there
                                                       LASPACK_SOLVERS;
#else                       // No valid linear solver package at compile time
                                                       INVALID_SOLVER_PACKAGE;
#endif



// ------------------------------------------------------------
// libMesh functions
void libMesh::init (int &argc, char** & argv)
{
  // should _not_ be initialized already.
  assert (!libMesh::initialized());

  // The following line is an optimization when simultaneous
  // C and C++ style access to output streams is not required.
  // The amount of benefit which occurs is probably implementation
  // defined, and may be nothing.  On the other hand, I have seen
  // some IO tests where IO peformance improves by a factor of two.
  std::ios_base::sync_with_stdio(false);
  
  // Build a command-line parser.
  command_line.reset(new GetPot (argc, argv));
  
#if defined(HAVE_PETSC)

  // Allow the user to bypass PETSc initialization
  if (!libMesh::on_command_line ("--disable-petsc"))      
    PetscInitialize (&argc, &argv, NULL, NULL);
  // Allow the user to bypass MPI initialization
  else if (!libMesh::on_command_line ("--disable-mpi"))
    MPI_Init (&argc, &argv);
   
  // Get the number of processors and our local processor id
  if (!libMesh::on_command_line ("--disable-petsc") &&
      !libMesh::on_command_line ("--disable-mpi"))
    {
      MPI_Comm_rank (PETSC_COMM_WORLD, &libMeshPrivateData::_processor_id);
      MPI_Comm_size (PETSC_COMM_WORLD, &libMeshPrivateData::_n_processors);
    }
  
#elif defined(HAVE_MPI)

  // Allow the user to bypass MPI initialization
  if (!libMesh::on_command_line ("--disable-mpi"))
    {
      MPI_Init (&argc, &argv);

      // Get the number of processors and our local processor id  
      MPI_Comm_rank (MPI_COMM_WORLD, &libMeshPrivateData::_processor_id);
      MPI_Comm_size (MPI_COMM_WORLD, &libMeshPrivateData::_n_processors);
    }
  
#else

  // No PETSC or MPI, can only be uniprocessor
  assert (libMeshPrivateData::_n_processors == 1);
  assert (libMeshPrivateData::_processor_id == 0);

#endif
  
  // Could we have gotten bad values from the above calls?
  assert (libMeshPrivateData::_n_processors >  0);
  assert (libMeshPrivateData::_processor_id >= 0);

  // Re-parse the command-line arguments.  Note that PETSc and MPI
  // initialization above may have removed command line arguments
  // that are not relevant to this application in the above calls.
  // We don't want a false-positive by detecting those arguments.
  command_line->parse_command_line (argc, argv);
    
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
  assert (libMesh::initialized());
  assert (!libMesh::closed());
}



int libMesh::close ()
{

#if defined(HAVE_PETSC)

  // Allow the user to bypass PETSc finalization
  if (!libMesh::on_command_line ("--disable-petsc"))      
    PetscFinalize ();
  // Allow the user to bypass MPI finalization
  else if (!libMesh::on_command_line ("--disable-mpi"))
    MPI_Finalize ();

#elif defined(HAVE_MPI)

  // Allow the user to bypass MPI finalization
  if (!libMesh::on_command_line ("--disable-mpi"))
    MPI_Finalize();
  
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
      
#if !defined(ENABLE_REFERENCE_COUNTING) || defined(NDEBUG)

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
  

  // Return the number of outstanding objects.
  // This is equivalent to return 0 if all of
  // the reference counted objects have been
  // deleted.
  return static_cast<int>(ReferenceCounter::n_objects());
}



bool libMesh::on_command_line (const std::string& arg)
{
  // Make sure the command line parser is ready for use
  assert (command_line.get() != NULL);
  
  return command_line->search (arg);
}



SolverPackage libMesh::default_solver_package ()
{
  assert (libMesh::initialized());
  
  static bool called = false;

  // Check the command line.  Since the command line is
  // unchanging it is sufficient to do this only once.
  if (!called)
    {
      called = true;

#ifdef HAVE_PETSC
      if (libMesh::on_command_line ("--use-petsc"))
	libMeshPrivateData::_solver_package = PETSC_SOLVERS;
#endif
      
#ifdef HAVE_LASPACK
      if (libMesh::on_command_line ("--use-laspack"  ) ||
	  libMesh::on_command_line ("--disable-petsc"))
	libMeshPrivateData::_solver_package = LASPACK_SOLVERS;
#endif
      
    }
  
  
  return libMeshPrivateData::_solver_package;  
}
