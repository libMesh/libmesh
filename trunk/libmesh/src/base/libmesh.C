// $Id: libmesh.C,v 1.17 2003-09-25 21:46:55 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
#include <math.h>


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
PerfLog            libMesh::log ("libMesh", (ENABLE_PERFORMANCE_LOGGING == 1));

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
#if   defined(HAVE_PETSC)
                                                       PETSC_SOLVERS;
#elif defined(HAVE_LASPACK)
                                                       LASPACK_SOLVERS;
#else
                                                       INVALID_SOLVER_PACKAGE;
#endif



// ------------------------------------------------------------
// libMesh functions
void libMesh::init (int &argc, char** & argv)
{
  // should _not_ be initialized already.
  assert (!libMesh::initialized());
  
#if defined(HAVE_PETSC)

  PetscInitialize (&argc, &argv, NULL, NULL);
 
  // Get the number of processors and our local processor id  
  MPI_Comm_rank (PETSC_COMM_WORLD, &libMeshPrivateData::_processor_id);
  MPI_Comm_size (PETSC_COMM_WORLD, &libMeshPrivateData::_n_processors);

#elif defined(HAVE_MPI)

  MPI_Init (&argc, &argv);

  // Get the number of processors and our local processor id  
  MPI_Comm_rank (MPI_COMM_WORLD, &libMeshPrivateData::_processor_id);
  MPI_Comm_size (MPI_COMM_WORLD, &libMeshPrivateData::_n_processors);

#else

  // No PETSC or MPI, can only be uniprocessor
  assert (libMeshPrivateData::_n_processors == 1);
  assert (libMeshPrivateData::_processor_id == 0);

#endif

  // Could we have gotten bad values from the above calls?
  assert (libMeshPrivateData::_n_processors >  0);
  assert (libMeshPrivateData::_processor_id >= 0);
    
  // Parse the command-line arguments
  command_line.reset(new GetPot (argc, argv));
  
  // redirect std::cout to nothing on all
  // other processors unless explicitly told
  // not to via the --keep-cout command-line argument.
  if (libMesh::processor_id() != 0)
    if (!libMesh::on_command_line("--keep-cout"))
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

  PetscFinalize ();

#elif defined(HAVE_MPI)

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
  assert (libMesh::initialized());
  
  return command_line->search (arg);
}



SolverPackage libMesh::default_solver_package ()
{
  assert (libMesh::initialized());
  
  static bool called = false;

  // Check the command line.  Since the command line is
  // unchangingt it is sufficient to do this only once.
  if (!called)
    {
      called = true;

#ifdef HAVE_PETSC
      if (libMesh::on_command_line ("--use-petsc"))
	libMeshPrivateData::_solver_package = PETSC_SOLVERS;
#endif

#ifdef HAVE_LASPACK
      if (libMesh::on_command_line ("--use-laspack"))
	libMeshPrivateData::_solver_package = LASPACK_SOLVERS;
#endif
      
    }


  return libMeshPrivateData::_solver_package;  
}
