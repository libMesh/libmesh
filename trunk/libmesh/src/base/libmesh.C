// $Id: libmesh.C,v 1.10 2003-04-07 18:34:48 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include <fstream>
#include <math.h>


// Local includes
#include "libmesh.h"
#include "reference_counter.h"


#if defined(HAVE_PETSC)
#  if  !defined(USE_COMPLEX_NUMBERS)
/**
 * Petsc include files.  PETSc with complex numbers 
 * is actually C++.
 */
namespace Petsc
{
  extern "C"
  {
#   include <petsc.h>
  }
}
using namespace Petsc;
#  else
#    include <petsc.h>
#  endif
#elif defined(HAVE_MPI)
namespace Mpi
{
  extern "C"
  {
#    include <mpi.h>
  }
}
using namespace Mpi;
#endif



// ------------------------------------------------------------
// libMesh static member initializations
//const Real       libMesh::pi        = 3.14159265358979323846264;
const Real       libMesh::pi        (acos(-1.));

#if defined(USE_COMPLEX_NUMBERS)
const Number     libMesh::imaginary (0., 1.);
#endif


#if defined(USE_COMPLEX_NUMBERS)
const Number     libMesh::zero      (0., 0.);
#else
const Number     libMesh::zero      (0.);
#endif



PerfLog      libMesh::log ("libMesh",
#ifdef ENABLE_PERFORMANCE_LOGGING
			   true
#else
			   false
#endif
			   );


std::ostream *libMesh::_out            = NULL;
std::ostream *libMesh::_err            = NULL;
bool          libMesh::_is_initialized = false;




// ------------------------------------------------------------
// libMesh member functions
#if defined(HAVE_MPI)
void libMesh::init (int & argc, char** & argv)
#else
void libMesh::init (int &     , char** &     )
#endif
{
#if defined(HAVE_PETSC)

  PetscInitialize (&argc, &argv, NULL, NULL);
 
  // Get the number of processors and our local processor id  
  MPI_Comm_rank (PETSC_COMM_WORLD, &_processor_id);
  MPI_Comm_size (PETSC_COMM_WORLD, &_n_processors);

#elif defined(HAVE_MPI)

  MPI_Init (&argc, &argv);

  // Get the number of processors and our local processor id  
  MPI_Comm_rank (MPI_COMM_WORLD, &_processor_id);
  MPI_Comm_size (MPI_COMM_WORLD, &_n_processors);

#else

  // No PETSC or MPI, can only be uniprocessor
  assert (_n_processors == 1);
  assert (_processor_id == 0);

#endif

  // Could we have gotten bad values from the above calls?
  assert (_n_processors >  0);
  assert (_processor_id >= 0);


  // Open the output streams.
  if (processor_id() == 0)
    {
      _out = &std::cout;
      _err = &std::cerr;
    }
  else
    {
      _out = new std::ostream(NULL);
      _err = new std::ostream(NULL);
    }

  
  // The library is now ready for use
  _is_initialized = true;

  // Make sure these work.  Library methods
  // depend on these being implemented properly,
  // so this is a good time to test them!
  assert (initialized());
  assert (!closed());
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


  // Close the output streams.
  {
    if (processor_id() != 0)
      {
	delete _out;
	delete _err;
      }
    
    _out = NULL;
    _err = NULL;
  }

  
  // Set the initialized() flag to false
  _is_initialized = false;
  

  // Return the number of outstanding objects.
  // This is equivalent to return 0 if all of
  // the reference counted objects have been
  // deleted.
  return static_cast<int>(ReferenceCounter::n_objects());
}
