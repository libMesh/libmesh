// $Id: libmesh.C,v 1.2 2003-02-17 04:05:48 benkirk Exp $

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
PerfLog      libMesh::log ("libMesh",
#ifdef ENABLE_PERFORMANCE_LOGGING
			   true
#else
			   false
#endif
			   );

bool         libMesh::_is_initialized = false;
int          libMesh::_n_processors   = 1;
int          libMesh::_processor_id   = 0;

// ------------------------------------------------------------
// libMesh member functions
void libMesh::init (int & argc, char** & argv)
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
  // That does not good if we are leaking memory!
  ReferenceCounter::print_info ();


  _is_initialized = false;
  
  // Return 0 if everything has gone well.
  return 0;
}
