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


// LibMesh includes
#include "nemesis_io.h"
#include "nemesis_io_helper.h"
#include "parallel_mesh.h"

// ------------------------------------------------------------
// Nemesis_IO class members
Nemesis_IO::Nemesis_IO (ParallelMesh& mesh) :
  MeshInput<ParallelMesh> (mesh, /*is_parallel_format=*/true),
  //MeshOutput<ParallelMesh> (mesh, /*is_parallel_format=*/true)
  _verbose (false)
{
}


Nemesis_IO::~Nemesis_IO ()
{
}



void Nemesis_IO::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;

#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
  // Set the verbose flag in the helper object
  // as well.
  nemhelper.verbose(_verbose);
#endif
}



#if defined(LIBMESH_HAVE_EXODUS_API) && defined(LIBMESH_HAVE_NEMESIS_API)
void Nemesis_IO::read (const std::string& base_filename)
{
  // Construct a filename string for this processor.
  //
  // FIXME: This assumes you are reading in a mesh on exactly the
  // same number of processors it was written out on!!
  // This should be generalized at some point...
  std::ostringstream file_oss;

  file_oss << base_filename;
    
  if ( libMesh::n_processors() > 1 )
    file_oss << '.' << libMesh::n_processors() << '.' << libMesh::processor_id();

  // In the 1 processor case, assume the base_filename is the filename...
  
  std::cout << "Opening file: " << file_oss.str() << std::endl;
}

#else

void Nemesis_IO::read (const std::string& )
{
  std::cerr <<  "ERROR, Nemesis API is not defined!" << std::endl;
  libmesh_error();
}

#endif
