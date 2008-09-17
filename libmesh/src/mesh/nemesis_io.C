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
#include "parallel.h"

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
  if (_verbose)
    {
      std::cout << "[" << libMesh::processor_id() << "] ";
      std::cout << "Reading Nemesis file on processor: " << libMesh::processor_id() << std::endl;
    }
  
  // Construct a filename string for this processor.
  //
  // FIXME: This assumes you are reading in a mesh on exactly the
  // same number of processors it was written out on!!
  // This should be generalized at some point...
  std::ostringstream file_oss;

  file_oss << base_filename;
    
  // In the 1 processor case, assume the base_filename is the filename...
  if ( libMesh::n_processors() > 1 )
    file_oss << '.' << libMesh::n_processors() << '.' << libMesh::processor_id();

  std::cout << "Opening file: " << file_oss.str() << std::endl;

  // Make it a little easier to access the ExodusII helper object...
  ExodusII_IO_Helper& ex2helper = nemhelper.ex2helper;
  
  // Open the Exodus file
  ex2helper.open(file_oss.str().c_str());

  // Only split files will have additional Nemesis information
  if ( libMesh::n_processors() > 1 )
    {
      // Get global information: number of nodes, elems, blocks, nodesets and sidesets
      nemhelper.get_init_global();

      // Get "load balance" information.  This includes the number of internal & border
      // nodes and elements as well as the number of communication maps.
      nemhelper.get_loadbal_param();

      // Communicate n_internal_elems and n_border_elems to all processors.  This
      // determines the global element numbering.  
      // std::vector<int> num_internal_elems_all (libMesh::n_processors());
      // Parallel::allgather(nemhelper.num_internal_elems, num_internal_elems_all);
      //       if (_verbose)
      // 	{
      // 	  std::cout << "[" << libMesh::processor_id() << "] ";
      // 	  for (unsigned int j=0; j<libMesh::n_processors(); ++j)
      // 	    std::cout << "num_internal_elems_all["
      // 		      <<	j << "]="
      // 		      << num_internal_elems_all[j]
      // 		      << " ";
      // 	  std::cout << std::endl;
      // 	}

      // The call to get_loadbal_param() gets 7 pieces of information.  We allgather
      // these now across all processors to determine some global numberings...
      std::vector<int> all_loadbal_data ( 7 );
      all_loadbal_data[0] = nemhelper.num_internal_nodes;
      all_loadbal_data[1] = nemhelper.num_border_nodes;
      all_loadbal_data[2] = nemhelper.num_external_nodes;
      all_loadbal_data[3] = nemhelper.num_internal_elems;
      all_loadbal_data[4] = nemhelper.num_border_elems;
      all_loadbal_data[5] = nemhelper.num_node_cmaps;
      all_loadbal_data[6] = nemhelper.num_elem_cmaps;
      Parallel::allgather(all_loadbal_data);
      
      if (_verbose)
	{
	  // Ex: print out the number of internal nodes on all processors
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  for (unsigned int j=0,c=0; c<libMesh::n_processors(); j+=7)
	    std::cout << "num_internal_nodes["
		      << c++ << "]="
		      << all_loadbal_data[j]
		      << " ";
	  std::cout << std::endl;
	}

      
      // Assertion: The sum of the border and internal elements on all processors
      // should equal nemhelper.num_elems_global
#ifndef NDEBUG
      int sum_internal_elems=0, sum_border_elems=0;
      for (unsigned int j=3,c=0; c<libMesh::n_processors(); j+=7,++c)
	sum_internal_elems += all_loadbal_data[j];

      for (unsigned int j=4,c=0; c<libMesh::n_processors(); j+=7,++c)
	sum_border_elems += all_loadbal_data[j];

      if (_verbose)
	{
	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "sum_internal_elems=" << sum_internal_elems << std::endl;

	  std::cout << "[" << libMesh::processor_id() << "] ";
	  std::cout << "sum_border_elems=" << sum_border_elems << std::endl;
	}

      libmesh_assert(sum_internal_elems+sum_border_elems == nemhelper.num_elems_global);
#endif

      // Compute my_elem_offset, the amount by which to offset the local elem numbering
      // on my processor.
      unsigned int my_elem_offset = 0;
      for (unsigned int i=0; i<libMesh::processor_id(); ++i)
	my_elem_offset += (all_loadbal_data[7*i + 3]+  // num_internal_elems, proc i
			   all_loadbal_data[7*i + 4]); // num_border_elems, proc i
      std::cout << "[" << libMesh::processor_id() << "] ";
      std::cout << "my_elem_offset=" << my_elem_offset << std::endl;
      
      // Get a reference to the ParallelMesh.  
      ParallelMesh& mesh = this->mesh();

      // Local information: Read the standard Exodus header
      ex2helper.read_header();
      ex2helper.print_header();

      // Be sure number of dimensions is equal to the number of dimensions in the mesh supplied.
      libmesh_assert(static_cast<unsigned int>(ex2helper.num_dim) == mesh.mesh_dimension());

      // Read nodes from the exodus file: this fills the ex2helper.x,y,z arrays.
      ex2helper.read_nodes();
    } // end if ( libMesh::n_processors() > 1 )

  else
    {
      // Do something different when running in serial??
      libmesh_error();
    }

  
  
}

#else

void Nemesis_IO::read (const std::string& )
{
  std::cerr <<  "ERROR, Nemesis API is not defined!" << std::endl;
  libmesh_error();
}

#endif
