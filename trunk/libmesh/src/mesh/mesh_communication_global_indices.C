// $Id: mesh_communication.C 2373 2007-11-08 23:04:25Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "libmesh_config.h"
#include "libmesh_common.h"
#include "libmesh_logging.h"
#include "mesh_base.h"
#include "mesh_tools.h"
#include "mesh_communication.h"
#include "parallel.h"
#include "parallel_sort.h"
#include "elem.h"
#ifdef HAVE_LIBHILBERT
#  include "hilbert.h"
#endif
  
#ifdef HAVE_LIBHILBERT
namespace {
  std::ostream& operator << (std::ostream& os, const Hilbert::BitVecType& t)
  {
    int r = t.rackCount()-1;
    
    while (r >= 0)
      os << "_" << t.racks()[r--];
    
    return os;
  }
}
#endif // #ifdef HAVE_LIBHILBERT



// ------------------------------------------------------------
// MeshCommunication class members
void MeshCommunication::find_global_indices (MeshBase& mesh) const
{
#ifdef HAVE_LIBHILBERT

  START_LOG ("find_global_indices()", "MeshCommunication");

  // This method determines partition-agnostic global indices
  // for nodes and elements.

  // Algorithm:
  // (1) compute the Hilbert key for each local node/element
  // (2) perform a parallel sort of the Hilbert key
  // (3) get the min/max value on each processor.
  // (4) determine the position in the global ranking for
  //     each local object

  // Global bounding box
  MeshTools::BoundingBox bbox =
    MeshTools::bounding_box (mesh);

  const Hilbert::inttype max_inttype = static_cast<Hilbert::inttype>(-1);
  const unsigned int sizeof_inttype = sizeof(Hilbert::inttype);

  //-------------------------------------------------------------
  // (1) compute Hilbert keys
  std::vector<Hilbert::BitVecType>
    node_keys, elem_keys;
  
  {
    CFixBitVec icoords[3];
    Hilbert::BitVecType hilbert_index(3*sizeof(double)*sizeof(Hilbert::inttype));
    //const unsigned int n_racks = hilbert_index.rackCount();
    node_keys.reserve (mesh.n_local_nodes());  
    elem_keys.reserve (mesh.n_local_elem()); 
    
    // Nodes first
    {
      MeshBase::const_node_iterator       it  = mesh.local_nodes_begin();
      const MeshBase::const_node_iterator end = mesh.local_nodes_end();

      for (; it != end; ++it)
	{
	  const Node* node = (*it);
	  assert (node != NULL);

	  const Point& p = *node;

	  const double // (x,y,z) in [0,1]^3
	    x = ((bbox.first(0) == bbox.second(0)) ? 0. :
		 (p(0)-bbox.first(0))/(bbox.second(0)-bbox.first(0))),
	  
	    y = ((bbox.first(1) == bbox.second(1)) ? 0. :
		 (p(1)-bbox.first(1))/(bbox.second(1)-bbox.first(1))),
	  
	    z = ((bbox.first(2) == bbox.second(2)) ? 0. :
		 (p(2)-bbox.first(2))/(bbox.second(2)-bbox.first(2)));
	
	  // (iccords) in [0,max_inttype]^3
	  icoords[0] = static_cast<Hilbert::inttype>(x*max_inttype);
	  icoords[1] = static_cast<Hilbert::inttype>(y*max_inttype);
	  icoords[2] = static_cast<Hilbert::inttype>(z*max_inttype);

	  // Compute the Hilbert Index
	  Hilbert::coordsToIndex (icoords, 8*sizeof_inttype, 3, hilbert_index);

	  node_keys.push_back(hilbert_index);
	}
    }
    
    // Elements next
    {
      MeshBase::const_element_iterator       it  = mesh.local_elements_begin();
      const MeshBase::const_element_iterator end = mesh.local_elements_end();

      for (; it != end; ++it)
	{
	  const Elem* elem = (*it);
	  assert (elem != NULL);

	  const Point p = elem->centroid();
	
	  const double // (x,y,z) in [0,1]^3
	    x = ((bbox.first(0) == bbox.second(0)) ? 0. :
		 (p(0)-bbox.first(0))/(bbox.second(0)-bbox.first(0))),
	  
	    y = ((bbox.first(1) == bbox.second(1)) ? 0. :
		 (p(1)-bbox.first(1))/(bbox.second(1)-bbox.first(1))),
	  
	    z = ((bbox.first(2) == bbox.second(2)) ? 0. :
		 (p(2)-bbox.first(2))/(bbox.second(2)-bbox.first(2)));
	
	  // (iccords) in [0,max_inttype]^3
	  icoords[0] = static_cast<Hilbert::inttype>(x*max_inttype);
	  icoords[1] = static_cast<Hilbert::inttype>(y*max_inttype);
	  icoords[2] = static_cast<Hilbert::inttype>(z*max_inttype);

	  // Compute the Hilbert Index
	  Hilbert::coordsToIndex (icoords, 8*sizeof_inttype, 3, hilbert_index);

	  elem_keys.push_back(hilbert_index);
	}
    }    
  } // done computing Hilbert keys
  
 
  //-------------------------------------------------------------
  // (2) parallel sort the Hilbert keys
  {
    Parallel::Sort<Hilbert::BitVecType> elem_sorter (elem_keys);
    elem_sorter.sort();

    for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
      {
	MPI_Barrier (libMesh::COMM_WORLD);

	if (libMesh::processor_id() == pid)
	  {
	    std::cerr << "PID [" << libMesh::processor_id() << "] size="
		      << elem_sorter.bin().size() << std::endl;

	    for (unsigned int i=0; i<elem_sorter.bin().size(); i++)
	      std::cerr << elem_sorter.bin()[i] << '\n';

	    std::cerr << std::endl;
	  }
	sleep (2);
	MPI_Barrier (libMesh::COMM_WORLD);
      }

    Parallel::Sort<Hilbert::BitVecType> node_sorter (node_keys);
    node_sorter.sort();
  }

  STOP_LOG ("find_global_indices()", "MeshCommunication");

#endif // HAVE_LIBHILBERT
}
