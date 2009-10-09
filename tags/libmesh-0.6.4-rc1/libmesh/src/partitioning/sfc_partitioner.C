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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "libmesh_config.h"
#include "mesh_base.h"
#include "sfc_partitioner.h"
#include "libmesh_logging.h"
#include "elem.h"

#ifdef LIBMESH_HAVE_SFCURVES
  namespace Sfc {
    extern "C" {
#     include "sfcurves.h"
    }
  }
#else
#  include "linear_partitioner.h"
#endif


// ------------------------------------------------------------
// SFCPartitioner implementation
void SFCPartitioner::_do_partition (MeshBase& mesh,
				    const unsigned int n)
{
  
  libmesh_assert (n > 0);

  // Check for an easy return
  if (n == 1)
    {
      this->single_partition (mesh);
      return;
    }

// What to do if the sfcurves library IS NOT present
#ifndef LIBMESH_HAVE_SFCURVES

  libmesh_here();
  std::cerr << "ERROR: The library has been built without"    << std::endl
	    << "Space Filling Curve support.  Using a linear" << std::endl
	    << "partitioner instead!" << std::endl;

  LinearPartitioner lp;

  lp.partition (mesh, n);
  
// What to do if the sfcurves library IS present
#else

  START_LOG("sfc_partition()", "SFCPartitioner");

  const unsigned int n_active_elem = mesh.n_active_elem();
  const unsigned int n_elem        = mesh.n_elem();
  
  // the forward_map maps the active element id
  // into a contiguous block of indices
  std::vector<unsigned int>
    forward_map (n_elem, libMesh::invalid_uint);

  // the reverse_map maps the contiguous ids back
  // to active elements
  std::vector<Elem*> reverse_map (n_active_elem, NULL);
  
  int size = static_cast<int>(n_active_elem);
  std::vector<double> x      (size);
  std::vector<double> y      (size);
  std::vector<double> z      (size);
  std::vector<int>    table  (size);


  // We need to map the active element ids into a
  // contiguous range.
  {
//     active_elem_iterator       elem_it (mesh.elements_begin());
//     const active_elem_iterator elem_end(mesh.elements_end());

    MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
    const MeshBase::element_iterator elem_end = mesh.active_elements_end(); 

    unsigned int el_num = 0;

    for (; elem_it != elem_end; ++elem_it)
      {
	libmesh_assert ((*elem_it)->id() < forward_map.size());
	libmesh_assert (el_num           < reverse_map.size());
	
	forward_map[(*elem_it)->id()] = el_num;
	reverse_map[el_num]           = *elem_it;
	el_num++;
      }
    libmesh_assert (el_num == n_active_elem);
  }
  

  // Get the centroid for each active element
  {
//     const_active_elem_iterator       elem_it (mesh.const_elements_begin());
//     const const_active_elem_iterator elem_end(mesh.const_elements_end());

    MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
    const MeshBase::element_iterator elem_end = mesh.active_elements_end(); 
    
    for (; elem_it != elem_end; ++elem_it)
      {
	const Elem* elem = *elem_it;
	
	libmesh_assert (elem->id() < forward_map.size());
	
	const Point p = elem->centroid();

	x[forward_map[elem->id()]] = p(0);
	y[forward_map[elem->id()]] = p(1);
	z[forward_map[elem->id()]] = p(2);
      }
  }

  // build the space-filling curve
  if (_sfc_type == "Hilbert")
    Sfc::hilbert (&x[0], &y[0], &z[0], &size, &table[0]);
  
  else if (_sfc_type == "Morton")
    Sfc::morton  (&x[0], &y[0], &z[0], &size, &table[0]);
  
  else
    {
      libmesh_here();
      std::cerr << "ERROR: Unknown type: " << _sfc_type << std::endl
		<< " Valid types are"                   << std::endl
		<< "  \"Hilbert\""                      << std::endl
		<< "  \"Morton\""                       << std::endl
		<< " "                                  << std::endl
		<< "Proceeding with a Hilbert curve."   << std::endl;
      
      Sfc::hilbert (&x[0], &y[0], &z[0], &size, &table[0]);
    }

  
  // Assign the partitioning to the active elements
  {
//      {
//        std::ofstream out ("sfc.dat");
//        out << "variables=x,y,z" << std::endl;
//        out << "zone f=point" << std::endl;
    
//        for (unsigned int i=0; i<n_active_elem; i++)
//  	out << x[i] << " "
//  	    << y[i] << " "
//  	    << z[i] << std::endl;
//      }      
    
    const unsigned int blksize = (n_active_elem+n-1)/n; 

    for (unsigned int i=0; i<n_active_elem; i++)
      {
	libmesh_assert (static_cast<unsigned int>(table[i]-1) < reverse_map.size());
	  
	Elem* elem = reverse_map[table[i]-1];

	elem->processor_id() = i/blksize;
      }
  }
  
  STOP_LOG("sfc_partition()", "SFCPartitioner");  

#endif
  
}
