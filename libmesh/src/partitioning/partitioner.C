// $Id: partitioner.C,v 1.8 2004-11-14 18:51:59 jwpeterson Exp $

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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "partitioner.h"
#include "mesh_base.h"
#include "elem.h"


// ------------------------------------------------------------
// Partitioner implementation


void Partitioner::partition (MeshBase& mesh,
			     const unsigned int n)
{
  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n;

  // Call the partitioning function
  this->_do_partition(mesh,n);
}





void Partitioner::repartition (MeshBase& mesh,
			       const unsigned int n)
{
  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n;
  
  // Call the partitioning function
  this->_do_repartition(mesh,n);
}





void Partitioner::single_partition (MeshBase& mesh)
{
  // Loop over all the elements and assign them to processor 0.
//   elem_iterator       elem_it (mesh.elements_begin());
//   const elem_iterator elem_end(mesh.elements_end());

  MeshBase::element_iterator       elem_it  = mesh.elements_begin();
  const MeshBase::element_iterator elem_end = mesh.elements_end(); 

  for ( ; elem_it != elem_end; ++elem_it)
    (*elem_it)->set_processor_id() = 0;
}
