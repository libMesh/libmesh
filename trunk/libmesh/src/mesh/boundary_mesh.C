// $Id: boundary_mesh.C,v 1.14 2005-02-22 22:17:39 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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

// Local includes
#include "boundary_mesh.h"



// ------------------------------------------------------------
// BoundaryMesh class member functions
BoundaryMesh::BoundaryMesh(unsigned int d) :
  Mesh(d)
{
}



BoundaryMesh::~BoundaryMesh()
{
  this->clear();
}



void BoundaryMesh::clear()
{
  // Reset the number of subdomains and partitions
  _n_sbd   = 1;
  _n_parts = 1;
  
  //   // Clear the elements data structure
  //   for (unsigned int e=0; e<_elements.size(); e++)
  //     if (_elements[e] != NULL)
  //       {
  // 	delete _elements[e];
  // 	_elements[e] = NULL;
  //       }	    
  
  // Delete all the elements.
  MeshBase::element_iterator it        = this->elements_begin();
  const MeshBase::element_iterator end = this->elements_end();

  for (; it != end; ++it)
    this->delete_elem(*it);

  // Note: in the future, we can't assume that there will be
  // a vector of elements to clear!
  // _elements.clear();
  
  // Don't delete the nodes here... They are simply pointers
  // to the nodes in \p MeshBase that will be deleted by another
  // class.

  // Note: in the future, we can't assume that there will be
  // a vector of nodes to clear!
  // _nodes.clear();
}







