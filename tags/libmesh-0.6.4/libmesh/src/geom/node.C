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




// Local includes
#include "node.h"




// ------------------------------------------------------------
// Node class static member initialization
//const unsigned int Node::invalid_id = libMesh::invalid_uint;


bool Node::operator==(const DofObject& rhs) const
{
  // Cast rhs to a Node*
  const Node* rhs_node = dynamic_cast<const Node*>(&rhs);

  // If we can't cast to a Node* then rhs must be an Elem
  if(rhs_node == NULL)
    return false;

  // Explicitly calling the operator== defined in Point
  return this->Point::operator==(*rhs_node);
}



#ifdef LIBMESH_HAVE_MPI
MPI_Datatype Node::PackedNode::create_mpi_datatype ()
{
  MPI_Datatype packed_node_type;
  MPI_Datatype types[] = { MPI_UNSIGNED, MPI_REAL };
  int blocklengths[] = { 2, 3 };
  MPI_Aint displs[2];
  
  // create a Packed node and get the addresses of the elements.
  // this will properly handle id/pid getting padded, for example,
  // in which case id and x may not be 2*sizeof(unsigned int) apart.
  Node::PackedNode pn;
  
  MPI_Address (&pn.id, &displs[0]);
  MPI_Address (&pn.x,  &displs[1]);
  displs[1] -= displs[0];
  displs[0] = 0;
  
#if MPI_VERSION > 1
  MPI_Type_create_struct (2, blocklengths, displs, types, &packed_node_type);
#else
  MPI_Type_struct (2, blocklengths, displs, types, &packed_node_type);
#endif // #if MPI_VERSION > 1

  return packed_node_type;
}
#endif // #ifdef LIBMESH_HAVE_MPI
