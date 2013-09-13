// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/tree.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tree class method

// constructor
template <unsigned int N>
Tree<N>::Tree (const MeshBase& m,
	       const unsigned int target_bin_size,
	       const Trees::BuildType bt) :
  TreeBase(m),
  root(m,target_bin_size),
  build_type(bt)
{
  // Set the root node bounding box equal to the bounding
  // box for the entire domain.
  root.set_bounding_box (MeshTools::bounding_box(mesh));


  if (build_type == Trees::NODES)
    {
      // Add all the nodes to the root node.  It will
      // automagically build the tree for us.
      MeshBase::const_node_iterator       it  = mesh.nodes_begin();
      const MeshBase::const_node_iterator end = mesh.nodes_end();

      for (; it != end; ++it)
	root.insert (*it);

      // Now the tree contains the nodes.
      // However, we want element pointers, so here we
      // convert between the two.
      std::vector<std::vector<const Elem*> > nodes_to_elem;

      MeshTools::build_nodes_to_elem_map (mesh, nodes_to_elem);
      root.transform_nodes_to_elements (nodes_to_elem);
    }

  else if (build_type == Trees::ELEMENTS)
    {
      // Add all active elements to the root node.  It will
      // automatically build the tree for us.
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();


      for (; it != end; ++it)
	root.insert (*it);
    }
}



template <unsigned int N>
const Elem* Tree<N>::find_element(const Point& p) const
{
  return root.find_element(p);
}


// ------------------------------------------------------------
// Explicit Instantiations
template class Tree<2>;
template class Tree<4>;
template class Tree<8>;

} // namespace libMesh
