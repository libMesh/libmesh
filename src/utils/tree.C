// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
Tree<N>::Tree (const MeshBase & m,
               unsigned int target_bin_size,
               const Trees::BuildType bt) :
  TreeBase(m),
  root(m,target_bin_size),
  build_type(bt)
{
  // Set the root node bounding box equal to the bounding
  // box for the entire domain.
  root.set_bounding_box (MeshTools::create_bounding_box(mesh));

  if (build_type == Trees::NODES)
    {
      // Add all the nodes to the root node.  It will
      // automagically build the tree for us.
      for (const auto & node : mesh.node_ptr_range())
        {
#ifndef NDEBUG
          bool node_was_inserted =
#endif
            root.insert (node);
          libmesh_assert(node_was_inserted);
        }

      // Now the tree contains the nodes.
      // However, we want element pointers, so here we
      // convert between the two.
      std::vector<std::vector<const Elem *>> nodes_to_elem;

      MeshTools::build_nodes_to_elem_map (mesh, nodes_to_elem);
      root.transform_nodes_to_elements (nodes_to_elem);
    }

  else if (build_type == Trees::ELEMENTS)
    {
      // Add all active elements to the root node.  It will
      // automatically build the tree for us.
      for (const auto & elem : mesh.active_element_ptr_range())
        {
#ifndef NDEBUG
          bool elem_was_inserted =
#endif
            root.insert (elem);
          libmesh_assert(elem_was_inserted);
        }
    }

  else if (build_type == Trees::LOCAL_ELEMENTS)
    {
      // Add all active, local elements to the root node.  It will
      // automatically build the tree for us.
      for (const auto & elem : mesh.active_local_element_ptr_range())
        {
#ifndef NDEBUG
          bool elem_was_inserted =
#endif
            root.insert (elem);
          libmesh_assert(elem_was_inserted);
        }
    }

  else
    libmesh_error_msg("Unknown build_type = " << build_type);
}



// copy-constructor is not implemented
template <unsigned int N>
Tree<N>::Tree (const Tree<N> & other_tree) :
  TreeBase   (other_tree),
  root       (other_tree.root),
  build_type (other_tree.build_type)
{
  libmesh_not_implemented();
}






template <unsigned int N>
void Tree<N>::print_nodes(std::ostream & my_out) const
{
  my_out << "Printing nodes...\n";
  root.print_nodes(my_out);
}



template <unsigned int N>
void Tree<N>::print_elements(std::ostream & my_out) const
{
  my_out << "Printing elements...\n";
  root.print_elements(my_out);
}



template <unsigned int N>
const Elem *
Tree<N>::find_element (const Point & p,
                       const std::set<subdomain_id_type> * allowed_subdomains,
                       Real relative_tol) const
{
  return root.find_element(p, allowed_subdomains, relative_tol);
}



template <unsigned int N>
const Elem *
Tree<N>::operator() (const Point & p,
                     const std::set<subdomain_id_type> * allowed_subdomains,
                     Real relative_tol) const
{
  return this->find_element(p, allowed_subdomains, relative_tol);
}


// ------------------------------------------------------------
// Explicit Instantiations
template class Tree<2>;
template class Tree<4>;
template class Tree<8>;

} // namespace libMesh
