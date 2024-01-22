// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_config.h"


// C++ includes
#include <sstream>

// Local includes
#include "libmesh/mesh_tet_interface.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/unstructured_mesh.h"

namespace libMesh
{

//----------------------------------------------------------------------
// MeshTetInterface class members
MeshTetInterface::MeshTetInterface (UnstructuredMesh & mesh) :
  _elem_type(TET4), _mesh(mesh)
{
}



unsigned MeshTetInterface::check_hull_integrity()
{
  // Check for easy return: if the Mesh is empty (i.e. if
  // somebody called triangulate_conformingDelaunayMesh on
  // a Mesh with no elements, then hull integrity check must
  // fail...
  if (_mesh.n_elem() == 0)
    return 3;

  for (auto & elem : this->_mesh.element_ptr_range())
    {
      // Check for proper element type
      if (elem->type() != TRI3)
        {
          //libmesh_error_msg("ERROR: Some of the elements in the original mesh were not TRI3!");
          return 1;
        }

      for (auto neigh : elem->neighbor_ptr_range())
        {
          if (neigh == nullptr)
            {
              // libmesh_error_msg("ERROR: Non-convex hull, cannot be tetrahedralized.");
              return 2;
            }
        }
    }

  // If we made it here, return success!
  return 0;
}



void MeshTetInterface::process_hull_integrity_result(unsigned result)
{
  if (result != 0)
    {
      libMesh::err << "Error! Conforming Delaunay mesh tetrahedralization requires a convex hull." << std::endl;

      if (result==1)
        {
          libMesh::err << "Non-TRI3 elements were found in the input Mesh.  ";
          libMesh::err << "A constrained Delaunay triangulation requires a convex hull of TRI3 elements." << std::endl;
        }

      libmesh_error_msg("Consider calling TetGenMeshInterface::pointset_convexhull() followed by Mesh::find_neighbors() first.");
    }
}



void MeshTetInterface::delete_2D_hull_elements()
{
  for (auto & elem : this->_mesh.element_ptr_range())
    {
      // Check for proper element type. Yes, we legally delete elements while
      // iterating over them because no entries from the underlying container
      // are actually erased.
      if (elem->type() == TRI3)
        _mesh.delete_elem(elem);
    }

  // We just removed any boundary info associated with hull element
  // edges, so let's update the boundary id caches.
  this->_mesh.get_boundary_info().regenerate_id_sets();
}



} // namespace libMesh
