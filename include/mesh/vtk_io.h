// The libMesh Finite Element Library.
// Copyright (C) 2002-2015 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_VTK_IO_H
#define LIBMESH_VTK_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"

#ifdef LIBMESH_HAVE_VTK
#include "vtkType.h"
#endif

// C++ includes
#include <cstddef>
#include <map>

// Forward declarations

class vtkUnstructuredGrid;
class vtkPoints;
class vtkCellArray;

namespace libMesh
{

class MeshBase;
class MeshData;

/**
 * This class implements reading and writing meshes in the VTK format.
 * Format description:
 * cf. <a href="http://www.vtk.org/">VTK home page</a>.
 *
 * This class will not have any functionality unless VTK is detected
 * during configure and hence LIBMESH_HAVE_VTK is defined.
 *
 * \author Wout Ruijter
 * \author John W. Peterson
 * \date 2007
 */
class VTKIO : public MeshInput<MeshBase>,
              public MeshOutput<MeshBase>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  explicit
  VTKIO (MeshBase& mesh, MeshData* mesh_data=NULL);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  explicit
  VTKIO (const MeshBase& mesh, MeshData* mesh_data=NULL);

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string&,
                                 const std::vector<Number>&,
                                 const std::vector<std::string>&) libmesh_override;

  /**
   * This method implements reading a mesh from a specified file
   * in VTK format.
   */
  virtual void read (const std::string&) libmesh_override;

  /**
   * Output the mesh without solutions to a .pvtu file
   */
  virtual void write (const std::string&) libmesh_override;

  /**
   * Get a pointer to the VTK datastructure
   */
  vtkUnstructuredGrid* get_vtk_grid();

  /**
   * Setter for compression flag
   */
  void set_compression(bool b);

private:
#ifdef LIBMESH_HAVE_VTK
  /**
   * Map libMesh element types to VTK element types
   */
  vtkIdType get_elem_type(ElemType type);
#endif

  /**
   * write the nodes from the mesh into a vtkUnstructuredGrid
   */
  void nodes_to_vtk();

  /**
   * write the cells from the mesh into a vtkUnstructuredGrid
   */
  void cells_to_vtk();

  /**
   * write the system vectors to vtk
   */
  void system_vectors_to_vtk(const EquationSystems& es, vtkUnstructuredGrid*& grid);

  /**
   * pointer to the VTK grid
   */
  vtkUnstructuredGrid* _vtk_grid;

  /**
   * A pointer to the MeshData object you would like to use.
   * with this VTKIO object.  Can be NULL.
   */
  MeshData* _mesh_data;

  /**
   * Flag to indicate whether the output should be compressed
   */
  bool _compress;

  /**
   * maps global node id to node id of partition
   */
  std::map<dof_id_type, dof_id_type> _local_node_map;
};



} // namespace libMesh


#endif // LIBMESH_VTK_IO_H
