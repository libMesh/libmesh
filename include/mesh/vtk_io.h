// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "vtkSmartPointer.h"
#endif

// C++ includes
#include <cstddef>
#include <map>

// Forward declarations
class vtkUnstructuredGrid;

namespace libMesh
{

class MeshBase;

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
  VTKIO (MeshBase & mesh);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  explicit
  VTKIO (const MeshBase & mesh);

  /**
   * Bring in base class functionality for name resolution and to
   * avoid warnings about hidden overloaded virtual functions.
   */
  using MeshOutput<MeshBase>::write_nodal_data;

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   *
   * As with other Mesh IO classes, this interface is still still
   * "available" when !LIBMESH_HAVE_VTK, however, it will throw a
   * runtime error.
   */
  virtual void write_nodal_data (const std::string &,
                                 const std::vector<Number> &,
                                 const std::vector<std::string> &) libmesh_override;

  /**
   * This method implements reading a mesh from a specified file
   * in VTK format.
   *
   * As with other Mesh IO classes, this interface is still still
   * "available" when !LIBMESH_HAVE_VTK, however, it will throw a
   * runtime error.
   */
  virtual void read (const std::string &) libmesh_override;

  /**
   * Output the mesh without solutions to a .pvtu file.
   *
   * As with other Mesh IO classes, this interface is still still
   * "available" when !LIBMESH_HAVE_VTK, however, it will throw a
   * runtime error.
   */
  virtual void write (const std::string &) libmesh_override;

#ifdef LIBMESH_HAVE_VTK

  /**
   * Setter for compression flag
   */
  void set_compression(bool b);

  /**
   * Get a pointer to the VTK unstructed grid data structure.
   */
  vtkUnstructuredGrid * get_vtk_grid();

private:
  /**
   * write the nodes from the mesh into a vtkUnstructuredGrid
   */
  void nodes_to_vtk();

  /**
   * write the cells from the mesh into a vtkUnstructuredGrid
   */
  void cells_to_vtk();

  /**
   * Write the system vectors to vtk
   *
   * This function is not currently used by anything, so it is commented
   * out, and may eventually be removed entirely.
   */
  // void system_vectors_to_vtk(const EquationSystems & es,
  //                            vtkUnstructuredGrid * & grid);

  /**
   * pointer to the VTK grid. the vtkSmartPointer will automatically
   * initialize the value to null and keep track of reference
   * counting.
   */
  vtkSmartPointer<vtkUnstructuredGrid> _vtk_grid;

  /**
   * Flag to indicate whether the output should be compressed
   */
  bool _compress;

  /**
   * maps global node id to node id of partition
   */
  std::map<dof_id_type, dof_id_type> _local_node_map;


  /**
   * Helper object that holds a map from VTK to libMesh element types
   * and vice-versa.
   */
  struct ElementMaps
  {
    // Associate libmesh_type with vtk_type (searchable in both directions).
    void associate(ElemType libmesh_type, vtkIdType vtk_type)
    {
      writing_map[libmesh_type] = vtk_type;
      reading_map[vtk_type] = libmesh_type;
    }

    // Find an entry in the writing map, or throw an error.
    vtkIdType find(ElemType libmesh_type)
    {
      std::map<ElemType, vtkIdType>::iterator it = writing_map.find(libmesh_type);

      if (it == writing_map.end())
        libmesh_error_msg("Element type " << libmesh_type << " not available in VTK.");

      return it->second;
    }

    // Find an entry in the reading map, or throw an error.
    ElemType find(vtkIdType vtk_type)
    {
      std::map<vtkIdType, ElemType>::iterator it = reading_map.find(vtk_type);

      if (it == reading_map.end())
        libmesh_error_msg("Element type " << vtk_type << " not available in libMesh.");

      return it->second;
    }

    std::map<ElemType, vtkIdType> writing_map;
    std::map<vtkIdType, ElemType> reading_map;
  };

  /**
   * ElementMaps object that is built statically and used by
   * all instances of this class.
   */
  static ElementMaps _element_maps;

  /**
   * Static function used to construct the _element_maps struct.
   */
  static ElementMaps build_element_maps();

#endif
};



} // namespace libMesh


#endif // LIBMESH_VTK_IO_H
