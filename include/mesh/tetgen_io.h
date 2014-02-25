// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_TETGEN_IO_H
#define LIBMESH_TETGEN_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"

// C++ includes
#include <cstddef>
#include <map>

namespace libMesh
{

// Forward declarations
class MeshBase;
class MeshData;


/**
 * This class implements reading and writing meshes in the TetGen format.
 * Format description:
 * cf. <a href="http://tetgen.berlios.de/">TetGen home page</a>.
 *
 * @author Benjamin S. Kirk, 2004
 */

// ------------------------------------------------------------
// TetGenIO class definition
class TetGenIO : public MeshInput<MeshBase>,
                 public MeshOutput<MeshBase>
{
public:

  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  explicit
  TetGenIO (MeshBase& mesh, MeshData* mesh_data=NULL);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  explicit
  TetGenIO (const MeshBase& mesh, MeshData* mesh_data=NULL);

  /**
   * This method implements reading a mesh from a specified file
   * in TetGen format.
   */
  virtual void read (const std::string& );

  /**
   * This method implements writing a mesh to a specified ".poly" file.
   * ".poly" files defines so called Piecewise Linear Complex (PLC).
   */
  virtual void write (const std::string& );

  /**
   * Data structure to hold node attributes read in from file.
   * What you do with these is up to you!
   */
  std::vector<std::vector<Real> > node_attributes;

  /**
   * Data structure to hold element attributes read in from file.
   * What you do with these is up to you!
   */
  std::vector<std::vector<Real> > element_attributes;

private:


  //-------------------------------------------------------------
  // read support methods

  /**
   * Reads a mesh (nodes & elements) from the file
   * provided through \p node_stream and ele_stream.
   */
  void read_nodes_and_elem (std::istream& node_stream,
                            std::istream& ele_stream);

  /**
   * Method reads nodes from \p node_stream and stores them in
   * vector<Node*> \p nodes in the order they come in.
   * The original node labels are being stored in the
   * map \p _assign_nodes in order to assign the elements to
   * the right nodes later.  In addition, provided it is
   * active, the \p MeshData gets to know the node id from
   * the file, too.
   */
  void node_in (std::istream& node_stream);

  /**
   * Method reads elements and stores them in
   * vector<Elem*> \p elements in the same order as they
   * come in. Within \p TetGenMeshInterface, element labels are
   * ignored, but \p MeshData takes care of such things
   * (if active).
   */
  void element_in (std::istream& ele_stream);

  //-------------------------------------------------------------
  // local data

  /**
   * stores new positions of nodes. Used when reading.
   */
  std::map<dof_id_type,dof_id_type> _assign_nodes;

  /**
   * total number of nodes. Primarily used when reading.
   */
  dof_id_type _num_nodes;

  /**
   * total number of elements. Primarily used when reading.
   */
  dof_id_type _num_elements;

  /**
   * A pointer to the MeshData object you would like to use.
   * with this TetGenIO object.  Can be NULL.
   */
  MeshData* _mesh_data;
};



// ------------------------------------------------------------
// TetGenIO inline members
inline
TetGenIO::TetGenIO (MeshBase& mesh, MeshData* mesh_data) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _mesh_data(mesh_data)
{
}



inline
TetGenIO::TetGenIO (const MeshBase& mesh, MeshData* mesh_data) :
  MeshOutput<MeshBase>(mesh),
  _mesh_data(mesh_data)
{
}


} // namespace libMesh


#endif // LIBMESH_TETGEN_IO_H
