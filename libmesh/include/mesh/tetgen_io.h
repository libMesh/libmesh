// $Id: tetgen_io.h,v 1.1 2004-03-19 19:16:52 benkirk Exp $

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



#ifndef __tetgen_io_h__
#define __tetgen_io_h__

// Local includes
#include "libmesh_common.h"
#include "mesh_io.h"


/**
 * This class implements reading and writing meshes in the TetGen format.
 *
 * @author Benjamin S. Kirk, 2004
 */

// ------------------------------------------------------------
// TetGenIO class definition
class TetGenIO : public MeshIO
{
 public:
  
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  TetGenIO (Mesh&);
  
  /**
   * This method implements reading a mesh from a specified file
   * in TetGen format.
   */
  virtual void read (const std::string& );

  
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
  std::map<unsigned int,unsigned int> _assign_nodes; 

  /**
   * total number of nodes. Primarily used when reading.
   */
  unsigned int _num_nodes;

  /**
   * total number of elements. Primarily used when reading.
   */
  unsigned int _num_elements;
 
};



// ------------------------------------------------------------
// TetGenIO inline members
inline
TetGenIO::TetGenIO (Mesh& mesh) :
  MeshIO (mesh)
{
}



#endif // #define __tetgen_io_h__
