// $Id: mesh_tetgen_support.h,v 1.1 2003-12-08 19:14:41 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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


#ifndef __mesh_tetgen_support_h__
#define __mesh_tetgen_support_h__

// C++ includes
#include <vector>
#include <map>



// Local includes
#include "mesh_base.h"



// Forward Declarations
class MeshData;


/**
 * Class \p TetGenMeshInterface provides an interface
 * for reading a mesh from a file in TetGen file format.
 */

class TetGenMeshInterface
{
public:

  /**
   * Constructor.  Takes the data relevant 
   * for reading mesh or reading data.
   */
  TetGenMeshInterface(std::vector<Node*>& nodes,
		      std::vector<Elem*>& elements,
		      MeshData& md);

  /**
   * Reads a mesh (nodes & elements) from the file
   * provided through \p node_stream and ele_stream.
   */
  void read (std::istream& node_stream, std::istream& ele_stream);

  /**
   * Destructor.
   */
  ~TetGenMeshInterface();

protected:


  //-------------------------------------------------------------
  // read support methods

  /**
   * Method reads nodes from \p node_stream and stores them in
   * vector<Node*> \p nodes in the order they come in.
   * The original node labels are being stored in the
   * map \p _assign_nodes in order to assign the elements to
   * the right nodes later.  In addition, provided it is 
   * active, the \p MeshData gets to know the node id from
   * the file, too.
   */
  void node_in(std::istream& node_stream);

  /**
   * Method reads elements and stores them in
   * vector<Elem*> \p elements in the same order as they
   * come in. Within \p TetGenMeshInterface, element labels are
   * ignored, but \p MeshData takes care of such things
   * (if active).
   */
  void element_in(std::istream& ele_stream);

  //-------------------------------------------------------------
  // local data
  
  /**
   * vector holding the nodes.
   */
  std::vector<Node*>& _nodes;    

  /**
   * vector holding the elements.
   */
  std::vector<Elem*>& _elements; 

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

  /**
   * writable reference to the class that
   * handles foreign node/element ids
   */
  MeshData& _mesh_data;

};

#endif





