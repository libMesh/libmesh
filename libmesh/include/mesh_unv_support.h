// $Id: mesh_unv_support.h,v 1.12 2003-07-12 16:56:39 ddreyer Exp $

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


#ifndef __mesh_unv_support_h__
#define __mesh_unv_support_h__

// C++ includes
#include <fstream>
#include <vector>
#include <set>



// Local includes
#include "mesh_base.h"



// Forward Declarations
class MeshData;


/**
 * Class \p UnvMeshInterface provides an interface
 * for reading a mesh (datasets 2411 and 2412)
 * from a file in I-deas Universal file format.
 *
 * @author: Tammo Kaschner
 */
class UnvMeshInterface
{
public:

  /**
   * Constructor.  Takes the data relevant 
   * for reading/writing mesh or reading data.
   * Note that for simplicity, the node and element
   * vectors have to be readable even when only
   * read access is needed.
   */
  UnvMeshInterface(std::vector<Node*>& nodes,
		   std::vector<Elem*>& elements,
		   MeshData& md);

  /**
   * Reads a mesh (nodes & elements) from the file
   * provided through \p in_stream, gives progress 
   * reports if desired using \p verbose.
   */
  void read (std::istream& in_stream,
	     const bool verbose = false);

  /**
   * Writes a mesh (nodes & elements) to the file
   * provided through \p out_stream.
   */
  void write (std::ostream& out_stream);
  
  /**
   * Destructor.
   */
  ~UnvMeshInterface();

protected:


  //-------------------------------------------------------------
  // read support methods
  /**
   * When reading, the node related data is
   * buffered in a temporary file prior to the actual
   * read process.  This helps e.g. in counting nodes
   * beforehand for pre-allocation.
   */
  void buffer_nodes (std::istream& physical_file,
		     std::fstream& temp_file);
  
  /**
   * Method reads nodes from \p tempfile and stores them in
   * vector<Node*> \p nodes in the order they come in.
   * The original node labels are being stored in the
   * map \p _assign_nodes in order to assign the elements to
   * the right nodes later.  In addition, provided it is 
   * active, the \p MeshData gets to know the node id from
   * the Universal file, too.
   */
  void node_in(std::fstream& temp_file);

  /**
   * When reading, the element related data is
   * buffered in a temporary file prior to the actual
   * read process.  This helps e.g. in counting elements
   * beforehand for pre-allocation.
   */
  void buffer_elements (std::istream& physical_file,
			std::fstream& temp_file);

  /**
   * Method reads elements and stores them in
   * vector<Elem*> \p elements in the same order as they
   * come in. Within \p UnvMeshInterface, element labels are
   * ignored, but \p MeshData takes care of such things
   * (if active).
   */
  void element_in(std::fstream& temp_file);


  //-------------------------------------------------------------
  // write support methods
  /**
   * Outputs nodes to the file \p out_file.
   * For this to work, the \p MeshData of the current
   * \p Mesh has to be active.  Do not use this directly,
   * but through the proper write method.
   */
  void node_out(std::ostream& out_file);

  /**
   * Outputs the element data to the file \p out_file.
   * For this to work, the \p MeshData of the current
   * \p Mesh has to be active. Do not use this directly,
   * but through the proper write method.
   */
  void element_out(std::ostream& out_file);


  //-------------------------------------------------------------
  // misc helpers
  /**
   * Method for setting the stream pointer
   * to the position of the dataset with label
   * \p ds_num in the temporary_file \p temp_file.
   */
  void set_stream_pointer(std::fstream& temp_file,
			  const std::string& ds_num);

  /**
   * Method for converting exponential notation
   * from "D" to "e", for example
   * \p 3.141592654D+00 \p --> \p 3.141592654e+00
   * in order to make it readable for C++.
   */
  std::string& D_to_e(std::string& number);


  //-------------------------------------------------------------
  // local data
  /**
   * vector holding the nodes.  Either used during writing
   * or reading.
   */
  std::vector<Node*>& _nodes;    

  /**
   * vector holding the elements.  Either used during writing
   * or reading.
   */
  std::vector<Elem*>& _elements; 

  /**
   * stores new positions of nodes.  Used when reading.
   */
  std::map<unsigned int,unsigned int> _assign_nodes; 

  /**
   * stores positions of datasets in the stream.  Used when reading.
   */
  std::map<std::string,std::streampos> _ds_position;

  /**
   * total number of nodes.  Primarily used when reading.
   */
  unsigned int _num_nodes;

  /**
   * total number of elements.  Primarily used when reading.
   */
  unsigned int _num_elements;

  /**
   * writable reference to the class that
   * handles foreign node/element ids
   */
  MeshData& _mesh_data;

  /**
   * labels for the node dataset
   */
  static const std::string _label_dataset_nodes;

  /**
   * labels for the element dataset
   */
  static const std::string _label_dataset_elements;

  /**
   * whether we need to convert notation of exponentials.
   * Used when reading.
   */
  bool _need_D_to_e;

};

#endif





