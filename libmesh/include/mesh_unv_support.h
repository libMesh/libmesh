// $Id: mesh_unv_support.h,v 1.7 2003-05-14 11:54:36 ddreyer Exp $

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

// Local includes
#include "mesh.h"



/**
 * Class UnvInterface provides an Interface to
 * a stream object that contains a mesh in
 * I-deas UNV format.
 *
 * @author: Tammo Kaschner
 */
class UnvInterface
{
public:

  /**
   * Constructor.
   * When you create an interface the file will be
   * read and the nodes and elements will be stored
   * in the vectors specified in the arguments
   */
  UnvInterface(std::istream& _in,
	       std::vector<Node*>&  _nodes,
	       std::vector<Elem*>& _elements,
	       BoundaryInfo& boundary_info);
  
  /**
   * Destructor.
   */
  virtual ~UnvInterface();

protected:
  
  /**
   * Scans the file, locates datasets and stores
   * their positions.
   */
  void init();

  /**
   * Scans a dataset and extracts important
   * information needed for further processing.
   * Currently recognizes nodes and elements.
   */
  void scan_dataset(std::string ds_num);

  /**
   * Method for setting the stream pointer
   * to the position of a dataset in temporary_file.
   * Argument is the dataset number, which must
   * be written as a string such as "2411".
   */
  void set_stream_pointer(std::string ds_num);
  
  /**
   * Method reads nodes from the virtual file and stores them in
   * vector<Node> nodes in the order they come in.
   * The original node labels are being stored in the
   * map _assign_nodes in order to assign the elements to
   * the right nodes later.  In addition, if there is
   * a \p BoundaryData, this \p BoundaryData gets to know 
   * the node id in the Universal file, too.
   */
  void node_in();

  /**
   * Method reads elements and stores them in
   * vector<Elem*> elements in the same order as they
   * come in. Within UnvInterface, element labels are
   * ignored, but \p BoundaryData takes care of such things.
   *
   * Obviously, nodes are re-numbered.
   */
  void element_in();

  /**
   * Method for converting exponential notation
   * from "D" to "e".
   */
  std::string& D_to_e(std::string& number);


  /**
   * input stream, physical file
   */
  std::istream& _phys_file;     

  /**
   * vector holding the nodes
   */
  std::vector<Node*>& _nodes;    

  /**
   * vector holding the elements
   */
  std::vector<Elem*>& _elements; 

  /**
   * stores new positions of nodes
   */
  std::map<unsigned int,unsigned int> _assign_nodes; 

  /**
   * stores positions of datasets in the stream
   */
  std::map<std::string,std::streampos> _ds_position;

  /**
   * total number of nodes
   */
  unsigned int _num_nodes;

  /**
   * total number of elements
   */
  unsigned int _num_elements;


private:

  /**
   * writable reference to the BoundaryInfo
   */
  BoundaryInfo& _boundary_info;

  /**
   * labels for the node dataset
   */
  std::string _label_dataset_nodes;

  /**
   * labels for the element dataset
   */
  std::string _label_dataset_elms;
    
  /**
   * temporary file, simplifies node & element conversion
   */
  std::fstream _temporary_file;

  /**
   * Name of temporary file
   */
  char* _temporary_file_name;

  /**
   * whether we need to convert notation of exponentials
   */
  bool _need_D_to_e;

};

#endif





