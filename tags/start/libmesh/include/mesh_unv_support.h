// $Id: mesh_unv_support.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>


// Local includes
#include "mesh.h"
#include "elem.h"
#include "point.h"

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
   * read and the vertices and elements will be stored
   * in the vectors specified in the arguments
   */
  UnvInterface(std::istream& _in,
	       std::vector<Point>& _vertices,
	       std::vector<Elem*>& _elements);
  
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
   * Method for reading nodes (dataset "2411").
   */
  void node_in();

  /**
   * Method for reading elements (dataset "2412").
   */
  void element_in();


  /**
   * Method for converting exponential notation
   * from "D" to "e".
   */
  std::string& D_to_e(std::string& number);


  // References
  std::istream& phys_file;                // input stream, physical file
  std::vector<Point>& vertices;           // vector holding the nodes
  std::vector<Elem*>& elements;           // vector holding the elements

  // stores new positions of nodes
  std::map<unsigned int,unsigned int> assign_nodes; 

  // stores positions of datasets in the stream
  std::map<std::string,std::streampos> ds_position;

  // store the total number of nodes and elements
  unsigned int num_nodes,
               num_elements;


private:

  // labels of important datasets
  std::string label_dataset_nodes,                  // Nodes
              label_dataset_elms;                   // Elements

    
  // temporary buffer
  std::fstream temporary_file;
  char* temporary_file_name;


  // whether we need to convert notation of exponentials
  bool need_D_to_e;
};

#endif





