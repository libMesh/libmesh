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



#ifndef LIBMESH_UNV_IO_H
#define LIBMESH_UNV_IO_H


// Local includes
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"

// C++ inludes
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace libMesh
{

// Forward declarations
class MeshBase;

/**
 * The \p UNVIO class implements the Ideas \p UNV universal
 * file format.  This class enables both reading and writing
 * \p UNV files.
 *
 * Author history
 * \author Tammo Kaschner
 * \author Daniel Dreyer
 * \author Benjamin S. Kirk
 * \author John W. Peterson
 * \date 2003, 2004, 2014
 */
class UNVIO : public MeshInput<MeshBase>,
              public MeshOutput<MeshBase>
{

public:

  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  UNVIO (MeshBase & mesh);

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  UNVIO (const MeshBase & mesh);

  /**
   * Destructor.
   */
  virtual ~UNVIO ();

  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string &) libmesh_override;

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string &) libmesh_override;

  /**
   * Set the flag indicationg if we should be verbose.
   */
  bool & verbose ();

  /**
   * Read a UNV data file containing a dataset of type "2414".
   * For more info, see http://tinyurl.com/htcf6zm
   */
  void read_dataset(std::string file_name);

  /**
   * @returns a pointer the values associated with the node \p node,
   * as read in by the read_dataset() method.  If no values exist for
   * the node in question, a libmesh_nullptr is returned instead.  It
   * is up to the user to check the return value before using it.
   */
  const std::vector<Number> * get_data (Node * node) const;

private:


  /**
   * The actual implementation of the read function.
   * The public read interface simply decides which
   * type of stream to pass the implementation.
   */
  void read_implementation (std::istream & in_stream);

  /**
   * The actual implementation of the write function.
   * The public write interface simply decides which
   * type of stream to pass the implementation.
   */
  void write_implementation (std::ostream & out_stream);

  //-------------------------------------------------------------
  // read support methods

  /**
   * Read nodes from file.
   */
  void nodes_in (std::istream & in_file);

  /**
   * Method reads elements and stores them in
   * \p std::vector<Elem *> \p _elements in the same order as they
   * come in. Within \p UNVIO, element labels are
   * ignored.
   */
  void elements_in (std::istream & in_file);

  /**
   * Reads the "groups" section of the file. The format of the groups section is described here:
   * http://www.sdrl.uc.edu/universal-file-formats-for-modal-analysis-testing-1/file-format-storehouse/unv_2467.htm
   */
  void groups_in(std::istream & in_file);

  //-------------------------------------------------------------
  // write support methods

  /**
   * Outputs nodes to the file \p out_file.  Do not use this directly,
   * but through the proper write method.
   */
  void nodes_out (std::ostream & out_file);

  /**
   * Outputs the element data to the file \p out_file. Do not use this
   * directly, but through the proper write method.
   */
  void elements_out (std::ostream & out_file);

  /**
   * Returns the maximum geometric element dimension encountered while
   * reading the Mesh.  Only valid after the elements have been read
   * in and the elems_of_dimension array has been populated.
   */
  unsigned char max_elem_dimension_seen ();

  /**
   * Replaces "1.1111D+00" with "1.1111e+00" if necessary.  Returns
   * true if the replacement occurs, false otherwise.  This function
   * only needs to be called once per stream, one can assume that if
   * one number needs rewriting, they all do.
   */
  bool need_D_to_e (std::string & number);

  //-------------------------------------------------------------
  // local data

  /**
   * should be be verbose?
   */
  bool _verbose;

  /**
   * Maps UNV node IDs to libMesh Node*s. Used when reading. Even if the
   * libMesh Mesh is renumbered, this map should continue to be valid.
   */
  std::map<dof_id_type, Node *> _unv_node_id_to_libmesh_node_ptr;

  /**
   * label for the node dataset
   */
  static const std::string _nodes_dataset_label;

  /**
   * label for the element dataset
   */
  static const std::string _elements_dataset_label;

  /**
   * label for the groups dataset
   */
  static const std::string _groups_dataset_label;

  /**
   * Map UNV element IDs to libmesh element IDs.
   */
  std::map<unsigned, unsigned> _unv_elem_id_to_libmesh_elem_id;

  /**
   * Map from libMesh Node* to data at that node, as read in by the
   * read_dataset() function.
   */
  std::map<Node *, std::vector<Number> > _node_data;
};



} // namespace libMesh


#endif // LIBMESH_UNV_IO_H
