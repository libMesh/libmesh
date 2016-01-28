// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MESH_INPUT_H
#define LIBMESH_MESH_INPUT_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_base.h"

// C++ includes
#include <cstddef>
#include <istream>
#include <string>
#include <vector>

namespace libMesh
{



/**
 * This class defines an abstract interface for \p Mesh input.
 * Specific classes derived from this class actually implement
 * reading various mesh formats.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
template <class MT>
class MeshInput
{
protected:

  /**
   * Default constructor. Will set the _obj to NULL, effectively
   * rendering this object useless.
   */
  explicit
  MeshInput (bool is_parallel_format = false);

  /**
   * Constructor.  Takes a writeable reference to an object.
   * This is the constructor required to read an object.
   */
  explicit
  MeshInput (MT &, const bool is_parallel_format = false);

public:

  /**
   * Destructor.
   */
  virtual ~MeshInput ();

  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string &) = 0;


protected:

  /**
   * Returns the object as a writeable reference.
   */
  MT & mesh ();

  /**
   * Sets the number of partitions in the mesh.  Typically this gets
   * done by the partitioner, but some parallel file formats begin
   * "pre-partitioned".
   */
  void set_n_partitions (unsigned int n_parts) { this->mesh().set_n_partitions() = n_parts; }

  /**
   * A vector of bools describing what dimension elements
   * have been encountered when reading a mesh.
   */
  std::vector<bool> elems_of_dimension;

  /**
   * Reads input from \p in, skipping all the lines
   * that start with the character \p comment_start.
   */
  void skip_comment_lines (std::istream & in,
                           const char comment_start);


private:


  /**
   * A pointer to a non-const object object.
   * This allows us to read the object from file.
   */
  MT * _obj;

  /**
   * Flag specifying whether this format is parallel-capable.
   * If this is false (default) I/O is only permitted when the mesh
   * has been serialized.
   */
  const bool _is_parallel_format;
};



// ------------------------------------------------------------
// MeshInput inline members
template <class MT>
inline
MeshInput<MT>::MeshInput (const bool is_parallel_format) :
  elems_of_dimension(),
  _obj (libmesh_nullptr),
  _is_parallel_format(is_parallel_format)
{
}



template <class MT>
inline
MeshInput<MT>::MeshInput (MT & obj, const bool is_parallel_format) :
  elems_of_dimension(),
  _obj (&obj),
  _is_parallel_format(is_parallel_format)
{
  if (!_is_parallel_format && !this->mesh().is_serial())
    {
      if (this->mesh().processor_id() == 0)
        {
          libmesh_do_once(libMesh::out <<
                          "Warning:  This MeshOutput subclass only supports meshes which have been serialized!"
                          << std::endl;);
        }
    }
}



template <class MT>
inline
MeshInput<MT>::~MeshInput ()
{
}



template <class MT>
inline
MT & MeshInput<MT>::mesh ()
{
  if (_obj == libmesh_nullptr)
    libmesh_error_msg("ERROR: _obj should not be NULL!");
  return *_obj;
}



template <class MT>
void MeshInput<MT>::skip_comment_lines (std::istream & in,
                                        const char comment_start)
{
  char c, line[256];

  while (in.get(c), c==comment_start)
    in.getline (line, 255);

  // put back first character of
  // first non-comment line
  in.putback (c);
}


} // namespace libMesh


#endif // LIBMESH_MESH_INPUT_H
