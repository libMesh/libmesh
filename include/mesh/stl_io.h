// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_STL_IO_H
#define LIBMESH_STL_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_output.h"
#include "libmesh/mesh_input.h"

// C++ includes
#include <map>

namespace libMesh
{

// Forward declarations
class MeshBase;
enum ElemType : int;

/**
 * This class implements reading and writing triangle meshes in the
 * STL format.
 *
 * \author Roy H. Stogner
 * \date 2024
 */
class STLIO : public MeshInput<MeshBase>,
              public MeshOutput<MeshBase>
{
public:

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  explicit
  STLIO (const MeshBase &);

  /**
   * Constructor.  Takes a writable reference to a mesh object.
   * This constructor is required to let us read in a mesh.
   */
  explicit
  STLIO (MeshBase &);

  /**
   * This method implements writing a mesh to a specified file.
   * We currently only support ASCII writes.
   */
  virtual void write (const std::string &) override;

  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string & mesh_file) override;

  /**
   * This method implements reading a mesh from a specified ASCII
   * input stream.
   */
  virtual void read_ascii (std::istream & input);

  /**
   * This method implements reading a mesh from a specified binary
   * input stream.
   *
   * If the size in bytes is known a priori then it can be passed in
   * to allow for additional error checking.
   */
  virtual void read_binary (std::istream & input,
                            std::size_t input_size = 0);

  /**
   * This method gets a name after a set or an ASCII read
   */
  const std::string & name () { return _name; }

  /**
   * This method sets a name to write
   */
  virtual void set_name (const std::string & name) { _name = name; }

  /**
   * Flag indicating whether or not to subdivide second order
   * elements when writing.  This option is not yet supported.
   */
  bool subdivide_second_order() { return _subdivide_second_order; }

  void set_subdivide_second_order(bool subdivide)
  { _subdivide_second_order = subdivide; }

private:

  /**
   * Helper to open possibly-zipped files
   */
  std::unique_ptr<std::istream> open_file(const std::string & filename);

  /**
   * Flag to subdivide second order elements
   */
  bool _subdivide_second_order;

  std::string _name;
};

} // namespace libMesh


#endif // LIBMESH_STL_IO_H
