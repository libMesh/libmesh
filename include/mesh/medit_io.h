// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MEDIT_IO_H
#define LIBMESH_MEDIT_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_output.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward declarations:
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;

/**
 * This class implements writing meshes in the mesh format
 * used by the MEdit visualization tool developed in the Gamma Project
 * at INRIA Roquencourt.
 * For a full description of the mesh format and to obtain the
 * MEdit software see the
 * <a href="http://www-rocq1.inria.fr/gamma/medit/medit.html">MEdit home page</a>.
 *
 * \author Florian Prill
 * \date 2004
 */
class MEDITIO : public MeshOutput<MeshBase>
{
public:

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  explicit
  MEDITIO (const MeshBase &);

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * and the desired scalar index for mesh colouring.
   * MEdit seems to understand only one scalar value.
   */
  MEDITIO (const MeshBase &, unsigned int c);

  /**
   * This method implements writing a mesh to a specified ".mesh" file.
   */
  virtual void write (const std::string &) override;

  /**
   * Bring in base class functionality for name resolution and to
   * avoid warnings about hidden overloaded virtual functions.
   */
  using MeshOutput<MeshBase>::write_nodal_data;

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string &,
                                 const std::vector<Number> &,
                                 const std::vector<std::string> &) override;

  /**
   * Flag indicating whether or not to write a binary file
   */
  bool & binary ();

private:

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write an ASCII file.
   */
  void write_ascii (const std::string &,
                    const std::vector<Number> * = nullptr,
                    const std::vector<std::string> * = nullptr);

  /**
   * Flag to write binary data.
   */
  bool _binary;

  unsigned int scalar_idx;
};



// ------------------------------------------------------------
// medit_io inline members
inline
MEDITIO::MEDITIO (const MeshBase & mesh_in) :
  MeshOutput<MeshBase> (mesh_in),
  _binary (false),
  scalar_idx (0)
{
}

inline
MEDITIO::MEDITIO (const MeshBase & mesh_in, unsigned int c) :
  MeshOutput<MeshBase> (mesh_in),
  _binary    (false),
  scalar_idx (c)
{
}


inline
bool & MEDITIO::binary ()
{
  return _binary;
}


} // namespace libMesh


#endif // LIBMESH_MEDIT_IO_H
