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

#ifndef LIBMESH_OFF_IO_H
#define LIBMESH_OFF_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"

namespace libMesh
{

// Forward declarations
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;

/**
 * This class is responsible for reading an unstructured, triangulated
 * surface in the standard OFF OOGL format.
 *
 * \author John W. Peterson
 * \date 2004
 * \brief Reads OOF OOGL triangulated surface files.
 */
class OFFIO : public MeshInput<MeshBase>
{
public:
  /**
   *  Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements.
   */
  explicit
  OFFIO (MeshBase &);

  /**
   * Reads in an OFF OOGL data file based on the string
   * you pass it.
   */
  virtual void read (const std::string & name) override;

private:
  /**
   * Implementation of the read() function.  This function
   * is called by the public interface function and implements
   * reading the file.
   */
  void read_stream (std::istream & in);
};


// ------------------------------------------------------------
// OFFIO inline members
inline
OFFIO::OFFIO (MeshBase & mesh_in) :
  MeshInput<MeshBase> (mesh_in)
{}


} // namespace libMesh


#endif // LIBMESH_OFF_IO_H
