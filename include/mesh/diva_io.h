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


#ifndef LIBMESH_DIVA_IO_H
#define LIBMESH_DIVA_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_output.h"

// C++ Includes   -----------------------------------

namespace libMesh
{

// Forward declarations
class MeshBase;

/**
 * This class implements writing meshes in the Diva format.
 * This is a scientific visualization program created by Kelly Gaither.
 * More information on Diva can be found here:
 * http://www.erc.msstate.edu/simcenter/docs/diva/
 *
 * \author John W. Peterson
 * \date 2004
 */
class DivaIO : public MeshOutput<MeshBase>
{
public:
  /**
   * Note that only writing diva files is supported since Diva is
   * not a mesh generator.
   */
  //TODO:[JWP] figure out how to implement a const write method!
  explicit
  DivaIO (const MeshBase &);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string &) libmesh_override;

private:

  /**
   * The actual implementation of writing the diva file.  This
   * file is called by the public interface file after it
   * constructs an ofstream.
   */
  void write_stream(std::ostream & out);
};




// ------------------------------------------------------------
// DivaIO inline members
inline
DivaIO::DivaIO (const MeshBase & mesh_in) :
  MeshOutput<MeshBase>  (mesh_in)
{}


} // namespace libMesh

#endif // LIBMESH_DIVA_IO_H
