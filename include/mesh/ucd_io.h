// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_UCD_IO_H
#define LIBMESH_UCD_IO_H

// C++ includes

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"
#include "libmesh/boundary_info.h"

namespace libMesh
{

// Forward declarations
class MeshBase;



/**
 * This class implements reading & writing meshes in the AVS's UCD format.
 *
 * @author Benjamin S. Kirk, 2004
 */

// ------------------------------------------------------------
// UCDIO class definition
class UCDIO : public MeshInput<MeshBase>,
	      public MeshOutput<MeshBase>
{
 public:

  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  explicit
  UCDIO (MeshBase&);

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  explicit
  UCDIO (const MeshBase&);

  /**
   * This method implements reading a mesh from a specified file
   * in UCD format.
   */
  virtual void read (const std::string& );

  /**
   * This method implements writing a mesh to a specified file
   * in UCD format.
   */
  virtual void write (const std::string& );

  /**
   * This method implements writing a mesh and solution to a specified file
   * in UCD format. This is internally called by MeshOutput::write_equation_systems
   */
  virtual void write_nodal_data(const std::string& fname, 
				const std::vector<Number>&soln, 
				const std::vector<std::string>& names);
  

 private:

  /**
   * The actual implementation of the read function.
   * The public read interface simply decides which
   * type of stream to pass the implementation.
   */
  void read_implementation (std::istream& in_stream);

  /**
   * The actual implementation of the write function.
   * The public write interface simply decides which
   * type of stream to pass the implementation.
   */
  void write_implementation (std::ostream& out_stream);

  /**
   * Write UCD format header
   */
  void write_header(std::ostream& out, const MeshBase& mesh,
		    dof_id_type n_elems, unsigned int n_vars );

  /**
   * Write node information
   */
  void write_nodes(std::ostream& out, const MeshBase& mesh);

  /**
   * Write element information
   */
  void write_interior_elems(std::ostream& out, const MeshBase& mesh);

  /**
   * Writes all nodal solution variables
   */
  void write_soln(std::ostream& out, const MeshBase& mesh,
		  const std::vector<std::string>& names,
		  const std::vector<Number>&soln);

 };



// ------------------------------------------------------------
// UCDIO inline members
inline
UCDIO::UCDIO (MeshBase& mesh) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh)
{
}



inline
UCDIO::UCDIO (const MeshBase& mesh) :
  MeshOutput<MeshBase> (mesh)
{
}


} // namespace libMesh


#endif // LIBMESH_UCD_IO_H
