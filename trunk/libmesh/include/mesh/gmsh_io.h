// $Id: gmsh_io.h,v 1.1 2004-07-13 21:48:41 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __gmsh_io_h__
#define __gmsh_io_h__

// Local includes
#include "libmesh_common.h"
#include "mesh_io.h"

// Forward declarations
class MeshBase;



/**
 * This class implements writing meshes in the Gmsh format.
 * For a full description of the Gmsh format and to obtain the
 * GMSH software see
 * <a href="http://http://www.geuz.org/gmsh/">the Gmsh home page</a>
 *
 * @author John W. Peterson, 2004
 */

// ------------------------------------------------------------
// GMVIO class definition
class GmshIO : public MeshIO<MeshBase>
{
 public:

  /**
   * Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements via the read() command.
   */
  GmshIO (MeshBase& mesh);
  
  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  GmshIO (const MeshBase& mesh);

  /**
   * Reads in a Gmsh data file based on the string
   * you pass it.
   */
  virtual void read (const std::string& name);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string& name);
  

private:
  /**
   * Implementation of the read() function.  This function
   * is called by the public interface function and implements
   * reading the file.
   */
  virtual void read_stream (std::istream& in);

  /**
   * This method implements writing a mesh to a
   * specified file.  This will write an ASCII file.
   */
  virtual void write_stream (std::ostream& out);

};



// ------------------------------------------------------------
// GmshIO inline members
inline
GmshIO::GmshIO (const MeshBase& mesh) :
  MeshIO<MeshBase> (mesh)
{
}


inline
GmshIO::GmshIO (MeshBase& mesh) :
  MeshIO<MeshBase>  (mesh)
{}





#endif // #define __gmsh_io_h__
