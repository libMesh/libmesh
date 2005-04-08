// $Id: gmsh_io.h,v 1.4 2005-04-08 20:39:24 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_input.h"
#include "mesh_output.h"

// Forward declarations
class MeshBase;



/**
 * This class implements writing meshes in the Gmsh format.
 * For a full description of the Gmsh format and to obtain the
 * GMSH software see
 * <a href="http://http://www.geuz.org/gmsh/">the Gmsh home page</a>
 *
 * @author John W. Peterson, 2004
 * @author Martin Lüthi, 2005 (support for reading meshes and writing results)
 */

// ------------------------------------------------------------
// GMVIO class definition
class GmshIO : public MeshInput<MeshBase>,
               public MeshOutput<MeshBase>
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
   * Reads in a mesh in the Gmsh *.msh format
   * from the ASCII file given by name.
   */
  virtual void read (const std::string& name);

  /**
   * This method implements writing a mesh to a specified file
   * in the Gmsh *.msh format.
   */
  virtual void write (const std::string& name);
  
  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string&,
				 const std::vector<Number>&,
				 const std::vector<std::string>&);

  /**
   * Flag indicating whether or not to write a binary file.  While binary
   * files may end up being smaller than equivalent ASCII files, they will
   * almost certainly take longer to write.  The reason for this is that
   * the ostream::write() function which is used to write "binary" data to
   * streams, only takes a pointer to char as its first argument.  This means
   * if you want to write anything other than a buffer of chars, you first
   * have to use a strange memcpy hack to get the data into the desired format.
   * See the templated to_binary_stream() function below.
   */
  bool & binary ();


private:
  /**
   * Implementation of the read() function.  This function
   * is called by the public interface function and implements
   * reading the file.
   */
  virtual void read_mesh (std::istream& in);

  /**
   * This method implements writing a mesh to a
   * specified file.  This will write an ASCII *.msh file.
   */
  virtual void write_mesh (std::ostream& out);

  /**
   * This method implements writing a mesh with nodal data to a specified file
   * where the nodal data and variable names are optionally provided.  This
   * will write an ASCII or binary *.pos file, depending on the binary flag.
   */
  void write_post (const std::string&,
                   const std::vector<Number>* = NULL,
                   const std::vector<std::string>* = NULL);

  /**
   * Flag to write binary data.
   */
  bool _binary;
 


};



// ------------------------------------------------------------
// GmshIO inline members
inline
GmshIO::GmshIO (const MeshBase& mesh) :
  MeshOutput<MeshBase> (mesh),
  _binary        (false)
{
}


inline
GmshIO::GmshIO (MeshBase& mesh) :
  MeshInput<MeshBase>  (mesh),
  MeshOutput<MeshBase> (mesh),
  _binary (false)
{}

inline
bool & GmshIO::binary ()
{
  return _binary;
}


#endif // #define __gmsh_io_h__
