// $Id: exodusII_io.h,v 1.2 2004-11-17 07:52:17 benkirk Exp $

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



#ifndef __exodusII_io_h__
#define __exodusiI_io_h__


// C++ inludes

// Local includes
#include "mesh_input.h"

// Forward declarations
class MeshBase;


/**
 * The \p ExodusII_IO class implements reading meshes in the
 * \p ExodusII file format from Sandia National Labs.  Due to
 * licensing restrictions the external libraries required to
 * read such files cannot be redistributed with \p libMesh.
 * However, if you have these libraries, this class will allow
 * you to use the API to read ExodusII, or \p .exd, files.
 *
 * @author Benjamin Kirk, John Peterson, 2004.
 */

// ------------------------------------------------------------
// ExodusII_IO class definition
class ExodusII_IO : public MeshInput<MeshBase>
{

 public:
  
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  ExodusII_IO (MeshBase& mesh);
  
  /**
   * Destructor.
   */
  virtual ~ExodusII_IO ();
  
  /**
   * This method implements reading a mesh from a specified file.
   * Open the file named \p name and read the mesh in Sandia National Lab's
   * ExodusII format. This is the method to use for reading in meshes generated
   * by cubit.  Works in 2D for \p TRIs, \p TRI6s, \p QUAD s, and \p QUAD9s.
   * Works in 3D for \p TET4s, \p TET10s, \p HEX8s, and \p HEX27s.
   */
  virtual void read (const std::string& name);

  /**
   * Set the flag indicationg if we should be verbose.
   */
  bool & verbose ();

  
 private:
  

  //-------------------------------------------------------------
  // local data

  /**
   * should be be verbose?
   */
  bool _verbose;
};



// ------------------------------------------------------------
// MeshIO inline members
inline
ExodusII_IO::ExodusII_IO (MeshBase& mesh) :
  MeshInput<MeshBase> (mesh),
  _verbose (false)
{
}



inline
ExodusII_IO::~ExodusII_IO ()
{
}



inline
bool & ExodusII_IO::verbose ()
{
  return _verbose;
}



#endif // #define __exodusII_io_h__
