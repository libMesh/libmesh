// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __postscript_io_h__
#define __postscript_io_h__

// Local includes
#include "libmesh_common.h"
#include "mesh_output.h"

// Forward declarations
class MeshBase;


/**
 * This class implements writing 2D meshes in Postscript.  It borrows
 * several ideas from, and is a more simple-minded version of, the
 * DataOutBase::write_eps() function from Deal II.  Only output is
 * supported here, and only the Mesh (none of the data) is written.
 * The main use I imagined for this class is creating nice Mesh
 * images for publications, since I didn't find/don't know of a free
 * visualization program which would do this.
 *
 * @author John W. Peterson, 2008
 */
class PostscriptIO : public MeshOutput<MeshBase>
{
 public:
  /**
   * Constructor.
   */
  PostscriptIO (const MeshBase& mesh);

  /**
   * Destructor.
   */
  virtual ~PostscriptIO ();

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string& );
  
  
};

#endif // #ifndef __postscript_io_h__
