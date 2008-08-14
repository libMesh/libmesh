// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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

#ifndef __off_io_h__
#define __off_io_h__

// Local includes
#include "libmesh_common.h"
#include "mesh_input.h"

// Forward declarations
class MeshBase;



/**
 * This class is repsonsible for reading an unstructured,
 * triangulated surface in the
 * standard OFF OOGL format.
 */

// ------------------------------------------------------------
// OFFIO class definition
class OFFIO : public MeshInput<MeshBase>
{
public:
  /**
   *  Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements.
   */
  OFFIO (MeshBase&);

  /**
   * Reads in an OFF OOGL data file based on the string
   * you pass it.
   */
  virtual void read (const std::string& name);

private:
  /**
   * Implementation of the read() function.  This function
   * is called by the public interface function and implements
   * reading the file.
   */
  virtual void read_stream (std::istream& in);
};


// ------------------------------------------------------------
// OFFIO inline members
inline
OFFIO::OFFIO (MeshBase& mesh) :
  MeshInput<MeshBase> (mesh)
{}



#endif
