// $Id: mesh_input.h,v 1.1 2004-11-17 07:52:17 benkirk Exp $

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



#ifndef __mesh_input_h__
#define __mesh_input_h__


// C++ inludes
#include <istream>
#include <string>

// Local includes
#include "libmesh_common.h"




/**
 * This class defines an abstract interface for \p Mesh input.
 * Specific classes derived from this class actually implement
 * reading various mesh formats.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 * \version $Revision: 1.1 $
 */

// ------------------------------------------------------------
// MeshInput class definition
template <class MT>
class MeshInput
{
 protected:

  /**
   * Default constructor. Will set the _obj to NULL, effectively
   * rendering this object useless.
   */
  MeshInput ();
  
  /**
   * Constructor.  Takes a writeable reference to an object.
   * This is the constructor required to read an object.
   */
  MeshInput (MT&);
  
 public:

  /**
   * Destructor.
   */
  virtual ~MeshInput ();
  
  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string&) = 0;

  
 protected:
  
  /**
   * Returns the object as a writeable reference.
   */
  MT& mesh ();
  
  /**
   * Reads input from \p in, skipping all the lines
   * that start with the character \p comment_start.
   */
  void skip_comment_lines (std::istream& in,
			   const char comment_start);

  
 private:
  

  /**
   * A pointer to a non-const object object.
   * This allows us to read the object from file.
   */ 
  MT* _obj;
};



// ------------------------------------------------------------
// MeshInput inline members
template <class MT>
inline
MeshInput<MT>::MeshInput () :
  _obj (NULL)
{
}



template <class MT>
inline
MeshInput<MT>::MeshInput (MT& obj) :
  _obj (&obj)
{
}



template <class MT>
inline
MeshInput<MT>::~MeshInput ()
{
}



template <class MT>
inline
MT& MeshInput<MT>::mesh ()
{
  if (_obj == NULL) error();
  return *_obj;
}



template <class MT>
void MeshInput<MT>::skip_comment_lines (std::istream &in,
					const char comment_start)
{    
  char c, line[256];
  
  while (in.get(c), c==comment_start) 
    in.getline (line, 255);
  
  // put back first character of
  // first non-comment line
  in.putback (c);
}



#endif // #define __mesh_io_h__
