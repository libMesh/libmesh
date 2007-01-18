
// $Id: patch.h,v 1.1 2007-01-18 22:24:47 roystgnr Exp $

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



#ifndef __patch_h__
#define __patch_h__

// C++ includes
#include <vector>
#include <string>
#include <set>

// Forward Declarations
class Elem;



/**
 * This class implements useful utility functions for a patch of
 * elements
 *
 * @author Roy H. Stogner, 2007.
 * based on code by Varis Carey, Benjamin S. Kirk, 2004.
 */
class Patch : public std::set<const Elem *>
{
public:

  /**
   * Constructor.
   */
  Patch() {}
  
  /**
   * Destructor.  
   */
  ~Patch() {}

  /**
   * This function finds all elements which touch the current patch at
   * a face, and adds them to the patch.
   */
  void add_face_neighbors();

  /**
   * This function finds all elements on the current processor which
   * touch the current patch at a face, and adds them to the patch.
   */
  void add_local_face_neighbors();

  /**
   * This function finds all elements which touch the current patch at
   * a vertex, and adds them to the patch.
   */
  void add_vertex_neighbors();

  /**
   * This function finds all elements on the current processor which
   * touch the current patch at a vertex, and adds them to the patch.
   */
  void add_local_vertex_neighbors();

  /**
   * Pointer to Member Function typedef
   */
  typedef void (Patch::*PMF)();

  /**
   * Erases any elements in the current patch, then builds a new patch
   * containing element \p elem by repeated addition of neighbors on
   * the current processor.  This procedure is repeated until the
   * number of elements meets or exceeds \p target_patch_size, or
   * until the patch has no more local neighbors.
   */
  void build_around_element(const Elem* elem,
			    const unsigned int target_patch_size = 10,
			    PMF patchtype = &Patch::add_local_vertex_neighbors);
};


#endif // #define __patch_h__
