// $Id: mesh_communication.h,v 1.4 2003-09-16 15:59:31 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __mesh_communication_h__
#define __mesh_communication_h__



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "mesh_common.h"


// Forward declarations
class MeshBase;
class Mesh;
class BoundaryInfo;



/**
 * This is the \p MeshCommunication class.  It handles all the details
 * of communicating mesh information from one processor to another.  All
 * parallelization of the \p Mesh data structures is done via this class.
 *
 * @author Benjamin S. Kirk, 2003
 */


// ------------------------------------------------------------
// MeshCommunication class definition
class MeshCommunication
{
public:

  /**
   * Constructor.
   */
  MeshCommunication () {}

  /**
   * Destructor.
   */
  ~MeshCommunication () {}

  /**
   * Clears all data structures and returns to a pristine state.
   */
  void clear ();

  /**
   * This method takes a mesh (which is assumed to reside on
   * processor 0) and distributes it to all the other processors.
   * It also distributes any boundary information the mesh has
   * associated with it.
   */
  void distribute (Mesh& );

  /**
   * Finds all the processors that may contain
   * elements that neighbor my elements.  This list
   * is guaranteed to include all processors that border
   * any of my elements, but may include additional ones as
   * well.  This method computes bounding boxes for the
   * elements on each processor and checks for overlaps.
   */
  void find_neighboring_processors() {}

  
private:

  /**
   *
   */
  void distribute_mesh (MeshBase& );

  /**
   *
   */
  void distribute_bcs (MeshBase&, BoundaryInfo&);
};



//--------------------------------------------------------------
// MeshCommunication inline members


#endif // end #ifndef __mesh_communication_h__ 
