// $Id: metis_partitioner.h,v 1.3 2003-08-26 22:58:45 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __metis_partitioner_h__
#define __metis_partitioner_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "partitioner.h"



/**
 * The \p MetisPartitioner uses the Metis graph partitioner
 * to partition the elements.
 */

// ------------------------------------------------------------
// MetisPartitioner class definition
class MetisPartitioner : public Partitioner
{
 public:

  /**
   * Constructor.  Requires a \p MeshBase.
   */
  MetisPartitioner (MeshBase& mesh) : Partitioner(mesh) {}

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void partition (const unsigned int n = libMesh::n_processors());

private:
};




#endif // #define __sfc_partitioner_h__
