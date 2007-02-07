// $Id: hilbert_sfc_partitioner.h,v 1.4 2005-02-22 22:17:34 jwpeterson Exp $

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



#ifndef __hilbert_sfc_partitioner_h__
#define __hilbert_sfc_partitioner_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "sfc_partitioner.h"



/**
 * The \p HilbertSFCPartitioner uses a Hilbert space
 * filling curve to partition the elements.
 */

// ------------------------------------------------------------
// HilbertSFCLinearPartitioner class definition
class HilbertSFCPartitioner : public SFCPartitioner
{
 public:

  /**
   * Constructor.
   */
  HilbertSFCPartitioner ()
  {
    this->set_sfc_type ("Hilbert");
  }

protected:
  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase& mesh,
			      const unsigned int n)
  { SFCPartitioner::_do_partition (mesh, n); }


  
private:
  
};


#endif // #define __hilbert_sfc_partitioner_h__
