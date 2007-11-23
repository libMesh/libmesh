// $Id: linear_partitioner.h,v 1.5 2007-10-21 20:48:44 benkirk Exp $

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



#ifndef __linear_partitioner_h__
#define __linear_partitioner_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "partitioner.h"



/**
 * The \p LinearPartitioner is the simplest of all possible partitioners.
 * It takes the element list and splits it into equal-sized chunks assigned
 * to each processor.  Health Warning: THIS PARTITIONER COULD BE ARBITRARILY
 * BAD!!
 */

// ------------------------------------------------------------
// LinearPartitioner class definition
class LinearPartitioner : public Partitioner
{
 public:

  /**
   * Constructor.
   */
  LinearPartitioner () {}

protected:
  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase& mesh,
			      const unsigned int n);

private:
  
};




#endif  // #define __linear_partitioner_h__
