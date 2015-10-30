// $Id: partitioner.h,v 1.3 2004-05-11 20:29:06 jwpeterson Exp $

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



#ifndef __partitioner_h__
#define __partitioner_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "libmesh.h"

// Forward Declarations
class MeshBase;



/**
 * The \p Partitioner class provides a uniform interface for
 * partitioning algorithms.  It takes a reference to a \p MeshBase
 * object as input, which it will partition into a number of
 * subdomains. 
 */

// ------------------------------------------------------------
// Partitioner class definition
class Partitioner
{
 public:

  /**
   * Constructor.
   */
  Partitioner () {}
  
  /**
   * Destructor. Virtual so that we can derive from this class.
   */
  virtual ~Partitioner() {}

  /**
   * Partition the \p MeshBase into \p n parts.  If the
   * user does not specify a number of pieces into which the
   * mesh should be partitioned, then the default behavior
   * of the partitioner is to partition according to the number
   * of processors defined in libMesh::n_processors().
   * The partitioner currently does not modify the subdomain_id
   * of each element.  This number is reserved for things like
   * material properties, etc.
   */
  void partition (MeshBase& mesh,
		  const unsigned int n=libMesh::n_processors());

  /**
   * Repartitions the \p MeshBase into \p n parts.  This
   * is required since some partitoning algorithms can repartition
   * more efficiently than computing a new partitioning from scratch.
   * The default behavior is to simply call this->partition(n)
  */
  void repartition (MeshBase& mesh,
		    const unsigned int n=libMesh::n_processors());

  
protected:

  /**
   * Trivially "partitions" the mesh for one processor.
   * Simply loops through the elements and assigns all of them
   * to processor 0.  Is is provided as a separate function
   * so that derived classes may use it without reimplementing it.
   */
  void single_partition (MeshBase& mesh);

  /**
   * This is the actual partitioning method which must be overloaded
   * in derived classes.  It is called via the public partition()
   * method above by the user.
   */
  virtual void _do_partition(MeshBase& mesh,
			     const unsigned int n) = 0;

  /**
   * This is the actual re-partitioning method which can be overloaded
   * in derived classes.  Note that the default behavior is to simply
   * call the partition function.
   */
  virtual void _do_repartition (MeshBase& mesh,
				const unsigned int n) { this->_do_partition (mesh, n); }
};




#endif
