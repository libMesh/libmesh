// $Id: partitioner.h,v 1.2 2003-07-25 20:58:24 benkirk Exp $

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



#ifndef __partitioner_h__
#define __partitioner_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------


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
   * Constructor.  Requires a \p MeshBase.
   */
  Partitioner (MeshBase& mesh) : _mesh(mesh) {}
  
  /**
   * Destructor. Virtual so that we can derive from this class.
   */
  virtual ~Partitioner() {}

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void partition (const unsigned int n) = 0;

  /**
   * Repartitions the \p MeshBase into \p n subdomains.  This
   * is required since some partitoning algorithms can repartition
   * more efficiently than computing a new partitioning from scratch.
   * IF this is not the case then your derived class should implement
   * \p repartition() to simply call \p this->partition().
   */
  virtual void repartition (const unsigned int n) = 0;

  
protected:

  /**
   * Trivially "partitions" the mesh for one processor.
   * Simply loops through the elements and assigns all of them
   * to processor 0.  Is is provided as a separate function
   * so that derived classes may use it without reimplementing it.
   */
  void single_partition ();
  
  /**
   * A reference to the \p MeshBase we will partition.
   */
  MeshBase& _mesh;
};




#endif
