// $Id: centroid_partitioner.h,v 1.1 2003-08-22 19:59:50 jwpeterson Exp $

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


#ifndef __centroid_partitioner_h__
#define __centroid_partitioner_h__

// C++ includes
#include <utility> // pair
#include <vector>

// Local includes
#include "partitioner.h"
#include "point.h"


// Forward declarations
class Elem;

/**
 * The centroid partitioner partitions simply based on the
 * locations of element centroids.  You must define what
 * you mean by "less than" for the list of element centroids, e.g.
 * if you only care about distance in the z-direction, you would
 * define "less than" differently than if you cared about radial
 * distance.
 *
 * @author John W. Peterson, 2003
 */



// CentroidPartitioner class definition
class CentroidPartitioner : public Partitioner
{
public:

  
  /**
   * A typedef which is reserved only for use within
   * this class.  If X is chosen, then centroid locations
   * will be sorted according to 
   */
  enum CentroidSortMethod {X=0,
			   Y,
			   Z,
			   RADIAL,
			   INVALID_METHOD};

  /**
   * Constructor.  Requires a reference to a MeshBase
   * and constructs the base class.
   */
  CentroidPartitioner (MeshBase& mesh, CentroidSortMethod sm);

  /**
   * Partitions the mesh into n subdomains.  This is
   * a required interface for the class.
   */
  virtual void partition (const unsigned int n);

  /**
   * Repartitions the mesh.  This is a required interface
   * for the class.  At this point, we'll just
   * call this->partition().  In the future, we can
   * optimize it maybe.
   */
  virtual void repartition(const unsigned int n);
  
private:

  /**
   * Computes a list of element centroids for the mesh.
   * This list will be kept around in case a repartition
   * is desired.
   */
  void compute_centroids();

  /**
   * Partition the list of centroids based on the
   * x-coordinate of the centroid.  This function provides
   * a function which may be passed to the std::sort
   * routine for sorting the elements by centroid.
   */
  static bool sort_x(const std::pair<Point, Elem*>& lhs,
		     const std::pair<Point, Elem*>& rhs);

  /**
   * Partition the list of centroids based on the
   * y-coordinate of the centroid.  This function provides
   * a function which may be passed to the std::sort
   * routine for sorting the elements by centroid.
   */
  static bool sort_y(const std::pair<Point, Elem*>& lhs,
		     const std::pair<Point, Elem*>& rhs);

  /**
   * Partition the list of centroids based on the
   * z-coordinate of the centroid.  This function provides
   * a function which may be passed to the std::sort
   * routine for sorting the elements by centroid.
   */
  static bool sort_z(const std::pair<Point, Elem*>& lhs,
		     const std::pair<Point, Elem*>& rhs);

  
  /**
   * Store a flag which tells which type of
   * sort method we are using.
   */
  const CentroidSortMethod _sort_method;


  
  /**
   * Vector which holds pairs of centroids and
   * their respective element pointers.
   */
  std::vector<std::pair<Point, Elem*> > _elem_centroids;
};


#endif
