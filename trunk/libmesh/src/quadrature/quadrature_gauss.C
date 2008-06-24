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


#include "quadrature_gauss.h"

// See the files:

// quadrature_gauss_1D.C
// quadrature_gauss_2D.C
// quadrature_gauss_3D.C

// for implementation of specific element types.


void QGauss::dunavant_rule(const Real rule_data[][4],
			   const unsigned int n_pts)
{
  // The input data array has 4 columns.  The first 3 are the permutation points.
  // The last column is the weights for a given set of permutation points.  A zero
  // in two of the first 3 columns implies the point is a 1-permutation (centroid).
  // A zero in one of the first 3 columns implies the point is a 3-permutation.
  // Otherwise each point is assumed to be a 6-permutation.

  // Always insert into the points & weights vector relative to the offset 
  unsigned int offset=0;

  
  for (unsigned int p=0; p<n_pts; ++p)
    {

      // There must always be a non-zero entry to start the row
      libmesh_assert(rule_data[p][0] != 0.0);
      
      // A zero weight may imply you did not set up the raw data correctly
      libmesh_assert(rule_data[p][3] != 0.0);

      // What kind of point is this?
      // One non-zero entry in first 3 cols   ? 1-perm (centroid) point = 1
      // Two non-zero entries in first 3 cols ? 3-perm point            = 3
      // Three non-zero entries               ? 6-perm point            = 6
      unsigned int pointtype=1;
      
      if (rule_data[p][1]!=0.0)      
	{
	  if (rule_data[p][2]!=0.0)
	    pointtype = 6;
	  else
	    pointtype = 3;
	}
      
      switch (pointtype)
	{
	case 1:
	  {
	    // Be sure we have enough space to insert this point
	    libmesh_assert(offset + 0 < _points.size());
	    
	    // The point has only a single permutation (the centroid!)
	    _points[offset  + 0] = Point(rule_data[p][0], rule_data[p][0]);

	    // The weight is always the last entry in the row.
	    _weights[offset + 0] = rule_data[p][3];

	    offset += 1;
	    break;
	  }

	case 3:
	  {
	    // Be sure we have enough space to insert these points
	    libmesh_assert(offset + 2 < _points.size());
	    
	    // Here it's understood the second entry is to be used twice, and
	    // thus there are three possible permutations.
	    _points[offset + 0] = Point(rule_data[p][0], rule_data[p][1]);
	    _points[offset + 1] = Point(rule_data[p][1], rule_data[p][0]);
	    _points[offset + 2] = Point(rule_data[p][1], rule_data[p][1]);

	    for (unsigned int j=0; j<3; ++j)
	      _weights[offset + j] = rule_data[p][3];
	    
	    offset += 3;
	    break;
	  }

	case 6:
	  {
	    // Be sure we have enough space to insert these points
	    libmesh_assert(offset + 5 < _points.size());
	    
	    // Three individual entries with six permutations.
	    _points[offset + 0] = Point(rule_data[p][0], rule_data[p][1]);
	    _points[offset + 1] = Point(rule_data[p][0], rule_data[p][2]);
	    _points[offset + 2] = Point(rule_data[p][1], rule_data[p][0]);
	    _points[offset + 3] = Point(rule_data[p][1], rule_data[p][2]);
	    _points[offset + 4] = Point(rule_data[p][2], rule_data[p][0]);
	    _points[offset + 5] = Point(rule_data[p][2], rule_data[p][1]);

	    for (unsigned int j=0; j<6; ++j)
	      _weights[offset + j] = rule_data[p][3];
	    
	    offset += 6;
	    break;
	  }

	default:
	  {
	    std::cerr << "Don't know what to do with this many permutation points!" << std::endl;
	    libmesh_error();
	  }
	}
    }
}
