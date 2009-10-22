// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



// Local includes
#include "quadrature_simpson.h"



void QSimpson::init_2D(const ElemType _type,
                       unsigned int)
{
#if LIBMESH_DIM > 1
  
  //-----------------------------------------------------------------------
  // 2D quadrature rules
  switch (_type)
    {


      //---------------------------------------------
      // Quadrilateral quadrature rules
    case QUAD4:
    case QUAD8:
    case QUAD9:
      {
	// We compute the 2D quadrature rule as a tensor
	// product of the 1D quadrature rule.
	QSimpson q1D(1);
	q1D.init(EDGE2);
	tensor_product_quad( q1D );
	return;
      }

	    
      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRI6:
      {
	// I'm not sure if you would call this Simpson's
	// rule for triangles.  What it *Really* is is
	// four trapezoidal rules combined to give a six
	// point rule.  The points lie at the nodal locations
	// of the TRI6, so you can get diagonal element
	// stiffness matrix entries for quadratic elements.
	// This rule should be able to integrate a little
	// better than linears exactly.
	
 	_points.resize(6);
 	_weights.resize(6);
	
 	_points[0](0) = 0.;
 	_points[0](1) = 0.;
	
 	_points[1](0) = 1.;
 	_points[1](1) = 0.;
	
 	_points[2](0) = 0.;
 	_points[2](1) = 1.;
	
 	_points[3](0) = 0.5;
 	_points[3](1) = 0.;
	
 	_points[4](0) = 0.;
 	_points[4](1) = 0.5;
	
 	_points[5](0) = 0.5;
 	_points[5](1) = 0.5;
	
 	_weights[0] = 0.041666666666666666666666666667; // 1./24.
 	_weights[1] = 0.041666666666666666666666666667; // 1./24.
 	_weights[2] = 0.041666666666666666666666666667; // 1./24.
 	_weights[3] = 0.125;                            // 1./8.
 	_weights[4] = 0.125;                            // 1./8.
 	_weights[5] = 0.125;                            // 1./8.
	
	return;
      }

      
      //---------------------------------------------
      // Unsupported type
    default:
      {
	std::cerr << "Element type not supported!:" << _type << std::endl;
	libmesh_error();
      }
    }

  libmesh_error();

  return;

#endif
}
