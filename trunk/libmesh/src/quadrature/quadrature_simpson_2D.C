// $Id: quadrature_simpson_2D.C,v 1.1 2003-01-20 16:31:45 jwpeterson Exp $

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



// Local includes
#include "quadrature_simpson.h"



void QSimpson::init_2D(const ElemType _type)
{
#if DIM > 1
  
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

	tensor_product_quad( &q1D );
	
	return;
      };

	    
      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRI6:
      {
	// TODO:[JWP] Fix Triangle rule for Simpson quadrature
	
// 	_points.resize(3);
// 	_weights.resize(3);
	
// 	_points[0](0) = 0.;
// 	_points[0](1) = 0.;
	
// 	_points[1](0) = 1.;
// 	_points[1](1) = 0.;
	
// 	_points[2](0) = 0.;
// 	_points[2](1) = 1.;
	
	
// 	_weights[0] = 1./6.;
// 	_weights[1] = 1./6.;
// 	_weights[2] = 1./6.;
	
	return;
      };

      
      //---------------------------------------------
      // Unsupported type
    default:
      {
	std::cerr << "Element type not supported!:" << _type << std::endl;
	error();
      }
    };

  error();

  return;

#endif
};




void QSimpson::init_2D(const ElemType _type,
		       const unsigned int side)
{
#if DIM > 1
  
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
	QSimpson q1D(1);
	q1D.init(EDGE2);
	side_rule_quad(&q1D, side);
	return;
      } 

	    
      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRI6:
      {
	QSimpson q1D(1);
	q1D.init(EDGE2);
	side_rule_tri(&q1D, side);
	return;
      }
      
    default:
      {
 	std::cerr << "Element type not supported!:" << _type << std::endl;
 	error();
      }
    };
  
  //  error();
  
  return;

#endif
};

