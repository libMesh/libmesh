// $Id: quadrature_trap_3D.C,v 1.3 2003-01-21 19:24:38 benkirk Exp $

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
#include "quadrature_trap.h"





void QTrap::init_3D(const ElemType _type)
{
#if DIM == 3
  
  //-----------------------------------------------------------------------
  // 3D quadrature rules
  switch (_type)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
	// We compute the 3D quadrature rule as a tensor
	// product of the 1D quadrature rule.
	QTrap q1D(1);
	q1D.init(EDGE2);

	tensor_product_hex( &q1D );
	
	return;
      }


      
      //---------------------------------------------
      // Tetrahedral quadrature rules
    case TET4:
    case TET10:
      {
	_points.resize(4);
	_weights.resize(4);
	
	_points[0](0) = 0.;
	_points[0](1) = 0.;
	_points[0](2) = 0.;
	
	_points[1](0) = 1.;
	_points[1](1) = 0.;
	_points[1](2) = 0.;
	
	_points[2](0) = 0.;
	_points[2](1) = 1.;
	_points[2](2) = 0.;
	
	_points[3](0) = 0.;
	_points[3](1) = 0.;
	_points[3](2) = 1.;
	
	
	
	_weights[0] = .0416666666666666666666666666666666666666666667;
	_weights[1] = _weights[0];
	_weights[2] = _weights[0];
	_weights[3] = _weights[0];
	
	return;
      };
      
      
      
      //---------------------------------------------
      // Prism quadrature rules
    case PRISM6:
    case PRISM18:
      {
	// We compute the 3D quadrature rule as a tensor
	// product of the 1D quadrature rule and a 2D
	// triangle quadrature rule
	    
	QTrap q1D(1);
	QTrap q2D(2);

	// Initialize 
	q1D.init(EDGE2);
	q2D.init(TRI3);

	tensor_product_prism(&q1D, &q2D);
	
	return;
      };

      
      //---------------------------------------------
      // Unsupported type
    default:
      {
	std::cerr << "ERROR: Unsupported type: " << _type << std::endl;
	error();
      }
    };

  error();

  return;
  
#endif
};







void QTrap::init_3D(const ElemType _type,
		    const unsigned int side)
{
#if DIM == 3
  
  //-----------------------------------------------------------------------
  // 3D quadrature rules
  switch (_type)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
	QTrap q2D(2);
	q2D.init(QUAD4);
	side_rule_hex(&q2D, side);
	return;
      }

      
	
      //---------------------------------------------
      // Tetrahedral quadrature rules
    case TET4:
    case TET10:
      {
	QTrap q2D(2);
	q2D.init(TRI3);
	side_rule_tet(&q2D, side);
	return;
      }

      
	    
      //---------------------------------------------
      // Prism quadrature rules
    case PRISM6:
    case PRISM18:
      {
	QTrap q2D(2);

	// There is no need to initialize q2D,
	// the side could be either a quad or a tri.
	side_rule_prism(&q2D, side);
	return;
      };


      
      //---------------------------------------------
      // Unsupported type
    default:
      {
	std::cerr << "ERROR: Unsupported type: " << _type << std::endl;
	error();
      }
    };

  error();
  
  return;
  
#endif
};

