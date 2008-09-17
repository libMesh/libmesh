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
#include "quadrature_monomial.h"
#include "quadrature_gauss.h"


void QMonomial::init_3D(const ElemType _type,
			unsigned int p)
{

  switch (_type)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
	switch(_order + 2*p)
	  {

	    // By default: construct and use a Gauss quadrature rule
	  default:
	    {
	      // Break out and fall down into the default: case for the
	      // outer switch statement.
	      break;
	    }
	    
	  } // end switch(_order + 2*p)
      } // end case HEX8/20/27

      
      // By default: construct and use a Gauss quadrature rule
    default:
      {
	QGauss gauss_rule(3, _order);
	gauss_rule.init(_type, p);

	// Swap points and weights with the about-to-be destroyed rule.
	_points.swap (gauss_rule.get_points() );
	_weights.swap(gauss_rule.get_weights());

	return;
      }
    } // end switch (_type)
}
