// $Id: quadrature_clough_2D.C,v 1.4 2007-10-21 20:48:53 benkirk Exp $

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



// Local includes
#include "quadrature_clough.h"
#include "quadrature_gauss.h"


void QClough::init_2D(const ElemType _type,
                      unsigned int p)
{
#if DIM > 1
  QGauss gauss_rule(2, _order);
  gauss_rule.init(TRI6, p);
  
  //-----------------------------------------------------------------------
  // 2D quadrature rules
  switch (_type)
    {

      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRI6:
      {
	std::vector<Point> &gausspoints = gauss_rule.get_points();
	std::vector<Real> &gaussweights = gauss_rule.get_weights();
	unsigned int numgausspts = gausspoints.size();
	_points.resize(numgausspts*3);
	_weights.resize(numgausspts*3);
        for (unsigned int i = 0; i != numgausspts; ++i)
          {
	    _points[3*i](0) = gausspoints[i](0) +
			    gausspoints[i](1) / 3.;
	    _points[3*i](1) = gausspoints[i](1) / 3.;
	    _points[3*i+1](0) = gausspoints[i](1) / 3.;
	    _points[3*i+1](1) = gausspoints[i](0) +
			    gausspoints[i](1) / 3.;
	    _points[3*i+2](0) = 1./3. +
			    gausspoints[i](0) * 2./3. -
			    gausspoints[i](1) / 3.;
	    _points[3*i+2](1) = 1./3. - 
			    gausspoints[i](0) / 3. +
			    gausspoints[i](1) * 2./3.;
	    _weights[3*i] = gaussweights[i] / 3.;
	    _weights[3*i+1] = _weights[3*i];
	    _weights[3*i+2] = _weights[3*i];
          }
	return;
      }

	    
      //---------------------------------------------
      // Unsupported type
    default:
      {
	std::cerr << "Element type not supported!:" << _type << std::endl;
	error();
      }
    }

  error();

  return;

#endif
}
