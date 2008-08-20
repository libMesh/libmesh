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


// C++ includes

// Local includes
#include "quadrature.h"

void QBase::init(const ElemType t,
                 unsigned int p)
{
  // check to see if we have already
  // done the work for this quadrature rule
  if (t == _type && p == _p_level)
    return;
  else
    {
      _type = t;
      _p_level = p;
    }
    
  
  
  switch(_dim)
    {
    case 0:
      this->init_0D(_type,_p_level);

      return;
      
    case 1:
      this->init_1D(_type,_p_level);

      return;
      
    case 2:
      this->init_2D(_type,_p_level);

      return;

    case 3:
      this->init_3D(_type,_p_level);

      return;

    default:
      libmesh_error();
    }
}



void QBase::init_0D(const ElemType,
                    unsigned int)
{
  _points.resize(1);
  _weights.resize(1);
  _points[0] = Point(0.);
  _weights[0] = 1.0;
}



void QBase::scale(std::pair<Real, Real> old_range,
		  std::pair<Real, Real> new_range)
{
  // Make sure we are in 1D
  libmesh_assert(_dim == 1);
  
  // Make sure that we have sane ranges 
  libmesh_assert(new_range.second > new_range.first);
  libmesh_assert(old_range.second > old_range.first);

  // Make sure there are some points
  libmesh_assert(_points.size() > 0);

  // We're mapping from old_range -> new_range 
  for (unsigned int i=0; i<_points.size(); i++)
    {
      _points[i](0) =
	(_points[i](0) - old_range.first) *
	(new_range.second - new_range.first) /
	(old_range.second - old_range.first) +
	new_range.first;
    }

  // Compute the scale factor and scale the weights
  const Real scfact = (new_range.second - new_range.first) /
                      (old_range.second - old_range.first);

  for (unsigned int i=0; i<_points.size(); i++)
    _weights[i] *= scfact;
}




void QBase::tensor_product_quad(const QBase& q1D)
{
  
  const unsigned int n_points = q1D.n_points();
  
  _points.resize(n_points * n_points);
  
  _weights.resize(n_points * n_points);
  
  unsigned int qp=0;
  
  for (unsigned int j=0; j<n_points; j++)
    for (unsigned int i=0; i<n_points; i++)
      {
	_points[qp](0) = q1D.qp(i)(0);
	_points[qp](1) = q1D.qp(j)(0);
	
	_weights[qp] = q1D.w(i)*q1D.w(j);
	
	qp++;
      }
}





void QBase::tensor_product_hex(const QBase& q1D)
{
  const unsigned int n_points = q1D.n_points();
  
  _points.resize(n_points * n_points * n_points);
  
  _weights.resize(n_points * n_points * n_points);
  
  unsigned int qp=0;
  
  for (unsigned int k=0; k<n_points; k++)
    for (unsigned int j=0; j<n_points; j++)
      for (unsigned int i=0; i<n_points; i++)
	{
	  _points[qp](0) = q1D.qp(i)(0);
	  _points[qp](1) = q1D.qp(j)(0);
	  _points[qp](2) = q1D.qp(k)(0);
	  
	  _weights[qp] = q1D.w(i) * q1D.w(j) * q1D.w(k);
	  
	  qp++;
	}
}




void QBase::tensor_product_prism(const QBase& q1D, const QBase& q2D)
{
  const unsigned int n_points1D = q1D.n_points();
  const unsigned int n_points2D = q2D.n_points();
  
  _points.resize  (n_points1D * n_points2D);
  _weights.resize (n_points1D * n_points2D);

  unsigned int qp=0;
  
  for (unsigned int j=0; j<n_points1D; j++)
    for (unsigned int i=0; i<n_points2D; i++)
      {
	_points[qp](0) = q2D.qp(i)(0);
	_points[qp](1) = q2D.qp(i)(1);
	_points[qp](2) = q1D.qp(j)(0);

	_weights[qp] = q2D.w(i) * q1D.w(j);

	qp++;
      }
  
}




std::ostream& operator << (std::ostream& os, const QBase& q)
{
  q.print_info(os);
  return os;
}
