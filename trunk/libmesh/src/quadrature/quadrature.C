// $Id: quadrature.C,v 1.7 2003-02-06 06:02:42 jwpeterson Exp $

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


#include <numeric>

#include "quadrature.h"

void QBase::init(const ElemType t)
{
  // check to see if we have already
  // done the work for this quadrature rule
  if (t == _type)
    return;
  else
    _type = t;
    
  
  
  switch(_dim)
    {
    case 1:
      init_1D(_type);

      return;
      
    case 2:
      init_2D(_type);

      return;

    case 3:
      init_3D(_type);

      return;

    default:
      error();
    };

  error();
  
  return;
};



void QBase::init(const ElemType t,
		 const unsigned int s)
{
  static unsigned int side = static_cast<unsigned int>(-1);

  
  // check to see if we have already
  // done the work for this quadrature rule
  if ((t == _type) &&
      (s == side)) 
    {
      //here();
      return;
    }
  else
    {
      _type = t;
      side = s;
    }


  switch (_dim)
    {
    case 2:
      init_2D(_type,side);

      return;

    case 3:
      init_3D(_type,side);

      return;


    default:
      error();

    };

  error();

  return;
};




void QBase::scale(std::pair<Real, Real> old_range,
		  std::pair<Real, Real> new_range)
{
  // Make sure we are in 1D
  assert(_dim == 1);
  
  // Make sure that we have sane ranges 
  assert(new_range.second > new_range.first);
  assert(old_range.second > old_range.first);

  // Make sure there are some points
  assert(_points.size() > 0);

  // We're mapping from old_range -> new_range 
  for (unsigned int i=0; i<_points.size(); i++)
    {
      //_points[i](0) = (_points[i](0)-y1) * (x2-x1) / (y2-y1) + x1;
      _points[i](0) = (_points[i](0) - old_range.first) *
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




void QBase::tensor_product_quad(QBase* q1D)
{
  
  const unsigned int n_points = q1D->n_points();
  
  _points.resize(n_points * n_points);
  
  _weights.resize(n_points * n_points);
  
  unsigned int qp=0;
  
  for (unsigned int j=0; j<n_points; j++)
    for (unsigned int i=0; i<n_points; i++)
      {
	_points[qp](0) = q1D->qp(i)(0);
	_points[qp](1) = q1D->qp(j)(0);
	
	_weights[qp] = q1D->w(i)*q1D->w(j);
	
	qp++;
      };

}




void QBase::tensor_product_tri(QBase* gauss1D, QBase* jac1D)
{
  
  // Both rules should be of the same order
  assert(gauss1D->n_points() == jac1D->n_points());

  // Save the number of points as a convenient variable
  const unsigned int n_points = gauss1D->n_points();
  
  // Both rules should be between x=0 and x=1
  assert(gauss1D->qp(0)(0) >= 0.0); assert(gauss1D->qp(n_points-1)(0) <= 1.0);
  assert(jac1D->qp(0)(0)   >= 0.0); assert(jac1D->qp(n_points-1)(0) <= 1.0);

  // Resize the points and weights vectors
  _points.resize(n_points * n_points);
  _weights.resize(n_points * n_points);

  // Compute the conical product
  unsigned int gp = 0;
  for (unsigned int i=0; i<n_points; i++)
    for (unsigned int j=0; j<n_points; j++)
      {
	_points[gp](0) = jac1D->qp(j)(0);                          //s[j];
	_points[gp](1) = gauss1D->qp(i)(0) * (1.-jac1D->qp(j)(0)); //r[i]*(1.-s[j]);
	_weights[gp]   = gauss1D->w(i) * jac1D->w(j);              //A[i]*B[j];
	gp++;
      }
}




void QBase::tensor_product_hex(QBase* q1D)
{
  const unsigned int n_points = q1D->n_points();
  
  _points.resize(n_points * n_points * n_points);
  
  _weights.resize(n_points * n_points * n_points);
  
  unsigned int qp=0;
  
  for (unsigned int k=0; k<n_points; k++)
    for (unsigned int j=0; j<n_points; j++)
      for (unsigned int i=0; i<n_points; i++)
	{
	  _points[qp](0) = q1D->qp(i)(0);
	  _points[qp](1) = q1D->qp(j)(0);
	  _points[qp](2) = q1D->qp(k)(0);
	  
	  _weights[qp] = q1D->w(i) * q1D->w(j) * q1D->w(k);
	  
	  qp++;
	};
}




void QBase::tensor_product_prism(QBase* q1D, QBase* q2D)
{
  const unsigned int n_points1D = q1D->n_points();
  const unsigned int n_points2D = q2D->n_points();
  
  _points.resize  (n_points1D * n_points2D);
  _weights.resize (n_points1D * n_points2D);

  unsigned int qp=0;

  for (unsigned int j=0; j<n_points1D; j++)
    for (unsigned int i=0; i<n_points2D; i++)
      {
	_points[qp](0) = q2D->qp(i)(0);
	_points[qp](1) = q2D->qp(i)(1);
	_points[qp](2) = q1D->qp(j)(0);

	_weights[qp] = q2D->w(i) * q1D->w(j);

	qp++;
      };
  
}




void QBase::tensor_product_tet(QBase* gauss1D, QBase* jacA1D, QBase* jacB1D)
{
  here();

  // All rules should be of the same order
  assert(gauss1D->n_points() == jacA1D->n_points());
  assert(jacA1D->n_points()  == jacB1D->n_points());
  
  // Save the number of points as a convenient variable
  const unsigned int n_points = gauss1D->n_points();
  
  // All rules should be between x=0 and x=1
  assert(gauss1D->qp(0)(0) >= 0.0); assert(gauss1D->qp(n_points-1)(0) <= 1.0);
  assert(jacA1D->qp(0)(0)  >= 0.0); assert(jacA1D->qp(n_points-1)(0)  <= 1.0);
  assert(jacB1D->qp(0)(0)  >= 0.0); assert(jacB1D->qp(n_points-1)(0)  <= 1.0);

  // Resize the points and weights vectors
  _points.resize(n_points * n_points * n_points);
  _weights.resize(n_points * n_points * n_points);

  // Compute the conical product
  unsigned int gp = 0;
  for (unsigned int i=0; i<n_points; i++)
    for (unsigned int j=0; j<n_points; j++)
      for (unsigned int k=0; k<n_points; k++)
      {
	_points[gp](0) = jacB1D->qp(k)(0);                                                  //t[k];
	_points[gp](1) = jacA1D->qp(j)(0)  * (1.-jacB1D->qp(k)(0));                         //s[j]*(1.-t[k]);
	_points[gp](2) = gauss1D->qp(i)(0) * (1.-jacA1D->qp(j)(0)) * (1.-jacB1D->qp(k)(0)); //r[i]*(1.-s[j])*(1.-t[k]);
	_weights[gp]   = gauss1D->w(i)     * jacA1D->w(j)          * jacB1D->w(k);          //A[i]*B[j]*C[k];
	gp++;
      }
  
  std::cout << "Sum of weights is: "
	    << std::accumulate(_weights.begin(), _weights.end(), 0., std::plus<Real>()) << std::endl;
}



void QBase::side_rule_quad(QBase* q1D, unsigned int side)
{

  const unsigned int n_points = q1D->n_points();
  
  _points.resize(n_points);
  _weights.resize(n_points);

  switch (side)
    {
    case 0:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = q1D->qp(qp)(0);
	    _points[qp](1) = -1.;

	    _weights[qp] = q1D->w(qp);
	  }
	return;
      }

    case 1:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = 1.;
	    _points[qp](1) = q1D->qp(qp)(0);

	    _weights[qp] = q1D->w(qp);
	  }
	return;
      }



    case 2:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = -q1D->qp(qp)(0);
	    _points[qp](1) = 1.;

	    _weights[qp] = q1D->w(qp);
	  }
	return;
      }

    case 3:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = -1.;
	    _points[qp](1) = -q1D->qp(qp)(0);

	    _weights[qp] = q1D->w(qp);
	  }
	return;
      }


    default:
      error();
    }
  
}




void QBase::side_rule_tri(QBase* q1D, unsigned int side)
{

  const unsigned int n_points = q1D->n_points();
  
  _points.resize(n_points);
  _weights.resize(n_points);


  switch (side)
    {
    case 0:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = .5*(q1D->qp(qp)(0) + 1.);
	    _points[qp](1) = 0.;

	    _weights[qp] = q1D->w(qp);
	  }
	return;
      }

    case 1:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = .5*(q1D->qp(n_points - (qp +1))(0) + 1.);
	    _points[qp](1) = 1. - _points[qp](0);

	    _weights[qp] = q1D->w(qp);
	  }
	return;
      }

    case 2:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = 0.;
	    _points[qp](1) = .5*(q1D->qp(n_points - (qp +1))(0) + 1.);

	    _weights[qp] = q1D->w(qp);
	  }
	return;
      }


    default:
      error();
    };
  
}




void QBase::side_rule_hex(QBase* q2D, unsigned int side)
{
  const unsigned int n_points = q2D->n_points();
  
  _points.resize(n_points);
  _weights.resize(n_points);
  
  switch (side)
    {

      // the face at z = -1
    case 0:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = q2D->qp(qp)(1);
	    _points[qp](1) = q2D->qp(qp)(0);
	    _points[qp](2) = -1.;

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }

      // the face at y = -1
    case 1:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = q2D->qp(qp)(0);
	    _points[qp](1) = -1.;
	    _points[qp](2) = q2D->qp(qp)(1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }



      // the face at x = 1
    case 2:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = 1.;
	    _points[qp](1) = q2D->qp(qp)(0);
	    _points[qp](2) = q2D->qp(qp)(1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }
      

      // the face at y = 1
    case 3:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = -q2D->qp(qp)(0);
	    _points[qp](1) = 1.;
	    _points[qp](2) = q2D->qp(qp)(1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }



      // the face at x = -1
    case 4:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = -1.;
	    _points[qp](1) = -q2D->qp(qp)(0);
	    _points[qp](2) = q2D->qp(qp)(1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }
      

      // the face at z = 1
    case 5:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = q2D->qp(qp)(0);
	    _points[qp](1) = q2D->qp(qp)(1);
	    _points[qp](2) = 1.;

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }



    default:
      error();
    }
      
}





void QBase::side_rule_tet(QBase* q2D, unsigned int side)
{
  const unsigned int n_points = q2D->n_points();
  
  _points.resize(n_points);

  _weights.resize(n_points);


  switch(side)
    {

      // the face at z = 0
    case 0:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = q2D->qp(qp)(0);
	    _points[qp](1) = q2D->qp(qp)(1);
	    _points[qp](2) = 0.;

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }


      // the face at y = 0
    case 1:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = q2D->qp(qp)(0);
	    _points[qp](1) = 0.;
	    _points[qp](2) = q2D->qp(qp)(1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }



      // the face at z = 1 - x - y
    case 2:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = q2D->qp(qp)(0);
	    _points[qp](1) = q2D->qp(qp)(1);
	    _points[qp](2) = 1. - _points[qp](0) - _points[qp](1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }

      
      // the face at x = 0
    case 3:
      {
	for (unsigned int qp=0; qp<n_points; qp++)
	  {
	    _points[qp](0) = 0.;
	    _points[qp](1) = q2D->qp(qp)(0);
	    _points[qp](2) = q2D->qp(qp)(1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }



    default:
      error();
    };
  
}




void QBase::side_rule_prism(QBase* q2D, unsigned int side)
{
  switch(side)
    {

      // the face at z = -1
    case 0:
      {
	q2D->init(TRI3);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) = q2D->qp(qp)(0);
	    _points[qp](1) = q2D->qp(qp)(1);
	    _points[qp](2) = -1.;

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }

      // the face at y = 0
    case 1:
      {
	q2D->init(QUAD4);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) = .5*(q2D->qp(qp)(0) + 1.);
	    _points[qp](1) = 0.;
	    _points[qp](2) = q2D->qp(qp)(1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }

      // the face at x + y = 1
    case 2:
      {
	q2D->init(QUAD4);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) = .5*(q2D->qp(qp)(0) + 1.);
	    _points[qp](1) = 1. - _points[qp](0);
	    _points[qp](2) = q2D->qp(qp)(1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }


      // the face at x = 0
    case 3:
      {
	q2D->init(QUAD4);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) = 0.;
	    _points[qp](1) = .5*(q2D->qp(qp)(0) + 1.);
	    _points[qp](2) = q2D->qp(qp)(1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }
      

      // the face at z = 1
    case 4:
      {
	q2D->init(TRI3);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) = q2D->qp(qp)(0);
	    _points[qp](1) = q2D->qp(qp)(1);
	    _points[qp](2) = 1.;

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }
    default:
      error();
    };
}




void QBase::side_rule_pyramid(QBase* q2D, unsigned int side)
{
  switch(side)
    {
      // the face at z = 0
    case 0:
      {
	q2D->init(QUAD4);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) = q2D->qp(qp)(0);
	    _points[qp](1) = q2D->qp(qp)(1);
	    _points[qp](2) = 0.;

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }

      // the face at z = 1 + y
    case 1:
      {
	q2D->init(TRI3);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	const Real a = -1.;
	const Real b =  1.;

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) =  a*q2D->qp(qp)(0) + b*q2D->qp(qp)(1);
	    _points[qp](1) = -b*q2D->qp(qp)(0) + a*q2D->qp(qp)(1);
	    _points[qp](2) = 1. + _points[qp](1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }


      // the face at z = 1 - x
    case 2:
      {
	q2D->init(TRI3);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	const Real a = 1.;
	const Real b = 1.;

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) =  a*q2D->qp(qp)(0) + b*q2D->qp(qp)(1);
	    _points[qp](1) = -b*q2D->qp(qp)(0) + a*q2D->qp(qp)(1);
	    _points[qp](2) = 1. - _points[qp](0);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }


      // the face at z = 1 - y
    case 3:
      {
	q2D->init(TRI3);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	const Real a = 1.;
	const Real b = -1.;

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) =  a*q2D->qp(qp)(0) + b*q2D->qp(qp)(1);
	    _points[qp](1) = -b*q2D->qp(qp)(0) + a*q2D->qp(qp)(1);
	    _points[qp](2) = 1. - _points[qp](1);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }
      

      // the face at z = 1 + x
    case 4:
      {
	q2D->init(TRI3);

	_points.resize(q2D->n_points());
	_weights.resize(q2D->n_points());

	const Real a = -1.;
	const Real b = -1.;

	for (unsigned int qp=0; qp<q2D->n_points(); qp++)
	  {
	    _points[qp](0) =  a*q2D->qp(qp)(0) + b*q2D->qp(qp)(1);
	    _points[qp](1) = -b*q2D->qp(qp)(0) + a*q2D->qp(qp)(1);
	    _points[qp](2) = 1. + _points[qp](0);

	    _weights[qp] = q2D->w(qp);
	  }

	return;
      }



    default:
      error();
    };
  

  

  
      
}
