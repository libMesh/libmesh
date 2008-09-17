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


#include "quadrature_conical.h"
#include "quadrature_gauss.h"
#include "quadrature_jacobi.h"

// See also the files:
// quadrature_conical_2D.C
// quadrature_conical_3D.C
// for additional implementation.




// Constructor
QConical::QConical(const unsigned int d,
		   const Order o) : QBase(d,o)
{
}



// Destructor
QConical::~QConical()
{
}



// Builds and scales a Gauss rule and a Jacobi rule.
// Then combines them to compute points and weights
// of a 2D conical product rule.
void QConical::conical_product_2D(unsigned int p)
{
  // Be sure the underlying rule object was built with the same dimension as the
  // rule we are about to construct.
  libmesh_assert (this->get_dim() == 2);
  
  QGauss  gauss1D(1,static_cast<Order>(_order+2*p));
  QJacobi jac1D(1,static_cast<Order>(_order+2*p),1,0);
	      
  // The Gauss rule needs to be scaled to [0,1]
  std::pair<Real, Real> old_range(-1.0L, 1.0L);
  std::pair<Real, Real> new_range( 0.0L, 1.0L);
  gauss1D.scale(old_range,
		new_range);

  // Now construct the points and weights for the conical product rule.

  // Both rules should have the same number of points.
  libmesh_assert(gauss1D.n_points() == jac1D.n_points());

  // Save the number of points as a convenient variable
  const unsigned int n_points = gauss1D.n_points();
  
  // Both rules should be between x=0 and x=1
  libmesh_assert(gauss1D.qp(0)(0) >= 0.0); libmesh_assert(gauss1D.qp(n_points-1)(0) <= 1.0);
  libmesh_assert(jac1D.qp(0)(0)   >= 0.0); libmesh_assert(jac1D.qp(n_points-1)(0) <= 1.0);

  // Resize the points and weights vectors
  _points.resize(n_points * n_points);
  _weights.resize(n_points * n_points);

  // Compute the conical product
  unsigned int gp = 0;
  for (unsigned int i=0; i<n_points; i++)
    for (unsigned int j=0; j<n_points; j++)
      {
	_points[gp](0) = jac1D.qp(j)(0);                          //s[j];
	_points[gp](1) = gauss1D.qp(i)(0) * (1.-jac1D.qp(j)(0)); //r[i]*(1.-s[j]);
	_weights[gp]   = gauss1D.w(i) * jac1D.w(j);              //A[i]*B[j];
	gp++;
      }
}




// Builds and scales a Gauss rule and a Jacobi rule.
// Then combines them to compute points and weights
// of a 3D conical product rule.
void QConical::conical_product_3D(unsigned int p)
{
  // Be sure the underlying rule object was built with the same dimension as the
  // rule we are about to construct.
  libmesh_assert (this->get_dim() == 3);
  
  QGauss  gauss1D(1,static_cast<Order>(_order+2*p));
  QJacobi jacA1D(1,static_cast<Order>(_order+2*p),1,0);
  QJacobi jacB1D(1,static_cast<Order>(_order+2*p),2,0);

  // The Gauss rule needs to be scaled to [0,1]
  std::pair<Real, Real> old_range(-1.0L, 1.0L);
  std::pair<Real, Real> new_range( 0.0L, 1.0L);
  gauss1D.scale(old_range,
		new_range);

  // Now construct the points and weights for the conical product rule.

  // All rules should have the same number of points
  libmesh_assert(gauss1D.n_points() == jacA1D.n_points());
  libmesh_assert(jacA1D.n_points()  == jacB1D.n_points());
  
  // Save the number of points as a convenient variable
  const unsigned int n_points = gauss1D.n_points();
  
  // All rules should be between x=0 and x=1
  libmesh_assert(gauss1D.qp(0)(0) >= 0.0); libmesh_assert(gauss1D.qp(n_points-1)(0) <= 1.0);
  libmesh_assert(jacA1D.qp(0)(0)  >= 0.0); libmesh_assert(jacA1D.qp(n_points-1)(0)  <= 1.0);
  libmesh_assert(jacB1D.qp(0)(0)  >= 0.0); libmesh_assert(jacB1D.qp(n_points-1)(0)  <= 1.0);

  // Resize the points and weights vectors
  _points.resize(n_points * n_points * n_points);
  _weights.resize(n_points * n_points * n_points);

  // Compute the conical product
  unsigned int gp = 0;
  for (unsigned int i=0; i<n_points; i++)
    for (unsigned int j=0; j<n_points; j++)
      for (unsigned int k=0; k<n_points; k++)
      {
	_points[gp](0) = jacB1D.qp(k)(0);                                                  //t[k];
	_points[gp](1) = jacA1D.qp(j)(0)  * (1.-jacB1D.qp(k)(0));                         //s[j]*(1.-t[k]);
	_points[gp](2) = gauss1D.qp(i)(0) * (1.-jacA1D.qp(j)(0)) * (1.-jacB1D.qp(k)(0)); //r[i]*(1.-s[j])*(1.-t[k]);
	_weights[gp]   = gauss1D.w(i)     * jacA1D.w(j)          * jacB1D.w(k);          //A[i]*B[j]*C[k];
	gp++;
      }
}
