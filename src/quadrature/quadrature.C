// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/elem.h"
#include "libmesh/quadrature.h"

namespace libMesh
{

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
      libmesh_error_msg("Invalid dimension _dim = " << _dim);
    }
}



void QBase::init (const Elem & elem,
                  const std::vector<Real> & /* vertex_distance_func */,
                  unsigned int p_level)
{
  // dispatch generic implementation
  this->init(elem.type(), p_level);
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
  libmesh_assert_equal_to (_dim, 1);

  Real
    h_new = new_range.second - new_range.first,
    h_old = old_range.second - old_range.first;

  // Make sure that we have sane ranges
  libmesh_assert_greater (h_new, 0.);
  libmesh_assert_greater (h_old, 0.);

  // Make sure there are some points
  libmesh_assert_greater (_points.size(), 0);

  // Compute the scale factor
  Real scfact = h_new/h_old;

  // We're mapping from old_range -> new_range
  for (std::size_t i=0; i<_points.size(); i++)
    {
      _points[i](0) = new_range.first +
        (_points[i](0) - old_range.first) * scfact;

      // Scale the weights
      _weights[i] *= scfact;
    }
}




void QBase::tensor_product_quad(const QBase & q1D)
{

  const unsigned int np = q1D.n_points();

  _points.resize(np * np);

  _weights.resize(np * np);

  unsigned int q=0;

  for (unsigned int j=0; j<np; j++)
    for (unsigned int i=0; i<np; i++)
      {
        _points[q](0) = q1D.qp(i)(0);
        _points[q](1) = q1D.qp(j)(0);

        _weights[q] = q1D.w(i)*q1D.w(j);

        q++;
      }
}





void QBase::tensor_product_hex(const QBase & q1D)
{
  const unsigned int np = q1D.n_points();

  _points.resize(np * np * np);

  _weights.resize(np * np * np);

  unsigned int q=0;

  for (unsigned int k=0; k<np; k++)
    for (unsigned int j=0; j<np; j++)
      for (unsigned int i=0; i<np; i++)
        {
          _points[q](0) = q1D.qp(i)(0);
          _points[q](1) = q1D.qp(j)(0);
          _points[q](2) = q1D.qp(k)(0);

          _weights[q] = q1D.w(i) * q1D.w(j) * q1D.w(k);

          q++;
        }
}




void QBase::tensor_product_prism(const QBase & q1D, const QBase & q2D)
{
  const unsigned int n_points1D = q1D.n_points();
  const unsigned int n_points2D = q2D.n_points();

  _points.resize  (n_points1D * n_points2D);
  _weights.resize (n_points1D * n_points2D);

  unsigned int q=0;

  for (unsigned int j=0; j<n_points1D; j++)
    for (unsigned int i=0; i<n_points2D; i++)
      {
        _points[q](0) = q2D.qp(i)(0);
        _points[q](1) = q2D.qp(i)(1);
        _points[q](2) = q1D.qp(j)(0);

        _weights[q] = q2D.w(i) * q1D.w(j);

        q++;
      }

}




std::ostream & operator << (std::ostream & os, const QBase & q)
{
  q.print_info(os);
  return os;
}

} // namespace libMesh
