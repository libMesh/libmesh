// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// libMesh includes
#include "libmesh/quadrature_monomial.h"
#include "libmesh/enum_quadrature_type.h"

namespace libMesh
{


// See the files:
// quadrature_monomial_2D.C
// quadrature_monomial_3D.C
// for implementation of specific element types.

QuadratureType QMonomial::type() const
{
  return QMONOMIAL;
}

void QMonomial::wissmann_rule(const Real rule_data[][3],
                              const unsigned int n_pts)
{
  for (unsigned int i=0, c=0; i<n_pts; ++i)
    {
      _points[c]  = Point( rule_data[i][0], rule_data[i][1] );
      _weights[c++] = rule_data[i][2];

      // This may be an (x1,x2) -> (-x1,x2) point, in which case
      // we will also generate the mirror point using the same weight.
      if (rule_data[i][0] != static_cast<Real>(0.0))
        {
          _points[c]  = Point( -rule_data[i][0], rule_data[i][1] );
          _weights[c++] = rule_data[i][2];
        }
    }
}



void QMonomial::stroud_rule(const Real rule_data[][3],
                            const unsigned int * rule_symmetry,
                            const unsigned int n_pts)
{
  for (unsigned int i=0, c=0; i<n_pts; ++i)
    {
      const Real
        x=rule_data[i][0],
        y=rule_data[i][1],
        wt=rule_data[i][2];

      switch(rule_symmetry[i])
        {
        case 0: // Single point (no symmetry)
          {
            _points[c]  = Point( x, y);
            _weights[c++] = wt;

            break;
          }
        case 1: // Fully-symmetric (x,y)
          {
            _points[c]    = Point( x, y);
            _weights[c++] = wt;

            _points[c]    = Point(-x, y);
            _weights[c++] = wt;

            _points[c]    = Point( x,-y);
            _weights[c++] = wt;

            _points[c]    = Point(-x,-y);
            _weights[c++] = wt;

            _points[c]    = Point( y, x);
            _weights[c++] = wt;

            _points[c]    = Point(-y, x);
            _weights[c++] = wt;

            _points[c]    = Point( y,-x);
            _weights[c++] = wt;

            _points[c]    = Point(-y,-x);
            _weights[c++] = wt;

            break;
          }
        case 2: // Fully-symmetric (x,x)
          {
            _points[c]    = Point( x, x);
            _weights[c++] = wt;

            _points[c]    = Point(-x, x);
            _weights[c++] = wt;

            _points[c]    = Point( x,-x);
            _weights[c++] = wt;

            _points[c]    = Point(-x,-x);
            _weights[c++] = wt;

            break;
          }
        case 3: // Fully-symmetric (x,0)
          {
            libmesh_assert_equal_to (y, 0.0);

            _points[c]    = Point( x,0.);
            _weights[c++] = wt;

            _points[c]    = Point(-x,0.);
            _weights[c++] = wt;

            _points[c]    = Point(0., x);
            _weights[c++] = wt;

            _points[c]    = Point(0.,-x);
            _weights[c++] = wt;

            break;
          }
        case 4: // Rotational invariant
          {
            _points[c]    = Point( x, y);
            _weights[c++] = wt;

            _points[c]    = Point(-x,-y);
            _weights[c++] = wt;

            _points[c]    = Point(-y, x);
            _weights[c++] = wt;

            _points[c]    = Point( y,-x);
            _weights[c++] = wt;

            break;
          }
        case 5: // Partial symmetry (Wissman's rules)
          {
            libmesh_assert_not_equal_to (x, 0.0);

            _points[c]    = Point( x, y);
            _weights[c++] = wt;

            _points[c]    = Point(-x, y);
            _weights[c++] = wt;

            break;
          }
        case 6: // Rectangular symmetry
          {
            _points[c]    = Point( x, y);
            _weights[c++] = wt;

            _points[c]    = Point(-x, y);
            _weights[c++] = wt;

            _points[c]    = Point(-x,-y);
            _weights[c++] = wt;

            _points[c]    = Point( x,-y);
            _weights[c++] = wt;

            break;
          }
        case 7: // Central symmetry
          {
            libmesh_assert_equal_to (x, 0.0);
            libmesh_assert_not_equal_to (y, 0.0);

            _points[c]    = Point(0., y);
            _weights[c++] = wt;

            _points[c]    = Point(0.,-y);
            _weights[c++] = wt;

            break;
          }
        default:
          libmesh_error_msg("Unknown symmetry!");
        } // end switch(rule_symmetry[i])
    }
}




void QMonomial::kim_rule(const Real rule_data[][4],
                         const unsigned int * rule_id,
                         const unsigned int n_pts)
{
  for (unsigned int i=0, c=0; i<n_pts; ++i)
    {
      const Real
        x=rule_data[i][0],
        y=rule_data[i][1],
        z=rule_data[i][2],
        wt=rule_data[i][3];

      switch(rule_id[i])
        {
        case 0: // (0,0,0) 1 permutation
          {
            _points[c]  = Point( x, y, z);    _weights[c++] = wt;

            break;
          }
        case 1: //  (x,0,0) 6 permutations
          {
            _points[c] = Point( x, 0., 0.);    _weights[c++] = wt;
            _points[c] = Point(0., -x, 0.);    _weights[c++] = wt;
            _points[c] = Point(-x, 0., 0.);    _weights[c++] = wt;
            _points[c] = Point(0.,  x, 0.);    _weights[c++] = wt;
            _points[c] = Point(0., 0., -x);    _weights[c++] = wt;
            _points[c] = Point(0., 0.,  x);    _weights[c++] = wt;

            break;
          }
        case 2: // (x,x,0) 12 permutations
          {
            _points[c] = Point( x,  x, 0.);    _weights[c++] = wt;
            _points[c] = Point( x, -x, 0.);    _weights[c++] = wt;
            _points[c] = Point(-x, -x, 0.);    _weights[c++] = wt;
            _points[c] = Point(-x,  x, 0.);    _weights[c++] = wt;
            _points[c] = Point( x, 0., -x);    _weights[c++] = wt;
            _points[c] = Point( x, 0.,  x);    _weights[c++] = wt;
            _points[c] = Point(0.,  x, -x);    _weights[c++] = wt;
            _points[c] = Point(0.,  x,  x);    _weights[c++] = wt;
            _points[c] = Point(0., -x, -x);    _weights[c++] = wt;
            _points[c] = Point(-x, 0., -x);    _weights[c++] = wt;
            _points[c] = Point(0., -x,  x);    _weights[c++] = wt;
            _points[c] = Point(-x, 0.,  x);    _weights[c++] = wt;

            break;
          }
        case 3: // (x,y,0) 24 permutations
          {
            _points[c] = Point( x,  y, 0.);    _weights[c++] = wt;
            _points[c] = Point( y, -x, 0.);    _weights[c++] = wt;
            _points[c] = Point(-x, -y, 0.);    _weights[c++] = wt;
            _points[c] = Point(-y,  x, 0.);    _weights[c++] = wt;
            _points[c] = Point( x, 0., -y);    _weights[c++] = wt;
            _points[c] = Point( x, -y, 0.);    _weights[c++] = wt;
            _points[c] = Point( x, 0.,  y);    _weights[c++] = wt;
            _points[c] = Point(0.,  y, -x);    _weights[c++] = wt;
            _points[c] = Point(-x,  y, 0.);    _weights[c++] = wt;
            _points[c] = Point(0.,  y,  x);    _weights[c++] = wt;
            _points[c] = Point( y, 0., -x);    _weights[c++] = wt;
            _points[c] = Point(0., -y, -x);    _weights[c++] = wt;
            _points[c] = Point(-y, 0., -x);    _weights[c++] = wt;
            _points[c] = Point( y,  x, 0.);    _weights[c++] = wt;
            _points[c] = Point(-y, -x, 0.);    _weights[c++] = wt;
            _points[c] = Point( y, 0.,  x);    _weights[c++] = wt;
            _points[c] = Point(0., -y,  x);    _weights[c++] = wt;
            _points[c] = Point(-y, 0.,  x);    _weights[c++] = wt;
            _points[c] = Point(-x, 0.,  y);    _weights[c++] = wt;
            _points[c] = Point(0., -x, -y);    _weights[c++] = wt;
            _points[c] = Point(0., -x,  y);    _weights[c++] = wt;
            _points[c] = Point(-x, 0., -y);    _weights[c++] = wt;
            _points[c] = Point(0.,  x,  y);    _weights[c++] = wt;
            _points[c] = Point(0.,  x, -y);    _weights[c++] = wt;

            break;
          }
        case 4: // (x,x,x) 8 permutations
          {
            _points[c] = Point( x,  x,  x);    _weights[c++] = wt;
            _points[c] = Point( x, -x,  x);    _weights[c++] = wt;
            _points[c] = Point(-x, -x,  x);    _weights[c++] = wt;
            _points[c] = Point(-x,  x,  x);    _weights[c++] = wt;
            _points[c] = Point( x,  x, -x);    _weights[c++] = wt;
            _points[c] = Point( x, -x, -x);    _weights[c++] = wt;
            _points[c] = Point(-x,  x, -x);    _weights[c++] = wt;
            _points[c] = Point(-x, -x, -x);    _weights[c++] = wt;

            break;
          }
        case 5: // (x,x,z) 24 permutations
          {
            _points[c] = Point( x,  x,  z);    _weights[c++] = wt;
            _points[c] = Point( x, -x,  z);    _weights[c++] = wt;
            _points[c] = Point(-x, -x,  z);    _weights[c++] = wt;
            _points[c] = Point(-x,  x,  z);    _weights[c++] = wt;
            _points[c] = Point( x,  z, -x);    _weights[c++] = wt;
            _points[c] = Point( x, -x, -z);    _weights[c++] = wt;
            _points[c] = Point( x, -z,  x);    _weights[c++] = wt;
            _points[c] = Point( z,  x, -x);    _weights[c++] = wt;
            _points[c] = Point(-x,  x, -z);    _weights[c++] = wt;
            _points[c] = Point(-z,  x,  x);    _weights[c++] = wt;
            _points[c] = Point( x, -z, -x);    _weights[c++] = wt;
            _points[c] = Point(-z, -x, -x);    _weights[c++] = wt;
            _points[c] = Point(-x,  z, -x);    _weights[c++] = wt;
            _points[c] = Point( x,  x, -z);    _weights[c++] = wt;
            _points[c] = Point(-x, -x, -z);    _weights[c++] = wt;
            _points[c] = Point( x,  z,  x);    _weights[c++] = wt;
            _points[c] = Point( z, -x,  x);    _weights[c++] = wt;
            _points[c] = Point(-x, -z,  x);    _weights[c++] = wt;
            _points[c] = Point(-x,  z,  x);    _weights[c++] = wt;
            _points[c] = Point( z, -x, -x);    _weights[c++] = wt;
            _points[c] = Point(-z, -x,  x);    _weights[c++] = wt;
            _points[c] = Point(-x, -z, -x);    _weights[c++] = wt;
            _points[c] = Point( z,  x,  x);    _weights[c++] = wt;
            _points[c] = Point(-z,  x, -x);    _weights[c++] = wt;

            break;
          }
        case 6: // (x,y,z) 24 permutations
          {
            _points[c] = Point( x,  y,  z);    _weights[c++] = wt;
            _points[c] = Point( y, -x,  z);    _weights[c++] = wt;
            _points[c] = Point(-x, -y,  z);    _weights[c++] = wt;
            _points[c] = Point(-y,  x,  z);    _weights[c++] = wt;
            _points[c] = Point( x,  z, -y);    _weights[c++] = wt;
            _points[c] = Point( x, -y, -z);    _weights[c++] = wt;
            _points[c] = Point( x, -z,  y);    _weights[c++] = wt;
            _points[c] = Point( z,  y, -x);    _weights[c++] = wt;
            _points[c] = Point(-x,  y, -z);    _weights[c++] = wt;
            _points[c] = Point(-z,  y,  x);    _weights[c++] = wt;
            _points[c] = Point( y, -z, -x);    _weights[c++] = wt;
            _points[c] = Point(-z, -y, -x);    _weights[c++] = wt;
            _points[c] = Point(-y,  z, -x);    _weights[c++] = wt;
            _points[c] = Point( y,  x, -z);    _weights[c++] = wt;
            _points[c] = Point(-y, -x, -z);    _weights[c++] = wt;
            _points[c] = Point( y,  z,  x);    _weights[c++] = wt;
            _points[c] = Point( z, -y,  x);    _weights[c++] = wt;
            _points[c] = Point(-y, -z,  x);    _weights[c++] = wt;
            _points[c] = Point(-x,  z,  y);    _weights[c++] = wt;
            _points[c] = Point( z, -x, -y);    _weights[c++] = wt;
            _points[c] = Point(-z, -x,  y);    _weights[c++] = wt;
            _points[c] = Point(-x, -z, -y);    _weights[c++] = wt;
            _points[c] = Point( z,  x,  y);    _weights[c++] = wt;
            _points[c] = Point(-z,  x, -y);    _weights[c++] = wt;

            break;
          }
        default:
          libmesh_error_msg("Unknown rule ID: " << rule_id[i] << "!");
        } // end switch(rule_id[i])
    }
}

} // namespace libMesh
