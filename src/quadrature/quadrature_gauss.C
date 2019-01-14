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
#include "libmesh/quadrature_gauss.h"
#include "libmesh/enum_quadrature_type.h"

namespace libMesh
{

// See the files:
// quadrature_gauss_1D.C
// quadrature_gauss_2D.C
// quadrature_gauss_3D.C
// for implementation of specific element types.


QuadratureType QGauss::type() const
{
  return QGAUSS;
}

void QGauss::keast_rule(const Real rule_data[][4],
                        const unsigned int n_pts)
{
  // Like the Dunavant rule, the input data should have 4 columns.  These columns
  // have the following format and implied permutations (w=weight).
  // {a, 0, 0, w} = 1-permutation  (a,a,a)
  // {a, b, 0, w} = 4-permutation  (a,b,b), (b,a,b), (b,b,a), (b,b,b)
  // {a, 0, b, w} = 6-permutation  (a,a,b), (a,b,b), (b,b,a), (b,a,b), (b,a,a), (a,b,a)
  // {a, b, c, w} = 12-permutation (a,a,b), (a,a,c), (b,a,a), (c,a,a), (a,b,a), (a,c,a)
  //                               (a,b,c), (a,c,b), (b,a,c), (b,c,a), (c,a,b), (c,b,a)

  // Always insert into the points & weights vector relative to the offset
  unsigned int offset=0;


  for (unsigned int p=0; p<n_pts; ++p)
    {

      // There must always be a non-zero entry to start the row
      libmesh_assert_not_equal_to (rule_data[p][0], static_cast<Real>(0.0));

      // A zero weight may imply you did not set up the raw data correctly
      libmesh_assert_not_equal_to (rule_data[p][3], static_cast<Real>(0.0));

      // What kind of point is this?
      // One non-zero entry in first 3 cols   ? 1-perm (centroid) point = 1
      // Two non-zero entries in first 3 cols ? 3-perm point            = 3
      // Three non-zero entries               ? 6-perm point            = 6
      unsigned int pointtype=1;

      if (rule_data[p][1] != static_cast<Real>(0.0))
        {
          if (rule_data[p][2] != static_cast<Real>(0.0))
            pointtype = 12;
          else
            pointtype = 4;
        }
      else
        {
          // The second entry is zero.  What about the third?
          if (rule_data[p][2] != static_cast<Real>(0.0))
            pointtype = 6;
        }


      switch (pointtype)
        {
        case 1:
          {
            // Be sure we have enough space to insert this point
            libmesh_assert_less (offset + 0, _points.size());

            const Real a = rule_data[p][0];

            // The point has only a single permutation (the centroid!)
            _points[offset  + 0] = Point(a,a,a);

            // The weight is always the last entry in the row.
            _weights[offset + 0] = rule_data[p][3];

            offset += pointtype;
            break;
          }

        case 4:
          {
            // Be sure we have enough space to insert these points
            libmesh_assert_less (offset + 3, _points.size());

            const Real a  = rule_data[p][0];
            const Real b  = rule_data[p][1];
            const Real wt = rule_data[p][3];

            // Here it's understood the second entry is to be used twice, and
            // thus there are three possible permutations.
            _points[offset + 0] = Point(a,b,b);
            _points[offset + 1] = Point(b,a,b);
            _points[offset + 2] = Point(b,b,a);
            _points[offset + 3] = Point(b,b,b);

            for (unsigned int j=0; j<pointtype; ++j)
              _weights[offset + j] = wt;

            offset += pointtype;
            break;
          }

        case 6:
          {
            // Be sure we have enough space to insert these points
            libmesh_assert_less (offset + 5, _points.size());

            const Real a  = rule_data[p][0];
            const Real b  = rule_data[p][2];
            const Real wt = rule_data[p][3];

            // Three individual entries with six permutations.
            _points[offset + 0] = Point(a,a,b);
            _points[offset + 1] = Point(a,b,b);
            _points[offset + 2] = Point(b,b,a);
            _points[offset + 3] = Point(b,a,b);
            _points[offset + 4] = Point(b,a,a);
            _points[offset + 5] = Point(a,b,a);

            for (unsigned int j=0; j<pointtype; ++j)
              _weights[offset + j] = wt;

            offset += pointtype;
            break;
          }


        case 12:
          {
            // Be sure we have enough space to insert these points
            libmesh_assert_less (offset + 11, _points.size());

            const Real a  = rule_data[p][0];
            const Real b  = rule_data[p][1];
            const Real c  = rule_data[p][2];
            const Real wt = rule_data[p][3];

            // Three individual entries with six permutations.
            _points[offset + 0] = Point(a,a,b);  _points[offset + 6]  = Point(a,b,c);
            _points[offset + 1] = Point(a,a,c); _points[offset + 7]  = Point(a,c,b);
            _points[offset + 2] = Point(b,a,a); _points[offset + 8]  = Point(b,a,c);
            _points[offset + 3] = Point(c,a,a); _points[offset + 9]  = Point(b,c,a);
            _points[offset + 4] = Point(a,b,a); _points[offset + 10] = Point(c,a,b);
            _points[offset + 5] = Point(a,c,a); _points[offset + 11] = Point(c,b,a);

            for (unsigned int j=0; j<pointtype; ++j)
              _weights[offset + j] = wt;

            offset += pointtype;
            break;
          }

        default:
          libmesh_error_msg("Don't know what to do with this many permutation points!");
        }

    }

}


// A number of different rules for triangles can be described by
// permutations of the following types of points:
// I:   "1"-permutation, (1/3,1/3)  (single point only)
// II:   3-permutation, (a,a,1-2a)
// III:  6-permutation, (a,b,1-a-b)
// The weights for a given set of permutations are all the same.
void QGauss::dunavant_rule2(const Real * wts,
                            const Real * a,
                            const Real * b,
                            const unsigned int * permutation_ids,
                            unsigned int n_wts)
{
  // Figure out how many total points by summing up the entries
  // in the permutation_ids array, and resize the _points and _weights
  // vectors appropriately.
  unsigned int total_pts = 0;
  for (unsigned int p=0; p<n_wts; ++p)
    total_pts += permutation_ids[p];

  // Resize point and weight vectors appropriately.
  _points.resize(total_pts);
  _weights.resize(total_pts);

  // Always insert into the points & weights vector relative to the offset
  unsigned int offset=0;

  for (unsigned int p=0; p<n_wts; ++p)
    {
      switch (permutation_ids[p])
        {
        case 1:
          {
            // The point has only a single permutation (the centroid!)
            // So we don't even need to look in the a or b arrays.
            _points [offset  + 0] = Point(Real(1)/3, Real(1)/3);
            _weights[offset + 0] = wts[p];

            offset += 1;
            break;
          }


        case 3:
          {
            // For this type of rule, don't need to look in the b array.
            _points[offset + 0] = Point(a[p],         a[p]);         // (a,a)
            _points[offset + 1] = Point(a[p],         1-2*a[p]); // (a,1-2a)
            _points[offset + 2] = Point(1-2*a[p], a[p]);         // (1-2a,a)

            for (unsigned int j=0; j<3; ++j)
              _weights[offset + j] = wts[p];

            offset += 3;
            break;
          }

        case 6:
          {
            // This type of point uses all 3 arrays...
            _points[offset + 0] = Point(a[p], b[p]);
            _points[offset + 1] = Point(b[p], a[p]);
            _points[offset + 2] = Point(a[p], 1-a[p]-b[p]);
            _points[offset + 3] = Point(1-a[p]-b[p], a[p]);
            _points[offset + 4] = Point(b[p], 1-a[p]-b[p]);
            _points[offset + 5] = Point(1-a[p]-b[p], b[p]);

            for (unsigned int j=0; j<6; ++j)
              _weights[offset + j] = wts[p];

            offset += 6;
            break;
          }

        default:
          libmesh_error_msg("Unknown permutation id: " << permutation_ids[p] << "!");
        }
    }

}


void QGauss::dunavant_rule(const Real rule_data[][4],
                           const unsigned int n_pts)
{
  // The input data array has 4 columns.  The first 3 are the permutation points.
  // The last column is the weights for a given set of permutation points.  A zero
  // in two of the first 3 columns implies the point is a 1-permutation (centroid).
  // A zero in one of the first 3 columns implies the point is a 3-permutation.
  // Otherwise each point is assumed to be a 6-permutation.

  // Always insert into the points & weights vector relative to the offset
  unsigned int offset=0;


  for (unsigned int p=0; p<n_pts; ++p)
    {

      // There must always be a non-zero entry to start the row
      libmesh_assert_not_equal_to ( rule_data[p][0], static_cast<Real>(0.0) );

      // A zero weight may imply you did not set up the raw data correctly
      libmesh_assert_not_equal_to ( rule_data[p][3], static_cast<Real>(0.0) );

      // What kind of point is this?
      // One non-zero entry in first 3 cols   ? 1-perm (centroid) point = 1
      // Two non-zero entries in first 3 cols ? 3-perm point            = 3
      // Three non-zero entries               ? 6-perm point            = 6
      unsigned int pointtype=1;

      if (rule_data[p][1] != static_cast<Real>(0.0))
        {
          if (rule_data[p][2] != static_cast<Real>(0.0))
            pointtype = 6;
          else
            pointtype = 3;
        }

      switch (pointtype)
        {
        case 1:
          {
            // Be sure we have enough space to insert this point
            libmesh_assert_less (offset + 0, _points.size());

            // The point has only a single permutation (the centroid!)
            _points[offset  + 0] = Point(rule_data[p][0], rule_data[p][0]);

            // The weight is always the last entry in the row.
            _weights[offset + 0] = rule_data[p][3];

            offset += 1;
            break;
          }

        case 3:
          {
            // Be sure we have enough space to insert these points
            libmesh_assert_less (offset + 2, _points.size());

            // Here it's understood the second entry is to be used twice, and
            // thus there are three possible permutations.
            _points[offset + 0] = Point(rule_data[p][0], rule_data[p][1]);
            _points[offset + 1] = Point(rule_data[p][1], rule_data[p][0]);
            _points[offset + 2] = Point(rule_data[p][1], rule_data[p][1]);

            for (unsigned int j=0; j<3; ++j)
              _weights[offset + j] = rule_data[p][3];

            offset += 3;
            break;
          }

        case 6:
          {
            // Be sure we have enough space to insert these points
            libmesh_assert_less (offset + 5, _points.size());

            // Three individual entries with six permutations.
            _points[offset + 0] = Point(rule_data[p][0], rule_data[p][1]);
            _points[offset + 1] = Point(rule_data[p][0], rule_data[p][2]);
            _points[offset + 2] = Point(rule_data[p][1], rule_data[p][0]);
            _points[offset + 3] = Point(rule_data[p][1], rule_data[p][2]);
            _points[offset + 4] = Point(rule_data[p][2], rule_data[p][0]);
            _points[offset + 5] = Point(rule_data[p][2], rule_data[p][1]);

            for (unsigned int j=0; j<6; ++j)
              _weights[offset + j] = rule_data[p][3];

            offset += 6;
            break;
          }

        default:
          libmesh_error_msg("Don't know what to do with this many permutation points!");
        }
    }
}

} // namespace libMesh
