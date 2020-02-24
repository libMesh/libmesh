// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_JACOBI_POLYNOMIALS_H
#define LIBMESH_JACOBI_POLYNOMIALS_H

// libMesh includes
#include "libmesh/libmesh_common.h" // Real

// C++ includes
#include <utility> // std::swap

namespace libMesh
{

namespace JacobiPolynomials
{

/**
 * The Jacobi polynomial value and derivative formulas are based on
 * the corresponding Wikipedia article [0]. Note that we have shifted
 * the indices used in the recurrence relation, otherwise this is
 * identical to the recurrence relation given in the article. When
 * alpha = beta = 0, the Jacobi polynomials reduce to the Legendre
 * polynomials.
 *
 * [0]: https://en.wikipedia.org/wiki/Jacobi_polynomials
 */
inline
Real value(unsigned n, unsigned alpha, unsigned beta, Real x)
{
  if (n == 0)
    return 1.;

  // Compute constants independent of loop index.
  unsigned int ab = alpha + beta;
  unsigned int a2 = alpha * alpha;
  unsigned int b2 = beta * beta;

  Real p0 = 1;
  Real p1 = (alpha + 1) + (ab + 2) * 0.5 * (x - 1);

  unsigned int i = 1;
  while (i < n)
    {
      // Note: we swap before updating p1, so p0 and p1 appear in
      // opposite positions than is usual in the update formula.
      std::swap(p0, p1);
      p1 = (((2*i + ab + 1) *
        ((2*i + ab + 2) * (2*i + ab) * x + a2 - b2)) * p0
        - 2 * (i + alpha) * (i + beta) * (2*i + ab + 2) * p1) /
       (2 * (i + 1) * (i + 1 + ab) * (2*i + ab));

      ++i;
    }
  return p1;
}

inline
Real deriv(unsigned n, unsigned alpha, unsigned beta, Real x)
{
  // Call value() with elevated (alpha, beta) and decremented n.
  return n == 0 ? 0 : 0.5 * (1 + alpha + beta + n) * value(n-1, alpha+1, beta+1, x);
}

} // namespace JacobiPolynomials
} // namespace libMesh

#endif // LIBMESH_JACOBI_POLYNOMIALS_H
