// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_FE_LAGRANGE_GLL_SHAPE_1D_H
#define LIBMESH_FE_LAGRANGE_GLL_SHAPE_1D_H

// Local includes
#include "libmesh/enum_order.h"
#include "libmesh/int_range.h"
#include "libmesh/quadrature_gauss_lobatto.h"

// C++ includes
#include <array>
#include <vector>

// The one-dimensional nodal basis on the Gauss-Lobatto points, which the tensor product
// elements build their shape functions from.
//
// The interpolation points of this basis are the points of the Gauss-Lobatto rule of the
// same count, read from QGaussLobatto so that the two carry the same values. A basis and a
// quadrature rule that share their points are collocated: every shape function is one at
// its own point and zero at the others, so the mass matrix a collocated rule assembles is
// diagonal, for any mapping and any element shape. That property is exact only while the
// two agree bit for bit, which is why the points are shared rather than recomputed here.

namespace libMesh
{

/**
 * The one-dimensional nodal basis on \p n Gauss-Lobatto points, of degree n-1.
 */
struct GLLBasis1D
{
  /// The interpolation points, in ascending order, the first and last being -1 and 1
  std::vector<Real> nodes;

  /// The barycentric weights, the i'th being the reciprocal of \f$\prod_{j \neq i}(\xi_i - \xi_j)\f$
  std::vector<Real> weights;

  /// \p deriv[k][i] is the derivative of the i'th shape function at the k'th point
  std::vector<std::vector<Real>> deriv;

  /// \p second[k][i] is the second derivative of the i'th shape function at the k'th point
  std::vector<std::vector<Real>> second;
};


/**
 * \returns The one-dimensional nodal basis on \p n Gauss-Lobatto points.
 */
inline GLLBasis1D build_gll_basis_1D(const unsigned int n)
{
  GLLBasis1D b;

  b.nodes = QGaussLobatto::points_1D(n);

  b.weights.resize(n);
  for (const auto i : make_range(n))
    {
      Real product = 1.;
      for (const auto j : make_range(n))
        if (j != i)
          product *= b.nodes[i] - b.nodes[j];
      b.weights[i] = 1. / product;
    }

  // Away from its own point the derivative of a shape function is the ratio of two
  // barycentric weights over the distance between the points. At its own point it is minus
  // the sum of the others in its row, since the shape functions sum to one everywhere and
  // so their derivatives sum to zero. That identity gives the two endpoint entries,
  // -p(p+1)/4 and p(p+1)/4, without a case of their own.
  b.deriv.assign(n, std::vector<Real>(n, 0.));
  for (const auto k : make_range(n))
    {
      Real row_sum = 0.;
      for (const auto i : make_range(n))
        if (i != k)
          {
            b.deriv[k][i] = (b.weights[i] / b.weights[k]) / (b.nodes[k] - b.nodes[i]);
            row_sum += b.deriv[k][i];
          }
      b.deriv[k][k] = -row_sum;
    }

  // The derivative of a polynomial of degree p is one of degree p-1, which this basis
  // interpolates exactly, so applying the first derivative twice gives the second.
  b.second.assign(n, std::vector<Real>(n, 0.));
  for (const auto k : make_range(n))
    for (const auto i : make_range(n))
      {
        Real sum = 0.;
        for (const auto m : make_range(n))
          sum += b.deriv[k][m] * b.deriv[m][i];
        b.second[k][i] = sum;
      }

  return b;
}

/**
 * \returns The one-dimensional basis of degree \p order, built once and shared by all
 * callers.
 *
 * The bases are built together on first use, which is a few thousand operations, so that
 * no lock is needed to reach one.
 */
inline const GLLBasis1D & fe_lagrange_gll_1D_basis(const Order order)
{
  static const std::array<GLLBasis1D, QGaussLobatto::max_points_1D + 1> bases = []()
    {
      std::array<GLLBasis1D, QGaussLobatto::max_points_1D + 1> b;
      for (unsigned int n = 2; n <= QGaussLobatto::max_points_1D; ++n)
        b[n] = build_gll_basis_1D(n);
      return b;
    }();

  const unsigned int n = static_cast<unsigned int>(order) + 1;

  libmesh_error_msg_if(n < 2 || n > QGaussLobatto::max_points_1D,
                       "No Gauss-Lobatto nodal basis of order " << order << " is available; "
                       "the tabulated Gauss-Lobatto points run to "
                       << QGaussLobatto::max_points_1D << ", so the orders run from 1 to "
                       << QGaussLobatto::max_points_1D - 1 << ".");

  return bases[n];
}


/**
 * \returns The index of the point of \p b that sits at \p xi, or \p invalid_uint when \p xi
 * is not one of them.
 *
 * The comparison is exact, and it holds for every point of a collocated quadrature rule
 * because the basis and the rule read the same values from QGaussLobatto.
 */
inline unsigned int fe_lagrange_gll_1D_point(const GLLBasis1D & b, const Real xi)
{
  for (const auto i : index_range(b.nodes))
    if (b.nodes[i] == xi)
      return cast_int<unsigned int>(i);

  return invalid_uint;
}

/**
 * \returns The value of the \p i'th one-dimensional Gauss-Lobatto nodal shape function of
 * degree \p order at \p xi.
 *
 * At one of the basis's own points the value is one or zero, which is the property a
 * collocated quadrature rule rests on. Elsewhere the barycentric form gives it, whose
 * numerator and denominator are the terms \f$w_j/(\xi - \xi_j)\f$ that the shape functions
 * share.
 */
inline Real fe_lagrange_gll_1D_shape(const Order order,
                                     const unsigned int i,
                                     const Real xi)
{
  const GLLBasis1D & b = fe_lagrange_gll_1D_basis(order);

  libmesh_assert_less (i, b.nodes.size());

  const unsigned int k = fe_lagrange_gll_1D_point(b, xi);
  if (k != invalid_uint)
    return (i == k) ? 1. : 0.;

  Real numerator = 0., denominator = 0.;

  for (const auto j : index_range(b.nodes))
    {
      const Real term = b.weights[j] / (xi - b.nodes[j]);

      denominator += term;

      if (j == i)
        numerator = term;
    }

  return numerator / denominator;
}

/**
 * \returns The derivative of the \p i'th one-dimensional Gauss-Lobatto nodal shape function
 * of degree \p order at \p xi.
 *
 * At one of the basis's own points this is a row of the differentiation matrix, which the
 * basis carries. Elsewhere it follows from differentiating the barycentric form, as
 * \f$\ell_i(\xi)\left[S_2/S_1 - 1/(\xi - \xi_i)\right]\f$ with
 * \f$S_m = \sum_j w_j/(\xi-\xi_j)^m\f$.
 */
inline Real fe_lagrange_gll_1D_shape_deriv(const Order order,
                                           const unsigned int i,
                                           const unsigned int libmesh_dbg_var(j),
                                           const Real xi)
{
  // only d()/dxi in 1D!
  libmesh_assert_equal_to (j, 0);

  const GLLBasis1D & b = fe_lagrange_gll_1D_basis(order);

  libmesh_assert_less (i, b.nodes.size());

  const unsigned int k = fe_lagrange_gll_1D_point(b, xi);
  if (k != invalid_uint)
    return b.deriv[k][i];

  Real s1 = 0., s2 = 0., numerator = 0.;

  for (const auto m : index_range(b.nodes))
    {
      const Real term = b.weights[m] / (xi - b.nodes[m]);

      s1 += term;
      s2 += term / (xi - b.nodes[m]);

      if (m == i)
        numerator = term;
    }

  const Real value = numerator / s1;

  return value * (s2 / s1 - 1. / (xi - b.nodes[i]));
}

/**
 * \returns The second derivative of the \p i'th one-dimensional Gauss-Lobatto nodal shape
 * function of degree \p order at \p xi.
 *
 * At one of the basis's own points this is a row of the second differentiation matrix.
 * Elsewhere it follows from differentiating the barycentric form twice, as
 * \f$2\ell_i(\xi)\left[1/d^2 - S_2/(d S_1) - S_3/S_1 + S_2^2/S_1^2\right]\f$ with
 * \f$d = \xi - \xi_i\f$ and \f$S_m = \sum_j w_j/(\xi-\xi_j)^m\f$.
 */
inline Real fe_lagrange_gll_1D_shape_second_deriv(const Order order,
                                                  const unsigned int i,
                                                  const unsigned int libmesh_dbg_var(j),
                                                  const Real xi)
{
  // only d^2()/dxi^2 in 1D!
  libmesh_assert_equal_to (j, 0);

  const GLLBasis1D & b = fe_lagrange_gll_1D_basis(order);

  libmesh_assert_less (i, b.nodes.size());

  const unsigned int k = fe_lagrange_gll_1D_point(b, xi);
  if (k != invalid_uint)
    return b.second[k][i];

  Real s1 = 0., s2 = 0., s3 = 0., numerator = 0.;

  for (const auto m : index_range(b.nodes))
    {
      const Real d = xi - b.nodes[m];
      const Real term = b.weights[m] / d;

      s1 += term;
      s2 += term / d;
      s3 += term / (d * d);

      if (m == i)
        numerator = term;
    }

  const Real value = numerator / s1;
  const Real d = xi - b.nodes[i];

  return 2. * value * (1. / (d * d) - s2 / (d * s1) - s3 / s1 + (s2 * s2) / (s1 * s1));
}

} // namespace libMesh

#endif // LIBMESH_FE_LAGRANGE_GLL_SHAPE_1D_H
