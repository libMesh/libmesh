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
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_order.h"
#include "libmesh/int_range.h"
#include "libmesh/quadrature_gauss_lobatto.h"

// C++ includes
#include <algorithm>
#include <array>
#include <iterator>
#include <vector>

// The one-dimensional nodal basis on the Gauss-Lobatto points. The tensor product element
// types build their shape functions from it.
//
// The n points are -1, 1, and the roots of the derivative of the Legendre polynomial of
// degree n-1, which is why they are also called the Gauss-Lobatto-Legendre (GLL) points.
// The basis reads them from QGaussLobatto, so it and the rule with the same number of
// points use identical values. A basis and a quadrature rule that share their points are
// collocated: each shape function is one at its own point and zero at the others. The mass
// matrix assembled with a collocated rule is therefore diagonal for any mapping and element
// shape. That requires the points to match exactly, which is why the basis takes them from
// the rule.

namespace libMesh
{

/**
 * The one-dimensional nodal basis of degree n-1 on n Gauss-Lobatto points.
 */
struct GLLBasis1D
{
  /// The interpolation points in ascending order. The first is -1 and the last is 1.
  std::vector<Real> nodes;

  /// The barycentric weights. The i'th is the reciprocal of \f$\prod_{j \neq i}(\xi_i - \xi_j)\f$
  std::vector<Real> weights;

  /// \p deriv[k][i] is the derivative of the i'th shape function at the k'th point
  std::vector<std::vector<Real>> deriv;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /// \p second[k][i] is the second derivative of the i'th shape function at the k'th point
  std::vector<std::vector<Real>> second;
#endif
};


/**
 * The number of points in the largest Gauss-Lobatto rule that QGaussLobatto tabulates,
 * which it uses for orders FORTYSECOND and FORTYTHIRD. The highest available basis degree
 * is one less than this.
 */
constexpr unsigned int fe_lagrange_gll_max_points_1D = 23;


/**
 * \returns The one-dimensional nodal basis on \p n Gauss-Lobatto points.
 */
inline GLLBasis1D build_gll_basis_1D(const unsigned int n)
{
  GLLBasis1D b;

  // QGaussLobatto chooses its number of points from the order: a rule of n points is exact
  // through degree 2n-3
  QGaussLobatto rule(1, static_cast<Order>(2*n - 3));
  rule.init(EDGE2);
  libmesh_assert_equal_to(rule.n_points(), n);

  b.nodes.resize(n);
  for (const auto i : make_range(n))
    b.nodes[i] = rule.qp(i)(0);

  b.weights.resize(n);
  for (const auto i : make_range(n))
    {
      Real product = 1.;
      for (const auto j : make_range(n))
        if (j != i)
          product *= b.nodes[i] - b.nodes[j];
      b.weights[i] = 1. / product;
    }

  // Entry (k,i) is the derivative of shape function i at point k. Off the diagonal it is
  // (w_i/w_k)/(x_k - x_i). The shape functions sum to one, so each row sums to zero, and
  // the diagonal entry is minus the sum of the rest of its row. That also gives the
  // endpoint entries -p(p+1)/4 and p(p+1)/4, so they don't need a special case.
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

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  // Off the diagonal, entry (k,i) of the second differentiation matrix is
  // 2 D_ki (D_kk - 1/(x_k - x_i)), with D the first differentiation matrix. Each row sums to
  // zero for the same reason as in D, which gives the diagonal.
  b.second.assign(n, std::vector<Real>(n, 0.));
  for (const auto k : make_range(n))
    {
      Real row_sum = 0.;
      for (const auto i : make_range(n))
        if (i != k)
          {
            b.second[k][i] =
              2. * b.deriv[k][i] * (b.deriv[k][k] - 1. / (b.nodes[k] - b.nodes[i]));
            row_sum += b.second[k][i];
          }
      b.second[k][k] = -row_sum;
    }
#endif

  return b;
}

/**
 * \returns The one-dimensional basis of degree \p order.
 *
 * Every available basis is built in one static initialization on the first call. C++ makes
 * that initialization thread safe, so later calls read the bases without a lock. Each basis
 * of n points costs \f$O(n^2)\f$ operations, and building all of them takes well under a
 * millisecond.
 */
inline const GLLBasis1D & fe_lagrange_gll_1D_basis(const Order order)
{
  static const std::array<GLLBasis1D, fe_lagrange_gll_max_points_1D + 1> bases = []()
    {
      std::array<GLLBasis1D, fe_lagrange_gll_max_points_1D + 1> b;
      for (const auto n : make_range(2u, fe_lagrange_gll_max_points_1D + 1))
        b[n] = build_gll_basis_1D(n);
      return b;
    }();

  const unsigned int n = static_cast<unsigned int>(order) + 1;

  libmesh_error_msg_if(n < 2 || n > fe_lagrange_gll_max_points_1D,
                       "The Gauss-Lobatto nodal basis supports orders 1 through "
                       << fe_lagrange_gll_max_points_1D - 1 << ", but order " << order
                       << " was requested.");

  return bases[n];
}


/**
 * \returns The index of the point of \p b equal to \p xi, or \p invalid_uint if there is
 * none.
 *
 * The points are in ascending order, so a binary search finds the candidate. The comparison
 * is exact. The points of a collocated quadrature rule still match, because the basis and
 * the rule read the same values from QGaussLobatto.
 */
inline unsigned int fe_lagrange_gll_1D_point(const GLLBasis1D & b, const Real xi)
{
  const auto it = std::lower_bound(b.nodes.begin(), b.nodes.end(), xi);

  if (it != b.nodes.end() && *it == xi)
    return cast_int<unsigned int>(std::distance(b.nodes.begin(), it));

  return invalid_uint;
}

/**
 * \returns The value of the \p i'th one-dimensional Gauss-Lobatto nodal shape function of
 * degree \p order at \p xi.
 *
 * At a basis point the value is exactly one or zero, which a collocated quadrature rule
 * relies on. Elsewhere it is evaluated with the barycentric formula
 * \f$\ell_i(\xi) = \frac{w_i/(\xi - \xi_i)}{\sum_j w_j/(\xi - \xi_j)}\f$.
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
 * \returns The barycentric interpolant at \p xi of the values \p f[k] at the points of
 * \p b, for \p xi away from those points.
 *
 * This is \f$\sum_k \ell_k(\xi) f_k\f$ in the second barycentric form, which stays
 * accurate as \p xi approaches a point.
 */
template <typename Column>
inline Real fe_lagrange_gll_1D_interpolate(const GLLBasis1D & b,
                                           const Real xi,
                                           const Column & f)
{
  Real numerator = 0., denominator = 0.;

  for (const auto k : index_range(b.nodes))
    {
      const Real term = b.weights[k] / (xi - b.nodes[k]);

      numerator += term * f(k);
      denominator += term;
    }

  return numerator / denominator;
}

/**
 * \returns The derivative of the \p i'th one-dimensional Gauss-Lobatto nodal shape function
 * of degree \p order at \p xi.
 *
 * At a basis point it is read from the differentiation matrix. The derivative is a
 * polynomial of degree p-1, which the basis represents exactly, so elsewhere it is the
 * interpolant of column \p i of that matrix. Differentiating the barycentric formula
 * directly would instead cancel catastrophically as \p xi approaches a point.
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

  return fe_lagrange_gll_1D_interpolate(b, xi, [&b, i](const std::size_t m)
                                        { return b.deriv[m][i]; });
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
/**
 * \returns The second derivative of the \p i'th one-dimensional Gauss-Lobatto nodal shape
 * function of degree \p order at \p xi.
 *
 * At a basis point it is read from the second differentiation matrix. Elsewhere it is the
 * interpolant of column \p i of that matrix, for the same reason as the first derivative.
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

  return fe_lagrange_gll_1D_interpolate(b, xi, [&b, i](const std::size_t m)
                                        { return b.second[m][i]; });
}
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh

#endif // LIBMESH_FE_LAGRANGE_GLL_SHAPE_1D_H
