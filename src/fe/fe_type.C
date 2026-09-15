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

// Local includes
#include "libmesh/fe_type.h"
#include "libmesh/quadrature_clough.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_gauss_lobatto.h"

// C++ Includes
#include <memory>


namespace libMesh
{

// ---------------------------------------
// FEType class members

std::unique_ptr<QBase>
FEType::default_quadrature_rule (const unsigned int dim,
                                 const int extraorder) const
{
  // Clough elements have at least piecewise cubic functions
  if (family == CLOUGH)
    {
      Order o = static_cast<Order>(std::max(static_cast<unsigned int>(this->default_quadrature_order()),
                                            static_cast<unsigned int>(7 + extraorder)));
      return std::make_unique<QClough>(dim, o);
    }

  if (family == SUBDIVISION)
    return std::make_unique<QGauss>(dim, static_cast<Order>(1 + extraorder));

  // The Gauss-Lobatto nodal bases interpolate at the points of the Gauss-Lobatto rule of
  // their own order, so that rule is the one they are collocated with: it makes their mass
  // matrix diagonal and asks nothing to evaluate their shape functions. A rule of n points
  // is exact through degree 2n-3, so the order that yields the p+1 points a basis of degree
  // p interpolates at is 2p-1. An extraorder carries past that, which buys accuracy on the
  // integrands a Gauss-Lobatto rule of p+1 points does not reach and gives up collocation.
  //
  // A p refined element keeps the pairing: QBase raises the order of a rule by twice the
  // p_level it is initialized with, and 2p-1 rises by two for each degree, so the rule on an
  // element of level l carries the p+l+1 points that the basis of degree p+l interpolates at.
  if (family == L2_LAGRANGE_GLL)
    {
      const int p = static_cast<int>(order.get_order());
      return std::make_unique<QGaussLobatto>(dim, static_cast<Order>(2*p - 1 + extraorder));
    }

  return std::make_unique<QGauss>(dim, static_cast<Order>(this->default_quadrature_order() + extraorder));
}


std::unique_ptr<QBase>
FEType::unweighted_quadrature_rule (const unsigned int dim,
                                    const int extraorder) const
{
  // Clough elements have at least piecewise cubic functions
  if (family == CLOUGH)
    {
      Order o = static_cast<Order>(std::max(static_cast<unsigned int>(this->unweighted_quadrature_order()),
                                            static_cast<unsigned int>(3 + extraorder)));
      return std::make_unique<QClough>(dim, o);
    }

  if (family == SUBDIVISION)
    return std::make_unique<QGauss>(dim, static_cast<Order>(1 + extraorder));

  return std::make_unique<QGauss>(dim, static_cast<Order>(this->unweighted_quadrature_order() + extraorder));
}

} // namespace libMesh
