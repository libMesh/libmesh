// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/quadrature_nodal.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/quadrature_simpson.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

void QNodal::init_2D(const ElemType, unsigned int)
{
#if LIBMESH_DIM > 1

  switch (_type)
    {

    case QUAD4:
    case QUADSHELL4:
    case TRI3:
    case TRISHELL3:
      {
        QTrap rule(/*dim=*/2, /*ignored*/_order);
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());
        return;
      }

    case QUAD8:
    case QUADSHELL8:
      {
        // A rule with 8 points which is exact for linears, and
        // naturally produces a lumped approximation to the mass
        // matrix. The quadrature points are numbered the same way as
        // the reference element nodes.
        _points =
          {
            Point(-1,-1), Point(+1,-1), Point(+1,+1), Point(-1,+1),
            Point(0.,-1), Point(+1,0.), Point(0.,+1), Point(-1,0.)
          };

        // The weights for the Quad8 nodal quadrature rule are
        // obtained from the following specific steps. Other
        // "serendipity" type rules are obtained similarly.
        //
        // 1.) Due to the symmetry of the bi-unit square domain, we
        // first note that there are only two "classes" of node in the
        // Quad8: vertices (with associated weight wv) and edges (with
        // associated weight we).
        //
        // 2.) In order for such a nodal quadrature rule to be exact
        // for constants, the weights must sum to the area of the
        // reference element, and therefore we must have:
        // 4*wv + 4*we = 4 --> wv + we = 1
        //
        // 3.) Due to symmetry (once again), such a rule is then
        // automatically exact for the linear polynomials "x" and "y",
        // regardless of the values of wv and we.
        //
        // 4.) We therefore still have two unknowns and one equation.
        // Attempting to additionally make the nodal quadrature rule
        // exact for the quadratic polynomials "x^2" and "y^2" leads
        // to a rule with negative weights, namely:
        // wv = -1/3
        // we = 4/3
        // Note: these are the same values one obtains by integrating
        // the corresponding Quad8 Lagrange basis functions over the
        // reference element.
        //
        // Since the weights appear on the diagonal of the nodal
        // quadrature rule's approximation to the mass matrix, rules
        // with negative weights yield an indefinite mass matrix,
        // i.e. one with both positive and negative eigenvalues. Rules
        // with negative weights can also produce a negative integral
        // approximation to a strictly positive integrand, which may
        // be unacceptable in some situations.
        //
        // 5.) Requiring all positive weights therefore restricts the
        // nodal quadrature rule to only be exact for linears. But
        // there is still one free parameter remaining, and thus we
        // need some other criterion in order to complete the
        // description of the rule.
        //
        // Here, we have decided to choose the quadrature weights in
        // such a way that the resulting nodal quadrature mass matrix
        // approximation is "proportional" to the true mass matrix's
        // diagonal entries for the reference element. We therefore
        // pose the following constrained optimization problem:
        //
        //     { min_{wv, we} |diag(M) - C*diag(W)|^2
        // (O) {
        //     { subject to wv + we = 1
        //
        // where:
        // * M is the true mass matrix
        // * W is the nodal quadrature approximation to the mass
        //   matrix. In this particular case:
        //   diag(W) = [wv,wv,wv,wv,we,we,we,we]
        // * C = tr(M) / vol(E) is the ratio between the trace of the
        //   true mass matrix and the volume of the reference
        //   element. For all Lagrange finite elements, we have C<1.
        //
        // 6.) The optimization problem (O) is solved directly by
        // substituting the algebraic constraint into the functional
        // which is to be minimized, then setting the derivative with
        // respect to the remaining parameter equal to zero and
        // solving. In the Quad8 and Hex20 cases, there is only one
        // free parameter, while in the Prism15 case there are two
        // free parameters, so a 2x2 system of linear equations must
        // be solved.
        Real wv = Real(12) / 79;
        Real we = Real(67) / 79;

        _weights = {wv, wv, wv, wv, we, we, we, we};

        return;
      }

    case QUAD9:
    case TRI6:
      {
        QSimpson rule(/*dim=*/2, /*ignored*/_order);
        rule.init(_type, /*ignored*/_p_level);
        _points.swap (rule.get_points());
        _weights.swap(rule.get_weights());
        return;
      }

    case TRI7:
      {
        // We can't exactly represent cubics with only seven nodes,
        // but with w_i = integral(phi_i) for Lagrange shape functions
        // phi_i, we not only get exact integrals of every Lagrange
        // shape function, including the cubic bubble, we also get
        // exact integrals of the rest of P^3 too.
         _points =
          {
            Point(0.,0.), Point(+1,0.), Point(0.,+1), Point(.5,0.),
            Point(.5,.5), Point(0.,.5), Point(1/Real(3),1/Real(3))
          };

        Real wv = Real(1)/15;
        Real we = Real(1)/40;
        _weights = {wv, wv, wv, we, we, we, Real(9)/40};
        return;
      }

    default:
      libmesh_error_msg("Element type not supported!:" << Utility::enum_to_string(_type));
    }
#endif
}

} // namespace libMesh
