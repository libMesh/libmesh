// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/edge_edge4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

// Edge4 class static member initializations
const int Edge4::num_nodes;
const int Edge4::num_sides;
const int Edge4::num_edges;
const int Edge4::num_children;
const int Edge4::nodes_per_side;
const int Edge4::nodes_per_edge;

#ifdef LIBMESH_ENABLE_AMR

const Real Edge4::_embedding_matrix[Edge4::num_children][Edge4::num_nodes][Edge4::num_nodes] =
  {
    // embedding matrix for child 0

    {
      // 0       1       2        3        // Shape function index
      {1.0,      0.0,    0.0,     0.0},    // left,  xi = -1
      {-0.0625, -0.0625, 0.5625,  0.5625}, // right, xi = 0
      {0.3125,   0.0625, 0.9375, -0.3125}, // middle left,  xi = -2/3
      {0.0,      0.0,    1.0,     0.0}     // middle right, xi = -1/3
    },

    // embedding matrix for child 1
    {
      // 0       1        2        3        // Shape function index
      {-0.0625, -0.0625,  0.5625,  0.5625}, // left,  xi = 0
      {0.0,      1.0,     0.0,     0.0},    // right, xi = 1
      {0.0,      0.0,     0.0,     1.0},    // middle left,  xi = 1/3
      {0.0625,   0.3125, -0.3125,  0.9375}  // middle right, xi = 2/3
    }
  };

#endif

bool Edge4::is_vertex(const unsigned int i) const
{
  return (i==0) || (i==1);
}

bool Edge4::is_edge(const unsigned int i) const
{
  return (i==2) || (i==3);
}

bool Edge4::is_face(const unsigned int ) const
{
  return false;
}

bool Edge4::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, 2);
  libmesh_assert_less (n, Edge4::num_nodes);
  return (s == n);
}

bool Edge4::is_node_on_edge(const unsigned int,
                            const unsigned int libmesh_dbg_var(e)) const
{
  libmesh_assert_equal_to (e, 0);
  return true;
}



bool Edge4::has_affine_map() const
{
  Point v = this->point(1) - this->point(0);
  if (!v.relative_fuzzy_equals
      ((this->point(2) - this->point(0))*3, affine_tol))
    return false;
  if (!v.relative_fuzzy_equals
      ((this->point(3) - this->point(0))*1.5, affine_tol))
    return false;
  return true;
}



bool Edge4::has_invertible_map(Real tol) const
{
  // At the moment this only makes sense for Lagrange elements
  libmesh_assert_equal_to(this->mapping_type(), LAGRANGE_MAP);

  // dx/dxi = a*xi^2 + b*xi + c,
  // where a, b, and c are vector quantities that depend on the
  // nodal positions as follows:
  const Point & x0 = this->point(0);
  const Point & x1 = this->point(1);
  const Point & x2 = this->point(2);
  const Point & x3 = this->point(3);

  Point a = 27./16 * (x1 - x0 + 3.*(x2 - x3));
  Point b = 18./16 * (x0 + x1 - x2 - x3);
  Point c = 1./16  * (x0 - x1 + 27.*(x3 - x2));

  // Normalized midpoint value of Jacobian. If c==0, then dx/dxi(0) =
  // 0 and the Jacobian is either zero on the whole element, or
  // changes sign within the element. Either way, the element is not
  // invertible. Note: use <= tol, so that if someone passes 0 for tol
  // it still works.
  Real c_norm = c.norm();
  if (c_norm <= tol)
    return false;

  // Coefficients of the quadratic scalar function
  // j(xi) := dot(c.unit(), dx/dxi)
  //        = alpha*xi^2 + beta*xi + gamma
  Real alpha = (a * c) / c_norm;
  Real beta = (b * c) / c_norm;
  Real gamma = c_norm;

  // If alpha and beta are both (approximately) zero, then the
  // Jacobian is actually constant but it's not zero (as we already
  // checked) so the element is invertible!
  if ((std::abs(alpha) <= tol) && (std::abs(beta) <= tol))
    return true;

  // If alpha is approximately zero, but beta is not, then j(xi) is
  // actually linear and the Jacobian changes sign at the point xi =
  // -gamma/beta, so we need to check whether that's in the
  // element.
  if (std::abs(alpha) <= tol)
    {
      // Debugging:
      // libMesh::out << "alpha = " << alpha << ", std::abs(alpha) <= tol" << std::endl;

      Real xi_0 = -gamma / beta;
      return ((xi_0 < -1.) || (xi_0 > 1.));
    }

  // If alpha is not (approximately) zero, then j(xi) is quadratic
  // and we need to solve for the roots.
  Real sqrt_term = beta*beta - 4*alpha*gamma;

  // Debugging:
  // libMesh::out << "sqrt_term = " << sqrt_term << std::endl;

  // If the term under the square root is negative, then the roots are
  // imaginary and the element is invertible.
  if (sqrt_term < 0)
    return true;

  sqrt_term = std::sqrt(sqrt_term);
  Real
    xi_1 = 0.5 * (-beta + sqrt_term) / alpha,
    xi_2 = 0.5 * (-beta - sqrt_term) / alpha;

  // libMesh::out << "xi_1 = " << xi_1 << std::endl;
  // libMesh::out << "xi_2 = " << xi_2 << std::endl;

  // If a root is outside the reference element, it is "OK".
  bool
    xi_1_ok = ((xi_1 < -1.) || (xi_1 > 1.)),
    xi_2_ok = ((xi_2 < -1.) || (xi_2 > 1.));

  // If both roots are OK, the element is invertible.
  return xi_1_ok && xi_2_ok;
}



Order Edge4::default_order() const
{
  return THIRD;
}



void Edge4::connectivity(const unsigned int sc,
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less_equal (sc, 2);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage
  conn.resize(2);

  switch(iop)
    {
    case TECPLOT:
      {
        switch (sc)
          {
          case 0:
            conn[0] = this->node_id(0)+1;
            conn[1] = this->node_id(2)+1;
            return;

          case 1:
            conn[0] = this->node_id(2)+1;
            conn[1] = this->node_id(3)+1;
            return;

          case 2:
            conn[0] = this->node_id(3)+1;
            conn[1] = this->node_id(1)+1;
            return;

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }

      }

    case VTK:
      {

        switch (sc)
          {
          case 0:
            conn[0] = this->node_id(0);
            conn[1] = this->node_id(2);
            return;

          case 1:
            conn[0] = this->node_id(2);
            conn[1] = this->node_id(3);
            return;

          case 2:
            conn[0] = this->node_id(3);
            conn[1] = this->node_id(1);
            return;

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



BoundingBox Edge4::loose_bounding_box () const
{
  // This might be a curved line through 2-space or 3-space, in which
  // case the full bounding box can be larger than the bounding box of
  // just the nodes.
  //
  // FIXME - I haven't yet proven the formula below to be correct for
  // cubics - RHS
  Point pmin, pmax;

  for (unsigned d=0; d<LIBMESH_DIM; ++d)
    {
      Real center = (this->point(2)(d) + this->point(3)(d))/2;
      Real hd = std::max(std::abs(center - this->point(0)(d)),
                         std::abs(center - this->point(1)(d)));

      pmin(d) = center - hd;
      pmax(d) = center + hd;
    }

  return BoundingBox(pmin, pmax);
}



dof_id_type Edge4::key () const
{
  return this->compute_key(this->node_id(0),
                           this->node_id(1),
                           this->node_id(2),
                           this->node_id(3));
}



Real Edge4::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = point(0),
    x1 = point(1),
    x2 = point(2),
    x3 = point(3);

  // Construct constant data vectors.
  // \vec{x}_{\xi}  = \vec{a1}*xi**2 + \vec{b1}*xi + \vec{c1}
  // This is copy-pasted directly from the output of a Python script.
  Point
    a1 = -27*x0/16 + 27*x1/16 + 81*x2/16 - 81*x3/16,
    b1 = 9*x0/8 + 9*x1/8 - 9*x2/8 - 9*x3/8,
    c1 = x0/16 - x1/16 - 27*x2/16 + 27*x3/16;

  // 4 point quadrature, 7th-order accurate
  const unsigned int N = 4;
  const Real q[N] = {-std::sqrt(525 + 70*std::sqrt(30.)) / 35,
                     -std::sqrt(525 - 70*std::sqrt(30.)) / 35,
                     std::sqrt(525 - 70*std::sqrt(30.)) / 35,
                     std::sqrt(525 + 70*std::sqrt(30.)) / 35};
  const Real w[N] = {(18 - std::sqrt(30.)) / 36,
                     (18 + std::sqrt(30.)) / 36,
                     (18 + std::sqrt(30.)) / 36,
                     (18 - std::sqrt(30.)) / 36};

  Real vol=0.;
  for (unsigned int i=0; i<N; ++i)
    vol += w[i] * (q[i]*q[i]*a1 + q[i]*b1 + c1).norm();

  return vol;
}

} // namespace libMesh
