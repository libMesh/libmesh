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
#include "libmesh/edge_edge3.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

// Edge3 class static member initializations
const int Edge3::num_nodes;
const int Edge3::num_sides;
const int Edge3::num_edges;
const int Edge3::num_children;
const int Edge3::nodes_per_side;
const int Edge3::nodes_per_edge;

#ifdef LIBMESH_ENABLE_AMR

const float Edge3::_embedding_matrix[2][3][3] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2
      {1.0, 0.0, 0.0}, // left
      {0.0, 0.0, 1.0}, // right
      {0.375,-0.125,0.75} // middle
    },

    // embedding matrix for child 1
    {
      // 0    1    2
      {0.0, 0.0, 1.0}, // left
      {0.0, 1.0, 0.0},  // right
      {-0.125,0.375,0.75} // middle
    }
  };

#endif

bool Edge3::is_vertex(const unsigned int i) const
{
  if (i < 2)
    return true;
  return false;
}

bool Edge3::is_edge(const unsigned int i) const
{
  if (i < 2)
    return false;
  return true;
}

bool Edge3::is_face(const unsigned int ) const
{
  return false;
}

bool Edge3::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, 2);
  libmesh_assert_less (n, Edge3::num_nodes);
  return (s == n);
}

bool Edge3::is_node_on_edge(const unsigned int,
                            const unsigned int libmesh_dbg_var(e)) const
{
  libmesh_assert_equal_to (e, 0);
  return true;
}



bool Edge3::has_affine_map() const
{
  Point v = this->point(1) - this->point(0);
  return (v.relative_fuzzy_equals
          ((this->point(2) - this->point(0))*2));
}



bool Edge3::has_invertible_map(Real tol) const
{
  // At the moment this only makes sense for Lagrange elements
  libmesh_assert_equal_to(this->mapping_type(), LAGRANGE_MAP);

  // The "Jacobian vector" (dx/dxi, dy/dxi, dz/dxi) is:
  // j(xi) := a*xi + b, where
  Point a = this->point(0) + this->point(1) - 2 * this->point(2);
  Point b = .5 * (this->point(1) - this->point(0));

  // Now we solve for the point xi_m where j(xi_m) \cdot j(0) = 0.
  // If this occurs somewhere on the reference element, then the
  // element is not invertible.
  // j(xi_m) . j(0) = 0
  // <=>
  // (a*xi_m + b) . b = 0
  // <=>
  // (a.b)*xi_m + b.b = 0
  // <=>
  // xi_m = -(b.b) / (a.b)

  // 1.) If b.b==0, then the endpoints of the Edge3 are at the same
  //     location, and the map is therefore not invertible.
  Real b_norm2 = b.norm_sq();
  if (b_norm2 <= tol*tol)
    return false;

  // 2.) If a.b==0, but b != 0 (see above), then the element is
  //     invertible but we don't want to divide by zero in the
  //     formula, so simply return true.
  Real ab = a * b;
  if (std::abs(ab) <= tol*tol)
    return true;

  Real xi_m = -b_norm2 / ab;
  return (xi_m < -1.) || (xi_m > 1.);
}



Order Edge3::default_order() const
{
  return SECOND;
}



void Edge3::connectivity(const unsigned int sc,
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less_equal (sc, 1);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage
  conn.resize(2);

  switch (iop)
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
            conn[1] = this->node_id(1)+1;
            return;

          default:
            libmesh_error_msg("Invalid sc = " << sc);
          }
      }


    case VTK:
      {
        conn.resize(3);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        return;

        /*
          switch (sc)
          {
          case 0:
          conn[0] = this->node_id(0);
          conn[1] = this->node_id(2);

          return;

          case 1:
          conn[0] = this->node_id(2);
          conn[1] = this->node_id(1);

          return;

          default:
          libmesh_error_msg("Invalid sc = " << sc);
          }
        */
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



std::pair<unsigned short int, unsigned short int>
Edge3::second_order_child_vertex (const unsigned int) const
{
  return std::pair<unsigned short int, unsigned short int>(0,0);
}



Real Edge3::volume () const
{
  // Finding the (exact) length of a general quadratic element
  // is a surprisingly complicated formula.
  Point A = this->point(0) + this->point(1) - 2*this->point(2);
  Point B = (this->point(1) - this->point(0))/2;

  const Real a = A.norm_sq();
  const Real c = B.norm_sq();

  // Degenerate straight line case
  if (a == 0.)
    return 2. * std::sqrt(c);

  const Real b = -2.*std::abs(A*B);

  const Real sqrt_term1 = a - b + c;
  const Real sqrt_term2 = a + b + c;

  // Fall back on straight line case instead of computing nan
  // Note: b can be positive or negative so we have to check both cases.
  if (sqrt_term1 < 0. || sqrt_term2 < 0.)
    return 2. * std::sqrt(c);

  const Real r1 = std::sqrt(sqrt_term1);
  const Real r2 = std::sqrt(sqrt_term2);
  const Real rsum = r1 + r2;
  const Real root_a = std::sqrt(a);
  const Real b_over_root_a = b / root_a;

  // Pre-compute the denominator of the log term. If it's zero, fall
  // back on the straight line case.
  const Real log_term_denom = -root_a - 0.5*b_over_root_a + r2;
  if (log_term_denom == 0. || rsum == 0.)
    return 2. * std::sqrt(c);

  Real log1p_arg = 2*(root_a - b/rsum) / log_term_denom;

  return
    0.5*(rsum + b_over_root_a*b_over_root_a/rsum +
        (c-0.25*b_over_root_a*b_over_root_a)*std::log1p(log1p_arg)/root_a);
}



BoundingBox Edge3::loose_bounding_box () const
{
  // This might be a curved line through 2-space or 3-space, in which
  // case the full bounding box can be larger than the bounding box of
  // just the nodes.
  Point pmin, pmax;

  for (unsigned d=0; d<LIBMESH_DIM; ++d)
    {
      Real center = this->point(2)(d);
      Real hd = std::max(std::abs(center - this->point(0)(d)),
                         std::abs(center - this->point(1)(d)));

      pmin(d) = center - hd;
      pmax(d) = center + hd;
    }

  return BoundingBox(pmin, pmax);
}



dof_id_type Edge3::key () const
{
  return this->compute_key(this->node_id(2));
}


} // namespace libMesh
