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
#include "libmesh/side.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{




// ------------------------------------------------------------
// Quad4 class static member initialization
const int Quad4::num_nodes;
const int Quad4::num_sides;
const int Quad4::num_children;
const int Quad4::nodes_per_side;

const unsigned int Quad4::side_nodes_map[Quad4::num_sides][Quad4::nodes_per_side] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 3}, // Side 2
    {3, 0}  // Side 3
  };


#ifdef LIBMESH_ENABLE_AMR

const Real Quad4::_embedding_matrix[Quad4::num_children][Quad4::num_nodes][Quad4::num_nodes] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2    3
      {1.0, 0.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0, 0.0}, // 1
      {.25, .25, .25, .25}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    },

    // embedding matrix for child 1
    {
      // 0    1    2    3
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 1.0, 0.0, 0.0}, // 1
      {0.0, 0.5, 0.5, 0.0}, // 2
      {.25, .25, .25, .25}  // 3
    },

    // embedding matrix for child 2
    {
      // 0    1    2    3
      {0.5, 0.0, 0.0, 0.5}, // 0
      {.25, .25, .25, .25}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.0, 0.0, 0.0, 1.0}  // 3
    },

    // embedding matrix for child 3
    {
      // 0    1    2    3
      {.25, .25, .25, .25}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.0, 0.0, 1.0, 0.0}, // 2
      {0.0, 0.0, 0.5, 0.5}  // 3
    }
  };

#endif





// ------------------------------------------------------------
// Quad4 class member functions

bool Quad4::is_vertex(const unsigned int) const
{
  return true;
}

bool Quad4::is_edge(const unsigned int) const
{
  return false;
}

bool Quad4::is_face(const unsigned int) const
{
  return false;
}

bool Quad4::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Quad4::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
Quad4::nodes_on_edge(const unsigned int e) const
{
  return nodes_on_side(e);
}

bool Quad4::has_affine_map() const
{
  Point v = this->point(3) - this->point(0);
  return (v.relative_fuzzy_equals(this->point(2) - this->point(1), affine_tol));
}



bool Quad4::has_invertible_map(Real tol) const
{
  // At the moment this only makes sense for Lagrange elements
  libmesh_assert_equal_to(this->mapping_type(), LAGRANGE_MAP);

  // Side vectors
  Point s0 = this->point(1) - this->point(0);
  Point s1 = this->point(2) - this->point(1);
  Point s2 = this->point(3) - this->point(2);
  Point s3 = this->point(0) - this->point(3);

  // Cross products of side vectors
  Point v1 = s3.cross(s0);
  Point v2 = s2.cross(s0);
  Point v3 = s3.cross(s1);

  // A unit vector in the direction of:
  // f(xi, eta) = (v1 + xi*v2 + eta*v3)
  // at the midpoint (xi, eta) = (1/2, 1/2) of the element. (Note that
  // we are using the [0,1]^2 reference element definition for the
  // Quad4 instead of the [-1,1]^2 reference element that is typically
  // used for FEM calculations.) We use this as a "reference" vector
  // and compare the sign of dot(n,f) at each vertex.
  Point n = v1 + Real(.5) * (v2 + v3);
  Real norm_n = n.norm();

  // If the Jacobian vector at the midpoint of the element is zero,
  // then it must be either zero for the entire element, or change
  // sign at the vertices. Either way the element is non-invertible.
  if (norm_n <= tol)
    return false;

  n /= norm_n;

  // Debugging
  // std::cout << "n=" << n << std::endl;

  // Compute scalar quantity n * (v1 + xi*v2 + eta*v3) at each
  // vertex. If it is non-zero and has the same sign at each
  // vertex, the the element is invertible, otherwise it is not.
  std::array<Real, 4> vertex_vals;
  unsigned int ctr = 0;
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<2; ++j)
      vertex_vals[ctr++] = n * (v1 + Real(i)*v2 + Real(j)*v3);

  // Debugging:
  // std::cout << "Vertex values: ";
  // for (const auto & val : vertex_vals)
  //   std::cout << val << " ";
  // std::cout << std::endl;

  auto result = std::minmax_element(vertex_vals.begin(), vertex_vals.end());
  Real min_vertex = *(result.first);
  Real max_vertex = *(result.second);

  // Debugging
  // std::cout << "min_vertex=" << min_vertex << std::endl;
  // std::cout << "max_vertex=" << max_vertex << std::endl;

  // If max and min are both on the same side of 0, we are invertible, otherwise we are not.
  return ((max_vertex > 0 && min_vertex > 0) || (max_vertex < 0 && min_vertex < 0));
}



Order Quad4::default_order() const
{
  return FIRST;
}



std::unique_ptr<Elem> Quad4::build_side_ptr (const unsigned int i,
                                             bool proxy)
{
  return this->simple_build_side_ptr<Edge2, Quad4>(i, proxy);
}



void Quad4::build_side_ptr (std::unique_ptr<Elem> & side,
                            const unsigned int i)
{
  this->simple_build_side_ptr<Quad4>(side, i, EDGE2);
}



void Quad4::connectivity(const unsigned int libmesh_dbg_var(sf),
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage.
  conn.resize(4);

  switch (iop)
    {
    case TECPLOT:
      {
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(3)+1;
        return;
      }

    case VTK:
      {
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        conn[3] = this->node_id(3);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



Point Quad4::true_centroid () const
{
  // Convenient references to our points
  const Point
    &x0 = point(0), &x1 = point(1),
    &x2 = point(2), &x3 = point(3);

  // Construct "dx/d(xi)" and "dx/d(eta)" vectors which are columns of the Jacobian.
  // \vec{x}_{\xi}  = \vec{a1}*eta + \vec{b1}
  // \vec{x}_{\eta} = \vec{a2}*xi  + \vec{b2}
  // This is copy-pasted directly from the output of a Python script. Note: we are off
  // by a factor of (1/4) here, but since the final result should have the same error
  // in the numerator and denominator, it should cancel out while saving us some math
  // operations.
  Point
    a1 = x0 - x1 + x2 - x3,
    b1 = -x0 + x1 + x2 - x3,
    a2 = a1,
    b2 = -x0 - x1 + x2 + x3;

  // Use 2x2 quadrature to compute the integral of each basis function
  // (as defined on the [-1,1]^2 reference domain). We use a 4-point
  // rule, which is exact for bi-cubics. The weights for this rule are
  // all equal to 1.
  const Real q[2] = {-std::sqrt(3.)/3, std::sqrt(3.)/3.};

  // Nodal areas
  Real A0 = 0., A1 = 0., A2 = 0., A3 = 0.;

  for (const auto & xi : q)
    for (const auto & eta : q)
      {
        Real jxw = cross_norm(eta*a1 + b1, xi*a2 + b2);

        A0 += jxw * (1-xi) * (1-eta); // 4 * phi_0
        A1 += jxw * (1+xi) * (1-eta); // 4 * phi_1
        A2 += jxw * (1+xi) * (1+eta); // 4 * phi_2
        A3 += jxw * (1-xi) * (1+eta); // 4 * phi_3
      }

  // Compute centroid
  return Point(x0(0)*A0 + x1(0)*A1 + x2(0)*A2 + x3(0)*A3,
               x0(1)*A0 + x1(1)*A1 + x2(1)*A2 + x3(1)*A3,
               x0(2)*A0 + x1(2)*A1 + x2(2)*A2 + x3(2)*A3) / (A0 + A1 + A2 + A3);
}



Real Quad4::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = point(0), x1 = point(1),
    x2 = point(2), x3 = point(3);

  // Construct constant data vectors.
  // \vec{x}_{\xi}  = \vec{a1}*eta + \vec{b1}
  // \vec{x}_{\eta} = \vec{a2}*xi  + \vec{b2}
  // This is copy-pasted directly from the output of a Python script.
  Point
    a1 = x0/4 - x1/4 + x2/4 - x3/4,
    b1 = -x0/4 + x1/4 + x2/4 - x3/4,
    a2 = a1,
    b2 = -x0/4 - x1/4 + x2/4 + x3/4;

  // Check for quick return for parallelogram QUAD4.
  if (a1.relative_fuzzy_equals(Point(0,0,0)))
    return 4. * b1.cross(b2).norm();

  // Otherwise, use 2x2 quadrature to approximate the surface area.

  // 4-point rule, exact for bi-cubics.  The weights for this rule are
  // all equal to 1.
  const Real q[2] = {-std::sqrt(3.)/3, std::sqrt(3.)/3.};

  Real vol=0.;
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<2; ++j)
      vol += cross_norm(q[j]*a1 + b1,
                        q[i]*a2 + b2);

  return vol;
}

BoundingBox
Quad4::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}


void Quad4::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 4);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap4nodes(0,1,2,3);
      swap4neighbors(0,1,2,3);
    }
}


ElemType Quad4::side_type (const unsigned int libmesh_dbg_var(s)) const
{
  libmesh_assert_less (s, 4);
  return EDGE2;
}

} // namespace libMesh
