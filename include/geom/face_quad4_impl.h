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

// Local includes
#include "libmesh/side.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

// ------------------------------------------------------------
// Quad4 class member functions

template <typename RealType>
bool Quad4Templ<RealType>::is_vertex(const unsigned int) const
{
  return true;
}

template <typename RealType>
bool Quad4Templ<RealType>::is_edge(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Quad4Templ<RealType>::is_face(const unsigned int) const
{
  return false;
}

template <typename RealType>
bool Quad4Templ<RealType>::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

template <typename RealType>
std::vector<unsigned>
Quad4Templ<RealType>::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, this->n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

template <typename RealType>
bool Quad4Templ<RealType>::has_affine_map() const
{
  Point v = this->point(3) - this->point(0);
  return (v.relative_fuzzy_equals(this->point(2) - this->point(1)));
}



template <typename RealType>
Order Quad4Templ<RealType>::default_order() const
{
  return FIRST;
}



template <typename RealType>
std::unique_ptr<ElemTempl<RealType>> Quad4Templ<RealType>::build_side_ptr (const unsigned int i,
                                             bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return libmesh_make_unique<Side<Edge2,Quad4>>(this,i);

  else
    {
      std::unique_ptr<Elem> edge = libmesh_make_unique<Edge2>();
      edge->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (auto n : edge->node_index_range())
        edge->set_node(n) = this->node_ptr(Quad4::side_nodes_map[i][n]);

      return edge;
    }
}



template <typename RealType>
void Quad4Templ<RealType>::build_side_ptr (std::unique_ptr<Elem> & side,
                            const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  if (!side.get() || side->type() != EDGE2)
    side = this->build_side_ptr(i, false);
  else
    {
      side->subdomain_id() = this->subdomain_id();

      for (auto n : side->node_index_range())
        side->set_node(n) = this->node_ptr(Quad4::side_nodes_map[i][n]);
    }
}



template <typename RealType>
void Quad4Templ<RealType>::connectivity(const unsigned int libmesh_dbg_var(sf),
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



template <typename RealType>
RealType Quad4Templ<RealType>::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = this->point(0), x1 = this->point(1),
    x2 = this->point(2), x3 = this->point(3);

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

  RealType vol=0.;
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<2; ++j)
      vol += cross_norm(q[j]*a1 + b1,
                        q[i]*a2 + b2);

  return vol;
}

template <typename RealType>
BoundingBoxTempl<RealType>
Quad4Templ<RealType>::loose_bounding_box () const
{
  return Elem::loose_bounding_box();
}


} // namespace libMesh
