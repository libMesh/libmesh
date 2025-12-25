// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes cont'd
#include "libmesh/face_inf_quad.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/face_inf_quad4.h"
#include "libmesh/enum_elem_quality.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfQuad class static member initializations
const int InfQuad::num_sides;
const int InfQuad::num_children;

// Note: we can omit initialization of the third entry of each row because
// static variables are automatically zero-initialized.
const Real InfQuad::_master_points[6][3] =
  {
    {-1, 0},
    {1, 0},
    {-1, 1},
    {1, 1},
    {0, 0},
    {0, 1}
  };

const unsigned int InfQuad::adjacent_sides_map[/*num_vertices*/4][/*max_adjacent_sides*/2] =
  {
    {0,  2}, // Sides adjacent to node 0
    {0,  1}, // Sides adjacent to node 1
    {2, 99}, // Sides adjacent to node 2
    {1, 99}  // Sides adjacent to node 3
  };

// ------------------------------------------------------------
// InfQuad class member functions
dof_id_type InfQuad::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  // The order of the node ids does not matter, they are sorted by the
  // compute_key() function.
  return this->compute_key(this->node_id(InfQuad4::side_nodes_map[s][0]),
                           this->node_id(InfQuad4::side_nodes_map[s][1]));
}



dof_id_type InfQuad::low_order_key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  // The order of the node ids does not matter, they are sorted by the
  // compute_key() function.
  return this->compute_key(this->node_id(InfQuad4::side_nodes_map[s][0]),
                           this->node_id(InfQuad4::side_nodes_map[s][1]));
}



unsigned int InfQuad::local_side_node(unsigned int side,
                                      unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, 2);

  return InfQuad4::side_nodes_map[side][side_node];
}



unsigned int InfQuad::local_edge_node(unsigned int edge,
                                      unsigned int edge_node) const
{
  return local_side_node(edge, edge_node);
}



std::unique_ptr<Elem> InfQuad::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  // Return value
  std::unique_ptr<Elem> edge;

  switch (i)
    {
    case 0: // base face
      {
        edge = std::make_unique<Edge2>();
        break;
      }

    case 1: // adjacent to another infinite element
    case 2: // adjacent to another infinite element
      {
        edge = std::make_unique<InfEdge2>();
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : edge->node_index_range())
    edge->set_node(n, this->node_ptr(InfQuad4::side_nodes_map[i][n]));

  return edge;
}



void InfQuad::side_ptr (std::unique_ptr<Elem> & side,
                        const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
      // the base face
    case 0:
      {
        if (!side.get() || side->type() != EDGE2)
          {
            side = this->side_ptr(i);
            return;
          }
        break;
      }

      // connecting to another infinite element
    case 1:
    case 2:
      {
        if (!side.get() || side->type() != INFEDGE2)
          {
            side = this->side_ptr(i);
            return;
          }
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  side->subdomain_id() = this->subdomain_id();

  // Set the nodes
  for (auto n : side->node_index_range())
    side->set_node(n, this->node_ptr(InfQuad4::side_nodes_map[i][n]));
}


bool InfQuad::is_child_on_side(const unsigned int c,
                               const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  return (s == 0 ||
          (c == 0 && s == 2) ||
          (c == 1 && s == 1));
}


bool
InfQuad::is_flipped() const
{
  return (
#if LIBMESH_DIM > 2
          // Don't bother outside the XY plane
          this->point(0)(2) && this->point(1)(2) &&
          this->point(2)(2) && this->point(3)(2) &&
#endif
          ((this->point(1)(0)-this->point(0)(0))*
           (this->point(2)(1)-this->point(0)(1)) >
           (this->point(2)(0)-this->point(0)(0))*
           (this->point(1)(1)-this->point(0)(1))));
}


std::vector<unsigned int>
InfQuad::edges_adjacent_to_node(const unsigned int n) const
{
  libmesh_assert_less(n, this->n_nodes());

  // For vertices, we use the adjacent_sides_map, otherwise node
  // 4 is on side 0 and node 5 is not any any side.
  if (this->is_vertex(n))
    {
      auto trim = (n < 2) ? 0 : 1;
      return {std::begin(adjacent_sides_map[n]), std::end(adjacent_sides_map[n]) - trim};
    }
  else if (n == 4)
    return {0};
  else
    return {};
}


Real InfQuad::quality (const ElemQuality q) const
{
  return Elem::quality(q); // Not implemented
}




std::pair<Real, Real> InfQuad::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {

    case ASPECT_RATIO:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case SKEW:
      bounds.first  = 0.;
      bounds.second = 0.5;
      break;

    case TAPER:
      bounds.first  = 0.;
      bounds.second = 0.7;
      break;

    case WARP:
      bounds.first  = 0.9;
      bounds.second = 1.;
      break;

    case STRETCH:
      bounds.first  = 0.25;
      bounds.second = 1.;
      break;

    case MIN_ANGLE:
      bounds.first  = 45.;
      bounds.second = 90.;
      break;

    case MAX_ANGLE:
      bounds.first  = 90.;
      bounds.second = 135.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case JACOBIAN:
    case SCALED_JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.;
      break;

    case SHEAR:
    case SHAPE:
    case SIZE:
      bounds.first  = 0.3;
      bounds.second = 1.;
      break;

    case DISTORTION:
      bounds.first  = 0.6;
      bounds.second = 1.;
      break;

    default:
      libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}



bool InfQuad::on_reference_element(const Point & p,
                                   const Real eps) const
{
  const Real & xi = p(0);
  const Real & eta = p(1);

  // The reference infquad is [-1,1]^2.
  return ((xi  >= -1.-eps) &&
          (xi  <=  1.+eps) &&
          (eta >= -1.-eps) &&
          (eta <=  1.+eps));
}



} // namespace libMesh



#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
