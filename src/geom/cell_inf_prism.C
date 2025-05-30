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

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/cell_inf_prism.h"
#include "libmesh/cell_inf_prism6.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_inf_quad4.h"
#include "libmesh/fe_type.h"
#include "libmesh/fe_interface.h"
#include "libmesh/inf_fe_map.h"
#include "libmesh/enum_order.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfPrism class static member initializations
const int InfPrism::num_sides;
const int InfPrism::num_edges;
const int InfPrism::num_children;

const Real InfPrism::_master_points[12][3] =
  {
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1},
    {1, 0, 1},
    {0, 1, 1},
    {.5, 0, 0},
    {.5, .5, 0},
    {0, .5, 0},
    {.5, 0, 1},
    {.5, .5, 1},
    {0, .5, 1}
  };

const unsigned int InfPrism::edge_sides_map[6][2] =
  {
    {0, 1}, // Edge 0
    {0, 2}, // Edge 1
    {0, 3}, // Edge 2
    {1, 3}, // Edge 3
    {1, 2}, // Edge 4
    {2, 3}  // Edge 5
  };

const unsigned int InfPrism::adjacent_edges_map[/*num_vertices*/6][/*max_adjacent_edges*/3] =
{
  {0,  2,  3}, // Edges adjacent to node 0
  {0,  1,  4}, // Edges adjacent to node 1
  {1,  2,  5}, // Edges adjacent to node 2
  {3, 99, 99}, // Edges adjacent to node 3
  {4, 99, 99}, // Edges adjacent to node 4
  {5, 99, 99}, // Edges adjacent to node 5
};

// ------------------------------------------------------------
// InfPrism class member functions
dof_id_type InfPrism::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0: // the triangular face at z=-1, base face
      return this->compute_key (this->node_id(InfPrism6::side_nodes_map[s][0]),
                                this->node_id(InfPrism6::side_nodes_map[s][1]),
                                this->node_id(InfPrism6::side_nodes_map[s][2]));

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      return this->compute_key (this->node_id(InfPrism6::side_nodes_map[s][0]),
                                this->node_id(InfPrism6::side_nodes_map[s][1]),
                                this->node_id(InfPrism6::side_nodes_map[s][2]),
                                this->node_id(InfPrism6::side_nodes_map[s][3]));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }
}



dof_id_type InfPrism::low_order_key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0: // the triangular face at z=-1, base face
      return this->compute_key (this->node_id(InfPrism6::side_nodes_map[s][0]),
                                this->node_id(InfPrism6::side_nodes_map[s][1]),
                                this->node_id(InfPrism6::side_nodes_map[s][2]));

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      return this->compute_key (this->node_id(InfPrism6::side_nodes_map[s][0]),
                                this->node_id(InfPrism6::side_nodes_map[s][1]),
                                this->node_id(InfPrism6::side_nodes_map[s][2]),
                                this->node_id(InfPrism6::side_nodes_map[s][3]));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }
}



unsigned int InfPrism::local_side_node(unsigned int side,
                                       unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  // Never more than 4 nodes per side.
  libmesh_assert_less(side_node, InfPrism6::nodes_per_side);

  // Some sides have 3 nodes.
  libmesh_assert(side != 0 || side_node < 3);

  return InfPrism6::side_nodes_map[side][side_node];
}



unsigned int InfPrism::local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const
{
  libmesh_assert_less(edge, this->n_edges());
  libmesh_assert_less(edge_node, InfPrism6::nodes_per_edge);

  return InfPrism6::edge_nodes_map[edge][edge_node];
}



std::unique_ptr<Elem> InfPrism::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;

  switch (i)
    {
    case 0: // the triangular face at z=-1, base face
      {
        face = std::make_unique<Tri3>();
        break;
      }

    case 1: // the quad face at y=0
    case 2: // the other quad face
    case 3: // the quad face at x=0
      {
        face = std::make_unique<InfQuad4>();
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(InfPrism6::side_nodes_map[i][n]));

  return face;
}



void InfPrism::side_ptr (std::unique_ptr<Elem> & side,
                         const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  switch (i)
    {
      // the base face
    case 0:
      {
        if (!side.get() || side->type() != TRI3)
          {
            side = this->side_ptr(i);
            return;
          }
        break;
      }

      // connecting to another infinite element
    case 1:
    case 2:
    case 3:
    case 4:
      {
        if (!side.get() || side->type() != INFQUAD4)
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
    side->set_node(n, this->node_ptr(InfPrism6::side_nodes_map[i][n]));
}



bool InfPrism::is_child_on_side(const unsigned int c,
                                const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  return (s == 0 || c+1 == s || c == s%3);
}



bool InfPrism::is_edge_on_side (const unsigned int e,
                                const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (edge_sides_map[e][0] == s || edge_sides_map[e][1] == s);
}

bool InfPrism::contains_point (const Point & p, Real tol) const
{
  // For infinite elements with linear base interpolation:
  // make use of the fact that infinite elements do not
  // live inside the envelope.  Use a fast scheme to
  // check whether point \p p is inside or outside
  // our relevant part of the envelope.  Note that
  // this is not exclusive: only when the distance is less,
  // we are safe.  Otherwise, we cannot say anything. The
  // envelope may be non-spherical, the physical point may lie
  // inside the envelope, outside the envelope, or even inside
  // this infinite element.  Therefore if this fails,
  // fall back to the inverse_map()
  const Point my_origin (this->origin());

  // determine the minimal distance of the base from the origin
  // use norm_sq(), it is faster than norm() and produces
  // the same behavior. Shift my_origin to center first:
  Point pt0_o(this->point(0) - my_origin);
  Point pt1_o(this->point(1) - my_origin);
  Point pt2_o(this->point(2) - my_origin);
  const Real tmp_min_distance_sq = std::min(pt0_o.norm_sq(),
                                            std::min(pt1_o.norm_sq(),
                                                     pt2_o.norm_sq()));

  // For a coarse grid, it is important to account for the fact
  // that the sides are not spherical, thus the closest point
  // can be closer than all edges.
  // This is an estimator using Pythagoras:
  const Real min_distance_sq = tmp_min_distance_sq
                              - .5*std::max((point(0)-point(1)).norm_sq(),
                                             std::max((point(0)-point(2)).norm_sq(),
                                                      (point(1)-point(2)).norm_sq()));

  // work with 1% allowable deviation.  We can still fall
  // back to the InfFE::inverse_map()
  const Real conservative_p_dist_sq = 1.01 * (Point(p - my_origin).norm_sq());

  if (conservative_p_dist_sq < min_distance_sq)
    {
      // the physical point is definitely not contained in the element:
      return false;
    }

  // this captures the case that the point is not (almost) in the direction of the element.:
  // first, project the problem onto the unit sphere:
  Point p_o(p - my_origin);
  pt0_o /= pt0_o.norm();
  pt1_o /= pt1_o.norm();
  pt2_o /= pt2_o.norm();
  p_o /= p_o.norm();

  // now, check if it is in the projected face; by comparing the distance of
  // any point in the element to \p p with the largest distance between this point
  // to any other point in the element.
  if ((p_o - pt0_o).norm_sq() > std::max((pt0_o - pt1_o).norm_sq(), (pt0_o - pt2_o).norm_sq()) ||
      (p_o - pt1_o).norm_sq() > std::max((pt1_o - pt2_o).norm_sq(), (pt1_o - pt0_o).norm_sq()) ||
      (p_o - pt2_o).norm_sq() > std::max((pt2_o - pt0_o).norm_sq(), (pt2_o - pt1_o).norm_sq()) )
    {
      // the physical point is definitely not contained in the element
      return false;
    }


  // It seems that \p p is at least close to this element.
  //
  // Declare a basic FEType.  Will use default in the base,
  // and something else (not important) in radial direction.
  FEType fe_type(default_order());

  const Point mapped_point = InfFEMap::inverse_map(dim(), this, p,
                                                   tol, false);

  return this->on_reference_element(mapped_point, tol);
}

std::vector<unsigned>
InfPrism::sides_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_sides_map[e]), std::end(edge_sides_map[e])};
}


bool
InfPrism::is_flipped() const
{
  return (triple_product(this->point(1)-this->point(0),
                         this->point(2)-this->point(0),
                         this->point(3)-this->point(0)) < 0);
}

std::vector<unsigned int>
InfPrism::edges_adjacent_to_node(const unsigned int n) const
{
  libmesh_assert_less(n, this->n_nodes());

  // For vertices, we use the InfPrism::adjacent_edges_map with
  // appropriate "trimming" based on whether the vertices are "at
  // infinity" or not.  Otherwise each of the mid-edge nodes is
  // adjacent only to the edge it is on, and the remaining nodes are
  // not adjacent to any edge.
  if (this->is_vertex(n))
    {
      auto trim = (n < 3) ? 0 : 2;
      return {std::begin(adjacent_edges_map[n]), std::end(adjacent_edges_map[n]) - trim};
    }
  else if (this->is_edge(n))
    return {n - this->n_vertices()};

  libmesh_assert(this->is_face(n) || this->is_internal(n));
  return {};
}



bool InfPrism::on_reference_element(const Point & p,
                                    const Real eps) const
{
  const Real & xi = p(0);
  const Real & eta = p(1);
  const Real & zeta = p(2);

  // inside the reference triangle with zeta in [-1,1]
  return ((xi   >=  0.-eps) &&
          (eta  >=  0.-eps) &&
          (zeta >= -1.-eps) &&
          (zeta <=  1.+eps) &&
          ((xi + eta) <= 1.+eps));
}




} // namespace libMesh



#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
