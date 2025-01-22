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
#include "libmesh/face_tri.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/enum_elem_quality.h"

// C++ includes
#include <array>

namespace libMesh
{


// ------------------------------------------------------------
// Tri class static member initializations
const int Tri::num_sides;
const int Tri::num_children;

// Note: we can omit initialization of the third entry of each row because
// static variables are automatically zero-initialized.
const Real Tri::_master_points[6][3] =
  {
    {0, 0},
    {1, 0},
    {0, 1},
    {0.5, 0},
    {0.5, 0.5},
    {0, 0.5}
  };

const unsigned int Tri::adjacent_sides_map[/*num_vertices*/3][/*n_adjacent_sides*/2] =
  {
    {0, 2},  // Sides adjacent to node 0
    {0, 1},  // Sides adjacent to node 1
    {1, 2}   // Sides adjacent to node 2
  };



// ------------------------------------------------------------
// Tri class member functions
dof_id_type Tri::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  return this->compute_key(this->node_id(Tri3::side_nodes_map[s][0]),
                           this->node_id(Tri3::side_nodes_map[s][1]));
}



dof_id_type Tri::low_order_key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  return this->compute_key(this->node_id(Tri3::side_nodes_map[s][0]),
                           this->node_id(Tri3::side_nodes_map[s][1]));
}



unsigned int Tri::local_side_node(unsigned int side,
                                  unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, Tri3::nodes_per_side);

  return Tri3::side_nodes_map[side][side_node];
}



unsigned int Tri::local_edge_node(unsigned int edge,
                                  unsigned int edge_node) const
{
  return local_side_node(edge, edge_node);
}



dof_id_type Tri::key () const
{
  return this->compute_key(this->node_id(0),
                           this->node_id(1),
                           this->node_id(2));
}



std::unique_ptr<Elem> Tri::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> edge = std::make_unique<Edge2>();

  for (auto n : edge->node_index_range())
    edge->set_node(n) = this->node_ptr(Tri3::side_nodes_map[i][n]);

  return edge;
}



void Tri::side_ptr (std::unique_ptr<Elem> & side,
                    const unsigned int i)
{
  this->simple_side_ptr<Tri,Tri3>(side, i, EDGE2);
}



bool Tri::is_child_on_side(const unsigned int c,
                           const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  return (c == s || c == (s+1)%3);
}


bool Tri::is_flipped() const
{
  return (
#if LIBMESH_DIM > 2
          // Don't bother outside the XY plane
          !this->point(0)(2) && !this->point(1)(2) &&
          !this->point(2)(2) &&
#endif
          ((this->point(1)(0)-this->point(0)(0))*
           (this->point(2)(1)-this->point(0)(1)) <
           (this->point(2)(0)-this->point(0)(0))*
           (this->point(1)(1)-this->point(0)(1))));
}


std::vector<unsigned int>
Tri::edges_adjacent_to_node(const unsigned int n) const
{
  libmesh_assert_less(n, this->n_nodes());

  // For vertices, we use the Tri::adjacent_sides_map, otherwise each
  // of the mid-edge nodes is adjacent only to the edge it is on, and the
  // center node is not adjacent to any edge.
  if (this->is_vertex(n))
    return {std::begin(adjacent_sides_map[n]), std::end(adjacent_sides_map[n])};
  else if (this->is_edge(n))
    return {n - this->n_vertices()};

  libmesh_assert(this->is_face(n));
  return {};
}


Real Tri::quality (const ElemQuality q) const
{
  switch (q)
    {
    case ASPECT_RATIO:
      {
        // Aspect Ratio definition from Ansys Theory Manual.
        // Reference: Ansys, Inc. Theory Reference, Ansys Release 9.0, 2004 (Chapter: 13.7.3)

        // Compute midpoint positions along each edge
        Point m[3] = {
          0.5 * (this->point(0) + this->point(1)),  // side opposite vertex 2
          0.5 * (this->point(1) + this->point(2)),  // side opposite vertex 0
          0.5 * (this->point(2) + this->point(0))}; // side opposite vertex 1

        // opposite[i] is the side index which is "opposite" vertex i
        static const unsigned int opposite[3] = {1, 2, 0};

        // other[i] is the side index which is _not_ i and _not_ opposite[i]
        static const unsigned int other[3] = {2, 0, 1};

        // Input is vertex index, i = 0, 1, 2
        auto vertex_aspect_ratio = [&](unsigned int i) -> Real
        {
          // Compute vectors:
          // v0: (vertex v, opposite midpoint)
          // v1: (midpoint[i], last midpoint)
          Point v0 = m[opposite[i]] - this->point(i);
          Point v1 = m[other[i]] - m[i];

          // Compute the length of the midlines
          Real v0_norm = v0.norm();
          Real v1_norm = v1.norm();

          // Instead of dividing by zero in the next step, just return
          // 0.  The optimal aspect ratio is 1.0, and "high" aspect
          // ratios are bad, but an aspect ratio of 0 should also be
          // considered bad.
          if (v0_norm == 0. || v1_norm == 0.)
            return 0.;

          // Compute sine of the angle between v0, v1.
          Real sin_theta = cross_norm(v0, v1) / v0_norm / v1_norm;
          Real v0s = v0_norm*sin_theta;
          Real v1s = v1_norm*sin_theta;

          // Determine the min, max of each midline length and its
          // projection.
          auto [min0, max0] = std::minmax(v0_norm, v1s);
          auto [min1, max1] = std::minmax(v0s, v1_norm);

          // Return the max of the two quotients
          return std::max(max0/min0, max1/min1);
        };

        return std::max(std::max(vertex_aspect_ratio(0), vertex_aspect_ratio(1)), vertex_aspect_ratio(2)) / std::sqrt(3);
      }

      /**
       * Source: Netgen, meshtool.cpp, TriangleQualityInst
       */
    case DISTORTION:
    case STRETCH:
      {
        const Point & p1 = this->point(0);
        const Point & p2 = this->point(1);
        const Point & p3 = this->point(2);

        Point v1 = p2 - p1;
        Point v2 = p3 - p1;
        Point v3 = p3 - p2;
        const Real l1 = v1.norm();
        const Real l2 = v2.norm();
        const Real l3 = v3.norm();

        // if one length is 0, quality is quite bad!
        if ((l1 <=0.) || (l2 <= 0.) || (l3 <= 0.))
          return 0.;

        const Real s1 = std::sin(std::acos(v1*v2/l1/l2)/2.);
        v1 *= -1;
        const Real s2 = std::sin(std::acos(v1*v3/l1/l3)/2.);
        const Real s3 = std::sin(std::acos(v2*v3/l2/l3)/2.);

        return 8. * s1 * s2 * s3;

      }

      // From: P. Knupp, "Algebraic mesh quality metrics for
      // unstructured initial meshes," Finite Elements in Analysis
      // and Design 39, 2003, p. 217-241, Section 3.2.
    case SHAPE:
      {
        // Unlike Quads, the Tri SHAPE metric is independent of the
        // node at which it is computed, we choose to compute it for
        // node 0.

        // The nodal Jacobian matrix A is a 3x2 matrix, hence we
        // represent it by a std:array with 6 entries.
        Point
          d01 = point(1) - point(0),
          d02 = point(2) - point(0);

        std::array<Real, 6> A =
          {{d01(0), d02(0),
            d01(1), d02(1),
            d01(2), d02(2)}};

        // Compute metric tensor entries, T = A^T * A.
        // This is a symmetric 2x2 matrix so we only
        // compute one of the off-diagonal entries.
        // As in the paper, we define lambda_ij := T_ij.
        Real
          lambda11 = A[0]*A[0] + A[2]*A[2] + A[4]*A[4],
          lambda12 = A[0]*A[1] + A[2]*A[3] + A[4]*A[5],
          lambda22 = A[1]*A[1] + A[3]*A[3] + A[5]*A[5];

        // Compute the denominator of the metric. If it is exactly
        // zero then return 0 (lowest quality) for this metric.
        Real den = lambda11 + lambda22 - lambda12;
        if (den == 0.0)
          return 0.;

        // Compute the nodal area
        Real alpha = std::sqrt(lambda11 * lambda22 - lambda12 * lambda12);

        // Finally, compute and return the metric.
        return std::sqrt(3) * alpha / den;
      }

    default:
      return Elem::quality(q);
    }

  // We won't get here.
  return Elem::quality(q);
}






std::pair<Real, Real> Tri::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {
      // A recent copy of the cubit manual [0] does not list bounds
      // for EDGE_LENGTH_RATIO or ASPECT_RATIO quality metrics, so we
      // have arbitrarily adopted the same values used for Quads here.
      // I'm open to suggestions of other appropriate values.
      //
      // [0]: https://cubit.sandia.gov/files/cubit/16.08/help_manual/WebHelp/mesh_generation/mesh_quality_assessment/triangular_metrics.htm
    case EDGE_LENGTH_RATIO:
    case ASPECT_RATIO:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case MAX_ANGLE:
      bounds.first  = 60.;
      bounds.second = 90.;
      break;

    case MIN_ANGLE:
      bounds.first  = 30.;
      bounds.second = 60.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 1.3;
      break;

    case JACOBIAN:
    case SCALED_JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.155;
      break;

    case SIZE:
    case SHAPE:
      bounds.first  = 0.25;
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

} // namespace libMesh
