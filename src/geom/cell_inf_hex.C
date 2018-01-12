// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// C++ includes
#include <algorithm> // for std::min, std::max

// Local includes cont'd
#include "libmesh/cell_inf_hex.h"
#include "libmesh/cell_inf_hex8.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_inf_quad4.h"
#include "libmesh/fe_type.h"
#include "libmesh/fe_interface.h"

namespace libMesh
{



// ------------------------------------------------------------
// InfHex class static member initializations


// We need to require C++11...
const Real InfHex::_master_points[18][3] =
  {
    {-1, -1, 0},
    {1, -1, 0},
    {1, 1, 0},
    {-1, 1, 0},
    {-1, -1, 1},
    {1, -1, 1},
    {1, 1, 1},
    {-1, 1, 1},
    {0, -1, 0},
    {1, 0, 0},
    {0, 1, 0},
    {-1, 0, 0},
    {0, -1, 1},
    {1, 0, 1},
    {0, 1, 1},
    {-1, 0, 1},
    {0, 0, 0},
    {0, 0, 1}
  };






// ------------------------------------------------------------
// InfHex class member functions
dof_id_type InfHex::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  // The order of the node ids does not matter, they are sorted by the
  // compute_key() function.
  return this->compute_key(this->node_id(InfHex8::side_nodes_map[s][0]),
                           this->node_id(InfHex8::side_nodes_map[s][1]),
                           this->node_id(InfHex8::side_nodes_map[s][2]),
                           this->node_id(InfHex8::side_nodes_map[s][3]));
}



unsigned int InfHex::which_node_am_i(unsigned int side,
                                     unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, 4);

  return InfHex8::side_nodes_map[side][side_node];
}



std::unique_ptr<Elem> InfHex::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  // Return value
  std::unique_ptr<Elem> face;

  // Think of a unit cube: (-1,1) x (-1,1) x (-1,1),
  // with (in general) the normals pointing outwards
  switch (i)
    {
      // the face at z = -1
      // the base, where the infinite element couples to conventional
      // elements
    case 0:
      {
        // Oops, here we are, claiming the normal of the face
        // elements point outwards -- and this is the exception:
        // For the side built from the base face,
        // the normal is pointing _into_ the element!
        // Why is that? - In agreement with build_side_ptr(),
        // which in turn _has_ to build the face in this
        // way as to enable the cool way \p InfFE re-uses \p FE.
        face = libmesh_make_unique<Quad4>();
        break;
      }

      // These faces connect to other infinite elements.
    case 1: // the face at y = -1
    case 2: // the face at x = 1
    case 3: // the face at y = 1
    case 4: // the face at x = -1
      {
        face = libmesh_make_unique<InfQuad4>();
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (unsigned n=0; n<face->n_nodes(); ++n)
    face->set_node(n) = this->node_ptr(InfHex8::side_nodes_map[i][n]);

  return face;
}



bool InfHex::is_child_on_side(const unsigned int c,
                              const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  return (s == 0 || c+1 == s || c == s%4);
}



bool InfHex::is_edge_on_side (const unsigned int e,
                              const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_sides());

  return (is_node_on_side(InfHex8::edge_nodes_map[e][0],s) &&
          is_node_on_side(InfHex8::edge_nodes_map[e][1],s));
}




Real InfHex::quality (const ElemQuality q) const
{
  switch (q)
    {

      /**
       * Compute the min/max diagonal ratio.
       * Source: CUBIT User's Manual.
       *
       * For infinite elements, we just only compute
       * the diagonal in the face...
       * Don't know whether this makes sense,
       * but should be a feasible way.
       */
    case DIAGONAL:
      {
        // Diagonal between node 0 and node 2
        const Real d02 = this->length(0,2);

        // Diagonal between node 1 and node 3
        const Real d13 = this->length(1,3);

        // Find the biggest and smallest diagonals
        const Real min = std::min(d02, d13);
        const Real max = std::max(d02, d13);

        libmesh_assert_not_equal_to (max, 0.0);

        return min / max;

        break;
      }

      /**
       * Minimum ratio of lengths derived from opposite edges.
       * Source: CUBIT User's Manual.
       *
       * For IFEMs, do this only for the base face...
       * Does this make sense?
       */
    case TAPER:
      {

        /**
         * Compute the side lengths.
         */
        const Real d01 = this->length(0,1);
        const Real d12 = this->length(1,2);
        const Real d23 = this->length(2,3);
        const Real d03 = this->length(0,3);

        std::vector<Real> edge_ratios(2);

        // Bottom
        edge_ratios[8] = std::min(d01, d23) / std::max(d01, d23);
        edge_ratios[9] = std::min(d03, d12) / std::max(d03, d12);

        return *(std::min_element(edge_ratios.begin(), edge_ratios.end())) ;

        break;
      }


      /**
       * Minimum edge length divided by max diagonal length.
       * Source: CUBIT User's Manual.
       *
       * And again, we mess around a bit, for the IFEMs...
       * Do this only for the base.
       */
    case STRETCH:
      {
        /**
         * Should this be a sqrt2, when we do this for the base only?
         */
        const Real sqrt3 = 1.73205080756888;

        /**
         * Compute the maximum diagonal in the base.
         */
        const Real d02 = this->length(0,2);
        const Real d13 = this->length(1,3);
        const Real max_diag = std::max(d02, d13);

        libmesh_assert_not_equal_to ( max_diag, 0.0 );

        /**
         * Compute the minimum edge length in the base.
         */
        std::vector<Real> edges(4);
        edges[0]  = this->length(0,1);
        edges[1]  = this->length(1,2);
        edges[2]  = this->length(2,3);
        edges[3]  = this->length(0,3);

        const Real min_edge = *(std::min_element(edges.begin(), edges.end()));
        return sqrt3 * min_edge / max_diag ;
      }


      /**
       * I don't know what to do for this metric.
       * Maybe the base class knows...
       */
    default:
      return Elem::quality(q);
    }
}



std::pair<Real, Real> InfHex::qual_bounds (const ElemQuality) const
{
  libmesh_not_implemented();

  std::pair<Real, Real> bounds;

  /*
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

    case SHEAR:
    case SHAPE:
    bounds.first  = 0.3;
    bounds.second = 1.;
    break;

    case CONDITION:
    bounds.first  = 1.;
    bounds.second = 8.;
    break;

    case JACOBIAN:
    bounds.first  = 0.5;
    bounds.second = 1.;
    break;

    case DISTORTION:
    bounds.first  = 0.6;
    bounds.second = 1.;
    break;

    case TAPER:
    bounds.first  = 0.;
    bounds.second = 0.4;
    break;

    case STRETCH:
    bounds.first  = 0.25;
    bounds.second = 1.;
    break;

    case DIAGONAL:
    bounds.first  = 0.65;
    bounds.second = 1.;
    break;

    case SIZE:
    bounds.first  = 0.5;
    bounds.second = 1.;
    break;

    default:
    libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
    bounds.first  = -1;
    bounds.second = -1;
    }
  */

  return bounds;
}





const unsigned short int InfHex::_second_order_adjacent_vertices[8][2] =
  {
    { 0,  1}, // vertices adjacent to node 8
    { 1,  2}, // vertices adjacent to node 9
    { 2,  3}, // vertices adjacent to node 10
    { 0,  3}, // vertices adjacent to node 11

    { 4,  5}, // vertices adjacent to node 12
    { 5,  6}, // vertices adjacent to node 13
    { 6,  7}, // vertices adjacent to node 14
    { 4,  7}  // vertices adjacent to node 15
  };



const unsigned short int InfHex::_second_order_vertex_child_number[18] =
  {
    99,99,99,99,99,99,99,99, // Vertices
    0,1,2,0,                 // Edges
    0,1,2,0,0,               // Faces
    0                        // Interior
  };



const unsigned short int InfHex::_second_order_vertex_child_index[18] =
  {
    99,99,99,99,99,99,99,99, // Vertices
    1,2,3,3,                 // Edges
    5,6,7,7,2,               // Faces
    6                        // Interior
  };


bool InfHex::contains_point (const Point & p, Real tol) const
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
  // fall back to the FEInterface::inverse_map()
  const Point my_origin (this->origin());

  // determine the minimal distance of the base from the origin
  // Use norm_sq() instead of norm(), it is faster
  Point pt0_o(this->point(0) - my_origin);
  Point pt1_o(this->point(1) - my_origin);
  Point pt2_o(this->point(2) - my_origin);
  Point pt3_o(this->point(3) - my_origin);
  const Real min_distance_sq = std::min(pt0_o.norm_sq(),
                                        std::min(pt1_o.norm_sq(),
                                                 std::min(pt2_o.norm_sq(),
                                                          pt3_o.norm_sq())));

  // work with 1% allowable deviation.  We can still fall
  // back to the InfFE::inverse_map()
  const Real conservative_p_dist_sq = 1.01 * (Point(p - my_origin).norm_sq());



  if (conservative_p_dist_sq < min_distance_sq)
    {
      // the physical point is definitely not contained in the element
      return false;
    }

  // this captures the case that the point is not (almost) in the direction of the element.:
  // first, project the problem onto the unit sphere:
  Point p_o(p - my_origin);
  pt0_o /= pt0_o.norm();
  pt1_o /= pt1_o.norm();
  pt2_o /= pt2_o.norm();
  pt3_o /= pt3_o.norm();
  p_o /= p_o.norm();


  // now, check if it is in the projected face; using that the diagonal contains
  // the largest distance between points in it
  Real max_h = std::max((pt0_o - pt2_o).norm_sq(),
                        (pt1_o - pt2_o).norm_sq())*1.01;

  if ((p_o - pt0_o).norm_sq() > max_h ||
      (p_o - pt1_o).norm_sq() > max_h ||
      (p_o - pt2_o).norm_sq() > max_h ||
      (p_o - pt3_o).norm_sq() > max_h )
    {
      // the physical point is definitely not contained in the element
      return false;
    }

  // Declare a basic FEType.  Will use default in the base,
  // and something else (not important) in radial direction.
  FEType fe_type(default_order());

  const Point mapped_point = FEInterface::inverse_map(dim(),
                                                      fe_type,
                                                      this,
                                                      p,
                                                      tol,
                                                      false);

  return FEInterface::on_reference_element(mapped_point, this->type(), tol);
}

} // namespace libMesh




#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
