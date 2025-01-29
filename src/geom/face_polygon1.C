// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/face_polygon1.h"

#include "libmesh/edge_edge2.h"
#include "libmesh/enum_order.h"
#include "libmesh/side.h"
#include "libmesh/tensor_value.h"

// C++ headers
#include <numeric> // std::iota

namespace libMesh
{


Polygon1::Polygon1 (const unsigned int num_sides, Elem * p) :
  Polygon(num_sides, num_sides, p)
{
  // A default triangulation is better than nothing
  for (int i : make_range(num_sides-2))
    this->_triangulation.push_back({0, i+1, i+2});
}


unsigned int Polygon1::opposite_node(const unsigned int node_in,
                                     const unsigned int side_in) const
{
  const auto ns = this->n_sides();
  if (ns % 2)
    libmesh_error();

  libmesh_assert_less (node_in, ns);
  libmesh_assert_less (node_in, this->n_nodes());
  libmesh_assert_less (side_in, this->n_sides());
  libmesh_assert(this->is_node_on_side(node_in, side_in));

  if (node_in == side_in)
    return (node_in + ns/2 + 1) % ns;

  libmesh_assert(node_in == side_in + 1);
  return (node_in + ns/2 - 1) % ns;
}



// ------------------------------------------------------------
// Polygon1 class member functions

bool Polygon1::is_vertex(const unsigned int libmesh_dbg_var(i)) const
{
  libmesh_assert (i < this->n_nodes());

  return true;
}

bool Polygon1::is_edge(const unsigned int libmesh_dbg_var(i)) const
{
  libmesh_assert (i < this->n_nodes());

  return false;
}

bool Polygon1::is_face(const unsigned int libmesh_dbg_var(i)) const
{
  libmesh_assert (i < this->n_nodes());

  return false;
}

bool Polygon1::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  const auto ns = this->n_sides();
  libmesh_assert_less (s, ns);
  libmesh_assert_less (n, this->n_nodes());

  return ((n % ns) == s) ||
         ((n < ns) && ((s+1)%ns) == n);
}

std::vector<unsigned int>
Polygon1::nodes_on_side(const unsigned int s) const
{
  const auto ns = this->n_sides();
  libmesh_assert(!(this->n_nodes() % ns));
  libmesh_assert_less(s, ns);

  std::vector<unsigned int> returnval(2);
  returnval[0] = s;
  returnval[1] = (s+1)%ns;

  return returnval;
}

std::vector<unsigned int>
Polygon1::nodes_on_edge(const unsigned int e) const
{
  return this->nodes_on_side(e);
}

bool Polygon1::has_affine_map() const
{
  const unsigned int ns = this->n_sides();

  // Get a good tolerance to use
  Real perimeter_l1 = 0;
  for (auto s : make_range(ns))
    perimeter_l1 += (this->point((s+1)%ns) - this->point(s)).l1_norm();
  const Real tol = perimeter_l1 * affine_tol;

  // Find the map we expect to have
  const Point veci = this->point(1) - this->point(0);
  const Point vec12 = this->point(2) - this->point(1);

  // Exterior angle
  const Real theta = 2 * libMesh::pi / ns;
  const Real costheta = cos(theta);
  const Real sintheta = sin(theta);

  RealTensor map;
  map(0, 0) = veci(0);
  map(1, 0) = veci(1);
  map(2, 0) = veci(2);
  map(0, 1) = (vec12(0) - veci(0)*costheta)/sintheta;
  map(1, 1) = (vec12(1) - veci(1)*costheta)/sintheta;
  map(2, 1) = (vec12(2) - veci(2)*costheta)/sintheta;

  libmesh_assert_less((map * this->master_point(1) -
                       (this->point(1) - this->point(0))).l1_norm(),
                      tol);
  libmesh_assert_less((map * this->master_point(2) -
                       (this->point(2) - this->point(0))).l1_norm(),
                      tol);

  for (auto i : make_range(3u, ns))
    if ((map * this->master_point(i) -
                     (this->point(i) - this->point(0))).l1_norm() >
        tol)
      return false;

  return true;
}



Order Polygon1::default_order() const
{
  return FIRST;
}



std::unique_ptr<Elem> Polygon1::build_side_ptr (const unsigned int i,
                                                bool proxy)
{
  const auto ns = this->n_sides();
  libmesh_assert_less (i, ns);

  std::unique_ptr<Elem> sidep;
  if (proxy)
    {
#ifdef LIBMESH_ENABLE_DEPRECATED
      sidep = std::make_unique<Side<Edge2,Polygon1>>(this,i);
      libmesh_deprecated();
#else
      libmesh_error();
#endif
    }
  else
    {
      sidep = std::make_unique<Edge2>(this);
      sidep->set_node(0) = this->node_ptr(i);
      sidep->set_node(1) = this->node_ptr((i+1)%ns);
    }

#ifdef LIBMESH_ENABLE_DEPRECATED
  if (!proxy) // proxy sides used to leave parent() set
#endif
    sidep->set_parent(nullptr);
  sidep->set_interior_parent(this);

  sidep->set_mapping_type(this->mapping_type());
  sidep->subdomain_id() = this->subdomain_id();
#ifdef LIBMESH_ENABLE_AMR
  sidep->set_p_level(this->p_level());
#endif

  return sidep;
}



void Polygon1::build_side_ptr (std::unique_ptr<Elem> & side,
                               const unsigned int i)
{
  const auto ns = this->n_sides();
  libmesh_assert_less (i, ns);

  if (!side.get() || side->type() != EDGE2)
    {
      side = this->build_side_ptr(i, false);
    }
  else
    {
      side->subdomain_id() = this->subdomain_id();
      side->set_mapping_type(this->mapping_type());
#ifdef LIBMESH_ENABLE_AMR
      side->set_p_level(this->p_level());
#endif

      side->set_node(0) = this->node_ptr(i);
      side->set_node(1) = this->node_ptr((i+1)%ns);
    }
}



void Polygon1::connectivity(const unsigned int /*sf*/,
                            const IOPackage /*iop*/,
                            std::vector<dof_id_type> & /*conn*/) const
{
  libmesh_not_implemented();
}



Real Polygon1::volume () const
{
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::volume();

  // We use a triangulation to calculate here; assert that it's as
  // consistent as possible.
  libmesh_assert_equal_to (this->_triangulation.size(),
                           this->n_nodes() - 2);

  Real double_area = 0;
  for (const auto & triangle : this->_triangulation)
    {
      Point v01 = this->point(triangle[1]) -
                  this->point(triangle[0]);
      Point v02 = this->point(triangle[2]) -
                  this->point(triangle[0]);
      double_area += std::sqrt(cross_norm_sq(v01, v02));
    }

  return double_area/2;
}



Point Polygon1::true_centroid () const
{
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::true_centroid();

  // We use a triangulation to calculate here; assert that it's as
  // consistent as possible.
  libmesh_assert_equal_to (this->_triangulation.size(),
                           this->n_nodes() - 2);

  Real double_area = 0;
  Point double_area_weighted_centroid;
  for (const auto & triangle : this->_triangulation)
    {
      Point v01 = this->point(triangle[1]) -
                  this->point(triangle[0]);
      Point v02 = this->point(triangle[2]) -
                  this->point(triangle[0]);

      const Real  double_tri_area = std::sqrt(cross_norm_sq(v01, v02));
      const Point tri_centroid = (this->point(triangle[0]) +
                                  this->point(triangle[1]) +
                                  this->point(triangle[2]))/3;

      double_area += double_tri_area;

      double_area_weighted_centroid += double_tri_area * tri_centroid;
    }

  return double_area_weighted_centroid / double_area;
}



std::pair<unsigned short int, unsigned short int>
Polygon1::second_order_child_vertex (const unsigned int /*n*/) const
{
  libmesh_not_implemented();
  return std::pair<unsigned short int, unsigned short int> (0,0);
}


void Polygon1::permute(unsigned int perm_num)
{
  const auto ns = this->n_sides();
  libmesh_assert_less (perm_num, ns);

  // This is mostly for unit testing so I'll just make it O(N).
  for (unsigned int p : make_range(perm_num))
    {
      libmesh_ignore(p);

      Node * tempnode = this->node_ptr(0);
      Elem * tempneigh = this->neighbor_ptr(0);
      for (unsigned int s : make_range(ns-1))
        {
          this->set_node(s) = this->node_ptr((s+1)%ns);
          this->set_neighbor(s, this->neighbor_ptr(s+1));
        }
      this->set_node(ns-1) = tempnode;
      this->set_neighbor(ns-1, tempneigh);

      // Keep the same triangulation, just with permuted node order,
      // so we can really expect this to act like the *same* polygon
      for (auto & triangle : this->_triangulation)
        for (auto i : make_range(3))
          triangle[i] = (triangle[i]+1)%ns;
    }
}


void Polygon1::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  const auto ns = this->n_sides();

  // Reorder nodes in such a way as to keep side ns-1 the same
  for (auto s : make_range(0u, ns/2))
    {
      swap2nodes(s,ns-1-s);
      swap2neighbors(s,ns-2-s);
      swap2boundarysides(s,ns-2-s,boundary_info);
      swap2boundaryedges(s,ns-2-s,boundary_info);
    }
}



ElemType Polygon1::side_type (const unsigned int libmesh_dbg_var(s)) const
{
  libmesh_assert_less (s, this->n_sides());
  return EDGE2;
}



void Polygon1::retriangulate()
{
  this->_triangulation.clear();

  // Start with the full polygon
  std::vector<int> remaining_nodes(this->n_nodes());
  std::iota(remaining_nodes.begin(), remaining_nodes.end(), 0);

  const auto ns = this->n_sides();

  const Point vavg = this->vertex_average();

  // Find out what plane we're in, if we're in 3D
#if LIBMESH_DIM > 2
  Point plane_normal;
  for (auto i : make_range(ns))
    {
      const Point vi     = this->point(i) - vavg;
      const Point viplus = this->point((i+1)%ns) - vavg;
      plane_normal += vi.cross(viplus);
    }
  plane_normal = plane_normal.unit();
#endif

  // Greedy heuristic: find the node with the most acutely convex
  // angle and strip it out as a triangle.
  while (remaining_nodes.size() > 2)
    {
      Real min_cos_angle = 1;
      int best_vertex = -1;
      for (auto n : index_range(remaining_nodes))
        {
          const Point & pn = this->point(remaining_nodes[n]);
          const Point & pnext = this->point(remaining_nodes[(n+1)%ns]);
          const Point & pprev = this->point(remaining_nodes[(n+ns-1)%ns]);
          const Point vprev = (pn - pprev).unit();
          const Point vnext = (pnext - pn).unit();

          const Real sign_check = (vprev.cross(vnext)) * plane_normal;
          if (sign_check <= 0)
            continue;

          const Real cos_angle = vprev * vnext;
          if (cos_angle < min_cos_angle)
            {
              min_cos_angle = cos_angle;
              best_vertex = n;
            }
        }

      libmesh_assert(best_vertex >= 0);

      this->_triangulation.push_back({remaining_nodes[(best_vertex+ns-1)%ns],
                                      remaining_nodes[best_vertex],
                                      remaining_nodes[(best_vertex+1)%ns]});
      remaining_nodes.erase(remaining_nodes.begin()+best_vertex);
    }
}



} // namespace libMesh
