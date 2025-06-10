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

#include "libmesh/cell_polyhedron.h"

// Local includes
#include "libmesh/face_polygon.h"
#include "libmesh/enum_elem_quality.h"
#include "libmesh/hashword.h"

// C++ includes
#include <array>
#include <unordered_map>
#include <unordered_set>


namespace libMesh
{

// ------------------------------------------------------------
// Polyhedron class static member initializations
const int Polyhedron::num_children;

// ------------------------------------------------------------
// Polyhedron class member functions


Polyhedron::Polyhedron (const std::vector<std::shared_ptr<Polygon>> & sides,
                        Elem * p) :
    Cell(/* unused here */ 0, sides.size(), p, nullptr, nullptr),
    _elemlinks_data(sides.size()+2), // neighbors + parent + interior_parent
    _nodelinks_data(0), // We'll have to resize *this* later too!
    _sidelinks_data(sides.size())
  {
    // Set our sides, and while we're at it figure out our node maps
    // and side normal directions and edge lookup table, and count our
    // sides' nodes.  If we have internal nodes too then the subclass
    // will append those afterward.
    unsigned int nn = 0;
    std::unordered_map<Node *, unsigned int> local_node_number;
    std::unordered_set<std::pair<const Node *, const Node *>,
                       libMesh::hash> edges_seen;
    std::unique_ptr<const Elem> edge;
    for (unsigned int s : index_range(sides))
      {
        libmesh_assert(sides[s].get());
        auto & side_tuple = _sidelinks_data[s];
        std::get<0>(side_tuple) = sides[s];

        Polygon & side = *sides[s]; // not const, for writeable nodes
        for (auto n : make_range(side.n_nodes()))
          {
            Node * node = side.node_ptr(n);
            if (auto it = local_node_number.find(node);
                it != local_node_number.end())
              {
                std::get<2>(side_tuple).push_back(it->second);
              }
            else
              {
                std::get<2>(side_tuple).push_back(nn);
                local_node_number[node] = nn++;
                _nodelinks_data.push_back(node);
              }
          }

        for (unsigned int e : make_range(side.n_edges()))
          {
            side.build_edge_ptr(edge, e);
            auto edge_vertices = std::make_pair(edge->node_ptr(0), edge->node_ptr(1));
            if (edge_vertices.first > edge_vertices.second)
              std::swap(edge_vertices.first, edge_vertices.second);

            if (!edges_seen.count(edge_vertices))
              {
                edges_seen.insert(edge_vertices);
                _edge_lookup.emplace_back(s, e);
              }
          }
      }

    // Do the manual initialization that Elem::Elem and Cell::Cell
    // couldn't, now that we've resized both our vectors.  No need to
    // manually set nullptr, though, since std::vector does that.
    this->_elemlinks = _elemlinks_data.data();
    this->_nodes = _nodelinks_data.data();
    this->_elemlinks[0] = p;

    libmesh_assert_equal_to(nn, this->n_nodes());

    // Figure out the orientation of our sides, now that we've got our
    // nodes organized enough to find our center.  The algorithm below
    // only works for convex polyhedra, but that's all we're
    // supporting for now.
    Point center;
    for (auto n : make_range(nn))
      center.add (this->point(n));
    center /= static_cast<Real>(nn);

    for (unsigned int s : index_range(sides))
      {
        const Polygon & side = *sides[s];
        const Point x_i = side.point(0);
        const Point n_i =
          (side.point(1) - side.point(0)).cross
          (side.point(0) - side.point(side.n_sides()-1)).unit();

        bool & inward_normal = std::get<1>(_sidelinks_data[s]);
        inward_normal = (n_i * (center - x_i) > TOLERANCE);
      }

    // We're betting a lot on "our polyhedra are all convex", so let's
    // check that if we have time.
#ifdef DEBUG
    for (unsigned int s : index_range(sides))
      {
        const Polygon & side = *sides[s];
        const Point x_i = side.point(0);
        const bool inward_normal = std::get<1>(this->_sidelinks_data[s]);

        const Point n_i =
          (side.point(1) - side.point(0)).cross
          (side.point(0) - side.point(side.n_sides()-1)).unit() *
          (inward_normal ? -1 : 1);

        for (const Point & node : this->node_ref_range())
          {
            const Point d_n = node - x_i;
            if (d_n * n_i > TOLERANCE * d_n.norm())
              libmesh_not_implemented_msg
                ("Cannot create a non-convex polyhedron");
          }
      }
#endif

    // Is this likely to ever be used?  We may do refinement with
    // polyhedra but it's probably not going to have a hierarchy...
    if (p)
      {
        this->subdomain_id() = p->subdomain_id();
        this->processor_id() = p->processor_id();
        _map_type = p->mapping_type();
        _map_data = p->mapping_data();

#ifdef LIBMESH_ENABLE_AMR
        this->set_p_level(p->p_level());
#endif
      }

    // Make sure the interior parent isn't undefined
    this->set_interior_parent(nullptr);
  }



Point Polyhedron::master_point (const unsigned int i) const
{
  return this->point(i);
}



bool Polyhedron::convex()
{
  for (unsigned int s : make_range(this->n_sides()))
    {
      const Polygon & side = *std::get<0>(this->_sidelinks_data[s]);
      const Point x_i = side.point(0);
      const bool inward_normal = std::get<1>(this->_sidelinks_data[s]);

      const Point n_i =
        (side.point(1) - side.point(0)).cross
        (side.point(0) - side.point(side.n_sides()-1)).unit() *
        (inward_normal ? -1 : 1);

      for (const Point & node : this->node_ref_range())
        {
          const Point d_n = node - x_i;
          if (d_n * n_i > TOLERANCE * d_n.norm())
            return false;
        }
    }
  return true;
}



bool Polyhedron::on_reference_element(const Point & p,
                                      const Real eps) const
{
  const unsigned int ns = this->n_sides();

  // Check that the point is on the same side of all the faces by
  // testing whether:
  //
  // n_i.(p - x_i) <= 0
  //
  // for each i, where:
  //   n_i is the outward normal of face i,
  //   x_i is a point on face i.

  for (auto i : make_range(ns))
    {
      const Polygon & face = *std::get<0>(this->_sidelinks_data[i]);
      const bool inward_normal = std::get<1>(this->_sidelinks_data[i]);

      const Point x_i = face.point(0);

      const Point n_i =
        (face.point(1) - face.point(0)).cross
        (face.point(0) - face.point(face.n_sides()-1)).unit() *
        (inward_normal ? -1 : 1);

  // This only works for polyhedra with flat sides.
#ifdef DEBUG
      for (auto j : make_range(face.n_sides()-1))
        {
          const Point x_j = face.point(j+1);
          const Point d_j = x_j - x_i;
          if (std::abs(d_j * n_i) > eps * d_j.norm())
            libmesh_not_implemented_msg
              ("Polyhedra with non-flat sides are not fully supported.");
        }
#endif

      if (n_i * (p - x_i) > eps)
        return false;
    }

  return true;
}



dof_id_type Polyhedron::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  const Polygon & face = *std::get<0>(this->_sidelinks_data[s]);

  return face.key();
}



dof_id_type Polyhedron::low_order_key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  const Polygon & face = *std::get<0>(this->_sidelinks_data[s]);

  const unsigned int nv = face.n_vertices();
  std::vector<dof_id_type> vertex_ids(nv);
  for (unsigned int v : make_range(nv))
    vertex_ids[v] = face.node_id(v);

  return Utility::hashword(vertex_ids);
}



unsigned int Polyhedron::local_side_node(unsigned int side,
                                         unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());

  const std::vector<unsigned int> & node_map =
    std::get<2>(this->_sidelinks_data[side]);
  libmesh_assert_less (side_node, node_map.size());

  return node_map[side_node];
}



unsigned int Polyhedron::local_edge_node(unsigned int edge,
                                         unsigned int edge_node) const
{
  libmesh_assert_less (edge, this->n_edges());
  libmesh_assert_less (edge, _edge_lookup.size());

  auto [side, edge_of_side] = _edge_lookup[edge];

  const Polygon & face = *std::get<0>(this->_sidelinks_data[side]);

  const std::vector<unsigned int> & node_map =
    std::get<2>(this->_sidelinks_data[side]);

  return node_map[face.local_edge_node(edge_of_side, edge_node)];
}



dof_id_type Polyhedron::key () const
{
  std::vector<dof_id_type> node_ids;
  for (const auto & n : this->node_ref_range())
    node_ids.push_back(n.id());

  return Utility::hashword(node_ids);
}



std::unique_ptr<Elem> Polyhedron::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  Polygon & face = *std::get<0>(this->_sidelinks_data[i]);
  std::unique_ptr<Elem> face_copy = face.disconnected_clone();
  for (auto n : face.node_index_range())
    face_copy->set_node(n, face.node_ptr(n));

  return face_copy;
}



void Polyhedron::side_ptr (std::unique_ptr<Elem> & side,
                           const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  // Polyhedra are irregular enough that we're not even going to try
  // and bother optimizing heap access here.
  side = this->side_ptr(i);
}



std::unique_ptr<Elem> Polyhedron::build_side_ptr (const unsigned int i)
{
  auto returnval = this->side_ptr(i);
  returnval->set_interior_parent(this);
  returnval->set_mapping_type(this->mapping_type());
  returnval->subdomain_id() = this->subdomain_id();
#ifdef LIBMESH_ENABLE_AMR
  returnval->set_p_level(this->p_level());
#endif
  return returnval;
}



void Polyhedron::build_side_ptr (std::unique_ptr<Elem> & side,
                                 const unsigned int i)
{
  this->side_ptr(side, i);
  side->set_interior_parent(this);
  side->set_mapping_type(this->mapping_type());
  side->subdomain_id() = this->subdomain_id();
#ifdef LIBMESH_ENABLE_AMR
  side->set_p_level(this->p_level());
#endif
}



std::unique_ptr<Elem> Polyhedron::build_edge_ptr (const unsigned int i)
{
  auto [s, se] = _edge_lookup[i];
  Polygon & face = *std::get<0>(_sidelinks_data[s]);
  return face.build_edge_ptr(se);
}



void Polyhedron::build_edge_ptr (std::unique_ptr<Elem> & elem,
                                 const unsigned int i)
{
  auto [s, se] = _edge_lookup[i];
  Polygon & face = *std::get<0>(_sidelinks_data[s]);
  face.build_edge_ptr(elem, se);
}



bool Polyhedron::is_child_on_side(const unsigned int /*c*/,
                                  const unsigned int /*s*/) const
{
  libmesh_not_implemented();
  return false;
}



unsigned int Polyhedron::opposite_side(const unsigned int /* side_in */) const
{
  // This is too ambiguous in general.
  libmesh_not_implemented();
  return libMesh::invalid_uint;
}



unsigned int Polyhedron::opposite_node(const unsigned int /* n */,
                                       const unsigned int /* s */) const
{
  // This is too ambiguous in general.
  libmesh_not_implemented();
  return libMesh::invalid_uint;
}



bool Polyhedron::is_flipped() const
{
  if (this->_triangulation.empty())
    return false;

  auto & tet = this->_triangulation[0];

  const Point v01 = this->point(tet[1]) - this->point(tet[0]);
  const Point v02 = this->point(tet[2]) - this->point(tet[0]);
  const Point v03 = this->point(tet[3]) - this->point(tet[0]);

  return (triple_product(v01, v02, v03) < 0);
}


std::vector<unsigned int>
Polyhedron::edges_adjacent_to_node(const unsigned int n) const
{
  libmesh_assert_less(n, this->n_nodes());

  // For mid-edge or mid-face nodes, the subclass had better have
  // overridden this.
  libmesh_assert_less(n, this->n_vertices());

  const unsigned int ne = this->n_edges();
  libmesh_assert_equal_to(ne, _edge_lookup.size());

  std::vector<unsigned int> adjacent_edges;

  unsigned int next_edge = 0;

  // Look for any adjacent edges on each side.  Make use of the fact
  // that we number our edges in order as we encounter them from each
  // side.
  for (auto t : index_range(_sidelinks_data))
    {
      const Polygon & face = *std::get<0>(_sidelinks_data[t]);
      const std::vector<unsigned int> & node_map =
        std::get<2>(_sidelinks_data[t]);

      while (_edge_lookup[next_edge].first < t)
        {
          ++next_edge;
          if (next_edge == ne)
            return adjacent_edges;
        }

      // If we haven't seen the next edge on this or an earlier side
      // then we might as well go on to the next.
      if (_edge_lookup[next_edge].first > t)
        continue;

      const unsigned int fnv = face.n_vertices();
      libmesh_assert_equal_to(fnv, face.n_edges());
      libmesh_assert_equal_to(fnv, face.n_sides());
      libmesh_assert_less_equal(fnv, node_map.size());

      // Polygon faces have one edge per vertex
      for (auto v : make_range(fnv))
        {
          libmesh_assert_equal_to (_edge_lookup[next_edge].first, t);

          if (_edge_lookup[next_edge].second > v)
            continue;

          while (_edge_lookup[next_edge].first == t &&
                 _edge_lookup[next_edge].second < v)
            {
              ++next_edge;
              if (next_edge == ne)
                return adjacent_edges;
            }

          if (_edge_lookup[next_edge].first > t)
            break;

          const unsigned int vn = node_map[v];
          const unsigned int vnp = node_map[(v+1)%fnv];
          libmesh_assert_less(vn, this->n_vertices());
          libmesh_assert_less(vnp, this->n_vertices());
          if (vn == n || vnp == n)
            adjacent_edges.push_back(next_edge);
        }
    }

  return adjacent_edges;
}


std::pair<Real, Real> Polyhedron::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {
    case EDGE_LENGTH_RATIO:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case MIN_ANGLE:
      bounds.first  = 30.;
      bounds.second = 180.;
      break;

    case MAX_ANGLE:
      bounds.first  = 60.;
      bounds.second = 180.;
      break;

    case JACOBIAN:
    case SCALED_JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.;
      break;

    default:
      libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}



std::vector<std::shared_ptr<Polygon>>
Polyhedron::side_clones() const
{
  const auto ns = this->n_sides();

  libmesh_assert_equal_to(ns, _sidelinks_data.size());

  std::vector<std::shared_ptr<Polygon>> cloned_sides(ns);

  for (auto i : make_range(ns))
    {
      const Polygon & face = *std::get<0>(this->_sidelinks_data[i]);

      Elem * clone = face.disconnected_clone().release();
      Polygon * polygon_clone = cast_ptr<Polygon *>(clone);
      cloned_sides[i] = std::shared_ptr<Polygon>(polygon_clone);

      // We can't actually use a *disconnected* clone to reconstruct
      // links between sides, so we'll temporarily give the clone our
      // own nodes; user code that typically replaces the usual
      // nullptr with permanent nodes will then instead place our
      // nodes with permanent nodes.
      for (auto n : make_range(face.n_nodes()))
        cloned_sides[i]->set_node
          (n, const_cast<Node *>(face.node_ptr(n)));
    }

  return cloned_sides;
}



bool Polyhedron::side_has_edge_nodes(unsigned int s,
                                     unsigned int min_node,
                                     unsigned int max_node) const
{
  const Polygon & face = *std::get<0>(_sidelinks_data[s]);
  const std::vector<unsigned int> & node_map =
    std::get<2>(this->_sidelinks_data[s]);

  for (unsigned int e : make_range(face.n_sides()))
    {
      std::vector<unsigned int> nodes_on_edge =
        face.nodes_on_side(e);
      libmesh_assert_equal_to(nodes_on_edge.size(), 2);
      nodes_on_edge[0] = node_map[nodes_on_edge[0]];
      nodes_on_edge[1] = node_map[nodes_on_edge[1]];
      if ((nodes_on_edge[0] == min_node) &&
          (nodes_on_edge[1] == max_node))
        return true;
      if ((nodes_on_edge[1] == min_node) &&
          (nodes_on_edge[0] == max_node))
        return true;
    }

  return false;
}



std::vector<unsigned int>
Polyhedron::sides_on_edge(const unsigned int e) const
{
  std::vector<unsigned int> returnval(2);
  auto [s1, s1e] = _edge_lookup[e];
  returnval[0] = s1;

  const Polygon & face1 = *std::get<0>(_sidelinks_data[s1]);
  const std::vector<unsigned int> & node_map =
    std::get<2>(this->_sidelinks_data[s1]);

  std::vector<unsigned int> nodes_on_edge =
    face1.nodes_on_side(s1e);
  libmesh_assert_equal_to(nodes_on_edge.size(), 2);
  nodes_on_edge[0] = node_map[nodes_on_edge[0]];
  nodes_on_edge[1] = node_map[nodes_on_edge[1]];

  if (nodes_on_edge[0] > nodes_on_edge[1])
    std::swap(nodes_on_edge[0], nodes_on_edge[1]);

  for (unsigned int s2 : make_range(this->n_sides()))
    {
      if (s2 == s1)
        continue;

      if (this->side_has_edge_nodes(s2, nodes_on_edge[0],
                                    nodes_on_edge[1]))
        {
          returnval[1] = s2;
          return returnval;
        }
    }

  libmesh_error();

  return returnval;
}



bool Polyhedron::is_edge_on_side(const unsigned int e,
                                 const unsigned int s) const
{
  auto [s1, s1e] = _edge_lookup[e];

  // Did we get lucky with our cache?
  if (s1 == s)
    return true;

  const Polygon & face1 = *std::get<0>(_sidelinks_data[s1]);
  const std::vector<unsigned int> & node_map =
    std::get<2>(this->_sidelinks_data[s1]);
  std::vector<unsigned int> nodes_on_edge1 =
    face1.nodes_on_side(s1e);
  libmesh_assert_equal_to(nodes_on_edge1.size(), 2);

  nodes_on_edge1[0] = node_map[nodes_on_edge1[0]];
  nodes_on_edge1[1] = node_map[nodes_on_edge1[1]];
  if (nodes_on_edge1[0] > nodes_on_edge1[1])
    std::swap(nodes_on_edge1[0], nodes_on_edge1[1]);

  return this->side_has_edge_nodes(s,
                                   nodes_on_edge1[0],
                                   nodes_on_edge1[1]);
}



std::array<Point, 4> Polyhedron::master_subelement (unsigned int i) const
{
  libmesh_assert_less(i, this->_triangulation.size());

  const auto & tet = this->_triangulation[i];

  return { this->master_point(tet[0]),
           this->master_point(tet[1]),
           this->master_point(tet[2]),
           this->master_point(tet[3]) };
}


std::tuple<unsigned int, Real, Real, Real>
Polyhedron::subelement_coordinates (const Point & p, Real tol) const
{
  std::tuple<unsigned int, Real, Real, Real> returnval =
    {libMesh::invalid_uint, -1, -1, -1};

  Real best_bad_coord = -1;

  for (auto s : make_range(this->n_subelements()))
    {
      const std::array<Point, 4> subtet =
        this->master_subelement(s);

      // Find barycentric coordinates in subelem
      const Point v0 = p - subtet[0];
      // const Point v1 = p - subtet[1];

      const Point v01 = subtet[1] - subtet[0];
      const Point v02 = subtet[2] - subtet[0];
      const Point v03 = subtet[3] - subtet[0];

      // const Point v12 = subtet[2] - subtet[1];
      // const Point v13 = subtet[3] - subtet[1];

      // const Real tp0 = triple_product(v1, v13, v12);
      const Real tp1 = triple_product(v0, v02, v03);
      const Real tp2 = triple_product(v0, v03, v01);
      const Real tp3 = triple_product(v0, v01, v02);

      const Real six_vol = triple_product(v01, v02, v03);

      const Real xi   = tp1 / six_vol;
      const Real eta  = tp2 / six_vol;
      const Real zeta = tp3 / six_vol;

      if (xi>=0 && eta>=0 && zeta>=0 && xi+eta+zeta<=1)
        return { s, xi, eta, zeta };

      const Real my_best_bad_coord =
        std::min(std::min(std::min(xi, eta), zeta), 1-xi-eta-zeta);

      if (my_best_bad_coord > best_bad_coord)
        {
          best_bad_coord = my_best_bad_coord;
          returnval = { s, xi, eta, zeta };
        }
    }

  if (best_bad_coord > -tol)
    return returnval;

  return {libMesh::invalid_uint, -1, -1, -1};
}


} // namespace libMesh
