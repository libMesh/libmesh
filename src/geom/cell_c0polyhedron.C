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
#include "libmesh/cell_c0polyhedron.h"

#include "libmesh/enum_order.h"
#include "libmesh/face_polygon.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/tensor_value.h"

// C++ headers
#include <numeric> // std::iota

namespace libMesh
{

// ------------------------------------------------------------
// C0Polyhedron class member functions

C0Polyhedron::C0Polyhedron
  (const std::vector<std::shared_ptr<Polygon>> & sides, Elem * p) :
  Polyhedron(sides, p)
{
  this->retriangulate();
}



C0Polyhedron::~C0Polyhedron() = default;



std::unique_ptr<Elem> C0Polyhedron::disconnected_clone() const
{
  auto sides = this->side_clones();

  std::unique_ptr<Elem> returnval =
    std::make_unique<C0Polyhedron>(sides);

  returnval->set_id() = this->id();
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (this->valid_unique_id())
    returnval->set_unique_id(this->unique_id());
#endif

  const auto n_elem_ints = this->n_extra_integers();
  returnval->add_extra_integers(n_elem_ints);
  for (unsigned int i = 0; i != n_elem_ints; ++i)
    returnval->set_extra_integer(i, this->get_extra_integer(i));

  returnval->inherit_data_from(*this);

  return returnval;
}



bool C0Polyhedron::is_vertex(const unsigned int libmesh_dbg_var(i)) const
{
  libmesh_assert (i < this->n_nodes());

  return true;
}



bool C0Polyhedron::is_edge(const unsigned int libmesh_dbg_var(i)) const
{
  libmesh_assert_less (i, this->n_nodes());

  return false;
}



bool C0Polyhedron::is_face(const unsigned int libmesh_dbg_var(i)) const
{
  libmesh_assert_less (i, this->n_nodes());

  return false;
}



bool C0Polyhedron::is_node_on_side(const unsigned int n,
                                   const unsigned int s) const
{
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (s, this->n_sides());

  const std::vector<unsigned int> & node_map =
    std::get<2>(_sidelinks_data[s]);

  const auto it = std::find(node_map.begin(), node_map.end(), n);

  return (it != node_map.end());
}

std::vector<unsigned int>
C0Polyhedron::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, this->n_sides());
  return std::get<2>(_sidelinks_data[s]);
}

std::vector<unsigned int>
C0Polyhedron::nodes_on_edge(const unsigned int e) const
{
  auto [s, se] = _edge_lookup[e];
  const Polygon & face = *std::get<0>(_sidelinks_data[s]);
  const std::vector<unsigned int> & node_map =
    std::get<2>(_sidelinks_data[s]);
  std::vector<unsigned int> nodes_on_edge =
    face.nodes_on_side(se);
  for (auto i : index_range(nodes_on_edge))
    nodes_on_edge[i] = node_map[nodes_on_edge[i]];
  return nodes_on_edge;
}



bool C0Polyhedron::is_node_on_edge(const unsigned int n,
                                   const unsigned int e) const
{
  libmesh_assert_less(e, _edge_lookup.size());
  auto [s, se] = _edge_lookup[e];

  const Polygon & face = *std::get<0>(_sidelinks_data[s]);
  const std::vector<unsigned int> & node_map =
    std::get<2>(_sidelinks_data[s]);
  std::vector<unsigned int> nodes_on_edge =
    face.nodes_on_side(se);
  for (auto noe : nodes_on_edge)
    if (node_map[noe] == n)
      return true;

  return false;
}



void C0Polyhedron::connectivity(const unsigned int /*sf*/,
                                const IOPackage /*iop*/,
                                std::vector<dof_id_type> & /*conn*/) const
{
  libmesh_not_implemented();
}



Real C0Polyhedron::volume () const
{
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::volume();

  // We use a triangulation to calculate here.

  Real six_vol = 0;
  for (const auto & subtet : this->_triangulation)
    {
      const Point p0 = this->point(subtet[0]);
      const Point p1 = this->point(subtet[1]);
      const Point p2 = this->point(subtet[2]);
      const Point p3 = this->point(subtet[3]);

      const Point v01 = p1 - p0;
      const Point v02 = p2 - p0;
      const Point v03 = p3 - p0;

      six_vol += triple_product(v01, v02, v03);
    }

  return six_vol * (1./6.);
}



Point C0Polyhedron::true_centroid () const
{
  // This specialization is good for Lagrange mappings only
  if (this->mapping_type() != LAGRANGE_MAP)
    return this->Elem::true_centroid();

  Real six_vol = 0;
  Point six_vol_weighted_centroid;
  for (const auto & subtet : this->_triangulation)
    {
      const Point p0 = this->point(subtet[0]);
      const Point p1 = this->point(subtet[1]);
      const Point p2 = this->point(subtet[2]);
      const Point p3 = this->point(subtet[3]);

      const Point v01 = p1 - p0;
      const Point v02 = p2 - p0;
      const Point v03 = p3 - p0;

      const Real six_tet_vol = triple_product(v01, v02, v03);

      const Point tet_centroid = (p0 + p1 + p2 + p3) * 0.25;

      six_vol += six_tet_vol;

      six_vol_weighted_centroid += six_tet_vol * tet_centroid;
    }

  return six_vol_weighted_centroid / six_vol;
}



std::pair<unsigned short int, unsigned short int>
C0Polyhedron::second_order_child_vertex (const unsigned int /*n*/) const
{
  libmesh_not_implemented();
  return std::pair<unsigned short int, unsigned short int> (0,0);
}



ElemType C0Polyhedron::side_type (const unsigned int libmesh_dbg_var(s)) const
{
  libmesh_assert_less (s, this->n_sides());
  return C0POLYGON;
}



Point
C0Polyhedron::side_vertex_average_normal(const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());
  libmesh_assert_equal_to(this->mapping_type(), LAGRANGE_MAP);
  const auto poly_side_ptr = this->side_ptr(s);
  const auto n_side_edges = poly_side_ptr->n_sides();

  // At the side vertex average, things simplify a bit
  // We get the side "plane" normal at all vertices, then average them
  Point normal;
  Point current_edge = poly_side_ptr->point(1) - poly_side_ptr->point(0);
  for (auto i : make_range(n_side_edges))
  {
    const Point next_edge = poly_side_ptr->point((i + 2) % n_side_edges) -
                            poly_side_ptr->point((i + 1) % n_side_edges);
    const Point normal_at_vertex = current_edge.cross(next_edge);
    normal += normal_at_vertex;
    // Note: the sides are planar, we don't need to test them all
    if (normal.norm_sq() > TOLERANCE)
      break;
    current_edge = next_edge;
  }
  bool outward_normal = std::get<1>(_sidelinks_data[s]);
  return (outward_normal ? 1. : -1.) * normal.unit();
}



void C0Polyhedron::retriangulate()
{
  this->_triangulation.clear();

  // We should already have a triangulation for each side.  We'll turn
  // that into a triangulation of the entire polyhedral surface
  // (stored as a mesh, so we can use `find_neighbors()`, then turn
  // that into a tetrahedralization by shrinking the volume one tet at
  // a time, which should work fine for convex polyhedra.

  Parallel::Communicator comm_self; // the default communicator is 1 rank
  ReplicatedMesh surface(comm_self);

  // poly_node_to_local_id[poly_node] is the local id of
  // poly_node in the Polyhedron, which is also the global id of the
  // corresponding node in the surface mesh.
  std::unordered_map<Node *, unsigned int> poly_node_to_id;

  for (unsigned int s : make_range(this->n_sides()))
    {
      const auto & [side, inward_normal, node_map] = this->_sidelinks_data[s];

      for (auto t : make_range(side->n_subtriangles()))
        {
          Elem * tri = surface.add_elem(Elem::build(TRI3));

          const std::array<int, 3> subtri = side->subtriangle(t);

          for (int i : make_range(3))
            {
              const int side_id = subtri[i];
              const Node * poly_node = side->node_ptr(side_id);

              libmesh_assert_less(side_id, node_map.size());
              const unsigned int local_id = node_map[side_id];

              Node * surf_node = surface.query_node_ptr(local_id);
              if (surf_node)
                libmesh_assert_equal_to(*(const Point*)poly_node,
                                        *(const Point*)surf_node);
              else
                surf_node = surface.add_point(*poly_node, local_id);

              // Get a consistent orientation for surface triangles,
              // facing their zeta directions outward
              const int tri_node = inward_normal ? i : 2-i;
              tri->set_node(tri_node, surf_node);
            }
        }
    }

  surface.allow_renumbering(false);
  surface.prepare_for_use();

  // We should have a watertight surface with consistent triangle
  // orientations.  That's expensive to check.
#ifdef DEBUG
  auto verify_surface = [& surface] ()
    {
      for (const Elem * elem : surface.element_ptr_range())
        {
          for (auto s : make_range(3))
            {
              const Elem * neigh = elem->neighbor_ptr(s);
              libmesh_assert(neigh);
              libmesh_assert_equal_to(neigh,
                                      surface.elem_ptr(neigh->id()));
              const unsigned int ns = neigh->which_neighbor_am_i(elem);
              libmesh_assert_less(ns, 3);
              libmesh_assert_equal_to(elem->node_ptr(s),
                                      neigh->node_ptr((ns+1)%3));
              libmesh_assert_equal_to(elem->node_ptr((s+1)%3),
                                      neigh->node_ptr(ns));
            }
        }
    };

  verify_surface();
#endif

  // We'll have to edit this as we change the surface elements, but we
  // have a method to initialize it easily.
  std::vector<std::vector<dof_id_type>> nodes_to_elem_vec_map;
  MeshTools::build_nodes_to_elem_map(surface, nodes_to_elem_vec_map);

  // We'll be inserting and deleting entries heavily, so we'll use
  // sets rather than vectors.  We want to get the same results in
  // parallel, so we'll use element ids rather than Elem pointers
  std::vector<std::set<dof_id_type>> nodes_to_elem_map;
  for (auto i : index_range(nodes_to_elem_vec_map))
    nodes_to_elem_map.emplace_back
      (nodes_to_elem_vec_map[i].begin(),
       nodes_to_elem_vec_map[i].end());

  // Now start meshing the volume enclosed by surface, one tet at a
  // time, with a greedy heuristic: find the vertex node with the most
  // acutely convex (solid) angle and strip it out as tetrahedra.

  // We'll want a vector of surrounding nodes for multiple uses,
  // sometimes with a similarly-sorted vector of surrounding elements
  // to go with it.
  auto surroundings_of =
    [&nodes_to_elem_map, & surface]
    (const Node & node,
     std::vector<Elem *> * surrounding_elems)
    {
      const std::set<dof_id_type> & elems_by_node =
        nodes_to_elem_map[node.id()];

      const unsigned int n_surrounding = elems_by_node.size();
      libmesh_assert_greater_equal(n_surrounding, 3);

      if (surrounding_elems)
        {
          libmesh_assert(surrounding_elems->empty());
          surrounding_elems->resize(n_surrounding);
        }

      std::vector<Node *> surrounding_nodes(n_surrounding);

      Elem * elem = surface.elem_ptr(*elems_by_node.begin());
      for (auto i : make_range(n_surrounding))
        {
          const unsigned int n = elem->get_node_index(&node);
          libmesh_assert_not_equal_to(n, invalid_uint);
          Node * next_node = elem->node_ptr((n+1)%3);
          surrounding_nodes[i] = next_node;
          if (surrounding_elems)
            (*surrounding_elems)[i] = elem;
          elem = elem->neighbor_ptr((n+2)%3);
          libmesh_assert(elem);
          libmesh_assert_equal_to(elem, surface.elem_ptr(elem->id()));

          // We should have a manifold here, but verifying that is
          // expensive
#ifdef DEBUG
          libmesh_assert_equal_to
            (std::count(surrounding_nodes.begin(),
                        surrounding_nodes.end(), next_node),
             1);
#endif
        }

      // We should have finished a loop
      libmesh_assert_equal_to
        (elem, surface.elem_ptr(*elems_by_node.begin()));

      return surrounding_nodes;
    };

  auto geometry_at = [&surroundings_of](const Node & node)
    {
      const std::vector<Node *> surrounding_nodes =
        surroundings_of(node, nullptr);

      // Now sum up solid angles from tetrahedra created from the
      // trivial triangulation of the surrounding nodes loop
      Real total_solid_angle = 0;
      const int n_surrounding =
        cast_int<int>(surrounding_nodes.size());

      for (auto n : make_range(n_surrounding-2))
        {
          const Point
            v01 = static_cast<Point>(*surrounding_nodes[n]) - node,
            v02 = static_cast<Point>(*surrounding_nodes[n+1]) - node,
            v03 = static_cast<Point>(*surrounding_nodes[n+2]) - node;

          total_solid_angle += solid_angle(v01, v02, v03);
        }

      return std::make_pair(n_surrounding, total_solid_angle);
    };

  // We'll keep track of solid angles and node valences so we don't
  // waste time recomputing them when they haven't changed.  We need
  // to be able to search efficiently for the smallest angles of each
  // valence, but also search efficiently for a particular node to
  // remove and reinsert it when its connectivity changes.
  //
  // Since C++11 multimap has guaranteed that pairs with matching keys
  // are kept in insertion order, so we can use Node * for values even
  // in parallel.
  typedef std::multimap<std::pair<int, Real>, Node*> node_map_type;
  node_map_type nodes_by_geometry;
  std::map<Node *, node_map_type::iterator> node_index;

  for (auto node : surface.node_ptr_range())
    node_index[node] =
      nodes_by_geometry.emplace(geometry_at(*node), node);

  // In 3D, this will require nested loops: an outer loop to remove
  // each vertex, and an inner loop to remove multiple tetrahedra in
  // cases where the vertex has more than 3 neighboring triangles.

  // We'll be done when there are only three "unremoved" nodes left,
  // so they don't actually enclose any volume.
  for (auto i : make_range(nodes_by_geometry.size()-3))
    {
      auto geometry_it = nodes_by_geometry.begin();
      auto geometry_key = geometry_it->first;
      auto [valence, angle] = geometry_key;
      Node * node = geometry_it->second;
      libmesh_ignore(i);

      // If our lowest-valence nodes are all points of non-convexity,
      // skip to a higher valence.
      while (angle > 2*pi-TOLERANCE)
        {
          geometry_it =
            nodes_by_geometry.upper_bound
              (std::make_pair(valence, Real(100)));
          libmesh_assert(geometry_it != nodes_by_geometry.end());

          std::tie(geometry_key, node) = *geometry_it;
          std::tie(valence, angle) = geometry_key;
        }

      std::vector<Elem *> surrounding_elems;
      std::vector<Node *> surrounding_nodes =
        surroundings_of(*node, &surrounding_elems);

      const std::size_t n_surrounding = surrounding_nodes.size();

      // As we separate surrounding nodes from our center node, we'll
      // be marking them as nullptr; we still need to be able to find
      // predecessor and successor nodes in order afterward.
      auto find_valid_nodes_around =
        [n_surrounding, & surrounding_nodes]
        (unsigned int j)
      {
        unsigned int jnext = (j+1)%n_surrounding;
        while (!surrounding_nodes[jnext])
          jnext = (jnext+1)%n_surrounding;

        unsigned int jprev = (j+n_surrounding-1)%n_surrounding;
        while (!surrounding_nodes[jprev])
          jprev = (jprev+n_surrounding-1)%n_surrounding;

        return std::make_pair(jprev, jnext);
      };

      // We may have too many surrounding nodes to handle with
      // just one tet.  In that case we'll keep a cache of the
      // element qualities that we'd get by making a tet with the
      // edge from the center node to each surrounding node, so we
      // can build the best tets first.
      //
      // In the case where we just have 3 nodes, we'll just pretend
      // they all have the same positive quality, so we can still
      // search this vector.
      std::vector<Real> local_tet_quality(n_surrounding, 1);

      // From our center node with N surrounding nodes we can make N-2
      // tetrahedra.  The first N-3 each replace two surface tets with
      // two new surface tets.
      //
      // My first idea was to greedily pick nodes with the smallest
      // local (solid) angles to get the best quality.  This works in
      // 2D, but such a node can give a pancake tet in 3D.
      //
      // My second idea was to greedily pick nodes with the highest
      // prospective tet quality.  This works for the first tets, but
      // can leave a pancake tet behind.
      //
      // My third idea is to try to fix the lowest quality tets first,
      // by picking cases where they have higher quality neighbors,
      // and creating those neighbors so as to change them.

      auto find_new_tet_nodes =
        [& local_tet_quality, & find_valid_nodes_around]
        ()
      {
        unsigned int jbest = 0;
        auto [jminus, jplus] = find_valid_nodes_around(jbest);
        Real qneighbest = std::min(local_tet_quality[jminus],
                                   local_tet_quality[jplus]);
        for (auto j : make_range(std::size_t(1),
                                 local_tet_quality.size()))
          {
            // We don't want to build a bad tet
            if (local_tet_quality[j] <= 0)
              continue;

            std::tie(jminus, jplus) = find_valid_nodes_around(j);
            Real qneighj = std::min(local_tet_quality[jminus],
                                    local_tet_quality[jplus]);

            // We don't want to build a tet that can't fix a neighbor
            // if we can build one that can.
            if (qneighbest <= 0 &&
                qneighj > 0)
              continue;

            // We want to try for the best possible fix.
            if ((local_tet_quality[j] - qneighj) >
                (local_tet_quality[jbest] - qneighj))
              {
                jbest = j;
                qneighbest = qneighj;
              }
          }

        libmesh_error_msg_if
          (local_tet_quality[jbest] <= 0,
           "Cannot build non-singular non-inverted tet");

        std::tie(jminus, jplus) = find_valid_nodes_around(jbest);

        return std::make_tuple(jbest, jminus, jplus);
      };

      if (n_surrounding > 3)
        {
          // We'll be searching local_tet_quality even after tet
          // extraction disconnects us from some nodes; when we do we
          // don't want to get one.
          constexpr Real far_node = -1e6;

          // Vectors from the center node to each of its surrounding
          // nodes are helpful for calculating prospective tet
          // quality.
          std::vector<Point> v0s(n_surrounding);
          for (auto j : make_range(n_surrounding))
            v0s[j] = *(Point *)surrounding_nodes[j] - *node;

          // Find the tet quality we'd potentially get from each
          // possible choice of tet
          auto local_tet_quality_of =
            [& surrounding_nodes, & v0s, & find_valid_nodes_around]
            (unsigned int j)
          {
            auto [jminus, jplus] = find_valid_nodes_around(j);

            // Anything proportional to the ratio of volume to
            // total-edge-length-cubed should peak for perfect tets
            // but hit 0 for pancakes and slivers.

            const Real total_len =
              v0s[j].norm() + v0s[jminus].norm() + v0s[jplus].norm() +
              (*(Point *)surrounding_nodes[jplus] -
               *(Point *)surrounding_nodes[j]).norm() +
              (*(Point *)surrounding_nodes[j] -
               *(Point *)surrounding_nodes[jminus]).norm() +
              (*(Point *)surrounding_nodes[jminus] -
               *(Point *)surrounding_nodes[jplus]).norm();

            // Orientation here is tricky.  Think of the triple
            // product as (v1 cross v2) dot v3, with right hand rule.
            const Real six_vol =
              triple_product(v0s[jminus], v0s[jplus], v0s[j]);

            return six_vol / (total_len * total_len * total_len);
          };

          for (auto j : make_range(n_surrounding))
            local_tet_quality[j] = local_tet_quality_of(j);

          // If we have N surrounding nodes, we can make N tets and
          // that'll bring us back to the 3-surrounding-node case to
          // finish.
          for (auto t : make_range(n_surrounding-3))
            {
              libmesh_ignore(t);

              auto [jbest, jminus, jplus] = find_new_tet_nodes();

              // Turn these four nodes into a tet
              Node * nbest  = surrounding_nodes[jbest],
                   * nminus = surrounding_nodes[jminus],
                   * nplus  = surrounding_nodes[jplus];
              this->add_tet(nminus->id(), nbest->id(), nplus->id(),
                            node->id());

              // Replace the old two triangles around these nodes with the
              // proper two new triangles.
              Elem * oldtri1 = surrounding_elems[jminus],
                   * oldtri2 = surrounding_elems[jbest],
                   * newtri1 = surface.add_elem(Elem::build(TRI3)),
                   * newtri2 = surface.add_elem(Elem::build(TRI3));

              const unsigned int c1 = oldtri1->get_node_index(node),
                                 c2 = oldtri2->get_node_index(node);

              newtri1->set_node(0, node);
              newtri1->set_node(1, nminus);
              newtri1->set_node(2, nplus);

              surrounding_elems[jminus] = newtri1;

              Elem * neigh10 = oldtri1->neighbor_ptr(c1);
              Elem * neigh12 = oldtri2->neighbor_ptr((c2+2)%3);
              newtri1->set_neighbor(0, neigh10);
              neigh10->set_neighbor(neigh10->which_neighbor_am_i(oldtri1), newtri1);
              newtri1->set_neighbor(1, newtri2);
              newtri1->set_neighbor(2, neigh12);
              neigh12->set_neighbor(neigh12->which_neighbor_am_i(oldtri2), newtri1);

              newtri2->set_node(0, nplus);
              newtri2->set_node(1, nminus);
              newtri2->set_node(2, nbest);

              Elem * neigh21 = oldtri1->neighbor_ptr((c1+1)%3);
              Elem * neigh22 = oldtri2->neighbor_ptr((c2+1)%3);
              newtri2->set_neighbor(0, newtri1);
              newtri2->set_neighbor(1, neigh21);
              neigh21->set_neighbor(neigh21->which_neighbor_am_i(oldtri1), newtri2);
              newtri2->set_neighbor(2, neigh22);
              neigh22->set_neighbor(neigh22->which_neighbor_am_i(oldtri2), newtri2);

              for (unsigned int p : make_range(3))
                {
                  nodes_to_elem_map[oldtri1->node_id(p)].erase(oldtri1->id());
                  nodes_to_elem_map[oldtri2->node_id(p)].erase(oldtri2->id());
                  nodes_to_elem_map[newtri1->node_id(p)].insert(newtri1->id());
                  nodes_to_elem_map[newtri2->node_id(p)].insert(newtri2->id());
                }

              // We've finished replacing the old triangles.
              surface.delete_elem(oldtri1);
              surface.delete_elem(oldtri2);

              // The solid angle for the far node should now stay
              // unchanged until we're out of this inner loop; let's
              // recalculate it here, and then we'll be done with it.
              Node * & nbestref = surrounding_nodes[jbest];
              nodes_by_geometry.erase(node_index[nbestref]);
              node_index[nbestref] =
                nodes_by_geometry.emplace(geometry_at(*nbestref), nbestref);

              // The far node is no longer sharing an edge with our center
              // node.  Make sure we don't use it again with the center
              // node.
              local_tet_quality[jbest] = far_node;
              nbestref = nullptr;

              // The potential tet qualities using the side nodes have
              // changed now that they're directly connected to each
              // other.
              local_tet_quality[jminus] =
                local_tet_quality_of(jminus);

              local_tet_quality[jplus] =
                local_tet_quality_of(jplus);
            }
        }

      // Now we should have just 3 surrounding nodes, with which to
      // make one tetrahedron.  Put them in a counterclockwise
      // (looking from outside) orientation, not the "best, clockwise,
      // counterclockwise" we get from the lambda.
      auto [j2, j1, j3] = find_new_tet_nodes();

      // Turn these last four nodes into a tet
      Node * n1 = surrounding_nodes[j1],
           * n2 = surrounding_nodes[j2],
           * n3 = surrounding_nodes[j3];
      this->add_tet(n1->id(), n2->id(), n3->id(), node->id());

      // Replace the three surface triangles of that tet with the new
      // fourth triangle.
      Elem * oldtri1 = surrounding_elems[j1],
           * oldtri2 = surrounding_elems[j2],
           * oldtri3 = surrounding_elems[j3],
           * newtri = surface.add_elem(Elem::build(TRI3));

      const unsigned int c1 = oldtri1->get_node_index(node),
                         c2 = oldtri2->get_node_index(node),
                         c3 = oldtri3->get_node_index(node);

      newtri->set_node(0, n1);
      newtri->set_node(1, n2);
      newtri->set_node(2, n3);
      Elem * neigh0 = oldtri1->neighbor_ptr((c1+1)%3);
      newtri->set_neighbor(0, neigh0);
      neigh0->set_neighbor(neigh0->which_neighbor_am_i(oldtri1), newtri);
      Elem * neigh1 = oldtri2->neighbor_ptr((c2+1)%3);
      newtri->set_neighbor(1, neigh1);
      neigh1->set_neighbor(neigh1->which_neighbor_am_i(oldtri2), newtri);
      Elem * neigh2 = oldtri3->neighbor_ptr((c3+1)%3);
      newtri->set_neighbor(2, neigh2);
      neigh2->set_neighbor(neigh2->which_neighbor_am_i(oldtri3), newtri);

      for (unsigned int p : make_range(3))
        {
          nodes_to_elem_map[oldtri1->node_id(p)].erase(oldtri1->id());
          nodes_to_elem_map[oldtri2->node_id(p)].erase(oldtri2->id());
          nodes_to_elem_map[oldtri3->node_id(p)].erase(oldtri3->id());
          nodes_to_elem_map[newtri->node_id(p)].insert(newtri->id());
        }

      // We've finished replacing the old triangles.
      surface.delete_elem(oldtri1);
      surface.delete_elem(oldtri2);
      surface.delete_elem(oldtri3);

      // We should have used up all our surrounding nodes now, and we
      // shouldn't have messed up our surface in the process, and our
      // center node should no longer be part of the surface.
#ifdef DEBUG
      surrounding_nodes[j1] = nullptr;
      surrounding_nodes[j2] = nullptr;
      surrounding_nodes[j3] = nullptr;

      for (auto ltq : surrounding_nodes)
        libmesh_assert(!ltq);

      if (surface.n_elem() > 3)
        verify_surface();

      for (const Elem * elem : surface.element_ptr_range())
        for (auto p : make_range(3))
          libmesh_assert_not_equal_to
            (elem->node_ptr(p), node);
#endif

      // We've used up our center node, so it's not something we can
      // eliminate again.
      nodes_by_geometry.erase(geometry_it);
    }

  // At this point our surface should just have two triangles left.
  libmesh_assert_equal_to(surface.n_elem(), 2);
}


void C0Polyhedron::add_tet(int n1,
                           int n2,
                           int n3,
                           int n4)
{
#ifndef NDEBUG
  const auto nn = this->n_nodes();
  libmesh_assert_less(n1, nn);
  libmesh_assert_less(n2, nn);
  libmesh_assert_less(n3, nn);
  libmesh_assert_less(n4, nn);

  const Point v12 = this->point(n2) - this->point(n1);
  const Point v13 = this->point(n3) - this->point(n1);
  const Point v14 = this->point(n4) - this->point(n1);
  const Real six_vol = triple_product(v12, v13, v14);
  libmesh_assert_greater(six_vol, Real(0));
#endif

  this->_triangulation.push_back({n1, n2, n3, n4});
}


} // namespace libMesh
