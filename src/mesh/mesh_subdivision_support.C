// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// C++ includes

// Local includes
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_subdivision_support.h"
#include "libmesh/boundary_info.h"

namespace libMesh
{


void MeshTools::Subdivision::find_one_ring(const Tri3Subdivision * elem,
                                           std::vector<const Node *> & nodes)
{
  libmesh_assert(elem->is_subdivision_updated());
  libmesh_assert(elem->get_ordered_node(0));

  unsigned int valence = elem->get_ordered_valence(0);
  nodes.resize(valence + 6);

  // The first three vertices in the patch are the ones from the element triangle
  nodes[0]       = elem->get_ordered_node(0);
  nodes[1]       = elem->get_ordered_node(1);
  nodes[valence] = elem->get_ordered_node(2);

  const unsigned int nn0 = elem->local_node_number(nodes[0]->id());

  const Tri3Subdivision * nb = dynamic_cast<const Tri3Subdivision *>(elem->neighbor_ptr(nn0));
  libmesh_assert(nb);

  unsigned int j, i = 1;

  do
    {
      ++i;
      j = nb->local_node_number(nodes[0]->id());
      nodes[i] = nb->node_ptr(next[j]);
      nb = static_cast<const Tri3Subdivision *>(nb->neighbor_ptr(j));
    } while (nb != elem);

  /* for nodes connected with N (= valence[0]) */
  nb = static_cast<const Tri3Subdivision *>(elem->neighbor_ptr(next[nn0]));
  j = nb->local_node_number(nodes[1]->id());
  nodes[valence+1] = nb->node_ptr(next[j]);

  nb = static_cast<const Tri3Subdivision *>(nb->neighbor_ptr(next[j]));
  j = nb->local_node_number(nodes[valence+1]->id());
  nodes[valence+4] = nb->node_ptr(next[j]);

  nb = static_cast<const Tri3Subdivision *>(nb->neighbor_ptr(next[j]));
  j = nb->local_node_number(nodes[valence+4]->id());
  nodes[valence+5] = nb->node_ptr(next[j]);

  /* for nodes connected with 1 */
  nb = static_cast<const Tri3Subdivision *>(elem->neighbor_ptr(next[nn0]));
  j = nb->local_node_number(nodes[1]->id());
  // nodes[valence+1] has been determined already

  nb = static_cast<const Tri3Subdivision *>(nb->neighbor_ptr(j));
  j = nb->local_node_number(nodes[1]->id());
  nodes[valence+2] = nb->node_ptr(next[j]);

  nb = static_cast<const Tri3Subdivision *>(nb->neighbor_ptr(j));
  j = nb->local_node_number(nodes[1]->id());
  nodes[valence+3] = nb->node_ptr(next[j]);

  return;
}


void MeshTools::Subdivision::all_subdivision(MeshBase & mesh)
{
  std::vector<Elem *> new_elements;
  new_elements.reserve(mesh.n_elem());
  const bool mesh_has_boundary_data =
    (mesh.get_boundary_info().n_boundary_ids() > 0);

  std::vector<Elem *> new_boundary_elements;
  std::vector<short int> new_boundary_sides;
  std::vector<boundary_id_type> new_boundary_ids;

  // Container to catch ids handed back from BoundaryInfo
  std::vector<boundary_id_type> ids;

  MeshBase::const_element_iterator       el     = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();
  for (; el != end_el; ++el)
    {
      const Elem * elem = *el;
      libmesh_assert_equal_to(elem->type(), TRI3);

      Elem * tri = new Tri3Subdivision;
      tri->set_id(elem->id());
      tri->subdomain_id() = elem->subdomain_id();
      tri->set_node(0) = (*el)->node_ptr(0);
      tri->set_node(1) = (*el)->node_ptr(1);
      tri->set_node(2) = (*el)->node_ptr(2);

      if (mesh_has_boundary_data)
        {
          for (unsigned short side = 0; side < elem->n_sides(); ++side)
            {
              mesh.get_boundary_info().boundary_ids(elem, side, ids);

              for (std::size_t id=0; id<ids.size(); ++id)
                {
                  // add the boundary id to the list of new boundary ids
                  new_boundary_ids.push_back(ids[id]);
                  new_boundary_elements.push_back(tri);
                  new_boundary_sides.push_back(side);
                }
            }

          // remove the original element from the BoundaryInfo structure
          mesh.get_boundary_info().remove(elem);
        }

      new_elements.push_back(tri);
      mesh.insert_elem(tri);
    }
  mesh.prepare_for_use();

  if (mesh_has_boundary_data)
    {
      // If the old mesh had boundary data, the new mesh better have some too.
      libmesh_assert_greater(new_boundary_elements.size(), 0);

      // We should also be sure that the lengths of the new boundary data vectors
      // are all the same.
      libmesh_assert_equal_to(new_boundary_sides.size(), new_boundary_elements.size());
      libmesh_assert_equal_to(new_boundary_sides.size(), new_boundary_ids.size());

      // Add the new boundary info to the mesh.
      for (std::size_t s = 0; s < new_boundary_elements.size(); ++s)
        mesh.get_boundary_info().add_side(new_boundary_elements[s],
                                          new_boundary_sides[s],
                                          new_boundary_ids[s]);
    }

  mesh.prepare_for_use();
}


void MeshTools::Subdivision::prepare_subdivision_mesh(MeshBase & mesh, bool ghosted)
{
  mesh.prepare_for_use();

  // convert all mesh elements to subdivision elements
  all_subdivision(mesh);

  if (!ghosted)
    {
      // add the ghost elements for the boundaries
      add_boundary_ghosts(mesh);
    }
  else
    {
      // This assumes that the mesh already has the ghosts. Only tagging them is required here.
      tag_boundary_ghosts(mesh);
    }

  mesh.prepare_for_use();

  std::vector<std::vector<const Elem *> > nodes_to_elem_map;
  MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

  // compute the node valences
  MeshBase::const_node_iterator       nd     = mesh.nodes_begin();
  const MeshBase::const_node_iterator end_nd = mesh.nodes_end();
  for (; nd != end_nd; ++nd)
    {
      Node * node = *nd;
      std::vector<const Node *> neighbors;
      MeshTools::find_nodal_neighbors(mesh, *node, nodes_to_elem_map, neighbors);
      const unsigned int valence =
        cast_int<unsigned int>(neighbors.size());
      libmesh_assert_greater(valence, 1);
      node->set_valence(valence);
    }

  MeshBase::const_element_iterator       el     = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();
  for (; el != end_el; ++el)
    {
      Tri3Subdivision * elem = dynamic_cast<Tri3Subdivision *>(*el);
      libmesh_assert(elem);
      if (!elem->is_ghost())
        elem->prepare_subdivision_properties();
    }
}


void MeshTools::Subdivision::tag_boundary_ghosts(MeshBase & mesh)
{
  MeshBase::element_iterator       el     = mesh.elements_begin();
  const MeshBase::element_iterator end_el = mesh.elements_end();
  for (; el != end_el; ++el)
    {
      Elem * elem = *el;
      libmesh_assert_equal_to(elem->type(), TRI3SUBDIVISION);

      Tri3Subdivision * sd_elem = static_cast<Tri3Subdivision *>(elem);
      for (unsigned int i = 0; i < elem->n_sides(); ++i)
        {
          if (elem->neighbor_ptr(i) == libmesh_nullptr)
            {
              sd_elem->set_ghost(true);
              // set all other neighbors to ghosts as well
              if (elem->neighbor_ptr(next[i]))
                {
                  Tri3Subdivision * nb = static_cast<Tri3Subdivision *>(elem->neighbor_ptr(next[i]));
                  nb->set_ghost(true);
                }
              if (elem->neighbor_ptr(prev[i]))
                {
                  Tri3Subdivision * nb = static_cast<Tri3Subdivision *>(elem->neighbor_ptr(prev[i]));
                  nb->set_ghost(true);
                }
            }
        }
    }
}


void MeshTools::Subdivision::add_boundary_ghosts(MeshBase & mesh)
{
  static const Real tol = 1e-5;

  // add the mirrored ghost elements (without using iterators, because the mesh is modified in the course)
  std::vector<Tri3Subdivision *> ghost_elems;
  std::vector<Node *> ghost_nodes;
  const unsigned int n_elem = mesh.n_elem();
  for (unsigned int eid = 0; eid < n_elem; ++eid)
    {
      Elem * elem = mesh.elem_ptr(eid);
      libmesh_assert_equal_to(elem->type(), TRI3SUBDIVISION);

      // If the triangle happens to be in a corner (two boundary
      // edges), we perform a counter-clockwise loop by mirroring the
      // previous triangle until we come back to the original
      // triangle.  This prevents degenerated triangles in the mesh
      // corners and guarantees that the node in the middle of the
      // loop is of valence=6.
      for (unsigned int i = 0; i < elem->n_sides(); ++i)
        {
          libmesh_assert_not_equal_to(elem->neighbor_ptr(i), elem);

          if (elem->neighbor_ptr(i) == libmesh_nullptr &&
              elem->neighbor_ptr(next[i]) == libmesh_nullptr)
            {
              Elem * nelem = elem;
              unsigned int k = i;
              for (unsigned int l=0;l<4;l++)
                {
                  // this is the vertex to be mirrored
                  Point point = nelem->point(k) + nelem->point(next[k]) - nelem->point(prev[k]);

                  // Check if the proposed vertex doesn't coincide
                  // with one of the existing vertices.  This is
                  // necessary because for some triangulations, it can
                  // happen that two mirrored ghost vertices coincide,
                  // which would then lead to a zero size ghost
                  // element below.
                  Node * node = libmesh_nullptr;
                  for (std::size_t j = 0; j < ghost_nodes.size(); ++j)
                    {
                      if ((*ghost_nodes[j] - point).norm() < tol * (elem->point(k) - point).norm())
                        {
                          node = ghost_nodes[j];
                          break;
                        }
                    }

                  // add the new vertex only if no other is nearby
                  if (node == libmesh_nullptr)
                    {
                      node = mesh.add_point(point);
                      ghost_nodes.push_back(node);
                    }

                  Tri3Subdivision * newelem = new Tri3Subdivision();

                  // add the first new ghost element to the list just as in the non-corner case
                  if (l == 0)
                    ghost_elems.push_back(newelem);

                  newelem->set_node(0) = nelem->node_ptr(next[k]);
                  newelem->set_node(1) = nelem->node_ptr(k);
                  newelem->set_node(2) = node;
                  newelem->set_neighbor(0, nelem);
                  newelem->set_ghost(true);
                  if (l>0)
                    newelem->set_neighbor(2, libmesh_nullptr);
                  nelem->set_neighbor(k, newelem);

                  mesh.add_elem(newelem);
                  mesh.get_boundary_info().add_node(nelem->node_ptr(k), 1);
                  mesh.get_boundary_info().add_node(nelem->node_ptr(next[k]), 1);
                  mesh.get_boundary_info().add_node(nelem->node_ptr(prev[k]), 1);
                  mesh.get_boundary_info().add_node(node, 1);

                  nelem = newelem;
                  k = 2 ;
                }

              Tri3Subdivision * newelem = new Tri3Subdivision();

              newelem->set_node(0) = elem->node_ptr(next[i]);
              newelem->set_node(1) = nelem->node_ptr(2);
              newelem->set_node(2) = elem->node_ptr(prev[i]);
              newelem->set_neighbor(0, nelem);
              nelem->set_neighbor(2, newelem);
              newelem->set_ghost(true);
              newelem->set_neighbor(2, elem);
              elem->set_neighbor(next[i],newelem);

              mesh.add_elem(newelem);

              break;
            }
        }

      for (unsigned int i = 0; i < elem->n_sides(); ++i)
        {
          libmesh_assert_not_equal_to(elem->neighbor_ptr(i), elem);
          if (elem->neighbor_ptr(i) == libmesh_nullptr)
            {
              // this is the vertex to be mirrored
              Point point = elem->point(i) + elem->point(next[i]) - elem->point(prev[i]);

              // Check if the proposed vertex doesn't coincide with
              // one of the existing vertices.  This is necessary
              // because for some triangulations, it can happen that
              // two mirrored ghost vertices coincide, which would
              // then lead to a zero size ghost element below.
              Node * node = libmesh_nullptr;
              for (std::size_t j = 0; j < ghost_nodes.size(); ++j)
                {
                  if ((*ghost_nodes[j] - point).norm() < tol * (elem->point(i) - point).norm())
                    {
                      node = ghost_nodes[j];
                      break;
                    }
                }

              // add the new vertex only if no other is nearby
              if (node == libmesh_nullptr)
                {
                  node = mesh.add_point(point);
                  ghost_nodes.push_back(node);
                }

              Tri3Subdivision * newelem = new Tri3Subdivision();
              ghost_elems.push_back(newelem);

              newelem->set_node(0) = elem->node_ptr(next[i]);
              newelem->set_node(1) = elem->node_ptr(i);
              newelem->set_node(2) = node;
              newelem->set_neighbor(0, elem);
              newelem->set_ghost(true);
              elem->set_neighbor(i, newelem);

              mesh.add_elem(newelem);
              mesh.get_boundary_info().add_node(elem->node_ptr(i), 1);
              mesh.get_boundary_info().add_node(elem->node_ptr(next[i]), 1);
              mesh.get_boundary_info().add_node(elem->node_ptr(prev[i]), 1);
              mesh.get_boundary_info().add_node(node, 1);
            }
        }
    }

  // add the missing ghost elements (connecting new ghost nodes)
  std::vector<Tri3Subdivision *> missing_ghost_elems;
  std::vector<Tri3Subdivision *>::iterator       ghost_el     = ghost_elems.begin();
  const std::vector<Tri3Subdivision *>::iterator end_ghost_el = ghost_elems.end();
  for (; ghost_el != end_ghost_el; ++ghost_el)
    {
      Tri3Subdivision * elem = *ghost_el;
      libmesh_assert(elem->is_ghost());

      for (unsigned int i = 0; i < elem->n_sides(); ++i)
        {
          if (elem->neighbor_ptr(i) == libmesh_nullptr &&
              elem->neighbor_ptr(prev[i]) != libmesh_nullptr)
            {
              // go around counter-clockwise
              Tri3Subdivision * nb1 = static_cast<Tri3Subdivision *>(elem->neighbor_ptr(prev[i]));
              Tri3Subdivision * nb2 = nb1;
              unsigned int j = i;
              unsigned int n_nb = 0;
              while (nb1 != libmesh_nullptr && nb1->id() != elem->id())
                {
                  j = nb1->local_node_number(elem->node_id(i));
                  nb2 = nb1;
                  nb1 = static_cast<Tri3Subdivision *>(nb1->neighbor_ptr(prev[j]));
                  libmesh_assert(nb1 == libmesh_nullptr || nb1->id() != nb2->id());
                  n_nb++;
                }

              libmesh_assert_not_equal_to(nb2->id(), elem->id());

              // Above, we merged coinciding ghost vertices. Therefore, we need
              // to exclude the case where there is no ghost element to add between
              // these two (identical) ghost nodes.
              if (elem->node_ptr(next[i])->id() == nb2->node_ptr(prev[j])->id())
                break;

              // If the number of already present neighbors is less than 4, we add another extra element
              // so that the node in the middle of the loop ends up being of valence=6.
              // This case usually happens when the middle node corresponds to a corner of the original mesh,
              // and the extra element below prevents degenerated triangles in the mesh corners.
              if (n_nb < 4)
                {
                  // this is the vertex to be mirrored
                  Point point = nb2->point(j) + nb2->point(prev[j]) - nb2->point(next[j]);

                  // Check if the proposed vertex doesn't coincide with one of the existing vertices.
                  // This is necessary because for some triangulations, it can happen that two mirrored
                  // ghost vertices coincide, which would then lead to a zero size ghost element below.
                  Node * node = libmesh_nullptr;
                  for (std::size_t k = 0; k < ghost_nodes.size(); ++k)
                    {
                      if ((*ghost_nodes[k] - point).norm() < tol * (nb2->point(j) - point).norm())
                        {
                          node = ghost_nodes[k];
                          break;
                        }
                    }

                  // add the new vertex only if no other is nearby
                  if (node == libmesh_nullptr)
                    {
                      node = mesh.add_point(point);
                      ghost_nodes.push_back(node);
                    }

                  Tri3Subdivision * newelem = new Tri3Subdivision();

                  newelem->set_node(0) = nb2->node_ptr(j);
                  newelem->set_node(1) = nb2->node_ptr(prev[j]);
                  newelem->set_node(2) = node;
                  newelem->set_neighbor(0, nb2);
                  newelem->set_neighbor(1, libmesh_nullptr);
                  newelem->set_ghost(true);
                  nb2->set_neighbor(prev[j], newelem);

                  mesh.add_elem(newelem);
                  mesh.get_boundary_info().add_node(nb2->node_ptr(j), 1);
                  mesh.get_boundary_info().add_node(nb2->node_ptr(prev[j]), 1);
                  mesh.get_boundary_info().add_node(node, 1);

                  nb2 = newelem;
                  j = nb2->local_node_number(elem->node_id(i));
                }

              Tri3Subdivision * newelem = new Tri3Subdivision();
              newelem->set_node(0) = elem->node_ptr(next[i]);
              newelem->set_node(1) = elem->node_ptr(i);
              newelem->set_node(2) = nb2->node_ptr(prev[j]);
              newelem->set_neighbor(0, elem);
              newelem->set_neighbor(1, nb2);
              newelem->set_neighbor(2, libmesh_nullptr);
              newelem->set_ghost(true);

              elem->set_neighbor(i, newelem);
              nb2->set_neighbor(prev[j], newelem);

              missing_ghost_elems.push_back(newelem);
              break;
            }
        } // end side loop
    } // end ghost element loop

  // add the missing ghost elements to the mesh
  std::vector<Tri3Subdivision *>::iterator       missing_el     = missing_ghost_elems.begin();
  const std::vector<Tri3Subdivision *>::iterator end_missing_el = missing_ghost_elems.end();
  for (; missing_el != end_missing_el; ++missing_el)
    mesh.add_elem(*missing_el);
}

} // namespace libMesh
