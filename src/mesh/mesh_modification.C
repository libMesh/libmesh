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



// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::acos()
#include <algorithm>
#include <limits>
#include <map>

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/function_base.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_tet10.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel.h"
#include "libmesh/remote_elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/unstructured_mesh.h"


namespace libMesh
{


// ------------------------------------------------------------
// MeshTools::Modification functions for mesh modification
void MeshTools::Modification::distort (MeshBase & mesh,
                                       const Real factor,
                                       const bool perturb_boundary)
{
  libmesh_assert (mesh.n_nodes());
  libmesh_assert (mesh.n_elem());
  libmesh_assert ((factor >= 0.) && (factor <= 1.));

  LOG_SCOPE("distort()", "MeshTools::Modification");

  // If we are not perturbing boundary nodes, make a
  // quickly-searchable list of node ids we can check against.
  std::unordered_set<dof_id_type> boundary_node_ids;
  if (!perturb_boundary)
    boundary_node_ids = MeshTools::find_boundary_nodes (mesh);

  // Now calculate the minimum distance to
  // neighboring nodes for each node.
  // hmin holds these distances.
  std::vector<float> hmin (mesh.max_node_id(),
                           std::numeric_limits<float>::max());

  for (const auto & elem : mesh.active_element_ptr_range())
    for (auto & n : elem->node_ref_range())
      hmin[n.id()] = std::min(hmin[n.id()],
                              static_cast<float>(elem->hmin()));

  // Now actually move the nodes
  {
    const unsigned int seed = 123456;

    // seed the random number generator.
    // We'll loop from 1 to n_nodes on every processor, even those
    // that don't have a particular node, so that the pseudorandom
    // numbers will be the same everywhere.
    std::srand(seed);

    // If the node is on the boundary or
    // the node is not used by any element (hmin[n]<1.e20)
    // then we should not move it.
    // [Note: Testing for (in)equality might be wrong
    // (different types, namely float and double)]
    for (auto n : IntRange<unsigned int>(0, mesh.max_node_id()))
      if ((perturb_boundary || !boundary_node_ids.count(n)) && hmin[n] < 1.e20)
        {
          // the direction, random but unit normalized
          Point dir (static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX),
                     (mesh.mesh_dimension() > 1) ? static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX) : 0.,
                     ((mesh.mesh_dimension() == 3) ? static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX) : 0.));

          dir(0) = (dir(0)-.5)*2.;
#if LIBMESH_DIM > 1
          if (mesh.mesh_dimension() > 1)
            dir(1) = (dir(1)-.5)*2.;
#endif
#if LIBMESH_DIM > 2
          if (mesh.mesh_dimension() == 3)
            dir(2) = (dir(2)-.5)*2.;
#endif

          dir = dir.unit();

          Node * node = mesh.query_node_ptr(n);
          if (!node)
            continue;

          (*node)(0) += dir(0)*factor*hmin[n];
#if LIBMESH_DIM > 1
          if (mesh.mesh_dimension() > 1)
            (*node)(1) += dir(1)*factor*hmin[n];
#endif
#if LIBMESH_DIM > 2
          if (mesh.mesh_dimension() == 3)
            (*node)(2) += dir(2)*factor*hmin[n];
#endif
        }
  }
}



void MeshTools::Modification::translate (MeshBase & mesh,
                                         const Real xt,
                                         const Real yt,
                                         const Real zt)
{
  const Point p(xt, yt, zt);

  for (auto & node : mesh.node_ptr_range())
    *node += p;
}


// void MeshTools::Modification::rotate2D (MeshBase & mesh,
//                                         const Real alpha)
// {
//   libmesh_assert_not_equal_to (mesh.mesh_dimension(), 1);

//   const Real pi = std::acos(-1);
//   const Real  a = alpha/180.*pi;
//   for (unsigned int n=0; n<mesh.n_nodes(); n++)
//     {
//       const Point p = mesh.node_ref(n);
//       const Real  x = p(0);
//       const Real  y = p(1);
//       const Real  z = p(2);
//       mesh.node_ref(n) = Point(std::cos(a)*x - std::sin(a)*y,
//                                std::sin(a)*x + std::cos(a)*y,
//                                z);
//     }

// }



void MeshTools::Modification::rotate (MeshBase & mesh,
                                      const Real phi,
                                      const Real theta,
                                      const Real psi)
{
#if LIBMESH_DIM == 3
  const Real  p = -phi/180.*libMesh::pi;
  const Real  t = -theta/180.*libMesh::pi;
  const Real  s = -psi/180.*libMesh::pi;
  const Real sp = std::sin(p), cp = std::cos(p);
  const Real st = std::sin(t), ct = std::cos(t);
  const Real ss = std::sin(s), cs = std::cos(s);

  // We follow the convention described at http://mathworld.wolfram.com/EulerAngles.html
  // (equations 6-14 give the entries of the composite transformation matrix).
  // The rotations are performed sequentially about the z, x, and z axes, in that order.
  // A positive angle yields a counter-clockwise rotation about the axis in question.
  for (auto & node : mesh.node_ptr_range())
    {
      const Point pt = *node;
      const Real  x  = pt(0);
      const Real  y  = pt(1);
      const Real  z  = pt(2);
      *node = Point(( cp*cs-sp*ct*ss)*x + ( sp*cs+cp*ct*ss)*y + (st*ss)*z,
                    (-cp*ss-sp*ct*cs)*x + (-sp*ss+cp*ct*cs)*y + (st*cs)*z,
                    ( sp*st)*x          + (-cp*st)*y          + (ct)*z   );
    }
#else
  libmesh_ignore(mesh, phi, theta, psi);
  libmesh_error_msg("MeshTools::Modification::rotate() requires libMesh to be compiled with LIBMESH_DIM==3");
#endif
}


void MeshTools::Modification::scale (MeshBase & mesh,
                                     const Real xs,
                                     const Real ys,
                                     const Real zs)
{
  const Real x_scale = xs;
  Real y_scale       = ys;
  Real z_scale       = zs;

  if (ys == 0.)
    {
      libmesh_assert_equal_to (zs, 0.);

      y_scale = z_scale = x_scale;
    }

  // Scale the x coordinate in all dimensions
  for (auto & node : mesh.node_ptr_range())
    (*node)(0) *= x_scale;

  // Only scale the y coordinate in 2 and 3D
  if (LIBMESH_DIM < 2)
    return;

  for (auto & node : mesh.node_ptr_range())
    (*node)(1) *= y_scale;

  // Only scale the z coordinate in 3D
  if (LIBMESH_DIM < 3)
    return;

  for (auto & node : mesh.node_ptr_range())
    (*node)(2) *= z_scale;
}


void MeshTools::Modification::smooth (MeshBase & mesh,
                                      const unsigned int n_iterations,
                                      const Real power)
{
  /**
   * This implementation assumes every element "side" has only 2 nodes.
   */
  libmesh_assert_equal_to (mesh.mesh_dimension(), 2);

  /*
   * Create a quickly-searchable list of boundary nodes.
   */
  std::unordered_set<dof_id_type> boundary_node_ids =
    MeshTools::find_boundary_nodes (mesh);

  for (unsigned int iter=0; iter<n_iterations; iter++)
    {
      /*
       * loop over the mesh refinement level
       */
      unsigned int n_levels = MeshTools::n_levels(mesh);
      for (unsigned int refinement_level=0; refinement_level != n_levels;
           refinement_level++)
        {
          // initialize the storage (have to do it on every level to get empty vectors
          std::vector<Point> new_positions;
          std::vector<Real>   weight;
          new_positions.resize(mesh.n_nodes());
          weight.resize(mesh.n_nodes());

          {
            // Loop over the elements to calculate new node positions
            for (const auto & elem : as_range(mesh.level_elements_begin(refinement_level),
                                              mesh.level_elements_end(refinement_level)))
              {
                /*
                 * We relax all nodes on level 0 first
                 * If the element is refined (level > 0), we interpolate the
                 * parents nodes with help of the embedding matrix
                 */
                if (refinement_level == 0)
                  {
                    for (auto s : elem->side_index_range())
                      {
                        /*
                         * Only operate on sides which are on the
                         * boundary or for which the current element's
                         * id is greater than its neighbor's.
                         * Sides get only built once.
                         */
                        if ((elem->neighbor_ptr(s) != nullptr) &&
                            (elem->id() > elem->neighbor_ptr(s)->id()))
                          {
                            std::unique_ptr<const Elem> side(elem->build_side_ptr(s));

                            const Node & node0 = side->node_ref(0);
                            const Node & node1 = side->node_ref(1);

                            Real node_weight = 1.;
                            // calculate the weight of the nodes
                            if (power > 0)
                              {
                                Point diff = node0-node1;
                                node_weight = std::pow(diff.norm(), power);
                              }

                            const dof_id_type id0 = node0.id(), id1 = node1.id();
                            new_positions[id0].add_scaled( node1, node_weight );
                            new_positions[id1].add_scaled( node0, node_weight );
                            weight[id0] += node_weight;
                            weight[id1] += node_weight;
                          }
                      } // element neighbor loop
                  }
#ifdef LIBMESH_ENABLE_AMR
                else   // refinement_level > 0
                  {
                    /*
                     * Find the positions of the hanging nodes of refined elements.
                     * We do this by calculating their position based on the parent
                     * (one level less refined) element, and the embedding matrix
                     */

                    const Elem * parent = elem->parent();

                    /*
                     * find out which child I am
                     */
                    unsigned int c = parent->which_child_am_i(elem);
                    /*
                     *loop over the childs (that is, the current elements) nodes
                     */
                    for (auto nc : elem->node_index_range())
                      {
                        /*
                         * the new position of the node
                         */
                        Point point;
                        for (auto n : parent->node_index_range())
                          {
                            /*
                             * The value from the embedding matrix
                             */
                            const float em_val = parent->embedding_matrix(c,nc,n);

                            if (em_val != 0.)
                              point.add_scaled (parent->point(n), em_val);
                          }

                        const dof_id_type id = elem->node_ptr(nc)->id();
                        new_positions[id] = point;
                        weight[id] = 1.;
                      }
                  } // if element refinement_level
#endif // #ifdef LIBMESH_ENABLE_AMR

              } // element loop

            /*
             * finally reposition the vertex nodes
             */
            for (auto nid : IntRange<unsigned int>(0, mesh.n_nodes()))
              if (!boundary_node_ids.count(nid) && weight[nid] > 0.)
                mesh.node_ref(nid) = new_positions[nid]/weight[nid];
          }

          // Now handle the additional second_order nodes by calculating
          // their position based on the vertex positions
          // we do a second loop over the level elements
          for (auto & elem : as_range(mesh.level_elements_begin(refinement_level),
                                      mesh.level_elements_end(refinement_level)))
            {
              const unsigned int son_begin = elem->n_vertices();
              const unsigned int son_end   = elem->n_nodes();
              for (unsigned int n=son_begin; n<son_end; n++)
                {
                  const unsigned int n_adjacent_vertices =
                    elem->n_second_order_adjacent_vertices(n);

                  Point point;
                  for (unsigned int v=0; v<n_adjacent_vertices; v++)
                    point.add(elem->point( elem->second_order_adjacent_vertex(n,v) ));

                  const dof_id_type id = elem->node_ptr(n)->id();
                  mesh.node_ref(id) = point/n_adjacent_vertices;
                }
            }
        } // refinement_level loop
    } // end iteration
}


void MeshTools::Modification::change_boundary_id (MeshBase & mesh,
                                                  const boundary_id_type old_id,
                                                  const boundary_id_type new_id)
{
  if (old_id == new_id)
    {
      // If the IDs are the same, this is a no-op.
      return;
    }

  // A reference to the Mesh's BoundaryInfo object, for convenience.
  BoundaryInfo & bi = mesh.get_boundary_info();

  {
    // Temporary vector to hold ids
    std::vector<boundary_id_type> bndry_ids;

    // build_node_list returns a vector of (node, bc) tuples.
    for (const auto & t : bi.build_node_list())
      if (std::get<1>(t) == old_id)
        {
          // Get the node in question
          const Node * node = mesh.node_ptr(std::get<0>(t));

          // Get all the current IDs for this node.
          bi.boundary_ids(node, bndry_ids);

          // Update the IDs accordingly
          std::replace(bndry_ids.begin(), bndry_ids.end(), old_id, new_id);

          // Remove all traces of that node from the BoundaryInfo object
          bi.remove(node);

          // Add it back with the updated IDs
          bi.add_node(node, bndry_ids);
        }
  }

  {
    // Temporary vector to hold ids
    std::vector<boundary_id_type> bndry_ids;

    // build_edge_list returns a vector of (elem, side, bc) tuples.
    for (const auto & t : bi.build_edge_list())
      if (std::get<2>(t) == old_id)
        {
          // Get the elem in question
          const Elem * elem = mesh.elem_ptr(std::get<0>(t));

          // The edge of the elem in question
          unsigned short int edge = std::get<1>(t);

          // Get all the current IDs for the edge in question.
          bi.edge_boundary_ids(elem, edge, bndry_ids);

          // Update the IDs accordingly
          std::replace(bndry_ids.begin(), bndry_ids.end(), old_id, new_id);

          // Remove all traces of that edge from the BoundaryInfo object
          bi.remove_edge(elem, edge);

          // Add it back with the updated IDs
          bi.add_edge(elem, edge, bndry_ids);
        }
  }

  {
    // Temporary vector to hold ids
    std::vector<boundary_id_type> bndry_ids;

    // build_shellface_list returns a vector of (elem, side, bc) tuples.
    for (const auto & t : bi.build_shellface_list())
      if (std::get<2>(t) == old_id)
        {
          // Get the elem in question
          const Elem * elem = mesh.elem_ptr(std::get<0>(t));

          // The shellface of the elem in question
          unsigned short int shellface = std::get<1>(t);

          // Get all the current IDs for the shellface in question.
          bi.shellface_boundary_ids(elem, shellface, bndry_ids);

          // Update the IDs accordingly
          std::replace(bndry_ids.begin(), bndry_ids.end(), old_id, new_id);

          // Remove all traces of that shellface from the BoundaryInfo object
          bi.remove_shellface(elem, shellface);

          // Add it back with the updated IDs
          bi.add_shellface(elem, shellface, bndry_ids);
        }
  }

  {
    // Temporary vector to hold ids
    std::vector<boundary_id_type> bndry_ids;

    // build_side_list returns a vector of (elem, side, bc) tuples.
    for (const auto & t : bi.build_side_list())
      if (std::get<2>(t) == old_id)
        {
          // Get the elem in question
          const Elem * elem = mesh.elem_ptr(std::get<0>(t));

          // The side of the elem in question
          unsigned short int side = std::get<1>(t);

          // Get all the current IDs for the side in question.
          bi.boundary_ids(elem, side, bndry_ids);

          // Update the IDs accordingly
          std::replace(bndry_ids.begin(), bndry_ids.end(), old_id, new_id);

          // Remove all traces of that side from the BoundaryInfo object
          bi.remove_side(elem, side);

          // Add it back with the updated IDs
          bi.add_side(elem, side, bndry_ids);
        }
  }

  // Remove any remaining references to the old_id so it does not show
  // up in lists of boundary ids, etc.
  bi.remove_id(old_id);
}



void MeshTools::Modification::change_subdomain_id (MeshBase & mesh,
                                                   const subdomain_id_type old_id,
                                                   const subdomain_id_type new_id)
{
  if (old_id == new_id)
    {
      // If the IDs are the same, this is a no-op.
      return;
    }

  for (auto & elem : mesh.element_ptr_range())
    {
      if (elem->subdomain_id() == old_id)
        elem->subdomain_id() = new_id;
    }
}


} // namespace libMesh
