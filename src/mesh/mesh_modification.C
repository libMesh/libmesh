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



// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::acos()
#include <algorithm>
#include <limits>
#include <map>
#include <array>

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
#include "libmesh/elem_side_builder.h"
#include "libmesh/tensor_value.h"

namespace
{
using namespace libMesh;

bool split_first_diagonal(const Elem * elem,
                          unsigned int diag_1_node_1,
                          unsigned int diag_1_node_2,
                          unsigned int diag_2_node_1,
                          unsigned int diag_2_node_2)
{
  return ((elem->node_id(diag_1_node_1) > elem->node_id(diag_2_node_1) &&
           elem->node_id(diag_1_node_1) > elem->node_id(diag_2_node_2)) ||
          (elem->node_id(diag_1_node_2) > elem->node_id(diag_2_node_1) &&
           elem->node_id(diag_1_node_2) > elem->node_id(diag_2_node_2)));
}


// Return the local index of the vertex on \p elem with the highest
// node id.
unsigned int highest_vertex_on(const Elem * elem)
{
  unsigned int highest_n = 0;
  dof_id_type highest_n_id = elem->node_id(0);
  for (auto n : make_range(1u, elem->n_vertices()))
    {
      const dof_id_type n_id = elem->node_id(n);
      if (n_id > highest_n_id)
        {
          highest_n = n;
          highest_n_id = n_id;
        }
    }

  return highest_n;
}


static const std::array<std::array<unsigned int, 3>, 8> opposing_nodes =
{{ {2,5,7},{3,4,6},{0,5,7},{1,4,6},{1,3,6},{0,2,7},{1,3,4},{0,2,5} }};


// Find the highest id on these side nodes of this element
std::pair<unsigned int, unsigned int>
split_diagonal(const Elem * elem,
               const std::vector<unsigned int> & nodes_on_side)
{
  libmesh_assert_equal_to(elem->type(), HEX8);

  unsigned int highest_n = nodes_on_side.front();
  dof_id_type highest_n_id = elem->node_id(nodes_on_side.front());
  for (auto n : nodes_on_side)
    {
      const dof_id_type n_id = elem->node_id(n);
      if (n_id > highest_n_id)
        {
          highest_n = n;
          highest_n_id = n_id;
        }
    }

  for (auto n : nodes_on_side)
    {
      for (auto n2 : opposing_nodes[highest_n])
        if (n2 == n)
          return std::make_pair(highest_n, n2);
    }

  libmesh_error();

  return std::make_pair(libMesh::invalid_uint, libMesh::invalid_uint);
}


// Reconstruct a C++20 feature in C++14
template <typename T>
struct reversion_wrapper { T& iterable; };

template <typename T>
auto begin (reversion_wrapper<T> w) {return std::rbegin(w.iterable);}

template <typename T>
auto end (reversion_wrapper<T> w) {return std::rend(w.iterable);}

template <typename T>
reversion_wrapper<T> reverse(T&& iterable) {return {iterable};}

}


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
    for (auto n : make_range(mesh.max_node_id()))
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

  // We haven't changed any topology, but just changing geometry could
  // have invalidated a point locator.
  mesh.clear_point_locator();
}



void MeshTools::Modification::permute_elements(MeshBase & mesh)
{
  LOG_SCOPE("permute_elements()", "MeshTools::Modification");

  // We don't yet support doing permute() on a parent element, which
  // would require us to consistently permute all its children and
  // give them different local child numbers.
  unsigned int n_levels = MeshTools::n_levels(mesh);
  if (n_levels > 1)
    libmesh_error();

  const unsigned int seed = 123456;

  // seed the random number generator.
  // We'll loop from 1 to max_elem_id on every processor, even those
  // that don't have a particular element, so that the pseudorandom
  // numbers will be the same everywhere.
  std::srand(seed);


  for (auto e_id : make_range(mesh.max_elem_id()))
    {
      int my_rand = std::rand();

      Elem * elem = mesh.query_elem_ptr(e_id);

      if (!elem)
        continue;

      const unsigned int max_permutation = elem->n_permutations();
      if (!max_permutation)
        continue;

      const unsigned int perm = my_rand % max_permutation;

      elem->permute(perm);
    }
}


void MeshTools::Modification::orient_elements(MeshBase & mesh)
{
  LOG_SCOPE("orient_elements()", "MeshTools::Modification");

  // We don't yet support doing orient() on a parent element, which
  // would require us to consistently orient all its children and
  // give them different local child numbers.
  unsigned int n_levels = MeshTools::n_levels(mesh);
  if (n_levels > 1)
    libmesh_not_implemented_msg("orient_elements() does not support refined meshes");

  BoundaryInfo & boundary_info = mesh.get_boundary_info();
  for (auto elem : mesh.element_ptr_range())
    elem->orient(&boundary_info);
}



void MeshTools::Modification::redistribute (MeshBase & mesh,
                                            const FunctionBase<Real> & mapfunc)
{
  libmesh_assert (mesh.n_nodes());
  libmesh_assert (mesh.n_elem());

  LOG_SCOPE("redistribute()", "MeshTools::Modification");

  DenseVector<Real> output_vec(LIBMESH_DIM);

  // FIXME - we should thread this later.
  std::unique_ptr<FunctionBase<Real>> myfunc = mapfunc.clone();

  for (auto & node : mesh.node_ptr_range())
    {
      (*myfunc)(*node, output_vec);

      (*node)(0) = output_vec(0);
#if LIBMESH_DIM > 1
      (*node)(1) = output_vec(1);
#endif
#if LIBMESH_DIM > 2
      (*node)(2) = output_vec(2);
#endif
    }

  // We haven't changed any topology, but just changing geometry could
  // have invalidated a point locator.
  mesh.clear_point_locator();
}



void MeshTools::Modification::translate (MeshBase & mesh,
                                         const Real xt,
                                         const Real yt,
                                         const Real zt)
{
  const Point p(xt, yt, zt);

  for (auto & node : mesh.node_ptr_range())
    *node += p;

  // We haven't changed any topology, but just changing geometry could
  // have invalidated a point locator.
  mesh.clear_point_locator();
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



RealTensorValue
MeshTools::Modification::rotate (MeshBase & mesh,
                                 const Real phi,
                                 const Real theta,
                                 const Real psi)
{
  // We won't change any topology, but just changing geometry could
  // invalidate a point locator.
  mesh.clear_point_locator();

#if LIBMESH_DIM == 3
  const auto R = RealTensorValue::intrinsic_rotation_matrix(phi, theta, psi);

  if (theta)
    mesh.set_spatial_dimension(3);

  for (auto & node : mesh.node_ptr_range())
    {
      Point & pt = *node;
      pt = R * pt;
    }

  return R;

#else
  libmesh_ignore(mesh, phi, theta, psi);
  libmesh_error_msg("MeshTools::Modification::rotate() requires libMesh to be compiled with LIBMESH_DIM==3");
  // We'll never get here
  return RealTensorValue();
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

  // We haven't changed any topology, but just changing geometry could
  // have invalidated a point locator.
  mesh.clear_point_locator();
}



void MeshTools::Modification::all_tri (MeshBase & mesh)
{
  if (!mesh.is_replicated() && !mesh.is_prepared())
    mesh.prepare_for_use();

  // The number of elements in the original mesh before any additions
  // or deletions.
  const dof_id_type n_orig_elem = mesh.n_elem();
  const dof_id_type max_orig_id = mesh.max_elem_id();

  // We store pointers to the newly created elements in a vector
  // until they are ready to be added to the mesh.  This is because
  // adding new elements on the fly can cause reallocation and invalidation
  // of existing mesh element_iterators.
  std::vector<std::unique_ptr<Elem>> new_elements;

  unsigned int max_subelems = 1;  // in 1D nothing needs to change
  if (mesh.mesh_dimension() == 2) // in 2D quads can split into 2 tris
    max_subelems = 2;
  if (mesh.mesh_dimension() == 3) // in 3D hexes can split into 6 tets
    max_subelems = 6;

  new_elements.reserve (max_subelems*n_orig_elem);

  // If the original mesh has *side* boundary data, we carry that over
  // to the new mesh with triangular elements.  We currently only
  // support bringing over side-based BCs to the all-tri mesh, but
  // that could probably be extended to node and edge-based BCs as
  // well.
  const bool mesh_has_boundary_data = (mesh.get_boundary_info().n_boundary_conds() > 0);

  // Temporary vectors to store the new boundary element pointers, side numbers, and boundary ids
  std::vector<Elem *> new_bndry_elements;
  std::vector<unsigned short int> new_bndry_sides;
  std::vector<boundary_id_type> new_bndry_ids;

  // We may need to add new points if we run into a 1.5th order
  // element; if we do that on a DistributedMesh in a ghost element then
  // we will need to fix their ids / unique_ids
  bool added_new_ghost_point = false;

  // Iterate over the elements, splitting:
  // QUADs into pairs of conforming triangles
  // PYRAMIDs into pairs of conforming tets,
  // PRISMs into triplets of conforming tets, and
  // HEXs into quintets or sextets of conforming tets.
  // We split on the shortest diagonal to give us better
  // triangle quality in 2D, and we split based on node ids
  // to guarantee consistency in 3D.

  // FIXME: This algorithm does not work on refined grids!
  {
#ifdef LIBMESH_ENABLE_UNIQUE_ID
    unique_id_type max_unique_id = mesh.parallel_max_unique_id();
#endif

    // For avoiding extraneous allocation when building side elements
    std::unique_ptr<const Elem> elem_side, subside_elem;

    for (auto & elem : mesh.element_ptr_range())
      {
        const ElemType etype = elem->type();

        // all_tri currently only works on coarse meshes
        if (elem->parent())
          libmesh_not_implemented_msg("Cannot convert a refined element into simplices\n");

        // The new elements we will split the quad into.
        // In 3D we may need as many as 6 tets per hex
        std::array<std::unique_ptr<Elem>, 6> subelem {};

        switch (etype)
          {
          case QUAD4:
            {
              subelem[0] = Elem::build(TRI3);
              subelem[1] = Elem::build(TRI3);

              // Check for possible edge swap
              if ((elem->point(0) - elem->point(2)).norm() <
                  (elem->point(1) - elem->point(3)).norm())
                {
                  subelem[0]->set_node(0, elem->node_ptr(0));
                  subelem[0]->set_node(1, elem->node_ptr(1));
                  subelem[0]->set_node(2, elem->node_ptr(2));

                  subelem[1]->set_node(0, elem->node_ptr(0));
                  subelem[1]->set_node(1, elem->node_ptr(2));
                  subelem[1]->set_node(2, elem->node_ptr(3));
                }

              else
                {
                  subelem[0]->set_node(0, elem->node_ptr(0));
                  subelem[0]->set_node(1, elem->node_ptr(1));
                  subelem[0]->set_node(2, elem->node_ptr(3));

                  subelem[1]->set_node(0, elem->node_ptr(1));
                  subelem[1]->set_node(1, elem->node_ptr(2));
                  subelem[1]->set_node(2, elem->node_ptr(3));
                }


              break;
            }

          case QUAD8:
            {
              if (elem->processor_id() != mesh.processor_id())
                added_new_ghost_point = true;

              subelem[0] = Elem::build(TRI6);
              subelem[1] = Elem::build(TRI6);

              // Add a new node at the center (vertex average) of the element.
              Node * new_node = mesh.add_point((mesh.point(elem->node_id(0)) +
                                                mesh.point(elem->node_id(1)) +
                                                mesh.point(elem->node_id(2)) +
                                                mesh.point(elem->node_id(3)))/4,
                                               DofObject::invalid_id,
                                               elem->processor_id());

              // Check for possible edge swap
              if ((elem->point(0) - elem->point(2)).norm() <
                  (elem->point(1) - elem->point(3)).norm())
                {
                  subelem[0]->set_node(0, elem->node_ptr(0));
                  subelem[0]->set_node(1, elem->node_ptr(1));
                  subelem[0]->set_node(2, elem->node_ptr(2));
                  subelem[0]->set_node(3, elem->node_ptr(4));
                  subelem[0]->set_node(4, elem->node_ptr(5));
                  subelem[0]->set_node(5, new_node);

                  subelem[1]->set_node(0, elem->node_ptr(0));
                  subelem[1]->set_node(1, elem->node_ptr(2));
                  subelem[1]->set_node(2, elem->node_ptr(3));
                  subelem[1]->set_node(3, new_node);
                  subelem[1]->set_node(4, elem->node_ptr(6));
                  subelem[1]->set_node(5, elem->node_ptr(7));

                }

              else
                {
                  subelem[0]->set_node(0, elem->node_ptr(3));
                  subelem[0]->set_node(1, elem->node_ptr(0));
                  subelem[0]->set_node(2, elem->node_ptr(1));
                  subelem[0]->set_node(3, elem->node_ptr(7));
                  subelem[0]->set_node(4, elem->node_ptr(4));
                  subelem[0]->set_node(5, new_node);

                  subelem[1]->set_node(0, elem->node_ptr(1));
                  subelem[1]->set_node(1, elem->node_ptr(2));
                  subelem[1]->set_node(2, elem->node_ptr(3));
                  subelem[1]->set_node(3, elem->node_ptr(5));
                  subelem[1]->set_node(4, elem->node_ptr(6));
                  subelem[1]->set_node(5, new_node);
                }

              break;
            }

          case QUAD9:
            {
              subelem[0] = Elem::build(TRI6);
              subelem[1] = Elem::build(TRI6);

              // Check for possible edge swap
              if ((elem->point(0) - elem->point(2)).norm() <
                  (elem->point(1) - elem->point(3)).norm())
                {
                  subelem[0]->set_node(0, elem->node_ptr(0));
                  subelem[0]->set_node(1, elem->node_ptr(1));
                  subelem[0]->set_node(2, elem->node_ptr(2));
                  subelem[0]->set_node(3, elem->node_ptr(4));
                  subelem[0]->set_node(4, elem->node_ptr(5));
                  subelem[0]->set_node(5, elem->node_ptr(8));

                  subelem[1]->set_node(0, elem->node_ptr(0));
                  subelem[1]->set_node(1, elem->node_ptr(2));
                  subelem[1]->set_node(2, elem->node_ptr(3));
                  subelem[1]->set_node(3, elem->node_ptr(8));
                  subelem[1]->set_node(4, elem->node_ptr(6));
                  subelem[1]->set_node(5, elem->node_ptr(7));
                }

              else
                {
                  subelem[0]->set_node(0, elem->node_ptr(0));
                  subelem[0]->set_node(1, elem->node_ptr(1));
                  subelem[0]->set_node(2, elem->node_ptr(3));
                  subelem[0]->set_node(3, elem->node_ptr(4));
                  subelem[0]->set_node(4, elem->node_ptr(8));
                  subelem[0]->set_node(5, elem->node_ptr(7));

                  subelem[1]->set_node(0, elem->node_ptr(1));
                  subelem[1]->set_node(1, elem->node_ptr(2));
                  subelem[1]->set_node(2, elem->node_ptr(3));
                  subelem[1]->set_node(3, elem->node_ptr(5));
                  subelem[1]->set_node(4, elem->node_ptr(6));
                  subelem[1]->set_node(5, elem->node_ptr(8));
                }

              break;
            }

          case HEX8:
            {
              BoundaryInfo & boundary_info = mesh.get_boundary_info();

              // Hexes all split into six tetrahedra
              subelem[0] = Elem::build(TET4);
              subelem[1] = Elem::build(TET4);
              subelem[2] = Elem::build(TET4);
              subelem[3] = Elem::build(TET4);
              subelem[4] = Elem::build(TET4);
              subelem[5] = Elem::build(TET4);

              // On faces, we choose the node with the highest
              // global id, and we split on the diagonal which
              // includes that node.  This ensures that (even in
              // parallel, even on distributed meshes) the same
              // diagonal split will be chosen for elements on either
              // side of the same quad face.
              const unsigned int highest_n = highest_vertex_on(elem);

              // opposing_node[n] is the local node number of the node
              // on the farthest corner of a hex8 from local node n
              static const std::array<unsigned int, 8> opposing_node =
                {6, 7, 4, 5, 2, 3, 0, 1};

              static const std::vector<std::vector<unsigned int>> sides_opposing_highest =
                {{2,3,5},{3,4,5},{1,4,5},{1,2,5},{0,2,3},{0,3,4},{0,1,4},{0,1,2}};
              static const std::vector<std::vector<unsigned int>> nodes_neighboring_highest =
                {{1,3,4},{0,2,5},{1,3,6},{0,2,7},{0,5,7},{1,4,6},{2,5,7},{3,4,6}};

              // Start by looking in three directions away from the
              // highest-id node.  In each direction there will be two
              // different possibilities for the split depending on
              // how the opposing face nodes are numbered.
              //
              // This is tricky enough that I'm not going to worry
              // about manually keeping tets oriented; we'll just call
              // orient() on each as we go.

              unsigned int next_subelem = 0;
              for (auto side : sides_opposing_highest[highest_n])
                {
                  const std::vector<unsigned int> nodes_on_side =
                    elem->nodes_on_side(side);

                  auto [dn, dn2] = split_diagonal(elem, nodes_on_side);

                  unsigned int split_on_neighbor = false;
                  for (auto n : nodes_neighboring_highest[highest_n])
                    if (dn == n || dn2 == n)
                      {
                        split_on_neighbor = true;
                        break;
                      }

                  // Add one or two elements for each opposing side,
                  // depending on whether the diagonal split there
                  // connects to the neighboring diagonal split or
                  // not.
                  if (split_on_neighbor)
                    {
                      subelem[next_subelem]->set_node(0, elem->node_ptr(highest_n));
                      subelem[next_subelem]->set_node(1, elem->node_ptr(dn));
                      subelem[next_subelem]->set_node(2, elem->node_ptr(dn2));
                      for (auto n : nodes_on_side)
                        if (n != dn && n != dn2)
                          {
                            subelem[next_subelem]->set_node(3, elem->node_ptr(n));
                            break;
                          }
                      subelem[next_subelem]->orient(&boundary_info);
                      ++next_subelem;

                      subelem[next_subelem]->set_node(0, elem->node_ptr(highest_n));
                      subelem[next_subelem]->set_node(1, elem->node_ptr(dn));
                      subelem[next_subelem]->set_node(2, elem->node_ptr(dn2));
                      for (auto n : reverse(nodes_on_side))
                        if (n != dn && n != dn2)
                          {
                            subelem[next_subelem]->set_node(3, elem->node_ptr(n));
                            break;
                          }
                      subelem[next_subelem]->orient(&boundary_info);
                      ++next_subelem;
                    }
                  else
                    {
                      subelem[next_subelem]->set_node(0, elem->node_ptr(highest_n));
                      subelem[next_subelem]->set_node(1, elem->node_ptr(dn));
                      subelem[next_subelem]->set_node(2, elem->node_ptr(dn2));
                      for (auto n : nodes_on_side)
                        for (auto n2 : nodes_neighboring_highest[highest_n])
                          if (n == n2)
                            {
                              subelem[next_subelem]->set_node(3, elem->node_ptr(n));
                              goto break_both_loops;
                            }

                      break_both_loops:
                      subelem[next_subelem]->orient(&boundary_info);
                      ++next_subelem;
                    }
                }

              // At this point we've created between 3 and 6 tets.
              // What's left to do depends on how many.

              // If we just chopped off three vertices into three
              // tets, then the best way to split this hex would be
              // the symmetric five-split.  Chop off the opposing
              // vertex too, and then the remaining interior is our
              // final tet.
              if (next_subelem == 3)
                {
                  subelem[next_subelem]->set_node(0, elem->node_ptr(opposing_nodes[highest_n][0]));
                  subelem[next_subelem]->set_node(1, elem->node_ptr(opposing_nodes[highest_n][1]));
                  subelem[next_subelem]->set_node(2, elem->node_ptr(opposing_nodes[highest_n][2]));
                  subelem[next_subelem]->set_node(3, elem->node_ptr(opposing_node[highest_n]));
                  subelem[next_subelem]->orient(&boundary_info);
                  ++next_subelem;

                  subelem[next_subelem]->set_node(0, elem->node_ptr(opposing_nodes[highest_n][0]));
                  subelem[next_subelem]->set_node(1, elem->node_ptr(opposing_nodes[highest_n][1]));
                  subelem[next_subelem]->set_node(2, elem->node_ptr(opposing_nodes[highest_n][2]));
                  subelem[next_subelem]->set_node(3, elem->node_ptr(highest_n));
                  subelem[next_subelem]->orient(&boundary_info);
                  ++next_subelem;

                  // We don't need the 6th tet after all
                  subelem[next_subelem].reset();
                  ++next_subelem;
                }

              // If we just chopped off one (or two) vertices into
              // tets, then the remaining gap is best (or only) filled
              // by pairing another tet with each.
              if (next_subelem == 4 ||
                  next_subelem == 5)
                {
                  for (auto side : sides_opposing_highest[highest_n])
                    {
                      const std::vector<unsigned int> nodes_on_side =
                        elem->nodes_on_side(side);

                      auto [dn, dn2] = split_diagonal(elem, nodes_on_side);

                      unsigned int split_on_neighbor = false;
                      for (auto n : nodes_neighboring_highest[highest_n])
                        if (dn == n || dn2 == n)
                          {
                            split_on_neighbor = true;
                            break;
                          }

                      // The two !split_on_neighbor sides are where we
                      // need the two remaining tets
                      if (!split_on_neighbor)
                        {
                          subelem[next_subelem]->set_node(0, elem->node_ptr(highest_n));
                          subelem[next_subelem]->set_node(1, elem->node_ptr(dn));
                          subelem[next_subelem]->set_node(2, elem->node_ptr(dn2));
                          subelem[next_subelem]->set_node(3, elem->node_ptr(opposing_node[highest_n]));
                          subelem[next_subelem]->orient(&boundary_info);
                          ++next_subelem;
                        }
                    }
                }

              // Whether we got there by creating six tets from the
              // first for loop or by patching up the split afterward,
              // we should have considered six tets (possibly
              // including one deleted one...) at this point.
              libmesh_assert(next_subelem == 6);

              break;
            }

          case PRISM6:
            {
              // Prisms all split into three tetrahedra
              subelem[0] = Elem::build(TET4);
              subelem[1] = Elem::build(TET4);
              subelem[2] = Elem::build(TET4);

              // Triangular faces are not split.

              // On quad faces, we choose the node with the highest
              // global id, and we split on the diagonal which
              // includes that node.  This ensures that (even in
              // parallel, even on distributed meshes) the same
              // diagonal split will be chosen for elements on either
              // side of the same quad face.  It also ensures that we
              // always have a mix of "clockwise" and
              // "counterclockwise" split faces (two of one and one
              // of the other on each prism; this is useful since the
              // alternative all-clockwise or all-counterclockwise
              // face splittings can't be turned into tets without
              // adding more nodes

              // Split on 0-4 diagonal
              if (split_first_diagonal(elem, 0,4, 1,3))
                {
                  // Split on 0-5 diagonal
                  if (split_first_diagonal(elem, 0,5, 2,3))
                    {
                      // Split on 1-5 diagonal
                      if (split_first_diagonal(elem, 1,5, 2,4))
                        {
                          subelem[0]->set_node(0, elem->node_ptr(0));
                          subelem[0]->set_node(1, elem->node_ptr(4));
                          subelem[0]->set_node(2, elem->node_ptr(5));
                          subelem[0]->set_node(3, elem->node_ptr(3));

                          subelem[1]->set_node(0, elem->node_ptr(0));
                          subelem[1]->set_node(1, elem->node_ptr(4));
                          subelem[1]->set_node(2, elem->node_ptr(1));
                          subelem[1]->set_node(3, elem->node_ptr(5));

                          subelem[2]->set_node(0, elem->node_ptr(0));
                          subelem[2]->set_node(1, elem->node_ptr(1));
                          subelem[2]->set_node(2, elem->node_ptr(2));
                          subelem[2]->set_node(3, elem->node_ptr(5));
                        }
                      else // Split on 2-4 diagonal
                        {
                          libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                          subelem[0]->set_node(0, elem->node_ptr(0));
                          subelem[0]->set_node(1, elem->node_ptr(4));
                          subelem[0]->set_node(2, elem->node_ptr(5));
                          subelem[0]->set_node(3, elem->node_ptr(3));

                          subelem[1]->set_node(0, elem->node_ptr(0));
                          subelem[1]->set_node(1, elem->node_ptr(4));
                          subelem[1]->set_node(2, elem->node_ptr(2));
                          subelem[1]->set_node(3, elem->node_ptr(5));

                          subelem[2]->set_node(0, elem->node_ptr(0));
                          subelem[2]->set_node(1, elem->node_ptr(1));
                          subelem[2]->set_node(2, elem->node_ptr(2));
                          subelem[2]->set_node(3, elem->node_ptr(4));
                        }
                    }
                  else // Split on 2-3 diagonal
                    {
                      libmesh_assert (split_first_diagonal(elem, 2,3, 0,5));

                      // 0-4 and 2-3 split implies 2-4 split
                      libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                      subelem[0]->set_node(0, elem->node_ptr(0));
                      subelem[0]->set_node(1, elem->node_ptr(4));
                      subelem[0]->set_node(2, elem->node_ptr(2));
                      subelem[0]->set_node(3, elem->node_ptr(3));

                      subelem[1]->set_node(0, elem->node_ptr(3));
                      subelem[1]->set_node(1, elem->node_ptr(4));
                      subelem[1]->set_node(2, elem->node_ptr(2));
                      subelem[1]->set_node(3, elem->node_ptr(5));

                      subelem[2]->set_node(0, elem->node_ptr(0));
                      subelem[2]->set_node(1, elem->node_ptr(1));
                      subelem[2]->set_node(2, elem->node_ptr(2));
                      subelem[2]->set_node(3, elem->node_ptr(4));
                    }
                }
              else // Split on 1-3 diagonal
                {
                  libmesh_assert (split_first_diagonal(elem, 1,3, 0,4));

                  // Split on 0-5 diagonal
                  if (split_first_diagonal(elem, 0,5, 2,3))
                    {
                      // 1-3 and 0-5 split implies 1-5 split
                      libmesh_assert (split_first_diagonal(elem, 1,5, 2,4));

                      subelem[0]->set_node(0, elem->node_ptr(1));
                      subelem[0]->set_node(1, elem->node_ptr(3));
                      subelem[0]->set_node(2, elem->node_ptr(4));
                      subelem[0]->set_node(3, elem->node_ptr(5));

                      subelem[1]->set_node(0, elem->node_ptr(1));
                      subelem[1]->set_node(1, elem->node_ptr(0));
                      subelem[1]->set_node(2, elem->node_ptr(3));
                      subelem[1]->set_node(3, elem->node_ptr(5));

                      subelem[2]->set_node(0, elem->node_ptr(0));
                      subelem[2]->set_node(1, elem->node_ptr(1));
                      subelem[2]->set_node(2, elem->node_ptr(2));
                      subelem[2]->set_node(3, elem->node_ptr(5));
                    }
                  else // Split on 2-3 diagonal
                    {
                      libmesh_assert (split_first_diagonal(elem, 2,3, 0,5));

                      // Split on 1-5 diagonal
                      if (split_first_diagonal(elem, 1,5, 2,4))
                        {
                          subelem[0]->set_node(0, elem->node_ptr(0));
                          subelem[0]->set_node(1, elem->node_ptr(1));
                          subelem[0]->set_node(2, elem->node_ptr(2));
                          subelem[0]->set_node(3, elem->node_ptr(3));

                          subelem[1]->set_node(0, elem->node_ptr(3));
                          subelem[1]->set_node(1, elem->node_ptr(1));
                          subelem[1]->set_node(2, elem->node_ptr(2));
                          subelem[1]->set_node(3, elem->node_ptr(5));

                          subelem[2]->set_node(0, elem->node_ptr(1));
                          subelem[2]->set_node(1, elem->node_ptr(3));
                          subelem[2]->set_node(2, elem->node_ptr(4));
                          subelem[2]->set_node(3, elem->node_ptr(5));
                        }
                      else // Split on 2-4 diagonal
                        {
                          libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                          subelem[0]->set_node(0, elem->node_ptr(0));
                          subelem[0]->set_node(1, elem->node_ptr(1));
                          subelem[0]->set_node(2, elem->node_ptr(2));
                          subelem[0]->set_node(3, elem->node_ptr(3));

                          subelem[1]->set_node(0, elem->node_ptr(2));
                          subelem[1]->set_node(1, elem->node_ptr(3));
                          subelem[1]->set_node(2, elem->node_ptr(4));
                          subelem[1]->set_node(3, elem->node_ptr(5));

                          subelem[2]->set_node(0, elem->node_ptr(3));
                          subelem[2]->set_node(1, elem->node_ptr(1));
                          subelem[2]->set_node(2, elem->node_ptr(2));
                          subelem[2]->set_node(3, elem->node_ptr(4));
                        }
                    }
                }

              break;
            }

          case PRISM20:
          case PRISM21:
            libmesh_experimental(); // We should upgrade this to TET14...
            libmesh_fallthrough();
          case PRISM18:
            {
              subelem[0] = Elem::build(TET10);
              subelem[1] = Elem::build(TET10);
              subelem[2] = Elem::build(TET10);

              // Split on 0-4 diagonal
              if (split_first_diagonal(elem, 0,4, 1,3))
                {
                  // Split on 0-5 diagonal
                  if (split_first_diagonal(elem, 0,5, 2,3))
                    {
                      // Split on 1-5 diagonal
                      if (split_first_diagonal(elem, 1,5, 2,4))
                        {
                          subelem[0]->set_node(0, elem->node_ptr(0));
                          subelem[0]->set_node(1, elem->node_ptr(4));
                          subelem[0]->set_node(2, elem->node_ptr(5));
                          subelem[0]->set_node(3, elem->node_ptr(3));

                          subelem[0]->set_node(4, elem->node_ptr(15));
                          subelem[0]->set_node(5, elem->node_ptr(13));
                          subelem[0]->set_node(6, elem->node_ptr(17));
                          subelem[0]->set_node(7, elem->node_ptr(9));
                          subelem[0]->set_node(8, elem->node_ptr(12));
                          subelem[0]->set_node(9, elem->node_ptr(14));

                          subelem[1]->set_node(0, elem->node_ptr(0));
                          subelem[1]->set_node(1, elem->node_ptr(4));
                          subelem[1]->set_node(2, elem->node_ptr(1));
                          subelem[1]->set_node(3, elem->node_ptr(5));

                          subelem[1]->set_node(4, elem->node_ptr(15));
                          subelem[1]->set_node(5, elem->node_ptr(10));
                          subelem[1]->set_node(6, elem->node_ptr(6));
                          subelem[1]->set_node(7, elem->node_ptr(17));
                          subelem[1]->set_node(8, elem->node_ptr(13));
                          subelem[1]->set_node(9, elem->node_ptr(16));

                          subelem[2]->set_node(0, elem->node_ptr(0));
                          subelem[2]->set_node(1, elem->node_ptr(1));
                          subelem[2]->set_node(2, elem->node_ptr(2));
                          subelem[2]->set_node(3, elem->node_ptr(5));

                          subelem[2]->set_node(4, elem->node_ptr(6));
                          subelem[2]->set_node(5, elem->node_ptr(7));
                          subelem[2]->set_node(6, elem->node_ptr(8));
                          subelem[2]->set_node(7, elem->node_ptr(17));
                          subelem[2]->set_node(8, elem->node_ptr(16));
                          subelem[2]->set_node(9, elem->node_ptr(11));
                        }
                      else // Split on 2-4 diagonal
                        {
                          libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                          subelem[0]->set_node(0, elem->node_ptr(0));
                          subelem[0]->set_node(1, elem->node_ptr(4));
                          subelem[0]->set_node(2, elem->node_ptr(5));
                          subelem[0]->set_node(3, elem->node_ptr(3));

                          subelem[0]->set_node(4, elem->node_ptr(15));
                          subelem[0]->set_node(5, elem->node_ptr(13));
                          subelem[0]->set_node(6, elem->node_ptr(17));
                          subelem[0]->set_node(7, elem->node_ptr(9));
                          subelem[0]->set_node(8, elem->node_ptr(12));
                          subelem[0]->set_node(9, elem->node_ptr(14));

                          subelem[1]->set_node(0, elem->node_ptr(0));
                          subelem[1]->set_node(1, elem->node_ptr(4));
                          subelem[1]->set_node(2, elem->node_ptr(2));
                          subelem[1]->set_node(3, elem->node_ptr(5));

                          subelem[1]->set_node(4, elem->node_ptr(15));
                          subelem[1]->set_node(5, elem->node_ptr(16));
                          subelem[1]->set_node(6, elem->node_ptr(8));
                          subelem[1]->set_node(7, elem->node_ptr(17));
                          subelem[1]->set_node(8, elem->node_ptr(13));
                          subelem[1]->set_node(9, elem->node_ptr(11));

                          subelem[2]->set_node(0, elem->node_ptr(0));
                          subelem[2]->set_node(1, elem->node_ptr(1));
                          subelem[2]->set_node(2, elem->node_ptr(2));
                          subelem[2]->set_node(3, elem->node_ptr(4));

                          subelem[2]->set_node(4, elem->node_ptr(6));
                          subelem[2]->set_node(5, elem->node_ptr(7));
                          subelem[2]->set_node(6, elem->node_ptr(8));
                          subelem[2]->set_node(7, elem->node_ptr(15));
                          subelem[2]->set_node(8, elem->node_ptr(10));
                          subelem[2]->set_node(9, elem->node_ptr(16));
                        }
                    }
                  else // Split on 2-3 diagonal
                    {
                      libmesh_assert (split_first_diagonal(elem, 2,3, 0,5));

                      // 0-4 and 2-3 split implies 2-4 split
                      libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                      subelem[0]->set_node(0, elem->node_ptr(0));
                      subelem[0]->set_node(1, elem->node_ptr(4));
                      subelem[0]->set_node(2, elem->node_ptr(2));
                      subelem[0]->set_node(3, elem->node_ptr(3));

                      subelem[0]->set_node(4, elem->node_ptr(15));
                      subelem[0]->set_node(5, elem->node_ptr(16));
                      subelem[0]->set_node(6, elem->node_ptr(8));
                      subelem[0]->set_node(7, elem->node_ptr(9));
                      subelem[0]->set_node(8, elem->node_ptr(12));
                      subelem[0]->set_node(9, elem->node_ptr(17));

                      subelem[1]->set_node(0, elem->node_ptr(3));
                      subelem[1]->set_node(1, elem->node_ptr(4));
                      subelem[1]->set_node(2, elem->node_ptr(2));
                      subelem[1]->set_node(3, elem->node_ptr(5));

                      subelem[1]->set_node(4, elem->node_ptr(12));
                      subelem[1]->set_node(5, elem->node_ptr(16));
                      subelem[1]->set_node(6, elem->node_ptr(17));
                      subelem[1]->set_node(7, elem->node_ptr(14));
                      subelem[1]->set_node(8, elem->node_ptr(13));
                      subelem[1]->set_node(9, elem->node_ptr(11));

                      subelem[2]->set_node(0, elem->node_ptr(0));
                      subelem[2]->set_node(1, elem->node_ptr(1));
                      subelem[2]->set_node(2, elem->node_ptr(2));
                      subelem[2]->set_node(3, elem->node_ptr(4));

                      subelem[2]->set_node(4, elem->node_ptr(6));
                      subelem[2]->set_node(5, elem->node_ptr(7));
                      subelem[2]->set_node(6, elem->node_ptr(8));
                      subelem[2]->set_node(7, elem->node_ptr(15));
                      subelem[2]->set_node(8, elem->node_ptr(10));
                      subelem[2]->set_node(9, elem->node_ptr(16));
                    }
                }
              else // Split on 1-3 diagonal
                {
                  libmesh_assert (split_first_diagonal(elem, 1,3, 0,4));

                  // Split on 0-5 diagonal
                  if (split_first_diagonal(elem, 0,5, 2,3))
                    {
                      // 1-3 and 0-5 split implies 1-5 split
                      libmesh_assert (split_first_diagonal(elem, 1,5, 2,4));

                      subelem[0]->set_node(0, elem->node_ptr(1));
                      subelem[0]->set_node(1, elem->node_ptr(3));
                      subelem[0]->set_node(2, elem->node_ptr(4));
                      subelem[0]->set_node(3, elem->node_ptr(5));

                      subelem[0]->set_node(4, elem->node_ptr(15));
                      subelem[0]->set_node(5, elem->node_ptr(12));
                      subelem[0]->set_node(6, elem->node_ptr(10));
                      subelem[0]->set_node(7, elem->node_ptr(16));
                      subelem[0]->set_node(8, elem->node_ptr(14));
                      subelem[0]->set_node(9, elem->node_ptr(13));

                      subelem[1]->set_node(0, elem->node_ptr(1));
                      subelem[1]->set_node(1, elem->node_ptr(0));
                      subelem[1]->set_node(2, elem->node_ptr(3));
                      subelem[1]->set_node(3, elem->node_ptr(5));

                      subelem[1]->set_node(4, elem->node_ptr(6));
                      subelem[1]->set_node(5, elem->node_ptr(9));
                      subelem[1]->set_node(6, elem->node_ptr(15));
                      subelem[1]->set_node(7, elem->node_ptr(16));
                      subelem[1]->set_node(8, elem->node_ptr(17));
                      subelem[1]->set_node(9, elem->node_ptr(14));

                      subelem[2]->set_node(0, elem->node_ptr(0));
                      subelem[2]->set_node(1, elem->node_ptr(1));
                      subelem[2]->set_node(2, elem->node_ptr(2));
                      subelem[2]->set_node(3, elem->node_ptr(5));

                      subelem[2]->set_node(4, elem->node_ptr(6));
                      subelem[2]->set_node(5, elem->node_ptr(7));
                      subelem[2]->set_node(6, elem->node_ptr(8));
                      subelem[2]->set_node(7, elem->node_ptr(17));
                      subelem[2]->set_node(8, elem->node_ptr(16));
                      subelem[2]->set_node(9, elem->node_ptr(11));
                    }
                  else // Split on 2-3 diagonal
                    {
                      libmesh_assert (split_first_diagonal(elem, 2,3, 0,5));

                      // Split on 1-5 diagonal
                      if (split_first_diagonal(elem, 1,5, 2,4))
                        {
                          subelem[0]->set_node(0, elem->node_ptr(0));
                          subelem[0]->set_node(1, elem->node_ptr(1));
                          subelem[0]->set_node(2, elem->node_ptr(2));
                          subelem[0]->set_node(3, elem->node_ptr(3));

                          subelem[0]->set_node(4, elem->node_ptr(6));
                          subelem[0]->set_node(5, elem->node_ptr(7));
                          subelem[0]->set_node(6, elem->node_ptr(8));
                          subelem[0]->set_node(7, elem->node_ptr(9));
                          subelem[0]->set_node(8, elem->node_ptr(15));
                          subelem[0]->set_node(9, elem->node_ptr(17));

                          subelem[1]->set_node(0, elem->node_ptr(3));
                          subelem[1]->set_node(1, elem->node_ptr(1));
                          subelem[1]->set_node(2, elem->node_ptr(2));
                          subelem[1]->set_node(3, elem->node_ptr(5));

                          subelem[1]->set_node(4, elem->node_ptr(15));
                          subelem[1]->set_node(5, elem->node_ptr(7));
                          subelem[1]->set_node(6, elem->node_ptr(17));
                          subelem[1]->set_node(7, elem->node_ptr(14));
                          subelem[1]->set_node(8, elem->node_ptr(16));
                          subelem[1]->set_node(9, elem->node_ptr(11));

                          subelem[2]->set_node(0, elem->node_ptr(1));
                          subelem[2]->set_node(1, elem->node_ptr(3));
                          subelem[2]->set_node(2, elem->node_ptr(4));
                          subelem[2]->set_node(3, elem->node_ptr(5));

                          subelem[2]->set_node(4, elem->node_ptr(15));
                          subelem[2]->set_node(5, elem->node_ptr(12));
                          subelem[2]->set_node(6, elem->node_ptr(10));
                          subelem[2]->set_node(7, elem->node_ptr(16));
                          subelem[2]->set_node(8, elem->node_ptr(14));
                          subelem[2]->set_node(9, elem->node_ptr(13));
                        }
                      else // Split on 2-4 diagonal
                        {
                          libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                          subelem[0]->set_node(0, elem->node_ptr(0));
                          subelem[0]->set_node(1, elem->node_ptr(1));
                          subelem[0]->set_node(2, elem->node_ptr(2));
                          subelem[0]->set_node(3, elem->node_ptr(3));

                          subelem[0]->set_node(4, elem->node_ptr(6));
                          subelem[0]->set_node(5, elem->node_ptr(7));
                          subelem[0]->set_node(6, elem->node_ptr(8));
                          subelem[0]->set_node(7, elem->node_ptr(9));
                          subelem[0]->set_node(8, elem->node_ptr(15));
                          subelem[0]->set_node(9, elem->node_ptr(17));

                          subelem[1]->set_node(0, elem->node_ptr(2));
                          subelem[1]->set_node(1, elem->node_ptr(3));
                          subelem[1]->set_node(2, elem->node_ptr(4));
                          subelem[1]->set_node(3, elem->node_ptr(5));

                          subelem[1]->set_node(4, elem->node_ptr(17));
                          subelem[1]->set_node(5, elem->node_ptr(12));
                          subelem[1]->set_node(6, elem->node_ptr(16));
                          subelem[1]->set_node(7, elem->node_ptr(11));
                          subelem[1]->set_node(8, elem->node_ptr(14));
                          subelem[1]->set_node(9, elem->node_ptr(13));

                          subelem[2]->set_node(0, elem->node_ptr(3));
                          subelem[2]->set_node(1, elem->node_ptr(1));
                          subelem[2]->set_node(2, elem->node_ptr(2));
                          subelem[2]->set_node(3, elem->node_ptr(4));

                          subelem[2]->set_node(4, elem->node_ptr(15));
                          subelem[2]->set_node(5, elem->node_ptr(7));
                          subelem[2]->set_node(6, elem->node_ptr(17));
                          subelem[2]->set_node(7, elem->node_ptr(12));
                          subelem[2]->set_node(8, elem->node_ptr(10));
                          subelem[2]->set_node(9, elem->node_ptr(16));
                        }
                    }
                }

              break;
            }

            // No need to split elements that are already simplicial:
          case EDGE2:
          case EDGE3:
          case EDGE4:
          case TRI3:
          case TRI6:
          case TRI7:
          case TET4:
          case TET10:
          case TET14:
          case INFEDGE2:
            // No way to split infinite quad/prism elements, so
            // hopefully no need to
          case INFQUAD4:
          case INFQUAD6:
          case INFPRISM6:
          case INFPRISM12:
            continue;
            // If we're left with an unimplemented hex we're probably
            // out of luck.  TODO: implement hexes
          default:
            {
              libMesh::err << "Error, encountered unimplemented element "
                           << Utility::enum_to_string<ElemType>(etype)
                           << " in MeshTools::Modification::all_tri()..."
                           << std::endl;
              libmesh_not_implemented();
            }
          } // end switch (etype)

        // Be sure the correct data is set for all subelems.
        const unsigned int nei = elem->n_extra_integers();
        for (unsigned int i=0; i != max_subelems; ++i)
          if (subelem[i]) {
            subelem[i]->processor_id() = elem->processor_id();
            subelem[i]->subdomain_id() = elem->subdomain_id();

            // Copy any extra element data.  Since the subelements
            // haven't been added to the mesh yet any allocation has
            // to be done manually.
            subelem[i]->add_extra_integers(nei);
            for (unsigned int ei=0; ei != nei; ++ei)
              subelem[ei]->set_extra_integer(ei, elem->get_extra_integer(ei));


            // Copy any mapping data.
            subelem[i]->set_mapping_type(elem->mapping_type());
            subelem[i]->set_mapping_data(elem->mapping_data());
          }

        // On a mesh with boundary data, we need to move that data to
        // the new elements.

        // On a mesh which is distributed, we need to move
        // remote_elem links to the new elements.
        bool mesh_is_serial = mesh.is_serial();

        if (mesh_has_boundary_data || !mesh_is_serial)
          {
            // Container to key boundary IDs handed back by the BoundaryInfo object.
            std::vector<boundary_id_type> bc_ids;

            for (auto sn : elem->side_index_range())
              {
                mesh.get_boundary_info().boundary_ids(elem, sn, bc_ids);

                if (bc_ids.empty() && elem->neighbor_ptr(sn) != remote_elem)
                  continue;

                // Make a sorted list of node ids for elem->side(sn)
                elem->build_side_ptr(elem_side, sn);
                std::vector<dof_id_type> elem_side_nodes(elem_side->n_nodes());
                for (unsigned int esn=0,
                     n_esn = cast_int<unsigned int>(elem_side_nodes.size());
                     esn != n_esn; ++esn)
                  elem_side_nodes[esn] = elem_side->node_id(esn);
                std::sort(elem_side_nodes.begin(), elem_side_nodes.end());

                for (unsigned int i=0; i != max_subelems; ++i)
                  if (subelem[i])
                    {
                      for (auto subside : subelem[i]->side_index_range())
                        {
                          subelem[i]->build_side_ptr(subside_elem, subside);

                          // Make a list of *vertex* node ids for this subside, see if they are all present
                          // in elem->side(sn).  Note 1: we can't just compare elem->key(sn) to
                          // subelem[i]->key(subside) in the Prism cases, since the new side is
                          // a different type.  Note 2: we only use vertex nodes since, in the future,
                          // a Hex20 or Prism15's QUAD8 face may be split into two Tri6 faces, and the
                          // original face will not contain the mid-edge node.
                          std::vector<dof_id_type> subside_nodes(subside_elem->n_vertices());
                          for (unsigned int ssn=0,
                               n_ssn = cast_int<unsigned int>(subside_nodes.size());
                               ssn != n_ssn; ++ssn)
                            subside_nodes[ssn] = subside_elem->node_id(ssn);
                          std::sort(subside_nodes.begin(), subside_nodes.end());

                          // std::includes returns true if every element of the second sorted range is
                          // contained in the first sorted range.
                          if (std::includes(elem_side_nodes.begin(), elem_side_nodes.end(),
                                            subside_nodes.begin(), subside_nodes.end()))
                            {
                              for (const auto & b_id : bc_ids)
                                if (b_id != BoundaryInfo::invalid_id)
                                  {
                                    new_bndry_ids.push_back(b_id);
                                    new_bndry_elements.push_back(subelem[i].get());
                                    new_bndry_sides.push_back(subside);
                                  }

                              // If the original element had a RemoteElem neighbor on side 'sn',
                              // then the subelem has one on side 'subside'.
                              if (elem->neighbor_ptr(sn) == remote_elem)
                                subelem[i]->set_neighbor(subside, const_cast<RemoteElem*>(remote_elem));
                            }
                        }
                    } // end for loop over subelem
              } // end for loop over sides

            // Remove the original element from the BoundaryInfo structure.
            mesh.get_boundary_info().remove(elem);

          } // end if (mesh_has_boundary_data)

        // Determine new IDs for the split elements which will be
        // the same on all processors, therefore keeping the Mesh
        // in sync.  Note: we offset the new IDs by max_orig_id to
        // avoid overwriting any of the original IDs.
        for (unsigned int i=0; i != max_subelems; ++i)
          if (subelem[i])
            {
              // Determine new IDs for the split elements which will be
              // the same on all processors, therefore keeping the Mesh
              // in sync.  Note: we offset the new IDs by the max of the
              // pre-existing ids to avoid conflicting with originals.
              subelem[i]->set_id( max_orig_id + 6*elem->id() + i );

#ifdef LIBMESH_ENABLE_UNIQUE_ID
              subelem[i]->set_unique_id(max_unique_id + max_subelems*elem->unique_id() + i);
#endif

              // Prepare to add the newly-created simplices
              new_elements.push_back(std::move(subelem[i]));
            }

        // Delete the original element
        mesh.delete_elem(elem);
      } // End for loop over elements
  } // end scope


  // Now, iterate over the new elements vector, and add them each to
  // the Mesh.
  for (auto & elem : new_elements)
    mesh.add_elem(std::move(elem));

  if (mesh_has_boundary_data)
    {
      // If the old mesh had boundary data, the new mesh better have
      // some.  However, we can't assert that the size of
      // new_bndry_elements vector is > 0, since we may not have split
      // any elements actually on the boundary.  We also can't assert
      // that the original number of boundary sides is equal to the
      // sum of the boundary sides currently in the mesh and the
      // newly-added boundary sides, since in 3D, we may have split a
      // boundary QUAD into two boundary TRIs.  Therefore, we won't be
      // too picky about the actual number of BCs, and just assert that
      // there are some, somewhere.
#ifndef NDEBUG
      bool nbe_nonempty = new_bndry_elements.size();
      mesh.comm().max(nbe_nonempty);
      libmesh_assert(nbe_nonempty ||
                     mesh.get_boundary_info().n_boundary_conds()>0);
#endif

      // We should also be sure that the lengths of the new boundary data vectors
      // are all the same.
      libmesh_assert_equal_to (new_bndry_elements.size(), new_bndry_sides.size());
      libmesh_assert_equal_to (new_bndry_sides.size(), new_bndry_ids.size());

      // Add the new boundary info to the mesh
      for (auto s : index_range(new_bndry_elements))
        mesh.get_boundary_info().add_side(new_bndry_elements[s],
                                          new_bndry_sides[s],
                                          new_bndry_ids[s]);
    }

  // In a DistributedMesh any newly added ghost node ids may be
  // inconsistent, and unique_ids of newly added ghost nodes remain
  // unset.
  // make_nodes_parallel_consistent() will fix all this.
  if (!mesh.is_serial())
    {
      mesh.comm().max(added_new_ghost_point);

      if (added_new_ghost_point)
        MeshCommunication().make_nodes_parallel_consistent (mesh);
    }

  // Prepare the newly created mesh for use.
  mesh.prepare_for_use();

  // Let the new_elements and new_bndry_elements vectors go out of scope.
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

  // For avoiding extraneous element side allocation
  ElemSideBuilder side_builder;

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
                            const Elem & side = side_builder(*elem, s);
                            const Node & node0 = side.node_ref(0);
                            const Node & node1 = side.node_ref(1);

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
                            const Real em_val = parent->embedding_matrix(c,nc,n);

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
            for (auto nid : make_range(mesh.n_nodes()))
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

  // We haven't changed any topology, but just changing geometry could
  // have invalidated a point locator.
  mesh.clear_point_locator();
}



#ifdef LIBMESH_ENABLE_AMR
void MeshTools::Modification::flatten(MeshBase & mesh)
{
  libmesh_assert(mesh.is_prepared() || mesh.is_replicated());

  // Algorithm:
  // .) For each active element in the mesh: construct a
  //    copy which is the same in every way *except* it is
  //    a level 0 element.  Store the pointers to these in
  //    a separate vector. Save any boundary information as well.
  //    Delete the active element from the mesh.
  // .) Loop over all (remaining) elements in the mesh, delete them.
  // .) Add the level-0 copies back to the mesh

  // Temporary storage for new element pointers
  std::vector<std::unique_ptr<Elem>> new_elements;

  // BoundaryInfo Storage for element ids, sides, and BC ids
  std::vector<Elem *>              saved_boundary_elements;
  std::vector<boundary_id_type>   saved_bc_ids;
  std::vector<unsigned short int> saved_bc_sides;

  // Container to catch boundary ids passed back by BoundaryInfo
  std::vector<boundary_id_type> bc_ids;

  // Reserve a reasonable amt. of space for each
  new_elements.reserve(mesh.n_active_elem());
  saved_boundary_elements.reserve(mesh.get_boundary_info().n_boundary_conds());
  saved_bc_ids.reserve(mesh.get_boundary_info().n_boundary_conds());
  saved_bc_sides.reserve(mesh.get_boundary_info().n_boundary_conds());

  for (auto & elem : mesh.active_element_ptr_range())
    {
      // Make a new element of the same type
      auto copy = Elem::build(elem->type());

      // Set node pointers (they still point to nodes in the original mesh)
      for (auto n : elem->node_index_range())
        copy->set_node(n, elem->node_ptr(n));

      // Copy over ids
      copy->processor_id() = elem->processor_id();
      copy->subdomain_id() = elem->subdomain_id();

      // Retain the original element's ID(s) as well, otherwise
      // the Mesh may try to create them for you...
      copy->set_id( elem->id() );
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      copy->set_unique_id(elem->unique_id());
#endif

      // This element could have boundary info or DistributedMesh
      // remote_elem links as well.  We need to save the (elem,
      // side, bc_id) triples and those links
      for (auto s : elem->side_index_range())
        {
          if (elem->neighbor_ptr(s) == remote_elem)
            copy->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));

          mesh.get_boundary_info().boundary_ids(elem, s, bc_ids);
          for (const auto & bc_id : bc_ids)
            if (bc_id != BoundaryInfo::invalid_id)
              {
                saved_boundary_elements.push_back(copy.get());
                saved_bc_ids.push_back(bc_id);
                saved_bc_sides.push_back(s);
              }
        }

      // Copy any extra element data.  Since the copy hasn't been
      // added to the mesh yet any allocation has to be done manually.
      const unsigned int nei = elem->n_extra_integers();
      copy->add_extra_integers(nei);
      for (unsigned int i=0; i != nei; ++i)
        copy->set_extra_integer(i, elem->get_extra_integer(i));

      // Copy any mapping data.
      copy->set_mapping_type(elem->mapping_type());
      copy->set_mapping_data(elem->mapping_data());

      // We're done with this element
      mesh.delete_elem(elem);

      // But save the copy
      new_elements.push_back(std::move(copy));
    }

  // Make sure we saved the same number of boundary conditions
  // in each vector.
  libmesh_assert_equal_to (saved_boundary_elements.size(), saved_bc_ids.size());
  libmesh_assert_equal_to (saved_bc_ids.size(), saved_bc_sides.size());

  // Loop again, delete any remaining elements
  for (auto & elem : mesh.element_ptr_range())
    mesh.delete_elem(elem);

  // Add the copied (now level-0) elements back to the mesh
  for (auto & new_elem : new_elements)
    {
      // Save the original ID, because the act of adding the Elem can
      // change new_elem's id!
      dof_id_type orig_id = new_elem->id();

      Elem * added_elem = mesh.add_elem(std::move(new_elem));

      // If the Elem, as it was re-added to the mesh, now has a
      // different ID (this is unlikely, so it's just an assert)
      // the boundary information will no longer be correct.
      libmesh_assert_equal_to (orig_id, added_elem->id());

      // Avoid compiler warnings in opt mode.
      libmesh_ignore(added_elem, orig_id);
    }

  // Finally, also add back the saved boundary information
  for (auto e : index_range(saved_boundary_elements))
    mesh.get_boundary_info().add_side(saved_boundary_elements[e],
                                      saved_bc_sides[e],
                                      saved_bc_ids[e]);

  // Trim unused and renumber nodes and elements
  mesh.prepare_for_use();
}
#endif // #ifdef LIBMESH_ENABLE_AMR



void MeshTools::Modification::change_boundary_id (MeshBase & mesh,
                                                  const boundary_id_type old_id,
                                                  const boundary_id_type new_id)
{
  // This is just a shim around the member implementation, now
  mesh.get_boundary_info().renumber_id(old_id, new_id);
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

  // We just invalidated mesh.get_subdomain_ids(), but it might not be
  // efficient to fix that here.
  mesh.set_hasnt_cached_elem_data();
}


} // namespace libMesh
