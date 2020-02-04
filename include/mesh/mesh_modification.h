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



#ifndef LIBMESH_MESH_MODIFICATION_H
#define LIBMESH_MESH_MODIFICATION_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/id_types.h" // for boundary_id_type, subdomain_id_type
#include "libmesh/mesh_base.h"
#include "libmesh/function_base.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/remote_elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/node.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_tet10.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/parallel.h"
#include "libmesh/mesh_communication.h"

#include <vector>

namespace libMesh
{


// ------------------------------------------------------------
// MeshTools::Modification namespace
namespace MeshTools
{
/**
 * Tools for \p Mesh modification.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
namespace Modification
{
/**
 * Randomly perturb the nodal locations.  This function will
 * move each node \p factor fraction of its minimum neighboring
 * node separation distance.  Nodes on the boundary are not moved
 * by default, however they may be by setting the flag
 * \p perturb_boundary true.
 */
void distort (MeshBase & mesh,
              const Real factor, const bool perturb_boundary=false);

/**
 * Deterministically perturb the nodal locations.  This function will
 * move each node from it's current x/y/z coordinates to a new x/y/z
 * coordinate given by the first LIBMESH_DIM components of the
 * specified function \p mapfunc
 *
 * Nodes on the boundary are also moved.
 *
 * Currently, non-vertex nodes are moved in the same way as vertex
 * nodes, according to (newx,newy,newz) = mapfunc(x,y,z).  This
 * behavior is often suboptimal for higher order geometries and may be
 * subject to change in future libMesh versions.
 */
template <typename RealType>
void redistribute (MeshBaseTempl<RealType> & mesh,
                                            const FunctionBase<RealType> & mapfunc)
{
  libmesh_assert (mesh.n_nodes());
  libmesh_assert (mesh.n_elem());

  LOG_SCOPE("redistribute()", "MeshTools::Modification");

  DenseVector<RealType> output_vec(LIBMESH_DIM);

  // FIXME - we should thread this later.
  std::unique_ptr<FunctionBase<RealType>> myfunc = mapfunc.clone();

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
}


/**
 * Translates the mesh.  The grid points are translated in the
 * \p x direction by \p xt, in the \p y direction by \p yt,
 * etc...
 */
void translate (MeshBase & mesh,
                const Real xt=0., const Real yt=0., const Real zt=0.);

//     /**
//      * Rotates the mesh in the xy plane. The rotation is
//      * counter-clock-wise (mathematical definition).
//      * The angle is in degrees (360 make a full circle)
//      */
//     void rotate2D (MeshBase & mesh,
//                    const Real alpha=0.);

/**
 * Rotates the mesh in 3D space.
 * Here the standard Euler angles are adopted
 * (http://mathworld.wolfram.com/EulerAngles.html)
 * The angles are in degrees (360 make a full circle)
 */
void rotate (MeshBase & mesh,
             const Real phi, const Real theta=0., const Real psi=0.);

/**
 * Scales the mesh.  The grid points are scaled in the
 * \p x direction by \p xs, in the \p y direction by \p ys,
 * etc...  If only \p xs is specified then the scaling is
 * assumed uniform in all directions.
 */
void scale (MeshBase & mesh,
            const Real xs, const Real ys=0., const Real zs=0.);


/**
 * Smooth the mesh with a simple Laplace smoothing algorithm.  The mesh is
 * smoothed \p n_iterations times.  If the parameter \p power is 0, each
 * node is moved to the average position of the neighboring connected
 * nodes. If \p power > 0, the node positions are weighted by their
 * distance.  The positions of higher order nodes, and nodes living in
 * refined elements, are calculated from the vertex positions of their
 * parent nodes.  Only works in 2D.
 *
 * \author Martin Luthi (luthi@gi.alaska.edu)
 * \date 2005
 */
void smooth(MeshBase &, unsigned int, Real);

#ifdef LIBMESH_ENABLE_AMR
/**
 * Removes all the refinement tree structure of Mesh, leaving
 * only the highest-level (most-refined) elements.  This is useful
 * when you want to write out a uniformly-refined grid to be treated later
 * as an initial mesh.
 *
 * \note Many functions in LibMesh assume a conforming (with no
 * hanging nodes) grid exists at some level, so you probably only want
 * to do this on meshes which have been uniformly refined.
 */
template <typename RealType>
void flatten(MeshBaseTempl<RealType> & mesh)
{
  typedef ElemTempl<RealType> Elem;
  typedef RemoteElemTempl<RealType> RemoteElem;

  // Algorithm:
  // .) For each active element in the mesh: construct a
  //    copy which is the same in every way *except* it is
  //    a level 0 element.  Store the pointers to these in
  //    a separate vector. Save any boundary information as well.
  //    Delete the active element from the mesh.
  // .) Loop over all (remaining) elements in the mesh, delete them.
  // .) Add the level-0 copies back to the mesh

  // Temporary storage for new element pointers
  std::vector<Elem *> new_elements;

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
      Elem * copy = Elem::build(elem->type()).release();

      // Set node pointers (they still point to nodes in the original mesh)
      for (auto n : elem->node_index_range())
        copy->set_node(n) = elem->node_ptr(n);

      // Copy over ids
      copy->processor_id() = elem->processor_id();
      copy->subdomain_id() = elem->subdomain_id();

      // Retain the original element's ID(s) as well, otherwise
      // the Mesh may try to create them for you...
      copy->set_id( elem->id() );
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      copy->set_unique_id() = elem->unique_id();
#endif

      // This element could have boundary info or DistributedMesh
      // remote_elem links as well.  We need to save the (elem,
      // side, bc_id) triples and those links
      for (auto s : elem->side_index_range())
        {
          if (elem->neighbor_ptr(s) == RemoteElem::get_instance())
            copy->set_neighbor(s, const_cast<RemoteElem *>(RemoteElem::get_instance()));

          mesh.get_boundary_info().boundary_ids(elem, s, bc_ids);
          for (const auto & bc_id : bc_ids)
            if (bc_id != BoundaryInfo::invalid_id)
              {
                saved_boundary_elements.push_back(copy);
                saved_bc_ids.push_back(bc_id);
                saved_bc_sides.push_back(s);
              }
        }


      // We're done with this element
      mesh.delete_elem(elem);

      // But save the copy
      new_elements.push_back(copy);
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

      Elem * added_elem = mesh.add_elem(new_elem);

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
  mesh.prepare_for_use(/*skip_renumber =*/ false);
}
#endif // #ifdef LIBMESH_ENABLE_AMR

/**
 * Finds any boundary ids that are currently old_id,
 * changes them to new_id
 */
void change_boundary_id (MeshBase & mesh,
                         const boundary_id_type old_id,
                         const boundary_id_type new_id);

/**
 * Finds any subdomain ids that are currently old_id,
 * changes them to new_id
 */
void change_subdomain_id (MeshBase & mesh,
                          const subdomain_id_type old_id,
                          const subdomain_id_type new_id);


namespace
{
template <typename RealType>
bool split_first_diagonal(const libMesh::ElemTempl<RealType> * elem,
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

}

/**
 * Converts the 2D quadrilateral elements of a Mesh into
 * triangular elements.
 *
 * \note Only works for 2D elements!  3D elements are ignored.
 * \note Probably won't do the right thing for meshes which
 * have been refined previously.
 */
template <typename RealType>
void all_tri (MeshBaseTempl<RealType> & mesh)
{
  typedef ElemTempl<RealType> Elem;
  typedef Tri3Templ<RealType> Tri3;
  typedef Tri6Templ<RealType> Tri6;
  typedef Tet4Templ<RealType> Tet4;
  typedef Tet10Templ<RealType> Tet10;
  typedef NodeTempl<RealType> Node;
  typedef BoundaryInfoTempl<RealType> BoundaryInfo;
  typedef RemoteElemTempl<RealType> RemoteElem;

  // The number of elements in the original mesh before any additions
  // or deletions.
  const dof_id_type n_orig_elem = mesh.n_elem();
  const dof_id_type max_orig_id = mesh.max_elem_id();

  // We store pointers to the newly created elements in a vector
  // until they are ready to be added to the mesh.  This is because
  // adding new elements on the fly can cause reallocation and invalidation
  // of existing mesh element_iterators.
  std::vector<Elem *> new_elements;

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

    for (auto & elem : mesh.element_ptr_range())
      {
        const ElemType etype = elem->type();

        // all_tri currently only works on coarse meshes
        libmesh_assert (!elem->parent());

        // The new elements we will split the quad into.
        // In 3D we may need as many as 6 tets per hex
        Elem * subelem[6];

        for (unsigned int i = 0; i != max_subelems; ++i)
          subelem[i] = nullptr;

        switch (etype)
          {
          case QUAD4:
            {
              subelem[0] = new Tri3;
              subelem[1] = new Tri3;

              // Check for possible edge swap
              if ((elem->point(0) - elem->point(2)).norm() <
                  (elem->point(1) - elem->point(3)).norm())
                {
                  subelem[0]->set_node(0) = elem->node_ptr(0);
                  subelem[0]->set_node(1) = elem->node_ptr(1);
                  subelem[0]->set_node(2) = elem->node_ptr(2);

                  subelem[1]->set_node(0) = elem->node_ptr(0);
                  subelem[1]->set_node(1) = elem->node_ptr(2);
                  subelem[1]->set_node(2) = elem->node_ptr(3);
                }

              else
                {
                  subelem[0]->set_node(0) = elem->node_ptr(0);
                  subelem[0]->set_node(1) = elem->node_ptr(1);
                  subelem[0]->set_node(2) = elem->node_ptr(3);

                  subelem[1]->set_node(0) = elem->node_ptr(1);
                  subelem[1]->set_node(1) = elem->node_ptr(2);
                  subelem[1]->set_node(2) = elem->node_ptr(3);
                }


              break;
            }

          case QUAD8:
            {
              if (elem->processor_id() != mesh.processor_id())
                added_new_ghost_point = true;

              subelem[0] = new Tri6;
              subelem[1] = new Tri6;

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
                  subelem[0]->set_node(0) = elem->node_ptr(0);
                  subelem[0]->set_node(1) = elem->node_ptr(1);
                  subelem[0]->set_node(2) = elem->node_ptr(2);
                  subelem[0]->set_node(3) = elem->node_ptr(4);
                  subelem[0]->set_node(4) = elem->node_ptr(5);
                  subelem[0]->set_node(5) = new_node;

                  subelem[1]->set_node(0) = elem->node_ptr(0);
                  subelem[1]->set_node(1) = elem->node_ptr(2);
                  subelem[1]->set_node(2) = elem->node_ptr(3);
                  subelem[1]->set_node(3) = new_node;
                  subelem[1]->set_node(4) = elem->node_ptr(6);
                  subelem[1]->set_node(5) = elem->node_ptr(7);

                }

              else
                {
                  subelem[0]->set_node(0) = elem->node_ptr(3);
                  subelem[0]->set_node(1) = elem->node_ptr(0);
                  subelem[0]->set_node(2) = elem->node_ptr(1);
                  subelem[0]->set_node(3) = elem->node_ptr(7);
                  subelem[0]->set_node(4) = elem->node_ptr(4);
                  subelem[0]->set_node(5) = new_node;

                  subelem[1]->set_node(0) = elem->node_ptr(1);
                  subelem[1]->set_node(1) = elem->node_ptr(2);
                  subelem[1]->set_node(2) = elem->node_ptr(3);
                  subelem[1]->set_node(3) = elem->node_ptr(5);
                  subelem[1]->set_node(4) = elem->node_ptr(6);
                  subelem[1]->set_node(5) = new_node;
                }

              break;
            }

          case QUAD9:
            {
              subelem[0] = new Tri6;
              subelem[1] = new Tri6;

              // Check for possible edge swap
              if ((elem->point(0) - elem->point(2)).norm() <
                  (elem->point(1) - elem->point(3)).norm())
                {
                  subelem[0]->set_node(0) = elem->node_ptr(0);
                  subelem[0]->set_node(1) = elem->node_ptr(1);
                  subelem[0]->set_node(2) = elem->node_ptr(2);
                  subelem[0]->set_node(3) = elem->node_ptr(4);
                  subelem[0]->set_node(4) = elem->node_ptr(5);
                  subelem[0]->set_node(5) = elem->node_ptr(8);

                  subelem[1]->set_node(0) = elem->node_ptr(0);
                  subelem[1]->set_node(1) = elem->node_ptr(2);
                  subelem[1]->set_node(2) = elem->node_ptr(3);
                  subelem[1]->set_node(3) = elem->node_ptr(8);
                  subelem[1]->set_node(4) = elem->node_ptr(6);
                  subelem[1]->set_node(5) = elem->node_ptr(7);
                }

              else
                {
                  subelem[0]->set_node(0) = elem->node_ptr(0);
                  subelem[0]->set_node(1) = elem->node_ptr(1);
                  subelem[0]->set_node(2) = elem->node_ptr(3);
                  subelem[0]->set_node(3) = elem->node_ptr(4);
                  subelem[0]->set_node(4) = elem->node_ptr(8);
                  subelem[0]->set_node(5) = elem->node_ptr(7);

                  subelem[1]->set_node(0) = elem->node_ptr(1);
                  subelem[1]->set_node(1) = elem->node_ptr(2);
                  subelem[1]->set_node(2) = elem->node_ptr(3);
                  subelem[1]->set_node(3) = elem->node_ptr(5);
                  subelem[1]->set_node(4) = elem->node_ptr(6);
                  subelem[1]->set_node(5) = elem->node_ptr(8);
                }

              break;
            }

          case PRISM6:
            {
              // Prisms all split into three tetrahedra
              subelem[0] = new Tet4;
              subelem[1] = new Tet4;
              subelem[2] = new Tet4;

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
                          subelem[0]->set_node(0) = elem->node_ptr(0);
                          subelem[0]->set_node(1) = elem->node_ptr(4);
                          subelem[0]->set_node(2) = elem->node_ptr(5);
                          subelem[0]->set_node(3) = elem->node_ptr(3);

                          subelem[1]->set_node(0) = elem->node_ptr(0);
                          subelem[1]->set_node(1) = elem->node_ptr(4);
                          subelem[1]->set_node(2) = elem->node_ptr(1);
                          subelem[1]->set_node(3) = elem->node_ptr(5);

                          subelem[2]->set_node(0) = elem->node_ptr(0);
                          subelem[2]->set_node(1) = elem->node_ptr(1);
                          subelem[2]->set_node(2) = elem->node_ptr(2);
                          subelem[2]->set_node(3) = elem->node_ptr(5);
                        }
                      else // Split on 2-4 diagonal
                        {
                          libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                          subelem[0]->set_node(0) = elem->node_ptr(0);
                          subelem[0]->set_node(1) = elem->node_ptr(4);
                          subelem[0]->set_node(2) = elem->node_ptr(5);
                          subelem[0]->set_node(3) = elem->node_ptr(3);

                          subelem[1]->set_node(0) = elem->node_ptr(0);
                          subelem[1]->set_node(1) = elem->node_ptr(4);
                          subelem[1]->set_node(2) = elem->node_ptr(2);
                          subelem[1]->set_node(3) = elem->node_ptr(5);

                          subelem[2]->set_node(0) = elem->node_ptr(0);
                          subelem[2]->set_node(1) = elem->node_ptr(1);
                          subelem[2]->set_node(2) = elem->node_ptr(2);
                          subelem[2]->set_node(3) = elem->node_ptr(4);
                        }
                    }
                  else // Split on 2-3 diagonal
                    {
                      libmesh_assert (split_first_diagonal(elem, 2,3, 0,5));

                      // 0-4 and 2-3 split implies 2-4 split
                      libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                      subelem[0]->set_node(0) = elem->node_ptr(0);
                      subelem[0]->set_node(1) = elem->node_ptr(4);
                      subelem[0]->set_node(2) = elem->node_ptr(2);
                      subelem[0]->set_node(3) = elem->node_ptr(3);

                      subelem[1]->set_node(0) = elem->node_ptr(3);
                      subelem[1]->set_node(1) = elem->node_ptr(4);
                      subelem[1]->set_node(2) = elem->node_ptr(2);
                      subelem[1]->set_node(3) = elem->node_ptr(5);

                      subelem[2]->set_node(0) = elem->node_ptr(0);
                      subelem[2]->set_node(1) = elem->node_ptr(1);
                      subelem[2]->set_node(2) = elem->node_ptr(2);
                      subelem[2]->set_node(3) = elem->node_ptr(4);
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

                      subelem[0]->set_node(0) = elem->node_ptr(1);
                      subelem[0]->set_node(1) = elem->node_ptr(3);
                      subelem[0]->set_node(2) = elem->node_ptr(4);
                      subelem[0]->set_node(3) = elem->node_ptr(5);

                      subelem[1]->set_node(0) = elem->node_ptr(1);
                      subelem[1]->set_node(1) = elem->node_ptr(0);
                      subelem[1]->set_node(2) = elem->node_ptr(3);
                      subelem[1]->set_node(3) = elem->node_ptr(5);

                      subelem[2]->set_node(0) = elem->node_ptr(0);
                      subelem[2]->set_node(1) = elem->node_ptr(1);
                      subelem[2]->set_node(2) = elem->node_ptr(2);
                      subelem[2]->set_node(3) = elem->node_ptr(5);
                    }
                  else // Split on 2-3 diagonal
                    {
                      libmesh_assert (split_first_diagonal(elem, 2,3, 0,5));

                      // Split on 1-5 diagonal
                      if (split_first_diagonal(elem, 1,5, 2,4))
                        {
                          subelem[0]->set_node(0) = elem->node_ptr(0);
                          subelem[0]->set_node(1) = elem->node_ptr(1);
                          subelem[0]->set_node(2) = elem->node_ptr(2);
                          subelem[0]->set_node(3) = elem->node_ptr(3);

                          subelem[1]->set_node(0) = elem->node_ptr(3);
                          subelem[1]->set_node(1) = elem->node_ptr(1);
                          subelem[1]->set_node(2) = elem->node_ptr(2);
                          subelem[1]->set_node(3) = elem->node_ptr(5);

                          subelem[2]->set_node(0) = elem->node_ptr(1);
                          subelem[2]->set_node(1) = elem->node_ptr(3);
                          subelem[2]->set_node(2) = elem->node_ptr(4);
                          subelem[2]->set_node(3) = elem->node_ptr(5);
                        }
                      else // Split on 2-4 diagonal
                        {
                          libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                          subelem[0]->set_node(0) = elem->node_ptr(0);
                          subelem[0]->set_node(1) = elem->node_ptr(1);
                          subelem[0]->set_node(2) = elem->node_ptr(2);
                          subelem[0]->set_node(3) = elem->node_ptr(3);

                          subelem[1]->set_node(0) = elem->node_ptr(2);
                          subelem[1]->set_node(1) = elem->node_ptr(3);
                          subelem[1]->set_node(2) = elem->node_ptr(4);
                          subelem[1]->set_node(3) = elem->node_ptr(5);

                          subelem[2]->set_node(0) = elem->node_ptr(3);
                          subelem[2]->set_node(1) = elem->node_ptr(1);
                          subelem[2]->set_node(2) = elem->node_ptr(2);
                          subelem[2]->set_node(3) = elem->node_ptr(4);
                        }
                    }
                }

              break;
            }

          case PRISM18:
            {
              subelem[0] = new Tet10;
              subelem[1] = new Tet10;
              subelem[2] = new Tet10;

              // Split on 0-4 diagonal
              if (split_first_diagonal(elem, 0,4, 1,3))
                {
                  // Split on 0-5 diagonal
                  if (split_first_diagonal(elem, 0,5, 2,3))
                    {
                      // Split on 1-5 diagonal
                      if (split_first_diagonal(elem, 1,5, 2,4))
                        {
                          subelem[0]->set_node(0) = elem->node_ptr(0);
                          subelem[0]->set_node(1) = elem->node_ptr(4);
                          subelem[0]->set_node(2) = elem->node_ptr(5);
                          subelem[0]->set_node(3) = elem->node_ptr(3);

                          subelem[0]->set_node(4) = elem->node_ptr(15);
                          subelem[0]->set_node(5) = elem->node_ptr(13);
                          subelem[0]->set_node(6) = elem->node_ptr(17);
                          subelem[0]->set_node(7) = elem->node_ptr(9);
                          subelem[0]->set_node(8) = elem->node_ptr(12);
                          subelem[0]->set_node(9) = elem->node_ptr(14);

                          subelem[1]->set_node(0) = elem->node_ptr(0);
                          subelem[1]->set_node(1) = elem->node_ptr(4);
                          subelem[1]->set_node(2) = elem->node_ptr(1);
                          subelem[1]->set_node(3) = elem->node_ptr(5);

                          subelem[1]->set_node(4) = elem->node_ptr(15);
                          subelem[1]->set_node(5) = elem->node_ptr(10);
                          subelem[1]->set_node(6) = elem->node_ptr(6);
                          subelem[1]->set_node(7) = elem->node_ptr(17);
                          subelem[1]->set_node(8) = elem->node_ptr(13);
                          subelem[1]->set_node(9) = elem->node_ptr(16);

                          subelem[2]->set_node(0) = elem->node_ptr(0);
                          subelem[2]->set_node(1) = elem->node_ptr(1);
                          subelem[2]->set_node(2) = elem->node_ptr(2);
                          subelem[2]->set_node(3) = elem->node_ptr(5);

                          subelem[2]->set_node(4) = elem->node_ptr(6);
                          subelem[2]->set_node(5) = elem->node_ptr(7);
                          subelem[2]->set_node(6) = elem->node_ptr(8);
                          subelem[2]->set_node(7) = elem->node_ptr(17);
                          subelem[2]->set_node(8) = elem->node_ptr(16);
                          subelem[2]->set_node(9) = elem->node_ptr(11);
                        }
                      else // Split on 2-4 diagonal
                        {
                          libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                          subelem[0]->set_node(0) = elem->node_ptr(0);
                          subelem[0]->set_node(1) = elem->node_ptr(4);
                          subelem[0]->set_node(2) = elem->node_ptr(5);
                          subelem[0]->set_node(3) = elem->node_ptr(3);

                          subelem[0]->set_node(4) = elem->node_ptr(15);
                          subelem[0]->set_node(5) = elem->node_ptr(13);
                          subelem[0]->set_node(6) = elem->node_ptr(17);
                          subelem[0]->set_node(7) = elem->node_ptr(9);
                          subelem[0]->set_node(8) = elem->node_ptr(12);
                          subelem[0]->set_node(9) = elem->node_ptr(14);

                          subelem[1]->set_node(0) = elem->node_ptr(0);
                          subelem[1]->set_node(1) = elem->node_ptr(4);
                          subelem[1]->set_node(2) = elem->node_ptr(2);
                          subelem[1]->set_node(3) = elem->node_ptr(5);

                          subelem[1]->set_node(4) = elem->node_ptr(15);
                          subelem[1]->set_node(5) = elem->node_ptr(16);
                          subelem[1]->set_node(6) = elem->node_ptr(8);
                          subelem[1]->set_node(7) = elem->node_ptr(17);
                          subelem[1]->set_node(8) = elem->node_ptr(13);
                          subelem[1]->set_node(9) = elem->node_ptr(11);

                          subelem[2]->set_node(0) = elem->node_ptr(0);
                          subelem[2]->set_node(1) = elem->node_ptr(1);
                          subelem[2]->set_node(2) = elem->node_ptr(2);
                          subelem[2]->set_node(3) = elem->node_ptr(4);

                          subelem[2]->set_node(4) = elem->node_ptr(6);
                          subelem[2]->set_node(5) = elem->node_ptr(7);
                          subelem[2]->set_node(6) = elem->node_ptr(8);
                          subelem[2]->set_node(7) = elem->node_ptr(15);
                          subelem[2]->set_node(8) = elem->node_ptr(10);
                          subelem[2]->set_node(9) = elem->node_ptr(16);
                        }
                    }
                  else // Split on 2-3 diagonal
                    {
                      libmesh_assert (split_first_diagonal(elem, 2,3, 0,5));

                      // 0-4 and 2-3 split implies 2-4 split
                      libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                      subelem[0]->set_node(0) = elem->node_ptr(0);
                      subelem[0]->set_node(1) = elem->node_ptr(4);
                      subelem[0]->set_node(2) = elem->node_ptr(2);
                      subelem[0]->set_node(3) = elem->node_ptr(3);

                      subelem[0]->set_node(4) = elem->node_ptr(15);
                      subelem[0]->set_node(5) = elem->node_ptr(16);
                      subelem[0]->set_node(6) = elem->node_ptr(8);
                      subelem[0]->set_node(7) = elem->node_ptr(9);
                      subelem[0]->set_node(8) = elem->node_ptr(12);
                      subelem[0]->set_node(9) = elem->node_ptr(17);

                      subelem[1]->set_node(0) = elem->node_ptr(3);
                      subelem[1]->set_node(1) = elem->node_ptr(4);
                      subelem[1]->set_node(2) = elem->node_ptr(2);
                      subelem[1]->set_node(3) = elem->node_ptr(5);

                      subelem[1]->set_node(4) = elem->node_ptr(12);
                      subelem[1]->set_node(5) = elem->node_ptr(16);
                      subelem[1]->set_node(6) = elem->node_ptr(17);
                      subelem[1]->set_node(7) = elem->node_ptr(14);
                      subelem[1]->set_node(8) = elem->node_ptr(13);
                      subelem[1]->set_node(9) = elem->node_ptr(11);

                      subelem[2]->set_node(0) = elem->node_ptr(0);
                      subelem[2]->set_node(1) = elem->node_ptr(1);
                      subelem[2]->set_node(2) = elem->node_ptr(2);
                      subelem[2]->set_node(3) = elem->node_ptr(4);

                      subelem[2]->set_node(4) = elem->node_ptr(6);
                      subelem[2]->set_node(5) = elem->node_ptr(7);
                      subelem[2]->set_node(6) = elem->node_ptr(8);
                      subelem[2]->set_node(7) = elem->node_ptr(15);
                      subelem[2]->set_node(8) = elem->node_ptr(10);
                      subelem[2]->set_node(9) = elem->node_ptr(16);
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

                      subelem[0]->set_node(0) = elem->node_ptr(1);
                      subelem[0]->set_node(1) = elem->node_ptr(3);
                      subelem[0]->set_node(2) = elem->node_ptr(4);
                      subelem[0]->set_node(3) = elem->node_ptr(5);

                      subelem[0]->set_node(4) = elem->node_ptr(15);
                      subelem[0]->set_node(5) = elem->node_ptr(12);
                      subelem[0]->set_node(6) = elem->node_ptr(10);
                      subelem[0]->set_node(7) = elem->node_ptr(16);
                      subelem[0]->set_node(8) = elem->node_ptr(14);
                      subelem[0]->set_node(9) = elem->node_ptr(13);

                      subelem[1]->set_node(0) = elem->node_ptr(1);
                      subelem[1]->set_node(1) = elem->node_ptr(0);
                      subelem[1]->set_node(2) = elem->node_ptr(3);
                      subelem[1]->set_node(3) = elem->node_ptr(5);

                      subelem[1]->set_node(4) = elem->node_ptr(6);
                      subelem[1]->set_node(5) = elem->node_ptr(9);
                      subelem[1]->set_node(6) = elem->node_ptr(15);
                      subelem[1]->set_node(7) = elem->node_ptr(16);
                      subelem[1]->set_node(8) = elem->node_ptr(17);
                      subelem[1]->set_node(9) = elem->node_ptr(14);

                      subelem[2]->set_node(0) = elem->node_ptr(0);
                      subelem[2]->set_node(1) = elem->node_ptr(1);
                      subelem[2]->set_node(2) = elem->node_ptr(2);
                      subelem[2]->set_node(3) = elem->node_ptr(5);

                      subelem[2]->set_node(4) = elem->node_ptr(6);
                      subelem[2]->set_node(5) = elem->node_ptr(7);
                      subelem[2]->set_node(6) = elem->node_ptr(8);
                      subelem[2]->set_node(7) = elem->node_ptr(17);
                      subelem[2]->set_node(8) = elem->node_ptr(16);
                      subelem[2]->set_node(9) = elem->node_ptr(11);
                    }
                  else // Split on 2-3 diagonal
                    {
                      libmesh_assert (split_first_diagonal(elem, 2,3, 0,5));

                      // Split on 1-5 diagonal
                      if (split_first_diagonal(elem, 1,5, 2,4))
                        {
                          subelem[0]->set_node(0) = elem->node_ptr(0);
                          subelem[0]->set_node(1) = elem->node_ptr(1);
                          subelem[0]->set_node(2) = elem->node_ptr(2);
                          subelem[0]->set_node(3) = elem->node_ptr(3);

                          subelem[0]->set_node(4) = elem->node_ptr(6);
                          subelem[0]->set_node(5) = elem->node_ptr(7);
                          subelem[0]->set_node(6) = elem->node_ptr(8);
                          subelem[0]->set_node(7) = elem->node_ptr(9);
                          subelem[0]->set_node(8) = elem->node_ptr(15);
                          subelem[0]->set_node(9) = elem->node_ptr(17);

                          subelem[1]->set_node(0) = elem->node_ptr(3);
                          subelem[1]->set_node(1) = elem->node_ptr(1);
                          subelem[1]->set_node(2) = elem->node_ptr(2);
                          subelem[1]->set_node(3) = elem->node_ptr(5);

                          subelem[1]->set_node(4) = elem->node_ptr(15);
                          subelem[1]->set_node(5) = elem->node_ptr(7);
                          subelem[1]->set_node(6) = elem->node_ptr(17);
                          subelem[1]->set_node(7) = elem->node_ptr(14);
                          subelem[1]->set_node(8) = elem->node_ptr(16);
                          subelem[1]->set_node(9) = elem->node_ptr(11);

                          subelem[2]->set_node(0) = elem->node_ptr(1);
                          subelem[2]->set_node(1) = elem->node_ptr(3);
                          subelem[2]->set_node(2) = elem->node_ptr(4);
                          subelem[2]->set_node(3) = elem->node_ptr(5);

                          subelem[2]->set_node(4) = elem->node_ptr(15);
                          subelem[2]->set_node(5) = elem->node_ptr(12);
                          subelem[2]->set_node(6) = elem->node_ptr(10);
                          subelem[2]->set_node(7) = elem->node_ptr(16);
                          subelem[2]->set_node(8) = elem->node_ptr(14);
                          subelem[2]->set_node(9) = elem->node_ptr(13);
                        }
                      else // Split on 2-4 diagonal
                        {
                          libmesh_assert (split_first_diagonal(elem, 2,4, 1,5));

                          subelem[0]->set_node(0) = elem->node_ptr(0);
                          subelem[0]->set_node(1) = elem->node_ptr(1);
                          subelem[0]->set_node(2) = elem->node_ptr(2);
                          subelem[0]->set_node(3) = elem->node_ptr(3);

                          subelem[0]->set_node(4) = elem->node_ptr(6);
                          subelem[0]->set_node(5) = elem->node_ptr(7);
                          subelem[0]->set_node(6) = elem->node_ptr(8);
                          subelem[0]->set_node(7) = elem->node_ptr(9);
                          subelem[0]->set_node(8) = elem->node_ptr(15);
                          subelem[0]->set_node(9) = elem->node_ptr(17);

                          subelem[1]->set_node(0) = elem->node_ptr(2);
                          subelem[1]->set_node(1) = elem->node_ptr(3);
                          subelem[1]->set_node(2) = elem->node_ptr(4);
                          subelem[1]->set_node(3) = elem->node_ptr(5);

                          subelem[1]->set_node(4) = elem->node_ptr(17);
                          subelem[1]->set_node(5) = elem->node_ptr(12);
                          subelem[1]->set_node(6) = elem->node_ptr(16);
                          subelem[1]->set_node(7) = elem->node_ptr(11);
                          subelem[1]->set_node(8) = elem->node_ptr(14);
                          subelem[1]->set_node(9) = elem->node_ptr(13);

                          subelem[2]->set_node(0) = elem->node_ptr(3);
                          subelem[2]->set_node(1) = elem->node_ptr(1);
                          subelem[2]->set_node(2) = elem->node_ptr(2);
                          subelem[2]->set_node(3) = elem->node_ptr(4);

                          subelem[2]->set_node(4) = elem->node_ptr(15);
                          subelem[2]->set_node(5) = elem->node_ptr(7);
                          subelem[2]->set_node(6) = elem->node_ptr(17);
                          subelem[2]->set_node(7) = elem->node_ptr(12);
                          subelem[2]->set_node(8) = elem->node_ptr(10);
                          subelem[2]->set_node(9) = elem->node_ptr(16);
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
          case TET4:
          case TET10:
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



        // Be sure the correct IDs are also set for all subelems.
        for (unsigned int i=0; i != max_subelems; ++i)
          if (subelem[i]) {
            subelem[i]->processor_id() = elem->processor_id();
            subelem[i]->subdomain_id() = elem->subdomain_id();
          }

        // On a mesh with boundary data, we need to move that data to
        // the new elements.

        // On a mesh which is distributed, we need to move
        // remote_elem links to the new elements.
        bool mesh_is_serial = mesh.is_serial();

        if (mesh_has_boundary_data || mesh_is_serial)
          {
            // Container to key boundary IDs handed back by the BoundaryInfo object.
            std::vector<boundary_id_type> bc_ids;

            for (auto sn : elem->side_index_range())
              {
                mesh.get_boundary_info().boundary_ids(elem, sn, bc_ids);
                for (const auto & b_id : bc_ids)
                  {
                    if (mesh_is_serial && b_id == BoundaryInfo::invalid_id)
                      continue;

                    // Make a sorted list of node ids for elem->side(sn)
                    std::unique_ptr<Elem> elem_side = elem->build_side_ptr(sn);
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
                              std::unique_ptr<Elem> subside_elem = subelem[i]->build_side_ptr(subside);

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
                                  if (b_id != BoundaryInfo::invalid_id)
                                    {
                                      new_bndry_ids.push_back(b_id);
                                      new_bndry_elements.push_back(subelem[i]);
                                      new_bndry_sides.push_back(subside);
                                    }

                                  // If the original element had a RemoteElem neighbor on side 'sn',
                                  // then the subelem has one on side 'subside'.
                                  if (elem->neighbor_ptr(sn) == RemoteElem::get_instance())
                                    subelem[i]->set_neighbor(subside, const_cast<RemoteElem*>(RemoteElem::get_instance()));
                                }
                            }
                        }
                  } // end for loop over boundary IDs
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
              subelem[i]->set_unique_id() = max_unique_id + 2*elem->unique_id() + i;
#endif

              // Prepare to add the newly-created simplices
              new_elements.push_back(subelem[i]);
            }

        // Delete the original element
        mesh.delete_elem(elem);
      } // End for loop over elements
  } // end scope


  // Now, iterate over the new elements vector, and add them each to
  // the Mesh.
  for (auto & elem : new_elements)
    mesh.add_elem(elem);

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
  mesh.prepare_for_use(/*skip_renumber =*/ false);

  // Let the new_elements and new_bndry_elements vectors go out of scope.
}

} // end namespace Meshtools::Modification
} // end namespace MeshTools

} // namespace libMesh


#endif // LIBMESH_MESH_MODIFICATION_H
