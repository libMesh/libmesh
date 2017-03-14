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
#include "libmesh/string_to_enum.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/partitioner.h"

namespace
{
bool split_first_diagonal(const libMesh::Elem * elem,
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

  // First find nodes on the boundary and flag them
  // so that we don't move them
  // on_boundary holds false (not on boundary) and true (on boundary)
  std::vector<bool> on_boundary (mesh.max_node_id(), false);

  if (!perturb_boundary) MeshTools::find_boundary_nodes (mesh, on_boundary);

  // Now calculate the minimum distance to
  // neighboring nodes for each node.
  // hmin holds these distances.
  std::vector<float> hmin (mesh.max_node_id(),
                           std::numeric_limits<float>::max());

  MeshBase::element_iterator       el  = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

  for (; el!=end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      hmin[(*el)->node_id(n)] = std::min(hmin[(*el)->node_id(n)],
                                         static_cast<float>((*el)->hmin()));


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
    for (unsigned int n=0; n<mesh.max_node_id(); n++)
      if (!on_boundary[n] && (hmin[n] < 1.e20) )
        {
          // the direction, random but unit normalized

          Point dir( static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX),
                     (mesh.mesh_dimension() > 1) ?
                     static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX)
                     : 0.,
                     ((mesh.mesh_dimension() == 3) ?
                      static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX)
                      : 0.)
                     );

          dir(0) = (dir(0)-.5)*2.;
          if (mesh.mesh_dimension() > 1)
            dir(1) = (dir(1)-.5)*2.;
          if (mesh.mesh_dimension() == 3)
            dir(2) = (dir(2)-.5)*2.;

          dir = dir.unit();

          Node * node = mesh.query_node_ptr(n);
          if (!node)
            continue;

          (*node)(0) += dir(0)*factor*hmin[n];
          if (mesh.mesh_dimension() > 1)
            (*node)(1) += dir(1)*factor*hmin[n];
          if (mesh.mesh_dimension() == 3)
            (*node)(2) += dir(2)*factor*hmin[n];
        }
  }
}



void MeshTools::Modification::redistribute (MeshBase & mesh,
                                            const FunctionBase<Real> & mapfunc)
{
  libmesh_assert (mesh.n_nodes());
  libmesh_assert (mesh.n_elem());

  LOG_SCOPE("redistribute()", "MeshTools::Modification");

  DenseVector<Real> output_vec(LIBMESH_DIM);

  // FIXME - we should thread this later.
  UniquePtr<FunctionBase<Real> > myfunc = mapfunc.clone();

  MeshBase::node_iterator       it  = mesh.nodes_begin();
  const MeshBase::node_iterator end = mesh.nodes_end();

  for (; it != end; ++it)
    {
      Node * node = *it;

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



void MeshTools::Modification::translate (MeshBase & mesh,
                                         const Real xt,
                                         const Real yt,
                                         const Real zt)
{
  const Point p(xt, yt, zt);

  const MeshBase::node_iterator nd_end = mesh.nodes_end();

  for (MeshBase::node_iterator nd = mesh.nodes_begin();
       nd != nd_end; ++nd)
    **nd += p;
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
  const MeshBase::node_iterator nd_end = mesh.nodes_end();

  for (MeshBase::node_iterator nd = mesh.nodes_begin();
       nd != nd_end; ++nd)
    {
      const Point pt = **nd;
      const Real  x  = pt(0);
      const Real  y  = pt(1);
      const Real  z  = pt(2);
      **nd = Point(( cp*cs-sp*ct*ss)*x + ( sp*cs+cp*ct*ss)*y + (st*ss)*z,
                   (-cp*ss-sp*ct*cs)*x + (-sp*ss+cp*ct*cs)*y + (st*cs)*z,
                   ( sp*st)*x          + (-cp*st)*y          + (ct)*z   );
    }
#else
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
  const MeshBase::node_iterator nd_end = mesh.nodes_end();

  for (MeshBase::node_iterator nd = mesh.nodes_begin();
       nd != nd_end; ++nd)
    (**nd)(0) *= x_scale;


  // Only scale the y coordinate in 2 and 3D
  if (LIBMESH_DIM < 2)
    return;

  for (MeshBase::node_iterator nd = mesh.nodes_begin();
       nd != nd_end; ++nd)
    (**nd)(1) *= y_scale;

  // Only scale the z coordinate in 3D
  if (LIBMESH_DIM < 3)
    return;

  for (MeshBase::node_iterator nd = mesh.nodes_begin();
       nd != nd_end; ++nd)
    (**nd)(2) *= z_scale;
}




// ------------------------------------------------------------
// UnstructuredMesh class member functions for mesh modification
void UnstructuredMesh::all_first_order ()
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  START_LOG("all_first_order()", "Mesh");

  /**
   * Prepare to identify (and then delete) a bunch of no-longer-used nodes.
   */
  std::vector<bool> node_touched_by_me(this->max_node_id(), false);

  /**
   * Loop over the high-ordered elements.
   * First make sure they _are_ indeed high-order, and then replace
   * them with an equivalent first-order element.
   */
  element_iterator endit = elements_end();
  for (element_iterator it = elements_begin();
       it != endit; ++it)
    {
      Elem * so_elem = *it;

      libmesh_assert(so_elem);

      /*
       * build the first-order equivalent, add to
       * the new_elements list.
       */
      Elem * lo_elem = Elem::build
        (Elem::first_order_equivalent_type
         (so_elem->type()), so_elem->parent()).release();

      for (unsigned int s=0; s != so_elem->n_sides(); ++s)
        if (so_elem->neighbor_ptr(s) == remote_elem)
          lo_elem->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));

#ifdef LIBMESH_ENABLE_AMR
      /*
       * Reset the parent links of any child elements
       */
      if (so_elem->has_children())
        for (unsigned int c=0; c != so_elem->n_children(); ++c)
          {
            so_elem->child_ptr(c)->set_parent(lo_elem);
            lo_elem->add_child(so_elem->child_ptr(c), c);
          }

      /*
       * Reset the child link of any parent element
       */
      if (so_elem->parent())
        {
          unsigned int c =
            so_elem->parent()->which_child_am_i(so_elem);
          lo_elem->parent()->replace_child(lo_elem, c);
        }

      /*
       * Copy as much data to the new element as makes sense
       */
      lo_elem->set_p_level(so_elem->p_level());
      lo_elem->set_refinement_flag(so_elem->refinement_flag());
      lo_elem->set_p_refinement_flag(so_elem->p_refinement_flag());
#endif

      libmesh_assert_equal_to (lo_elem->n_vertices(), so_elem->n_vertices());

      /*
       * By definition the vertices of the linear and
       * second order element are identically numbered.
       * transfer these.
       */
      for (unsigned int v=0; v < so_elem->n_vertices(); v++)
        {
          lo_elem->set_node(v) = so_elem->node_ptr(v);
          node_touched_by_me[lo_elem->node_id(v)] = true;
        }

      /*
       * find_neighbors relies on remote_elem neighbor links being
       * properly maintained.
       */
      for (unsigned short s=0; s<so_elem->n_sides(); s++)
        {
          if (so_elem->neighbor(s) == remote_elem)
            lo_elem->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));
        }

      /**
       * If the second order element had any boundary conditions they
       * should be transfered to the first-order element.  The old
       * boundary conditions will be removed from the BoundaryInfo
       * data structure by insert_elem.
       */
      this->get_boundary_info().copy_boundary_ids
        (this->get_boundary_info(), so_elem, lo_elem);

      /*
       * The new first-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the second-order element.
       */
      lo_elem->set_id(so_elem->id());
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      lo_elem->set_unique_id() = so_elem->unique_id();
#endif
      lo_elem->processor_id() = so_elem->processor_id();
      lo_elem->subdomain_id() = so_elem->subdomain_id();
      this->insert_elem(lo_elem);
    }

  const MeshBase::node_iterator nd_end = this->nodes_end();
  MeshBase::node_iterator nd = this->nodes_begin();
  while (nd != nd_end)
    {
      Node * the_node = *nd;
      ++nd;
      if (!node_touched_by_me[the_node->id()])
        this->delete_node(the_node);
    }

  // If crazy people applied boundary info to non-vertices and then
  // deleted those non-vertices, we should make sure their boundary id
  // caches are correct.
  this->get_boundary_info().regenerate_id_sets();

  STOP_LOG("all_first_order()", "Mesh");

  // On hanging nodes that used to also be second order nodes, we
  // might now have an invalid nodal processor_id()
  Partitioner::set_node_processor_ids(*this);

  // delete or renumber nodes, etc
  this->prepare_for_use(/*skip_renumber =*/ false);
}



void UnstructuredMesh::all_second_order (const bool full_ordered)
{
  // This function must be run on all processors at once
  parallel_object_only();

  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  /*
   * If the mesh is empty
   * then we have nothing to do
   */
  if (!this->n_elem())
    return;

  /*
   * If the mesh is already second order
   * then we have nothing to do.
   * We have to test for this in a round-about way to avoid
   * a bug on distributed parallel meshes with more processors
   * than elements.
   */
  bool already_second_order = false;
  if (this->elements_begin() != this->elements_end() &&
      (*(this->elements_begin()))->default_order() != FIRST)
    already_second_order = true;
  this->comm().max(already_second_order);
  if (already_second_order)
    return;

  START_LOG("all_second_order()", "Mesh");

  /*
   * this map helps in identifying second order
   * nodes.  Namely, a second-order node:
   * - edge node
   * - face node
   * - bubble node
   * is uniquely defined through a set of adjacent
   * vertices.  This set of adjacent vertices is
   * used to identify already added higher-order
   * nodes.  We are safe to use node id's since we
   * make sure that these are correctly numbered.
   */
  std::map<std::vector<dof_id_type>, Node *> adj_vertices_to_so_nodes;

  /*
   * for speed-up of the \p add_point() method, we
   * can reserve memory.  Guess the number of additional
   * nodes for different dimensions
   */
  switch (this->mesh_dimension())
    {
    case 1:
      /*
       * in 1D, there can only be order-increase from Edge2
       * to Edge3.  Something like 1/2 of n_nodes() have
       * to be added
       */
      this->reserve_nodes(static_cast<unsigned int>
                          (1.5*static_cast<double>(this->n_nodes())));
      break;

    case 2:
      /*
       * in 2D, either refine from Tri3 to Tri6 (double the nodes)
       * or from Quad4 to Quad8 (again, double) or Quad9 (2.25 that much)
       */
      this->reserve_nodes(static_cast<unsigned int>
                          (2*static_cast<double>(this->n_nodes())));
      break;


    case 3:
      /*
       * in 3D, either refine from Tet4 to Tet10 (factor = 2.5) up to
       * Hex8 to Hex27 (something  > 3).  Since in 3D there _are_ already
       * quite some nodes, and since we do not want to overburden the memory by
       * a too conservative guess, use the lower bound
       */
      this->reserve_nodes(static_cast<unsigned int>
                          (2.5*static_cast<double>(this->n_nodes())));
      break;

    default:
      // Hm?
      libmesh_error_msg("Unknown mesh dimension " << this->mesh_dimension());
    }



  /*
   * form a vector that will hold the node id's of
   * the vertices that are adjacent to the son-th
   * second-order node.  Pull this outside of the
   * loop so that silly compilers don't repeatedly
   * create and destroy the vector.
   */
  std::vector<dof_id_type> adjacent_vertices_ids;

  /**
   * Loop over the low-ordered elements in the _elements vector.
   * First make sure they _are_ indeed low-order, and then replace
   * them with an equivalent second-order element.  Don't
   * forget to delete the low-order element, or else it will leak!
   */
  element_iterator
    it = elements_begin(),
    endit = elements_end();

  for (; it != endit; ++it)
    {
      // the linear-order element
      Elem * lo_elem = *it;

      libmesh_assert(lo_elem);

      // make sure it is linear order
      if (lo_elem->default_order() != FIRST)
        libmesh_error_msg("ERROR: This is not a linear element: type=" << lo_elem->type());

      // this does _not_ work for refined elements
      libmesh_assert_equal_to (lo_elem->level (), 0);

      /*
       * build the second-order equivalent, add to
       * the new_elements list.  Note that this here
       * is the only point where \p full_ordered
       * is necessary.  The remaining code works well
       * for either type of seconrd-order equivalent, e.g.
       * Hex20 or Hex27, as equivalents for Hex8
       */
      Elem * so_elem =
        Elem::build (Elem::second_order_equivalent_type(lo_elem->type(),
                                                        full_ordered) ).release();

      libmesh_assert_equal_to (lo_elem->n_vertices(), so_elem->n_vertices());


      /*
       * By definition the vertices of the linear and
       * second order element are identically numbered.
       * transfer these.
       */
      for (unsigned int v=0; v < lo_elem->n_vertices(); v++)
        so_elem->set_node(v) = lo_elem->node_ptr(v);

      /*
       * Now handle the additional mid-side nodes.  This
       * is simply handled through a map that remembers
       * the already-added nodes.  This map maps the global
       * ids of the vertices (that uniquely define this
       * higher-order node) to the new node.
       * Notation: son = second-order node
       */
      const unsigned int son_begin = so_elem->n_vertices();
      const unsigned int son_end   = so_elem->n_nodes();


      for (unsigned int son=son_begin; son<son_end; son++)
        {
          const unsigned int n_adjacent_vertices =
            so_elem->n_second_order_adjacent_vertices(son);

          adjacent_vertices_ids.resize(n_adjacent_vertices);

          for (unsigned int v=0; v<n_adjacent_vertices; v++)
            adjacent_vertices_ids[v] =
              so_elem->node_id( so_elem->second_order_adjacent_vertex(son,v) );

          /*
           * \p adjacent_vertices_ids is now in order of the current
           * side.  sort it, so that comparisons  with the
           * \p adjacent_vertices_ids created through other elements'
           * sides can match
           */
          std::sort(adjacent_vertices_ids.begin(),
                    adjacent_vertices_ids.end());


          // does this set of vertices already has a mid-node added?
          std::pair<std::map<std::vector<dof_id_type>, Node *>::iterator,
                    std::map<std::vector<dof_id_type>, Node *>::iterator>
            pos = adj_vertices_to_so_nodes.equal_range (adjacent_vertices_ids);

          // no, not added yet
          if (pos.first == pos.second)
            {
              /*
               * for this set of vertices, there is no
               * second_order node yet.  Add it.
               *
               * compute the location of the new node as
               * the average over the adjacent vertices.
               */
              Point new_location = this->point(adjacent_vertices_ids[0]);
              for (unsigned int v=1; v<n_adjacent_vertices; v++)
                new_location += this->point(adjacent_vertices_ids[v]);

              new_location /= static_cast<Real>(n_adjacent_vertices);

              /* Add the new point to the mesh.
               * If we are on a serialized mesh, then we're doing this
               * all in sync, and the node processor_id will be
               * consistent between processors.
               * If we are on a distributed mesh, we can fix
               * inconsistent processor ids later, but only if every
               * processor gives new nodes a *locally* consistent
               * processor id, so we'll give the new node the
               * processor id of an adjacent element for now and then
               * we'll update that later if appropriate.
               */
              Node * so_node = this->add_point
                (new_location, DofObject::invalid_id,
                 lo_elem->processor_id());

              /*
               * insert the new node with its defining vertex
               * set into the map, and relocate pos to this
               * new entry, so that the so_elem can use
               * \p pos for inserting the node
               */
              adj_vertices_to_so_nodes.insert(pos.first,
                                              std::make_pair(adjacent_vertices_ids,
                                                             so_node));

              so_elem->set_node(son) = so_node;
            }
          // yes, already added.
          else
            {
              Node *so_node = pos.first->second;
              libmesh_assert(so_node);

              so_elem->set_node(son) = so_node;

              // We need to ensure that the processor who should own a
              // node *knows* they own the node.
              if (so_node->processor_id() > lo_elem->processor_id())
                so_node->processor_id() = lo_elem->processor_id();
            }
        }

      /*
       * find_neighbors relies on remote_elem neighbor links being
       * properly maintained.
       */
      for (unsigned short s=0; s<lo_elem->n_sides(); s++)
        {
          if (lo_elem->neighbor_ptr(s) == remote_elem)
            so_elem->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));
        }

      /**
       * If the linear element had any boundary conditions they
       * should be transfered to the second-order element.  The old
       * boundary conditions will be removed from the BoundaryInfo
       * data structure by insert_elem.
       *
       * Also, prepare_for_use() will reconstruct most of our neighbor
       * links, but if we have any remote_elem links in a distributed
       * mesh, they need to be preserved.  We do that in the same loop
       * here.
       */
      this->get_boundary_info().copy_boundary_ids
        (this->get_boundary_info(), lo_elem, so_elem);

      /*
       * The new second-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the first-order element.
       */
      so_elem->set_id(lo_elem->id());
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      so_elem->set_unique_id() = lo_elem->unique_id();
#endif
      so_elem->processor_id() = lo_elem->processor_id();
      so_elem->subdomain_id() = lo_elem->subdomain_id();
      this->insert_elem(so_elem);
    }

  // we can clear the map
  adj_vertices_to_so_nodes.clear();


  STOP_LOG("all_second_order()", "Mesh");

  // In a DistributedMesh our ghost node processor ids may be bad,
  // the ids of nodes touching remote elements may be inconsistent,
  // and unique_ids of newly added non-local nodes remain unset.
  // make_nodes_parallel_consistent() will fix all this.
  if (!this->is_serial())
    MeshCommunication().make_nodes_parallel_consistent (*this);

  // renumber nodes, elements etc
  this->prepare_for_use(/*skip_renumber =*/ false);
}






void MeshTools::Modification::all_tri (MeshBase & mesh)
{
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

    MeshBase::element_iterator       el  = mesh.elements_begin();
    const MeshBase::element_iterator end = mesh.elements_end();

    for (; el!=end; ++el)
      {
        Elem * elem = *el;

        const ElemType etype = elem->type();

        // all_tri currently only works on coarse meshes
        libmesh_assert (!elem->parent());

        // The new elements we will split the quad into.
        // In 3D we may need as many as 6 tets per hex
        Elem * subelem[6];

        for (unsigned int i = 0; i != max_subelems; ++i)
          subelem[i] = libmesh_nullptr;

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

            for (unsigned int sn=0; sn<elem->n_sides(); ++sn)
              {
                mesh.get_boundary_info().boundary_ids(*el, sn, bc_ids);
                for (std::vector<boundary_id_type>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
                  {
                    const boundary_id_type b_id = *id_it;

                    if (mesh_is_serial && b_id == BoundaryInfo::invalid_id)
                      continue;

                    // Make a sorted list of node ids for elem->side(sn)
                    UniquePtr<Elem> elem_side = elem->build_side_ptr(sn);
                    std::vector<dof_id_type> elem_side_nodes(elem_side->n_nodes());
                    for (std::size_t esn=0; esn<elem_side_nodes.size(); ++esn)
                      elem_side_nodes[esn] = elem_side->node_id(esn);
                    std::sort(elem_side_nodes.begin(), elem_side_nodes.end());

                    for (unsigned int i=0; i != max_subelems; ++i)
                      if (subelem[i])
                        {
                          for (unsigned int subside=0; subside < subelem[i]->n_sides(); ++subside)
                            {
                              UniquePtr<Elem> subside_elem = subelem[i]->build_side_ptr(subside);

                              // Make a list of *vertex* node ids for this subside, see if they are all present
                              // in elem->side(sn).  Note 1: we can't just compare elem->key(sn) to
                              // subelem[i]->key(subside) in the Prism cases, since the new side is
                              // a different type.  Note 2: we only use vertex nodes since, in the future,
                              // a Hex20 or Prism15's QUAD8 face may be split into two Tri6 faces, and the
                              // original face will not contain the mid-edge node.
                              std::vector<dof_id_type> subside_nodes(subside_elem->n_vertices());
                              for (std::size_t ssn=0; ssn<subside_nodes.size(); ++ssn)
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
                                  if (elem->neighbor_ptr(sn) == remote_elem)
                                    subelem[i]->set_neighbor(subside, const_cast<RemoteElem*>(remote_elem));
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
  {
    std::vector<Elem *>::iterator el        = new_elements.begin();
    const std::vector<Elem *>::iterator end = new_elements.end();
    for (; el != end; ++el)
      mesh.add_elem(*el);
  }

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
      for (std::size_t s=0; s<new_bndry_elements.size(); ++s)
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


void MeshTools::Modification::smooth (MeshBase & mesh,
                                      const unsigned int n_iterations,
                                      const Real power)
{
  /**
   * This implementation assumes every element "side" has only 2 nodes.
   */
  libmesh_assert_equal_to (mesh.mesh_dimension(), 2);

  /*
   * find the boundary nodes
   */
  std::vector<bool>  on_boundary;
  MeshTools::find_boundary_nodes(mesh, on_boundary);

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
            /*
             * Loop over the elements to calculate new node positions
             */
            MeshBase::element_iterator       el  = mesh.level_elements_begin(refinement_level);
            const MeshBase::element_iterator end = mesh.level_elements_end(refinement_level);

            for (; el != end; ++el)
              {
                /*
                 * Constant handle for the element
                 */
                const Elem * elem = *el;

                /*
                 * We relax all nodes on level 0 first
                 * If the element is refined (level > 0), we interpolate the
                 * parents nodes with help of the embedding matrix
                 */
                if (refinement_level == 0)
                  {
                    for (unsigned int s=0; s<elem->n_neighbors(); s++)
                      {
                        /*
                         * Only operate on sides which are on the
                         * boundary or for which the current element's
                         * id is greater than its neighbor's.
                         * Sides get only built once.
                         */
                        if ((elem->neighbor_ptr(s) != libmesh_nullptr) &&
                            (elem->id() > elem->neighbor_ptr(s)->id()))
                          {
                            UniquePtr<const Elem> side(elem->build_side_ptr(s));

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
                    for (unsigned int c=0; c < parent->n_children(); c++)
                      {
                        if (parent->child_ptr(c) == elem)
                          {
                            /*
                             *loop over the childs (that is, the current elements) nodes
                             */
                            for (unsigned int nc=0; nc < elem->n_nodes(); nc++)
                              {
                                /*
                                 * the new position of the node
                                 */
                                Point point;
                                for (unsigned int n=0; n<parent->n_nodes(); n++)
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

                          } // if parent->child == elem
                      } // for parent->n_children
                  } // if element refinement_level
#endif // #ifdef LIBMESH_ENABLE_AMR

              } // element loop

            /*
             * finally reposition the vertex nodes
             */
            for (unsigned int nid=0; nid<mesh.n_nodes(); ++nid)
              if (!on_boundary[nid] && weight[nid] > 0.)
                mesh.node_ref(nid) = new_positions[nid]/weight[nid];
          }

          {
            /*
             * Now handle the additional second_order nodes by calculating
             * their position based on the vertex postitions
             * we do a second loop over the level elements
             */
            MeshBase::element_iterator       el  = mesh.level_elements_begin(refinement_level);
            const MeshBase::element_iterator end = mesh.level_elements_end(refinement_level);

            for (; el != end; ++el)
              {
                /*
                 * Constant handle for the element
                 */
                const Elem * elem = *el;
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
          }

        } // refinement_level loop

    } // end iteration
}



#ifdef LIBMESH_ENABLE_AMR
void MeshTools::Modification::flatten(MeshBase & mesh)
{
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
  {
    MeshBase::element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::element_iterator end = mesh.active_elements_end();

    for (; it != end; ++it)
      {
        Elem * elem = *it;

        // Make a new element of the same type
        Elem * copy = Elem::build(elem->type()).release();

        // Set node pointers (they still point to nodes in the original mesh)
        for(unsigned int n=0; n<elem->n_nodes(); n++)
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
        for (unsigned short s=0; s<elem->n_sides(); s++)
          {
            if (elem->neighbor_ptr(s) == remote_elem)
              copy->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));

            mesh.get_boundary_info().boundary_ids(elem, s, bc_ids);
            for (std::vector<boundary_id_type>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
              {
                const boundary_id_type bc_id = *id_it;

                if (bc_id != BoundaryInfo::invalid_id)
                  {
                    saved_boundary_elements.push_back(copy);
                    saved_bc_ids.push_back(bc_id);
                    saved_bc_sides.push_back(s);
                  }
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
  }


  // Loop again, delete any remaining elements
  {
    MeshBase::element_iterator       it  = mesh.elements_begin();
    const MeshBase::element_iterator end = mesh.elements_end();

    for (; it != end; ++it)
      mesh.delete_elem( *it );
  }


  // Add the copied (now level-0) elements back to the mesh
  {
    for (std::vector<Elem *>::iterator it = new_elements.begin();
         it != new_elements.end();
         ++it)
      {
#ifndef NDEBUG
        dof_id_type orig_id = (*it)->id();

        // ugly mid-statement endif to avoid unused variable warnings
        Elem * added_elem =
#endif
          mesh.add_elem(*it);

#ifndef NDEBUG
        dof_id_type added_id = added_elem->id();
#endif

        // If the Elem, as it was re-added to the mesh, now has a
        // different ID (this is unlikely, so it's just an assert)
        // the boundary information will no longer be correct.
        libmesh_assert_equal_to (orig_id, added_id);
      }
  }

  // Finally, also add back the saved boundary information
  for (std::size_t e=0; e<saved_boundary_elements.size(); ++e)
    mesh.get_boundary_info().add_side(saved_boundary_elements[e],
                                      saved_bc_sides[e],
                                      saved_bc_ids[e]);

  // Trim unused and renumber nodes and elements
  mesh.prepare_for_use(/*skip_renumber =*/ false);
}
#endif // #ifdef LIBMESH_ENABLE_AMR



void MeshTools::Modification::change_boundary_id (MeshBase & mesh,
                                                  const boundary_id_type old_id,
                                                  const boundary_id_type new_id)
{
  if(old_id == new_id)
    {
      // If the IDs are the same, this is a no-op.
      return;
    }

  // A reference to the Mesh's BoundaryInfo object, for convenience.
  BoundaryInfo & bi = mesh.get_boundary_info();

  {
    // Build a list of all nodes that have boundary IDs
    std::vector<dof_id_type> node_list;
    std::vector<boundary_id_type> bc_id_list;
    bi.build_node_list (node_list, bc_id_list);

    // Temporary vector to hold ids
    std::vector<boundary_id_type> bndry_ids;

    // For each node with the old_id...
    for (std::size_t idx=0; idx<node_list.size(); ++idx)
      if (bc_id_list[idx] == old_id)
        {
          // Get the node in question
          const Node * node = mesh.node_ptr(node_list[idx]);

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
    // Build a list of all edges that have boundary IDs
    std::vector<dof_id_type> elem_list;
    std::vector<unsigned short int> edge_list;
    std::vector<boundary_id_type> bc_id_list;
    bi.build_edge_list (elem_list, edge_list, bc_id_list);

    // Temporary vector to hold ids
    std::vector<boundary_id_type> bndry_ids;

    // For each edge with the old_id...
    for (std::size_t idx=0; idx<elem_list.size(); ++idx)
      if (bc_id_list[idx] == old_id)
        {
          // Get the elem in question
          const Elem * elem = mesh.elem_ptr(elem_list[idx]);

          // The edge of the elem in question
          unsigned short int edge = edge_list[idx];

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
    // Build a list of all shell-faces that have boundary IDs
    std::vector<dof_id_type> elem_list;
    std::vector<unsigned short int> shellface_list;
    std::vector<boundary_id_type> bc_id_list;
    bi.build_shellface_list (elem_list, shellface_list, bc_id_list);

    // Temporary vector to hold ids
    std::vector<boundary_id_type> bndry_ids;

    // For each shellface with the old_id...
    for (std::size_t idx=0; idx<elem_list.size(); ++idx)
      if (bc_id_list[idx] == old_id)
        {
          // Get the elem in question
          const Elem * elem = mesh.elem_ptr(elem_list[idx]);

          // The shellface of the elem in question
          unsigned short int shellface = shellface_list[idx];

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
    // Build a list of all sides that have boundary IDs
    std::vector<dof_id_type> elem_list;
    std::vector<unsigned short int> side_list;
    std::vector<boundary_id_type> bc_id_list;
    bi.build_side_list (elem_list, side_list, bc_id_list);

    // Temporary vector to hold ids
    std::vector<boundary_id_type> bndry_ids;

    // For each side with the old_id...
    for (std::size_t idx=0; idx<elem_list.size(); ++idx)
      if (bc_id_list[idx] == old_id)
        {
          // Get the elem in question
          const Elem * elem = mesh.elem_ptr(elem_list[idx]);

          // The side of the elem in question
          unsigned short int side = side_list[idx];

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
  if(old_id == new_id)
    {
      // If the IDs are the same, this is a no-op.
      return;
    }

  MeshBase::element_iterator           el = mesh.elements_begin();
  const MeshBase::element_iterator end_el = mesh.elements_end();

  for (; el != end_el; ++el)
    {
      Elem * elem = *el;

      if (elem->subdomain_id() == old_id)
        elem->subdomain_id() = new_id;
    }
}


} // namespace libMesh
