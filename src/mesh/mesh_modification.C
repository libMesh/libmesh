// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/location_maps.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel.h"
#include "libmesh/remote_elem.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/unstructured_mesh.h"

namespace libMesh
{


// ------------------------------------------------------------
// MeshTools::Modification functions for mesh modification
void MeshTools::Modification::distort (MeshBase& mesh,
				       const Real factor,
				       const bool perturb_boundary)
{
  libmesh_assert (mesh.n_nodes());
  libmesh_assert (mesh.n_elem());
  libmesh_assert ((factor >= 0.) && (factor <= 1.));

  START_LOG("distort()", "MeshTools::Modification");



  // First find nodes on the boundary and flag them
  // so that we don't move them
  // on_boundary holds false (not on boundary) and true (on boundary)
  std::vector<bool> on_boundary (mesh.n_nodes(), false);

  if (!perturb_boundary) MeshTools::find_boundary_nodes (mesh, on_boundary);

  // Now calculate the minimum distance to
  // neighboring nodes for each node.
  // hmin holds these distances.
  std::vector<float> hmin (mesh.n_nodes(),
			   std::numeric_limits<float>::max());

  MeshBase::element_iterator       el  = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

  for (; el!=end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      hmin[(*el)->node(n)] = std::min(hmin[(*el)->node(n)],
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
    for (unsigned int n=0; n<mesh.n_nodes(); n++)
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

          Node *node = mesh.node_ptr(n);
          if (!node)
            continue;

          (*node)(0) += dir(0)*factor*hmin[n];
	  if (mesh.mesh_dimension() > 1)
            (*node)(1) += dir(1)*factor*hmin[n];
          if (mesh.mesh_dimension() == 3)
            (*node)(2) += dir(2)*factor*hmin[n];
	}
  }


  // All done
  STOP_LOG("distort()", "MeshTools::Modification");
}



void MeshTools::Modification::translate (MeshBase& mesh,
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


// void MeshTools::Modification::rotate2D (MeshBase& mesh,
//                                         const Real alpha)
// {
//   libmesh_assert_not_equal_to (mesh.mesh_dimension(), 1);

//   const Real pi = std::acos(-1);
//   const Real  a = alpha/180.*pi;
//   for (unsigned int n=0; n<mesh.n_nodes(); n++)
//     {
//       const Point p = mesh.node(n);
//       const Real  x = p(0);
//       const Real  y = p(1);
//       const Real  z = p(2);
//       mesh.node(n) = Point(std::cos(a)*x - std::sin(a)*y,
//                            std::sin(a)*x + std::cos(a)*y,
//                            z);
//     }

// }



void MeshTools::Modification::rotate (MeshBase& mesh,
				      const Real phi,
				      const Real theta,
				      const Real psi)
{
  libmesh_assert_not_equal_to (mesh.mesh_dimension(), 1);

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
}


void MeshTools::Modification::scale (MeshBase& mesh,
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
  if (mesh.spatial_dimension() < 2)
    return;

  for (MeshBase::node_iterator nd = mesh.nodes_begin();
       nd != nd_end; ++nd)
    (**nd)(1) *= y_scale;

  // Only scale the z coordinate in 3D
  if (mesh.spatial_dimension() < 3)
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
      Elem* so_elem = *it;

      libmesh_assert(so_elem);

      /*
       * build the first-order equivalent, add to
       * the new_elements list.
       */
      Elem* lo_elem = Elem::build
        (Elem::first_order_equivalent_type
          (so_elem->type()), so_elem->parent()).release();

      for (unsigned int s=0; s != so_elem->n_sides(); ++s)
        if (so_elem->neighbor(s) == remote_elem)
          lo_elem->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));

#ifdef LIBMESH_ENABLE_AMR
      /*
       * Reset the parent links of any child elements
       */
      if (so_elem->has_children())
        for (unsigned int c=0; c != so_elem->n_children(); ++c)
          {
            so_elem->child(c)->set_parent(lo_elem);
            lo_elem->add_child(so_elem->child(c), c);
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
	  lo_elem->set_node(v) = so_elem->get_node(v);
          node_touched_by_me[lo_elem->node(v)] = true;
        }

      /**
       * If the second order element had any boundary conditions they
       * should be transfered to the first-order element.  The old
       * boundary conditions will be removed from the BoundaryInfo
       * data structure by insert_elem.
       */
      libmesh_assert_equal_to (lo_elem->n_sides(), so_elem->n_sides());

      for (unsigned int s=0; s<so_elem->n_sides(); s++)
	{
	  const std::vector<boundary_id_type> boundary_ids =
	    this->boundary_info->raw_boundary_ids (so_elem, s);

	  this->boundary_info->add_side (lo_elem, s, boundary_ids);
	}

      /*
       * The new first-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the second-order element.
       */
      lo_elem->set_id(so_elem->id());
      lo_elem->processor_id() = so_elem->processor_id();
      lo_elem->subdomain_id() = so_elem->subdomain_id();
      this->insert_elem(lo_elem);
    }

  const MeshBase::node_iterator nd_end = this->nodes_end();
  MeshBase::node_iterator nd = this->nodes_begin();
    while (nd != nd_end)
    {
      Node *the_node = *nd;
      ++nd;
      if (!node_touched_by_me[the_node->id()])
        this->delete_node(the_node);
    }

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
  std::map<std::vector<dof_id_type>, Node*> adj_vertices_to_so_nodes;

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
      libmesh_error();
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
  const_element_iterator endit = elements_end();
  for (const_element_iterator it = elements_begin();
       it != endit; ++it)
    {
      // the linear-order element
      const Elem* lo_elem = *it;

      libmesh_assert(lo_elem);

      // make sure it is linear order
      if (lo_elem->default_order() != FIRST)
        {
	  libMesh::err << "ERROR: This is not a linear element: type="
		        << lo_elem->type() << std::endl;
	  libmesh_error();
	}

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
      Elem* so_elem =
	Elem::build (Elem::second_order_equivalent_type(lo_elem->type(),
							full_ordered) ).release();

      libmesh_assert_equal_to (lo_elem->n_vertices(), so_elem->n_vertices());


      /*
       * By definition the vertices of the linear and
       * second order element are identically numbered.
       * transfer these.
       */
      for (unsigned int v=0; v < lo_elem->n_vertices(); v++)
	so_elem->set_node(v) = lo_elem->get_node(v);

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
	      so_elem->node( so_elem->second_order_adjacent_vertex(son,v) );

	  /*
	   * \p adjacent_vertices_ids is now in order of the current
	   * side.  sort it, so that comparisons  with the
	   * \p adjacent_vertices_ids created through other elements'
	   * sides can match
	   */
	  std::sort(adjacent_vertices_ids.begin(),
		    adjacent_vertices_ids.end());


	  // does this set of vertices already has a mid-node added?
	  std::pair<std::map<std::vector<dof_id_type>, Node*>::iterator,
                    std::map<std::vector<dof_id_type>, Node*>::iterator>
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

	      /* Add the new point to the mesh, giving it a globally
               * well-defined processor id.
	       */
	      Node* so_node = this->add_point
                (new_location, DofObject::invalid_id,
                this->node(adjacent_vertices_ids[0]).processor_id());

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
	      libmesh_assert(pos.first->second);

	      so_elem->set_node(son) = pos.first->second;
	    }
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
      libmesh_assert_equal_to (lo_elem->n_sides(), so_elem->n_sides());

      for (unsigned int s=0; s<lo_elem->n_sides(); s++)
	{
	  const std::vector<boundary_id_type> boundary_ids =
	    this->boundary_info->raw_boundary_ids (lo_elem, s);

	  this->boundary_info->add_side (so_elem, s, boundary_ids);

	  if (lo_elem->neighbor(s) == remote_elem)
            so_elem->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));
	}

      /*
       * The new second-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the first-order element.
       */
      so_elem->set_id(lo_elem->id());
      so_elem->processor_id() = lo_elem->processor_id();
      so_elem->subdomain_id() = lo_elem->subdomain_id();
      this->insert_elem(so_elem);
    }

  // we can clear the map
  adj_vertices_to_so_nodes.clear();


  STOP_LOG("all_second_order()", "Mesh");

  // In a ParallelMesh our ghost node processor ids may be bad and
  // the ids of nodes touching remote elements may be inconsistent.
  // Fix them.
  if (!this->is_serial())
    {
      LocationMap<Node> loc_map;
      MeshCommunication().make_nodes_parallel_consistent
        (*this, loc_map);
    }

  // renumber nodes, elements etc
  this->prepare_for_use(/*skip_renumber =*/ false);
}






void MeshTools::Modification::all_tri (MeshBase& mesh)
{
  // The number of elements in the original mesh before any additions
  // or deletions.
  const dof_id_type n_orig_elem = mesh.n_elem();

  // We store pointers to the newly created elements in a vector
  // until they are ready to be added to the mesh.  This is because
  // adding new elements on the fly can cause reallocation and invalidation
  // of existing iterators.
  std::vector<Elem*> new_elements;
  new_elements.reserve (2*n_orig_elem);

  // If the original mesh has boundary data, we carry that over
  // to the new mesh with triangular elements.
  const bool mesh_has_boundary_data = (mesh.boundary_info->n_boundary_ids() > 0);

  // Temporary vectors to store the new boundary element pointers, side numbers, and boundary ids
  std::vector<Elem*> new_bndry_elements;
  std::vector<unsigned short int> new_bndry_sides;
  std::vector<boundary_id_type> new_bndry_ids;

  // Iterate over the elements, splitting QUADS into
  // pairs of conforming triangles.
  // FIXME: This algorithm does not work on refined grids!
  {
    MeshBase::element_iterator       el  = mesh.elements_begin();
    const MeshBase::element_iterator end = mesh.elements_end();

    for (; el!=end; ++el)
      {
	Elem* elem = *el;

	const ElemType etype = elem->type();

	// all_tri currently only works on coarse meshes
	libmesh_assert (!elem->parent());

	// We split the quads using the shorter of the two diagonals
	// to maintain the best angle properties.
	bool edge_swap = false;

	// True if we actually split the current element.
	bool split_elem = false;

	// The two new triangular elements we will split the quad into.
	Elem* tri0 = NULL;
	Elem* tri1 = NULL;


	switch (etype)
	  {
	  case QUAD4:
	    {
	      split_elem = true;

	      tri0 = new Tri3;
	      tri1 = new Tri3;

	      // Check for possible edge swap
	      if ((elem->point(0) - elem->point(2)).size() <
		  (elem->point(1) - elem->point(3)).size())
		{
		  tri0->set_node(0) = elem->get_node(0);
		  tri0->set_node(1) = elem->get_node(1);
		  tri0->set_node(2) = elem->get_node(2);

		  tri1->set_node(0) = elem->get_node(0);
		  tri1->set_node(1) = elem->get_node(2);
		  tri1->set_node(2) = elem->get_node(3);
		}

	      else
		{
		  edge_swap=true;

		  tri0->set_node(0) = elem->get_node(0);
		  tri0->set_node(1) = elem->get_node(1);
		  tri0->set_node(2) = elem->get_node(3);

		  tri1->set_node(0) = elem->get_node(1);
		  tri1->set_node(1) = elem->get_node(2);
		  tri1->set_node(2) = elem->get_node(3);
		}


	      break;
	    }

	  case QUAD8:
	    {
	      split_elem =  true;

	      tri0 = new Tri6;
	      tri1 = new Tri6;

	      Node* new_node = mesh.add_point( (mesh.node(elem->node(0)) +
						mesh.node(elem->node(1)) +
						mesh.node(elem->node(2)) +
						mesh.node(elem->node(3)) / 4)
					       );

	      // Check for possible edge swap
	      if ((elem->point(0) - elem->point(2)).size() <
		  (elem->point(1) - elem->point(3)).size())
		{
		  tri0->set_node(0) = elem->get_node(0);
		  tri0->set_node(1) = elem->get_node(1);
		  tri0->set_node(2) = elem->get_node(2);
		  tri0->set_node(3) = elem->get_node(4);
		  tri0->set_node(4) = elem->get_node(5);
		  tri0->set_node(5) = new_node;

		  tri1->set_node(0) = elem->get_node(0);
		  tri1->set_node(1) = elem->get_node(2);
		  tri1->set_node(2) = elem->get_node(3);
		  tri1->set_node(3) = new_node;
		  tri1->set_node(4) = elem->get_node(6);
		  tri1->set_node(5) = elem->get_node(7);

		}

	      else
		{
		  edge_swap=true;

		  tri0->set_node(0) = elem->get_node(3);
		  tri0->set_node(1) = elem->get_node(0);
		  tri0->set_node(2) = elem->get_node(1);
		  tri0->set_node(3) = elem->get_node(7);
		  tri0->set_node(4) = elem->get_node(4);
		  tri0->set_node(5) = new_node;

		  tri1->set_node(0) = elem->get_node(1);
		  tri1->set_node(1) = elem->get_node(2);
		  tri1->set_node(2) = elem->get_node(3);
		  tri1->set_node(3) = elem->get_node(5);
		  tri1->set_node(4) = elem->get_node(6);
		  tri1->set_node(5) = new_node;
		}

	      break;
	    }

	  case QUAD9:
	    {
	      split_elem =  true;

	      tri0 = new Tri6;
	      tri1 = new Tri6;

	      // Check for possible edge swap
	      if ((elem->point(0) - elem->point(2)).size() <
		  (elem->point(1) - elem->point(3)).size())
		{
		  tri0->set_node(0) = elem->get_node(0);
		  tri0->set_node(1) = elem->get_node(1);
		  tri0->set_node(2) = elem->get_node(2);
		  tri0->set_node(3) = elem->get_node(4);
		  tri0->set_node(4) = elem->get_node(5);
		  tri0->set_node(5) = elem->get_node(8);

		  tri1->set_node(0) = elem->get_node(0);
		  tri1->set_node(1) = elem->get_node(2);
		  tri1->set_node(2) = elem->get_node(3);
		  tri1->set_node(3) = elem->get_node(8);
		  tri1->set_node(4) = elem->get_node(6);
		  tri1->set_node(5) = elem->get_node(7);
		}

	      else
		{
		  edge_swap=true;

		  tri0->set_node(0) = elem->get_node(0);
		  tri0->set_node(1) = elem->get_node(1);
		  tri0->set_node(2) = elem->get_node(3);
		  tri0->set_node(3) = elem->get_node(4);
		  tri0->set_node(4) = elem->get_node(8);
		  tri0->set_node(5) = elem->get_node(7);

		  tri1->set_node(0) = elem->get_node(1);
		  tri1->set_node(1) = elem->get_node(2);
		  tri1->set_node(2) = elem->get_node(3);
		  tri1->set_node(3) = elem->get_node(5);
		  tri1->set_node(4) = elem->get_node(6);
		  tri1->set_node(5) = elem->get_node(8);
		}

	      break;
	    }
          // No need to split elements that are already triangles
          case TRI3:
          case TRI6:
            continue;
          // Try to ignore non-2D elements for now
	  default:
	    {
	      libMesh::err << "Warning, encountered non-2D element "
                            << Utility::enum_to_string<ElemType>(etype)
			    << " in MeshTools::Modification::all_tri(), hope that's OK..."
			    << std::endl;
	    }
	  } // end switch (etype)



	if (split_elem)
	  {
	    // Be sure the correct ID's are also set for tri0 and
	    // tri1.
            tri0->processor_id() = elem->processor_id();
            tri0->subdomain_id() = elem->subdomain_id();
            tri1->processor_id() = elem->processor_id();
            tri1->subdomain_id() = elem->subdomain_id();

	    if (mesh_has_boundary_data)
	      {
		for (unsigned int sn=0; sn<elem->n_sides(); ++sn)
		  {
                    const std::vector<boundary_id_type>& bc_ids = mesh.boundary_info->boundary_ids(*el, sn);
                    for (std::vector<boundary_id_type>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
                      {
                        const boundary_id_type b_id = *id_it;

		        if (b_id != BoundaryInfo::invalid_id)
		          {
			    // Add the boundary ID to the list of new boundary ids
			    new_bndry_ids.push_back(b_id);

			    // Convert the boundary side information of the old element to
			    // boundary side information for the new element.
			    if (!edge_swap)
			      {
			        switch (sn)
			          {
			          case 0:
				    {
				      // New boundary side is Tri 0, side 0
				      new_bndry_elements.push_back(tri0);
				      new_bndry_sides.push_back(0);
				      break;
				    }
			          case 1:
				    {
				      // New boundary side is Tri 0, side 1
				      new_bndry_elements.push_back(tri0);
				      new_bndry_sides.push_back(1);
				      break;
				    }
			          case 2:
				    {
				      // New boundary side is Tri 1, side 1
				      new_bndry_elements.push_back(tri1);
				      new_bndry_sides.push_back(1);
				      break;
				    }
			          case 3:
				    {
				      // New boundary side is Tri 1, side 2
				      new_bndry_elements.push_back(tri1);
				      new_bndry_sides.push_back(2);
				      break;
				    }

			          default:
				    {
				      libMesh::err << "Quad4/8/9 cannot have more than 4 sides." << std::endl;
				      libmesh_error();
				    }
			          }
			      }

			    else // edge_swap==true
			      {
			        switch (sn)
			          {
			          case 0:
				    {
				      // New boundary side is Tri 0, side 0
				      new_bndry_elements.push_back(tri0);
				      new_bndry_sides.push_back(0);
				      break;
				    }
			          case 1:
				    {
				      // New boundary side is Tri 1, side 0
				      new_bndry_elements.push_back(tri1);
				      new_bndry_sides.push_back(0);
				      break;
				    }
			          case 2:
				    {
				      // New boundary side is Tri 1, side 1
				      new_bndry_elements.push_back(tri1);
				      new_bndry_sides.push_back(1);
				      break;
				    }
			          case 3:
				    {
				      // New boundary side is Tri 0, side 2
				      new_bndry_elements.push_back(tri0);
				      new_bndry_sides.push_back(2);
				      break;
				    }

			          default:
				    {
				      libMesh::err << "Quad4/8/9 cannot have more than 4 sides." << std::endl;
				      libmesh_error();
				    }
			          }
			      } // end edge_swap==true
		          } // end if (b_id != BoundaryInfo::invalid_id)
                      } // end for loop over boundary IDs
		  } // end for loop over sides

		// Remove the original element from the BoundaryInfo structure.
		mesh.boundary_info->remove(elem);

	      } // end if (mesh_has_boundary_data)


	    // On a distributed mesh, we need to preserve remote_elem
	    // links, since prepare_for_use can't reconstruct them for
	    // us.
            for (unsigned int sn=0; sn<elem->n_sides(); ++sn)
              {
                if (elem->neighbor(sn) == remote_elem)
                  {
                    // Create a remote_elem link on one of the new
                    // elements corresponding to the link from the old
                    // element.
                    if (!edge_swap)
                      {
                        switch (sn)
                              {
                              case 0:
                                {
                                  // New remote side is Tri 0, side 0
                                  tri0->set_neighbor(0, const_cast<RemoteElem*>(remote_elem));
                                  break;
                                }
                              case 1:
                                {
                                  // New remote side is Tri 0, side 1
                                  tri0->set_neighbor(1, const_cast<RemoteElem*>(remote_elem));
                                  break;
                                }
                              case 2:
                                {
                                  // New remote side is Tri 1, side 1
                                  tri1->set_neighbor(1, const_cast<RemoteElem*>(remote_elem));
                                  break;
                                }
                              case 3:
                                {
                                  // New remote side is Tri 1, side 2
                                  tri1->set_neighbor(2, const_cast<RemoteElem*>(remote_elem));
                                  break;
                                }

                              default:
                                {
                                  libMesh::err << "Quad4/8/9 cannot have more than 4 sides." << std::endl;
                                  libmesh_error();
                                }
                              }
                          }

                        else // edge_swap==true
                          {
                            switch (sn)
                              {
                              case 0:
                                {
                                  // New remote side is Tri 0, side 0
                                  tri0->set_neighbor(0, const_cast<RemoteElem*>(remote_elem));
                                  break;
                                }
                              case 1:
                                {
                                  // New remote side is Tri 1, side 0
                                  tri1->set_neighbor(0, const_cast<RemoteElem*>(remote_elem));
                                  break;
                                }
                              case 2:
                                {
                                  // New remote side is Tri 1, side 1
                                  tri1->set_neighbor(1, const_cast<RemoteElem*>(remote_elem));
                                  break;
                                }
                              case 3:
                                {
                                  // New remote side is Tri 0, side 2
                                  tri0->set_neighbor(2, const_cast<RemoteElem*>(remote_elem));
                                  break;
                                }

                              default:
                                {
                                  libMesh::err << "Quad4/8/9 cannot have more than 4 sides." << std::endl;
                                  libmesh_error();
                                }
                              }
                          } // end edge_swap==true
                      } // end if (elem->neighbor(sn) == remote_elem)
              } // end for loop over sides

	    // Determine new IDs for the split elements which will be
	    // the same on all processors, therefore keeping the Mesh
	    // in sync.  Note: we offset the new IDs by n_orig_elem to
	    // avoid overwriting any of the original IDs, this assumes
	    // they were contiguously-numbered to begin with...
	    tri0->set_id( n_orig_elem + 2*elem->id() + 0 );
	    tri1->set_id( n_orig_elem + 2*elem->id() + 1 );

	    // Add the newly-created triangles to the temporary vector of new elements.
	    new_elements.push_back(tri0);
	    new_elements.push_back(tri1);

	    // Delete the original element
	    mesh.delete_elem(elem);
	  } // end if (split_elem)
      } // End for loop over elements
  } // end scope


  // Now, iterate over the new elements vector, and add them each to
  // the Mesh.
  {
    std::vector<Elem*>::iterator el        = new_elements.begin();
    const std::vector<Elem*>::iterator end = new_elements.end();
    for (; el != end; ++el)
      mesh.add_elem(*el);
  }

  if (mesh_has_boundary_data)
    {
      // By this time, we should have removed all of the original boundary sides
      // - except on a hybrid mesh, where we can't "start from a blank slate"! - RHS
      // libmesh_assert_equal_to (mesh.boundary_info->n_boundary_conds(), 0);

      // Clear the boundary info, to be sure and start from a blank slate.
      // mesh.boundary_info->clear();

      // If the old mesh had boundary data, the new mesh better have some.
      libmesh_assert_greater (new_bndry_elements.size(), 0);

      // We should also be sure that the lengths of the new boundary data vectors
      // are all the same.
      libmesh_assert_equal_to (new_bndry_elements.size(), new_bndry_sides.size());
      libmesh_assert_equal_to (new_bndry_sides.size(), new_bndry_ids.size());

      // Add the new boundary info to the mesh
      for (unsigned int s=0; s<new_bndry_elements.size(); ++s)
	mesh.boundary_info->add_side(new_bndry_elements[s],
				     new_bndry_sides[s],
				     new_bndry_ids[s]);
    }


  // Prepare the newly created mesh for use.
  mesh.prepare_for_use(/*skip_renumber =*/ false);

  // Let the new_elements and new_bndry_elements vectors go out of scope.
}


void MeshTools::Modification::smooth (MeshBase& mesh,
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
                const Elem* elem = *el;

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
                        if ((elem->neighbor(s) != NULL) &&
                            (elem->id() > elem->neighbor(s)->id()) )
                          {
                            AutoPtr<Elem> side(elem->build_side(s));

                            Node* node0 = side->get_node(0);
                            Node* node1 = side->get_node(1);

                            Real node_weight = 1.;
                            // calculate the weight of the nodes
                            if (power > 0)
                              {
                                Point diff = (*node0)-(*node1);
                                node_weight = std::pow( diff.size(), power );
                              }

                            const dof_id_type id0 = node0->id(), id1 = node1->id();
                            new_positions[id0].add_scaled( *node1, node_weight );
                            new_positions[id1].add_scaled( *node0, node_weight );
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

                    const Elem* parent = elem->parent();

                    /*
                     * find out which child I am
                     */
                    for (unsigned int c=0; c < parent->n_children(); c++)
                      {
                        if (parent->child(c) == elem)
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

                                const dof_id_type id = elem->get_node(nc)->id();
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
                mesh.node(nid) = new_positions[nid]/weight[nid];
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
                const Elem* elem = *el;
                const unsigned int son_begin = elem->n_vertices();
                const unsigned int son_end   = elem->n_nodes();
                for (unsigned int n=son_begin; n<son_end; n++)
                  {
                    const unsigned int n_adjacent_vertices =
                      elem->n_second_order_adjacent_vertices(n);

                    Point point;
                    for (unsigned int v=0; v<n_adjacent_vertices; v++)
                      point.add(elem->point( elem->second_order_adjacent_vertex(n,v) ));

                    const dof_id_type id = elem->get_node(n)->id();
                    mesh.node(id) = point/n_adjacent_vertices;
                  }
              }
          }

        } // refinement_level loop

    } // end iteration
}



#ifdef LIBMESH_ENABLE_AMR
void MeshTools::Modification::flatten(MeshBase& mesh)
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
  std::vector<Elem*> new_elements;

  // BoundaryInfo Storage for element ids, sides, and BC ids
  std::vector<Elem*>              saved_boundary_elements;
  std::vector<boundary_id_type>   saved_bc_ids;
  std::vector<unsigned short int> saved_bc_sides;

  // Reserve a reasonable amt. of space for each
  new_elements.reserve(mesh.n_active_elem());
  saved_boundary_elements.reserve(mesh.boundary_info->n_boundary_conds());
  saved_bc_ids.reserve(mesh.boundary_info->n_boundary_conds());
  saved_bc_sides.reserve(mesh.boundary_info->n_boundary_conds());
  {
    MeshBase::element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::element_iterator end = mesh.active_elements_end();

    for (; it != end; ++it)
      {
	Elem* elem = *it;

	// Make a new element of the same type
	Elem* copy = Elem::build(elem->type()).release();

	// Set node pointers (they still point to nodes in the original mesh)
	for(unsigned int n=0; n<elem->n_nodes(); n++)
	  copy->set_node(n) = elem->get_node(n);

	// Copy over ids
        copy->processor_id() = elem->processor_id();
        copy->subdomain_id() = elem->subdomain_id();

	// Retain the original element's ID as well, otherwise ParallelMesh will
	// try to create one for you...
	copy->set_id( elem->id() );

	// This element could have boundary info or ParallelMesh
	// remote_elem links as well.  We need to save the (elem,
	// side, bc_id) triples and those links
	for (unsigned int s=0; s<elem->n_sides(); s++)
	  {
            if (elem->neighbor(s) == remote_elem)
              copy->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));

	    const std::vector<boundary_id_type>& bc_ids = mesh.boundary_info->boundary_ids(elem,s);
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
    for (std::vector<Elem*>::iterator it = new_elements.begin();
	 it != new_elements.end();
	 ++it)
      {
	dof_id_type orig_id = (*it)->id();

	Elem* added_elem = mesh.add_elem(*it);

	dof_id_type added_id = added_elem->id();

	// If the Elem, as it was re-added to the mesh, now has a
	// different ID (this is unlikely, so it's just an assert)
	// the boundary information will no longer be correct.
	libmesh_assert_equal_to (orig_id, added_id);
      }
  }

  // Finally, also add back the saved boundary information
  for (unsigned int e=0; e<saved_boundary_elements.size(); ++e)
    mesh.boundary_info->add_side(saved_boundary_elements[e],
				 saved_bc_sides[e],
				 saved_bc_ids[e]);

  // Trim unused and renumber nodes and elements
  mesh.prepare_for_use(/*skip_renumber =*/ false);
}
#endif // #ifdef LIBMESH_ENABLE_AMR



void MeshTools::Modification::change_boundary_id (MeshBase& mesh,
                                                  const boundary_id_type old_id,
                                                  const boundary_id_type new_id)
{
  // Only level-0 elements store BCs.  Loop over them.
  MeshBase::element_iterator           el = mesh.level_elements_begin(0);
  const MeshBase::element_iterator end_el = mesh.level_elements_end(0);
  for (; el != end_el; ++el)
    {
      Elem *elem = *el;

      unsigned int n_nodes = elem->n_nodes();
      for (unsigned int n=0; n != n_nodes; ++n)
        {
          const std::vector<boundary_id_type>& old_ids = mesh.boundary_info->boundary_ids(elem->get_node(n));
          if (std::find(old_ids.begin(), old_ids.end(), old_id) != old_ids.end())
            {
              std::vector<boundary_id_type> new_ids(old_ids);
              std::replace(new_ids.begin(), new_ids.end(), old_id, new_id);
              mesh.boundary_info->remove(elem->get_node(n));
              mesh.boundary_info->add_node(elem->get_node(n), new_ids);
            }
        }

      unsigned int n_edges = elem->n_edges();
      for (unsigned int edge=0; edge != n_edges; ++edge)
        {
          const std::vector<boundary_id_type>& old_ids = mesh.boundary_info->edge_boundary_ids(elem, edge);
          if (std::find(old_ids.begin(), old_ids.end(), old_id) != old_ids.end())
            {
              std::vector<boundary_id_type> new_ids(old_ids);
              std::replace(new_ids.begin(), new_ids.end(), old_id, new_id);
              mesh.boundary_info->remove_edge(elem, edge);
              mesh.boundary_info->add_edge(elem, edge, new_ids);
            }
        }

      unsigned int n_sides = elem->n_sides();
      for (unsigned int s=0; s != n_sides; ++s)
        {
          const std::vector<boundary_id_type>& old_ids = mesh.boundary_info->boundary_ids(elem, s);
          if (std::find(old_ids.begin(), old_ids.end(), old_id) != old_ids.end())
            {
              std::vector<boundary_id_type> new_ids(old_ids);
              std::replace(new_ids.begin(), new_ids.end(), old_id, new_id);
              mesh.boundary_info->remove_side(elem, s);
              mesh.boundary_info->add_side(elem, s, new_ids);
            }
        }
    }
}



void MeshTools::Modification::change_subdomain_id (MeshBase& mesh,
						   const subdomain_id_type old_id,
						   const subdomain_id_type new_id)
{
  MeshBase::element_iterator           el = mesh.elements_begin();
  const MeshBase::element_iterator end_el = mesh.elements_end();

  for (; el != end_el; ++el)
    {
      Elem *elem = *el;

      if (elem->subdomain_id() == old_id)
	elem->subdomain_id() = new_id;
    }
}


} // namespace libMesh
