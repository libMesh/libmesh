// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include <cmath> // for std::acos()
#include <algorithm>
#include <map>

// Local includes
#include "mesh_tools.h"
#include "mesh_modification.h"
#include "unstructured_mesh.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "libmesh_logging.h"
#include "boundary_info.h"
#include "string_to_enum.h"


// ------------------------------------------------------------
// MeshTools::Modification functions for mesh modification
void MeshTools::Modification::distort (MeshBase& mesh,
				       const Real factor,
				       const bool perturb_boundary)
{
  assert (mesh.mesh_dimension() != 1);
  assert (mesh.n_nodes());
  assert (mesh.n_elem());
  assert ((factor >= 0.) && (factor <= 1.));

  START_LOG("distort()", "MeshTools::Modification");



  // First find nodes on the boundary and flag them
  // so that we don't move them
  // on_boundary holds false (not on boundary) and true (on boundary)
  std::vector<bool> on_boundary (mesh.n_nodes(), false);
  
  if (!perturb_boundary) MeshTools::find_boundary_nodes (mesh, on_boundary);

  // Now calculate the minimum distance to
  // neighboring nodes for each node.
  // hmin holds these distances.
  std::vector<float> hmin (mesh.n_nodes(), 1.e20);

  MeshBase::element_iterator       el  = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end(); 

  for (; el!=end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      hmin[(*el)->node(n)] = std::min(hmin[(*el)->node(n)],
                                      static_cast<float>((*el)->hmin()));
  
  
  // Now actually move the nodes
  {
    const unsigned int seed = 123456;
    
    // seed the random number generator
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
		     static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX),
		     ((mesh.mesh_dimension() == 3) ?
		      static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX) :
		      0.)
		     );
	  
	  dir(0) = (dir(0)-.5)*2.;
	  dir(1) = (dir(1)-.5)*2.;
	  if (mesh.mesh_dimension() == 3)
	    dir(2) = (dir(2)-.5)*2.;
	  
	  dir = dir.unit();

          mesh.node(n)(0) += dir(0)*factor*hmin[n];
          mesh.node(n)(1) += dir(1)*factor*hmin[n];
          if (mesh.mesh_dimension() == 3)
            mesh.node(n)(2) += dir(2)*factor*hmin[n];
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

  for (unsigned int n=0; n<mesh.n_nodes(); n++)
    mesh.node(n) += p;
}


// void MeshTools::Modification::rotate2D (MeshBase& mesh,
//                                         const Real alpha)
// {
//   assert (mesh.mesh_dimension() != 1);

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
  assert (mesh.mesh_dimension() != 1);

  const Real pi = std::acos(-1.);
  const Real  p = -phi/180.*pi;
  const Real  t = -theta/180.*pi;
  const Real  s = -psi/180.*pi;
  const Real sp = std::sin(p), cp = std::cos(p);
  const Real st = std::sin(t), ct = std::cos(t);
  const Real ss = std::sin(s), cs = std::cos(s);

  for (unsigned int n=0; n<mesh.n_nodes(); n++)
    {
      const Point p = mesh.node(n);
      const Real  x = p(0);
      const Real  y = p(1);
      const Real  z = p(2);
      mesh.node(n) = Point(( cp*cs-sp*ct*ss)*x + ( sp*cs+cp*ct*ss)*y + (st*ss)*z,
                           (-cp*ss-sp*ct*cs)*x + (-sp*st+cp*ct*cs)*y + (st*cs)*z,
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
      assert (zs == 0.);

      y_scale = z_scale = x_scale;
    }

  // Scale the x coordinate in all dimensions
  for (unsigned int n=0; n<mesh.n_nodes(); n++)
    mesh.node(n)(0) *= x_scale;


  // Only scale the y coordinate in 2 and 3D
  if (mesh.spatial_dimension() > 1)
    {

      for (unsigned int n=0; n<mesh.n_nodes(); n++)
	mesh.node(n)(1) *= y_scale;

      // Only scale the z coordinate in 3D
      if (mesh.spatial_dimension() == 3)
	{
	  for (unsigned int n=0; n<mesh.n_nodes(); n++)
	    mesh.node(n)(2) *= z_scale;
	}
    }
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
   * Loop over the high-ordered elements.
   * First make sure they _are_ indeed high-order, and then replace
   * them with an equivalent first-order element.
   */
  const_element_iterator endit = elements_end();
  for (const_element_iterator it = elements_begin();
       it != endit; ++it)
    {
      Elem* so_elem = *it;

      assert (so_elem != NULL);

      /*
       * build the first-order equivalent, add to
       * the new_elements list.
       */
      Elem *newparent = so_elem->parent();
      Elem* lo_elem = Elem::build
        (Elem::first_order_equivalent_type
          (so_elem->type()), newparent).release();

#ifdef ENABLE_AMR
      /*
       * Add this element to it's parent if it has one
       */
      if (newparent)
        newparent->add_child(lo_elem);

      /*
       * Reset the parent links of any child elements
       */
      if (so_elem->has_children())
        {
          for (unsigned int c=0; c != so_elem->n_children(); ++c)
            so_elem->child(c)->set_parent(lo_elem);
        }

      /*
       * Copy as much data to the new element as makes sense
       */
      lo_elem->set_p_level(so_elem->p_level());
      lo_elem->set_refinement_flag(so_elem->refinement_flag());
      lo_elem->set_p_refinement_flag(so_elem->p_refinement_flag());
#endif

      assert (lo_elem->n_vertices() == so_elem->n_vertices());

      /*
       * By definition the vertices of the linear and
       * second order element are identically numbered.
       * transfer these.
       */
      for (unsigned int v=0; v < so_elem->n_vertices(); v++)
	lo_elem->set_node(v) = so_elem->get_node(v);

      /**
       * If the second order element had any boundary conditions they
       * should be transfered to the first-order element.  The old
       * boundary conditions will be removed from the BoundaryInfo
       * data structure by insert_elem.
       */
      assert (lo_elem->n_sides() == so_elem->n_sides());
	
      for (unsigned int s=0; s<so_elem->n_sides(); s++)
	{
	  const short int boundary_id =
	    this->boundary_info->boundary_id (so_elem, s);
	    
	  if (boundary_id != this->boundary_info->invalid_id)
	    this->boundary_info->add_side (lo_elem, s, boundary_id);
	}

      /*
       * The new first-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the second-order element.
       */
      lo_elem->set_id(so_elem->id());
      this->insert_elem(lo_elem);
    }

  STOP_LOG("all_first_order()", "Mesh");

  // delete or renumber nodes, etc
  this->prepare_for_use();
}



void UnstructuredMesh::all_second_order (const bool full_ordered)
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and 
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  /*
   * If the mesh is empty or already second order
   * then we have nothing to do
   */
  if (!this->n_elem() ||
      (*(this->elements_begin()))->default_order() != FIRST)
    return;

  // does this work also in parallel?
  // assert (this->n_processors() == 1);

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
  std::map<std::vector<unsigned int>, Node*> adj_vertices_to_so_nodes;

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
      this->reserve_nodes(static_cast<unsigned int>(1.5*this->n_nodes()));
      break;

    case 2:
      /*
       * in 2D, either refine from Tri3 to Tri6 (double the nodes)
       * or from Quad4 to Quad8 (again, double) or Quad9 (2.25 that much)
       */
      this->reserve_nodes(static_cast<unsigned int>(2*this->n_nodes()));
      break;


    case 3:
      /*
       * in 3D, either refine from Tet4 to Tet10 (factor = 2.5) up to
       * Hex8 to Hex27 (something  > 3).  Since in 3D there _are_ already
       * quite some nodes, and since we do not want to overburden the memory by
       * a too conservative guess, use the lower bound
       */
      this->reserve_nodes(static_cast<unsigned int>(2.5*this->n_nodes()));
      break;
	
    default:
      // Hm?
      error();
  }



  /*
   * form a vector that will hold the node id's of
   * the vertices that are adjacent to the son-th
   * second-order node.  Pull this outside of the
   * loop so that silly compilers don't repeatedly
   * create and destroy the vector.
   */
  std::vector<unsigned int> adjacent_vertices_ids;


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
      Elem* lo_elem = *it;

      assert (lo_elem != NULL);

      // make sure it is linear order
      if (lo_elem->default_order() != FIRST)
        {	  
	  std::cerr << "ERROR: This is not a linear element: type=" 
		    << lo_elem->type() << std::endl;
	  error();
	}

      // this does _not_ work for refined elements
      assert (lo_elem->level () == 0);

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

      assert (lo_elem->n_vertices() == so_elem->n_vertices());


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
	  std::pair<std::map<std::vector<unsigned int>, Node*>::iterator,
                    std::map<std::vector<unsigned int>, Node*>::iterator>	    
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
	      
	      // add the new point to the mesh
	      Node* so_node = this->add_point (new_location);

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
	      assert (pos.first->second != NULL);
	      
	      so_elem->set_node(son) = pos.first->second;
	    }
	}


      /**
       * If the linear element had any boundary conditions they
       * should be transfered to the second-order element.  The old
       * boundary conditions will be removed from the BoundaryInfo
       * data structure by insert_elem.
       */
      assert (lo_elem->n_sides() == so_elem->n_sides());
	
      for (unsigned int s=0; s<lo_elem->n_sides(); s++)
	{
	  const short int boundary_id =
	    this->boundary_info->boundary_id (lo_elem, s);
	    
	  if (boundary_id != this->boundary_info->invalid_id)
	    this->boundary_info->add_side (so_elem, s, boundary_id);
	}

      /*
       * The new second-order element is ready.
       * Inserting it into the mesh will replace and delete
       * the first-order element.
       */
      so_elem->set_id(lo_elem->id());
      this->insert_elem(so_elem);
    }


  // we can clear the map
  adj_vertices_to_so_nodes.clear();


  STOP_LOG("all_second_order()", "Mesh");

  // renumber nodes, elements etc
  this->prepare_for_use();
}






void MeshTools::Modification::all_tri (MeshBase& mesh)
{
  // We store pointers to the newly created elements in a vector
  // until they are ready to be added to the mesh.  This is because
  // adding new elements on the fly can cause reallocation and invalidation
  // of existing iterators.
  std::vector<Elem*> new_elements;
  new_elements.reserve (2*mesh.n_active_elem());

  // If the original mesh has boundary data, we carry that over
  // to the new mesh with triangular elements.
  const bool mesh_has_boundary_data = (mesh.boundary_info->n_boundary_ids() > 0);

  // Temporary vectors to store the new boundary element pointers, side numbers, and boundary ids
  std::vector<Elem*> new_bndry_elements;
  std::vector<unsigned short int> new_bndry_sides;
  std::vector<short int> new_bndry_ids;
  
  // Iterate over the elements, splitting QUADS into
  // pairs of conforming triangles.
  {
    MeshBase::element_iterator       el  = mesh.elements_begin();
    const MeshBase::element_iterator end = mesh.elements_end(); 

    for (; el!=end; ++el)
      {
	const ElemType etype = (*el)->type();

	// all_tri currently only works on coarse meshes
	assert ((*el)->parent() == NULL);

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
	      if (((*el)->point(0) - (*el)->point(2)).size() <
		  ((*el)->point(1) - (*el)->point(3)).size())
		{	      
		  tri0->set_node(0) = (*el)->get_node(0);
		  tri0->set_node(1) = (*el)->get_node(1);
		  tri0->set_node(2) = (*el)->get_node(2);
		
		  tri1->set_node(0) = (*el)->get_node(0);
		  tri1->set_node(1) = (*el)->get_node(2);
		  tri1->set_node(2) = (*el)->get_node(3);
		}
	    
	      else
		{
		  edge_swap=true;
		  
		  tri0->set_node(0) = (*el)->get_node(0);
		  tri0->set_node(1) = (*el)->get_node(1);
		  tri0->set_node(2) = (*el)->get_node(3);
		
		  tri1->set_node(0) = (*el)->get_node(1);
		  tri1->set_node(1) = (*el)->get_node(2);
		  tri1->set_node(2) = (*el)->get_node(3);
		}

	      
	      break;
	    }
      
	  case QUAD8:
	    {
	      split_elem =  true;
	      
	      tri0 = new Tri6;
	      tri1 = new Tri6;
	  
	      Node* new_node = mesh.add_point((mesh.node((*el)->node(0)) +
					       mesh.node((*el)->node(1)) +
					       mesh.node((*el)->node(2)) +
					       mesh.node((*el)->node(3)))*.25
					       );
	  
	      // Check for possible edge swap
	      if (((*el)->point(0) - (*el)->point(2)).size() <
		  ((*el)->point(1) - (*el)->point(3)).size())
		{	      
		  tri0->set_node(0) = (*el)->get_node(0);
		  tri0->set_node(1) = (*el)->get_node(1);
		  tri0->set_node(2) = (*el)->get_node(2);
		  tri0->set_node(3) = (*el)->get_node(4);
		  tri0->set_node(4) = (*el)->get_node(5);
		  tri0->set_node(5) = new_node;
	      
		  tri1->set_node(0) = (*el)->get_node(0);
		  tri1->set_node(1) = (*el)->get_node(2);
		  tri1->set_node(2) = (*el)->get_node(3);
		  tri1->set_node(3) = new_node;
		  tri1->set_node(4) = (*el)->get_node(6);
		  tri1->set_node(5) = (*el)->get_node(7);

		}
	  
	      else
		{
		  edge_swap=true;
		  
		  tri0->set_node(0) = (*el)->get_node(3);
		  tri0->set_node(1) = (*el)->get_node(0);
		  tri0->set_node(2) = (*el)->get_node(1);
		  tri0->set_node(3) = (*el)->get_node(7);
		  tri0->set_node(4) = (*el)->get_node(4);
		  tri0->set_node(5) = new_node;
	      
		  tri1->set_node(0) = (*el)->get_node(1);
		  tri1->set_node(1) = (*el)->get_node(2);
		  tri1->set_node(2) = (*el)->get_node(3);
		  tri1->set_node(3) = (*el)->get_node(5);
		  tri1->set_node(4) = (*el)->get_node(6);
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
	      if (((*el)->point(0) - (*el)->point(2)).size() <
		  ((*el)->point(1) - (*el)->point(3)).size())
		{	      
		  tri0->set_node(0) = (*el)->get_node(0);
		  tri0->set_node(1) = (*el)->get_node(1);
		  tri0->set_node(2) = (*el)->get_node(2);
		  tri0->set_node(3) = (*el)->get_node(4);
		  tri0->set_node(4) = (*el)->get_node(5);
		  tri0->set_node(5) = (*el)->get_node(8);
	      
		  tri1->set_node(0) = (*el)->get_node(0);
		  tri1->set_node(1) = (*el)->get_node(2);
		  tri1->set_node(2) = (*el)->get_node(3);
		  tri1->set_node(3) = (*el)->get_node(8);
		  tri1->set_node(4) = (*el)->get_node(6);
		  tri1->set_node(5) = (*el)->get_node(7);
		}

	      else
		{
		  edge_swap=true;
		  
		  tri0->set_node(0) = (*el)->get_node(0);
		  tri0->set_node(1) = (*el)->get_node(1);
		  tri0->set_node(2) = (*el)->get_node(3);
		  tri0->set_node(3) = (*el)->get_node(4);
		  tri0->set_node(4) = (*el)->get_node(8);
		  tri0->set_node(5) = (*el)->get_node(7);
	      
		  tri1->set_node(0) = (*el)->get_node(1);
		  tri1->set_node(1) = (*el)->get_node(2);
		  tri1->set_node(2) = (*el)->get_node(3);
		  tri1->set_node(3) = (*el)->get_node(5);
		  tri1->set_node(4) = (*el)->get_node(6);
		  tri1->set_node(5) = (*el)->get_node(8);
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
	      std::cerr << "Warning, encountered non-2D element "
                        << Utility::enum_to_string<ElemType>(etype)
			<< " in MeshTools::Modification::all_tri(), hope that's OK..."
			<< std::endl;
	    }
	  } // end switch (etype)

	

	if (split_elem)
	  {
	    // If this mesh has boundary data, be sure the correct
	    // ID's are also set for tri0 and tri1.
	    if (mesh_has_boundary_data)
	      {
		for (unsigned int sn=0; sn<(*el)->n_sides(); ++sn)
		  {
		    short int b_id = mesh.boundary_info->boundary_id(*el, sn);

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
				  std::cerr << "Quad4/8/9 cannot have more than 4 sides." << std::endl;
				  error();
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
				  std::cerr << "Quad4/8/9 cannot have more than 4 sides." << std::endl;
				  error();
				}
			      }
			  } // end edge_swap==true
		      } // end if (b_id != BoundaryInfo::invalid_id)
		  } // end for loop over sides

		// Remove the original element from the BoundaryInfo structure.
		mesh.boundary_info->remove(*el);

	      } // end if (mesh_has_boundary_data)
	    
	    // Add the newly-created triangles to the temporary vector of new elements.
	    new_elements.push_back(tri0);
	    new_elements.push_back(tri1);

	    // Delete the original element
	    mesh.delete_elem(*el);
	  } // end if (split_elem)
      } // End for loop over elements
  } 

  
    // Now, iterate over the new elements vector, and add them each to
    // the Mesh.
  {
    std::vector<Elem*>::iterator el        = new_elements.begin();
    const std::vector<Elem*>::iterator end = new_elements.end();
    for (; el != end; ++el)
      {
	mesh.add_elem(*el);
      }
  }

  
  if (mesh_has_boundary_data)
    {
      // By this time, we should have removed all of the original boundary sides
      // - except on a hybrid mesh, where we can't "start from a blank slate"! - RHS
      // assert (mesh.boundary_info->n_boundary_conds()==0);

      // Clear the boundary info, to be sure and start from a blank slate.
      // mesh.boundary_info->clear();

      // If the old mesh had boundary data, the new mesh better have some.
      assert (new_bndry_elements.size() > 0);

      // We should also be sure that the lengths of the new boundary data vectors
      // are all the same.
      assert (new_bndry_elements.size() == new_bndry_sides.size());
      assert (new_bndry_sides.size()    == new_bndry_ids.size());

      // Add the new boundary info to the mesh
      for (unsigned int s=0; s<new_bndry_elements.size(); ++s)
	mesh.boundary_info->add_side(new_bndry_elements[s],
				     new_bndry_sides[s],
				     new_bndry_ids[s]);
    }
  

  // Prepare the newly created mesh for use.
  mesh.prepare_for_use();

  // Let the new_elements and new_bndry_elements vectors go out of scope.
}


void MeshTools::Modification::smooth (MeshBase& mesh,
                                      const unsigned int n_iterations,
                                      const Real power)
{
  /**
   * This implementation assumes every element "side" has only 2 nodes.
   */
  assert (mesh.mesh_dimension() == 2);
  
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

                            const unsigned int id0 = node0->id(), id1 = node1->id();
                            new_positions[id0].add_scaled( *node1, node_weight );
                            new_positions[id1].add_scaled( *node0, node_weight );
                            weight[id0] += node_weight;
                            weight[id1] += node_weight;
                          }
                      } // element neighbor loop
                  } 
#ifdef ENABLE_AMR
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

                                const unsigned int id = elem->get_node(nc)->id();
                                new_positions[id] = point;
                                weight[id] = 1.;
                              }
                    
                          } // if parent->child == elem
                      } // for parent->n_children
                  } // if element refinement_level
#endif // #ifdef ENABLE_AMR

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

                    const unsigned int id = elem->get_node(n)->id();
                    mesh.node(id) = point/n_adjacent_vertices;
                  }
              }
          }
      
        } // refinement_level loop

    } // end iteration
}



#ifdef ENABLE_AMR
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
  std::vector<short int>          saved_bc_ids;
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

	// Set node pointers (they point to nodes in the original mesh)
	// The copy's element ID will be set upon adding it to the mesh
	for(unsigned int n=0; n<elem->n_nodes(); n++)
	  copy->set_node(n) = elem->get_node(n);

	// This element could have boundary info as well.  We need
	// to save the (elem, side, bc_id) triples
	for (unsigned int s=0; s<elem->n_sides(); s++)
	    if (elem->neighbor(s) == NULL)
	      {
		short int bc_id = mesh.boundary_info->boundary_id (elem,s);

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
    assert (saved_boundary_elements.size() == saved_bc_ids.size());
    assert (saved_bc_ids.size()            == saved_bc_sides.size());
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
      mesh.add_elem(*it);
  }

  // Finally, also add back the saved boundary information
  for (unsigned int e=0; e<saved_boundary_elements.size(); ++e)
    mesh.boundary_info->add_side(saved_boundary_elements[e],
				 saved_bc_sides[e],
				 saved_bc_ids[e]);

  // Trim unused and renumber nodes and elements
  mesh.prepare_for_use();
}
#endif // #ifdef ENABLE_AMR
