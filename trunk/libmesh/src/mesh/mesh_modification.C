// $Id: mesh_modification.C,v 1.15 2005-08-26 21:00:35 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include <cmath> // for acos()
#include <algorithm>
#include <map>

// Local includes
#include "mesh_tools.h"
#include "mesh_modification.h"
#include "mesh.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "libmesh_logging.h"
#include "boundary_info.h"



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
    srand(seed);
    
    // If the node is on the boundary or
    // the node is not used by any element (hmin[n]<1.e20)
    // then we should not move it.
    // [Note: Testing for (in)equality might be wrong
    // (different types, namely float and double)]
    for (unsigned int n=0; n<mesh.n_nodes(); n++)
      if (!on_boundary[n] && (hmin[n] < 1.e20) )
	{
	  // the direction, random but unit normalized
	  
	  Point dir( static_cast<Real>(rand())/static_cast<Real>(RAND_MAX),
		     static_cast<Real>(rand())/static_cast<Real>(RAND_MAX),
		     ((mesh.mesh_dimension() == 3) ?
		      static_cast<Real>(rand())/static_cast<Real>(RAND_MAX) :
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

//   const Real pi = acos(-1);
//   const Real  a = alpha/180.*pi;
//   for (unsigned int n=0; n<mesh.n_nodes(); n++)
//     {
//       const Point p = mesh.node(n);
//       const Real  x = p(0);
//       const Real  y = p(1);
//       const Real  z = p(2);
//       mesh.node(n) = Point(cos(a)*x - sin(a)*y,
//                            sin(a)*x + cos(a)*y,
//                            z);
//     }

// }



void MeshTools::Modification::rotate (MeshBase& mesh,
				      const Real phi,
				      const Real theta,
				      const Real psi)
{
  assert (mesh.mesh_dimension() != 1);

  const Real pi = acos(-1.);
  const Real  p = -phi/180.*pi;
  const Real  t = -theta/180.*pi;
  const Real  s = -psi/180.*pi;
  const Real sp = sin(p), cp = cos(p);
  const Real st = sin(t), ct = cos(t);
  const Real ss = sin(s), cs = cos(s);

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
// Mesh class member functions for mesh modification
void Mesh::all_second_order (const bool full_ordered)
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and 
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  // does this work also in parallel?
  // assert (this->n_processors() == 1);

  START_LOG("all_second_order()", "MeshBase");

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
      this->_nodes.reserve(static_cast<unsigned int>(1.5*this->_nodes.size()));
      break;

    case 2:
      /*
       * in 2D, either refine from Tri3 to Tri6 (double the nodes)
       * or from Quad4 to Quad8 (again, double) or Quad9 (2.25 that much)
       */
      this->_nodes.reserve(static_cast<unsigned int>(2*this->_nodes.size()));
      break;


    case 3:
      /*
       * in 3D, either refine from Tet4 to Tet10 (factor = 2.5) up to
       * Hex8 to Hex27 (something  > 3).  Since in 3D there _are_ already
       * quite some nodes, and since we do not want to overburden the memory by
       * a too conservative guess, use the lower bound
       */
      this->_nodes.reserve(static_cast<unsigned int>(2.5*this->_nodes.size()));
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
  for (unsigned int e=0; e<_elements.size(); e++)
    {
      // the linear-order element
      Elem* lo_elem = _elements[e];

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
       * should be transfered to the second-order element, and then
       * removed from the BoundaryInfo data structure.
       */
      {
	assert (lo_elem->n_sides() == so_elem->n_sides());
	
	for (unsigned int s=0; s<lo_elem->n_sides(); s++)
	  {
	    const short int boundary_id =
	      this->boundary_info->boundary_id (lo_elem, s);
	    
	    if (boundary_id != this->boundary_info->invalid_id)
	      this->boundary_info->add_side (so_elem, s, boundary_id);
	  }
	
	/**
	 * We have taken any boundary conditions the low-ordered
	 * element may have had.  Since we are about to delete
	 * the low-ordered element, we should first un-associate
	 * any boundary conditions it has.
	 */
	this->boundary_info->remove (lo_elem);
      }

      
      /*
       * The new second-order element is ready.
       * Delete the linear element and replace it with
       * the second-order element.
       */
      delete lo_elem;
      _elements[e] = so_elem;
    }


  // we can clear the map
  adj_vertices_to_so_nodes.clear();


  STOP_LOG("all_second_order()", "MeshBase");

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

  // Iterate over the elements, splitting QUADS into
  // pairs of conforming triangles.
  {
    MeshBase::element_iterator       el  = mesh.elements_begin();
    const MeshBase::element_iterator end = mesh.elements_end(); 

    for (; el!=end; ++el)
      {
	const ElemType etype = (*el)->type();

	switch (etype)
	  {
	  case QUAD4:
	    {
	      Elem* tri0 = new Tri3;
	      Elem* tri1 = new Tri3;
	    
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
		  tri0->set_node(0) = (*el)->get_node(0);
		  tri0->set_node(1) = (*el)->get_node(1);
		  tri0->set_node(2) = (*el)->get_node(3);
		
		  tri1->set_node(0) = (*el)->get_node(1);
		  tri1->set_node(1) = (*el)->get_node(2);
		  tri1->set_node(2) = (*el)->get_node(3);
		}

	      new_elements.push_back(tri0);
	      new_elements.push_back(tri1);

	      mesh.delete_elem(*el);
	      break;
	    }
      
	  case QUAD8:
	    {
	      Elem* tri0 = new Tri6;
	      Elem* tri1 = new Tri6;
	  
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
	  
	      new_elements.push_back(tri0);
	      new_elements.push_back(tri1);

	      mesh.delete_elem(*el);
	      break;
	    }
      
	  case QUAD9:
	    {
	      Elem* tri0 = new Tri6;
	      Elem* tri1 = new Tri6;

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
	    
	      new_elements.push_back(tri0);
	      new_elements.push_back(tri1);
	      
	      mesh.delete_elem(*el);
	      break;
	    }
	  
	  default:
	    {
	      // If not one of the QUAD* types, the Elem must
	      // be a TRI* type already, or a 3D element, so just leave it.
	    }
	  }
      }
  }

  {
    // Now, iterate over the new elements vector, and add them each to
    // the Mesh.
    std::vector<Elem*>::iterator el        = new_elements.begin();
    const std::vector<Elem*>::iterator end = new_elements.end();
    for (; el != end; ++el)
      {
	mesh.add_elem(*el);
      }
  }
  

  // Prepare the newly created mesh for use.
  mesh.prepare_for_use();

  // Let the new_elements vector go out of scope.
}
