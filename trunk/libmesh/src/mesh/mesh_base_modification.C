// $Id: mesh_base_modification.C,v 1.8 2003-09-02 18:02:43 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
#include <algorithm>
#include <map>


// Local includes
#include "mesh_base.h"
#include "libmesh.h"
#include "face_inf_quad4.h"
#include "face_inf_quad6.h"
#include "cell_inf_prism6.h"
#include "cell_inf_prism12.h"
#include "cell_inf_hex8.h"
#include "cell_inf_hex16.h"
#include "cell_inf_hex18.h"
#include "mesh_logging.h"


// ------------------------------------------------------------
// Mesh class member functions for mesh modification
void MeshBase::all_second_order (const bool full_ordered)
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and 
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
      renumber_nodes_and_elements ();

  // does this work also in parallel?
  assert (this->n_processors() == 1);

  START_LOG("all_second_order()", "MeshBase");

  /*
   * the vector holding the new second-order
   * replacement elements
   */
  std::vector<Elem*> new_elements;
  new_elements.reserve(this->n_elem());

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
   * iterate over all elements contained in the 
   * mesh
   */
  elem_iterator       old_elements_it  (this->elements_begin());
  const elem_iterator old_elements_end (this->elements_end());

  for (; old_elements_it != old_elements_end; ++old_elements_it)
    {
      // the linear-order element
      Elem* lo_elem = *old_elements_it;

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
      Elem* so_elem = Elem::build ( Elem::second_order_equivalent_type(lo_elem->type(), 
								       full_ordered) );
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
	  const unsigned int n_adjacent_vertices = so_elem->n_second_order_adjacent_vertices(son);

	  /*
	   * form a vector that will hold the node id's of
	   * the vertices that are adjacent to the son-th
	   * second-order node
	   */
	  std::vector<unsigned int> adjacent_vertices_ids;
	  adjacent_vertices_ids.resize(n_adjacent_vertices);
	  for (unsigned int v=0; v<n_adjacent_vertices; v++)
	      adjacent_vertices_ids[v] = so_elem->node( so_elem->second_order_adjacent_vertex(son,v) );

	  /*
	   * \p adjacent_vertices_ids is now in order of the current
	   * side.  sort it, so that comparisons  with the 
	   * \p adjacent_vertices_ids created through other elements' 
	   * sides can match
	   */
	  std::sort(adjacent_vertices_ids.begin(), adjacent_vertices_ids.end());


	  // does this set of vertices already has a mid-node added?
	  std::map<std::vector<unsigned int>, Node*>::const_iterator pos =  
	      adj_vertices_to_so_nodes.find(adjacent_vertices_ids);

	  if (pos == adj_vertices_to_so_nodes.end())
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
	      adj_vertices_to_so_nodes.insert(std::make_pair(adjacent_vertices_ids, so_node));

	      so_elem->set_node(son) = so_node;
	    }
	  else
	      so_elem->set_node(son) = pos->second;

	}


      /*
       * The new second-order element is ready.
       * Add it to the new_elements vector
       */
      new_elements.push_back(so_elem);

    }


  // we can clear the map
  adj_vertices_to_so_nodes.clear();


  /*
   * the \p _elements vector has to be replaced
   * by the \p new_elements vector.  Delete the
   * old element, then put the new element in place.
   */
  {
    for (unsigned int e=0; e<_elements.size(); e++)
      {
	assert (_elements[e] != NULL);
	delete _elements[e];
	_elements[e] = new_elements[e];
      }

    // now we can safely clear our local \p new_elements vector
    new_elements.clear();
  }


  STOP_LOG("all_second_order()", "MeshBase");

  // renumber nodes, elements etc
  this->prepare_for_use();
}








#ifdef ENABLE_INFINITE_ELEMENTS


const Point MeshBase::build_inf_elem(bool be_verbose)
{
  START_LOG("build_inf_elem()", "MeshBase");

  // determine origin automatically,
  // works only if the mesh has no symmetry planes.
  std::pair<Point, Point> b_box = bounding_box();
  Point origin = (b_box.first+b_box.second)/2.;
    
  if (be_verbose)
    {
#ifdef DEBUG
      std::cout << " Determined origin for Infinite Elements:" 
		<< std::endl
		<< "  ";
      origin.write_unformatted(std::cout);
      std::cout << std::endl;
#endif
    }

  build_inf_elem(origin, false, false, false, be_verbose);

  STOP_LOG("build_inf_elem()", "MeshBase");

  /*
   * when finished with building the Ifems,
   * it remains to prepare the mesh for use:
   * find neighbors (again), partition (if needed)...
   */
  this->prepare_for_use ();

  return origin;
}







const Point MeshBase::build_inf_elem (const InfElemOriginValue& origin_x,
				      const InfElemOriginValue& origin_y,
				      const InfElemOriginValue& origin_z,
				      const bool x_sym,
				      const bool y_sym,
				      const bool z_sym,
				      const bool be_verbose,
				      std::vector<const Node*>* inner_boundary_nodes)
{
  START_LOG("build_inf_elem()", "MeshBase");

  /*
   * first determine the origin of the 
   * infinite elements.  For this, the
   * origin defaults to the given values,
   * and may be overridden when the user
   * provided values
   */
  Point origin(origin_x.second, origin_y.second, origin_z.second);

  /*
   * when only _one_ of the origin coordinates is _not_
   * given, we have to determine it on our own
   */
  if ( !origin_x.first || !origin_y.first || !origin_z.first)
    {
      // determine origin
      std::pair<Point, Point> b_box = this->bounding_box();
      const Point auto_origin = (b_box.first+b_box.second)/2.;

      // override default values, if necessary
      if (!origin_x.first)
	  origin(0) = auto_origin(0);
      if (!origin_y.first)
	  origin(1) = auto_origin(1);
      if (!origin_z.first)
	  origin(2) = auto_origin(2);

      if (be_verbose)
        {
	  std::cout << " Origin for Infinite Elements:" << std::endl;

	  if (!origin_x.first)
	      std::cout << "  determined x-coordinate" << std::endl;
	  if (!origin_y.first)
	      std::cout << "  determined y-coordinate" << std::endl;
	  if (!origin_z.first)
	      std::cout << "  determined z-coordinate" << std::endl;

	  std::cout << "  coordinates: ";
	  origin.write_unformatted(std::cout);
	  std::cout << std::endl;
	}
    }

  else if (be_verbose)

    {
      std::cout << " Origin for Infinite Elements:" << std::endl;
      std::cout << "  coordinates: ";
      origin.write_unformatted(std::cout);
      std::cout << std::endl;
    }



  /*
   * Now that we have the origin,
   * check if the user provided an
   * \p inner_boundary_nodes.  If so,
   * we pass a std::set to the
   * actual implementation of the
   * build_inf_elem(), so that we can
   * convert this to the Node* vector
   */
  if (inner_boundary_nodes != NULL)
    {
      /*
       * note that the std::set that we will get
       * from build_inf_elem() uses the index of
       * the element in this->_elements vector,
       * and the second entry is the side index
       * for this element.  Therefore, we do _not_
       * need to renumber nodes and elements
       * prior to building the infinite elements.
       *
       * However, note that this method here uses
       * node id's... Do we need to renumber?
       */


      // Form the list of faces of elements which finally
      // will tell us which nodes should receive boundary
      // conditions (to form the std::vector<const Node*>)
      std::set< std::pair<unsigned int,
	                  unsigned int> > inner_faces;


      // build infinite elements
      this->build_inf_elem(origin, 
			   x_sym, y_sym, z_sym, 
			   be_verbose,
			   &inner_faces);

      if (be_verbose)
        {
	  this->print_info();
	  std::cout << "Data pre-processing:" << std::endl
		    << " convert the <int,int> list to a Node* list..."
		    << std::endl;
	}

      /*
       * First use a std::vector<unsigned int> that holds
       * the global node numbers.  Then sort this vector,
       * so that it can be made unique (no multiple occurence
       * of a node), and then finally insert the Node* in
       * the vector inner_boundary_nodes.
       *
       * Reserve memory for the vector<unsigned int> with 
       * 4 times the size of the number of elements in the 
       * std::set. This is a good bet for Quad4 face elements.  
       * For higher-order elements, this probably _has_ to lead
       * to additional allocations...
       * Practice has to show how this affects performance.
       */
      std::vector<unsigned int> inner_boundary_node_numbers;
      inner_boundary_node_numbers.reserve(4*inner_faces.size());

      /*
       * Now transform the set of pairs to a list of (possibly
       * duplicate) global node numbers.
       */
      std::set< std::pair<unsigned int,unsigned int> > :: iterator face_it = inner_faces.begin();
      const std::set< std::pair<unsigned int,unsigned int> > :: iterator face_end = inner_faces.end();
      for(; face_it!=face_end;++face_it)
        {
	  std::pair<unsigned int,unsigned int> p=*face_it;
	
	  // build a full-ordered side element to get _all_ the base nodes
	  AutoPtr<Elem> side( this->elem(p.first)->build_side(p.second) ); 

	  // insert all the node numbers in inner_boundary_node_numbers
	  for (unsigned int n=0; n< side->n_nodes(); n++)
	      inner_boundary_node_numbers.push_back(side->node(n));
	}


      /*
       * inner_boundary_node_numbers now still holds multiple entries of
       * node numbers.  So first sort, then unique the vector.
       * Note that \p std::unique only puts the new ones in
       * front, while to leftovers are @e not deleted.  Instead,
       * it returns a pointer to the end of the unique range.
       */
      const unsigned int ibn_size_before = inner_boundary_node_numbers.size();
      std::sort (inner_boundary_node_numbers.begin(), inner_boundary_node_numbers.end());
      std::vector<unsigned int>::iterator unique_end = 
	  std::unique (inner_boundary_node_numbers.begin(), inner_boundary_node_numbers.end());

      const int unique_size = std::distance(inner_boundary_node_numbers.begin(), unique_end);
      assert (unique_size <= static_cast<const int>(ibn_size_before));

      /*
       * Finally, create const Node* in the inner_boundary_nodes
       * vector.  Reserve, not resize (otherwise, the push_back
       * would append the interesting nodes, while NULL-nodes
       * live in the resize'd area...
       */
      inner_boundary_nodes->reserve (unique_size);
      inner_boundary_nodes->clear();


      std::vector<unsigned int>::iterator pos_it = inner_boundary_node_numbers.begin();
      for (; pos_it != unique_end; ++pos_it)
        {
	  const Node& node = this->node(*pos_it);
	  inner_boundary_nodes->push_back(&node);
	}

      if (be_verbose)
	  std::cout << "  finished identifying " << unique_size 
		    << " target nodes." << std::endl;  
    }

  else

    {
      // simply build the infinite elements
      this->build_inf_elem(origin, x_sym, y_sym, z_sym, be_verbose);
    }


  STOP_LOG("build_inf_elem()", "MeshBase");

  /*
   * when finished with building the Ifems,
   * it remains to prepare the mesh for use:
   * find neighbors again, partition (if needed)...
   */
  this->prepare_for_use ();

  return origin;
}







void MeshBase::build_inf_elem(const Point& origin, 
			      const bool x_sym, 
			      const bool y_sym, 
			      const bool z_sym,
			      const bool be_verbose,
			      std::set< std::pair<unsigned int,
			                          unsigned int> >* inner_faces)
{

  if (be_verbose)
    {
#ifdef DEBUG
      std::cout << " Building Infinite Elements:" << std::endl;
      std::cout << "  updating element neighbor tables..." << std::endl;
#else
      std::cout << " Verbose mode disabled in non-debug mode." << std::endl;
#endif
    }


  PAUSE_LOG("build_inf_elem()", "MeshBase");

  this->find_neighbors();	// update elem->neighbor() tables

  RESTART_LOG("build_inf_elem()", "MeshBase");

  // pairs: (first: element number, second: side number)
  std::set< std::pair<unsigned int,unsigned int> > faces,ofaces;
  std::set< std::pair<unsigned int,unsigned int> > :: iterator face_it;
	
  std::set<unsigned int> onodes;
  std::set<unsigned int> :: iterator on_it;
	
  Real max_r=0.;
  unsigned int max_r_node;

#ifdef DEBUG
  if (be_verbose)
    {
      std::cout << "  collecting boundary sides";
      if (x_sym || y_sym || z_sym)
	std::cout << ", skipping sides in symmetry planes..." << std::endl;
      else
	std::cout << "..." << std::endl;
    }
#endif

  /**
   * Iterate through all elements and sides, collect indices of all active
   * boundary sides in the faces set. Skip sides which lie in symmetry planes.
   * Later, sides of the inner boundary will be sorted out.
   * Can't use iterators here since the element
   * index (e) is explicitly used later...
   */
  for(unsigned int e=0;e<this->n_elem();e++)
    {
      if (!(this->elem(e)->active()))
   	continue;
      
      for (unsigned int s=0; s<this->elem(e)->n_neighbors(); s++)
	{
	  if (this->elem(e)->neighbor(s) != NULL)
	    continue;	 // check if elem(e) is on the boundary
	  
	  /* 
	   * note that it is safe to use the Elem::side() method, 
	   * which gives a non-full-ordered element 
	   */
	  AutoPtr<Elem> side(this->elem(e)->side(s));
		
	  bool sym_side=false;		
			
	  bool on_x_sym=true;
	  bool on_y_sym=true;
	  bool on_z_sym=true;
			
		
	  /*
	   * check whether the nodes are on the symmetry planes;
	   * therefore sufficient to use a non-full-ordered side element
	   */
	  for(unsigned int n=0;n<side->n_nodes();n++)
	    {
	
	      Point dist_from_origin=this->point(side->node(n))-origin;

	      if(x_sym)
		if( fabs(dist_from_origin(0)) > 1.e-6 )
		  on_x_sym=false;

	      if(y_sym)
		if( fabs(dist_from_origin(1)) > 1.e-6 )
		  on_y_sym=false;

	      if(z_sym)
		if( fabs(dist_from_origin(2)) > 1.e-6 )
		  on_z_sym=false;

	      //find the node most distant from origin

	      Real r=dist_from_origin.size();
	      if (r>max_r)
		{
		  max_r=r;max_r_node=side->node(n);
		}

	    }

	  sym_side=(x_sym&&on_x_sym)||(y_sym&&on_y_sym)||(z_sym&&on_z_sym);
				
	  std::pair<unsigned int,unsigned int> p(e,s);
	  
	  if (!sym_side)
	    faces.insert(p);					
					    

	}	// sides
    }   // elems

			
  /* 
   *	If a boundary side has one node on the outer boundary,
   *	all points of this side are on the outer boundary.
   *
   *	Start with the node most distant from origin, which has
   *	to be on the outer boundary, then recursively find all
   *	sides and nodes connected to it. Found sides are moved
   *  from faces to ofaces, nodes are collected in onodes.  
   *
   *  Here, the search is done iterative, because, depending on
   *  the mesh, a very high level of recursion might be necessary.
   */		
	 
  onodes.insert(max_r_node);	
		
  face_it = faces.begin();
  unsigned int facesfound=0;
	
  do {
			
    std::pair<unsigned int,unsigned int> p;
    p=*face_it;

    /*
     * This has to be a full-ordered side element,
     * since we need the correct n_nodes,
     */
    AutoPtr<Elem> side(elem(p.first)->build_side(p.second));
		
    bool found=false;
    for(unsigned int sn=0; sn<side->n_nodes(); sn++)
      if(onodes.count(side->node(sn))) {found=true;break;}
    	
    	
    /* If a new oface is found, include it's nodes in onodes */		
    	
    if(found)		
      {
	for(unsigned int sn=0;sn<side->n_nodes();sn++)
	  onodes.insert(side->node(sn));
    				
	ofaces.insert(p);
	face_it++;			// iteration is done here
	faces.erase(p);
    		
	facesfound++;
      }
    else face_it++;			// iteration is done here


    /* If at least one new oface was found in this cycle, 
     * do another search cycle. */

    if(facesfound>0&&face_it==faces.end()){	
      facesfound=0;
      face_it=faces.begin();
    }
  }
  while(face_it!=faces.end());
	
	
#ifdef DEBUG
  if (be_verbose)
    std::cout << "  found " 
	      << faces.size() 
	      << " inner and " 
	      << ofaces.size() 
	      << " outer boundary faces" 
	      << std::endl;
#endif	
	
  /*
   * When the user provided a non-null pointer to
   * inner_faces, that implies he wants to have
   * this std::set.  For now, simply copy the data.
   */
  if (inner_faces != NULL)
      *inner_faces = faces;

  /*
   * free memory, clear our local variable, no need
   * for it any more.
   */
  faces.clear();


  // outer_nodes maps onodes to their duplicates

  std::map<unsigned int, Node *> outer_nodes;


  /**
   * for each boundary node, add an outer_node with 
   * double distance from origin.
   */


  for(on_it=onodes.begin();on_it!=onodes.end();++on_it)
    {
      Point p=Point(this->point(*on_it))*2.-origin;
      outer_nodes[*on_it]=this->add_point(p);
    }


#ifdef DEBUG
  // for verbose, remember n_elem
  unsigned int _n_conventional_elem = this->n_elem();
#endif


  /**
   * build Elems based on boundary side type
   */
  for(face_it=ofaces.begin();face_it!=ofaces.end();++face_it)
    {
	
      std::pair<unsigned int,unsigned int> p=*face_it;
	
      // build a full-ordered side element to get the base nodes
      AutoPtr<Elem> side(elem(p.first)->build_side(p.second)); 

      bool is_higher_order_elem=false;


      /*
       * create cell depending on side type, assign nodes,
       * use braces to force scope
       */
      {
        Elem* el;
	switch(side->type())
	{
	  // 3D infinite elements
	  // TRIs					
	  case TRI3:	
	    el=new InfPrism6;
	    break;
					 		
	  case TRI6: 
	    el=new InfPrism12; 
	    is_higher_order_elem=true;
	    break;
							
	  // QUADs					
	  case QUAD4: 
	    el=new InfHex8;
	    break;
							
	  case QUAD8: 
	    el=new InfHex16;
	    is_higher_order_elem=true;
	    break;
							
	  case QUAD9: 
	    el=new InfHex18;
	    /*
	     * the method of assigning nodes (which follows below)
	     * omits in the case of QUAD9 the bubble node; therefore
	     * we assign these by hand here
	     */
	    el->set_node(16) = side->get_node(8);
	    el->set_node(17) = outer_nodes[side->node(8)];
	    is_higher_order_elem=true;
	    break;

	  // 2D infinite elements		 
	  case EDGE2:
	    el=new InfQuad4;
	    break;

	  case EDGE3:
	    el=new InfQuad6;
	    el->set_node(4) = side->get_node(2);
	    break;

	  // 1D infinite elements not supported
	  default: 
	    std::cout << "MeshBase::build_inf_elem(Point, bool, bool, bool, bool): invalid face element" 
		      << std::endl;
	    continue;
	}


	/*
	 * assign vertices to the new infinite element
	 */
	const unsigned int n_base_vertices = side->n_vertices();
	for(unsigned int i=0; i<n_base_vertices; i++)
	  {
	    el->set_node(i                ) = side->get_node(i);
	    el->set_node(i+n_base_vertices) = outer_nodes[side->node(i)];
	  }


	/*
	 * when this is a higher order element, 
	 * assign also the nodes in between
	 */
	if (is_higher_order_elem)
          {
	    /* 
	     * n_safe_base_nodes is the number of nodes in \p side
	     * that may be safely assigned using below for loop.
	     * Actually, n_safe_base_nodes is _identical_ with el->n_vertices(),
	     * since for QUAD9, the 9th node was already assigned above
	     */
	    const unsigned int n_safe_base_nodes   = el->n_vertices();

	    for(unsigned int i=n_base_vertices; i<n_safe_base_nodes; i++)
	      {
	        el->set_node(i+n_base_vertices)   = side->get_node(i);
	        el->set_node(i+n_safe_base_nodes) = outer_nodes[side->node(i)];
	      }
	  }
			

 	/*
	 * add infinite element to mesh
	 */
	this->add_elem(el);

      }

    }


#ifdef DEBUG
  if (be_verbose)
    std::cout << "  added "
	      << this->n_elem()-_n_conventional_elem
	      << " infinite elements and "
	      << onodes.size() 
	      << " nodes to the mesh"
	      << std::endl
	      << std::endl;
#endif

}



#endif /* ENABLE_INFINITE_ELEMENTS */



