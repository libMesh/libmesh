// $Id: mesh_base.C,v 1.41 2003-06-24 05:33:51 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



// library configuration
#include "mesh_config.h"

// C++ includes
#include <algorithm>
#include <sstream>
#include <math.h>
#include <set>
#include <map>

#if   defined(HAVE_HASH_MAP)
# include <hash_map>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#endif

// Local includes
#include "mesh_base.h"
#include "libmesh.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_inf_quad4.h"
#include "face_inf_quad6.h"
#include "cell_inf_prism6.h"
#include "cell_inf_prism12.h"
#include "cell_inf_hex8.h"
#include "cell_inf_hex16.h"
#include "cell_inf_hex18.h"
#include "petsc_matrix.h"
#include "mesh_logging.h"
#include "metis_partitioner.h"




// ------------------------------------------------------------
// MeshBase class member functions
MeshBase::MeshBase (unsigned int d) :
#ifdef ENABLE_AMR
  mesh_refinement    (*this),
#endif
  boundary_info      (*this),
  data               (*this),
  mesh_communication (*this),
  _n_sbd             (1),
  _dim               (d),
  _is_prepared       (false)
{
  assert (DIM <= 3);
  assert (DIM >= _dim);
  assert (libMesh::initialized());
}



MeshBase::MeshBase (const MeshBase& other_mesh) :
#ifdef ENABLE_AMR
  mesh_refinement    (*this),
#endif
  boundary_info      (*this),
  data               (*this),
  mesh_communication (*this),
  _nodes             (other_mesh._nodes),
  _elements          (other_mesh._elements),
  _n_sbd             (other_mesh._n_sbd),
  _dim               (other_mesh._dim),
  _is_prepared       (other_mesh._is_prepared)

{
}



MeshBase::~MeshBase()
{
  this->clear();

  assert (!libMesh::closed());
}



void MeshBase::prepare_for_use ()
{
  // Renumber the nodes and elements so that they in
  // contiguous blocks.
  this->renumber_nodes_and_elements();
  
  // Let all the elements find their neighbors
  this->find_neighbors();

  // Partition the mesh.
  this->partition();
  
  // The mesh is now prepared for use.
  _is_prepared = true;
}



Node* MeshBase::add_point (const Point& p,
			   const unsigned int num)
{  
  START_LOG("add_point()", "MeshBase");

  if (num >= this->n_nodes())
    {
      _nodes.push_back (Node::build(p, this->n_nodes()));

      STOP_LOG("add_point()", "MeshBase");
  
      return _nodes.back();
    }
  
  else
    {
      assert (num < this->n_nodes());
      assert (this->node_ptr(num)       != NULL);
      assert (this->node_ptr(num)->id() != Node::invalid_id);
      
      this->node(num) = p;
      
      assert (this->node(num).id() == num);
      
      STOP_LOG("add_point()", "MeshBase");
  
      return this->node_ptr(num);
    }

  
  // We'll never get here...
  error();
  return NULL;
}



void MeshBase::add_elem (Elem* e, const unsigned int n)
{
  START_LOG("add_elem()", "MeshBase");

  if (n >= _elements.size())
    {
      if (e != NULL)
	e->set_id (_elements.size());
      
      _elements.push_back(e);
    }
  
  else
    {
      assert (n < _elements.size());

      if (e != NULL)
	e->set_id (n);
      
      _elements[n] = e;
    }

  STOP_LOG("add_elem()", "MeshBase");
}



unsigned int MeshBase::n_active_elem () const
{
  unsigned int num=0;

  const_active_elem_iterator       el (this->elements_begin());
  const const_active_elem_iterator end(this->elements_end()); 

  for (; el!=end; ++el)
    num++;

  return num;
}



void MeshBase::clear ()
{
  // Clear other data structures
#ifdef ENABLE_AMR
  
  mesh_refinement.clear();
  
#endif

  boundary_info.clear();

  mesh_communication.clear();

  
  // Reset the number of subdomains and the
  // number of processors
  _n_sbd  = 1;

  // Clear the elements data structure
  {
    for (unsigned int e=0; e<_elements.size(); e++)
      if (_elements[e] != NULL)
	{
	  delete _elements[e];
	  _elements[e] = NULL;
	}
    
    _elements.clear();
  }

  // clear the nodes data structure
  {
    for (unsigned int n=0; n<_nodes.size(); n++)
      if (_nodes[n] != NULL)
	{
	  delete _nodes[n];
	  _nodes[n] = NULL;
	}
    
    _nodes.clear();
  }

  // Reset the _is_prepared flag
  _is_prepared = false;
}



unsigned int MeshBase::n_sub_elem () const
{
  unsigned int ne=0;

  const_elem_iterator       el (this->elements_begin());
  const const_elem_iterator end(this->elements_end());
  
  for (; el!=end; ++el)
    ne += (*el)->n_sub_elem(); 

  return ne;
}



unsigned int MeshBase::n_active_sub_elem () const
{
  unsigned int ne=0;

  const_active_elem_iterator       el (this->elements_begin());
  const const_active_elem_iterator end(this->elements_end());
  
  for (; el!=end; ++el)
    ne += (*el)->n_sub_elem(); 

  return ne;
}



std::vector<ElemType> MeshBase::elem_types() const
{
  std::vector<ElemType> et;

  assert (n_elem());
	
  const_elem_iterator       el (this->elements_begin());
  const const_elem_iterator end(this->elements_end());
  
  /**
   * Automatically get the first type
   */
  et.push_back((*el)->type());  ++el;

  /**
   * Loop over the rest of the elements.
   * If the current element type isn't in the
   * vector, insert it.
   */
  for (; el != end; ++el)
    {
      if (!std::count(et.begin(), et.end(), (*el)->type()))
	{
	  et.push_back((*el)->type());
	}
    }
  
  return et;
}



unsigned int MeshBase::n_elem_of_type(const ElemType type) const
{
  unsigned int cnt=0;

  const_type_elem_iterator       el (this->elements_begin(), type);
  const const_type_elem_iterator end(this->elements_end(),   type);

  for (; el!=end; ++el)
    cnt++;
  
  return cnt;
}



unsigned int MeshBase::n_active_elem_of_type(const ElemType type) const
{
  unsigned int cnt=0;

  const_active_type_elem_iterator       el (this->elements_begin(), type);
  const const_active_type_elem_iterator end(this->elements_end(),   type);
  
  for (; el!=end; ++el)
    cnt++;
    
  return cnt;
}



unsigned int MeshBase::total_weight() const
{
  unsigned int weight=0;

  const_elem_iterator       el (this->elements_begin());
  const const_elem_iterator end(this->elements_end());

  for ( ; el != end; ++el)
    weight += (*el)->n_nodes();
  
  return weight;
}



std::string MeshBase::get_info() const
{
  std::ostringstream out;

  out << " Mesh Information:"     << std::endl
      << "  mesh_dimension()="    << this->mesh_dimension()    << std::endl
      << "  spatial_dimension()=" << this->spatial_dimension() << std::endl
      << "  n_nodes()="           << this->n_nodes()           << std::endl
      << "  n_elem()="            << this->n_elem()            << std::endl
#ifdef ENABLE_AMR
      << "   n_active_elem()="    << this->n_active_elem()     << std::endl
#endif
      << "  n_subdomains()="      << this->n_subdomains()      << std::endl
      << "  n_processors()="      << this->n_processors()      << std::endl
      << "  processor_id()="      << this->processor_id()      << std::endl;

  return out.str();
}


void MeshBase::print_info() const
{
  std::cout << this->get_info()
	    << std::endl;
}



void MeshBase::skip_comment_lines (std::istream &in,
				   const char comment_start)
{    
  char c;
  while (in.get(c), c==comment_start) 
    {
      char line[256];
      in.get (line, 255, '\n'); // ignore rest of line, at most 256 chars
      in.get (c);               // ignore '\n' at end of line.
    }
  
  // put back first character of
  // first non-comment line
  in.putback (c);
}



void MeshBase::find_neighbors()
{
  assert(this->n_nodes() != 0);
  assert(this->n_elem()  != 0);

  
  if (_dim == 1)
    error();


  START_LOG("find_neighbors()", "MeshBase");
  
  //TODO [BSK]: This should be removed later?!
  for (unsigned int e=0; e<this->n_elem(); e++)
    for (unsigned int s=0; s<this->elem(e)->n_neighbors(); s++)
      this->elem(e)->set_neighbor(s,NULL);

  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  {
    // data structures -- Use the hash_multimap if available
    typedef unsigned int                    key_type;
    typedef std::pair<Elem*, unsigned char> val_type;
    typedef std::pair<key_type, val_type>   key_val_pair;
    
#if   defined(HAVE_HASH_MAP)    
    typedef std::hash_multimap<key_type, val_type>       map_type;    
#elif defined(HAVE_EXT_HASH_MAP)
# if  __GNUC__ >= 3
    typedef __gnu_cxx::hash_multimap<key_type, val_type> map_type;
# else
    DIE A HORRIBLE DEATH
# endif
#else
    typedef std::multimap<key_type, val_type>            map_type;
#endif
    
    // A map from side keys to corresponding elements & side numbers  
    map_type side_to_elem_map;
  
    elem_iterator       el (this->elements_begin());
    const elem_iterator end(this->elements_end());
    
    for (; el != end; ++el)
      {
	Elem* element = *el;
	
	for (unsigned int ms=0; ms<element->n_neighbors(); ms++)
	  {
	  next_side:
	    
	    if (element->neighbor(ms) == NULL)
	      {
		// Get the key for the side of this element
		const unsigned int key = element->key(ms);
		
		// Look for elements that have an identical side key
		std::pair <map_type::iterator, map_type::iterator>
		  bounds = side_to_elem_map.equal_range(key);
		
		// May be multiple keys, check all the possible
		// elements which _might_ be neighbors.
		if (bounds.first != bounds.second)
		  {
		    // Get the side for this element
		    const AutoPtr<Elem> my_side(element->side(ms));

		    // Look at all the entries with an equivalent key
		    while (bounds.first != bounds.second)
		      {
			// Get the potential element
			Elem* neighbor = bounds.first->second.first;
			
			// Get the side for the neighboring element
			const unsigned int ns = bounds.first->second.second;
			const AutoPtr<Elem> their_side(neighbor->side(ns));
			
			// If found a match wbounds.firsth my side
			if (*my_side == *their_side) 
			  {
			    // So we are neighbors. Tell the other element 
			    element->set_neighbor (ms,neighbor);
			    neighbor->set_neighbor(ns,element);
			    
			    side_to_elem_map.erase (bounds.first);
			    
			    // get out of this nested crap
			    goto next_side; 
			  }

			++bounds.first;
		      }
		  }
		    
		// didn't find a match...
		// Build the map entry for this element
		key_val_pair kvp;
		
		kvp.first         = key;
		kvp.second.first  = element;
		kvp.second.second = ms;
		
		// use the lower bound as a hint for
		// where to put it.
#if defined(HAVE_HASH_MAP) || defined(HAVE_EXT_HASH_MAP)
		side_to_elem_map.insert (kvp);
#else
		side_to_elem_map.insert (bounds.first,kvp);
#endif
	      }
	  }
      }
  }

  
  
#ifdef ENABLE_AMR

  /**
   * Here we look at all of the child elements.
   * If a child element has a NULL neighbor it
   * is either because it is on the boundary
   * or because its neighbor is at a different
   * level.  In the latter case we must get the
   * neighbor from the parent.
   *
   * Furthermore, that neighbor better be active,
   * otherwise we missed a child somewhere.
   */
  not_level_elem_iterator el (this->elements_begin(), 0);
  not_level_elem_iterator end(this->elements_end(),   0);

  for (; el != end; ++el)
    {
      Elem* elem = *el;
      
      assert (elem->parent() != NULL);
      
      for (unsigned int s=0; s < elem->n_neighbors(); s++)
	if (elem->neighbor(s) == NULL)
	  {	    
	    elem->set_neighbor(s, elem->parent()->neighbor(s));
	    
#ifdef DEBUG	    
	    if (elem->neighbor(s) != NULL)
	      if (!elem->neighbor(s)->active())
		{
		  std::cerr << "I'm confused..." << std::endl;
		  this->write_gmv("bad_mesh.gmv");
		  error();
		}
#endif
	  }
    }
  
#endif

  STOP_LOG("find_neighbors()", "MeshBase");
}



#ifdef ENABLE_INFINITE_ELEMENTS


const Point MeshBase::build_inf_elem(bool be_verbose)
{
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


  this->find_neighbors();	// update elem->neighbor() tables


  START_LOG("build_inf_elem()", "MeshBase");

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


  STOP_LOG("build_inf_elem()", "MeshBase");

}


#endif // ifdef ENABLE_INFINITE_ELEMENTS




void MeshBase::build_nodes_to_elem_map (std::vector<std::vector<unsigned int> >&
					nodes_to_elem_map) const
{
  nodes_to_elem_map.resize (this->n_nodes());

  const_elem_iterator       el (this->elements_begin());
  const const_elem_iterator end(this->elements_end());

  for (; el != end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      {
	assert ((*el)->node(n) < nodes_to_elem_map.size());
	assert ((*el)->id()    < this->n_elem());
	
	nodes_to_elem_map[(*el)->node(n)].push_back((*el)->id());
      }
}



void MeshBase::all_tri ()
{
  assert (this->mesh_dimension() == 2);
	  
  std::vector<Elem*> new_elements;
  new_elements.reserve (2*this->n_active_elem());

  active_elem_iterator el (this->elements_begin());
  active_elem_iterator end(this->elements_end());

  for (; el!=end; ++el)
    {
      if ((*el)->type() == QUAD4)
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
	  
	  delete *el; //_elements[e];
	}
      
      else if ((*el)->type() == QUAD8)
	{
	  Elem* tri0 = new Tri6;
	  Elem* tri1 = new Tri6;
	  
	  Node* new_node = add_point((node((*el)->node(0)) +
				      node((*el)->node(1)) +
				      node((*el)->node(2)) +
				      node((*el)->node(3)))*.25
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
	  
	  delete *el; //_elements[e];
	}
      
      else if ((*el)->type() == QUAD9)
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

	  delete *el; //_elements[e];
	}
      else
	new_elements.push_back(*el);
    }
  
  _elements = new_elements;

  this->prepare_for_use();
}



void MeshBase::partition (const unsigned int n_sbdmns)
{
  assert (n_sbdmns > 0);

  MetisPartitioner mp (*this);

  mp.partition(n_sbdmns);
}



void MeshBase::renumber_nodes_and_elements ()
{
  START_LOG("renumber_nodes_and_elem()", "MeshBase");

  std::vector<Elem*> new_elem;
  std::vector<Node*> new_nodes;

  // Reserve space in the new containers.
  new_elem.reserve  (_elements.size());
  new_nodes.reserve (_nodes.size());
  
  // Begin by setting all node and element ids
  // to an invalid value.
  {
    for (unsigned int e=0; e<_elements.size(); e++)
      if (_elements[e] != NULL)
	_elements[e]->invalidate_id();

    for (unsigned int n=0; n<_nodes.size(); n++)
      if (_nodes[n] != NULL)	
	{
	  _nodes[n]->invalidate_id();
	  _nodes[n]->invalidate_processor_id();
	}
  }


  // Renumber the elements and the nodes.
  { 
    unsigned int next_free_elem = 0;
    unsigned int next_free_node = 0;

    
    // If there is only one processor, loop over the nodes and elements
    // and set their ids based on where they lie in their vectors.
    //if (this->n_processors() == 1)
    if (true)
      {
	elem_iterator       el    (this->elements_begin());
	const elem_iterator end_el(this->elements_end());
	
	for (; el != end_el; ++el)
	  {
	    // this element should _not_ have been numbered already
	    assert ((*el)->id() == Elem::invalid_id);
	    
	    (*el)->set_id(next_free_elem++);

	    // Add the element to the new list
	    new_elem.push_back(*el);
	  }


	node_iterator       nd    (this->nodes_begin());
	const node_iterator end_nd(this->nodes_end());

	//TODO:[BSK] This will not produce any inactive nodes to be deleted in the last step.  Might want to change this later.
	for (; nd != end_nd; ++nd)
	  {
	    // this node should _not_ have been numbered already
	    assert ((*nd)->id() == Node::invalid_id);
	    
	    (*nd)->set_id(next_free_node++);

	    //assert (this->processor_id() == 0);
	    
	    (*nd)->set_processor_id(0);
	     
	    // Add the node to the new list
	    new_nodes.push_back(*nd);
	  }	       
      }
	  

  
//     // Otherwise renumber the elements to be in contiguous blocks
//     // on the processors.  The nodes are numbered in the order they
//     // are encountered on the element.
//     else
//       {    
// 	for (unsigned int proc_id=0;
// 	     proc_id<this->n_processors(); proc_id++)
// 	  {
// 	    // Loop over the elements on the processor proc_id
// 	    pid_elem_iterator       el    (this->elements_begin(), proc_id);
// 	    const pid_elem_iterator end_el(this->elements_end(),   proc_id);
	    
// 	    for (; el != end_el; ++el)
// 	      {
// 		// this element should _not_ have been numbered already
// 		assert ((*el)->id() == Elem::invalid_id);
		
// 		(*el)->set_id(next_free_elem++);
		
// 		// Add the element to the new list
// 		new_elem.push_back(*el);
		
// 		// Number the nodes on the element that
// 		// have not been numbered already.
// 		for (unsigned int n=0; n<(*el)->n_nodes(); n++)
// 		  if ((*el)->get_node(n)->id() == Node::invalid_id)
// 		    {
// 		      (*el)->get_node(n)->set_id(next_free_node++);
		      
// 		      // Add the node to the new list
// 		      new_nodes.push_back((*el)->get_node(n));
		      
// 		      // How could this fail?  Only if something is WRONG!
// 		      assert ((*el)->get_node(n)->processor_id() ==
// 			      Node::invalid_processor_id);
		      
// 		      (*el)->get_node(n)->set_processor_id((*el)->processor_id());
// 		    }
// 	      }
// 	  }
//       }

    // This could only fail if we did something seriously WRONG!
    assert (new_elem.size()  == next_free_elem);
    assert (new_nodes.size() == next_free_node);
  }


//   // Delete the inactive nodes
//   for (unsigned int n=0; n<_nodes.size(); n++)
//     if (_nodes[n] != NULL)
//       if (!_nodes[n]->active())
// 	{
// 	  delete _nodes[n];
// 	  _nodes[n] = NULL;
// 	}
    
  // Finally, reassign the _nodes and _elem vectors
  _elements = new_elem;
  _nodes    = new_nodes;
  
  STOP_LOG("renumber_nodes_and_elem()", "MeshBase");
}




void MeshBase::find_boundary_nodes(std::vector<bool>& on_boundary)
{
  // Resize the vector which holds boundary nodes and fill with zeros.
  on_boundary.resize(this->n_nodes());
  std::fill(on_boundary.begin(),
	    on_boundary.end(),
	    false);

  // Loop over elements, find those on boundary, and
  // mark them as 1 in on_boundary.
  active_elem_iterator       el (this->elements_begin());
  const active_elem_iterator end(this->elements_end());
  
  for (; el != end; ++el)
    for (unsigned int s=0; s<(*el)->n_neighbors(); s++)
      if ((*el)->neighbor(s) == NULL) // on the boundary
	{
	  const AutoPtr<Elem> side((*el)->build_side(s));
	  
	  for (unsigned int n=0; n<side->n_nodes(); n++)
	    on_boundary[side->node(n)] = true;
	}
}




void MeshBase::distort (const Real factor,
			const bool perturb_boundary)
{
  assert (mesh_dimension() != 1);
  assert (n_nodes());
  assert (n_elem());
  assert ((factor >= 0.) && (factor <= 1.));

  START_LOG("distort()", "MeshBase");



  // First find nodes on the boundary and flag them
  // so that we don't move them
  // on_boundary holds 0's (not on boundary) and 1's (on boundary)
  std::vector<bool> on_boundary(this->n_nodes());
  
  if (!perturb_boundary)
    this->find_boundary_nodes(on_boundary);

  // Now calculate the minimum distance to
  // neighboring nodes for each node.
  // hmin holds these distances.
  std::vector<float> hmin(this->n_nodes(), 1.e20);
  
  active_elem_iterator       el (this->elements_begin());
  const active_elem_iterator end(this->elements_end());
  
  for (; el!=end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      hmin[(*el)->node(n)] = std::min(hmin[(*el)->node(n)],
				      static_cast<float>((*el)->hmin()));		

  
  // Now actually move the nodes
  {
    const unsigned int seed = 123456;
    
    // seed the random number generator
    srand(seed);
    
    for (unsigned int n=0; this->n_nodes(); n++)
      if (!on_boundary[n])
	{
	  // the direction, random but unit normalized
	  
	  Point dir( ((Real) rand())/((Real) RAND_MAX),
		     ((Real) rand())/((Real) RAND_MAX),
		     ((this->mesh_dimension() == 3) ?
		      ((Real) rand())/((Real) RAND_MAX) :
		      0.)
		     );
	  
	  dir(0) = (dir(0)-.5)*2.;
	  dir(1) = (dir(1)-.5)*2.;

	  if (this->mesh_dimension() == 3)
	    dir(2) = (dir(2)-.5)*2.;
	  
	  dir = dir.unit();

	  // if hmin[n]=1.e20 then the node is not
	  // used by any element.  We should not
	  // move it.
	  if (hmin[n] != 1.e20)
	    {
	      this->node(n)(0) += dir(0)*factor*hmin[n];
	      this->node(n)(1) += dir(1)*factor*hmin[n];
	      
	      if (this->mesh_dimension() == 3)
		this->node(n)(2) += dir(2)*factor*hmin[n];
	    }
	}
  }


  // All done  
  STOP_LOG("distort()", "MeshBase");
}



void MeshBase::translate (const Real xt,
			  const Real yt,
			  const Real zt)
{
  const Point p(xt, yt, zt);

  for (unsigned int n=0; n<this->n_nodes(); n++)
    this->node(n) += p;
}



void MeshBase::rotate (const Real,
		       const Real,
		       const Real)
{
  error();
}



void MeshBase::scale (const Real xs,
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
  for (unsigned int n=0; n<this->n_nodes(); n++)
    this->node(n)(0) *= x_scale;


  // Only scale the y coordinate in 2 and 3D
  if (this->spatial_dimension() > 1)
    {

      for (unsigned int n=0; n<this->n_nodes(); n++)
	this->node(n)(1) *= y_scale;

      // Only scale the z coordinate in 3D
      if (this->spatial_dimension() == 3)
	{
	  for (unsigned int n=0; n<this->n_nodes(); n++)
	    this->node(n)(2) *= z_scale;
	}
    }
}



std::pair<Point, Point> 
MeshBase::bounding_box() const
{
  // processor bounding box with no arguments
  // computes the global bounding box
  return this->processor_bounding_box();
}



Sphere
MeshBase::bounding_sphere() const
{
  std::pair<Point, Point> bbox = this->bounding_box();

  const Real  diag = (bbox.second - bbox.first).size();
  const Point cent = (bbox.second + bbox.first)/2.;

  return Sphere (cent, .5*diag);
}



std::pair<Point, Point> 
MeshBase::processor_bounding_box (const unsigned int pid) const
{
  assert (this->n_nodes() != 0);

  Point min(1.e30,   1.e30,  1.e30);
  Point max(-1.e30, -1.e30, -1.e30);

  // By default no processor is specified and we compute
  // the bounding box for the whole domain.
  if (pid == static_cast<unsigned int>(-1))
    {
      for (unsigned int n=0; n<this->n_nodes(); n++)
	for (unsigned int i=0; i<this->spatial_dimension(); i++)
	  {
	    min(i) = std::min(min(i), this->point(n)(i));
	    max(i) = std::max(max(i), this->point(n)(i));
	  }      
    }
  // if a specific processor id is specified then we need
  // to only consider those elements living on that processor
  else
    {
      const_pid_elem_iterator       el (this->elements_begin(), pid);
      const const_pid_elem_iterator end(this->elements_end(),   pid);

      for (; el != end; ++el)
	for (unsigned int n=0; n<(*el)->n_nodes(); n++)
	    for (unsigned int i=0; i<this->spatial_dimension(); i++)
	      {
		min(i) = std::min(min(i), this->point((*el)->node(n))(i));
		max(i) = std::max(max(i), this->point((*el)->node(n))(i));
	      }      
    }

  const std::pair<Point, Point> ret_val(min, max);

  return ret_val;  
}



Sphere
MeshBase::processor_bounding_sphere (const unsigned int pid) const
{
  std::pair<Point, Point> bbox = this->processor_bounding_box(pid);

  const Real  diag = (bbox.second - bbox.first).size();
  const Point cent = (bbox.second + bbox.first)/2.;

  return Sphere (cent, .5*diag);
}



std::pair<Point, Point> 
MeshBase::subdomain_bounding_box (const unsigned int sid) const
{
  assert (this->n_nodes() != 0);

  Point min( 1.e30,  1.e30,  1.e30);
  Point max(-1.e30, -1.e30, -1.e30);

  // By default no subdomain is specified and we compute
  // the bounding box for the whole domain.
  if (sid == static_cast<unsigned int>(-1))
    {
      for (unsigned int n=0; n<this->n_nodes(); n++)
	for (unsigned int i=0; i<this->spatial_dimension(); i++)
	  {
	    min(i) = std::min(min(i), this->point(n)(i));
	    max(i) = std::max(max(i), this->point(n)(i));
	  }      
    }

  // if a specific subdomain id is specified then we need
  // to only consider those elements living on that subdomain
  else
    {
      for (unsigned int e=0; e<this->n_elem(); e++)
	if (this->elem(e)->subdomain_id() == sid)
	  for (unsigned int n=0; n<this->elem(e)->n_nodes(); n++)
	    for (unsigned int i=0; i<this->spatial_dimension(); i++)
	      {
		min(i) = std::min(min(i), this->point(this->elem(e)->node(n))(i));
		max(i) = std::max(max(i), this->point(this->elem(e)->node(n))(i));
	      }      
    }

  const std::pair<Point, Point> ret_val(min, max);

  return ret_val;  
}



Sphere MeshBase::subdomain_bounding_sphere (const unsigned int sid) const
{
  std::pair<Point, Point> bbox = this->subdomain_bounding_box(sid);

  const Real  diag = (bbox.second - bbox.first).size();
  const Point cent = (bbox.second + bbox.first)/2.;

  return Sphere (cent, .5*diag);
}



void MeshBase::build_L_graph (PetscMatrix<Number>& conn) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: This fuctionality requires PETSC support!"
	    << std::endl;
  error();

#else

  // Initialize the connectivity matrix.
  {
    conn.init(this->n_nodes(),
	      this->n_nodes(),
	      this->n_nodes(),
	      this->n_nodes());

    // be sure the diagonals are all 0 so that ++ works.
    for (unsigned int n=0; n<this->n_nodes(); n++)
      conn.set(n,n,0.);

  }
  
  switch (this->mesh_dimension())
    {
    case 1:
      {
	std::cerr << "ERROR:  The connectivity graph doesn't make much sense"
		  << std::endl
		  << " in 1D!"
		  << std::endl;
	error();
      }

      
      // Create the graph for a 2D mesh.  Do this by looking
      // at element edges.
    case 2:
      {
	const_active_elem_iterator       el (this->elements_begin());
	const const_active_elem_iterator end(this->elements_end());
	
	for (; el != end; ++el)
	  for (unsigned int s=0; s<(*el)->n_neighbors(); s++)
	    if (((*el)->neighbor(s) == NULL) ||
		((*el)->id() > (*el)->neighbor(s)->id()))
	      {
		AutoPtr<Elem> side((*el)->build_side(s));
		
		const unsigned int n0 = side->node(0);
		const unsigned int n1 = side->node(1);
		
		conn.set(n0,n0, conn(n0,n0) + 1.);		  
		conn.set(n0,n1, -1.);
		
		conn.set(n1,n1, conn(n1,n1) + 1.);
		conn.set(n1,n0, -1.);
	      }
	
	// All done.
	break;
      }



      // Create the graph for a 3D mesh.  Do this by looking
      // at element faces, then at the edges of the face.
    case 3:
      {
	const_active_elem_iterator       el (this->elements_begin());
	const const_active_elem_iterator end(this->elements_end());

	for (; el != end; ++el)
	  for (unsigned int f=0; f<(*el)->n_neighbors(); f++) // Loop over faces
	    if (((*el)->neighbor(f) == NULL) ||
		((*el)->id() > (*el)->neighbor(f)->id()))
	      {
		AutoPtr<Elem> face((*el)->build_side(f));
		
		for (unsigned int s=0; s<face->n_neighbors(); s++) // Loop over face's edges
		  {
		    AutoPtr<Elem> side(face->build_side(s));
		    
		    const unsigned int n0 = side->node(0);
		    const unsigned int n1 = side->node(1);
		    
		    // If this is the first time we've seen this edge
		    if (conn(n0,n1) == 0.)
		      {
			assert (conn(n1,n0) == 0.);
			
			conn.set(n0,n0, conn(n0,n0) + 1.);		  
			conn.set(n0,n1, -1.);
			
			conn.set(n1,n1, conn(n1,n1) + 1.);
			conn.set(n1,n0, -1.);
		      }
		  }
	      }

	// All done
	break;
      }
      

    default:
      // what?
      error();
    }


  // OK, now the matrix is built.  Close it
  // and return.
  conn.close();
  
  return;

#endif
}



void MeshBase::build_script_L_graph (PetscMatrix<Number>& conn) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: This fuctionality requires PETSC support!"
	    << std::endl;
  error();

#else

  // Inefficient at the moment.  We build an L
  // matrix and use it to create the script L matrix
  PetscMatrix<Number> l_conn;

  build_L_graph (l_conn);

  conn.init(l_conn.m(),
	    l_conn.n(),
	    l_conn.m(),
	    l_conn.n());

  
  switch (mesh_dimension())
    {
    case 1:
      {
	std::cerr << "ERROR:  The connectivity graph doesn't make much sense"
		  << std::endl
		  << " in 1D!"
		  << std::endl;
	error();
      }

      
      // Create the graph for a 2D mesh.  Do this by looking
      // at element edges.
    case 2:
      {
	const_active_elem_iterator       el (this->elements_begin());
	const const_active_elem_iterator end(this->elements_end());
	
	for (; el != end; ++el)
	  for (unsigned int s=0; s<(*el)->n_neighbors(); s++)
	    if (((*el)->neighbor(s) == NULL) ||
	        ((*el) > (*el)->neighbor(s)))
	      {
		AutoPtr<Elem> side((*el)->build_side(s));
		
		const unsigned int n0 = side->node(0);
		const unsigned int n1 = side->node(1);
		
		conn.set(n0,n0, 1.);
		conn.set(n1,n1, 1.);
		
#ifdef USE_COMPLEX_NUMBERS
		const Real prod_term = -1./sqrt(l_conn(n0,n0).real()*l_conn(n1,n1).real());
#else
		const Real prod_term = -1./sqrt(l_conn(n0,n0)*l_conn(n1,n1));
#endif
		
		conn.set(n0,n1, prod_term);
		conn.set(n1,n0, prod_term);
	      }

	// All done.
	break;
      }



      // Create the graph for a 3D mesh.  Do this by looking
      // at element faces, then at the edges of the face.
    case 3:
      {
	const_active_elem_iterator       el (this->elements_begin());
	const const_active_elem_iterator end(this->elements_end());
	
	for (; el != end; ++el)
	  for (unsigned int f=0; f<(*el)->n_neighbors(); f++) // Loop over faces
	    if (((*el)->neighbor(f) == NULL) ||
		((*el) > (*el)->neighbor(f)))
	      {
		AutoPtr<Elem> face((*el)->build_side(f));
		
		for (unsigned int s=0; s<face->n_neighbors(); s++) // Loop over face's edges
		  {
		    AutoPtr<Elem> side(face->build_side(s));
		    
		    const unsigned int n0 = side->node(0);
		    const unsigned int n1 = side->node(1);
		    
		    conn.set(n0,n0, 1.);
		    conn.set(n1,n1, 1.);
		    
#ifdef USE_COMPLEX_NUMBERS
		    const Real prod_term = -1./sqrt(l_conn(n0,n0).real()*l_conn(n1,n1).real());
#else
		    const Real prod_term = -1./sqrt(l_conn(n0,n0)*l_conn(n1,n1));
#endif
		    
		    conn.set(n0,n1, prod_term);
		    conn.set(n1,n0, prod_term);
		  }
	      }
	
	// All done
	break;
      }
      

    default:
      // what?
      error();
    }


  // OK, now the matrix is built.  Close it
  // and return.
  conn.close();
  
  return;  

#endif
}



void MeshBase::read(const std::string&)
{
  std::cerr << "ERROR:  You shouldn't be calling this" << std::endl
	    << " Use Mesh::read() instead." << std::endl;
  error();
}



void MeshBase::write(const std::string& name)
{
  START_LOG("write()", "MeshBase");
  
  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      this->write_tecplot (name);
    
    else if (name.rfind(".plt") < name.size())
      this->write_tecplot_binary (name);

    else if (name.rfind(".ucd") < name.size())
      this->write_ucd (name);

    else if (name.rfind(".gmv") < name.size())
      {
	if (n_subdomains() > 1)
	  this->write_gmv_binary(name, NULL, NULL, true);
	else
	  this->write_gmv_binary(name);
      }
  }

  STOP_LOG("write()", "MeshBase");
}



void MeshBase::write(const std::string& name,
		     std::vector<Number>& v,
		     std::vector<std::string>& vn)
{
  START_LOG("write()", "MeshBase");
  
  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      this->write_tecplot (name, &v, &vn);
    
    else if (name.rfind(".plt") < name.size())
      this->write_tecplot_binary (name, &v, &vn);
    
    else if (name.rfind(".gmv") < name.size())
      {
	if (this->n_subdomains() > 1)
	  this->write_gmv_binary(name, &v, &vn, true);
	else
	  this->write_gmv_binary(name, &v, &vn);
      }
  }

  STOP_LOG("write()", "MeshBase");
}



#ifdef USE_COMPLEX_NUMBERS

const char* MeshBase::complex_filename(const std::string& _n,
				       unsigned int r_o_c) const
{
  std::string loc=_n;
  
  if (r_o_c == 0)
    loc.append(".real");
  
  else
    loc.append(".imag");
  
  return loc.c_str();
}



void MeshBase::prepare_complex_data(const std::vector<Number>* source,
				    std::vector<Real>* real_part,
				    std::vector<Real>* imag_part) const
{
  real_part->resize(source->size());
  imag_part->resize(source->size());


  for (unsigned int i=0; i<source->size(); i++)
    {
      (*real_part)[i] = (*source)[i].real();
      (*imag_part)[i] = (*source)[i].imag();
    }
}

#endif // #ifdef USE_COMPLEX_NUMBERS



#ifdef ENABLE_AMR

void MeshBase::trim_unused_elements (std::set<unsigned int>& unused_elements)
{
  /**
   * Anything we clear in this routiune
   * will invalidate the unknowing boundary
   * mesh, so we need to clear it.  It must
   * be recreated before reuse.  
   */
  //boundary_info.boundary_mesh.clear();
  
  
  /**
   * Trim the unused elements
   */
  {
    // We don't Really need this in the
    // current implementation
    unused_elements.clear();

    // for the time being we make a copy
    // of the elements vector since the pointers
    // are relatively small.  Note that this is
    // not _necessary_, but it should be
    // less expensive than repeated calls
    // to std::vector<>::erase()    
    std::vector<Elem*> new_elements;
    
    new_elements.resize(this->n_elem());

    unsigned int ne=0;
    
    for (unsigned int e=0; e<this->n_elem(); e++)
      if (_elements[e] != NULL)
	new_elements[ne++] = _elements[e]; 

    new_elements.resize(ne);
    
    _elements = new_elements;

    /**
     * Explicitly check that it worked
     * in DEBUG mode.
     */
#ifdef DEBUG

    for (unsigned int e=0; e<this->n_elem(); e++)
      assert (_elements[e] != NULL);
    
#endif
    
  }
}

#endif // #ifdef ENABLE_AMR




