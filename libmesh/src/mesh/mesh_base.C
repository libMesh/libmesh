// $Id: mesh_base.C,v 1.5 2003-01-24 17:24:44 jwpeterson Exp $

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



// C++ includes
#include <algorithm>
#include <sstream>
#include <math.h>
#include <set>

// Local includes
#include "mesh_base.h"
#include "face_tri3.h"
#include "face_tri6.h"


#ifdef HAVE_SFCURVES
// prototype for SFC code
namespace sfc {
  extern "C" {
#include "sfcurves.h"
  }
}
#endif





// ------------------------------------------------------------
// MeshBase class member functions
MeshBase::MeshBase (unsigned int d,
		    unsigned int pid) :
  _n_sbd(1),
  _n_proc(1),
  _dim(d),
  _proc_id(pid),
#ifdef ENABLE_PERFORMANCE_LOGGING
  _perf_log("MeshBase", true)
#else
  _perf_log("MeshBase", false)
#endif
  
{
  assert (DIM <= 3);
  assert (DIM >= _dim);
};



MeshBase::MeshBase (const MeshBase& other_mesh) :
  _n_sbd(other_mesh._n_sbd),
  _n_proc(other_mesh._n_proc),
  _dim(other_mesh._dim),
  _proc_id(other_mesh._proc_id)

{

  _nodes = other_mesh._nodes;
  _elements = other_mesh._elements;
  
};


MeshBase::~MeshBase()
{
  clear();
};



Node* MeshBase::add_point(const Point& p,
			  const unsigned int num)
{  
  _perf_log.start_event("add_point()");

  if ((num == static_cast<unsigned int>(-1)) ||
      (num == n_nodes()))
    {
      _nodes.push_back(Node::build(p, n_nodes()));

      _perf_log.stop_event("add_point()");
  
      return node_ptr(n_nodes()-1);
    }
  else
    {
      assert (num < n_nodes() );
      assert (node_ptr(num) != NULL);
      assert (node_ptr(num)->id() != Node::invalid_id);
      
      node(num)          = p;
      node(num).set_id() = num;
      
      _perf_log.stop_event("add_point()");
  
      return node_ptr(num);
    }

  
  // We'll never get here...
  error();
  return NULL;
};



void MeshBase::add_elem(Elem* e, const unsigned int n)
{
  _perf_log.start_event("add_elem()");

  if ((n == static_cast<unsigned int>(-1)) ||
      (n == n_elem()))
    _elements.push_back(e);
  else
    {
      assert( n < n_elem() );

      _elements[n] = e;
    }

  _perf_log.stop_event("add_elem()");
};



unsigned int MeshBase::n_active_elem() const
{
  unsigned int num=0;

  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->active())
      num++;

  return num;
};



void MeshBase::clear()
{
  // Reset the number of subdomains and the
  // number of processors
  _n_sbd  = 1;
  _n_proc = 1;

  // Clear the elements data structure
  {
    for (unsigned int i=0; i<n_elem(); i++)
      if (elem(i) != NULL)
	delete elem(i);
    
    _elements.clear();
  };

  // clear the nodes data structure
  {
    for (unsigned int n=0; n<n_nodes(); n++)
      if (node_ptr(n) != NULL)
	delete node_ptr(n);
    
    _nodes.clear();
  };  
};



unsigned int MeshBase::n_sub_elem() const
{
  unsigned int ne=0;

  for (unsigned int i=0; i<n_elem(); i++)
    ne += elem(i)->n_sub_elem();

  return ne;
};



unsigned int MeshBase::n_active_sub_elem() const
{
  unsigned int ne=0;

  for (unsigned int i=0; i<n_elem(); i++)
    if (elem(i)->active())
      ne += elem(i)->n_sub_elem();

  return ne;
};



std::vector<ElemType> MeshBase::elem_types() const
{
  std::vector<ElemType> et;

  assert (n_elem());
	
  /**
   * Automatically get the first type
   */
  et.push_back(elem(0)->type()); 

  /**
   * Loop over the rest of the elements.
   * If the current element type isn't in the
   * vector, insert it.
   */
  for (unsigned int e=1; e<n_elem(); e++)
    {
      if (!std::count(et.begin(), et.end(), elem(e)->type()))
	{
	  et.push_back(elem(e)->type());
	}
    }
  
  return et;
};



unsigned int MeshBase::n_elem_of_type(const ElemType type) const
{
  unsigned int cnt=0;

  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->type() == type)
      cnt++;
  
  return cnt;
};



unsigned int MeshBase::n_active_elem_of_type(const ElemType type) const
{
  unsigned int cnt=0;

  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->active())
      if (elem(e)->type() == type)
	cnt++;
  
  return cnt;
};



unsigned int MeshBase::total_weight() const
{
  unsigned int weight=0;

  for (unsigned int e=0; e<n_elem(); e++)
    weight += elem(e)->n_nodes();
  
  return weight;
};



std::string MeshBase::get_info() const
{
  std::ostringstream out;

  out << " Mesh Information:" << std::endl
      << "  mesh_dimension()=" << mesh_dimension() << std::endl
      << "  spatial_dimension()=" << spatial_dimension() << std::endl
      << "  n_nodes()=" << n_nodes() << std::endl
      << "  n_elem()=" << n_elem() << std::endl
#ifdef ENABLE_AMR
      << "  n_active_elem()=" << n_active_elem() << std::endl
#endif
      << "  n_subdomains()=" << n_subdomains() << std::endl
      << "  n_processors()=" << n_processors() << std::endl
      << "  processor_id()=" << processor_id() << std::endl;

  return out.str();
};


void MeshBase::print_info() const
{
  std::cout << get_info()
	    << std::endl;
};



void MeshBase::skip_comment_lines (std::istream &in,
				   const char comment_start)
{    
  char c;
  while (in.get(c), c==comment_start) 
    {
      char line[256];
      in.get (line, 255, '\n'); // ignore rest of line, at most 256 chars
      in.get (c);               // ignore '\n' at end of line.
    };
  
  // put back first character of
  // first non-comment line
  in.putback (c);
};



void MeshBase::find_neighbors()
{
  assert(n_nodes() != 0);
  assert(n_elem() != 0);

  
  if (_dim == 1)
    error();


  _perf_log.start_event("find_neighbors()");
  
  // data structures
  typedef std::pair<unsigned int, Elem*> key_val_pair;
  std::multimap<unsigned int, Elem*>     side_to_elem;
  
  //TODO [BSK]: This should be removed later?!
  for (unsigned int e=0; e<n_elem(); e++)
    for (unsigned int s=0; s<elem(e)->n_sides(); s++)
      elem(e)->set_neighbor(s,NULL);

  unsigned int c0=0, c1=0;

  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  for (unsigned int e=0; e<n_elem(); e++)
    {
      Elem* element = elem(e);
      
      for (unsigned int ms=0; ms<element->n_sides(); ms++)
	{
	next_side:
	  
	  if (element->neighbor(ms) == NULL)
	    {
	      const AutoPtr<Elem> my_side(element->side(ms));
	      const unsigned int key      = my_side->key();

	      // Look for elements that have an identical side key
	      std::pair <std::multimap<unsigned int, Elem*>::iterator,
		std::multimap<unsigned int, Elem*>::iterator >
		bounds = side_to_elem.equal_range(key);

	      // If no side corresponding to the key was found...
	      if (bounds.first == bounds.second)
		{
		  c0++;
		  
		  // use the lower bound as a hint for
		  // where to put it.
		  side_to_elem.insert (bounds.first,
				       key_val_pair(key, element));	      
		}

	      // Otherwise may be multiple keys, check all the possible
	      // elements which _might_ be neighbors.  Be sure not to check
	      // the element of interest to avoid a false match!
	      else
		{
		  c1++;
		  for (std::multimap<unsigned int, Elem*>::iterator
			 it = bounds.first; it != bounds.second; ++it)
		    if (it->second != element)
		      {
			Elem* neighbor = it->second;
			
			// look at all their sides
			for (unsigned int ns=0; ns<neighbor->n_sides(); ns++) 
			  {
			    const AutoPtr<Elem> their_side(neighbor->side(ns));

			    // and find a match with my side
			    if (*my_side == *their_side) 
			      {
				// So we are neighbors. Tell the other
				// element to avoid duplicate searches.
				element->set_neighbor (ms,neighbor);
				neighbor->set_neighbor(ns,element);
				
				side_to_elem.erase (it);
				
				// get out of this nested crap
				goto next_side; 
			      };
			  };
		      };
		  
		  // didn't find a match...
		  side_to_elem.insert (bounds.first,
				       key_val_pair(key, element));
		};
	    };
	};
    };
  
  side_to_elem.clear();

  
  
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
  
  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->level() != 0)
      for (unsigned int s=0; s<elem(e)->n_sides(); s++)
	if (elem(e)->neighbor(s) == NULL)
	  {
	    elem(e)->set_neighbor(s,elem(e)->parent()->neighbor(s));

#ifdef DEBUG	    
	    if (elem(e)->neighbor(s) != NULL)
	      assert (elem(e)->neighbor(s)->active());
#endif
	  }
  
#endif

//   std::cout << "Counter 0=" << c0
// 	    << ", Counter 1=" << c1
// 	    << std::endl;
    
    

  _perf_log.stop_event("find_neighbors()");
};



#ifdef ENABLE_INFINITE_ELEMENTS


void MeshBase::build_inf_elem()
{
  // determine origin automatically,
  // works only if the mesh has no symmetry planes.
  std::pair<Point, Point> b_box = bounding_box();
    
  build_inf_elem( (b_box.first+b_box.second)/2. );	 /* untested */
}



void MeshBase::build_inf_elem(const Point& origin, 
			      const bool x_sym, 
			      const bool y_sym, 
			      const bool z_sym,
			      const bool be_verbose)
{
		
  if (be_verbose)
    std::cout << "Updating element neighbor tables..." << std::endl;
	
  find_neighbors();	// update elem->neighbor() tables


  _perf_log.start_event("build_inf_elem()");

  std::set< std::pair<unsigned int,unsigned int> > faces,ofaces;
  std::set< std::pair<unsigned int,unsigned int> > :: iterator face_it;
	
  std::set<unsigned int> onodes;
  std::set<unsigned int> :: iterator on_it;
	
  real max_r=0.;
  unsigned int max_r_node;

  if (be_verbose)
    {
      std::cout << "Collecting boundary sides";
      if (x_sym || y_sym || z_sym)
	std::cout << ", skipping sides in symmetry planes..." << std::endl;
      else
	std::cout << "..." << std::endl;
    }

  /**
   * Iterate through all elements and sides, collect indices of all active
   * boundary sides in the faces set. Skip sides which lie in symmetry planes.
   * Later, sides of the inner boundary will be sorted out.
   */

  for(unsigned int e=0;e<n_elem();e++)
    {
      if (!(elem(e)->active()))
	continue;

      for (unsigned int s=0; s<elem(e)->n_sides(); s++)
	{
	  if (elem(e)->neighbor(s) != NULL)
	    continue;	 // check if elem(e) is on the boundary
		
	  /* note that it is safe to use the Elem::side() method, 
	     which gives a non-full-ordered element */
	  AutoPtr<Elem> side(elem(e)->side(s));			
		
	  bool sym_side=false;		
			
	  bool on_x_sym=true;
	  bool on_y_sym=true;
	  bool on_z_sym=true;
			
	  // TODO:[HVDH] Find a better criterion based on mesh properties

		
	  for(unsigned int n=0;n<side->n_nodes();n++)
	    {
	
	      Point dist_from_origin=point(side->node(n))-origin;

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

	      real r=dist_from_origin.size();
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

    AutoPtr<Elem> side(elem(p.first)->side(p.second));
		
    bool found=false;
    for(unsigned int node=0;node<side->n_nodes();node++)
      if(onodes.count(side->node(node))){found=true;break;}
    	
    	
    /* If a new oface is found, include it's nodes in onodes */		
    	
    if(found)		
      {
	for(unsigned int n=0;n<side->n_nodes();n++)
	  onodes.insert(side->node(n));
    				
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
	
	
  if (be_verbose)
    std::cout << "Found " << faces.size() << " inner and " 
	      << ofaces.size() << " outer boundary faces." << std::endl;
	
	
  faces.clear();		//free memory


  // outer_nodes maps onodes to their duplicates

  std::map<unsigned int, Node *> outer_nodes;


  /**
   * for each boundary node, add an outer_node with 
   * double distance from origin.
   */

  for(on_it=onodes.begin();on_it!=onodes.end();on_it++)
    {
      Point p=Point(point(*on_it))*2-origin;
      outer_nodes[*on_it]=add_point(p);
    }


  // for verbose, remember n_elem
  unsigned int n_conventional_elem = n_elem();


  /**
   * build Elems based on boundary side type
   */

  for(face_it=ofaces.begin();face_it!=ofaces.end();face_it++)
    {
	
      std::pair<unsigned int,unsigned int> p=*face_it;
	
	
      AutoPtr<Elem> side(elem(p.first)->build_side(p.second));							
				
      // create cell depending on side type
				
      Elem* el;
      switch(side->n_nodes())
	{
	  // TRIs					
	case 3:	
	  el=Elem::build(INFPRISM6);
	  break;
					 		
	case 6: 
	  el=Elem::build(INFPRISM12);
	  break;
							
	  // QUADs					
	case 4: 
	  el=Elem::build(INFHEX8);
	  break;
							
	case 8: 
	  el=Elem::build(INFHEX16);
	  break;
							
	case 9: 
	  el=Elem::build(INFHEX18);
	  break;

	default: 
	  std::cout << "MeshBase::build_inf_elem(Point, bool, bool, bool, bool): invalid face element" 
		    << std::endl;
	  continue;
	}
			
      // assign nodes to cell
						
      for(unsigned int i=0;i<side->n_nodes();i++)
	{
	  el->set_node(i  )=side->get_node(i);
	  el->set_node(i+side->n_nodes())=outer_nodes[side->node(i)];
	}
			
      // add infinite element to mesh			
      add_elem(el);	
    }


  if (be_verbose)
    std::cout << "Added "
	      << n_elem()-n_conventional_elem
	      << " infinite elements to mesh."
	      << std::endl;


  _perf_log.stop_event("build_inf_elem()");

}


#endif // ifdef ENABLE_INFINITE_ELEMENTS




void MeshBase::build_nodes_to_elem_map (std::vector<std::vector<unsigned int> >&
					nodes_to_elem_map) const
{
  nodes_to_elem_map.resize (n_nodes());

  for (unsigned int e=0; e<n_elem(); e++)
    for (unsigned int n=0; n<elem(e)->n_nodes(); n++)
      nodes_to_elem_map[elem(e)->node(n)].push_back(e);
};



void MeshBase::all_tri ()
{
  assert (mesh_dimension() == 2);
	  
  std::vector<Elem*> new_elements;
  new_elements.reserve (2*n_elem());

  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->active())
      if (elem(e)->type() == QUAD4)
	{
	  Elem* tri0 = new Tri3;
	  Elem* tri1 = new Tri3;
	  
	  // Check for possible edge swap
	  if ((elem(e)->point(0) - elem(e)->point(2)).size() <
	      (elem(e)->point(1) - elem(e)->point(3)).size())
	    {	      
	      tri0->set_node(0) = elem(e)->get_node(0);
	      tri0->set_node(1) = elem(e)->get_node(1);
	      tri0->set_node(2) = elem(e)->get_node(2);
	      
	      tri1->set_node(0) = elem(e)->get_node(0);
	      tri1->set_node(1) = elem(e)->get_node(2);
	      tri1->set_node(2) = elem(e)->get_node(3);
	    }

	  else
	    {
	      tri0->set_node(0) = elem(e)->get_node(0);
	      tri0->set_node(1) = elem(e)->get_node(1);
	      tri0->set_node(2) = elem(e)->get_node(3);
	      
	      tri1->set_node(0) = elem(e)->get_node(1);
	      tri1->set_node(1) = elem(e)->get_node(2);
	      tri1->set_node(2) = elem(e)->get_node(3);
	    }
	  
	  new_elements.push_back(tri0);
	  new_elements.push_back(tri1);
	  
	  delete _elements[e];
	}
      
      else if (elem(e)->type() == QUAD8)
	{
	  Elem* tri0 = new Tri6;
	  Elem* tri1 = new Tri6;
	  
	  Node* new_node = add_point((node(elem(e)->node(0)) +
				      node(elem(e)->node(1)) +
				      node(elem(e)->node(2)) +
				      node(elem(e)->node(3)))*.25
				     );
	  
	  // Check for possible edge swap
	  if ((elem(e)->point(0) - elem(e)->point(2)).size() <
	      (elem(e)->point(1) - elem(e)->point(3)).size())
	    {	      
	      tri0->set_node(0) = elem(e)->get_node(0);
	      tri0->set_node(1) = elem(e)->get_node(1);
	      tri0->set_node(2) = elem(e)->get_node(2);
	      tri0->set_node(3) = elem(e)->get_node(4);
	      tri0->set_node(4) = elem(e)->get_node(5);
	      tri0->set_node(5) = new_node;
	      
	      tri1->set_node(0) = elem(e)->get_node(0);
	      tri1->set_node(1) = elem(e)->get_node(2);
	      tri1->set_node(2) = elem(e)->get_node(3);
	      tri1->set_node(3) = new_node;
	      tri1->set_node(4) = elem(e)->get_node(6);
	      tri1->set_node(5) = elem(e)->get_node(7);

	    }
	  
	  else
	    {
	      tri0->set_node(0) = elem(e)->get_node(3);
	      tri0->set_node(1) = elem(e)->get_node(0);
	      tri0->set_node(2) = elem(e)->get_node(1);
	      tri0->set_node(3) = elem(e)->get_node(7);
	      tri0->set_node(4) = elem(e)->get_node(4);
	      tri0->set_node(5) = new_node;
	      
	      tri1->set_node(0) = elem(e)->get_node(1);
	      tri1->set_node(1) = elem(e)->get_node(2);
	      tri1->set_node(2) = elem(e)->get_node(3);
	      tri1->set_node(3) = elem(e)->get_node(5);
	      tri1->set_node(4) = elem(e)->get_node(6);
	      tri1->set_node(5) = new_node;
	    }
	  
	  new_elements.push_back(tri0);
	  new_elements.push_back(tri1);
	  
	  delete _elements[e];
	}
      
      else if (elem(e)->type() == QUAD9)
	{
	  Elem* tri0 = new Tri6;
	  Elem* tri1 = new Tri6;

	  // Check for possible edge swap
	  if ((elem(e)->point(0) - elem(e)->point(2)).size() <
	      (elem(e)->point(1) - elem(e)->point(3)).size())
	    {	      
	      tri0->set_node(0) = elem(e)->get_node(0);
	      tri0->set_node(1) = elem(e)->get_node(1);
	      tri0->set_node(2) = elem(e)->get_node(2);
	      tri0->set_node(3) = elem(e)->get_node(4);
	      tri0->set_node(4) = elem(e)->get_node(5);
	      tri0->set_node(5) = elem(e)->get_node(8);
	      
	      tri1->set_node(0) = elem(e)->get_node(0);
	      tri1->set_node(1) = elem(e)->get_node(2);
	      tri1->set_node(2) = elem(e)->get_node(3);
	      tri1->set_node(3) = elem(e)->get_node(8);
	      tri1->set_node(4) = elem(e)->get_node(6);
	      tri1->set_node(5) = elem(e)->get_node(7);
	    }

	  else
	    {
	      tri0->set_node(0) = elem(e)->get_node(0);
	      tri0->set_node(1) = elem(e)->get_node(1);
	      tri0->set_node(2) = elem(e)->get_node(3);
	      tri0->set_node(3) = elem(e)->get_node(4);
	      tri0->set_node(4) = elem(e)->get_node(8);
	      tri0->set_node(5) = elem(e)->get_node(7);
	      
	      tri1->set_node(0) = elem(e)->get_node(1);
	      tri1->set_node(1) = elem(e)->get_node(2);
	      tri1->set_node(2) = elem(e)->get_node(3);
	      tri1->set_node(3) = elem(e)->get_node(5);
	      tri1->set_node(4) = elem(e)->get_node(6);
	      tri1->set_node(5) = elem(e)->get_node(8);
	    }
	  
	  new_elements.push_back(tri0);
	  new_elements.push_back(tri1);
	  
	  delete _elements[e];
	}
      else
	new_elements.push_back(elem(e));

  _elements = new_elements;

  find_neighbors();
};



void MeshBase::sfc_partition(const unsigned int n_sbdmns,
			     const std::string& type)
{
  // won't work without Bill's library!
#ifndef HAVE_SFCURVES

  {
    std::cerr << "ERROR:  Not compiled with space-filling curve" << std::endl
	      << " support.  Using linear partitioning instead" << std::endl
	      << " This partitioning could be arbitrarily bad!" << std::endl
	      <<  std::endl;
    
    set_n_subdomains() = n_sbdmns;
    set_n_processors() = n_sbdmns;

    // check for easy return
    if (n_sbdmns == 1)
      {
	for (unsigned int e=0; e<n_elem(); e++)
	  elem(e)->subdomain_id() = 
	    elem(e)->processor_id() = 0;
	
	return;
      };

    
    const unsigned int blksize = n_elem()/n_sbdmns; 
    
    for (unsigned int e=0; e<n_elem(); e++)
      elem(e)->subdomain_id() = 
	elem(e)->processor_id() = 
	(int) (e/blksize);
    
    return;
  }
  
#else

  assert (n_nodes() != 0);
  assert (n_sbdmns > 0);
  assert (n_sbdmns <= n_elem());
  

  set_n_subdomains() = n_sbdmns;
  set_n_processors() = n_sbdmns;

  // check for easy return
  if (n_sbdmns == 1)
    {
      for (unsigned int e=0; e<n_elem(); e++)
	elem(e)->subdomain_id() = 
	  elem(e)->processor_id() = 0;
      
      return;
    };
  
  _perf_log.start_event("sfc_partition()");
    
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<unsigned int> table;
  
  x.resize(n_elem(),0.);
  y.resize(n_elem(),0.);
  z.resize(n_elem(),0.);
  
  for (unsigned int e=0; e<n_elem(); e++)
    {
      const Point p = elem(e)->centroid();      
      
      x[e] = p(0);
      y[e] = p(1);

      if (_dim == 3)
	z[e] = p(2);
    };

  int size = static_cast<int>(n_elem());
  table.resize(size);
  
  if (type == "hilbert")
    sfc::hilbert(&x[0], &y[0], &z[0], &size, (int*) &table[0]);
  else if (type == "morton")
    sfc::morton(&x[0], &y[0], &z[0], &size, (int*) &table[0]);
  else
    error();

  const unsigned int wgt_per_proc = total_weight()/n_subdomains();
  unsigned int wgt = 0;
	
  for (unsigned int e=0; e<n_elem(); e++)
    {
      elem(table[e]-1)->subdomain_id() = 
	elem(table[e]-1)->processor_id() = 
	wgt/wgt_per_proc;

      wgt += elem(table[e]-1)->n_nodes();
    };
  
  _perf_log.stop_event("sfc_partiton()");

  return;
  
#endif
};



void MeshBase::distort(const real factor,
		       const bool perturb_boundary)
{
  assert (mesh_dimension() != 1);
  assert (n_nodes());
  assert (n_elem());
  assert ((factor >= 0.) && (factor <= 1.));

  _perf_log.start_event("distort()");

  std::vector<real>      hmin(n_nodes(), 1.e20);
  std::vector<short int> on_boundary(n_nodes(), 0);


  // First find nodes on the boundary and flag them
  // so that we don't move them
  if (!perturb_boundary)
    {
      for (unsigned int e=0; e<n_elem(); e++)
	if (elem(e)->active())
	  for (unsigned int s=0; s<elem(e)->n_sides(); s++)
	    if (elem(e)->neighbor(s) == NULL) // on the boundary
	      {
		const AutoPtr<Elem> side(elem(e)->build_side(s));
		
		for (unsigned int n=0; n<side->n_nodes(); n++)
		  on_boundary[side->node(n)] = 1;
	      };
    };


  // Now calculate the minimum distance to
  // neighboring nodes for each node
  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->active())
      for (unsigned int n=0; n<elem(e)->n_nodes(); n++)
	hmin[elem(e)->node(n)] = std::min(hmin[elem(e)->node(n)],
					  elem(e)->hmin());		

  
  // Now actually move the nodes
  {
    const int seed = 123456;
    
    // seed the random number generator
    srand(seed);
    
    for (unsigned int n=0; n<n_nodes(); n++)
      if (!on_boundary[n])
	{
	  // the direction, random but unit normalized
	  
	  Point dir( ((real) rand())/((real) RAND_MAX),
		     ((real) rand())/((real) RAND_MAX),
		     ((mesh_dimension() == 3) ?
		      ((real) rand())/((real) RAND_MAX) :
		      0.)
		     );
	  
	  dir(0) = (dir(0)-.5)*2.;
	  dir(1) = (dir(1)-.5)*2.;

	  if (spatial_dimension() == 3)
	    dir(2) = (dir(2)-.5)*2.;
	  
	  dir = dir.unit();

	  // if hmin[n]=1.e20 then the node is not
	  // used by any element.  We should not
	  // move it.
	  if (hmin[n] != 1.e20)
	    {
	      node(n)(0) += dir(0)*factor*hmin[n];
	      node(n)(1) += dir(1)*factor*hmin[n];
	      
	      if (spatial_dimension() == 3)
		node(n)(2) += dir(2)*factor*hmin[n];
	    }
	}
  };


  // All done  
  _perf_log.stop_event("distort()");

  return;
};



void MeshBase::translate (const real xt,
			  const real yt,
			  const real zt)
{
  const Point p(xt, yt, zt);

  for (unsigned int n=0; n<n_nodes(); n++)
    node(n) += p;
};



void MeshBase::rotate (const real,
		       const real,
		       const real)
{
  error();
};



void MeshBase::scale (const real xs,
		      const real ys,
		      const real zs)
{
  const real x_scale = xs;
  real y_scale       = ys;
  real z_scale       = zs;
  
  if (ys == 0.)
    {
      assert (zs == 0.);

      y_scale = z_scale = x_scale;
    };

  // Scale the x coordinate in all dimensions
  for (unsigned int n=0; n<n_nodes(); n++)
    node(n)(0) = node(n)(0)*x_scale;


  // Only scale the y coordinate in 2 and 3D
  if (spatial_dimension() > 1)
    {

      for (unsigned int n=0; n<n_nodes(); n++)
	node(n)(1) = node(n)(1)*y_scale;

      // Only scale the z coordinate in 3D
      if (spatial_dimension() == 3)
	{
	  for (unsigned int n=0; n<n_nodes(); n++)
	    node(n)(2) = node(n)(2)*z_scale;
	};
    };
};



std::pair<Point, Point> 
MeshBase::bounding_box() const
{
  // processor bounding box with no arguments
  // computes the global bounding box
  return processor_bounding_box();
};



Sphere
MeshBase::bounding_sphere() const
{
  std::pair<Point, Point> bbox = bounding_box();

  const real  diag = (bbox.second - bbox.first).size();
  const Point cent = (bbox.second + bbox.first)/2.;

  Sphere sphere (cent, .5*diag);

  return sphere;
};



std::pair<Point, Point> 
MeshBase::processor_bounding_box (const unsigned int pid) const
{
  assert (n_nodes() != 0);

  Point min(1.e30, 1.e30, 1.e30);
  Point max(-1.e30, -1.e30, -1.e30);

  // By default no processor is specified and we compute
  // the bounding box for the whole domain.
  if (pid == static_cast<unsigned int>(-1))
    {
      for (unsigned int n=0; n<n_nodes(); n++)
	for (unsigned int i=0; i<spatial_dimension(); i++)
	  {
	    min(i) = std::min(min(i), point(n)(i));
	    max(i) = std::max(max(i), point(n)(i));
	  };      
    }
  // if a specific processor id is specified then we need
  // to only consider those elements living on that processor
  else
    {
      for (unsigned int e=0; e<n_elem(); e++)
	if (elem(e)->processor_id() == pid)
	  for (unsigned int n=0; n<elem(e)->n_nodes(); n++)
	    for (unsigned int i=0; i<spatial_dimension(); i++)
	      {
		min(i) = std::min(min(i), point(elem(e)->node(n))(i));
		max(i) = std::max(max(i), point(elem(e)->node(n))(i));
	      };      
    };

  const std::pair<Point, Point> ret_val(min, max);

  return ret_val;  
};



Sphere
MeshBase::processor_bounding_sphere (const unsigned int pid) const
{
  std::pair<Point, Point> bbox = processor_bounding_box(pid);

  const real  diag = (bbox.second - bbox.first).size();
  const Point cent = (bbox.second + bbox.first)/2.;

  Sphere sphere (cent, .5*diag);

  return sphere;
};



std::pair<Point, Point> 
MeshBase::subdomain_bounding_box (const unsigned int sid) const
{
  assert (n_nodes() != 0);

  Point min(1.e30, 1.e30, 1.e30);
  Point max(-1.e30, -1.e30, -1.e30);

  // By default no subdomain is specified and we compute
  // the bounding box for the whole domain.
  if (sid == static_cast<unsigned int>(-1))
    {
      for (unsigned int n=0; n<n_nodes(); n++)
	for (unsigned int i=0; i<spatial_dimension(); i++)
	  {
	    min(i) = std::min(min(i), point(n)(i));
	    max(i) = std::max(max(i), point(n)(i));
	  };      
    }

  // if a specific subdomain id is specified then we need
  // to only consider those elements living on that subdomain
  else
    {
      for (unsigned int e=0; e<n_elem(); e++)
	if (elem(e)->subdomain_id() == sid)
	  for (unsigned int n=0; n<elem(e)->n_nodes(); n++)
	    for (unsigned int i=0; i<spatial_dimension(); i++)
	      {
		min(i) = std::min(min(i), point(elem(e)->node(n))(i));
		max(i) = std::max(max(i), point(elem(e)->node(n))(i));
	      };      
    };

  const std::pair<Point, Point> ret_val(min, max);

  return ret_val;  
};



Sphere MeshBase::subdomain_bounding_sphere (const unsigned int sid) const
{
  std::pair<Point, Point> bbox = subdomain_bounding_box(sid);

  const real  diag = (bbox.second - bbox.first).size();
  const Point cent = (bbox.second + bbox.first)/2.;

  Sphere sphere (cent, .5*diag);

  return sphere;
};



void MeshBase::read(const std::string&)
{
  std::cerr << "ERROR:  You shouldn't be calling this" << std::endl
	    << " Use Mesh::read() instead." << std::endl;
  error();
};



void MeshBase::write(const std::string& name)
{
  _perf_log.start_event("write()");
  
  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      write_tecplot (name);
    
    else if (name.rfind(".plt") < name.size())
      write_tecplot_binary (name);

    else if (name.rfind(".ucd") < name.size())
      write_ucd (name);

    else if (name.rfind(".gmv") < name.size())
      {
	if (n_subdomains() > 1)
	  write_gmv_binary(name, NULL, NULL, true);
	else
	  write_gmv_binary(name);
      }
  };

  _perf_log.stop_event("write()");
};



void MeshBase::write(const std::string& name,
		     std::vector<number>& v,
		     std::vector<std::string>& vn)
{
  _perf_log.start_event("write()");
  
  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      write_tecplot (name, &v, &vn);
    
    else if (name.rfind(".plt") < name.size())
      write_tecplot_binary (name, &v, &vn);
    
    else if (name.rfind(".gmv") < name.size())
      {
	if (n_subdomains() > 1)
	  write_gmv_binary(name, &v, &vn, true);
	else
	  write_gmv_binary(name, &v, &vn);
      }
  };

  _perf_log.stop_event("write()");
};



#ifdef USE_COMPLEX_NUMBERS

const char* MeshBase::complex_filename(const std::string& _n,
				       unsigned int r_o_c)
{
  std::string loc=_n;
  if (r_o_c == 0)
    loc.append(".real");
  else
    loc.append(".imag");
  return loc.c_str();
};


void MeshBase::prepare_complex_data(const std::vector<number>* source,
				    std::vector<real>* real_part,
				    std::vector<real>* imag_part)
{
  real_part->resize(source->size());
  imag_part->resize(source->size());


  for (unsigned int i=0; i< source->size(); i++)
    {
      (*real_part)[i] = (*source)[i].real();
      (*imag_part)[i] = (*source)[i].real();
    };
};


#endif




