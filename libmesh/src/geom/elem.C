// $Id: elem.C,v 1.20 2003-05-15 23:34:35 benkirk Exp $

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
#include <iterator>
#include <set>

// Local includes
#include "elem.h"
#include "mesh_base.h"
#include "fe_type.h"
#include "fe_interface.h"
#include "edge_edge2.h"
#include "edge_edge3.h"
#include "edge_inf_edge2.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad4.h"
#include "face_quad8.h"
#include "face_quad9.h"
#include "face_inf_quad.h"
#include "face_inf_quad4.h"
#include "face_inf_quad6.h"
#include "cell_tet.h"
#include "cell_tet4.h"
#include "cell_tet10.h"
#include "cell_hex.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_hex27.h"
#include "cell_inf_hex.h"
#include "cell_inf_hex8.h"
#include "cell_inf_hex16.h"
#include "cell_inf_hex18.h"
#include "cell_prism.h"
#include "cell_prism6.h"
#include "cell_prism15.h"
#include "cell_prism18.h"
#include "cell_inf_prism.h"
#include "cell_inf_prism6.h"
#include "cell_inf_prism12.h"
#include "cell_pyramid.h"
#include "cell_pyramid5.h"



// ------------------------------------------------------------
// Elem class member funcions
Elem* Elem::build(const ElemType type,
		  const Elem* p)
{
  switch (type)
    {
      // 1D elements 
    case EDGE2:
      {
	return new Edge2(p);
      }
    case EDGE3:
      {
	return new Edge3(p);
      }


      
      // 2D elements
    case TRI3:
      {
	return new Tri3(p);
      }
    case TRI6:
      {
	return new Tri6(p);
      }
    case QUAD4:
      {
	return new Quad4(p);
      }
    case QUAD8:
      {
	return new Quad8(p);
      }
    case QUAD9:
      {
	return new Quad9(p);
      }
   

      // 3D elements
    case TET4:
      {
	return new Tet4(p);
      }
    case TET10:
      {
	return new Tet10(p);
      }
    case HEX8:
      {
	return new Hex8(p);
      }
    case HEX20:
      {
	return new Hex20(p);
      }
    case HEX27:
      {
	return new Hex27(p);
      }
    case PRISM6:
      {
	return new Prism6(p);
      }
    case PRISM15:
      {
	return new Prism15(p);
      }
    case PRISM18:
      {
	return new Prism18(p);
      }
    case PYRAMID5:
      {
	return new Pyramid5(p);
      }



#ifdef ENABLE_INFINITE_ELEMENTS

      // 1D infinite elements
    case INFEDGE2:
      {
	return new InfEdge2(p);
      }


      // 2D infinite elements
    case INFQUAD4:
      {
	return new InfQuad4(p);
      }
    case INFQUAD6:
      {
	return new InfQuad6(p);
      }

   
    // 3D infinite elements
    case INFHEX8:
      {
	return new InfHex8(p);
      }
    case INFHEX16:
      {
	return new InfHex16(p);
      }
    case INFHEX18:
      {
	return new InfHex18(p);
      }
    case INFPRISM6:
      {
	return new InfPrism6(p);
      }
    case INFPRISM12:
      {
	return new InfPrism12(p);
      }

#endif
           
    default:
      {
	std::cerr << "ERROR: Undefined element type!." << std::endl;
	error();
      }
    }
    

  // This will never happen...  Look at the case-structure above.
  // error() will abort, so we won't get here...
  error();
  
  return NULL;
}



unsigned int Elem::which_neighbor_am_i (const Elem* e) const
{
  assert (e != NULL);
  
  for (unsigned int s=0; s<this->n_neighbors(); s++)
    if (this->neighbor(s) == e)
      return s;
    

  std::cerr << "ERROR:  Elements are not neighbors!" 
	    << std::endl;

  error();

  return static_cast<unsigned int>(-1);
}



void Elem::write_tecplot_connectivity(std::ostream& out) const
{
  assert (!out.bad());
  assert (_nodes != NULL);
  
  for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
    {
      std::vector<unsigned int> conn = tecplot_connectivity(sc);
      
      // Orignial code
      //for (unsigned int i=0; i<conn.size(); i++)
      // 	out << conn[i] << " ";
      
      // New code
      std::copy(conn.begin(),
 		conn.end(),
 		std::ostream_iterator<unsigned int>(out, " "));
      
      out << std::endl;
    }
}



unsigned int Elem::key() const
{
  const unsigned int nv = n_vertices();
  
  assert (nv != 0);

  std::vector<unsigned int> vec (nv, 0);

  for (unsigned int v=0; v<nv; v++)
    vec[v] = this->node(v);

  
  std::sort(vec.begin(), vec.end());
  
  unsigned int n       = vec[0];
  const unsigned int m = vec[0];
  
  for (unsigned int i=1; i<nv; i++)
    n = n^vec[i];
  

  return n + m*m;  
} 



Point Elem::centroid() const
{
  Point cp;

  for (unsigned int n=0; n<n_vertices(); n++)
    cp.add (this->point(n));

  return (cp /= static_cast<Real>(n_vertices()));    
}



Real Elem::hmin() const
{
  Real h_min=1.e30;

  for (unsigned int n_outer=0; n_outer<n_vertices(); n_outer++)
    for (unsigned int n_inner=0; n_inner<n_vertices(); n_inner++)
      if (n_outer != n_inner) // would create false 0
	{
	  const Point diff = (point(n_outer) - point(n_inner));
	  
	  h_min = std::min(h_min,diff.size());
	}

  return h_min;
}



Real Elem::hmax() const
{
  Real h_max=0;

  for (unsigned int n_outer=0; n_outer<n_vertices(); n_outer++)
    for (unsigned int n_inner=0; n_inner<n_vertices(); n_inner++)
      if (n_outer != n_inner) // will be 0, definately _not_ max
	{
	  const Point diff = (point(n_outer) - point(n_inner));
	  
	  h_max = std::max(h_max,diff.size());
	}

  return h_max;
}



Real Elem::length(const unsigned int n1, 
		  const unsigned int n2) const
{
  assert ( n1 < n_vertices() );
  assert ( n2 < n_vertices() );

  return (point(n1) - point(n2)).size();
}



bool Elem::operator == (const Elem& rhs) const
{
  assert (n_nodes());
  assert (rhs.n_nodes());

  // Elements can only be equal if they
  // contain the same number of nodes.
  if (n_nodes() == rhs.n_nodes())
    {
      // Create a set that contains our global
      // node numbers and those of our neighbor.
      // If the set is the same size as the number
      // of nodes in both elements then they must
      // be connected to the same nodes.     
      std::set<unsigned int> nodes_set;

      for (unsigned int n=0; n<this->n_nodes(); n++)
	{
	  nodes_set.insert(node(n));
	  nodes_set.insert(rhs.node(n));
	}

      // If this passes the elements are connected
      // to the same global nodes
      if (nodes_set.size() == n_nodes())
	return true;
    }

  // If we get here it is because the elements either
  // do not have the same number of nodes or they are
  // connected to different nodes.  Either way they
  // are not the same element
  return false;
}



void Elem::write_ucd_connectivity(std::ostream &out) const
{
  assert (out);
  assert (_nodes != NULL);

  // Original Code
  for (unsigned int i=0; i<n_nodes(); i++)
    out << node(i)+1 << "\t";

  out << std::endl;
}



Real Elem::quality (const ElemQuality q) const
{
  switch (q)
    {    
      /**
       * I don't know what to do for this metric. 
       */
    default:
      {
	here();

	std::cerr << "ERROR:  unknown quality metric: "
		  << q 
		  << std::endl
		  << "Cowardly returning 1."
		  << std::endl;

	return 1.;
      }
    }

    
    // Will never get here...
    error();
    return 0.;
}



bool Elem::contains_point (const Point& p) const
{
  // Declare a basic FEType.  Will ue a Lagrange
  // element by default.
  FEType fe_type(default_order());
  
  const Point mapped_point = FEInterface::inverse_map(dim(),
						      fe_type,
						      this,
						      p,
						      1.e-4,
						      false);

  return FEInterface::on_reference_element(mapped_point, this->type());
}



void Elem::nullify_neighbors ()
{
  // Tell any of my neighbors about my death...
  // Looks strange, huh?
  for (unsigned int n=0; n<this->n_neighbors(); n++)
    if (this->neighbor(n) != NULL)
      {
	Elem* neighbor = this->neighbor(n);

	// Note:  it is possible that I see the neighbor
	// (which is at a lower refinement level than me)
	// but they don't see me, so avoid that case.
	if (neighbor->level() == this->level())
	  {	
// 	    std::cout << "this=" << this
// 		      << ", neighbor=" << neighbor
// 		      << ", this->id()=" << this->id()
// 		      << std::endl;
	
	    const unsigned int w_n_a_i = neighbor->which_neighbor_am_i(this);
	    neighbor->set_neighbor(w_n_a_i, NULL);
	    this->set_neighbor(n, NULL);
	  }
      }
}



/**
 * The following functions only apply when
 * AMR is enabled and thus are not present
 * otherwise.
 */ 
#ifdef ENABLE_AMR

void Elem::refine (MeshBase& mesh)
{
  assert (this->refinement_flag() == Elem::REFINE);
  assert (this->active());
  assert (_children == NULL);

  // Create my children
  {
    _children = new Elem*[this->n_children()];

    for (unsigned int c=0; c<this->n_children(); c++)
      {
	_children[c] = Elem::build(this->type(), this);
	_children[c]->set_refinement_flag(Elem::JUST_REFINED);
      }
  }


  // Compute new nodal locations
  // and asssign nodes to children
  {
    std::vector<std::vector<Point> >  p(this->n_children());
    
    for (unsigned int c=0; c<this->n_children(); c++)
      p[c].resize(this->child(c)->n_nodes());
    

    // compute new nodal locations
    for (unsigned int c=0; c<this->n_children(); c++)
      for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
	for (unsigned int n=0; n<this->n_nodes(); n++)
	  if (this->embedding_matrix(c,nc,n) != 0.)
	    p[c][nc].add_scaled (this->point(n), this->embedding_matrix(c,nc,n));
    
    
    // assign nodes to children & add them to the mesh
    for (unsigned int c=0; c<this->n_children(); c++)
      {
	for (unsigned int nc=0; nc<this->child(c)->n_nodes(); nc++)
	  _children[c]->set_node(nc) = mesh.mesh_refinement.add_point(p[c][nc]);

	mesh.add_elem(this->child(c), mesh.mesh_refinement.new_element_number());
      }
  }


  
  // Possibly add boundary information
  for (unsigned int s=0; s<this->n_neighbors(); s++)
    if (this->neighbor(s) == NULL)
      {
	const short int id = mesh.boundary_info.boundary_id(this, s);
	
	if (id != mesh.boundary_info.invalid_id)
	  for (unsigned int sc=0; sc<this->n_children_per_side(s); sc++)
	    mesh.boundary_info.add_side(this->child(this->side_children_matrix(s,sc)), s, id);
      }

  
  // Un-set my refinement flag now
  this->set_refinement_flag(Elem::DO_NOTHING);
}



void Elem::coarsen()
{
  assert (this->refinement_flag() == Elem::COARSEN);
  assert (!this->active());

  // Delete the storage for my children
  delete [] _children;

  _children = NULL;

  this->set_refinement_flag(Elem::DO_NOTHING);
}

#endif // #ifdef ENABLE_AMR


