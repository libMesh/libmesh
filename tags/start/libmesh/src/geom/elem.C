// $Id: elem.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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

// Local includes
#include "elem.h"
#include "mesh.h"
#include "cell.h"
#include "fe_interface.h"

// Temporary 1D element includes
#include "edge_edge2.h"
#include "edge_edge3.h"
#include "edge_inf_edge2.h"

// Temporary 2D element includes
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad4.h"
#include "face_quad8.h"
#include "face_quad9.h"
#include "face_inf_quad4.h"
#include "face_inf_quad6.h"

// Temporary 3D element includes
#include "cell_tet.h"
#include "cell_tet4.h"
#include "cell_tet10.h"
#include "cell_hex.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_hex27.h"
#include "cell_inf_hex8.h"
#include "cell_inf_hex16.h"
#include "cell_inf_hex18.h"
#include "cell_prism.h"
#include "cell_prism6.h"
#include "cell_inf_prism6.h"
#include "cell_inf_prism12.h"
#include "cell_pyramid.h"
#include "cell_pyramid5.h"



// ------------------------------------------------------------
// Elem class member funcions
Elem* Elem::build(const ElemType type)
{
  switch (type)
    {
      // 1D elements
    case EDGE2:
      {
	return new Edge2;
      }
    case EDGE3:
      {
	return new Edge3;
      }


      
      // 2D elements
    case TRI3:
      {
	return new Tri3;
      }
    case TRI6:
      {
	return new Tri6;
      }
    case QUAD4:
      {
	return new Quad4;
      }
    case QUAD8:
      {
	return new Quad8;
      }
    case QUAD9:
      {
	return new Quad9;
      }
   

      // 3D elements
    case TET4:
      {
	return new Tet4;
      }
    case TET10:
      {
	return new Tet10;
      }
    case HEX8:
      {
	return new Hex8;
      }
    case HEX20:
      {
	return new Hex20;
      }
    case HEX27:
      {
	return new Hex27;
      }
    case PRISM6:
      {
	return new Prism6;
      }
    case PRISM18:
      {
	error();
	return new Prism6;
      }
    case PYRAMID5:
      {
	return new Pyramid5;
      }



#ifdef ENABLE_INFINITE_ELEMENTS

      // 1D infinite elements
    case INFEDGE2:
      {
	return new InfEdge2;
      }


      // 2D infinite elements
    case INFQUAD4:
      {
	return new InfQuad4;
      }
    case INFQUAD6:
      {
	return new InfQuad6;
      }

   
    // 3D infinite elements
    case INFHEX8:
      {
	return new InfHex8;
      }
    case INFHEX16:
      {
	return new InfHex16;
      }
    case INFHEX18:
      {
	return new InfHex18;
      }
    case INFPRISM6:
      {
	return new InfPrism6;
      }
    case INFPRISM12:
      {
	return new InfPrism12;
      }

#endif
           
    default:
      {
	std::cout << "Undefined element type!." << std::endl;
	error();
      };
    };
    

  // This will never happen...  Look at the case-structure above.
  // error() will abort, so we won't get here...
  error();
  
  return NULL;
};



unsigned int Elem::which_neighbor_am_i (const Elem* e) const
{
  assert (e != NULL);

  for (unsigned int s=0; s<n_sides(); s++)
    if (neighbor(s) == e)
      return s;

  std::cerr << "ERROR:  Elements are not neighbors!" 
	    << std::endl;

  error();

  return static_cast<unsigned int>(-1);
};



void Elem::write_tecplot_connectivity(std::ostream& out) const
{
  assert (out.good());
  assert (!_nodes.empty());
  
  for (unsigned int sc=0; sc<n_sub_elem(); sc++)
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
    };
};



unsigned int Elem::node (const unsigned int i) const
{
  assert (i < n_nodes());

  return _nodes[i];
};



unsigned int & Elem::node (const unsigned int i)
{
  assert (i < n_nodes());

  return _nodes[i];
};



unsigned int Elem::key() const
{
  assert (n_nodes());

  return key(_nodes);
};



unsigned int Elem::key(std::vector<unsigned int> vec) const
{
  std::sort(vec.begin(), vec.end());
  
  unsigned int n = vec[0];

  unsigned int m = vec[0];
  
  for (unsigned int i=1; i<vec.size(); i++)
    {
      n = n^vec[i];
    };

  return n + m*m;
};
 



Point Elem::centroid(const MeshBase& mesh) const
{
  Point cp;

  for (unsigned int n=0; n<n_vertices(); n++)
    cp += mesh.vertex(node(n));

  return cp/((real) n_vertices());    
};



real Elem::hmin(const MeshBase& mesh) const
{
  real h_min=1.e30;

  for (unsigned int n_outer=0; n_outer<n_vertices(); n_outer++)
    for (unsigned int n_inner=0; n_inner<n_vertices(); n_inner++)
      if (n_outer != n_inner)
	{
	  const Point diff = (mesh.vertex(node(n_outer)) -
			      mesh.vertex(node(n_inner))
			      );
	  
	  h_min = std::min(h_min,diff.size());
	};

  return h_min;
};



real Elem::hmax(const MeshBase& mesh) const
{
  real h_max=0;

  for (unsigned int n_outer=0; n_outer<n_vertices(); n_outer++)
    for (unsigned int n_inner=0; n_inner<n_vertices(); n_inner++)
      if (n_outer != n_inner)
	{
	  const Point diff = (mesh.vertex(node(n_outer)) -
			      mesh.vertex(node(n_inner))
			      );
	  
	  h_max = std::max(h_max,diff.size());
	};

  return h_max;
};



real Elem::length(const MeshBase& mesh,
		  const unsigned int n1, 
		  const unsigned int n2) const
{
  assert ( n1 < n_vertices() );
  assert ( n2 < n_vertices() );

  return (mesh.vertex(node(n1)) - mesh.vertex(node(n2))).size();
}



bool Elem::operator == (const Elem& rhs) const
{
  assert (n_nodes());
  
  std::vector<unsigned int> rnodes = rhs.get_nodes();
  std::vector<unsigned int> lnodes = get_nodes();

  // order the vectors
  std::sort(rnodes.begin(),rnodes.end());
  std::sort(lnodes.begin(),lnodes.end());

  if (rnodes == lnodes)
    return true;
 
  return false;
};



bool Elem::operator < (const Elem& rhs) const
{
  assert (n_nodes());

  if (n_nodes() < rhs.n_nodes())
    return true;
  
  else if (n_nodes() > rhs.n_nodes())
    return false;
  
  else // n_nodes() == rhs.n_nodes()
    for (unsigned int n=0; n<n_nodes(); n++)
      if ((node(n)) < (rhs.node(n)))
	return true;
    
  return false;
};    



void Elem::write_ucd_connectivity(std::ostream &out) const
{
  assert (out);
  assert (!_nodes.empty());

  // Original Code
  // for (unsigned int i=0; i<nodes.size(); i++)
  //     out << node(i)+1 << "\t";

  // New code
  std::transform(_nodes.begin(),
		 _nodes.end(),
		 std::ostream_iterator<unsigned int>(out, "\t"),
		 std::bind2nd(std::plus<unsigned int>(), 1));

  out << std::endl;
};



real Elem::quality (const MeshBase& mesh, const ElemQuality q) const
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
};



bool Elem::contains_point (const MeshBase& mesh,
			   const Point& p) const
{
  // Declare a basic FEType.  Will ue a Lagrange
  // element by default.
  FEType fe_type(default_order());
  
  const Point mapped_point = FEInterface::inverse_map(dim(),
						      fe_type,
						      mesh,
						      this,
						      p);

  return FEBase::on_reference_element(mapped_point, type());
};



/**
 * The following functions only apply when
 * AMR is enabled and thus are not present
 * otherwise.
 */ 
#ifdef ENABLE_AMR



unsigned int Elem::level() const
{
  // if I don't have a parent I was
  // created directly from file
  // or by the user, so I am a
  // level-0 element
  if (parent() == NULL)
    return 0;

  // otherwise we are at a level one
  // higher than our parent
  return (parent()->level() + 1);
};

#endif


