// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __edge_inf_edge2_h__
#define __edge_inf_edge2_h__


// Local includes
#include "libmesh_common.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#include "edge.h"


/**
 * The \p InfEdge2 is an infinte element in 1D composed of 2 nodes. 
 * It is numbered like this:
 *
   \verbatim
    INFEDGE2:

        o         closer to infinity
        | 1
        |
        |
        |
        o         base node
          0       
   \endverbatim
 */
// ------------------------------------------------------------
// InfEdge2 class definition
class InfEdge2 : public Edge
{
 public:

  /**
   * Constructor.  By default this element has no parent.
   */
  InfEdge2 (Elem* p=NULL) :
    Edge(InfEdge2::n_nodes(), p) {}

  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  InfEdge2(const unsigned int nn,
	   const unsigned int ns,
	   Elem* p) :
    Edge(nn, ns, p) {}

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
  
  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const;
  
  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
			       const unsigned int s) const;
  
  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge (i.e. "returns true" in 1D)
   */
  virtual bool is_node_on_edge(const unsigned int n,
			       const unsigned int e) const;
  
  /**
   * @returns \p INFEDGE2
   */
  ElemType type()  const { return INFEDGE2; }
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; }

  virtual void connectivity(const unsigned int se,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const;

//   const std::vector<unsigned int> tecplot_connectivity(const unsigned int se=0) const;
  
//   void vtk_connectivity(const unsigned int se,
// 			std::vector<unsigned int>*conn = NULL) const;
  
//   unsigned int vtk_element_type (const unsigned int) const
//   { return 3; }


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p true.  This is an infinite element. 
   */
  bool infinite () const { return true; }

  /**
   * @returns the origin of this infinite element.
   */
  Point origin () const;

#endif

  
protected:

  
#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int,
			 const unsigned int,
			 const unsigned int) const
  { libmesh_error(); return 0.; }

#endif
};




// ------------------------------------------------------------
// InfEdge2 class member functions
inline
Point InfEdge2::origin () const
{
  return ( this->point(0)*2. - this->point(1) );
}


#endif

#endif
