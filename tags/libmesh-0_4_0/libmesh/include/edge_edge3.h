// $Id: edge_edge3.h,v 1.11 2003-05-22 21:18:02 benkirk Exp $

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



#ifndef __edge3_h__
#define __edge3_h__


// Local includes
#include "mesh_common.h"
#include "edge.h"



/**
 * The \p Edge3 is an element in 1D composed of 3 nodes. It is numbered
 * like this:
 *
 * \verbatim
 *  EGDE3: o----o----o
 *         0    2    1
 * \endverbatim
 */

// ------------------------------------------------------------
// Edge3 class definition
class Edge3 : public Edge
{
 public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Edge3 (const Elem* p=NULL);
  
  /**
   * @returns 3
   */
  unsigned int n_nodes() const { return 3; }

  /**
   * @returns 2
   */
  unsigned int n_sub_elem() const { return 2; }

  /**
   * @returns \p EDGE3
   */
  ElemType type()  const { return EDGE3; }
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int se=0) const;

  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 3; }


#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p false.  This is a finite element. 
   */
  bool infinite () const { return false; }

#endif

  
protected:

  
#ifdef ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int,
			 const unsigned int,
			 const unsigned int) const
  { error(); return 0.; }
  
#endif
};




// ------------------------------------------------------------
// Edge3 class member functions
inline
Edge3::Edge3(const Elem* p) :
  Edge(Edge3::n_nodes(), p) 
{
}






#endif
