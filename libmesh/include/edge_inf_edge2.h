// $Id: edge_inf_edge2.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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



#ifndef __inf_edge2_h__
#define __inf_edge2_h__


// Local includes
#include "mesh_common.h"

#ifdef ENABLE_INFINITE_ELEMENTS

#include "edge.h"


/**
 * The \p InfEdge2 is an infinte element in 1D composed of 2 nodes. 
 * It is numbered like this:
 *
 * \verbatim
 *  INFEDGE2: o--------o
 *            0        1
 * \endverbatim
 */
// ------------------------------------------------------------
// InfEdge2 class definition
class InfEdge2 : public Edge
{
 public:

  /**
   * Constructor.  By default this element has no parent.
   */
  InfEdge2 (Edge* p=NULL);

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; };
  
  /**
   * @returns \p INFEDGE2
   */
  ElemType type()  const { return INFEDGE2; };
  
  /**
   * @returns FIRST
   */
  Order default_order() const { return FIRST; };
  
  const std::vector<unsigned int> tecplot_connectivity(const unsigned int se=0) const;
  
#ifdef ENABLE_AMR

  /**
   * Refine the element.
   */
  void refine(Mesh& mesh)
    { error(); return; };
  
  /**
   * Refine the element.
   */
  void coarsen() { error(); return; };

#endif
  
 private:
  
};




// ------------------------------------------------------------
// InfEdge2 class member functions
inline
InfEdge2::InfEdge2(Edge* p) :
  Edge(InfEdge2::n_nodes(), p) 
{
};

#endif

#endif
