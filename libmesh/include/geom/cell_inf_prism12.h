// $Id: cell_inf_prism12.h,v 1.6 2005-02-22 22:17:31 jwpeterson Exp $

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



#ifndef __cell_inf_prism12_h__
#define __cell_inf_prism12_h__

// C++ includes

// Local includes
#include "libmesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS

#include "cell_inf_prism.h"




/**
 * The \p InfPrism12 is an infinite element in 3D composed of 12 nodes.
 * It is numbered like this:
   \verbatim
   INFPRISM12:
            5      
            o      
            :      
            :      
            :   
     11 o   :   o 10
        :  2:   :   
        :   o   :        closer to infinity
        :  . .  :   
   3o   : . o9. :   o4
    |   :.  |  .:   |
    |   o   |   o   |
    |  . 8  |  7 .  |
    | .     |     . |
    |.      |      .|     base face
    o-------o-------o
    0       6       1
   \endverbatim
 */

// ------------------------------------------------------------
// InfPrism12 class definition
class InfPrism12 : public InfPrism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  InfPrism12  (const Elem* p=NULL);

  /**
   * @returns 12.  The \p InfPrism12 has 12 nodes.
   */
  unsigned int n_nodes() const { return 12; }
  
  /**
   * @returns \p INFPRISM12
   */
  ElemType     type () const   { return INFPRISM12; }
  
  /**
   * @returns 4
   */
  unsigned int n_sub_elem() const { return 4; }
  
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
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }
  
  /**
   * Returns a \p TRI6 built coincident with face 0, an \p INFQUAD6 
   * built coincident with faces 1 to 3.  Note that the \p AutoPtr<Elem>
   * takes care of freeing memory.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  virtual void connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const;

//   void tecplot_connectivity(const unsigned int sc,
// 			    std::vector<unsigned int>& conn) const;
  
//   void vtk_connectivity(const unsigned int,
// 			std::vector<unsigned int>*) const
//   { error(); }
  
//   unsigned int vtk_element_type (const unsigned int) const
//   { return 13; }

  /**
   * @returns 2 for all \p n
   */
  unsigned int n_second_order_adjacent_vertices (const unsigned int) const
      { return 2; }

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   * Note that \p n is counted as depicted above, \f$ 6 \le n < 12 \f$.
   */
  unsigned short int second_order_adjacent_vertex (const unsigned int n,
						   const unsigned int v) const;



protected:  

#ifdef ENABLE_AMR
  
  /**
   * Matrix used to create the elements children.
   */
  float embedding_matrix (const unsigned int i,
			  const unsigned int j,
			  const unsigned int k) const
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[4][12][12];
  
#endif


private:
  
  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[6][2];
  
};



// ------------------------------------------------------------
// InfPrism12 class member functions
inline
InfPrism12::InfPrism12(const Elem* p) :
  InfPrism(InfPrism12::n_nodes(), p) 
{
}


#endif  // ifdef ENABLE_INFINITE_ELEMENTS

#endif
