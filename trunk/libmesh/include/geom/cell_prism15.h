// $Id: cell_prism15.h,v 1.4 2004-11-15 22:09:11 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __cell_prism15_h__
#define __cell_prism15_h__

// C++ includes

// Local includes
#include "cell_prism.h"




/**
 * The \p Prism15 is an element in 3D composed of 15 nodes.
 * It is numbered like this:
   \verbatim
   PRISM15:
            5
            o
           /:\
          / : \
         /  :  \
        /   :   \
    14 o    :    o 13
      /     :     \ 
     /      :      \
    /       o 11    \
 3 /        :        \4
  o---------o---------o
  |         :12       |
  |         :         |
  |         :         |
  |         o         |
  |        .2.        |
  |       .   .       |
9 o      .     .      o 10
  |     .       .     |
  |  8 o         o 7  |
  |   .           .   |
  |  .             .  |
  | .               . |
  |.                 .|
  o---------o---------o
  0         6         1

  
   \endverbatim
 */

// ------------------------------------------------------------
// Prism class definition
class Prism15 : public Prism
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Prism15  (const Elem* p=NULL);
  
  /**
   * @returns \p PRISM15
   */
  ElemType     type () const   { return PRISM15; }

  /**
   * @returns 15
   */
  unsigned int n_nodes() const { return 15; }

  /**
   * @returns 1
   */
  unsigned int n_sub_elem() const { return 1; }
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; }
  
  /**
   * Builds a \p QUAD8 or \p TRI6 built coincident with face i.  
   * The \p AutoPtr<Elem> handles the memory aspect.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  virtual void connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const;

//   void tecplot_connectivity(const unsigned int sc,
// 			    std::vector<unsigned int>& conn) const;
  
  
//   void vtk_connectivity(const unsigned int sc,
// 			std::vector<unsigned int> *conn = NULL) const;
  
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
   * Note that \p n is counted as depicted above, \f$ 6 \le n < 15 \f$.
   */
  unsigned short int second_order_adjacent_vertex (const unsigned int n,
						   const unsigned int v) const;
  
  
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
// Prism15 class member functions
inline
Prism15::Prism15(const Elem* p) :
  Prism(Prism15::n_nodes(), p) 
{
}



#endif // #define __cell_prism15_h__

