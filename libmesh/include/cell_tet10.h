// $Id: cell_tet10.h,v 1.2 2003-01-20 16:31:21 jwpeterson Exp $

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



#ifndef __cell_tet10_h__
#define __cell_tet10_h__

// C++ includes

// Local includes
#include "cell_tet.h"




/**
 * The \p Tet10 is an element in 3D composed of 10 nodes.
 * It is numbered like this:
  \verbatim
              3
  TET10:      o
             /|\
            / | \
        7  /  |  \9
          o   |   o
         /    |8   \
        /     o     \
       /    6 |      \
    0 o-----o-|-------o 2
       \      |      /
        \     |     /
         \    |    /
        4 o   |   o 5
           \  |  /
            \ | /
             \|/
              o
              1
   \endverbatim
 */

// ------------------------------------------------------------
// Tet10 class definition
class Tet10 : public Tet
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Tet10  (Cell* p=NULL);
  
  /**
   * @returns \p TET10
   */
  ElemType     type ()   const { return TET10; };

  /**
   * @returns 10
   */
  unsigned int n_nodes() const { return 10; };

  /**
   * @returns 8
   */
  unsigned int n_sub_elem() const { return 8; };
  
  /**
   * @returns SECOND
   */
  Order default_order() const { return SECOND; };
  
  /**
   * Builds a TRI6 built coincident with face i.  This
   * method allocates memory, so be sure to delete
   * the returned pointer when it is no longer needed.
   */
  AutoPtr<Elem> build_side (const unsigned int i) const;

  const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const;
  
  
  void vtk_connectivity(const unsigned int sc,
			std::vector<unsigned int> *conn = NULL) const;
  
  unsigned int vtk_element_type (const unsigned int) const
  { return 10; };
  
#ifdef ENABLE_AMR

  /**
   * Refine the element.
   */
  void refine(Mesh&)
  { error(); return; };
  
  /**
   * Coarsen the element.
   */
  void coarsen()
  { error(); return; };

#endif
  
private:
  
};



// ------------------------------------------------------------
// Tet10 class member functions
inline
Tet10::Tet10(Cell* p) :
  Tet(Tet10::n_nodes(), p) 
{
};



#endif
