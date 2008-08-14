// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __quadrature_clough_h__
#define __quadrature_clough_h__

// Local includes
#include "quadrature.h"




/**
 * This class creates a gaussian quadrature rule duplicated for each
 * subelement of a Clough-Tocher divided macroelement.
 */

// ------------------------------------------------------------
// QClough class definition

class QClough : public QBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QClough (const unsigned int _dim,
	   const Order _order=INVALID_ORDER);

  /**
   * Destructor.
   */
  ~QClough();

  /**
   * @returns \p QCLOUGH
   */
  QuadratureType type() const { return QCLOUGH; }

 
 private:

  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_2D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);

};



// ------------------------------------------------------------
// QClough class members
inline
QClough::QClough(const unsigned int d,
	         const Order o) : QBase(d,o)
{
}




inline
QClough::~QClough()
{
}




#endif
