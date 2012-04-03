// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __quadrature_grid_h__
#define __quadrature_grid_h__

// Local includes
#include "quadrature.h"

namespace libMesh
{




/**
 * This class creates quadrature points on a uniform grid, Order
 * points on a side.
 */

// ------------------------------------------------------------
// QGrid class definition

class QGrid : public QBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QGrid (const unsigned int _dim,
	 const Order _order=INVALID_ORDER);

  /**
   * Destructor.
   */
  ~QGrid();

  /**
   * @returns \p QGRID
   */
  QuadratureType type() const { return QGRID; }


 private:

  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_2D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);

};



// ------------------------------------------------------------
// QGauss class members
inline
QGrid::QGrid(const unsigned int d,
	     const Order o) : QBase(d,o)
{
}




inline
QGrid::~QGrid()
{
}


} // namespace libMesh



#endif
