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



#ifndef LIBMESH_QUADRATURE_COMPOSITE_H
#define LIBMESH_QUADRATURE_COMPOSITE_H

// Local includes
#include "libmesh/quadrature.h"

// C++ includes

namespace libMesh
{




/**
 * This class implements generic composite quadrature rules.
 * Composite quadrature rules are constructed from any of the
 * supported rules by breaking an element into subelements and
 * applying the base rule on each subelement.
 */

// ------------------------------------------------------------
// QComposite class definition
template <class QSubCell>
class QComposite : public QSubCell
{  
 public:

  /**
   * Import base class data into namespace scope.
   */
  using QSubCell::_dim;
  
  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QComposite (const unsigned int _dim,
	      const Order _order=INVALID_ORDER);

  /**
   * Destructor.
   */
  ~QComposite();

  /**
   * @returns \p QCOMPOSITE
   */
  QuadratureType type() const { return QCOMPOSITE; }


 private:
  
  /**
   * Subcell quadrature object.
   */
  QSubCell _q_subcell;

};



// ------------------------------------------------------------
// QComposite class members
template <class QSubCell>
inline
QComposite<QSubCell>::QComposite(const unsigned int d,
				 const Order o) :
  QSubCell(d,o), // explicitly call base class constructor
  _q_subcell(d,o)
{
  // explicitly call the init function in 1D since the
  // other tensor-product rules require this one.
  // note that EDGE will not be used internally, however
  // if we called the function with INVALID_ELEM it would try to
  // be smart and return, thinking it had already done the work.
  if (_dim == 1)
    this->init(EDGE2);
}




template <class QSubCell>
inline
QComposite<QSubCell>::~QComposite()
{
}


} // namespace libMesh



#endif // LIBMESH_QUADRATURE_COMPOSITE_H
