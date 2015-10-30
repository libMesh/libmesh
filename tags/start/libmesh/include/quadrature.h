// $Id: quadrature.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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



#ifndef __quadrature_h__
#define __quadrature_h__

// C++ includes
#include <vector>

// Local includes
#include "mesh_common.h"
#include "point.h"
#include "elem_type.h"
#include "order.h"



/**
 * This is the \p QBase class.  It provides the basic functionality
 * from which various quadrature rules can be derived.  The class contains
 * \p dim dimensional points describing the quadrature locations
 * (referenced to a master object) and associated weights.
 *
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// QBase class definition

class QBase
{
protected:

  /**
   * Constructor. Protected to prevent instantiation
   * of this base class.
   */
  QBase (const unsigned int _dim,
	 const Order _order=INVALID_ORDER);
  
public:
  
  /**
   * Destructor.
   */
  virtual ~QBase() {};

  /**
   * @returns the number of points associated with the quadrature rule.
   */    
  unsigned int n_points() const
    { assert (!_points.empty()); return _points.size(); };

  /**
   * @returns a \p std::vector containing the quadrature point locations
   * on a reference object.
   */
  const std::vector<Point>& get_points()  const { return _points;  };

  /**
   * @returns a \p std::vector containing the quadrature weights.
   */
  const std::vector<real>& get_weights() const { return _weights; };

  /**
   * @returns the \f$ i^{th} \f$ quadrature point on the reference object.
   */
  Point qp(const unsigned int i) const
    { assert (i < _points.size()); return _points[i]; };

  /**
   * @returns the \f$ i^{th} \f$ quadrature weight.
   */
  real w(const unsigned int i) const
    { assert (i < _weights.size()); return _weights[i]; };

  /**
   * Initializes the data structures to contain a quadrature rule
   * for an object of type \p type.
   */
  virtual void init (const ElemType _type=INVALID_ELEM) = 0;
  
  /**
   * Initializes the data structures to contain a quadrature rule
   * for side \p side of object of type \p _type.
   */
  virtual void init (const ElemType _type,
		     const unsigned int side) = 0;

  /**
   * @returns the order of the quadrature rule.   
   */
  Order get_order() const { return _order; };
  
  /**
   * Prints information relevant to the quadrature rule.
   */
  void print_info() const;



protected:

  
  /**
   * The dimension
   */
  const unsigned int _dim;
  
  /**
   * The order of the quadrature rule.
   */ 
  const Order _order;

  /**
   * The type of element for which the current values have
   * been computed.
   */
  ElemType _type;

  /**
   * The reference element locations of the
   * quadrature points.
   */
  std::vector<Point> _points;

  /**
   * The value of the quadrature weights.
   */
  std::vector<real> _weights;
};



/**
 * This class implemenets specific orders of Gauss quadrature.
 * Gauss quadrature rules of order \p p have the property of
 * integrating polynomials of degree \p 2p-1 exactly.
 */

// ------------------------------------------------------------
// QGauss class definition

class QGauss : public QBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QGauss (const unsigned int _dim,
	  const Order _order=INVALID_ORDER);

  /**
   * Destructor.
   */
  ~QGauss();

  /**
   * Initializes the data structures to contain a quadrature rule
   * for an object of type \p _type.
   */
  void init(const ElemType _type=INVALID_ELEM);
  
  /**
   * Initializes the data structures to contain a quadrature rule
   * for side \p side of object of type \p _type.
   */
  void init (const ElemType _type,
	     const unsigned int side);
 private:

  void init_1D (const ElemType _type=INVALID_ELEM);
  void init_2D (const ElemType _type=INVALID_ELEM);
  void init_3D (const ElemType _type=INVALID_ELEM);
  
  void init_2D (const ElemType _type,
		const unsigned int side);
  void init_3D (const ElemType _type,
		const unsigned int side);

};



/**
 * This class implemenets trapezoidal quadratue.  These rules
 * sample at the corners and will integrate linears exactly.
 */

// ------------------------------------------------------------
// QTrap class definition

class QTrap : public QBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QTrap (const unsigned int _dim);

  /**
   * Destructor.
   */
  ~QTrap();

  /**
   * Initializes the data structures to contain a quadrature rule
   * for an object of type \p type.
   */
  void init(const ElemType _type=INVALID_ELEM);
  
  /**
   * Initializes the data structures to contain a quadrature rule
   * for side \p side of object of type \p _type.
   */
  void init (const ElemType _type,
	     const unsigned int side);
 private:

  void init_1D (const ElemType _type=INVALID_ELEM);
  void init_2D (const ElemType _type=INVALID_ELEM);
  void init_3D (const ElemType _type=INVALID_ELEM);
  
  void init_2D (const ElemType _type,
		const unsigned int side);
  void init_3D (const ElemType _type,
		const unsigned int side);

};



// ------------------------------------------------------------
// QBase class members

inline
QBase::QBase(const unsigned int d,
	     const Order o) :
  _dim(d),
  _order(o),
  _type(INVALID_ELEM)
{
};




inline
void QBase::print_info() const
{
  assert(!_points.empty());
  assert(!_weights.empty());

  std::cout << "N_Q_Points=" << n_points() << std::endl << std::endl;
  for (unsigned int qp=0; qp<n_points(); qp++)
    {
      std::cout << " Point " << qp << ":" << std::endl << "  ";
      _points[qp].print();
      std::cout << " Weight: " << std::endl;
      std::cout << "  w=" << _weights[qp] << std::endl << std::endl;
    };
};


// ------------------------------------------------------------
// QGauss class members
inline
QGauss::QGauss(const unsigned int d,
	       const Order o) : QBase(d,o)
{
  // explicitly call the init function in 1D since the
  // other tensor-product rules require this one.
  // note that EDGE will not be used internally, however
  // if we called the function with INVALID_ELEM it would try to
  // be smart and return, thinking it had already done the work.
  if (_dim == 1)
    init(EDGE2);
};




inline
QGauss::~QGauss()
{
};




// ------------------------------------------------------------
// QTrap class members
inline
QTrap::QTrap(const unsigned int d) : QBase(d,FIRST)
{
  // explicitly call the init function in 1D since the
  // other tensor-product rules require this one.
  // note that EDGE will not be used internally, however
  // if we called the function with INVALID_ELEM it would try to
  // be smart and return, thinking it had already done the work.
  if (_dim == 1)
    init(EDGE2);
};




inline
QTrap::~QTrap()
{
};



#endif
