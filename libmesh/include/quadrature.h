// $Id: quadrature.h,v 1.4 2003-01-21 19:24:35 benkirk Exp $

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
#include "enum_elem_type.h"
#include "enum_order.h"



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
  void init (const ElemType _type=INVALID_ELEM);
  
  /**
   * Initializes the data structures to contain a quadrature rule
   * for side \p side of object of type \p _type.
   */
  void init (const ElemType _type,
	     const unsigned int side);

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
   * Initializes the 1D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class. 
   */
  virtual void init_1D (const ElemType _type=INVALID_ELEM) = 0;

  /**
   * Initializes the 2D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   */
  virtual void init_2D (const ElemType _type=INVALID_ELEM) = 0;

  /**
   * Initializes the 3D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   */
  virtual void init_3D (const ElemType _type=INVALID_ELEM) = 0;


  
  /**
   * Initialize the 1D quadrature rule for a side (edge).
   */
  virtual void init_2D (const ElemType _type,
			const unsigned int side) = 0;
  /**
   * Initialize the 2D quadrature rule for a side (face).
   */
  virtual void init_3D (const ElemType _type,
			const unsigned int side) = 0;
  
  /**
   * Computes the tensor product of
   * two 1D rules and returns a 2D rule.
   * Used in the init_2D routines for
   * quadrilateral element types.
   */
  void tensor_product_quad (QBase* q1D);

  /**
   * Computes the tensor product of
   * three 1D rules and returns a 3D rule.
   * Used in the init_3D routines for
   * hexahedral element types.
   */
  void tensor_product_hex (QBase* q1D);
  
  /**
   * Computes the tensor product of
   * three 1D rules and returns a 3D rule.
   * Used in the init_3D routines for
   * hexahedral element types.
   */
  void tensor_product_prism (QBase* q1D, QBase* q2D);

  /**
   * Computes the quadrature rule for side
   * \p side of a quadrilateral element.
   * Used by the init_2D routines when passed
   * a side number.
   */
  void side_rule_quad (QBase* q1D, unsigned int side);

  /**
   * Computes the quadrature rule for side
   * \p side of a triangular element.
   * Used by the init_2D routines when passed
   * a side number.
   */
  void side_rule_tri (QBase* q1D, unsigned int side);

  /**
   * Computes the quadrature rule for side \p side
   * of a hexahedral element.  (The sides of hexes
   * are of course quads.)  Used by the init_3D
   * routines when passed a side number;
   */
  void side_rule_hex (QBase* q2D, unsigned int side);

  /**
   * Computes the quadrature rule for side \p side
   * of a tetrahedral element. (The sides of tets
   * are of course tris.)  Used by the init_3D routines
   * when passed a side number.
   */
  void side_rule_tet (QBase* q2D, unsigned int side);

  /**
   * Computes a quadrature rule for side \p side
   * of a prismatic element.  (The sides of prisms
   * can be either tris or quads.)  There is no need
   * to call init on q2D before calling this routine,
   * although it shouldn't hurt.
   */
  void side_rule_prism (QBase* q2D, unsigned int side);

  /**
   * Computes a quadrature rule for side \p side of
   * a pyramid.  (The sides of pyramids are either
   * tris or quads.)  There is no need to call init
   * on q2D before calling this routine, although
   * it shouldn't hurt.
   */
  void side_rule_pyramid (QBase* q2D, unsigned int side);
  
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

#endif
