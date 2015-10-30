// $Id: quadrature.h,v 1.2 2004-01-03 15:37:42 benkirk Exp $

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



#ifndef __quadrature_h__
#define __quadrature_h__

// C++ includes
#include <vector>
#include <utility>

// Local includes
#include "libmesh_common.h"
#include "reference_counted_object.h"
#include "point.h"
#include "enum_elem_type.h"
#include "enum_order.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"



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

class QBase : public ReferenceCountedObject<QBase>
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
  virtual ~QBase() {}

  /**
   * @returns the quadrature type in derived classes.
   */
  virtual QuadratureType type() const = 0;

  /**
   * Builds a specific quadrature rule, identified through the
   * \p QuadratureType.  An \p AutoPtr<QBase> is returned
   * to prevent a memory leak. This way the user need not
   * remember to delete the object.  Enables run-time decision of
   * the quadrature rule.
   */
  static AutoPtr<QBase> build(const QuadratureType _qt,
			      const unsigned int _dim,
			      const Order _order=INVALID_ORDER);

  /**
   * @returns the number of points associated with the quadrature rule.
   */    
  unsigned int n_points() const
    { assert (!_points.empty()); return _points.size(); }

  /**
   * @returns a \p std::vector containing the quadrature point locations
   * on a reference object.
   */
  const std::vector<Point>& get_points() const { return _points;  }

  /**
   * @returns a \p std::vector containing the quadrature point locations
   * on a reference object as a writeable reference.
   */
  std::vector<Point>& get_points() { return _points;  }

  /**
   * @returns a \p std::vector containing the quadrature weights.
   */
  const std::vector<Real>& get_weights() const { return _weights; }

  /**
   * @returns a \p std::vector containing the quadrature weights.
   */
  std::vector<Real>& get_weights() { return _weights; }

  /**
   * @returns the \f$ i^{th} \f$ quadrature point on the reference object.
   */
  Point qp(const unsigned int i) const
    { assert (i < _points.size()); return _points[i]; }

  /**
   * @returns the \f$ i^{th} \f$ quadrature weight.
   */
  Real w(const unsigned int i) const
    { assert (i < _weights.size()); return _weights[i]; }
  
  /**
   * Initializes the data structures to contain a quadrature rule
   * for an object of type \p type.  
   */
  void init (const ElemType _type=INVALID_ELEM);

  /**
   * @returns the order of the quadrature rule.   
   */
  Order get_order() const { return _order; }
  
  /**
   * Prints information relevant to the quadrature rule.
   */
  void print_info() const;



protected:

  
  /**
   * Initializes the 1D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   * It is assumed that derived quadrature rules will at least
   * define the init_1D function, therefore it is pure virtual.
   */
  virtual void init_1D (const ElemType _type=INVALID_ELEM) = 0;

  /**
   * Initializes the 2D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   * Should not be pure virtual since a derived quadrature rule
   * may only be defined in 1D.  If not redefined, gives an
   * error (when \p DEBUG defined) when called.
   */
  virtual void init_2D (const ElemType)
#ifndef DEBUG
  {}
#else
  {  
    std::cerr << "ERROR: Seems as if this quadrature rule" << std::endl
	      << " is not implemented for 2D." << std::endl;
    error();
  }
#endif

  /**
   * Initializes the 3D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   * Should not be pure virtual since a derived quadrature rule
   * may only be defined in 1D.  If not redefined, gives an
   * error (when \p DEBUG defined) when called.
   */
  virtual void init_3D (const ElemType)
#ifndef DEBUG
  {}
#else
  {  
    std::cerr << "ERROR: Seems as if this quadrature rule" << std::endl
	      << " is not implemented for 3D." << std::endl;
    error();
  }
#endif
  
  /**
   * Maps the points of a 1D interval quadrature rule (typically [-1,1])
   * to any other 1D interval (typically [0,1]) and scales the weights
   * accordingly.  The quadrature rule will be mapped from the
   * entries of old_range to the entries of new_range.
   */
  void scale(std::pair<Real, Real> old_range,
	     std::pair<Real, Real> new_range);
  
  /**
   * Computes the tensor product of
   * two 1D rules and returns a 2D rule.
   * Used in the init_2D routines for
   * quadrilateral element types.
   */
  void tensor_product_quad (QBase& q1D);

  /**
   * Computes the conical product of
   * two 1D rules to generate a (sub-optimal)
   * 2D rule for triangles.  Note that:
   * gauss1D = 1D Gauss rule
   * jacA1D  = 1D Jacobi-Gauss rule with (1-x) wt. funtion
   * Method can be found in:
   * Approximate Calculation of Multiple Integrals, Stroud, A. H.
   */
  void tensor_product_tri (QBase& gauss1D, QBase& jacA1D);
  
  /**
   * Computes the tensor product quadrature rule
   * [q1D x q1D x q1D] from the 1D rule q1D.
   * Used in the init_3D routines for
   * hexahedral element types.
   */
  void tensor_product_hex (QBase& q1D);
  
  /**
   * Computes the tensor product of
   * a 1D quadrature rule and a 2D
   * quadrature rule.
   * Used in the init_3D routines for
   * prismatic element types.
   */
  void tensor_product_prism (QBase& q1D, QBase& q2D);

  /**
   * Computes the conical product of
   * three 1D rules to generate a (sub-optimal)
   * 3D rule for tets.  Note that:
   * gauss1D = 1D Gauss rule
   * jacA1D  = 1D Jacobi-Gauss rule with (1-x) wt. funtion
   * jacB1D  = 1D Jacobi-Gauss rule with (1-x)^2 wt. function 
   * Method can be found in:
   * Approximate Calculation of Multiple Integrals, Stroud, A. H.
   */
  void tensor_product_tet (QBase& gauss1D, QBase& jacA1D, QBase& jacB1D);

  
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
  std::vector<Real> _weights;
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
}




inline
void QBase::print_info() const
{
  assert(!_points.empty());
  assert(!_weights.empty());

  std::cout << "N_Q_Points=" << this->n_points() << std::endl << std::endl;
  for (unsigned int qp=0; qp<this->n_points(); qp++)
    {
      std::cout << " Point " << qp << ":" << std::endl << "  ";
      _points[qp].print();
      std::cout << " Weight: " << std::endl;
      std::cout << "  w=" << _weights[qp] << std::endl << std::endl;
    }
}

#endif
