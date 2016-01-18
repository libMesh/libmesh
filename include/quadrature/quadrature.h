// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_QUADRATURE_H
#define LIBMESH_QUADRATURE_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/point.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/auto_ptr.h"

// C++ includes
#include <vector>
#include <string>
#include <utility>

namespace libMesh
{


// forward declarations
class Elem;

/**
 * This is the \p QBase class.  It provides the basic functionality
 * from which various quadrature rules can be derived.  The class contains
 * \p dim dimensional points describing the quadrature locations
 * (referenced to a master object) and associated weights.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 */
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
   * \p name string.  An \p UniquePtr<QBase> is returned
   * to prevent a memory leak. This way the user need not
   * remember to delete the object.  Enables run-time decision of
   * the quadrature rule.  The input parameter \p name
   * must be mappable through the \p Utility::string_to_enum<>()
   * function.
   */
  static UniquePtr<QBase> build (const std::string & name,
                                 const unsigned int _dim,
                                 const Order _order=INVALID_ORDER);

  /**
   * Builds a specific quadrature rule, identified through the
   * \p QuadratureType.  An \p UniquePtr<QBase> is returned
   * to prevent a memory leak. This way the user need not
   * remember to delete the object.  Enables run-time decision of
   * the quadrature rule.
   */
  static UniquePtr<QBase> build (const QuadratureType _qt,
                                 const unsigned int _dim,
                                 const Order _order=INVALID_ORDER);

  /**
   * @returns the current element type we're set up for
   */
  ElemType get_elem_type() const
  { return _type; }

  /**
   * @returns the current p refinement level we're initialized with
   */
  unsigned int get_p_level() const
  { return _p_level; }

  /**
   * @returns the number of points associated with the quadrature rule.
   */
  unsigned int n_points() const
  { libmesh_assert (!_points.empty());
    return cast_int<unsigned int>(_points.size()); }

  /**
   * @returns the dimension of the quadrature rule.
   */
  unsigned int get_dim() const { return _dim;  }

  /**
   * @returns a \p std::vector containing the quadrature point locations
   * on a reference object.
   */
  const std::vector<Point> & get_points() const { return _points;  }

  /**
   * @returns a \p std::vector containing the quadrature point locations
   * on a reference object as a writeable reference.
   */
  std::vector<Point> & get_points() { return _points;  }

  /**
   * @returns a \p std::vector containing the quadrature weights.
   */
  const std::vector<Real> & get_weights() const { return _weights; }

  /**
   * @returns a \p std::vector containing the quadrature weights.
   */
  std::vector<Real> & get_weights() { return _weights; }

  /**
   * @returns the \f$ i^{th} \f$ quadrature point on the reference object.
   */
  Point qp(const unsigned int i) const
  { libmesh_assert_less (i, _points.size()); return _points[i]; }

  /**
   * @returns the \f$ i^{th} \f$ quadrature weight.
   */
  Real w(const unsigned int i) const
  { libmesh_assert_less (i, _weights.size()); return _weights[i]; }

  /**
   * Initializes the data structures to contain a quadrature rule
   * for an object of type \p type.
   */
  virtual void init (const ElemType type=INVALID_ELEM,
                     unsigned int p_level=0);

  /**
   * Initializes the data structures for a specific, potentially cut
   * element.  The array \p vertex_distance_func contains vertex
   * values of a signed distance function that cuts the element.  This
   * interface is indended to be extended by derived classes that can
   * cut the element into subelements, for example, and constuct a
   * composite quadrature rule for the cut element.
   */
  virtual void init (const Elem & elem,
                     const std::vector<Real> & vertex_distance_func,
                     unsigned int p_level=0);

  /**
   * @returns the order of the quadrature rule.
   */
  Order get_order() const { return static_cast<Order>(_order + _p_level); }

  /**
   * Prints information relevant to the quadrature rule, by default to
   * libMesh::out.
   */
  void print_info(std::ostream & os=libMesh::out) const;

  /**
   * Maps the points of a 1D interval quadrature rule (typically [-1,1])
   * to any other 1D interval (typically [0,1]) and scales the weights
   * accordingly.  The quadrature rule will be mapped from the
   * entries of old_range to the entries of new_range.
   */
  void scale(std::pair<Real, Real> old_range,
             std::pair<Real, Real> new_range);

  /**
   * Same as above, but allows you to use the stream syntax.
   */
  friend std::ostream & operator << (std::ostream & os, const QBase & q);

  /**
   * Returns true if the shape functions need to be recalculated.
   *
   * This can happen if the number of points or their positions change.
   *
   * By default this will return false.
   */
  virtual bool shapes_need_reinit() { return false; }

  /**
   * Flag (default true) controlling the use of quadrature rules with negative
   * weights.  Set this to false to ONLY use (potentially) safer but more expensive
   * rules with all positive weights.
   *
   * Negative weights typically appear in Gaussian quadrature rules
   * over three-dimensional elements.  Rules with negative weights can
   * be unsuitable for some problems.  For example, it is possible for
   * a rule with negative weights to obtain a negative result when
   * integrating a positive function.
   *
   * A particular example: if rules with negative weights are not allowed,
   * a request for TET,THIRD (5 points) will return the TET,FIFTH (14 points)
   * rule instead, nearly tripling the computational effort required!
   */
  bool allow_rules_with_negative_weights;

protected:


  /**
   * Initializes the 0D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  Generally this
   * is just one point with weight 1.
   */
  virtual void init_0D (const ElemType type=INVALID_ELEM,
                        unsigned int p_level=0);

  /**
   * Initializes the 1D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   * It is assumed that derived quadrature rules will at least
   * define the init_1D function, therefore it is pure virtual.
   */
  virtual void init_1D (const ElemType type=INVALID_ELEM,
                        unsigned int p_level=0) = 0;

  /**
   * Initializes the 2D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   * Should not be pure virtual since a derived quadrature rule
   * may only be defined in 1D.  If not redefined, gives an
   * error (when \p DEBUG defined) when called.
   */
  virtual void init_2D (const ElemType,
                        unsigned int =0)
#ifndef DEBUG
  {}
#else
  {
    libmesh_error_msg("ERROR: Seems as if this quadrature rule \nis not implemented for 2D.");
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
  virtual void init_3D (const ElemType,
                        unsigned int =0)
#ifndef DEBUG
  {}
#else
  {
    libmesh_error_msg("ERROR: Seems as if this quadrature rule \nis not implemented for 3D.");
  }
#endif


  /**
   * Computes the tensor product of
   * two 1D rules and returns a 2D rule.
   * Used in the init_2D routines for
   * quadrilateral element types.
   */
  void tensor_product_quad (const QBase & q1D);

  /**
   * Computes the tensor product quadrature rule
   * [q1D x q1D x q1D] from the 1D rule q1D.
   * Used in the init_3D routines for
   * hexahedral element types.
   */
  void tensor_product_hex (const QBase & q1D);

  /**
   * Computes the tensor product of
   * a 1D quadrature rule and a 2D
   * quadrature rule.
   * Used in the init_3D routines for
   * prismatic element types.
   */
  void tensor_product_prism (const QBase & q1D, const QBase & q2D);



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
   * The p level of element for which the current values have
   * been computed.
   */
  unsigned int _p_level;

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
  allow_rules_with_negative_weights(true),
  _dim(d),
  _order(o),
  _type(INVALID_ELEM),
  _p_level(0)
{
}




inline
void QBase::print_info(std::ostream & os) const
{
  libmesh_assert(!_points.empty());
  libmesh_assert(!_weights.empty());

  Real summed_weights=0;
  os << "N_Q_Points=" << this->n_points() << std::endl << std::endl;
  for (unsigned int qpoint=0; qpoint<this->n_points(); qpoint++)
    {
      os << " Point " << qpoint << ":\n"
         << "  "
         << _points[qpoint]
         << "\n Weight:\n "
         << "  w=" << _weights[qpoint] << "\n" << std::endl;

      summed_weights += _weights[qpoint];
    }
  os << "Summed Weights: " << summed_weights << std::endl;
}

} // namespace libMesh

#endif // LIBMESH_QUADRATURE_H
