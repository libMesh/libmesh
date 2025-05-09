// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/enum_elem_type.h" // INVALID_ELEM
#include "libmesh/enum_order.h" // INVALID_ORDER

// C++ includes
#include <vector>
#include <string>
#include <utility>
#include <memory>

// We used to have additional arguments to the init_0D, init_1D, etc.
// classes, but they became unused after QBase members made them
// redundant, so they were deprecated and eventually removed.  Users'
// derived classes that need to compile against both old and new
// libMesh version can test this macro to determine what signature
// their overrides should have.
#define LIBMESH_QBASE_INIT_ARGUMENTS_REMOVED

namespace libMesh
{

// forward declarations
class Elem;
enum QuadratureType : int;

/**
 * The \p QBase class provides the basic functionality from which
 * various quadrature rules can be derived.  It computes and stores
 * the quadrature points (in reference element space) and associated
 * weights.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief Base class for all quadrature families and orders.
 */
class QBase : public ReferenceCountedObject<QBase>
{
protected:

  /**
   * Constructor. Protected to prevent instantiation of this base
   * class.  Use the build() method instead.
   */
  QBase (unsigned int dim,
         Order order=INVALID_ORDER);

public:

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this simple class.
   */
  QBase (const QBase &) = default;
  QBase (QBase &&) = default;
  QBase & operator= (const QBase &) = default;
  QBase & operator= (QBase &&) = default;
  virtual ~QBase() = default;

  /**
   * \returns A copy of this quadrature rule wrapped in a smart pointer.
   */
  virtual std::unique_ptr<QBase> clone() const;

  /**
   * \returns The quadrature type in derived classes.
   */
  virtual QuadratureType type() const = 0;

  /**
   * Builds a specific quadrature rule based on the \p name
   * string. This enables selection of the quadrature rule at
   * run-time.  The input parameter \p name must be mappable through
   * the \p Utility::string_to_enum<>() function.
   *
   * This function allocates memory, therefore a \p std::unique_ptr<QBase>
   * is returned so that the user does not accidentally leak it.
   */
  static std::unique_ptr<QBase> build (std::string_view name,
                                       const unsigned int dim,
                                       const Order order=INVALID_ORDER);

  /**
   * Builds a specific quadrature rule based on the QuadratureType.
   * This enables selection of the quadrature rule at run-time.
   *
   * This function allocates memory, therefore a \p std::unique_ptr<QBase>
   * is returned so that the user does not accidentally leak it.
   */
  static std::unique_ptr<QBase> build (const QuadratureType qt,
                                       const unsigned int dim,
                                       const Order order=INVALID_ORDER);

  /**
   * \returns The element type we're currently using.
   */
  ElemType get_elem_type() const { return _type; }

  /**
   * \returns The p-refinement level we're currently using.
   */
  unsigned int get_p_level() const { return _p_level; }

  /**
   * \returns The number of points associated with the quadrature rule.
   */
  unsigned int n_points() const
  {
    libmesh_assert (!_points.empty());
    return cast_int<unsigned int>(_points.size());
  }

  /**
   * Alias for n_points() to enable use in index_range
   *
   * \returns The number of points associated with the quadrature rule.
   */
  unsigned int size() const
  {
    return n_points();
  }

  /**
   * \returns The spatial dimension of the quadrature rule.
   */
  unsigned int get_dim() const { return _dim; }

  /**
   * \returns A \p std::vector containing the quadrature point locations
   * in reference element space.
   */
  const std::vector<Point> & get_points() const { return _points; }

  /**
   * \returns A \p std::vector containing the quadrature point locations
   * in reference element space as a writable reference.
   */
  std::vector<Point> & get_points() { return _points; }

  /**
   * \returns A constant reference to a \p std::vector containing the
   * quadrature weights.
   */
  const std::vector<Real> & get_weights() const { return _weights; }

  /**
   * \returns A writable references to a \p std::vector containing the
   * quadrature weights.
   */
  std::vector<Real> & get_weights() { return _weights; }

  /**
   * \returns The \f$ i^{th} \f$ quadrature point in reference element space.
   */
  Point qp(const unsigned int i) const
  {
    libmesh_assert_less (i, _points.size());
    return _points[i];
  }

  /**
   * \returns The \f$ i^{th} \f$ quadrature weight.
   */
  Real w(const unsigned int i) const
  {
    libmesh_assert_less (i, _weights.size());
    return _weights[i];
  }

  /**
   * Initializes the data structures for a quadrature rule for the
   * element \p e.  If \p p_level is specified it overrides the element
   * p_level() elevation to use.
   */
  virtual void init (const Elem & e,
                     unsigned int p_level=invalid_uint);

  /**
   * Initializes the data structures for a quadrature rule for an
   * element of type \p type.  Some types, such as Polygon subclasses,
   * might require more detailed element information and so might not
   * be compatible with this API.
   *
   * New code should use the Elem-based API for most use cases, but
   * some code may initialize quadrature rules with simple ElemType
   * values like triangles and edges for use in tensor product or
   * conical product constructions; this code can set the
   * \p simple_type_only flag to avoid being identified as deprecated.
   */
  virtual void init (const ElemType type=INVALID_ELEM,
                     unsigned int p_level=0,
                     bool simple_type_only=false);

  /**
   * Initializes the data structures for a quadrature rule based on
   * the element, element type, and p_level settings of \p other_rule
   */
  virtual void init (const QBase & other_rule);

  /**
   * Initializes the data structures for an element potentially "cut"
   * by a signed distance function.  The array \p vertex_distance_func
   * contains vertex values of the signed distance function.  If the
   * signed distance function changes sign on the vertices, then the
   * element is considered to be cut.) This interface can be extended
   * by derived classes in order to subdivide the element and construct
   * a composite quadrature rule.
   */
  virtual void init (const Elem & elem,
                     const std::vector<Real> & vertex_distance_func,
                     unsigned int p_level=0);

  /**
   * \returns The current "total" order of the quadrature rule which
   * can vary element by element, depending on the Elem::p_level(),
   * which gets passed to us during init().
   *
   * Each additional power of p increases the quadrature order
   * required to integrate the mass matrix by 2, hence the formula
   * below.
   *
   * \todo This function should also be used in all of the Order
   * switch statements in the rules themselves.
   */
  Order get_order() const { return static_cast<Order>(_order + 2 * _p_level); }

  /**
   * \returns The "base" order of the quadrature rule, independent of
   * element.
   *
   * This function should be used when comparing quadrature objects
   * independently of their last initialization.
   */
  Order get_base_order() const { return static_cast<Order>(_order); }

  /**
   * Prints information relevant to the quadrature rule, by default to
   * libMesh::out.
   */
  void print_info(std::ostream & os=libMesh::out) const;

  /**
   * Maps the points of a 1D quadrature rule defined by "old_range" to
   * another 1D interval defined by "new_range" and scales the weights
   * accordingly.
   */
  void scale(std::pair<Real, Real> old_range,
             std::pair<Real, Real> new_range);

  /**
   * Same as above, but allows you to use the stream syntax.
   */
  friend std::ostream & operator << (std::ostream & os, const QBase & q);

  /**
   * \returns \p true if the shape functions need to be recalculated,
   * \p false otherwise.
   *
   * This may be required if the number of quadrature points or their
   * position changes.
   */
  virtual bool shapes_need_reinit() { return false; }

  /**
   * Flag (default true) controlling the use of quadrature rules with
   * negative weights.  Set this to false to require rules with all
   * positive weights.
   *
   * Rules with negative weights can be unsuitable for some problems.
   * For example, it is possible for a rule with negative weights to
   * obtain a negative result when integrating a positive function.
   *
   * A particular example: if rules with negative weights are not allowed,
   * a request for TET,THIRD (5 points) will return the TET,FIFTH (14 points)
   * rule instead, nearly tripling the computational effort required!
   */
  bool allow_rules_with_negative_weights;

  /**
   * The flag's value defaults to false so that one does not accidentally
   * use a nodal quadrature rule on Pyramid elements, since evaluating the
   * inverse element Jacobian (e.g. dphi) is not well-defined at the
   * Pyramid apex because the element Jacobian is zero there.
   *
   * We do not want to completely prevent someone from using a nodal
   * quadrature rule on Pyramids, however, since there are legitimate
   * use cases (lumped mass matrix) so the flag can be set to true to
   * override this behavior.
   */
  bool allow_nodal_pyramid_quadrature;

protected:


  /**
   * Initializes the 0D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  Generally this
   * is just one point with weight 1.
   */
  virtual void init_0D ();

  /**
   * Initializes the 1D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   * It is assumed that derived quadrature rules will at least
   * define the init_1D function, therefore it is pure virtual.
   */
  virtual void init_1D () = 0;

  /**
   * Initializes the 2D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   * Should not be pure virtual since a derived quadrature rule
   * may only be defined in 1D.  If not overridden, throws an
   * error.
   */
  virtual void init_2D ();

  /**
   * Initializes the 3D quadrature rule by filling the points and
   * weights vectors with the appropriate values.  The order of
   * the rule will be defined by the implementing class.
   * Should not be pure virtual since a derived quadrature rule
   * may only be defined in 1D.  If not overridden, throws an
   * error.
   */
  virtual void init_3D ();

  /**
   * Constructs a 2D rule from the tensor product of \p q1D with
   * itself.  Used in the \p init_2D() routines for quadrilateral
   * element types.
   */
  void tensor_product_quad (const QBase & q1D);

  /**
   * Computes the tensor product quadrature rule [q1D x q1D x q1D]
   * from the 1D rule q1D.  Used in the init_3D routines for
   * hexahedral element types.
   */
  void tensor_product_hex (const QBase & q1D);

  /**
   * Computes the tensor product of a 1D quadrature rule and a 2D
   * quadrature rule.  Used in the init_3D routines for prismatic
   * element types.
   */
  void tensor_product_prism (const QBase & q1D, const QBase & q2D);

  /**
   * The spatial dimension of the quadrature rule.
   */
  unsigned int _dim;

  /**
   * The polynomial order which the quadrature rule is capable of
   * integrating exactly.
   */
  Order _order;

  /**
   * The type of element for which the current values have been
   * computed.
   */
  ElemType _type;

  /**
   * The element for which the current values were computed, or
   * nullptr if values were computed without a specific element.
   */
  const Elem * _elem;

  /**
   * The p-level of the element for which the current values have
   * been computed.
   */
  unsigned int _p_level;

  /**
   * The locations of the quadrature points in reference element
   * space.
   */
  std::vector<Point> _points;

  /**
   * The quadrature weights.  The order of the weights matches the
   * ordering of the _points vector.
   */
  std::vector<Real> _weights;
};

} // namespace libMesh

#endif // LIBMESH_QUADRATURE_H
