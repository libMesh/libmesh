// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FE_H
#define LIBMESH_FE_H

// Local includes
#include "libmesh/fe_base.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// forward declarations
class DofConstraints;
class DofMap;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
class InfFE;

#endif


/**
 * Most finite element types in libMesh are scalar-valued
 */
template <FEFamily T>
struct FEOutputType
{
  typedef Real type;
};


/**
 * Specialize for non-scalar-valued elements
 */
template<>
struct FEOutputType<LAGRANGE_VEC>
{
  typedef RealVectorValue type;
};

template<>
struct FEOutputType<NEDELEC_ONE>
{
  typedef RealVectorValue type;
};

template<>
struct FEOutputType<MONOMIAL_VEC>
{
  typedef RealVectorValue type;
};


/**
 * A specific instantiation of the \p FEBase class. This
 * class is templated, and specific template instantiations
 * will result in different Finite Element families. Full specialization
 * of the template for specific dimensions(\p Dim) and families
 * (\p T) provide support for specific finite element types.
 * The use of templates allows for compile-time optimization,
 * however it requires that the specific finite element family
 * and dimension is also known at compile time.  If this is
 * too restricting for your application you can use the
 * \p FEBase::build() member to create abstract (but still optimized)
 * finite elements.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 * \brief Template class which generates the different FE families and orders.
 */
template <unsigned int Dim, FEFamily T>
class FE : public FEGenericBase<typename FEOutputType<T>::type>
{
public:

  /**
   * Constructor.
   */
  explicit
  FE(const FEType & fet);

  typedef typename
  FEGenericBase<typename FEOutputType<T>::type>::OutputShape
  OutputShape;

  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method allows you to specify the dimension,
   * element type, and order directly.  This allows the method to
   * be static.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static OutputShape shape(const ElemType t,
                           const Order o,
                           const unsigned int i,
                           const Point & p);

  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method allows you to specify the dimension,
   * element type, and order directly.  This allows the method to
   * be static.
   *
   * On a p-refined element, \p o should be the base order of the
   * element if \p add_p_level is left \p true, or can be the base
   * order of the element if \p add_p_level is set to \p false.
   */
  static OutputShape shape(const Elem * elem,
                           const Order o,
                           const unsigned int i,
                           const Point & p,
                           const bool add_p_level = true);

  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method allows you to specify the dimension and
   * element type directly. The order is given by the FEType.
   * This allows the method to be static.
   *
   * On a p-refined element, \p o should be the base order of the
   * element if \p add_p_level is left \p true, or can be the base
   * order of the element if \p add_p_level is set to \p false.
   */
   static OutputShape shape(const FEType fet,
                            const Elem * elem,
                            const unsigned int i,
                            const Point & p,
                            const bool add_p_level = true);

  /**
   * Fills \p v with the values of the \f$ i^{th} \f$
   * shape function, evaluated at all points p.  You must specify
   * element order directly.  \p v should already be the appropriate
   * size.
   *
   * On a p-refined element, \p o should be the base order of the
   * element if \p add_p_level is left \p true, or can be the base
   * order of the element if \p add_p_level is set to \p false.
   */
  static void shapes(const Elem * elem,
                     const Order o,
                     const unsigned int i,
                     const std::vector<Point> & p,
                     std::vector<OutputShape> & v,
                     const bool add_p_level = true);


  /**
   * Fills \p v[i][qp] with the values of the \f$ i^{th} \f$
   * shape functions, evaluated at points qp in p.  You must specify
   * element order directly.  \p v should already be the appropriate
   * size.
   *
   * On a p-refined element, \p o should be the base order of the
   * element if \p add_p_level is left \p true, or can be the base
   * order of the element if \p add_p_level is set to \p false.
   */
  static void all_shapes(const Elem * elem,
                         const Order o,
                         const std::vector<Point> & p,
                         std::vector<std::vector<OutputShape> > & v,
                         const bool add_p_level = true);

  /**
   * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at point \p p.  This method allows you to
   * specify the dimension, element type, and order directly.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static OutputShape shape_deriv(const ElemType t,
                                 const Order o,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p);

  /**
   * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function.  You must specify element type, and order directly.
   *
   * On a p-refined element, \p o should be the base order of the
   * element if \p add_p_level is left \p true, or can be the base
   * order of the element if \p add_p_level is set to \p false.
   */
  static OutputShape shape_deriv(const Elem * elem,
                                 const Order o,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level = true);

  /**
   * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function.  You must specify element type, and order (via
   * FEType) directly.
   *
   * On a p-refined element, \p o should be the base order of the
   * element if \p add_p_level is left \p true, or can be the base
   * order of the element if \p add_p_level is set to \p false.
   */
  static OutputShape shape_deriv(const FEType fet,
                                 const Elem * elem,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p,
                                 const bool add_p_level = true);

  /**
   * Fills \p v with the \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function, evaluated at all points p.  You must specify
   * element order directly.  \p v should already be the appropriate
   * size.
   *
   * On a p-refined element, \p o should be the base order of the
   * element if \p add_p_level is left \p true, or can be the base
   * order of the element if \p add_p_level is set to \p false.
   */
  static void shape_derivs(const Elem * elem,
                           const Order o,
                           const unsigned int i,
                           const unsigned int j,
                           const std::vector<Point> & p,
                           std::vector<OutputShape> & v,
                           const bool add_p_level = true);


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * \returns The second \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at the point \p p.
   *
   * \note Cross-derivatives are indexed according to:
   * j = 0 ==> d^2 phi / dxi^2
   * j = 1 ==> d^2 phi / dxi deta
   * j = 2 ==> d^2 phi / deta^2
   * j = 3 ==> d^2 phi / dxi dzeta
   * j = 4 ==> d^2 phi / deta dzeta
   * j = 5 ==> d^2 phi / dzeta^2
   *
   * \note Computing second derivatives is not currently supported for
   * all element types: \f$ C^1 \f$ (Clough, Hermite and Subdivision),
   * Lagrange, Hierarchic, L2_Hierarchic, and Monomial are supported.
   * All other element types return an error when asked for second
   * derivatives.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static OutputShape shape_second_deriv(const ElemType t,
                                        const Order o,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p);

  /**
   * \returns The second \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at the point \p p.
   *
   * \note Cross-derivatives are indexed according to:
   * j = 0 ==> d^2 phi / dxi^2
   * j = 1 ==> d^2 phi / dxi deta
   * j = 2 ==> d^2 phi / deta^2
   * j = 3 ==> d^2 phi / dxi dzeta
   * j = 4 ==> d^2 phi / deta dzeta
   * j = 5 ==> d^2 phi / dzeta^2
   *
   * \note Computing second derivatives is not currently supported for
   * all element types: \f$ C^1 \f$ (Clough, Hermite and Subdivision),
   * Lagrange, Hierarchic, L2_Hierarchic, and Monomial are supported.
   * All other element types return an error when asked for second
   * derivatives.
   *
   * On a p-refined element, \p o should be the base order of the
   * element if \p add_p_level is left \p true, or can be the base
   * order of the element if \p add_p_level is set to \p false.
   */
  static OutputShape shape_second_deriv(const Elem * elem,
                                        const Order o,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level = true);

  /**
   * \returns The second \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at the point \p p.
   *
   * \note Cross-derivatives are indexed according to:
   * j = 0 ==> d^2 phi / dxi^2
   * j = 1 ==> d^2 phi / dxi deta
   * j = 2 ==> d^2 phi / deta^2
   * j = 3 ==> d^2 phi / dxi dzeta
   * j = 4 ==> d^2 phi / deta dzeta
   * j = 5 ==> d^2 phi / dzeta^2
   *
   * \note Computing second derivatives is not currently supported for
   * all element types: \f$ C^1 \f$ (Clough, Hermite and Subdivision),
   * Lagrange, Hierarchic, L2_Hierarchic, and Monomial are supported.
   * All other element types return an error when asked for second
   * derivatives.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static OutputShape shape_second_deriv(const FEType fet,
                                        const Elem * elem,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p,
                                        const bool add_p_level = true);

#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Build the nodal soln from the element soln.
   * This is the solution that will be plotted.
   *
   * On a p-refined element, \p o should be the base order of the element.
   */
  static void nodal_soln(const Elem * elem, const Order o,
                         const std::vector<Number> & elem_soln,
                         std::vector<Number> & nodal_soln);

  /**
   * \returns The number of shape functions associated with
   * this finite element.
   */
  virtual unsigned int n_shape_functions () const override;

  /**
   * \returns The number of shape functions associated with
   * a finite element of type \p t and approximation order \p o.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static unsigned int n_shape_functions (const ElemType t,
                                         const Order o)
  { return FE<Dim,T>::n_dofs (t,o); }

  /**
   * \returns The number of shape functions associated with this
   * finite element.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static unsigned int n_dofs(const ElemType t,
                             const Order o);

  /**
   * \returns The number of dofs at node \p n for a finite element
   * of type \p t and order \p o.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static unsigned int n_dofs_at_node(const ElemType t,
                                     const Order o,
                                     const unsigned int n);

  /**
   * \returns The number of dofs interior to the element,
   * not associated with any interior nodes.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static unsigned int n_dofs_per_elem(const ElemType t,
                                      const Order o);

  /**
   * \returns The continuity level of the finite element.
   */
  virtual FEContinuity get_continuity() const override;

  /**
   * \returns \p true if the finite element's higher order shape functions are
   * hierarchic
   */
  virtual bool is_hierarchic() const override;

  /**
   * Fills the vector di with the local degree of freedom indices
   * associated with side \p s of element \p elem
   *
   * On a p-refined element, \p o should be the base order of the element.
   */
  static void dofs_on_side(const Elem * const elem,
                           const Order o,
                           unsigned int s,
                           std::vector<unsigned int> & di);
  /**
   * Fills the vector di with the local degree of freedom indices
   * associated with edge \p e of element \p elem
   *
   * On a p-refined element, \p o should be the base order of the element.
   */
  static void dofs_on_edge(const Elem * const elem,
                           const Order o,
                           unsigned int e,
                           std::vector<unsigned int> & di);

  static Point inverse_map (const Elem * elem,
                            const Point & p,
                            const Real tolerance = TOLERANCE,
                            const bool secure = true) {
    // libmesh_deprecated(); // soon
    return FEMap::inverse_map(Dim, elem, p, tolerance, secure, secure);
  }

  static void inverse_map (const Elem * elem,
                           const std::vector<Point> & physical_points,
                           std::vector<Point> &       reference_points,
                           const Real tolerance = TOLERANCE,
                           const bool secure = true) {
    // libmesh_deprecated(); // soon
    FEMap::inverse_map(Dim, elem, physical_points, reference_points,
                       tolerance, secure, secure);
  }

  /**
   * This is at the core of this class. Use this for each
   * new element in the mesh.  Reinitializes all the physical
   * element-dependent data based on the current element
   * \p elem.  By default the shape functions and associated
   * data are computed at the quadrature points specified
   * by the quadrature rule \p qrule, but may be any points
   * specified on the reference element specified in the optional
   * argument \p pts.
   */
  virtual void reinit (const Elem * elem,
                       const std::vector<Point> * const pts = nullptr,
                       const std::vector<Real> * const weights = nullptr) override;

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of \p face.  The \p tolerance parameter is passed to
   * the involved call to \p inverse_map().  By default the shape
   * functions and associated data are computed at the quadrature
   * points specified by the quadrature rule \p qrule, but may be any
   * points specified on the reference \em side element specified in
   * the optional argument \p pts.
   */
  virtual void reinit (const Elem * elem,
                       const unsigned int side,
                       const Real tolerance = TOLERANCE,
                       const std::vector<Point> * const pts = nullptr,
                       const std::vector<Real> * const weights = nullptr) override;

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p edge.  The \p tolerance parameter is passed to the
   * involved call to \p inverse_map().  By default the shape
   * functions and associated data are computed at the quadrature
   * points specified by the quadrature rule \p qrule, but may be any
   * points specified on the reference \em side element specified in
   * the optional argument \p pts.
   */
  virtual void edge_reinit (const Elem * elem,
                            const unsigned int edge,
                            const Real tolerance = TOLERANCE,
                            const std::vector<Point> * const pts = nullptr,
                            const std::vector<Real> * const weights = nullptr) override;

  /**
   * Computes the reference space quadrature points on the side of
   * an element based on the side quadrature points.
   */
  virtual void side_map (const Elem * elem,
                         const Elem * side,
                         const unsigned int s,
                         const std::vector<Point> & reference_side_points,
                         std::vector<Point> &       reference_points) override;

  /**
   * Provides the class with the quadrature rule, which provides the
   * locations (on a reference element) where the shape functions are
   * to be calculated.
   */
  virtual void attach_quadrature_rule (QBase * q) override;

  /**
   * \returns The total number of quadrature points.  Call this
   * to get an upper bound for the \p for loop in your simulation
   * for matrix assembly of the current element.
   */
  virtual unsigned int n_quadrature_points () const override;

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Computes the constraint matrix contributions (for
   * non-conforming adapted meshes) corresponding to
   * variable number \p var_number, using element-specific
   * optimizations if possible.
   */
  static void compute_constraints (DofConstraints & constraints,
                                   DofMap & dof_map,
                                   const unsigned int variable_number,
                                   const Elem * elem);
#endif // #ifdef LIBMESH_ENABLE_AMR

  /**
   * \returns \p true when the shape functions (for
   * this \p FEFamily) depend on the particular
   * element, and therefore needs to be re-initialized
   * for each new element.  \p false otherwise.
   */
  virtual bool shapes_need_reinit() const override;

  static Point map (const Elem * elem,
                    const Point & reference_point) {
    // libmesh_deprecated(); // soon
    return FEMap::map(Dim, elem, reference_point);
  }

  static Point map_xi (const Elem * elem,
                       const Point & reference_point) {
    // libmesh_deprecated(); // soon
    return FEMap::map_deriv(Dim, elem, 0, reference_point);
  }

  static Point map_eta (const Elem * elem,
                        const Point & reference_point) {
    // libmesh_deprecated(); // soon
    return FEMap::map_deriv(Dim, elem, 1, reference_point);
  }

  static Point map_zeta (const Elem * elem,
                         const Point & reference_point) {
    // libmesh_deprecated(); // soon
    return FEMap::map_deriv(Dim, elem, 2, reference_point);
  }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  /**
   * make InfFE classes friends, so that these may access
   * the private \p map, map_xyz methods
   */
  template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
  friend class InfFE;
#endif

protected:

  /**
   * Update the various member data fields \p phi,
   * \p dphidxi, \p dphideta, \p dphidzeta, etc.
   * for the current element.  These data will be computed
   * at the points \p qp, which are generally (but need not be)
   * the quadrature points.
   */
  virtual void init_shape_functions(const std::vector<Point> & qp,
                                    const Elem * e);

  /**
   * Init \p dual_phi and potentially \p dual_dphi, \p dual_d2phi
   */
  void init_dual_shape_functions(unsigned int n_shapes, unsigned int n_qp);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * Initialize the data fields for the base of an
   * an infinite element.
   */
  virtual void init_base_shape_functions(const std::vector<Point> & qp,
                                         const Elem * e) override;

#endif

  /**
   * A default implementation for shapes
   */
  static void default_shapes (const Elem * elem,
                              const Order o,
                              const unsigned int i,
                              const std::vector<Point> & p,
                              std::vector<OutputShape> & v,
                              const bool add_p_level = true)
    {
      libmesh_assert_equal_to(p.size(), v.size());
      for (auto vi : index_range(v))
        v[vi] = FE<Dim,T>::shape (elem, o, i, p[vi], add_p_level);
    }

  /**
   * A default implementation for all_shapes
   */
  static void default_all_shapes (const Elem * elem,
                                  const Order o,
                                  const std::vector<Point> & p,
                                  std::vector<std::vector<OutputShape>> & v,
                                  const bool add_p_level = true)
    {
      for (auto i : index_range(v))
        {
          libmesh_assert_equal_to ( p.size(), v[i].size() );
          FE<Dim,T>::shapes (elem, o, i, p, v[i], add_p_level);
        }
    }

  /**
   * A default implementation for shape_derivs
   */
  static void default_shape_derivs (const Elem * elem,
                                    const Order o,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const std::vector<Point> & p,
                                    std::vector<OutputShape> & v,
                                    const bool add_p_level = true)
    {
      libmesh_assert_equal_to(p.size(), v.size());
      for (auto vi : index_range(v))
        v[vi] = FE<Dim,T>::shape_deriv (elem, o, i, j, p[vi], add_p_level);
    }

  /**
   * An array of the node locations on the last
   * element we computed on
   */
  std::vector<Point> cached_nodes;

  /**
   * The last side and last edge we did a reinit on
   */
  ElemType last_side;

  unsigned int last_edge;
};



/**
 * Clough-Tocher finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Roy Stogner
 * \date 2004
 */
template <unsigned int Dim>
class FEClough : public FE<Dim,CLOUGH>
{
public:

  /**
   * Constructor. Creates a hierarchic finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FEClough(const FEType & fet) :
    FE<Dim,CLOUGH> (fet)
  {}
};



/**
 * Hermite finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Roy Stogner
 * \date 2005
 */
template <unsigned int Dim>
class FEHermite : public FE<Dim,HERMITE>
{
public:

  /**
   * Constructor. Creates a hierarchic finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FEHermite(const FEType & fet) :
    FE<Dim,HERMITE> (fet)
  {}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * 1D hermite functions on unit interval
   */
  static Real hermite_raw_shape_second_deriv(const unsigned int basis_num,
                                             const Real xi);
#endif
  static Real hermite_raw_shape_deriv(const unsigned int basis_num,
                                      const Real xi);
  static Real hermite_raw_shape(const unsigned int basis_num,
                                const Real xi);
};



/**
 * Subdivision finite elements.
 *
 * Template specialization prototypes are needed for calling from
 * inside FESubdivision::init_shape_functions
 */
template <>
Real FE<2,SUBDIVISION>::shape(const Elem * elem,
                              const Order order,
                              const unsigned int i,
                              const Point & p,
                              const bool add_p_level);

template <>
Real FE<2,SUBDIVISION>::shape_deriv(const Elem * elem,
                                    const Order order,
                                    const unsigned int i,
                                    const unsigned int j,
                                    const Point & p,
                                    const bool add_p_level);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <>
Real FE<2,SUBDIVISION>::shape_second_deriv(const Elem * elem,
                                           const Order order,
                                           const unsigned int i,
                                           const unsigned int j,
                                           const Point & p,
                                           const bool add_p_level);

#endif

class FESubdivision : public FE<2,SUBDIVISION>
{
public:

  /**
   * Constructor. Creates a subdivision surface finite element.
   * Currently only supported for two-dimensional meshes in
   * three-dimensional space.
   */
  FESubdivision(const FEType & fet);

  /**
   * This is at the core of this class. Use this for each new
   * non-ghosted element in the mesh.  Reinitializes all the physical
   * element-dependent data based on the current element
   * \p elem.  By default the shape functions and associated
   * data are computed at the quadrature points specified
   * by the quadrature rule \p qrule, but may be any points
   * specified on the reference element specified in the optional
   * argument \p pts.
   */
  virtual void reinit (const Elem * elem,
                       const std::vector<Point> * const pts = nullptr,
                       const std::vector<Real> * const weights = nullptr) override;

  /**
   * This prevents some compilers being confused by partially
   * overriding this virtual function.
   */
  virtual void reinit (const Elem *,
                       const unsigned int,
                       const Real = TOLERANCE,
                       const std::vector<Point> * const = nullptr,
                       const std::vector<Real> * const = nullptr) override
  { libmesh_not_implemented(); }

  /**
   * Provides the class with the quadrature rule, which provides the
   * locations (on a reference element) where the shape functions are
   * to be calculated.
   */
  virtual void attach_quadrature_rule (QBase * q) override;

  /**
   * Update the various member data fields \p phi,
   * \p dphidxi, \p dphideta, \p dphidzeta, etc.
   * for the current element.  These data will be computed
   * at the points \p qp, which are generally (but need not be)
   * the quadrature points.
   */
  virtual void init_shape_functions(const std::vector<Point> & qp,
                                    const Elem * elem) override;

  /**
   * \returns The value of the \f$ i^{th} \f$ of the 12 quartic
   * box splines interpolating a regular Loop subdivision
   * element, evaluated at the barycentric coordinates \p v,
   * \p w.
   */
  static Real regular_shape(const unsigned int i,
                            const Real v,
                            const Real w);

  /**
   * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th}
   * \f$ of the 12 quartic box splines interpolating a regular
   * Loop subdivision element, evaluated at the barycentric
   * coordinates \p v, \p w.
   */
  static Real regular_shape_deriv(const unsigned int i,
                                  const unsigned int j,
                                  const Real v,
                                  const Real w);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * \returns The second \f$ j^{th} \f$ derivative of the
   * \f$ i^{th} \f$ of the 12 quartic box splines interpolating
   * a regular Loop subdivision element, evaluated at the
   * barycentric coordinates \p v, \p w.
   */
  static Real regular_shape_second_deriv(const unsigned int i,
                                         const unsigned int j,
                                         const Real v,
                                         const Real w);


#endif // LIBMESH_ENABLE_SECOND_DERIVATIVE
  /**
   * Fills the vector \p weights with the weight coefficients
   * of the Loop subdivision mask for evaluating the limit surface
   * at a node explicitly. The size of \p weights will be
   * 1 + \p valence, where \p valence is the number of neighbor
   * nodes of the node where the limit surface is to be
   * evaluated. The weight for the node itself is the first
   * element of \p weights.
   */
  static void loop_subdivision_mask(std::vector<Real> & weights,
                                    const unsigned int valence);


  /**
   * Builds the subdivision matrix \p A for the Loop scheme. The
   * size depends on the element's \p valence.
   */
  static void init_subdivision_matrix(DenseMatrix<Real> & A,
                                      unsigned int valence);
};



/**
 * Hierarchic finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 */
template <unsigned int Dim>
class FEHierarchic : public FE<Dim,HIERARCHIC>
{
public:

  /**
   * Constructor. Creates a hierarchic finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FEHierarchic(const FEType & fet) :
    FE<Dim,HIERARCHIC> (fet)
  {}
};



/**
 * Discontinuous Hierarchic finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Truman E. Ellis
 * \date 2011
 */
template <unsigned int Dim>
class FEL2Hierarchic : public FE<Dim,L2_HIERARCHIC>
{
public:

  /**
   * Constructor. Creates a hierarchic finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FEL2Hierarchic(const FEType & fet) :
    FE<Dim,L2_HIERARCHIC> (fet)
  {}
};



/**
 * Lagrange finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 */
template <unsigned int Dim>
class FELagrange : public FE<Dim,LAGRANGE>
{
public:

  /**
   * Constructor. Creates a Lagrange finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FELagrange(const FEType & fet) :
    FE<Dim,LAGRANGE> (fet)
  {}
};


/**
 * Discontinuous Lagrange finite elements.
 */
template <unsigned int Dim>
class FEL2Lagrange : public FE<Dim,L2_LAGRANGE>
{
public:

  /**
   * Constructor. Creates a discontinuous Lagrange finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FEL2Lagrange(const FEType & fet) :
    FE<Dim,L2_LAGRANGE> (fet)
  {}
};


/**
 * Monomial finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 */
template <unsigned int Dim>
class FEMonomial : public FE<Dim,MONOMIAL>
{
public:

  /**
   * Constructor. Creates a monomial finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FEMonomial(const FEType & fet) :
    FE<Dim,MONOMIAL> (fet)
  {}
};


/**
 * The FEScalar class is used for working with SCALAR variables.
 */
template <unsigned int Dim>
class FEScalar : public FE<Dim,SCALAR>
{
public:

  /**
   * Constructor. Creates a SCALAR finite element
   * which simply represents one or more
   * extra DOFs coupled to all other DOFs in
   * the system.
   */
  explicit
  FEScalar(const FEType & fet) :
    FE<Dim,SCALAR> (fet)
  {}
};


/**
 * XYZ finite elements.  These require specialization
 * because the shape functions are defined in terms of
 * physical XYZ coordinates rather than local coordinates.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 */
template <unsigned int Dim>
class FEXYZ : public FE<Dim,XYZ>
{
public:

  /**
   * Constructor. Creates a monomial finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FEXYZ(const FEType & fet) :
    FE<Dim,XYZ> (fet)
  {}

  /**
   * Explicitly call base class method.  This prevents some
   * compilers being confused by partially overriding this virtual function.
   * \note: pts need to be in reference space coordinates, not physical ones.
   */
  virtual void reinit (const Elem * elem,
                       const std::vector<Point> * const pts = nullptr,
                       const std::vector<Real> * const weights = nullptr) override
  { FE<Dim,XYZ>::reinit (elem, pts, weights); }

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of \p face.
   */
  virtual void reinit (const Elem * elem,
                       const unsigned int side,
                       const Real tolerance = TOLERANCE,
                       const std::vector<Point> * const pts = nullptr,
                       const std::vector<Real> * const weights = nullptr) override;


protected:

  /**
   * Update the various member data fields \p phi,
   * \p dphidxi, \p dphideta, \p dphidzeta, etc.
   * for the current element.  These data will be computed
   * at the points \p qp, which are generally (but need not be)
   * the quadrature points.
   */
  virtual void init_shape_functions(const std::vector<Point> & qp,
                                    const Elem * e) override;

  /**
   * After having updated the jacobian and the transformation
   * from local to global coordinates in \p FEAbstract::compute_map(),
   * the first derivatives of the shape functions are
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx, \p dphidy, and \p dphidz. This method
   * should rarely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected.
   */
  virtual void compute_shape_functions(const Elem * elem, const std::vector<Point> & qp) override;

  /**
   * Compute the map & shape functions for this face.
   */
  void compute_face_values (const Elem * elem,
                            const Elem * side,
                            const std::vector<Real> & weights);
};



/**
 * FELagrangeVec objects are used for working with vector-valued
 * finite elements
 *
 * \author Paul T. Bauman
 * \date 2013
 */
template <unsigned int Dim>
class FELagrangeVec : public FE<Dim,LAGRANGE_VEC>
{
public:

  /**
   * Constructor. Creates a vector Lagrange finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FELagrangeVec(const FEType & fet) :
    FE<Dim,LAGRANGE_VEC> (fet)
  {}
};



/**
 * FENedelecOne objects are used for working with vector-valued
 * Nedelec finite elements of the first kind.
 *
 * \author Paul T. Bauman
 * \date 2013
 */
template <unsigned int Dim>
class FENedelecOne : public FE<Dim,NEDELEC_ONE>
{
public:
  /**
   * Constructor. Creates a vector Lagrange finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FENedelecOne(const FEType & fet) :
    FE<Dim,NEDELEC_ONE> (fet)
  {}
};

/**
 * FEMonomialVec objects are used for working with vector-valued
 * discontinuous finite elements
 *
 * \author Alex D. Lindsay
 * \date 2019
 */
template <unsigned int Dim>
class FEMonomialVec : public FE<Dim,MONOMIAL_VEC>
{
public:

  /**
   * Constructor. Creates a vector Monomial finite element
   * to be used in dimension \p Dim.
   */
  explicit
  FEMonomialVec(const FEType & fet) :
    FE<Dim,MONOMIAL_VEC> (fet)
  {}
};



/**
 * Provide Typedefs for various element types.
 */
namespace FiniteElements
{
/**
 * Convenient definition for a 2D
 * Clough-Tocher finite element.
 */
typedef FEClough<2> FEClough2D;

/**
 * Convenient definition for a 1D
 * Hierarchic finite element.
 */
typedef FE<1,HIERARCHIC> FEHierarchic1D;

/**
 * Convenient definition for a 2D
 * Hierarchic finite element.
 */
typedef FE<2,HIERARCHIC> FEHierarchic2D;

/**
 * Convenient definition for a 3D
 * Hierarchic finite element.
 */
typedef FE<3,HIERARCHIC> FEHierarchic3D;


/**
 * Convenient definition for a 1D
 * Discontinuous Hierarchic finite element.
 */
typedef FE<1,L2_HIERARCHIC> FEL2Hierarchic1D;

/**
 * Convenient definition for a 2D
 * Discontinuous Hierarchic finite element.
 */
typedef FE<2,L2_HIERARCHIC> FEL2Hierarchic2D;

/**
 * Convenient definition for a 3D
 * Discontinuous Hierarchic finite element.
 */
typedef FE<3,L2_HIERARCHIC> FEL2Hierarchic3D;


/**
 * Convenient definition for a 1D
 * Lagrange finite element.
 */
typedef FE<1,LAGRANGE> FELagrange1D;

/**
 * Convenient definition for a 2D
 * Lagrange finite element.
 */
typedef FE<2,LAGRANGE> FELagrange2D;

/**
 * Convenient definition for a 3D
 * Lagrange finite element.
 */
typedef FE<3,LAGRANGE> FELagrange3D;


/**
 * Convenient definition for a 1D
 * Discontinuous Lagrange finite element.
 */
typedef FE<1,L2_LAGRANGE> FEL2Lagrange1D;

/**
 * Convenient definition for a 2D
 * Discontinuous Lagrange finite element.
 */
typedef FE<2,L2_LAGRANGE> FEL2Lagrange2D;

/**
 * Convenient definition for a 3D
 * Discontinuous Lagrange finite element.
 */
typedef FE<3,L2_LAGRANGE> FEL2Lagrange3D;


/**
 * Convenient definition for a 1D
 * Monomial finite element.
 */
typedef FE<1,MONOMIAL> FEMonomial1D;

/**
 * Convenient definition for a 2D
 * Monomial finite element.
 */
typedef FE<2,MONOMIAL> FEMonomial2D;

/**
 * Convenient definition for a 3D
 * Monomial finite element.
 */
typedef FE<3,MONOMIAL> FEMonomial3D;

}

} // namespace libMesh

#define LIBMESH_DEFAULT_VECTORIZED_FE(MyDim, MyType) \
template<>                                           \
void FE<MyDim,MyType>::all_shapes                    \
  (const Elem * elem,                                \
   const Order o,                                    \
   const std::vector<Point> & p,                     \
   std::vector<std::vector<OutputShape>> & v,        \
   const bool add_p_level)                           \
{                                                    \
  FE<MyDim,MyType>::default_all_shapes               \
    (elem,o,p,v,add_p_level);                        \
}                                                    \
                                                     \
template<>                                           \
void FE<MyDim,MyType>::shapes                        \
  (const Elem * elem,                                \
   const Order o,                                    \
   const unsigned int i,                             \
   const std::vector<Point> & p,                     \
   std::vector<OutputShape> & v,                     \
   const bool add_p_level)                           \
{                                                    \
  FE<MyDim,MyType>::default_shapes                   \
    (elem,o,i,p,v,add_p_level);                      \
}                                                    \
                                                     \
template<>                                           \
void FE<MyDim,MyType>::shape_derivs                  \
  (const Elem * elem,                                \
   const Order o,                                    \
   const unsigned int i,                             \
   const unsigned int j,                             \
   const std::vector<Point> & p,                     \
   std::vector<OutputShape> & v,                     \
   const bool add_p_level)                           \
{                                                    \
  FE<MyDim,MyType>::default_shape_derivs             \
    (elem,o,i,j,p,v,add_p_level);                    \
}


#endif // LIBMESH_FE_H
