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



#ifndef LIBMESH_FE_ABSTRACT_H
#define LIBMESH_FE_ABSTRACT_H

// Local includes
#include "libmesh/reference_counted_object.h"
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/fe_type.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/fe_map.h"

// C++ includes
#include <cstddef>
#include <vector>

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#include "libmesh/tensor_value.h"
#endif

namespace libMesh
{


// forward declarations
template <typename T> class DenseMatrix;
template <typename T> class DenseVector;
class BoundaryInfo;
class DofConstraints;
class DofMap;
class Elem;
class MeshBase;
template <typename T> class NumericVector;
class QBase;

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
class NodeConstraints;
#endif

#ifdef LIBMESH_ENABLE_PERIODIC
class PeriodicBoundaries;
class PointLocatorBase;
#endif

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
class InfFE;
#endif



/**
 * This class forms the foundation from which generic finite
 * elements may be derived.  In the current implementation the
 * templated derived class \p FE offers a wide variety of commonly
 * used finite element concepts.  Check there for details.
 * Use the \p FEAbstract::build() method to create an object of any of
 * the derived classes.
 * Note that the amount of virtual members is kept to a minimum,
 * and the sophisticated template scheme of \p FE is quite
 * likely to offer acceptably fast code.
 *
 * All calls to static members of the \p FE classes should be
 * requested through the \p FEInterface.  This interface class
 * offers sort-of runtime polymorphism for the templated finite
 * element classes.  Even internal library classes, like \p DofMap,
 * request the number of dof's through this interface class.
 * Note that this also enables the co-existence of various
 * element-based schemes.
 * This class is well 'at the heart' of the library, so
 * things in here should better remain unchanged.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 */
class FEAbstract : public ReferenceCountedObject<FEAbstract>
{
protected:

  /**
   * Constructor.  Optionally initializes required data
   * structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  FEAbstract (const unsigned int dim,
              const FEType & fet);

public:

  /**
   * Destructor.
   */
  virtual ~FEAbstract();

  /**
   * Builds a specific finite element type.  A \p
   * UniquePtr<FEAbstract> is returned to prevent a memory leak. This
   * way the user need not remember to delete the object.
   */
  static UniquePtr<FEAbstract> build (const unsigned int dim,
                                      const FEType & type);

  /**
   * This is at the core of this class. Use this for each
   * new element in the mesh.  Reinitializes the requested physical
   * element-dependent data based on the current element
   * \p elem. By default the element data are computed at the quadrature
   * points specified by the quadrature rule \p qrule, but any set
   * of points on the reference element may be specified in the optional
   * argument \p pts.
   *
   * Note that the FE classes decide which data to initialize based on
   * which accessor functions such as \p get_phi() or \p get_d2phi() have
   * been called, so all such accessors should be called before the first
   * \p reinit().
   */
  virtual void reinit (const Elem * elem,
                       const std::vector<Point> * const pts = libmesh_nullptr,
                       const std::vector<Real> * const weights = libmesh_nullptr) = 0;

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of the element \p elem.  The \p tolerance paremeter
   * is passed to the involved call to \p inverse_map().  By default the
   * element data are computed at the quadrature points specified by the
   * quadrature rule \p qrule, but any set of points on the reference
   * \em side element may be specified in the optional argument \p pts.
   */
  virtual void reinit (const Elem * elem,
                       const unsigned int side,
                       const Real tolerance = TOLERANCE,
                       const std::vector<Point> * const pts = libmesh_nullptr,
                       const std::vector<Real> * const weights = libmesh_nullptr) = 0;

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p edge of the element \p elem.  The \p tolerance paremeter
   * is passed to the involved call to \p inverse_map().  By default the
   * element data are computed at the quadrature points specified by the
   * quadrature rule \p qrule, but any set of points on the reference
   * \em edge element may be specified in the optional argument \p pts.
   */
  virtual void edge_reinit (const Elem * elem,
                            const unsigned int edge,
                            const Real tolerance = TOLERANCE,
                            const std::vector<Point> * pts = libmesh_nullptr,
                            const std::vector<Real> * weights = libmesh_nullptr) = 0;

  /**
   * Computes the reference space quadrature points on the side of
   * an element based on the side quadrature points.
   */
  virtual void side_map (const Elem * elem,
                         const Elem * side,
                         const unsigned int s,
                         const std::vector<Point> & reference_side_points,
                         std::vector<Point> &       reference_points) = 0;

  /**
   * @returns true if the point p is located on the reference element
   * for element type t, false otherwise.  Since we are doing floating
   * point comparisons here the parameter \p eps can be specified to
   * indicate a tolerance.  For example, \f$ x \le 1 \f$  becomes
   * \f$ x \le 1 + \epsilon \f$.
   */
  static bool on_reference_element(const Point & p,
                                   const ElemType t,
                                   const Real eps = TOLERANCE);
  /**
   * returns the reference space nodes coordinates
   * given the element type
   */
  static void get_refspace_nodes(const ElemType t,
                                 std::vector<Point> & nodes);


#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * Computes the nodal constraint contributions (for
   * non-conforming adapted meshes), using Lagrange geometry
   */
  static void compute_node_constraints (NodeConstraints & constraints,
                                        const Elem * elem);
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#ifdef LIBMESH_ENABLE_PERIODIC

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * Computes the node position constraint equation contributions (for
   * meshes with periodic boundary conditions)
   */
  static void compute_periodic_node_constraints (NodeConstraints & constraints,
                                                 const PeriodicBoundaries & boundaries,
                                                 const MeshBase & mesh,
                                                 const PointLocatorBase * point_locator,
                                                 const Elem * elem);
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#endif // LIBMESH_ENABLE_PERIODIC

  /**
   * @returns the \p xyz spatial locations of the quadrature
   * points on the element.
   */
  const std::vector<Point> & get_xyz() const
  { return this->_fe_map->get_xyz(); }

  /**
   * @returns the element Jacobian times the quadrature weight for
   * each quadrature point.
   */
  const std::vector<Real> & get_JxW() const
  { return this->_fe_map->get_JxW(); }

  /**
   * @returns the element tangents in xi-direction at the quadrature
   * points.
   */
  const std::vector<RealGradient> & get_dxyzdxi() const
  { return this->_fe_map->get_dxyzdxi(); }

  /**
   * @returns the element tangents in eta-direction at the quadrature
   * points.
   */
  const std::vector<RealGradient> & get_dxyzdeta() const
  { return this->_fe_map->get_dxyzdeta(); }

  /**
   * @returns the element tangents in zeta-direction at the quadrature
   * points.
   */
  const std::vector<RealGradient> & get_dxyzdzeta() const
  { return _fe_map->get_dxyzdzeta(); }

  /**
   * @returns the second partial derivatives in xi.
   */
  const std::vector<RealGradient> & get_d2xyzdxi2() const
  { return this->_fe_map->get_d2xyzdxi2(); }

  /**
   * @returns the second partial derivatives in eta.
   */
  const std::vector<RealGradient> & get_d2xyzdeta2() const
  { return this->_fe_map->get_d2xyzdeta2(); }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * @returns the second partial derivatives in zeta.
   */
  const std::vector<RealGradient> & get_d2xyzdzeta2() const
  { return this->_fe_map->get_d2xyzdzeta2(); }

#endif

  /**
   * @returns the second partial derivatives in xi-eta.
   */
  const std::vector<RealGradient> & get_d2xyzdxideta() const
  { return this->_fe_map->get_d2xyzdxideta(); }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * @returns the second partial derivatives in xi-zeta.
   */
  const std::vector<RealGradient> & get_d2xyzdxidzeta() const
  { return this->_fe_map->get_d2xyzdxidzeta(); }

  /**
   * @returns the second partial derivatives in eta-zeta.
   */
  const std::vector<RealGradient> & get_d2xyzdetadzeta() const
  { return this->_fe_map->get_d2xyzdetadzeta(); }

#endif

  /**
   * @returns the dxi/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dxidx() const
  { return this->_fe_map->get_dxidx(); }

  /**
   * @returns the dxi/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dxidy() const
  { return this->_fe_map->get_dxidy(); }

  /**
   * @returns the dxi/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dxidz() const
  { return this->_fe_map->get_dxidz(); }

  /**
   * @returns the deta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_detadx() const
  { return this->_fe_map->get_detadx(); }

  /**
   * @returns the deta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_detady() const
  { return this->_fe_map->get_detady(); }

  /**
   * @returns the deta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_detadz() const
  { return this->_fe_map->get_detadz(); }

  /**
   * @returns the dzeta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dzetadx() const
  { return this->_fe_map->get_dzetadx(); }

  /**
   * @returns the dzeta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dzetady() const
  { return this->_fe_map->get_dzetady(); }

  /**
   * @returns the dzeta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dzetadz() const
  { return this->_fe_map->get_dzetadz(); }

  /**
   * @returns the tangent vectors for face integration.
   */
  const std::vector<std::vector<Point> > & get_tangents() const
  { return this->_fe_map->get_tangents(); }

  /**
   * @returns the outward pointing normal vectors for face integration.
   */
  const std::vector<Point> & get_normals() const
  { return this->_fe_map->get_normals(); }

  /**
   * @returns the curvatures for use in face integration.
   */
  const std::vector<Real> & get_curvatures() const
  { return this->_fe_map->get_curvatures();}

  /**
   * Provides the class with the quadrature rule.  Implement
   * this in derived classes.
   */
  virtual void attach_quadrature_rule (QBase * q) = 0;

  /**
   * @returns the total number of approximation shape functions
   * for the current element.  Useful during matrix assembly.
   * Implement this in derived classes.
   */
  virtual unsigned int n_shape_functions () const = 0;

  /**
   * @returns the total number of quadrature points.  Useful
   * during matrix assembly.  Implement this in derived classes.
   */
  virtual unsigned int n_quadrature_points () const = 0;

  /**
   * @returns the element type that the current shape functions
   * have been calculated for.  Useful in determining when shape
   * functions must be recomputed.
   */
  ElemType get_type()  const { return elem_type; }

  /**
   * @returns the p refinement level that the current shape
   * functions have been calculated for.
   */
  unsigned int get_p_level() const { return _p_level; }

  /**
   * @returns the FE Type (approximation order and family) of the finite element.
   */
  FEType get_fe_type()  const { return fe_type; }

  /**
   * @returns the approximation order of the finite element.
   */
  Order get_order()  const { return static_cast<Order>(fe_type.order + _p_level); }

  /**
   * Sets the *base* FE order of the finite element.
   */
  void set_fe_order(int new_order) { fe_type.order = new_order; }

  /**
   * @returns the continuity level of the finite element.
   */
  virtual FEContinuity get_continuity() const = 0;

  /**
   * @returns true if the finite element's higher order shape functions are
   * hierarchic
   */
  virtual bool is_hierarchic() const = 0;

  /**
   * @returns the finite element family of this element.
   */
  FEFamily get_family()  const { return fe_type.family; }

  /**
   * @returns the mapping object
   */
  const FEMap & get_fe_map() const { return *_fe_map.get(); }

  /**
   * Prints the Jacobian times the weight for each quadrature point.
   */
  void print_JxW(std::ostream & os) const;

  /**
   * Prints the value of each shape function at each quadrature point.
   * Implement in derived class since this depends on whether the element
   * is vector-valued or not.
   */
  virtual void print_phi(std::ostream & os) const =0;

  /**
   * Prints the value of each shape function's derivative
   * at each quadrature point. Implement in derived class since this
   * depends on whether the element is vector-valued or not.
   */
  virtual void print_dphi(std::ostream & os) const =0;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Prints the value of each shape function's second derivatives
   * at each quadrature point. Implement in derived class since this
   * depends on whether the element is vector-valued or not.
   */
  virtual void print_d2phi(std::ostream & os) const =0;

#endif

  /**
   * Prints the spatial location of each quadrature point
   * (on the physical element).
   */
  void print_xyz(std::ostream & os) const;

  /**
   * Prints all the relevant information about the current element.
   */
  void print_info(std::ostream & os) const;

  /**
   * Same as above, but allows you to print to a stream.
   */
  friend std::ostream & operator << (std::ostream & os, const FEAbstract & fe);


protected:

  /**
   * After having updated the jacobian and the transformation
   * from local to global coordinates in \p FEMap::compute_map(),
   * the first derivatives of the shape functions are
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx, \p dphidy, and \p dphidz. This method
   * should rarely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected. This needs to be implemented in the
   * derived class since this function depends on whether
   * the shape functions are vector-valued or not.
   */
  virtual void compute_shape_functions(const Elem *, const std::vector<Point> & ) =0;

  UniquePtr<FEMap> _fe_map;


  /**
   * The dimensionality of the object
   */
  const unsigned int dim;

  /**
   * Have calculations with this object already been started?
   * Then all get_* functions should already have been called.
   */
  mutable bool calculations_started;

  /**
   * Should we calculate shape functions?
   */
  mutable bool calculate_phi;

  /**
   * Should we calculate shape function gradients?
   */
  mutable bool calculate_dphi;

  /**
   * Should we calculate shape function hessians?
   */
  mutable bool calculate_d2phi;

  /**
   * Should we calculate shape function curls?
   */
  mutable bool calculate_curl_phi;

  /**
   * Should we calculate shape function divergences?
   */
  mutable bool calculate_div_phi;

  /**
   * Should we calculate reference shape function gradients?
   */
  mutable bool calculate_dphiref;


  /**
   * The finite element type for this object.  Note that this
   * should be constant for the object.
   */
  FEType fe_type;

  /**
   * The element type the current data structures are
   * set up for.
   */
  ElemType elem_type;

  /**
   * The p refinement level the current data structures are
   * set up for.
   */
  unsigned int _p_level;

  /**
   * A pointer to the quadrature rule employed
   */
  QBase * qrule;

  /**
   * A flag indicating if current data structures
   * correspond to quadrature rule points
   */
  bool shapes_on_quadrature;

  /**
   * @returns \p true when the shape functions (for
   * this \p FEFamily) depend on the particular
   * element, and therefore needs to be re-initialized
   * for each new element.  \p false otherwise.
   */
  virtual bool shapes_need_reinit() const = 0;

};




// ------------------------------------------------------------
// FEAbstract class inline members
inline
FEAbstract::FEAbstract(const unsigned int d,
                       const FEType & fet) :
  _fe_map( FEMap::build(fet) ),
  dim(d),
  calculations_started(false),
  calculate_phi(false),
  calculate_dphi(false),
  calculate_d2phi(false),
  calculate_curl_phi(false),
  calculate_div_phi(false),
  calculate_dphiref(false),
  fe_type(fet),
  elem_type(INVALID_ELEM),
  _p_level(0),
  qrule(libmesh_nullptr),
  shapes_on_quadrature(false)
{
}


inline
FEAbstract::~FEAbstract()
{
}

} // namespace libMesh

#endif // LIBMESH_FE_ABSTRACT_H
