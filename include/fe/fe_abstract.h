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



#ifndef LIBMESH_FE_ABSTRACT_H
#define LIBMESH_FE_ABSTRACT_H

// Local includes
#include "libmesh/reference_counted_object.h"
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/fe_type.h"
#include "libmesh/fe_map.h"
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#include "libmesh/tensor_value.h"
#endif

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum ElemType : int;
}
#else
#include "libmesh/enum_elem_type.h"
#endif

// C++ includes
#include <cstddef>
#include <vector>
#include <memory>

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#define virtual_for_inffe virtual
#else
#define virtual_for_inffe
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
 * This class forms the foundation from which generic finite elements
 * may be derived.  In the current implementation, the templated
 * derived class \p FE offers a wide variety of commonly used finite
 * element concepts.  Check there for details.  Use the \p
 * FEAbstract::build() method to create an object of any of the
 * derived classes.
 *
 * \note In the present design, the number of virtual members is kept
 * to a minimum for performance reasons, although this is not based on
 * rigorous profiling.
 *
 * All calls to static members of the \p FE classes should be
 * requested through the \p FEInterface.  This interface class
 * approximates runtime polymorphism for the templated finite element
 * classes.  Even internal library classes, like \p DofMap, request
 * the number of DOFs through this interface class.  This approach
 * also enables the co-existence of various element-based schemes.
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
   * Builds a specific finite element type.
   *
   * \returns A std::unique_ptr<FEAbstract> to the FE object to prevent
   * memory leaks.
   */
  static std::unique_ptr<FEAbstract> build (const unsigned int dim,
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
   * \note The FE classes decide which data to initialize based on
   * which accessor functions such as \p get_phi() or \p get_d2phi()
   * have been called, so all such accessors should be called before
   * the first \p reinit().
   */
  virtual void reinit (const Elem * elem,
                       const std::vector<Point> * const pts = nullptr,
                       const std::vector<Real> * const weights = nullptr) = 0;

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of the element \p elem.  The \p tolerance parameter
   * is passed to the involved call to \p inverse_map().  By default the
   * element data are computed at the quadrature points specified by the
   * quadrature rule \p qrule, but any set of points on the reference
   * \em side element may be specified in the optional argument \p pts.
   */
  virtual void reinit (const Elem * elem,
                       const unsigned int side,
                       const Real tolerance = TOLERANCE,
                       const std::vector<Point> * const pts = nullptr,
                       const std::vector<Real> * const weights = nullptr) = 0;

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p edge of the element \p elem.  The \p tolerance parameter
   * is passed to the involved call to \p inverse_map().  By default the
   * element data are computed at the quadrature points specified by the
   * quadrature rule \p qrule, but any set of points on the reference
   * \em edge element may be specified in the optional argument \p pts.
   */
  virtual void edge_reinit (const Elem * elem,
                            const unsigned int edge,
                            const Real tolerance = TOLERANCE,
                            const std::vector<Point> * pts = nullptr,
                            const std::vector<Real> * weights = nullptr) = 0;

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
   * \returns \p true if the point p is located on the reference element
   * for element type t, false otherwise.  Since we are doing floating
   * point comparisons here the parameter \p eps can be specified to
   * indicate a tolerance.  For example, \f$ x \le 1 \f$  becomes
   * \f$ x \le 1 + \epsilon \f$.
   */
  static bool on_reference_element(const Point & p,
                                   const ElemType t,
                                   const Real eps = TOLERANCE);
  /**
   * \returns The reference space coordinates of \p nodes based on the
   * element type.
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
   * \returns the dimension of this FE
   */
  unsigned int get_dim() const
  { return dim; }

  /**
   * \returns nothing, but lets the FE know you're explicitly
   * prerequesting calculations.  This is useful when you only want
   * the FE for n_quadrature_points, n_dofs_on_side, or other methods
   * that don't require shape function calculations, but you don't
   * want libMesh "backwards compatibility" mode to assume you've made
   * no prerequests and need to calculate everything.
   */
  void get_nothing() const
  { calculate_nothing = true; }

  /**
   * \returns The \p xyz spatial locations of the quadrature
   * points on the element.
   *
   * It is overwritten by infinite elements since there
   * \p FEMap cannot be used to compute \p xyz.
   */
  virtual_for_inffe
  const std::vector<Point> & get_xyz() const
  { calculate_map = true; return this->_fe_map->get_xyz(); }


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  /**
   * This function is the variant of \p get_JxW() for \p InfFE.
   * Since J diverges there, a respectize decay-function must be
   * applied to obtain well-defined quantities.
   *
   * For FE, it is equivalent to the common \p get_JxW().
   */
  virtual const std::vector<Real> & get_JxWxdecay_sq () const
  { return get_JxW();}
#endif

  /**
   * \returns The element Jacobian times the quadrature weight for
   * each quadrature point.
   *
   * For \p InfFE, use \p get_JxWxdecay_sq() instead.
   */
  virtual_for_inffe
  const std::vector<Real> & get_JxW() const
  { calculate_map = true; return this->_fe_map->get_JxW(); }

  /**
   * \returns The element tangents in xi-direction at the quadrature
   * points.
   */
  virtual_for_inffe
  const std::vector<RealGradient> & get_dxyzdxi() const
  { calculate_map = true; return this->_fe_map->get_dxyzdxi(); }

  /**
   * \returns The element tangents in eta-direction at the quadrature
   * points.
   */
  virtual_for_inffe
  const std::vector<RealGradient> & get_dxyzdeta() const
  { calculate_map = true; return this->_fe_map->get_dxyzdeta(); }

  /**
   * \returns The element tangents in zeta-direction at the quadrature
   * points.
   */
  virtual_for_inffe
  const std::vector<RealGradient> & get_dxyzdzeta() const
  { return _fe_map->get_dxyzdzeta(); }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The second partial derivatives in xi.
   */
  virtual_for_inffe
  const std::vector<RealGradient> & get_d2xyzdxi2() const
  { calculate_map = true; return this->_fe_map->get_d2xyzdxi2(); }

  /**
   * \returns The second partial derivatives in eta.
   */
  virtual_for_inffe
  const std::vector<RealGradient> & get_d2xyzdeta2() const
  { calculate_map = true; return this->_fe_map->get_d2xyzdeta2(); }

  /**
   * \returns The second partial derivatives in zeta.
   */
  virtual_for_inffe
  const std::vector<RealGradient> & get_d2xyzdzeta2() const
  { calculate_map = true; return this->_fe_map->get_d2xyzdzeta2(); }

  /**
   * \returns The second partial derivatives in xi-eta.
   */
  virtual_for_inffe
  const std::vector<RealGradient> & get_d2xyzdxideta() const
  { calculate_map = true; return this->_fe_map->get_d2xyzdxideta(); }

  /**
   * \returns The second partial derivatives in xi-zeta.
   */
  virtual_for_inffe
  const std::vector<RealGradient> & get_d2xyzdxidzeta() const
  { calculate_map = true; return this->_fe_map->get_d2xyzdxidzeta(); }

  /**
   * \returns The second partial derivatives in eta-zeta.
   */
  virtual_for_inffe
  const std::vector<RealGradient> & get_d2xyzdetadzeta() const
  { calculate_map = true; return this->_fe_map->get_d2xyzdetadzeta(); }

#endif

  /**
   * \returns The dxi/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  virtual_for_inffe
  const std::vector<Real> & get_dxidx() const
  { calculate_map = true; return this->_fe_map->get_dxidx(); }

  /**
   * \returns The dxi/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  virtual_for_inffe
  const std::vector<Real> & get_dxidy() const
  { calculate_map = true; return this->_fe_map->get_dxidy(); }

  /**
   * \returns The dxi/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  virtual_for_inffe
  const std::vector<Real> & get_dxidz() const
  { calculate_map = true; return this->_fe_map->get_dxidz(); }

  /**
   * \returns The deta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  virtual_for_inffe
  const std::vector<Real> & get_detadx() const
  { calculate_map = true; return this->_fe_map->get_detadx(); }

  /**
   * \returns The deta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  virtual_for_inffe
  const std::vector<Real> & get_detady() const
  { calculate_map = true; return this->_fe_map->get_detady(); }

  /**
   * \returns The deta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  virtual_for_inffe
  const std::vector<Real> & get_detadz() const
  { calculate_map = true; return this->_fe_map->get_detadz(); }

  /**
   * \returns The dzeta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  virtual_for_inffe
  const std::vector<Real> & get_dzetadx() const
  { calculate_map = true; return this->_fe_map->get_dzetadx(); }

  /**
   * \returns The dzeta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  virtual_for_inffe
  const std::vector<Real> & get_dzetady() const
  { calculate_map = true; return this->_fe_map->get_dzetady(); }

  /**
   * \returns The dzeta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  virtual_for_inffe
  const std::vector<Real> & get_dzetadz() const
  { calculate_map = true; return this->_fe_map->get_dzetadz(); }

  /**
   * \returns The tangent vectors for face integration.
   */
  virtual_for_inffe
  const std::vector<std::vector<Point>> & get_tangents() const
  { calculate_map = true; return this->_fe_map->get_tangents(); }

  /**
   * \returns The outward pointing normal vectors for face integration.
   */
  virtual_for_inffe
  const std::vector<Point> & get_normals() const
  { calculate_map = true; return this->_fe_map->get_normals(); }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * \returns The curvatures for use in face integration.
   */
  virtual_for_inffe
  const std::vector<Real> & get_curvatures() const
  { calculate_map = true; return this->_fe_map->get_curvatures();}

#endif

  /**
   * Provides the class with the quadrature rule.  Implement
   * this in derived classes.
   */
  virtual void attach_quadrature_rule (QBase * q) = 0;

  /**
   * \returns The total number of approximation shape functions
   * for the current element.  Useful during matrix assembly.
   * Implement this in derived classes.
   */
  virtual unsigned int n_shape_functions () const = 0;

  /**
   * \returns The total number of quadrature points.  Useful
   * during matrix assembly.  Implement this in derived classes.
   */
  virtual unsigned int n_quadrature_points () const = 0;

  /**
   * \returns The element type that the current shape functions
   * have been calculated for.  Useful in determining when shape
   * functions must be recomputed.
   */
  ElemType get_type()  const { return elem_type; }

  /**
   * \returns The p refinement level that the current shape
   * functions have been calculated for.
   */
  unsigned int get_p_level() const { return _p_level; }

  /**
   * \returns The FE Type (approximation order and family) of the finite element.
   */
  FEType get_fe_type()  const { return fe_type; }

  /**
   * \returns The approximation order of the finite element.
   */
  Order get_order()  const { return static_cast<Order>(fe_type.order + _p_level); }

  /**
   * Sets the *base* FE order of the finite element.
   */
  void set_fe_order(int new_order) { fe_type.order = new_order; }

  /**
   * \returns The continuity level of the finite element.
   */
  virtual FEContinuity get_continuity() const = 0;

  /**
   * \returns \p true if the finite element's higher order shape functions are
   * hierarchic
   */
  virtual bool is_hierarchic() const = 0;

  /**
   * \returns The finite element family of this element.
   */
  FEFamily get_family()  const { return fe_type.family; }

  /**
   * \returns The mapping object
   *
   * \note for InfFE, this gives a useless object.
   */
  const FEMap & get_fe_map() const { return *_fe_map.get(); }
  FEMap & get_fe_map() { return *_fe_map.get(); }

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
  virtual void print_dual_phi(std::ostream & os) const =0;

  /**
   * Prints the value of each shape function's derivative
   * at each quadrature point. Implement in derived class since this
   * depends on whether the element is vector-valued or not.
   */
  virtual void print_dphi(std::ostream & os) const =0;
  virtual void print_dual_dphi(std::ostream & os) const =0;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Prints the value of each shape function's second derivatives
   * at each quadrature point. Implement in derived class since this
   * depends on whether the element is vector-valued or not.
   */
  virtual void print_d2phi(std::ostream & os) const =0;
  virtual void print_dual_d2phi(std::ostream & os) const =0;

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

  /**
   * request phi calculations
   */
  virtual void request_phi() const = 0;
  virtual void request_dual_phi() const = 0;

  /**
   * request dphi calculations
   */
  virtual void request_dphi() const = 0;
  virtual void request_dual_dphi() const = 0;

  /**
   * set calculate_dual as needed
   */
  void set_calculate_dual(const bool val){calculate_dual = val; }

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

  std::unique_ptr<FEMap> _fe_map;


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
   * Are we calculating dual basis?
   */
  mutable bool calculate_dual;

  /**
   * Are we potentially deliberately calculating nothing?
   */
  mutable bool calculate_nothing;

  /**
   * Are we calculating mapping functions?
   */
  mutable bool calculate_map;

  /**
   * Should we calculate shape functions?
   */
  mutable bool calculate_phi;

  /**
   * Should we calculate shape function gradients?
   */
  mutable bool calculate_dphi;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Should we calculate shape function hessians?
   */
  mutable bool calculate_d2phi;
#else
  // Otherwise several interfaces need to be redone.
  const bool calculate_d2phi=false;

#endif

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
   * The finite element type for this object.
   *
   * \note This should be constant for the object.
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
   * \returns \p true when the shape functions (for
   * this \p FEFamily) depend on the particular
   * element, and therefore needs to be re-initialized
   * for each new element.  \p false otherwise.
   */
  virtual bool shapes_need_reinit() const = 0;

};

} // namespace libMesh

#endif // LIBMESH_FE_ABSTRACT_H
