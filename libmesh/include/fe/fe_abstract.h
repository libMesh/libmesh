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



#ifndef __fe_abstract_h__
#define __fe_abstract_h__

// Local includes
#include "reference_counted_object.h"
#include "point.h"
#include "vector_value.h"
#include "enum_elem_type.h"
#include "fe_type.h"
#include "auto_ptr.h"

// C++ includes
#include <cstddef>
#include <vector>

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#include "tensor_value.h"
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
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// FEAbstract class definition
class FEAbstract : public ReferenceCountedObject<FEAbstract>
{
protected:

  /**
   * Constructor.  Optionally initializes required data
   * structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  FEAbstract (const unsigned int dim,
	      const FEType& fet);

public:

  /**
   * Destructor.
   */
  virtual ~FEAbstract();

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
  virtual void reinit (const Elem* elem,
		       const std::vector<Point>* const pts = NULL,
                       const std::vector<Real>* const weights = NULL) = 0;

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of the element \p elem.  The \p tolerance paremeter
   * is passed to the involved call to \p inverse_map().  By default the
   * element data are computed at the quadrature points specified by the
   * quadrature rule \p qrule, but any set of points on the reference
   * \em side element may be specified in the optional argument \p pts.
   */
  virtual void reinit (const Elem* elem,
		       const unsigned int side,
		       const Real tolerance = TOLERANCE,
                       const std::vector<Point>* const pts = NULL,
                       const std::vector<Real>* const weights = NULL) = 0;

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p edge of the element \p elem.  The \p tolerance paremeter
   * is passed to the involved call to \p inverse_map().  By default the
   * element data are computed at the quadrature points specified by the
   * quadrature rule \p qrule, but any set of points on the reference
   * \em edge element may be specified in the optional argument \p pts.
   */
  virtual void edge_reinit (const Elem* elem,
		            const unsigned int edge,
			    const Real tolerance = TOLERANCE,
                            const std::vector<Point>* pts = NULL,
                            const std::vector<Real>* weights = NULL) = 0;

  /**
   * Computes the reference space quadrature points on the side of
   * an element based on the side quadrature points.
   */
  virtual void side_map (const Elem* elem,
                         const Elem* side,
                         const unsigned int s,
                         const std::vector<Point>& reference_side_points,
                         std::vector<Point>&       reference_points) = 0;

  /**
   * @returns true if the point p is located on the reference element
   * for element type t, false otherwise.  Since we are doing floating
   * point comparisons here the parameter \p eps can be specified to
   * indicate a tolerance.  For example, \f$ x \le 1 \f$  becomes
   * \f$ x \le 1 + \epsilon \f$.
   */
  static bool on_reference_element(const Point& p,
				   const ElemType t,
				   const Real eps = TOLERANCE);
  /**
   * returns the reference space nodes coordinates
   * given the element type
   */
  static void get_refspace_nodes(const ElemType t,
	                         std::vector<Point>& nodes);


#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * Computes the nodal constraint contributions (for
   * non-conforming adapted meshes), using Lagrange geometry
   */
  static void compute_node_constraints (NodeConstraints &constraints,
                                        const Elem* elem);
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#ifdef LIBMESH_ENABLE_PERIODIC

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  /**
   * Computes the node position constraint equation contributions (for
   * meshes with periodic boundary conditions)
   */
  static void compute_periodic_node_constraints (NodeConstraints &constraints,
                                                 const PeriodicBoundaries &boundaries,
                                                 const MeshBase& mesh,
                                                 const PointLocatorBase* point_locator,
                                                 const Elem* elem);
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

#endif // LIBMESH_ENABLE_PERIODIC

  /**
   * @returns the \p xyz spatial locations of the quadrature
   * points on the element.
   */
  const std::vector<Point>& get_xyz() const
  { return xyz; }

  /**
   * @returns the element Jacobian times the quadrature weight for
   * each quadrature point.
   */
  const std::vector<Real>& get_JxW() const
  { return JxW; }

  /**
   * @returns the element tangents in xi-direction at the quadrature
   * points.
   */
  const std::vector<RealGradient>& get_dxyzdxi() const
  { return dxyzdxi_map; }

  /**
   * @returns the element tangents in eta-direction at the quadrature
   * points.
   */
  const std::vector<RealGradient>& get_dxyzdeta() const
  { return dxyzdeta_map; }

  /**
   * @returns the element tangents in zeta-direction at the quadrature
   * points.
   */
  const std::vector<RealGradient>& get_dxyzdzeta() const
  { return dxyzdzeta_map; }

  /**
   * @returns the second partial derivatives in xi.
   */
  const std::vector<RealGradient>& get_d2xyzdxi2() const
  { return d2xyzdxi2_map; }

  /**
   * @returns the second partial derivatives in eta.
   */
  const std::vector<RealGradient>& get_d2xyzdeta2() const
  { return d2xyzdeta2_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * @returns the second partial derivatives in zeta.
   */
  const std::vector<RealGradient>& get_d2xyzdzeta2() const
  { return d2xyzdzeta2_map; }

#endif

  /**
   * @returns the second partial derivatives in xi-eta.
   */
  const std::vector<RealGradient>& get_d2xyzdxideta() const
  { return d2xyzdxideta_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * @returns the second partial derivatives in xi-zeta.
   */
  const std::vector<RealGradient>& get_d2xyzdxidzeta() const
  { return d2xyzdxidzeta_map; }

  /**
   * @returns the second partial derivatives in eta-zeta.
   */
  const std::vector<RealGradient>& get_d2xyzdetadzeta() const
  { return d2xyzdetadzeta_map; }

#endif

  /**
   * @returns the dxi/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real>& get_dxidx() const
  { return dxidx_map; }

  /**
   * @returns the dxi/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real>& get_dxidy() const
  { return dxidy_map; }

  /**
   * @returns the dxi/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real>& get_dxidz() const
  { return dxidz_map; }

  /**
   * @returns the deta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real>& get_detadx() const
  { return detadx_map; }

  /**
   * @returns the deta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real>& get_detady() const
  { return detady_map; }

  /**
   * @returns the deta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real>& get_detadz() const
  { return detadz_map; }

  /**
   * @returns the dzeta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real>& get_dzetadx() const
  { return dzetadx_map; }

  /**
   * @returns the dzeta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real>& get_dzetady() const
  { return dzetady_map; }

  /**
   * @returns the dzeta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real>& get_dzetadz() const
  { return dzetadz_map; }

  /**
   * @returns the tangent vectors for face integration.
   */
  const std::vector<std::vector<Point> >& get_tangents() const
  { return tangents; }

  /**
   * @returns the normal vectors for face integration.
   */
  const std::vector<Point>& get_normals() const
  { return normals; }

  /**
   * @returns the curvatures for use in face integration.
   */
  const std::vector<Real>& get_curvatures() const
  { return curvatures;}

  /**
   * Provides the class with the quadrature rule.  Implement
   * this in derived classes.
   */
  virtual void attach_quadrature_rule (QBase* q) = 0;

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
   * Prints the Jacobian times the weight for each quadrature point.
   */
  void print_JxW(std::ostream& os) const;

  /**
   * Prints the value of each shape function at each quadrature point.
   * Implement in derived class since this depends on whether the element
   * is vector-valued or not.
   */
  virtual void print_phi(std::ostream& os) const =0;

  /**
   * Prints the value of each shape function's derivative
   * at each quadrature point. Implement in derived class since this 
   * depends on whether the element is vector-valued or not.
   */
  virtual void print_dphi(std::ostream& os) const =0;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Prints the value of each shape function's second derivatives
   * at each quadrature point. Implement in derived class since this 
   * depends on whether the element is vector-valued or not.
   */
  virtual void print_d2phi(std::ostream& os) const =0;

#endif

  /**
   * Prints the spatial location of each quadrature point
   * (on the physical element).
   */
  void print_xyz(std::ostream& os) const;

  /**
   * Prints all the relevant information about the current element.
   */
  void print_info(std::ostream& os) const;

  /**
   * Same as above, but allows you to print to a stream.
   */
  friend std::ostream& operator << (std::ostream& os, const FEAbstract& fe);


protected:


  /**
   * Compute the jacobian and some other additional
   * data fields. Takes the integration weights
   * as input, along with a pointer to the element.
   */
  virtual void compute_map(const std::vector<Real>& qw,
			   const Elem* e);

  /**
   * Compute the jacobian and some other additional
   * data fields. Takes the integration weights
   * as input, along with a pointer to the element.
   * The element is assumed to have a constant Jacobian
   */
  virtual void compute_affine_map(const std::vector<Real>& qw,
				  const Elem* e);

  /**
   * Compute the jacobian and some other additional
   * data fields at the single point with index p.
   */
  void compute_single_point_map(const std::vector<Real>& qw,
		                const Elem* e,
				unsigned int p);
  
  /**
   * A utility function for use by compute_*_map
   */
  void resize_map_vectors(unsigned int n_qp);

  /**
   * Same as compute_map, but for a side.  Useful for boundary integration.
   */
  void compute_face_map(const std::vector<Real>& qw,
			const Elem* side);

  /**
   * Same as before, but for an edge.  Useful for some projections.
   */
  void compute_edge_map(const std::vector<Real>& qw,
			const Elem* side);

  /**
   * After having updated the jacobian and the transformation
   * from local to global coordinates in \p FEAbstract::compute_map(),
   * the first derivatives of the shape functions are
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx, \p dphidy, and \p dphidz. This method
   * should rarely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected. This needs to be implemented in the
   * derived class since this function depends on whether
   * the shape functions are vector-valued or not.
   */
  virtual void compute_shape_functions(const Elem*) =0;

  /**
   * Used in \p FEAbstract::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the x value of the pth entry of the dxzydxi_map.
   */
  Real dxdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](0); }

  /**
   * Used in \p FEAbstract::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the y value of the pth entry of the dxzydxi_map.
   */
  Real dydxi_map(const unsigned int p) const   { return dxyzdxi_map[p](1); }

  /**
   * Used in \p FEAbstract::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the z value of the pth entry of the dxzydxi_map.
   */
  Real dzdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](2); }

  /**
   * Used in \p FEAbstract::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the x value of the pth entry of the dxzydeta_map.
   */
  Real dxdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](0); }

  /**
   * Used in \p FEAbstract::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the y value of the pth entry of the dxzydeta_map.
   */
  Real dydeta_map(const unsigned int p) const  { return dxyzdeta_map[p](1); }

  /**
   * Used in \p FEAbstract::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the z value of the pth entry of the dxzydeta_map.
   */
  Real dzdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](2); }

  /**
   * Used in \p FEAbstract::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the x value of the pth entry of the dxzydzeta_map.
   */
  Real dxdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](0); }

  /**
   * Used in \p FEAbstract::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the y value of the pth entry of the dxzydzeta_map.
   */
  Real dydzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](1); }

  /**
   * Used in \p FEAbstract::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the z value of the pth entry of the dxzydzeta_map.
   */
  Real dzdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](2); }





  /**
   * The dimensionality of the object
   */
  const unsigned int dim;

  /**
   * The spatial locations of the quadrature points
   */
  std::vector<Point> xyz;



  /**
   * Vector of parital derivatives:
   * d(x)/d(xi), d(y)/d(xi), d(z)/d(xi)
   */
  std::vector<RealGradient> dxyzdxi_map;

  /**
   * Vector of parital derivatives:
   * d(x)/d(eta), d(y)/d(eta), d(z)/d(eta)
   */
  std::vector<RealGradient> dxyzdeta_map;

  /**
   * Vector of parital derivatives:
   * d(x)/d(zeta), d(y)/d(zeta), d(z)/d(zeta)
   */
  std::vector<RealGradient> dxyzdzeta_map;

  /**
   * Vector of second partial derivatives in xi:
   * d^2(x)/d(xi)^2, d^2(y)/d(xi)^2, d^2(z)/d(xi)^2
   */
  std::vector<RealGradient> d2xyzdxi2_map;

  /**
   * Vector of mixed second partial derivatives in xi-eta:
   * d^2(x)/d(xi)d(eta) d^2(y)/d(xi)d(eta) d^2(z)/d(xi)d(eta)
   */
  std::vector<RealGradient> d2xyzdxideta_map;

  /**
   * Vector of second partial derivatives in eta:
   * d^2(x)/d(eta)^2
   */
  std::vector<RealGradient> d2xyzdeta2_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Vector of second partial derivatives in xi-zeta:
   * d^2(x)/d(xi)d(zeta), d^2(y)/d(xi)d(zeta), d^2(z)/d(xi)d(zeta)
   */
  std::vector<RealGradient> d2xyzdxidzeta_map;

  /**
   * Vector of mixed second partial derivatives in eta-zeta:
   * d^2(x)/d(eta)d(zeta) d^2(y)/d(eta)d(zeta) d^2(z)/d(eta)d(zeta)
   */
  std::vector<RealGradient> d2xyzdetadzeta_map;

  /**
   * Vector of second partial derivatives in zeta:
   * d^2(x)/d(zeta)^2
   */
  std::vector<RealGradient> d2xyzdzeta2_map;

#endif

  /**
   * Map for partial derivatives:
   * d(xi)/d(x). Needed for the Jacobian.
   */
  std::vector<Real>  dxidx_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(y). Needed for the Jacobian.
   */
  std::vector<Real>  dxidy_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(z). Needed for the Jacobian.
   */
  std::vector<Real>  dxidz_map;




  /**
   * Map for partial derivatives:
   * d(eta)/d(x). Needed for the Jacobian.
   */
  std::vector<Real>  detadx_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(y). Needed for the Jacobian.
   */
  std::vector<Real>  detady_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(z). Needed for the Jacobian.
   */
  std::vector<Real>  detadz_map;





  /**
   * Map for partial derivatives:
   * d(zeta)/d(x). Needed for the Jacobian.
   */
  std::vector<Real>  dzetadx_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(y). Needed for the Jacobian.
   */
  std::vector<Real>  dzetady_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(z). Needed for the Jacobian.
   */
  std::vector<Real>  dzetadz_map;

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
   * Shape function values.
   */
  std::vector<std::vector<Real> >   phi;

  /**
   * Shape function derivative values.
   */
  std::vector<std::vector<RealGradient> >  dphi;

  /**
   * Shape function derivatives in the xi direction.
   */
  std::vector<std::vector<Real> >   dphidxi;

  /**
   * Shape function derivatives in the eta direction.
   */
  std::vector<std::vector<Real> >   dphideta;

  /**
   * Shape function derivatives in the zeta direction.
   */
  std::vector<std::vector<Real> >   dphidzeta;

  /**
   * Shape function derivatives in the x direction.
   */
  std::vector<std::vector<Real> >   dphidx;

  /**
   * Shape function derivatives in the y direction.
   */
  std::vector<std::vector<Real> >   dphidy;

  /**
   * Shape function derivatives in the z direction.
   */
  std::vector<std::vector<Real> >   dphidz;


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Shape function second derivative values.
   */
  std::vector<std::vector<RealTensor> >  d2phi;

  /**
   * Shape function second derivatives in the xi direction.
   */
  std::vector<std::vector<Real> >   d2phidxi2;

  /**
   * Shape function second derivatives in the xi-eta direction.
   */
  std::vector<std::vector<Real> >   d2phidxideta;

  /**
   * Shape function second derivatives in the xi-zeta direction.
   */
  std::vector<std::vector<Real> >   d2phidxidzeta;

  /**
   * Shape function second derivatives in the eta direction.
   */
  std::vector<std::vector<Real> >   d2phideta2;

  /**
   * Shape function second derivatives in the eta-zeta direction.
   */
  std::vector<std::vector<Real> >   d2phidetadzeta;

  /**
   * Shape function second derivatives in the zeta direction.
   */
  std::vector<std::vector<Real> >   d2phidzeta2;

  /**
   * Shape function second derivatives in the x direction.
   */
  std::vector<std::vector<Real> >   d2phidx2;

  /**
   * Shape function second derivatives in the x-y direction.
   */
  std::vector<std::vector<Real> >   d2phidxdy;

  /**
   * Shape function second derivatives in the x-z direction.
   */
  std::vector<std::vector<Real> >   d2phidxdz;

  /**
   * Shape function second derivatives in the y direction.
   */
  std::vector<std::vector<Real> >   d2phidy2;

  /**
   * Shape function second derivatives in the y-z direction.
   */
  std::vector<std::vector<Real> >   d2phidydz;

  /**
   * Shape function second derivatives in the z direction.
   */
  std::vector<std::vector<Real> >   d2phidz2;

#endif




  /**
   * Map for the shape function phi.
   */
  std::vector<std::vector<Real> >   phi_map;

  /**
   * Map for the derivative, d(phi)/d(xi).
   */
  std::vector<std::vector<Real> >   dphidxi_map;

  /**
   * Map for the derivative, d(phi)/d(eta).
   */
  std::vector<std::vector<Real> >   dphideta_map;

  /**
   * Map for the derivative, d(phi)/d(zeta).
   */
  std::vector<std::vector<Real> >   dphidzeta_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Map for the second derivative, d^2(phi)/d(xi)^2.
   */
  std::vector<std::vector<Real> >   d2phidxi2_map;

  /**
   * Map for the second derivative, d^2(phi)/d(xi)d(eta).
   */
  std::vector<std::vector<Real> >   d2phidxideta_map;

  /**
   * Map for the second derivative, d^2(phi)/d(xi)d(zeta).
   */
  std::vector<std::vector<Real> >   d2phidxidzeta_map;

  /**
   * Map for the second derivative, d^2(phi)/d(eta)^2.
   */
  std::vector<std::vector<Real> >   d2phideta2_map;

  /**
   * Map for the second derivative, d^2(phi)/d(eta)d(zeta).
   */
  std::vector<std::vector<Real> >   d2phidetadzeta_map;

  /**
   * Map for the second derivative, d^2(phi)/d(zeta)^2.
   */
  std::vector<std::vector<Real> >   d2phidzeta2_map;

#endif

  /**
   * Map for the side shape functions, psi.
   */
  std::vector<std::vector<Real> >   psi_map;

  /**
   * Map for the derivative of the side functions,
   * d(psi)/d(xi).
   */
  std::vector<std::vector<Real> >   dpsidxi_map;

  /**
   * Map for the derivative of the side function,
   * d(psi)/d(eta).
   */
  std::vector<std::vector<Real> >   dpsideta_map;

  /**
   * Map for the second derivatives (in xi) of the
   * side shape functions.  Useful for computing
   * the curvature at the quadrature points.
   */
  std::vector<std::vector<Real> > d2psidxi2_map;

  /**
   * Map for the second (cross) derivatives in xi, eta
   * of the side shape functions.  Useful for
   * computing the curvature at the quadrature points.
   */
  std::vector<std::vector<Real> > d2psidxideta_map;

  /**
   * Map for the second derivatives (in eta) of the
   * side shape functions.  Useful for computing the
   * curvature at the quadrature points.
   */
  std::vector<std::vector<Real> > d2psideta2_map;

  /**
   * Tangent vectors on boundary at quadrature points.
   */
  std::vector<std::vector<Point> >  tangents;

  /**
   * Normal vectors on boundary at quadrature points
   */
  std::vector<Point>                normals;

  /**
   * The mean curvature (= one half the sum of the principal
   * curvatures) on the boundary at the quadrature points.
   * The mean curvature is a scalar value.
   */
  std::vector<Real>                 curvatures;

  /**
   * Jacobian*Weight values at quadrature points
   */
  std::vector<Real>                 JxW;

  /**
   * The finite element type for this object.  Note that this
   * should be constant for the object.
   */
  const FEType fe_type;

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
  QBase* qrule;

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
		       const FEType& fet) :
  dim(d),
  xyz(),
  dxyzdxi_map(),
  dxyzdeta_map(),
  dxyzdzeta_map(),
  d2xyzdxi2_map(),
  d2xyzdxideta_map(),
  d2xyzdeta2_map(),
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  d2xyzdxidzeta_map(),
  d2xyzdetadzeta_map(),
  d2xyzdzeta2_map(),
#endif
  dxidx_map(),
  dxidy_map(),
  dxidz_map(),
  detadx_map(),
  detady_map(),
  detadz_map(),
  dzetadx_map(),
  dzetady_map(),
  dzetadz_map(),
  calculations_started(false),
  calculate_phi(false),
  calculate_dphi(false),
  calculate_d2phi(false),
  phi_map(),
  dphidxi_map(),
  dphideta_map(),
  dphidzeta_map(),
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  d2phidxi2_map(),
  d2phidxideta_map(),
  d2phidxidzeta_map(),
  d2phideta2_map(),
  d2phidetadzeta_map(),
  d2phidzeta2_map(),
#endif
  psi_map(),
  dpsidxi_map(),
  dpsideta_map(),
  d2psidxi2_map(),
  d2psidxideta_map(),
  d2psideta2_map(),
  tangents(),
  normals(),
  curvatures(),
  JxW(),
  fe_type(fet),
  elem_type(INVALID_ELEM),
  _p_level(0),
  qrule(NULL),
  shapes_on_quadrature(false)
{
}


inline
FEAbstract::~FEAbstract()
{
}

} // namespace libMesh

#endif
