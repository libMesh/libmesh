// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FE_MAP_H
#define LIBMESH_FE_MAP_H

// libMesh includes
#include "libmesh/reference_counted_object.h"
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/fe_type.h"
#include "libmesh/fe_shim.h"

// C++ includes
#include <memory>
#include <type_traits>

namespace libMesh
{

// forward declarations
template <typename>
class ElemTempl;
template <typename>
class NodeTempl;

/**
 * Class contained in FE that encapsulates mapping (i.e. from physical
 * space to reference space and vice-versa) quantities and operations.
 *
 * \author Paul Bauman
 * \date 2012
 * \brief Computes finite element mapping function values, gradients, etc.
 */
template <typename RealType = Real>
class FEMapTempl
{
public:
  typedef FEMapTempl<RealType> FEMap;
  typedef NodeTempl<RealType> Node;
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef VectorValue<RealType> RealVectorValue;
  typedef TensorValue<RealType> RealTensor;
  typedef TensorValue<RealType> RealTensorValue;

  FEMapTempl(Real jtol = 0);
  virtual ~FEMapTempl(){}

  static std::unique_ptr<FEMap> build(FEType fe_type);

  static FEFamily map_fe_type(const Elem & elem);

  template<unsigned int Dim>
  void init_reference_to_physical_map(const std::vector<Point> & qp,
                                      const Elem * elem);

  /**
   * Method for casting real points into AD type points for use when initializing
   * shape functions based on a quadrature rule, which is always going to supply the same
   * reference points regardless of mesh dislacements
   */
  template<unsigned int Dim, typename T = RealType,
           typename std::enable_if<!std::is_same<T,Real>::value,int>::type = 0>
  void init_reference_to_physical_map(const std::vector<PointTempl<Real>> & qp,
                                      const Elem * elem);

  /**
   * Compute the jacobian and some other additional data fields at the
   * single point with index p.  Takes the integration weights as
   * input, along with a pointer to the element and a list of points
   * that contribute to the element.  Also takes a boolean flag
   * telling whether second derivatives should actually be computed.
   */
  void compute_single_point_map(const unsigned int dim,
                                const std::vector<Real> & qw,
                                const Elem * elem,
                                unsigned int p,
                                const std::vector<const Node *> & elem_nodes,
                                bool compute_second_derivatives);

  /**
   * Compute the jacobian and some other additional
   * data fields. Takes the integration weights
   * as input, along with a pointer to the element.
   * The element is assumed to have a constant Jacobian
   */
  virtual void compute_affine_map(const unsigned int dim,
                                  const std::vector<Real> & qw,
                                  const Elem * elem);

  /**
   * Assign a fake jacobian and some other additional data fields.
   * Takes the integration weights as input.  For use on non-element
   * evaluations.
   */
  virtual void compute_null_map(const unsigned int dim,
                                const std::vector<Real> & qw);


  /**
   * Compute the jacobian and some other additional
   * data fields. Takes the integration weights
   * as input, along with a pointer to the element.
   * Also takes a boolean parameter indicating whether second
   * derivatives need to be calculated, allowing us to potentially
   * skip unnecessary, expensive computations.
   */
  virtual void compute_map(const unsigned int dim,
                           const std::vector<Real> & qw,
                           const Elem * elem,
                           bool calculate_d2phi);

  /**
   * Same as compute_map, but for a side.  Useful for boundary integration.
   */
  virtual void compute_face_map(int dim,
                                const std::vector<Real> & qw,
                                const Elem * side);

  /**
   * Same as before, but for an edge.  Useful for some projections.
   */
  void compute_edge_map(int dim,
                        const std::vector<Real> & qw,
                        const Elem * side);

  /**
   * Initializes the reference to physical element map for a side.
   * This is used for boundary integration.
   */
  template<unsigned int Dim>
  void init_face_shape_functions(const std::vector<Point> & qp,
                                 const Elem * side);

  /**
   * Method for casting real points into AD type points for use when initializing
   * shape functions based on a quadrature rule, which is always going to supply the same
   * reference points regardless of mesh dislacements
   */
  template<unsigned int Dim, typename T = RealType,
           typename std::enable_if<!std::is_same<T,Real>::value,int>::type = 0>
  void init_face_shape_functions(const std::vector<PointTempl<Real>> & qp,
                                 const Elem * side);

  /**
   * Same as before, but for an edge. This is used for some projection
   * operators.
   */
  template<unsigned int Dim>
  void init_edge_shape_functions(const std::vector<Point> & qp,
                                 const Elem * edge);

 /**
   * Method for casting real points into AD type points for use when initializing
   * shape functions based on a quadrature rule, which is always going to supply the same
   * reference points regardless of mesh dislacements
   */
  template<unsigned int Dim, typename T = RealType,
           typename std::enable_if<!std::is_same<T,Real>::value,int>::type = 0>
  void init_edge_shape_functions(const std::vector<PointTempl<Real>> & qp,
                                 const Elem * edge);

  /**
   * \returns The location (in physical space) of the point
   * \p p located on the reference element.
   */
  static Point map (const unsigned int dim,
                    const Elem * elem,
                    const Point & reference_point);

  /**
   * \returns component \p j of d(xyz)/d(xi eta zeta) (in physical
   * space) of the point \p p located on the reference element.
   */
  static Point map_deriv (const unsigned int dim,
                          const Elem * elem,
                          const unsigned int j,
                          const Point & reference_point);

  /**
   * \returns The location (on the reference element) of the
   * point \p p located in physical space.  This function requires
   * inverting the (possibly nonlinear) transformation map, so
   * it is not trivial. The optional parameter \p tolerance defines
   * how close is "good enough."  The map inversion iteration
   * computes the sequence \f$ \{ p_n \} \f$, and the iteration is
   * terminated when \f$ \|p - p_n\| < \mbox{\texttt{tolerance}} \f$
   * The parameter secure (always assumed false in non-debug mode)
   * switches on integrity-checks on the mapped points.
   */
  static Point inverse_map (const unsigned int dim,
                            const Elem * elem,
                            const Point & p,
                            const Real tolerance = TOLERANCE,
                            const bool secure = true);

  /**
   * Takes a number points in physical space (in the \p
   * physical_points vector) and finds their location on the reference
   * element for the input element \p elem.  The values on the
   * reference element are returned in the vector \p
   * reference_points. The optional parameter \p tolerance defines how
   * close is "good enough."  The map inversion iteration computes the
   * sequence \f$ \{ p_n \} \f$, and the iteration is terminated when
   * \f$ \|p - p_n\| < \mbox{\texttt{tolerance}} \f$
   * The parameter secure (always assumed false in non-debug mode)
   * switches on integrity-checks on the mapped points.
   */
  static void inverse_map (unsigned int dim,
                           const Elem * elem,
                           const std::vector<Point> & physical_points,
                           std::vector<Point> &       reference_points,
                           const Real tolerance = TOLERANCE,
                           const bool secure = true);

  /**
   * \returns The \p xyz spatial locations of the quadrature
   * points on the element.
   */
  const std::vector<Point> & get_xyz() const
  { libmesh_assert(!calculations_started || calculate_xyz);
    calculate_xyz = true; return xyz; }

  /**
   * \returns The element Jacobian for each quadrature point.
   */
  const std::vector<RealType> & get_jacobian() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return jac; }

  /**
   * \returns The element Jacobian times the quadrature weight for
   * each quadrature point.
   */
  const std::vector<RealType> & get_JxW() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return JxW; }

  /**
   * \returns The element tangents in xi-direction at the quadrature
   * points.
   */
  const std::vector<RealVectorValue> & get_dxyzdxi() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxyzdxi_map; }

  /**
   * \returns The element tangents in eta-direction at the quadrature
   * points.
   */
  const std::vector<RealVectorValue> & get_dxyzdeta() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxyzdeta_map; }

  /**
   * \returns The element tangents in zeta-direction at the quadrature
   * points.
   */
  const std::vector<RealVectorValue> & get_dxyzdzeta() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxyzdzeta_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The second partial derivatives in xi.
   */
  const std::vector<RealVectorValue> & get_d2xyzdxi2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdxi2_map; }

  /**
   * \returns The second partial derivatives in eta.
   */
  const std::vector<RealVectorValue> & get_d2xyzdeta2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdeta2_map; }

  /**
   * \returns The second partial derivatives in zeta.
   */
  const std::vector<RealVectorValue> & get_d2xyzdzeta2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdzeta2_map; }

  /**
   * \returns The second partial derivatives in xi-eta.
   */
  const std::vector<RealVectorValue> & get_d2xyzdxideta() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdxideta_map; }

  /**
   * \returns The second partial derivatives in xi-zeta.
   */
  const std::vector<RealVectorValue> & get_d2xyzdxidzeta() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdxidzeta_map; }

  /**
   * \returns The second partial derivatives in eta-zeta.
   */
  const std::vector<RealVectorValue> & get_d2xyzdetadzeta() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdetadzeta_map; }

#endif

  /**
   * \returns The dxi/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<RealType> & get_dxidx() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxidx_map; }

  /**
   * \returns The dxi/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<RealType> & get_dxidy() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxidy_map; }

  /**
   * \returns The dxi/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<RealType> & get_dxidz() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxidz_map; }

  /**
   * \returns The deta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<RealType> & get_detadx() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return detadx_map; }

  /**
   * \returns The deta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<RealType> & get_detady() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return detady_map; }

  /**
   * \returns The deta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<RealType> & get_detadz() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return detadz_map; }

  /**
   * \returns The dzeta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<RealType> & get_dzetadx() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dzetadx_map; }

  /**
   * \returns The dzeta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<RealType> & get_dzetady() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dzetady_map; }

  /**
   * \returns The dzeta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<RealType> & get_dzetadz() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dzetadz_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Second derivatives of "xi" reference coordinate wrt physical coordinates.
   */
  const std::vector<std::vector<RealType>> & get_d2xidxyz2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xidxyz2_map; }

  /**
   * Second derivatives of "eta" reference coordinate wrt physical coordinates.
   */
  const std::vector<std::vector<RealType>> & get_d2etadxyz2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2etadxyz2_map; }

  /**
   * Second derivatives of "zeta" reference coordinate wrt physical coordinates.
   */
  const std::vector<std::vector<RealType>> & get_d2zetadxyz2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2zetadxyz2_map; }
#endif

  /**
   * \returns The reference to physical map for the side/edge
   */
  const std::vector<std::vector<RealType>> & get_psi() const
  { return psi_map; }

  /**
   * \returns The reference to physical map for the element
   */
  const std::vector<std::vector<RealType>> & get_phi_map() const
  { libmesh_assert(!calculations_started || calculate_xyz);
    calculate_xyz = true; return phi_map; }

  /**
   * \returns The reference to physical map derivative
   */
  const std::vector<std::vector<RealType>> & get_dphidxi_map() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphidxi_map; }

  /**
   * \returns The reference to physical map derivative
   */
  const std::vector<std::vector<RealType>> & get_dphideta_map() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphideta_map; }

  /**
   * \returns The reference to physical map derivative
   */
  const std::vector<std::vector<RealType>> & get_dphidzeta_map() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphidzeta_map; }

  /**
   * \returns The tangent vectors for face integration.
   */
  const std::vector<std::vector<Point>> & get_tangents() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return tangents; }

  /**
   * \returns The outward pointing normal vectors for face integration.
   */
  const std::vector<Point> & get_normals() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return normals; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The curvatures for use in face integration.
   */
  const std::vector<RealType> & get_curvatures() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return curvatures;}

#endif

  /**
   * Prints the Jacobian times the weight for each quadrature point.
   */
  void print_JxW(std::ostream & os) const;

  /**
   * Prints the spatial location of each quadrature point
   * (on the physical element).
   */
  void print_xyz(std::ostream & os) const;

  /* FIXME: PB: The public functions below break encapsulation! I'm not happy
     about it, but I don't have time to redo the infinite element routines.
     These are needed in inf_fe_boundary.C and inf_fe.C. A proper implementation
     would subclass this class and hide all the radial function business. */
  /**
   * \returns The reference to physical map for the side/edge
   */
  std::vector<std::vector<RealType>> & get_psi()
  { return psi_map; }

  /**
   * \returns The reference to physical map derivative for the side/edge
   */
  std::vector<std::vector<RealType>> & get_dpsidxi()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dpsidxi_map; }

  const std::vector<std::vector<RealType>> & get_dpsidxi() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dpsidxi_map; }

  /**
   * \returns The reference to physical map derivative for the side/edge
   */
  std::vector<std::vector<RealType>> & get_dpsideta()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dpsideta_map; }

  const std::vector<std::vector<RealType>> & get_dpsideta() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dpsideta_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The reference to physical map 2nd derivative for the side/edge
   */
  std::vector<std::vector<RealType>> & get_d2psidxi2()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2psidxi2_map; }

  /**
   * \returns const reference to physical map 2nd derivative for the side/edge
   */
  const std::vector<std::vector<RealType>> & get_d2psidxi2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2psidxi2_map; }

  /**
   * \returns The reference to physical map 2nd derivative for the side/edge
   */
  std::vector<std::vector<RealType>> & get_d2psidxideta()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2psidxideta_map; }

  /**
   * \returns const reference to physical map 2nd derivative for the side/edge
   */
  const std::vector<std::vector<RealType>> & get_d2psidxideta() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2psidxideta_map; }

  /**
   * \returns The reference to physical map 2nd derivative for the side/edge
   */
  std::vector<std::vector<RealType>> & get_d2psideta2()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2psideta2_map; }


  /**
   * \returns const reference to physical map 2nd derivative for the side/edge
   */
  const std::vector<std::vector<RealType>> & get_d2psideta2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2psideta2_map; }

#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The reference to physical map for the element
   */
  std::vector<std::vector<RealType>> & get_phi_map()
  { libmesh_assert(!calculations_started || calculate_xyz);
    calculate_xyz = true; return phi_map; }

  /**
   * \returns The reference to physical map derivative
   */
  std::vector<std::vector<RealType>> & get_dphidxi_map()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphidxi_map; }

  /**
   * \returns The reference to physical map derivative
   */
  std::vector<std::vector<RealType>> & get_dphideta_map()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphideta_map; }

  /**
   * \returns The reference to physical map derivative
   */
  std::vector<std::vector<RealType>> & get_dphidzeta_map()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphidzeta_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<RealType>> & get_d2phidxi2_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidxi2_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<RealType>> & get_d2phidxideta_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidxideta_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<RealType>> & get_d2phidxidzeta_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidxidzeta_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<RealType>> & get_d2phideta2_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phideta2_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<RealType>> & get_d2phidetadzeta_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidetadzeta_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<RealType>> & get_d2phidzeta2_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidzeta2_map; }
#endif

  /* FIXME: PB: This function breaks encapsulation! Needed in FE<>::reinit and
     InfFE<>::reinit. Not sure yet if the algorithm can be redone to avoid this. */
  /**
   * \returns Writable reference to the element Jacobian times
   * the quadrature weight for each quadrature point.
   */
  std::vector<RealType> & get_JxW()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return JxW; }

  /**
   * Set the Jacobian tolerance used for determining when the mapping fails. The mapping is
   * determined to fail if jac <= jacobian_tolerance.
   */
  void set_jacobian_tolerance(Real tol) { jacobian_tolerance = tol; }

protected:

  /**
   * Determine which values are to be calculated
   */
  void determine_calculations() {
    calculations_started = true;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    // Second derivative calculations currently have first derivative
    // calculations as a prerequisite
    if (calculate_d2xyz)
      calculate_dxyz = true;
#endif
  }

  /**
   * A utility function for use by compute_*_map
   */
  void resize_quadrature_map_vectors(const unsigned int dim, unsigned int n_qp);

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The x value of the pth entry of the dxzydxi_map.
   */
  RealType dxdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](0); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The y value of the pth entry of the dxzydxi_map.
   */
  RealType dydxi_map(const unsigned int p) const   { return dxyzdxi_map[p](1); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The z value of the pth entry of the dxzydxi_map.
   */
  RealType dzdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](2); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The x value of the pth entry of the dxzydeta_map.
   */
  RealType dxdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](0); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The y value of the pth entry of the dxzydeta_map.
   */
  RealType dydeta_map(const unsigned int p) const  { return dxyzdeta_map[p](1); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The z value of the pth entry of the dxzydeta_map.
   */
  RealType dzdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](2); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The x value of the pth entry of the dxzydzeta_map.
   */
  RealType dxdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](0); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The y value of the pth entry of the dxzydzeta_map.
   */
  RealType dydzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](1); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The z value of the pth entry of the dxzydzeta_map.
   */
  RealType dzdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](2); }

  /**
   * The spatial locations of the quadrature points
   */
  std::vector<Point> xyz;

  /**
   * Vector of partial derivatives:
   * d(x)/d(xi), d(y)/d(xi), d(z)/d(xi)
   */
  std::vector<RealVectorValue> dxyzdxi_map;

  /**
   * Vector of partial derivatives:
   * d(x)/d(eta), d(y)/d(eta), d(z)/d(eta)
   */
  std::vector<RealVectorValue> dxyzdeta_map;

  /**
   * Vector of partial derivatives:
   * d(x)/d(zeta), d(y)/d(zeta), d(z)/d(zeta)
   */
  std::vector<RealVectorValue> dxyzdzeta_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Vector of second partial derivatives in xi:
   * d^2(x)/d(xi)^2, d^2(y)/d(xi)^2, d^2(z)/d(xi)^2
   */
  std::vector<RealVectorValue> d2xyzdxi2_map;

  /**
   * Vector of mixed second partial derivatives in xi-eta:
   * d^2(x)/d(xi)d(eta) d^2(y)/d(xi)d(eta) d^2(z)/d(xi)d(eta)
   */
  std::vector<RealVectorValue> d2xyzdxideta_map;

  /**
   * Vector of second partial derivatives in eta:
   * d^2(x)/d(eta)^2
   */
  std::vector<RealVectorValue> d2xyzdeta2_map;

  /**
   * Vector of second partial derivatives in xi-zeta:
   * d^2(x)/d(xi)d(zeta), d^2(y)/d(xi)d(zeta), d^2(z)/d(xi)d(zeta)
   */
  std::vector<RealVectorValue> d2xyzdxidzeta_map;

  /**
   * Vector of mixed second partial derivatives in eta-zeta:
   * d^2(x)/d(eta)d(zeta) d^2(y)/d(eta)d(zeta) d^2(z)/d(eta)d(zeta)
   */
  std::vector<RealVectorValue> d2xyzdetadzeta_map;

  /**
   * Vector of second partial derivatives in zeta:
   * d^2(x)/d(zeta)^2
   */
  std::vector<RealVectorValue> d2xyzdzeta2_map;

#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Map for partial derivatives:
   * d(xi)/d(x). Needed for the Jacobian.
   */
  std::vector<RealType> dxidx_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(y). Needed for the Jacobian.
   */
  std::vector<RealType> dxidy_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(z). Needed for the Jacobian.
   */
  std::vector<RealType> dxidz_map;


  /**
   * Map for partial derivatives:
   * d(eta)/d(x). Needed for the Jacobian.
   */
  std::vector<RealType> detadx_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(y). Needed for the Jacobian.
   */
  std::vector<RealType> detady_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(z). Needed for the Jacobian.
   */
  std::vector<RealType> detadz_map;


  /**
   * Map for partial derivatives:
   * d(zeta)/d(x). Needed for the Jacobian.
   */
  std::vector<RealType> dzetadx_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(y). Needed for the Jacobian.
   */
  std::vector<RealType> dzetady_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(z). Needed for the Jacobian.
   */
  std::vector<RealType> dzetadz_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Second derivatives of "xi" reference coordinate wrt physical coordinates.
   * At each qp: (xi_{xx}, xi_{xy}, xi_{xz}, xi_{yy}, xi_{yz}, xi_{zz})
   */
  std::vector<std::vector<RealType>> d2xidxyz2_map;

  /**
   * Second derivatives of "eta" reference coordinate wrt physical coordinates.
   * At each qp: (eta_{xx}, eta_{xy}, eta_{xz}, eta_{yy}, eta_{yz}, eta_{zz})
   */
  std::vector<std::vector<RealType>> d2etadxyz2_map;

  /**
   * Second derivatives of "zeta" reference coordinate wrt physical coordinates.
   * At each qp: (zeta_{xx}, zeta_{xy}, zeta_{xz}, zeta_{yy}, zeta_{yz}, zeta_{zz})
   */
  std::vector<std::vector<RealType>> d2zetadxyz2_map;
#endif

  /**
   * Map for the shape function phi.
   */
  std::vector<std::vector<RealType>> phi_map;

  /**
   * Map for the derivative, d(phi)/d(xi).
   */
  std::vector<std::vector<RealType>> dphidxi_map;

  /**
   * Map for the derivative, d(phi)/d(eta).
   */
  std::vector<std::vector<RealType>> dphideta_map;

  /**
   * Map for the derivative, d(phi)/d(zeta).
   */
  std::vector<std::vector<RealType>> dphidzeta_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Map for the second derivative, d^2(phi)/d(xi)^2.
   */
  std::vector<std::vector<RealType>> d2phidxi2_map;

  /**
   * Map for the second derivative, d^2(phi)/d(xi)d(eta).
   */
  std::vector<std::vector<RealType>> d2phidxideta_map;

  /**
   * Map for the second derivative, d^2(phi)/d(xi)d(zeta).
   */
  std::vector<std::vector<RealType>> d2phidxidzeta_map;

  /**
   * Map for the second derivative, d^2(phi)/d(eta)^2.
   */
  std::vector<std::vector<RealType>> d2phideta2_map;

  /**
   * Map for the second derivative, d^2(phi)/d(eta)d(zeta).
   */
  std::vector<std::vector<RealType>> d2phidetadzeta_map;

  /**
   * Map for the second derivative, d^2(phi)/d(zeta)^2.
   */
  std::vector<std::vector<RealType>> d2phidzeta2_map;

#endif

  /**
   * Map for the side shape functions, psi.
   */
  std::vector<std::vector<RealType>> psi_map;

  /**
   * Map for the derivative of the side functions,
   * d(psi)/d(xi).
   */
  std::vector<std::vector<RealType>> dpsidxi_map;

  /**
   * Map for the derivative of the side function,
   * d(psi)/d(eta).
   */
  std::vector<std::vector<RealType>> dpsideta_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Map for the second derivatives (in xi) of the
   * side shape functions.  Useful for computing
   * the curvature at the quadrature points.
   */
  std::vector<std::vector<RealType>> d2psidxi2_map;

  /**
   * Map for the second (cross) derivatives in xi, eta
   * of the side shape functions.  Useful for
   * computing the curvature at the quadrature points.
   */
  std::vector<std::vector<RealType>> d2psidxideta_map;

  /**
   * Map for the second derivatives (in eta) of the
   * side shape functions.  Useful for computing the
   * curvature at the quadrature points.
   */
  std::vector<std::vector<RealType>> d2psideta2_map;

#endif

  /**
   * Tangent vectors on boundary at quadrature points.
   */
  std::vector<std::vector<Point>> tangents;

  /**
   * Normal vectors on boundary at quadrature points
   */
  std::vector<Point> normals;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * The mean curvature (= one half the sum of the principal
   * curvatures) on the boundary at the quadrature points.
   * The mean curvature is a scalar value.
   */
  std::vector<RealType> curvatures;

#endif

  /**
   * Jacobian values at quadrature points
   */
  std::vector<RealType> jac;

  /**
   * Jacobian*Weight values at quadrature points
   */
  std::vector<RealType> JxW;

  /**
   * Have calculations with this object already been started?
   * Then all get_* functions should already have been called.
   */
  mutable bool calculations_started;

  /**
   * Should we calculate physical point locations?
   */
  mutable bool calculate_xyz;

  /**
   * Should we calculate mapping gradients?
   */
  mutable bool calculate_dxyz;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Should we calculate mapping hessians?
   */
  mutable bool calculate_d2xyz;

#endif

  /**
   * FE classes should be able to reset calculations_started in a few
   * special cases.
   */
  template <unsigned int, FEFamily, typename>
  friend class FE;
  template <unsigned int, typename>
  friend struct FEMapShapeFuncShim;
  template <unsigned int, FEFamily, typename>
  friend struct FEReinitShim;
  template <unsigned int, FEFamily, typename>
  friend struct FEEdgeReinitShim;
  template <unsigned int, FEFamily, typename>
  friend struct FESideMapShim;

  /**
   * The Jacobian tolerance used for determining when the mapping fails. The mapping is
   * determined to fail if jac <= jacobian_tolerance. If not set by the user, this number
   * defaults to 0
   */
  Real jacobian_tolerance;

private:
  /**
   * A helper function used by FEMap::compute_single_point_map() to
   * compute second derivatives of the inverse map.
   */
  void compute_inverse_map_second_derivs(unsigned p);

  /**
   * Work vector for compute_affine_map()
   */
  std::vector<const Node *> _elem_nodes;
};

template <typename RealType>
template <unsigned int Dim>
void FEMapTempl<RealType>::init_face_shape_functions(const std::vector<Point> & point, const Elem * elem)
{
  FEMapShapeFuncShim<Dim, RealType>::init_face_shape_functions(*this, point, elem);
}

template <typename RealType>
template <unsigned int Dim>
void FEMapTempl<RealType>::init_edge_shape_functions(const std::vector<Point> & point, const Elem * elem)
{
  FEMapShapeFuncShim<Dim, RealType>::init_edge_shape_functions(*this, point, elem);
}

template <typename RealType>
template <unsigned int Dim, typename T,
          typename std::enable_if<!std::is_same<T,Real>::value,int>::type>
void FEMapTempl<RealType>::init_face_shape_functions(const std::vector<PointTempl<Real>> & qp,
                                                     const Elem * side)
{
  std::vector<Point> my_qp(qp.size());

  for (std::size_t point = 0; point < qp.size(); ++point)
    my_qp[point] = qp[point];

  this->init_face_shape_functions(my_qp, side);
}

template <typename RealType>
template <unsigned int Dim, typename T,
          typename std::enable_if<!std::is_same<T,Real>::value,int>::type>
void FEMapTempl<RealType>::init_edge_shape_functions(const std::vector<PointTempl<Real>> & qp,
                                                     const Elem * edge)
{
  std::vector<Point> my_qp(qp.size());

  for (std::size_t point = 0; point < qp.size(); ++point)
    my_qp[point] = qp[point];

  this->init_edge_shape_functions(my_qp, edge);
}

template <typename RealType>
template <unsigned int Dim, typename T,
          typename std::enable_if<!std::is_same<T,Real>::value,int>::type>
void FEMapTempl<RealType>::init_reference_to_physical_map(const std::vector<PointTempl<Real>> & qp,
                                                          const Elem * elem)
{
  std::vector<Point> my_qp(qp.size());

  for (std::size_t point = 0; point < qp.size(); ++point)
    my_qp[point] = qp[point];

  this->init_reference_to_physical_map(my_qp, elem);
}

typedef FEMapTempl<Real> FEMap;
}

#endif // LIBMESH_FE_MAP_H
