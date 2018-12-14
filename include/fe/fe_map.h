// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/auto_ptr.h" // deprecated

// C++ includes
#include <memory>

namespace libMesh
{

// forward declarations
class Elem;
class Node;

/**
 * Class contained in FE that encapsulates mapping (i.e. from physical
 * space to reference space and vice-versa) quantities and operations.
 *
 * \author Paul Bauman
 * \date 2012
 * \brief Computes finite element mapping function values, gradients, etc.
 */
class FEMap
{
public:

  FEMap();
  virtual ~FEMap(){}

  static std::unique_ptr<FEMap> build(FEType fe_type);

  template<unsigned int Dim>
  void init_reference_to_physical_map(const std::vector<Point> & qp,
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
   * Same as before, but for an edge. This is used for some projection
   * operators.
   */
  template<unsigned int Dim>
  void init_edge_shape_functions(const std::vector<Point> & qp,
                                 const Elem * edge);

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
  const std::vector<Real> & get_jacobian() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return jac; }

  /**
   * \returns The element Jacobian times the quadrature weight for
   * each quadrature point.
   */
  const std::vector<Real> & get_JxW() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return JxW; }

  /**
   * \returns The element tangents in xi-direction at the quadrature
   * points.
   */
  const std::vector<RealGradient> & get_dxyzdxi() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxyzdxi_map; }

  /**
   * \returns The element tangents in eta-direction at the quadrature
   * points.
   */
  const std::vector<RealGradient> & get_dxyzdeta() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxyzdeta_map; }

  /**
   * \returns The element tangents in zeta-direction at the quadrature
   * points.
   */
  const std::vector<RealGradient> & get_dxyzdzeta() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxyzdzeta_map; }

  /**
   * \returns The second partial derivatives in xi.
   */
  const std::vector<RealGradient> & get_d2xyzdxi2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdxi2_map; }

  /**
   * \returns The second partial derivatives in eta.
   */
  const std::vector<RealGradient> & get_d2xyzdeta2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdeta2_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The second partial derivatives in zeta.
   */
  const std::vector<RealGradient> & get_d2xyzdzeta2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdzeta2_map; }

#endif

  /**
   * \returns The second partial derivatives in xi-eta.
   */
  const std::vector<RealGradient> & get_d2xyzdxideta() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdxideta_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The second partial derivatives in xi-zeta.
   */
  const std::vector<RealGradient> & get_d2xyzdxidzeta() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdxidzeta_map; }

  /**
   * \returns The second partial derivatives in eta-zeta.
   */
  const std::vector<RealGradient> & get_d2xyzdetadzeta() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xyzdetadzeta_map; }

#endif

  /**
   * \returns The dxi/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dxidx() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxidx_map; }

  /**
   * \returns The dxi/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dxidy() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxidy_map; }

  /**
   * \returns The dxi/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dxidz() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dxidz_map; }

  /**
   * \returns The deta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_detadx() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return detadx_map; }

  /**
   * \returns The deta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_detady() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return detady_map; }

  /**
   * \returns The deta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_detadz() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return detadz_map; }

  /**
   * \returns The dzeta/dx entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dzetadx() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dzetadx_map; }

  /**
   * \returns The dzeta/dy entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dzetady() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dzetady_map; }

  /**
   * \returns The dzeta/dz entry in the transformation
   * matrix from physical to local coordinates.
   */
  const std::vector<Real> & get_dzetadz() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dzetadz_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Second derivatives of "xi" reference coordinate wrt physical coordinates.
   */
  const std::vector<std::vector<Real>> & get_d2xidxyz2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2xidxyz2_map; }

  /**
   * Second derivatives of "eta" reference coordinate wrt physical coordinates.
   */
  const std::vector<std::vector<Real>> & get_d2etadxyz2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2etadxyz2_map; }

  /**
   * Second derivatives of "zeta" reference coordinate wrt physical coordinates.
   */
  const std::vector<std::vector<Real>> & get_d2zetadxyz2() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2zetadxyz2_map; }
#endif

  /**
   * \returns The reference to physical map for the side/edge
   */
  const std::vector<std::vector<Real>> & get_psi() const
  { return psi_map; }

  /**
   * \returns The reference to physical map for the element
   */
  const std::vector<std::vector<Real>> & get_phi_map() const
  { libmesh_assert(!calculations_started || calculate_xyz);
    calculate_xyz = true; return phi_map; }

  /**
   * \returns The reference to physical map derivative
   */
  const std::vector<std::vector<Real>> & get_dphidxi_map() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphidxi_map; }

  /**
   * \returns The reference to physical map derivative
   */
  const std::vector<std::vector<Real>> & get_dphideta_map() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphideta_map; }

  /**
   * \returns The reference to physical map derivative
   */
  const std::vector<std::vector<Real>> & get_dphidzeta_map() const
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

  /**
   * \returns The curvatures for use in face integration.
   */
  const std::vector<Real> & get_curvatures() const
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return curvatures;}

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
  std::vector<std::vector<Real>> & get_psi()
  { return psi_map; }

  /**
   * \returns The reference to physical map derivative for the side/edge
   */
  std::vector<std::vector<Real>> & get_dpsidxi()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dpsidxi_map; }

  const std::vector<std::vector<Real>> & get_dpsidxi() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dpsidxi_map; }

  /**
   * \returns The reference to physical map derivative for the side/edge
   */
  std::vector<std::vector<Real>> & get_dpsideta()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dpsideta_map; }

  const std::vector<std::vector<Real>> & get_dpsideta() const
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dpsideta_map; }

  /**
   * \returns The reference to physical map 2nd derivative for the side/edge
   */
  std::vector<std::vector<Real>> & get_d2psidxi2()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2psidxi2_map; }

  /**
   * \returns The reference to physical map 2nd derivative for the side/edge
   */
  std::vector<std::vector<Real>> & get_d2psidxideta()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2psidxideta_map; }

  /**
   * \returns The reference to physical map 2nd derivative for the side/edge
   */
  std::vector<std::vector<Real>> & get_d2psideta2()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2psideta2_map; }

  /**
   * \returns The reference to physical map for the element
   */
  std::vector<std::vector<Real>> & get_phi_map()
  { libmesh_assert(!calculations_started || calculate_xyz);
    calculate_xyz = true; return phi_map; }

  /**
   * \returns The reference to physical map derivative
   */
  std::vector<std::vector<Real>> & get_dphidxi_map()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphidxi_map; }

  /**
   * \returns The reference to physical map derivative
   */
  std::vector<std::vector<Real>> & get_dphideta_map()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphideta_map; }

  /**
   * \returns The reference to physical map derivative
   */
  std::vector<std::vector<Real>> & get_dphidzeta_map()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return dphidzeta_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<Real>> & get_d2phidxi2_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidxi2_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<Real>> & get_d2phidxideta_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidxideta_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<Real>> & get_d2phidxidzeta_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidxidzeta_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<Real>> & get_d2phideta2_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phideta2_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<Real>> & get_d2phidetadzeta_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidetadzeta_map; }

  /**
   * \returns The reference to physical map 2nd derivative
   */
  std::vector<std::vector<Real>> & get_d2phidzeta2_map()
  { libmesh_assert(!calculations_started || calculate_d2xyz);
    calculate_d2xyz = true; return d2phidzeta2_map; }
#endif

  /* FIXME: PB: This function breaks encapsulation! Needed in FE<>::reinit and
     InfFE<>::reinit. Not sure yet if the algorithm can be redone to avoid this. */
  /**
   * \returns Writable reference to the element Jacobian times
   * the quadrature weight for each quadrature point.
   */
  std::vector<Real> & get_JxW()
  { libmesh_assert(!calculations_started || calculate_dxyz);
    calculate_dxyz = true; return JxW; }

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
  Real dxdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](0); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The y value of the pth entry of the dxzydxi_map.
   */
  Real dydxi_map(const unsigned int p) const   { return dxyzdxi_map[p](1); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The z value of the pth entry of the dxzydxi_map.
   */
  Real dzdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](2); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The x value of the pth entry of the dxzydeta_map.
   */
  Real dxdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](0); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The y value of the pth entry of the dxzydeta_map.
   */
  Real dydeta_map(const unsigned int p) const  { return dxyzdeta_map[p](1); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The z value of the pth entry of the dxzydeta_map.
   */
  Real dzdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](2); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The x value of the pth entry of the dxzydzeta_map.
   */
  Real dxdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](0); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The y value of the pth entry of the dxzydzeta_map.
   */
  Real dydzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](1); }

  /**
   * Used in \p FEMap::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   *
   * \returns The z value of the pth entry of the dxzydzeta_map.
   */
  Real dzdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](2); }

  /**
   * The spatial locations of the quadrature points
   */
  std::vector<Point> xyz;

  /**
   * Vector of partial derivatives:
   * d(x)/d(xi), d(y)/d(xi), d(z)/d(xi)
   */
  std::vector<RealGradient> dxyzdxi_map;

  /**
   * Vector of partial derivatives:
   * d(x)/d(eta), d(y)/d(eta), d(z)/d(eta)
   */
  std::vector<RealGradient> dxyzdeta_map;

  /**
   * Vector of partial derivatives:
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
  std::vector<Real> dxidx_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(y). Needed for the Jacobian.
   */
  std::vector<Real> dxidy_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(z). Needed for the Jacobian.
   */
  std::vector<Real> dxidz_map;


  /**
   * Map for partial derivatives:
   * d(eta)/d(x). Needed for the Jacobian.
   */
  std::vector<Real> detadx_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(y). Needed for the Jacobian.
   */
  std::vector<Real> detady_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(z). Needed for the Jacobian.
   */
  std::vector<Real> detadz_map;


  /**
   * Map for partial derivatives:
   * d(zeta)/d(x). Needed for the Jacobian.
   */
  std::vector<Real> dzetadx_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(y). Needed for the Jacobian.
   */
  std::vector<Real> dzetady_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(z). Needed for the Jacobian.
   */
  std::vector<Real> dzetadz_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Second derivatives of "xi" reference coordinate wrt physical coordinates.
   * At each qp: (xi_{xx}, xi_{xy}, xi_{xz}, xi_{yy}, xi_{yz}, xi_{zz})
   */
  std::vector<std::vector<Real>> d2xidxyz2_map;

  /**
   * Second derivatives of "eta" reference coordinate wrt physical coordinates.
   * At each qp: (eta_{xx}, eta_{xy}, eta_{xz}, eta_{yy}, eta_{yz}, eta_{zz})
   */
  std::vector<std::vector<Real>> d2etadxyz2_map;

  /**
   * Second derivatives of "zeta" reference coordinate wrt physical coordinates.
   * At each qp: (zeta_{xx}, zeta_{xy}, zeta_{xz}, zeta_{yy}, zeta_{yz}, zeta_{zz})
   */
  std::vector<std::vector<Real>> d2zetadxyz2_map;
#endif

  /**
   * Map for the shape function phi.
   */
  std::vector<std::vector<Real>> phi_map;

  /**
   * Map for the derivative, d(phi)/d(xi).
   */
  std::vector<std::vector<Real>> dphidxi_map;

  /**
   * Map for the derivative, d(phi)/d(eta).
   */
  std::vector<std::vector<Real>> dphideta_map;

  /**
   * Map for the derivative, d(phi)/d(zeta).
   */
  std::vector<std::vector<Real>> dphidzeta_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Map for the second derivative, d^2(phi)/d(xi)^2.
   */
  std::vector<std::vector<Real>> d2phidxi2_map;

  /**
   * Map for the second derivative, d^2(phi)/d(xi)d(eta).
   */
  std::vector<std::vector<Real>> d2phidxideta_map;

  /**
   * Map for the second derivative, d^2(phi)/d(xi)d(zeta).
   */
  std::vector<std::vector<Real>> d2phidxidzeta_map;

  /**
   * Map for the second derivative, d^2(phi)/d(eta)^2.
   */
  std::vector<std::vector<Real>> d2phideta2_map;

  /**
   * Map for the second derivative, d^2(phi)/d(eta)d(zeta).
   */
  std::vector<std::vector<Real>> d2phidetadzeta_map;

  /**
   * Map for the second derivative, d^2(phi)/d(zeta)^2.
   */
  std::vector<std::vector<Real>> d2phidzeta2_map;

#endif

  /**
   * Map for the side shape functions, psi.
   */
  std::vector<std::vector<Real>> psi_map;

  /**
   * Map for the derivative of the side functions,
   * d(psi)/d(xi).
   */
  std::vector<std::vector<Real>> dpsidxi_map;

  /**
   * Map for the derivative of the side function,
   * d(psi)/d(eta).
   */
  std::vector<std::vector<Real>> dpsideta_map;

  /**
   * Map for the second derivatives (in xi) of the
   * side shape functions.  Useful for computing
   * the curvature at the quadrature points.
   */
  std::vector<std::vector<Real>> d2psidxi2_map;

  /**
   * Map for the second (cross) derivatives in xi, eta
   * of the side shape functions.  Useful for
   * computing the curvature at the quadrature points.
   */
  std::vector<std::vector<Real>> d2psidxideta_map;

  /**
   * Map for the second derivatives (in eta) of the
   * side shape functions.  Useful for computing the
   * curvature at the quadrature points.
   */
  std::vector<std::vector<Real>> d2psideta2_map;

  /**
   * Tangent vectors on boundary at quadrature points.
   */
  std::vector<std::vector<Point>> tangents;

  /**
   * Normal vectors on boundary at quadrature points
   */
  std::vector<Point> normals;

  /**
   * The mean curvature (= one half the sum of the principal
   * curvatures) on the boundary at the quadrature points.
   * The mean curvature is a scalar value.
   */
  std::vector<Real> curvatures;

  /**
   * Jacobian values at quadrature points
   */
  std::vector<Real> jac;

  /**
   * Jacobian*Weight values at quadrature points
   */
  std::vector<Real> JxW;

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

  /**
   * Should we calculate mapping hessians?
   */
  mutable bool calculate_d2xyz;

  /**
   * FE classes should be able to reset calculations_started in a few
   * special cases.
   */
  template <unsigned int Dim, FEFamily T>
  friend class FE;

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

}

#endif // LIBMESH_FE_MAP_H
