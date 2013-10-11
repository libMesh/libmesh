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



#ifndef LIBMESH_FE_MAP_H
#define LIBMESH_FE_MAP_H

#include "libmesh/reference_counted_object.h"
#include "libmesh/point.h"
#include "libmesh/vector_value.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/fe_type.h"
#include "libmesh/auto_ptr.h"

namespace libMesh
{

// forward declarations
class Elem;

  class FEMap
  {
  public:

    FEMap();
    virtual ~FEMap(){}

    static AutoPtr<FEMap> build(FEType fe_type);

    template<unsigned int Dim>
    void init_reference_to_physical_map(const std::vector<Point>& qp,
                                        const Elem* elem);

    /**
     * Compute the jacobian and some other additional
     * data fields at the single point with index p.
     */
    void compute_single_point_map(const unsigned int dim,
                                  const std::vector<Real>& qw,
                                  const Elem* elem,
                                  unsigned int p);

    /**
     * Compute the jacobian and some other additional
     * data fields. Takes the integration weights
     * as input, along with a pointer to the element.
     * The element is assumed to have a constant Jacobian
     */
    virtual void compute_affine_map(const unsigned int dim,
                                    const std::vector<Real>& qw,
                                    const Elem* elem);

    /**
     * Compute the jacobian and some other additional
     * data fields. Takes the integration weights
     * as input, along with a pointer to the element.
     */
    virtual void compute_map(const unsigned int dim,
                             const std::vector<Real>& qw,
                             const Elem* elem);

    /**
     * Same as compute_map, but for a side.  Useful for boundary integration.
     */
    virtual void compute_face_map(int dim,
                                  const std::vector<Real>& qw,
				  const Elem* side);

    /**
     * Same as before, but for an edge.  Useful for some projections.
     */
    void compute_edge_map(int dim,
                          const std::vector<Real>& qw,
			  const Elem* side);

    /**
     * Initalizes the reference to physical element map for a side.
     * This is used for boundary integration.
     */
    template< unsigned int Dim>
    void init_face_shape_functions(const std::vector<Point>& qp,
                                   const Elem* side);

    /**
     * Same as before, but for an edge. This is used for some projection
     * operators.
     */
    template< unsigned int Dim>
    void init_edge_shape_functions(const std::vector<Point>& qp,
                                   const Elem* edge);

    /**
     * @returns the \p xyz spatial locations of the quadrature
     * points on the element.
     */
    const std::vector<Point>& get_xyz() const
    { return xyz; }

    /**
     * @returns the element Jacobian for each quadrature point.
     */
    const std::vector<Real>& get_jacobian() const
    { return jac; }

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
     * @returns the reference to physical map for the side/edge
     */
    const std::vector<std::vector<Real> >& get_psi() const
    { return psi_map; }

    /**
     * @returns the reference to physical map for the element
     */
    const std::vector<std::vector<Real> >& get_phi_map() const
    { return phi_map; }

    /**
     * @returns the reference to physical map derivative
     */
    const std::vector<std::vector<Real> >& get_dphidxi_map() const
    { return dphidxi_map; }

    /**
     * @returns the reference to physical map derivative
     */
    const std::vector<std::vector<Real> >& get_dphideta_map() const
    { return dphideta_map; }

    /**
     * @returns the reference to physical map derivative
     */
    const std::vector<std::vector<Real> >& get_dphidzeta_map() const
    { return dphidzeta_map; }

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
     * Prints the Jacobian times the weight for each quadrature point.
     */
    void print_JxW(std::ostream& os) const;

    /**
     * Prints the spatial location of each quadrature point
     * (on the physical element).
     */
    void print_xyz(std::ostream& os) const;

    /* FIXME: PB: The public functions below break encapsulation! I'm not happy
       about it, but I don't have time to redo the infinite element routines.
       These are needed in inf_fe_boundary.C and inf_fe.C. A proper implementation
       would subclass this class and hide all the radial function business. */
    /**
     * @returns the reference to physical map for the side/edge
     */
    std::vector<std::vector<Real> >& get_psi()
    { return psi_map; }

    /**
     * @returns the reference to physical map derivative for the side/edge
     */
    std::vector<std::vector<Real> >& get_dpsidxi()
    { return dpsidxi_map; }

    /**
     * @returns the reference to physical map derivative for the side/edge
     */
    std::vector<std::vector<Real> >& get_dpsideta()
    { return dpsideta_map; }

    /**
     * @returns the reference to physical map 2nd derivative for the side/edge
     */
    std::vector<std::vector<Real> >& get_d2psidxi2()
    { return d2psidxi2_map; }

    /**
     * @returns the reference to physical map 2nd derivative for the side/edge
     */
    std::vector<std::vector<Real> >& get_d2psidxideta()
    { return d2psidxideta_map; }

    /**
     * @returns the reference to physical map 2nd derivative for the side/edge
     */
    std::vector<std::vector<Real> >& get_d2psideta2()
    { return d2psideta2_map; }

    /**
     * @returns the reference to physical map for the element
     */
    std::vector<std::vector<Real> >& get_phi_map()
    { return phi_map; }

    /**
     * @returns the reference to physical map derivative
     */
    std::vector<std::vector<Real> >& get_dphidxi_map()
    { return dphidxi_map; }

    /**
     * @returns the reference to physical map derivative
     */
    std::vector<std::vector<Real> >& get_dphideta_map()
    { return dphideta_map; }

    /**
     * @returns the reference to physical map derivative
     */
    std::vector<std::vector<Real> >& get_dphidzeta_map()
    { return dphidzeta_map; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    /**
     * @returns the reference to physical map 2nd derivative
     */
    std::vector<std::vector<Real> >& get_d2phidxi2_map()
    { return d2phidxi2_map; }

    /**
     * @returns the reference to physical map 2nd derivative
     */
    std::vector<std::vector<Real> >& get_d2phidxideta_map()
    { return d2phidxideta_map; }

    /**
     * @returns the reference to physical map 2nd derivative
     */
    std::vector<std::vector<Real> >& get_d2phidxidzeta_map()
    { return d2phidxidzeta_map; }

    /**
     * @returns the reference to physical map 2nd derivative
     */
    std::vector<std::vector<Real> >& get_d2phideta2_map()
    { return d2phideta2_map; }

    /**
     * @returns the reference to physical map 2nd derivative
     */
    std::vector<std::vector<Real> >& get_d2phidetadzeta_map()
    { return d2phidetadzeta_map; }

    /**
     * @returns the reference to physical map 2nd derivative
     */
    std::vector<std::vector<Real> >& get_d2phidzeta2_map()
    { return d2phidzeta2_map; }
#endif

    /* FIXME: PB: This function breaks encapsulation! Needed in FE<>::reinit and
       InfFE<>::reinit. Not sure yet if the algorithm can be redone to avoid this. */
    /**
     * @returns writable reference to the element Jacobian times
     * the quadrature weight for each quadrature point.
     */
    std::vector<Real>& get_JxW()
    { return JxW; }

  protected:

    /**
     * A utility function for use by compute_*_map
     */
    void resize_quadrature_map_vectors(const unsigned int dim, unsigned int n_qp);

    /**
     * Used in \p FEMap::compute_map(), which should be
     * be usable in derived classes, and therefore protected.
     * Returns the x value of the pth entry of the dxzydxi_map.
     */
    Real dxdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](0); }

    /**
     * Used in \p FEMap::compute_map(), which should be
     * be usable in derived classes, and therefore protected.
     * Returns the y value of the pth entry of the dxzydxi_map.
     */
    Real dydxi_map(const unsigned int p) const   { return dxyzdxi_map[p](1); }

    /**
     * Used in \p FEMap::compute_map(), which should be
     * be usable in derived classes, and therefore protected.
     * Returns the z value of the pth entry of the dxzydxi_map.
     */
    Real dzdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](2); }

    /**
     * Used in \p FEMap::compute_map(), which should be
     * be usable in derived classes, and therefore protected.
     * Returns the x value of the pth entry of the dxzydeta_map.
     */
    Real dxdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](0); }

    /**
     * Used in \p FEMap::compute_map(), which should be
     * be usable in derived classes, and therefore protected.
     * Returns the y value of the pth entry of the dxzydeta_map.
     */
    Real dydeta_map(const unsigned int p) const  { return dxyzdeta_map[p](1); }

    /**
     * Used in \p FEMap::compute_map(), which should be
     * be usable in derived classes, and therefore protected.
     * Returns the z value of the pth entry of the dxzydeta_map.
     */
    Real dzdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](2); }

    /**
     * Used in \p FEMap::compute_map(), which should be
     * be usable in derived classes, and therefore protected.
     * Returns the x value of the pth entry of the dxzydzeta_map.
     */
    Real dxdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](0); }

    /**
     * Used in \p FEMap::compute_map(), which should be
     * be usable in derived classes, and therefore protected.
     * Returns the y value of the pth entry of the dxzydzeta_map.
     */
    Real dydzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](1); }

    /**
     * Used in \p FEMap::compute_map(), which should be
     * be usable in derived classes, and therefore protected.
     * Returns the z value of the pth entry of the dxzydzeta_map.
     */
    Real dzdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](2); }

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

    /**
     * Map for the shape function phi.
     */
    std::vector<std::vector<Real> > phi_map;

    /**
     * Map for the derivative, d(phi)/d(xi).
     */
    std::vector<std::vector<Real> > dphidxi_map;

    /**
     * Map for the derivative, d(phi)/d(eta).
     */
    std::vector<std::vector<Real> > dphideta_map;

    /**
     * Map for the derivative, d(phi)/d(zeta).
     */
    std::vector<std::vector<Real> > dphidzeta_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

    /**
     * Map for the second derivative, d^2(phi)/d(xi)^2.
     */
    std::vector<std::vector<Real> > d2phidxi2_map;

    /**
     * Map for the second derivative, d^2(phi)/d(xi)d(eta).
     */
    std::vector<std::vector<Real> > d2phidxideta_map;

    /**
     * Map for the second derivative, d^2(phi)/d(xi)d(zeta).
     */
    std::vector<std::vector<Real> > d2phidxidzeta_map;

    /**
     * Map for the second derivative, d^2(phi)/d(eta)^2.
     */
    std::vector<std::vector<Real> > d2phideta2_map;

    /**
     * Map for the second derivative, d^2(phi)/d(eta)d(zeta).
     */
    std::vector<std::vector<Real> > d2phidetadzeta_map;

    /**
     * Map for the second derivative, d^2(phi)/d(zeta)^2.
     */
    std::vector<std::vector<Real> > d2phidzeta2_map;

#endif

    /**
     * Map for the side shape functions, psi.
     */
    std::vector<std::vector<Real> > psi_map;

    /**
     * Map for the derivative of the side functions,
     * d(psi)/d(xi).
     */
    std::vector<std::vector<Real> > dpsidxi_map;

    /**
     * Map for the derivative of the side function,
     * d(psi)/d(eta).
     */
    std::vector<std::vector<Real> > dpsideta_map;

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
    std::vector<std::vector<Point> > tangents;

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
  };

}

#endif // LIBMESH_FE_MAP_H
