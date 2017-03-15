// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#include "libmesh/inf_fe.h"
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/inf_fe_macro.h"
#include "libmesh/libmesh_logging.h"

namespace libMesh
{



// ------------------------------------------------------------
// InfFE static class members concerned with coordinate
// mapping


template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Point InfFE<Dim,T_radial,T_map>::map (const Elem * inf_elem,
                                      const Point & reference_point)
{
  libmesh_assert(inf_elem);
  libmesh_assert_not_equal_to (Dim, 0);

  UniquePtr<Elem>      base_elem (Base::build_elem (inf_elem));

  const Order        radial_mapping_order (Radial::mapping_order());
  const Real         v                    (reference_point(Dim-1));

  // map in the base face
  Point base_point;
  switch (Dim)
    {
    case 1:
      base_point = inf_elem->point(0);
      break;
    case 2:
      base_point = FE<1,LAGRANGE>::map (base_elem.get(), reference_point);
      break;
    case 3:
      base_point = FE<2,LAGRANGE>::map (base_elem.get(), reference_point);
      break;
#ifdef DEBUG
    default:
      libmesh_error_msg("Unknown Dim = " << Dim);
#endif
    }


  // map in the outer node face not necessary. Simply
  // compute the outer_point = base_point + (base_point-origin)
  const Point outer_point (base_point * 2. - inf_elem->origin());

  Point p;

  // there are only two mapping shapes in radial direction
  p.add_scaled (base_point,  InfFE<Dim,INFINITE_MAP,T_map>::eval (v, radial_mapping_order, 0));
  p.add_scaled (outer_point, InfFE<Dim,INFINITE_MAP,T_map>::eval (v, radial_mapping_order, 1));

  return p;
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Point InfFE<Dim,T_radial,T_map>::inverse_map (const Elem * inf_elem,
                                              const Point & physical_point,
                                              const Real tolerance,
                                              const bool secure)
{
  libmesh_assert(inf_elem);
  libmesh_assert_greater_equal (tolerance, 0.);

  // Start logging the map inversion.
  LOG_SCOPE("inverse_map()", "InfFE");

  // 1.)
  // build a base element to do the map inversion in the base face
  UniquePtr<Elem> base_elem (Base::build_elem (inf_elem));

  const ElemType inf_elem_type = inf_elem->type();
  if (inf_elem_type != INFHEX8 &&
      inf_elem_type != INFPRISM6)
    libmesh_error_msg("ERROR: InfFE::inverse_map is currently implemented only for \ninfinite elments of type InfHex8 and InfPrism6.");

  // 2.)
  // just like in FE<Dim-1,LAGRANGE>::inverse_map(): compute
  // the local coordinates, but only in the base element.
  // The radial part can then be computed directly later on.

  // How much did the point on the reference
  // element change by in this Newton step?
  Real inverse_map_error = 0.;

  // The point on the reference element.  This is
  // the "initial guess" for Newton's method.  The
  // centroid seems like a good idea, but computing
  // it is a little more intensive than, say taking
  // the zero point.
  //
  // Convergence should be insensitive of this choice
  // for "good" elements.
  Point p; // the zero point.  No computation required

  // Now find the intersection of a plane represented by the base
  // element nodes and the line given by the origin of the infinite
  // element and the physical point.
  Point intersection;

  // the origin of the infinite lement
  const Point o = inf_elem->origin();

  switch (Dim)
    {
      // unnecessary for 1D
    case 1:
      break;

    case 2:
      libmesh_error_msg("ERROR: InfFE::inverse_map is not yet implemented in 2d");

    case 3:
      {
        // references to the nodal points of the base element
        const Point & p0 = base_elem->point(0);
        const Point & p1 = base_elem->point(1);
        const Point & p2 = base_elem->point(2);

        // a reference to the physical point
        const Point & fp = physical_point;

        // The intersection of the plane and the line is given by
        // can be computed solving a linear 3x3 system
        // a*({p1}-{p0})+b*({p2}-{p0})-c*({fp}-{o})={fp}-{p0}.
        const Real c_factor = -(p1(0)*fp(1)*p0(2)-p1(0)*fp(2)*p0(1)
                                +fp(0)*p1(2)*p0(1)-p0(0)*fp(1)*p1(2)
                                +p0(0)*fp(2)*p1(1)+p2(0)*fp(2)*p0(1)
                                -p2(0)*fp(1)*p0(2)-fp(0)*p2(2)*p0(1)
                                +fp(0)*p0(2)*p2(1)+p0(0)*fp(1)*p2(2)
                                -p0(0)*fp(2)*p2(1)-fp(0)*p0(2)*p1(1)
                                +p0(2)*p2(0)*p1(1)-p0(1)*p2(0)*p1(2)
                                -fp(0)*p1(2)*p2(1)+p2(1)*p0(0)*p1(2)
                                -p2(0)*fp(2)*p1(1)-p1(0)*fp(1)*p2(2)
                                +p2(2)*p1(0)*p0(1)+p1(0)*fp(2)*p2(1)
                                -p0(2)*p1(0)*p2(1)-p2(2)*p0(0)*p1(1)
                                +fp(0)*p2(2)*p1(1)+p2(0)*fp(1)*p1(2))/
          (fp(0)*p1(2)*p0(1)-p1(0)*fp(2)*p0(1)
           +p1(0)*fp(1)*p0(2)-p1(0)*o(1)*p0(2)
           +o(0)*p2(2)*p0(1)-p0(0)*fp(2)*p2(1)
           +p1(0)*o(1)*p2(2)+fp(0)*p0(2)*p2(1)
           -fp(0)*p1(2)*p2(1)-p0(0)*o(1)*p2(2)
           +p0(0)*fp(1)*p2(2)-o(0)*p0(2)*p2(1)
           +o(0)*p1(2)*p2(1)-p2(0)*fp(2)*p1(1)
           +fp(0)*p2(2)*p1(1)-p2(0)*fp(1)*p0(2)
           -o(2)*p0(0)*p1(1)-fp(0)*p0(2)*p1(1)
           +p0(0)*o(1)*p1(2)+p0(0)*fp(2)*p1(1)
           -p0(0)*fp(1)*p1(2)-o(0)*p1(2)*p0(1)
           -p2(0)*o(1)*p1(2)-o(2)*p2(0)*p0(1)
           -o(2)*p1(0)*p2(1)+o(2)*p0(0)*p2(1)
           -fp(0)*p2(2)*p0(1)+o(2)*p2(0)*p1(1)
           +p2(0)*o(1)*p0(2)+p2(0)*fp(1)*p1(2)
           +p2(0)*fp(2)*p0(1)-p1(0)*fp(1)*p2(2)
           +p1(0)*fp(2)*p2(1)-o(0)*p2(2)*p1(1)
           +o(2)*p1(0)*p0(1)+o(0)*p0(2)*p1(1));


        // Check whether the above system is ill-posed.  It should
        // only happen when \p physical_point is not in \p
        // inf_elem. In this case, any point that is not in the
        // element is a valid answer.
        if (libmesh_isinf(c_factor))
          return o;

        // Compute the intersection with
        // {intersection} = {fp} + c*({fp}-{o}).
        intersection.add_scaled(fp,1.+c_factor);
        intersection.add_scaled(o,-c_factor);

        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << Dim);
    }

  // The number of iterations in the map inversion process.
  unsigned int cnt = 0;

  // Newton iteration loop.
  do
    {
      // Increment in current iterate \p p, will be computed.
      // Automatically initialized to all zero.  Note that
      // in 3D, actually only the first two entries are
      // filled by the inverse map, and in 2D only the first.
      Point dp;

      // The form of the map and how we invert it depends
      // on the dimension that we are in.
      switch (Dim)
        {
          // 1D infinite element - no map inversion necessary
        case 1:
          break;

          // 2D infinite element - 1D map inversion
          //
          // In this iteration scheme only search for the local coordinate
          // in xi direction.  Once xi is determined, the radial coordinate eta is
          // uniquely determined, and there is no need to iterate in that direction.
        case 2:
          {
            // Where our current iterate \p p maps to.
            const Point physical_guess = FE<1,LAGRANGE>::map (base_elem.get(), p);

            // How far our current iterate is from the actual point.
            const Point delta = physical_point - physical_guess;

            const Point dxi = FE<1,LAGRANGE>::map_xi (base_elem.get(), p);

            // For details on Newton's method see fe_map.C
            const Real G = dxi*dxi;

            if (secure)
              libmesh_assert_greater (G, 0.);

            const Real Ginv = 1./G;

            const Real  dxidelta = dxi*delta;

            // compute only the first coordinate
            dp(0) = Ginv*dxidelta;

            break;
          }

          // 3D infinite element - 2D map inversion
          //
          // In this iteration scheme only search for the local coordinates
          // in xi and eta direction.  Once xi, eta are determined, the radial
          // coordinate zeta may directly computed.
        case 3:
          {
            // Where our current iterate \p p maps to.
            const Point physical_guess = FE<2,LAGRANGE>::map (base_elem.get(), p);

            // How far our current iterate is from the actual point.
            // const Point delta = physical_point - physical_guess;
            const Point delta = intersection - physical_guess;

            const Point dxi  = FE<2,LAGRANGE>::map_xi  (base_elem.get(), p);
            const Point deta = FE<2,LAGRANGE>::map_eta (base_elem.get(), p);

            // For details on Newton's method see fe_map.C
            const Real
              G11 = dxi*dxi,  G12 = dxi*deta,
              G21 = dxi*deta, G22 = deta*deta;

            const Real det = (G11*G22 - G12*G21);

            if (secure)
              {
                libmesh_assert_greater (det, 0.);
                libmesh_assert_greater (std::abs(det), 1.e-10);
              }

            const Real inv_det = 1./det;

            const Real
              Ginv11 =  G22*inv_det,
              Ginv12 = -G12*inv_det,

              Ginv21 = -G21*inv_det,
              Ginv22 =  G11*inv_det;

            const Real  dxidelta  = dxi*delta;
            const Real  detadelta = deta*delta;

            // compute only the first two coordinates.
            dp(0) = (Ginv11*dxidelta + Ginv12*detadelta);
            dp(1) = (Ginv21*dxidelta + Ginv22*detadelta);

            break;
          }

          // Some other dimension?
        default:
          libmesh_error_msg("Unknown Dim = " << Dim);
        } // end switch(Dim), dp now computed

      // determine the error in computing the local coordinates
      // in the base: ||P_n+1 - P_n||
      inverse_map_error = dp.norm();

      // P_n+1 = P_n + dp
      p.add (dp);

      // Increment the iteration count.
      cnt++;

      // Watch for divergence of Newton's
      // method.
      if (cnt > 10)
        {
          if (secure || !secure)
            {
              libmesh_here();
              {
                libMesh::err << "WARNING: Newton scheme has not converged in "
                             << cnt << " iterations:\n"
                             << "   physical_point="
                             << physical_point
                             << "   dp="
                             << dp
                             << "   p="
                             << p
                             << "   error=" << inverse_map_error
                             << std::endl;
              }
            }

          if (cnt > 20)
            libmesh_error_msg("ERROR: Newton scheme FAILED to converge in " << cnt << " iterations!");
        }
    }
  while (inverse_map_error > tolerance);

  // 4.
  //
  // Now that we have the local coordinates in the base,
  // we compute the radial distance with Newton iteration.

  // distance from the physical point to the ifem origin
  const Real fp_o_dist = Point(o-physical_point).norm();

  // the distance from the intersection on the
  // base to the origin
  const Real a_dist = Point(o-intersection).norm();

  // element coordinate in radial direction
  // here our first guess is 0.
  Real v = 0.;

  // We know the analytic answer for CARTESIAN,
  // other schemes do not yet have a better guess.
  if (T_map == CARTESIAN)
    v = 1.-2.*a_dist/fp_o_dist;

  // the order of the radial mapping
  const Order radial_mapping_order (Radial::mapping_order());

  unsigned int cnt2 = 0;
  inverse_map_error = 0.;

  // Newton iteration in 1-D
  do
    {
      if (v < -1.)
        {
          // in this case, physical_point is not in this element.
          // We therefore give back the best approximation:
          p(Dim-1) = -1;
          return p;
        }

      // the mapping in radial direction
      // note that we only have two mapping functions in
      // radial direction
      const Real r = a_dist * InfFE<Dim,INFINITE_MAP,T_map>::eval (v, radial_mapping_order, 0)
        + 2. * a_dist * InfFE<Dim,INFINITE_MAP,T_map>::eval (v, radial_mapping_order, 1);

      const Real dr = a_dist * InfFE<Dim,INFINITE_MAP,T_map>::eval_deriv (v, radial_mapping_order, 0)
        + 2. * a_dist * InfFE<Dim,INFINITE_MAP,T_map>::eval_deriv (v, radial_mapping_order, 1);

      const Real G = dr*dr;
      const Real Ginv = 1./G;

      const Real delta = fp_o_dist - r;
      const Real drdelta = dr*delta;

      Real dp = Ginv*drdelta;

      // update the radial coordinate
      v += dp;

      // note that v should be smaller than 1,
      // since radial mapping function tends to infinity
      if (v >= 1.)
        v = .9999;

      inverse_map_error = std::abs(dp);

      // increment iteration count
      cnt2 ++;
      if (cnt2 > 20)
        libmesh_error_msg("ERROR: 1D Newton scheme FAILED to converge");
    }
  while (inverse_map_error > tolerance);

  // Set the last coordinate of the Point to v.
  p(Dim-1) = v;

  // If we are in debug mode do a sanity check.  Make sure
  // the point \p p on the reference element actually does
  // map to the point \p physical_point within a tolerance.
#ifdef DEBUG
  const Point check = InfFE<Dim,T_radial,T_map>::map (inf_elem, p);
  const Point diff  = physical_point - check;

  if (diff.norm() > tolerance)
    libmesh_warning("WARNING:  diff is " << diff.norm());
#endif

  return p;
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::inverse_map (const Elem * elem,
                                             const std::vector<Point> & physical_points,
                                             std::vector<Point> &       reference_points,
                                             const Real tolerance,
                                             const bool secure)
{
  // The number of points to find the
  // inverse map of
  const std::size_t n_points = physical_points.size();

  // Resize the vector to hold the points
  // on the reference element
  reference_points.resize(n_points);

  // Find the coordinates on the reference
  // element of each point in physical space
  for (unsigned int p=0; p<n_points; p++)
    reference_points[p] =
      InfFE<Dim,T_radial,T_map>::inverse_map (elem, physical_points[p], tolerance, secure);
}




//--------------------------------------------------------------
// Explicit instantiations using the macro from inf_fe_macro.h
//INSTANTIATE_INF_FE(1,CARTESIAN);

//INSTANTIATE_INF_FE(2,CARTESIAN);

//INSTANTIATE_INF_FE(3,CARTESIAN);

INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, Point, map(const Elem *, const Point &));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, Point, map(const Elem *, const Point &));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, Point, map(const Elem *, const Point &));

INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, Point, inverse_map(const Elem *, const Point &, const Real, const bool));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, Point, inverse_map(const Elem *, const Point &, const Real, const bool));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, Point, inverse_map(const Elem *, const Point &, const Real, const bool));

INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, inverse_map(const Elem *, const std::vector<Point> &, std::vector<Point> &, const Real,  const bool));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, inverse_map(const Elem *, const std::vector<Point> &, std::vector<Point> &, const Real,  const bool));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, inverse_map(const Elem *, const std::vector<Point> &, std::vector<Point> &, const Real,  const bool));


} // namespace libMesh


#endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
