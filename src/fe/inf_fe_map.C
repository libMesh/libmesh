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

  std::unique_ptr<Elem>      base_elem (Base::build_elem (inf_elem));

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
    default:
#ifdef DEBUG
      libmesh_error_msg("Unknown Dim = " << Dim);
#endif
      break;
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
  libmesh_assert(Dim > 0);

  // Start logging the map inversion.
  LOG_SCOPE("inverse_map()", "InfFE");

  // The strategy is:
  // compute the intersection of the line
  // physical_point - origin with the base element,
  // find its internal coordinatels using FE<Dim-1,LAGRANGE>::inverse_map():
  // The radial part can then be computed directly later on.

  // 1.)
  // build a base element to do the map inversion in the base face
  std::unique_ptr<Elem> base_elem (Base::build_elem (inf_elem));

  // The point on the reference element (which we are looking for).
  // start with an invalid guess:
  Point p;
  p(Dim-1)=-2.;

  // 2.)
  // Now find the intersection of a plane represented by the base
  // element nodes and the line given by the origin of the infinite
  // element and the physical point.
  Point intersection;

  // the origin of the infinite element
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
          return p;

        // The intersection should always be between the origin and the physical point.
        // It can become positive close to the border, but should
        // be only very small.
        //  So as in the case above, we can be sufficiently sure here
        // that \p fp is not in this element:
        if (c_factor > 0.01)
          return p;

        // Compute the intersection with
        // {intersection} = {fp} + c*({fp}-{o}).
        intersection.add_scaled(fp,1.+c_factor);
        intersection.add_scaled(o,-c_factor);

        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << Dim);
    }

  // 3.)
  // Now we have the intersection-point (projection of physical point onto base-element).
  // Lets compute its internal coordinates (being p(0) and p(1)):
  p= FE<Dim-1,LAGRANGE>::inverse_map(base_elem.get(), intersection, tolerance, secure);

  // 4.
  // Now that we have the local coordinates in the base,
  // we compute the radial distance with Newton iteration.

  // distance from the physical point to the ifem origin
  const Real fp_o_dist = Point(o-physical_point).norm();

  // the distance from the intersection on the
  // base to the origin
  const Real a_dist = Point(o-intersection).norm();

  // element coordinate in radial direction

  // fp_o_dist is at infinity.
  if (libmesh_isinf(fp_o_dist))
    {
      p(Dim-1)=1;
      return p;
    }

  // when we are somewhere in this element:
  Real v = 0;

  if (T_map == CARTESIAN)
    v = 1.-2.*a_dist/fp_o_dist;
  else
    libmesh_not_implemented();

  // do not put the point back into the element, otherwise the contains_point-function
  // gives false positives!
  //if (v <= -1.-1e-5)
  //  v=-1.;
  //if (v >= 1.)
  //  v=1.-1e-5;

  p(Dim-1)=v;
#ifdef DEBUG
  // first check whether we are in the reference-element:
  if (-1.-1.e-5 < v && v < 1.)
    {
      const Point check = InfFE<Dim,T_radial,T_map>::map (inf_elem, p);
      const Point diff  = physical_point - check;

      if (diff.norm() > tolerance)
        libmesh_warning("WARNING:  diff is " << diff.norm());
    }
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
