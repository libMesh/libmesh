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



// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#include "libmesh/inf_fe_map.h"
#include "libmesh/inf_fe.h"
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/inf_fe_macro.h"
#include "libmesh/libmesh_logging.h"

#include "libmesh/type_tensor.h"

namespace libMesh
{

// ------------------------------------------------------------
// InfFE static class members concerned with coordinate
// mapping


Point InfFEMap::map (const unsigned int dim,
                     const Elem * inf_elem,
                     const Point & reference_point)
{
  libmesh_assert(inf_elem);
  libmesh_assert_not_equal_to (dim, 0);

  std::unique_ptr<Elem>      base_elem (InfFEBase::build_elem (inf_elem));

  const Order        radial_mapping_order (InfFERadial::mapping_order());
  const Real         v                    (reference_point(dim-1));

  // map in the base face
  Point base_point;
  switch (dim)
    {
    case 1:
      base_point = inf_elem->point(0);
      break;
    case 2:
    case 3:
      base_point = FEMap::map (dim-1, base_elem.get(), reference_point);
      break;
    default:
#ifdef DEBUG
      libmesh_error_msg("Unknown dim = " << dim);
#endif
      break;
    }

  // This is the same as the algorithm used below,
  // but is more explicit in the actual form in the end.
  //
  // NOTE: the form used below can be implemented to yield
  // more general/flexible mappings, but the current form is
  // used e.g. for \p inverse_map() and \p reinit() explicitly.
  return (base_point-inf_elem->origin())*2./(1.-v)+inf_elem->origin();

  // map in the outer node face not necessary. Simply
  // compute the outer_point = base_point + (base_point-origin)
  const Point outer_point (base_point * 2. - inf_elem->origin());

  Point p;

  // there are only two mapping shapes in radial direction
  p.add_scaled (base_point,  eval (v, radial_mapping_order, 0));
  p.add_scaled (outer_point, eval (v, radial_mapping_order, 1));

  return p;
}



Point InfFEMap::inverse_map (const unsigned int dim,
                             const Elem * inf_elem,
                             const Point & physical_point,
                             const Real tolerance,
                             const bool secure)
{
  libmesh_assert(inf_elem);
  libmesh_assert_greater_equal (tolerance, 0.);
  libmesh_assert(dim > 0);

  // Start logging the map inversion.
  LOG_SCOPE("inverse_map()", "InfFEMap");

  // The strategy is:
  // compute the intersection of the line
  // physical_point - origin with the base element,
  // find its internal coordinatels using FEMap::inverse_map():
  // The radial part can then be computed directly later on.

  // 1.)
  // build a base element to do the map inversion in the base face
  std::unique_ptr<Elem> base_elem (InfFEBase::build_elem (inf_elem));

  // The point on the reference element (which we are looking for).
  // start with an invalid guess:
  Point p;
  p(dim-1)=-2.;

  // 2.)
  // Now find the intersection of a plane represented by the base
  // element nodes and the line given by the origin of the infinite
  // element and the physical point.
  Point intersection;

  // the origin of the infinite element
  const Point o = inf_elem->origin();

  switch (dim)
    {
      // unnecessary for 1D
    case 1:
      break;

    case 2:
      libmesh_error_msg("ERROR: InfFE::inverse_map is not yet implemented in 2d");

    case 3:
      {
        const Point xi ( base_elem->point(1) - base_elem->point(0));
        const Point eta( base_elem->point(2) - base_elem->point(0));
        const Point zeta( physical_point - o);

        // normal vector of the base elements plane
        Point n(xi.cross(eta));
        Real c_factor = (base_elem->point(0) - o)*n/(zeta*n) - 1.;

        // Check whether the above system is ill-posed.  It should
        // only happen when \p physical_point is not in \p
        // inf_elem. In this case, any point that is not in the
        // element is a valid answer.
        if (libmesh_isinf(c_factor))
          return p;

        // Compute the intersection with
        // {intersection} = {physical_point} + c*({physical_point}-{o}).
        intersection.add_scaled(physical_point,1.+c_factor);
        intersection.add_scaled(o,-c_factor);

        // For non-planar elements, the intersecting point is obtained via Newton-iteration
        if (!base_elem->has_affine_map())
          {
            unsigned int iter_max = 20;

            // the number of shape functions needed for the base_elem
            unsigned int n_sf = FE<2,LAGRANGE>::n_shape_functions(base_elem->type(),base_elem->default_order());

            // guess base element coordinates: p=xi,eta,0
            // in first iteration, never run with 'secure' to avoid false warnings.
            Point ref_point= FEMap::inverse_map(dim-1, base_elem.get(), intersection,
                                                tolerance, false);

            // Newton iteration
            for(unsigned int it=0; it<iter_max; ++it)
              {
                // Get the shape function and derivative values at the reference coordinate
                // phi.size() == dphi.size()
                Point dxyz_dxi;
                Point dxyz_deta;

                Point intersection_guess;
                for(unsigned int i=0; i<n_sf; ++i)
                  {

                    intersection_guess += base_elem->node_ref(i) * FE<2,LAGRANGE>::shape(base_elem->type(),
                                                                                         base_elem->default_order(),
                                                                                         i,
                                                                                         ref_point);

                    dxyz_dxi += base_elem->node_ref(i) * FE<2,LAGRANGE>::shape_deriv(base_elem->type(),
                                                                                     base_elem->default_order(),
                                                                                     i,
                                                                                     0, // d()/dxi
                                                                                     ref_point);

                    dxyz_deta+= base_elem->node_ref(i) * FE<2,LAGRANGE>::shape_deriv(base_elem->type(),
                                                                                     base_elem->default_order(),
                                                                                     i,
                                                                                     1, // d()/deta
                                                                                     ref_point);
                  } // for i

                TypeVector<Real> F(physical_point + c_factor*(physical_point-o) - intersection_guess);

                TypeTensor<Real> J;
                J(0,0) = (physical_point-o)(0);
                J(0,1) = -dxyz_dxi(0);
                J(0,2) = -dxyz_deta(0);
                J(1,0) = (physical_point-o)(1);
                J(1,1) = -dxyz_dxi(1);
                J(1,2) = -dxyz_deta(1);
                J(2,0) = (physical_point-o)(2);
                J(2,1) = -dxyz_dxi(2);
                J(2,2) = -dxyz_deta(2);

                // delta will be the newton step
                TypeVector<Real> delta;
                J.solve(F,delta);

                // check for convergence
                Real tol = std::min( TOLERANCE*0.01, TOLERANCE*base_elem->hmax() );
                if ( delta.norm() < tol )
                  {
                    // newton solver converged, now make sure it converged to a point on the base_elem
                    if (base_elem->contains_point(intersection_guess,TOLERANCE*0.1))
                      {
                        intersection = intersection_guess;
                      }
                    //else
                    //  {
                    //     err<<" Error: inverse_map failed!"<<std::endl;
                    //  }
                    break; // break out of 'for it'
                  }
                else
                  {
                    c_factor     -= delta(0);
                    ref_point(0) -= delta(1);
                    ref_point(1) -= delta(2);
                  }

              }

          }
        break;
      }

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    }

#ifndef NDEBUG
  // In debug mode, the validity of base_elems id is checked.
  // Thus, we should set it to a valid  one:
  base_elem->set_id(inf_elem->id());
#endif

  // 3.)
  // Now we have the intersection-point (projection of physical point onto base-element).
  // Lets compute its internal coordinates (being p(0) and p(1)):
  p= FEMap::inverse_map(dim-1, base_elem.get(), intersection,
                        tolerance, secure);

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
      p(dim-1)=1;
      return p;
    }

  // when we are somewhere in this element:
  Real v = 0;

  // For now we're sticking with T_map == CARTESIAN
  // if (T_map == CARTESIAN)
  v = 1.-2.*a_dist/fp_o_dist;
  // else
  //   libmesh_not_implemented();


  p(dim-1)=v;
#ifdef DEBUG
  // first check whether we are in the reference-element:
  if ((-1.-1.e-5 < v && v < 1.) || secure)
    {
      const Point check = map (dim, inf_elem, p);
      const Point diff  = physical_point - check;

      if (diff.norm() > tolerance)
        libmesh_warning("WARNING:  diff is " << diff.norm());
    }
#endif

  return p;
}



void InfFEMap::inverse_map (const unsigned int dim,
                            const Elem * elem,
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
      inverse_map (dim, elem, physical_points[p], tolerance, secure);
}


} // namespace libMesh


#endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
