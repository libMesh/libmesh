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



// C++ includes


// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#include "libmesh/inf_fe.h"
#include "libmesh/inf_fe_macro.h"
#include "libmesh/quadrature.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/elem.h"

namespace libMesh
{

// Method for 2D, 3D -- see inf_fe_1D.C for a 1D version of this
template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
void InfFE<Dim,T_radial,T_base>::reinit(const Elem * inf_elem,
                                        const unsigned int s,
                                        const Real /*tolerance*/,
                                        const std::vector<Point> * const pts,
                                        const std::vector<Real> * const weights)
{
  if (weights != nullptr)
    libmesh_not_implemented_msg("ERROR: User-specified weights for infinite elements are not implemented!");

  if (pts != nullptr)
    libmesh_not_implemented_msg("ERROR: User-specified points for infinite elements are not implemented!");

  // We don't do this for 1D elements!
  libmesh_assert_not_equal_to (Dim, 1);

  libmesh_assert(inf_elem);
  libmesh_assert(qrule);

  // Build the side of interest
  const std::unique_ptr<const Elem> side(inf_elem->build_side_ptr(s));

  // set the element type
  elem_type = inf_elem->type();

  // eventually initialize radial quadrature rule
  bool radial_qrule_initialized = false;

  // if we are working on the base-side, the radial function is constant.
  // With this, we ensure that at least for base elements we reinitialize all quantities
  // when we enter for the first time.
  if (s == 0)
    current_fe_type.radial_order = 0;
  else
    /**
     * After the recent larger changes, this case was not tested.
     * It might work, but maybe it gives wrong results.
     */
    libmesh_not_implemented();

  if (current_fe_type.radial_order != fe_type.radial_order)
    {
      if (s > 0)
        {
          current_fe_type.radial_order = fe_type.radial_order;
          radial_qrule->init(EDGE2, inf_elem->p_level());
        }
      else
        {
          // build a new 0-dimensional quadrature-rule:
          radial_qrule=QBase::build(QGAUSS, 0, fe_type.radial_order);
          radial_qrule->init(NODEELEM, 0);

          //the base_qrule is set up with dim-1, but apparently we need dim, so we replace it:
          base_qrule=QBase::build(qrule->type(), side->dim(), qrule->get_order());

          unsigned int side_p_level = inf_elem->p_level();
          if (inf_elem->neighbor_ptr(s) != nullptr)
            side_p_level = std::max(side_p_level, inf_elem->neighbor_ptr(s)->p_level());
          base_qrule->init(side->type(), side_p_level);
        }
      radial_qrule_initialized = true;
    }

  // Initialize the face shape functions
  if (this->get_type() != inf_elem->type() ||
      base_fe->shapes_need_reinit()        ||
      radial_qrule_initialized)
    this->init_face_shape_functions (qrule->get_points(), side.get());

  // The reinit() function computes all what we want except for
  //  - normal, tangents: They are not considered
  // This is done below:
  compute_face_functions();
}



// Method for 2D, 3D -- see inf_fe_1D.C for a 1D version of this
template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
void InfFE<Dim,T_radial,T_base>::edge_reinit(const Elem *,
                                             const unsigned int,
                                             const Real,
                                             const std::vector<Point> * const pts,
                                             const std::vector<Real> * const /*weights*/)
{
  // We don't do this for 1D elements!
  //libmesh_assert_not_equal_to (Dim, 1);
  libmesh_not_implemented_msg("ERROR: Edge conditions for infinite elements not implemented!");

  if (pts != nullptr)
    libmesh_not_implemented_msg("ERROR: User-specified points for infinite elements not implemented!");
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
void InfFE<Dim,T_radial,T_base>::init_face_shape_functions(const std::vector<Point> &,
                                                           const Elem * inf_side)
{
  libmesh_assert(inf_side);

  // Currently, this makes only sense in 3-D!
  libmesh_assert_equal_to (Dim, 3);

  // Initialize the radial shape functions (in particular som)
  this->init_radial_shape_functions(inf_side);

  // Initialize the base shape functions
  if (inf_side->infinite())
    this->update_base_elem(inf_side);
  else
    // in this case, I need the 2D base
    this->update_base_elem(inf_side->interior_parent());

  // Initialize the base quadrature rule
  base_qrule->init(base_elem->type(), inf_side->p_level());

  // base_fe still corresponds to the (dim-1)-dimensional base of the InfFE object,
  // so update the fe_base.
  if (inf_side->infinite())
    {
      base_fe = FEBase::build(1, this->fe_type);
      base_fe->attach_quadrature_rule(base_qrule.get());
    }
  else
    {
      base_fe = FEBase::build(Dim-1, this->fe_type);
      base_fe->attach_quadrature_rule(base_qrule.get());
    }

  //before initializing, we should say what to compute:
  base_fe->_fe_map->get_xyz();
  base_fe->_fe_map->get_JxW();

  // initialize the shape functions on the base
  base_fe->init_base_shape_functions(base_fe->qrule->get_points(),
                                     base_elem.get());

  // the number of quadrature points
  const unsigned int n_radial_qp =
    cast_int<unsigned int>(som.size());
  const unsigned int n_base_qp   = base_qrule->n_points();
  const unsigned int n_total_qp  = n_radial_qp * n_base_qp;

#ifdef DEBUG
  // when evaluating the base side, there should be only one radial point.
  if (!inf_side->infinite())
    libmesh_assert_equal_to (n_radial_qp, 1);
#endif

  // the quadrature weights
  _total_qrule_weights.resize(n_total_qp);
  std::vector<Point> qp(n_total_qp);

  // quadrature rule weights
  if (Dim < 3)
    {
      // the quadrature points must be assembled differently for lower dims.
      libmesh_not_implemented();
    }
  {
    const std::vector<Real> & radial_qw = radial_qrule->get_weights();
    const std::vector<Real> & base_qw   = base_qrule->get_weights();
    const std::vector<Point> & radial_qp = radial_qrule->get_points();
    const std::vector<Point> & base_qp   = base_qrule->get_points();

    libmesh_assert_equal_to (radial_qw.size(), n_radial_qp);
    libmesh_assert_equal_to (base_qw.size(), n_base_qp);

    for (unsigned int rp=0; rp<n_radial_qp; rp++)
      for (unsigned int bp=0; bp<n_base_qp; bp++)
        {
          _total_qrule_weights[bp + rp*n_base_qp] = radial_qw[rp] * base_qw[bp];
          // initialize the quadrature-points for the 2D side element
          // - either the base element or it has a 1D base + radial direction.
          if (inf_side->infinite())
            qp[bp + rp*n_base_qp]=Point(base_qp[bp](0),
                                        0.,
                                        radial_qp[rp](0));
          else
            qp[bp + rp*n_base_qp]=Point(base_qp[bp](0),
                                        base_qp[bp](1),
                                        -1.);
        }
  }

  this->reinit(inf_side->interior_parent(), &qp);

}

template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
void InfFE<Dim,T_radial,T_base>::compute_face_functions()
{

  const unsigned int n_qp = cast_int<unsigned int>(_total_qrule_weights.size());
  if (calculate_dxyz)
    {
      this->normals.resize(n_qp);

      if (Dim > 1)
        {
          this->tangents.resize(n_qp);
          for (unsigned int p=0; p<n_qp; ++p)
            this->tangents[p].resize(LIBMESH_DIM-1);
        }
      else
        {
          libMesh::err << "tangents have no sense in 1-dimensional elements!"<<std::endl;
          libmesh_error_msg("Exiting...");
        }
    }

  // If we have no quadrature points, there's nothing else to do
  if (!n_qp)
    return;

  switch(Dim)
    {
    case 1:
    case 2:
      {
        libmesh_not_implemented();
        break;
      }
    case 3:
      {
        for (unsigned int p=0; p<n_qp; ++p)
          {

            if (calculate_dxyz)
              {
                //
                // seeking dxyzdx, dxyzdeta means to compute
                //        / dx/dxi    dy/dxi   dz/dxi \.
                // J^-1= |                             |
                //       \ dx/deta  dy/deta   dz/deta /.
                // which is the psudo-inverse of J, i.e.
                //
                // J^-1 = (J^T J)^-1 J^T
                //
                // where J^T T is the 2x2 matrix 'g' used to compute the
                // Jacobian determinant; thus
                //
                // J^-1 = ________1________   / g22  -g21 \  / dxi/dx  dxi/dy   dxi/dz \.
                //        g11*g22 - g21*g12   \-g12  g11  /  \ deta/dx deta/dy deta/dz /.
                const std::vector<Real> & base_dxidx = base_fe->get_dxidx();
                const std::vector<Real> & base_dxidy = base_fe->get_dxidy();
                const std::vector<Real> & base_dxidz = base_fe->get_dxidz();
                const std::vector<Real> & base_detadx = base_fe->get_detadx();
                const std::vector<Real> & base_detady = base_fe->get_detady();
                const std::vector<Real> & base_detadz = base_fe->get_detadz();

                const Real g11 = (base_dxidx[p]*base_dxidx[p] +
                                  base_dxidy[p]*base_dxidy[p] +
                                  base_dxidz[p]*base_dxidz[p]);
                const Real g12 = (base_dxidx[p]*base_detadx[p] +
                                  base_dxidy[p]*base_detady[p] +
                                  base_dxidz[p]*base_detadz[p]);
                const Real g21 = g12;
                const Real g22 = (base_detadx[p]*base_detadx[p] +
                                  base_detady[p]*base_detady[p] +
                                  base_detadz[p]*base_detadz[p]);

                // det is scaled by r^6
                const Real det = (g11*g22 - g12*g21);

                // scaled by r^-3
                Point dxyzdxi_map((g22*base_dxidx[p]-g21*base_detadx[p])/det,
                                  (g22*base_dxidy[p]-g21*base_detady[p])/det,
                                  (g22*base_dxidz[p]-g21*base_detadz[p])/det);

                Point dxyzdeta_map((g11*base_detadx[p] - g12*base_dxidx[p])/det,
                                   (g11*base_detady[p] - g12*base_dxidy[p])/det,
                                   (g11*base_detadz[p] - g12*base_dxidz[p])/det);
                // scaled by r^-2

                this->tangents[p][0] = dxyzdxi_map.unit();

                this->tangents[p][1] = (dxyzdeta_map - (dxyzdeta_map*tangents[p][0])*tangents[p][0] ).unit();

                this->normals[p]     = tangents[p][0].cross(tangents[p][1]).unit();
                // recompute JxW using the 2D Jacobian:
                this->JxWxdecay[p] = _total_qrule_weights[p]/std::sqrt(det);
                this->JxW[p] = JxWxdecay[p];
              }


          }
        break;
      }
    default:
      libmesh_error_msg("Unsupported dim = " << dim);
    }

}



// Explicit instantiations - doesn't make sense in 1D, but as
// long as we only return errors, we are fine... ;-)
//#include "libmesh/inf_fe_instantiate_1D.h"
//#include "libmesh/inf_fe_instantiate_2D.h"
//#include "libmesh/inf_fe_instantiate_3D.h"
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, reinit(const Elem *, const unsigned int, const Real, const std::vector<Point> * const, const std::vector<Real> * const));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, reinit(const Elem *, const unsigned int, const Real, const std::vector<Point> * const, const std::vector<Real> * const));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, reinit(const Elem *, const unsigned int, const Real, const std::vector<Point> * const, const std::vector<Real> * const));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, edge_reinit(const Elem *, const unsigned int, const Real, const std::vector<Point> * const, const std::vector<Real> * const));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, edge_reinit(const Elem *, const unsigned int, const Real, const std::vector<Point> * const, const std::vector<Real> * const));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, edge_reinit(const Elem *, const unsigned int, const Real, const std::vector<Point> * const, const std::vector<Real> * const));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, init_face_shape_functions(const std::vector<Point> &, const Elem *));
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, init_face_shape_functions(const std::vector<Point> &, const Elem *));
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, init_face_shape_functions(const std::vector<Point> &, const Elem *));
INSTANTIATE_INF_FE_MBRF(1, CARTESIAN, void, compute_face_functions());
INSTANTIATE_INF_FE_MBRF(2, CARTESIAN, void, compute_face_functions());
INSTANTIATE_INF_FE_MBRF(3, CARTESIAN, void, compute_face_functions());

} // namespace libMesh

#endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
