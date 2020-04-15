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

#include "libmesh/inf_fe.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/int_range.h"
#include "libmesh/auto_ptr.h"

namespace libMesh
{



// Constructor
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
InfFE<Dim,T_radial,T_map>::InfFE (const FEType & fet) :
  FEBase       (Dim, fet),

  _n_total_approx_sf (0),
  _n_total_qp        (0),

  // initialize the current_fe_type to all the same
  // values as \p fet (since the FE families and coordinate
  // map type should not change), but use an invalid order
  // for the radial part (since this is the only order
  // that may change!).
  // the data structures like \p phi etc are not initialized
  // through the constructor, but through reinit()
  current_fe_type (FEType(fet.order,
                          fet.family,
                          INVALID_ORDER,
                          fet.radial_family,
                          fet.inf_map))

{
  // Sanity checks
  libmesh_assert_equal_to (T_radial, fe_type.radial_family);
  libmesh_assert_equal_to (T_map, fe_type.inf_map);

  // build the base_fe object
  if (Dim != 1)
    base_fe = FEBase::build(Dim-1, fet);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::attach_quadrature_rule (QBase * q)
{
  libmesh_assert(q);
  libmesh_assert(base_fe);

  const Order base_int_order   = q->get_order();
  const Order radial_int_order = static_cast<Order>(2 * (static_cast<unsigned int>(fe_type.radial_order.get_order()) + 1) +2);
  const unsigned int qrule_dim = q->get_dim();

  if (Dim != 1)
    {
      // build a Dim-1 quadrature rule of the type that we received
      base_qrule = QBase::build(q->type(), qrule_dim-1, base_int_order);
      base_fe->attach_quadrature_rule(base_qrule.get());
    }

  // in radial direction, always use Gauss quadrature
  radial_qrule = libmesh_make_unique<QGauss>(1, radial_int_order);

  // Maybe helpful to store the QBase *
  // with which we initialized our own quadrature rules.
  // Used e.g. in \p InfFE::reinit(elem,side)
  qrule = q;
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
void InfFE<Dim,T_radial,T_base>::update_base_elem (const Elem * inf_elem)
{
  base_elem.reset(InfFEBase::build_elem(inf_elem));
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::reinit(const Elem * inf_elem,
                                       const std::vector<Point> * const pts,
                                       const std::vector<Real> * const weights)
{
  libmesh_assert(base_fe.get());
  libmesh_assert(inf_elem);

  // I don't understand infinite elements well enough to risk
  // calculating too little.  :-(  RHS
  this->calculate_phi = this->calculate_dphi = this->calculate_dphiref = true;
  this->get_xyz();
  this->determine_calculations();
  base_fe->calculate_phi = base_fe->calculate_dphi = base_fe->calculate_dphiref = true;
  base_fe->get_xyz();
  base_fe->determine_calculations();

  if (pts == nullptr)
    {
      libmesh_assert(base_fe->qrule);
      libmesh_assert_equal_to (base_fe->qrule, base_qrule.get());
      libmesh_assert(radial_qrule.get());

      bool init_shape_functions_required = false;

      // init the radial data fields only when the radial order changes
      if (current_fe_type.radial_order != fe_type.radial_order)
        {
          current_fe_type.radial_order = fe_type.radial_order;

          // Watch out: this call to QBase->init() only works for
          // current_fe_type = const!   To allow variable Order,
          // the init() of QBase has to be modified...
          radial_qrule->init(EDGE2);

          // initialize the radial shape functions
          this->init_radial_shape_functions(inf_elem);

          init_shape_functions_required=true;
        }

      bool update_base_elem_required=true;

      // update the type in accordance to the current cell
      // and reinit if the cell type has changed or (as in
      // the case of the hierarchics) the shape functions
      // depend on the particular element and need a reinit
      if ((Dim != 1) &&
          ((this->get_type() != inf_elem->type())  ||
           (base_fe->shapes_need_reinit())))
        {
          // store the new element type, update base_elem
          // here.  Through \p update_base_elem_required,
          // remember whether it has to be updated (see below).
          elem_type = inf_elem->type();
          this->update_base_elem(inf_elem);

          update_base_elem_required=false;

          // initialize the base quadrature rule for the new element
          base_qrule->init(base_elem->type());

          init_shape_functions_required=true;

        }

      // computing the reference-to-physical map and coordinates works
      // only, if we have the current base_elem stored.
      // This happens when fe_type is const,
      // the inf_elem->type remains the same.  Then we have to
      // update the base elem _here_.
      if (update_base_elem_required)
        this->update_base_elem(inf_elem);

      // initialize the shape functions in the base
      base_fe->init_base_shape_functions(base_fe->qrule->get_points(),
                                         base_elem.get());

      // compute the shape functions and map functions of base_fe
      // before using them later in compute_shape_functions.
      base_fe->_fe_map->compute_map (base_fe->dim, base_fe->qrule->get_weights(),
                                     base_elem.get(), base_fe->calculate_d2phi);
      base_fe->compute_shape_functions(base_elem.get(), base_fe->qrule->get_points());

      // when either the radial or base part change,
      // we have to init the whole fields
      if (init_shape_functions_required)
        this->init_shape_functions (radial_qrule->get_points(),
                                    base_fe->qrule->get_points(),
                                    inf_elem);

      // Compute the shape functions and the derivatives
      // at all quadrature points.
      this->compute_shape_functions (inf_elem,
                                     base_fe->qrule->get_points(),
                                     radial_qrule->get_points()
                                     /* weights are computed insid the function*/
                                     );
    }

  else // if pts != nullptr
    {
      // update the elem_type
      elem_type = inf_elem->type();

      // We'll assume that pts is a tensor product mesh of points.
      // That will handle the pts.size()==1 case that we care about
      // right now, and it will generalize a bit, and it won't break
      // the assumptions elsewhere in InfFE.
      std::vector<Point> radial_pts;
      if (pts->size() > 0)
        {
          Real radius = (*pts)[0](Dim-1);
          radial_pts.push_back(radius);
          unsigned int n_radial_pts=1;
          unsigned int n_angular_pts=1;
          for (auto p : IntRange<std::size_t>(1, pts->size()))
            {
              radius = (*pts)[p](Dim-1);
              // check for repetition of radius: The max. allowed distance is somewhat arbitrary
              // but the given value should not produce false positives...
              if (std::abs(radial_pts[p-n_radial_pts*n_angular_pts](0) - radius)< 1e-4)
                {
                  // should the next one be in the next segment?
                  if (p+1 == n_radial_pts*(n_angular_pts+1))
                    ++n_angular_pts;
                }
              // if none yet repeated:
              else if (n_angular_pts == 1)
                {
                  radial_pts.push_back(radius);
                  ++n_radial_pts;
                }
              // if there was repetition but this does not, the assumed
              // format does not work:
              else
                {
                  libmesh_error_msg("We assumed that the "<<pts->size()
                                  <<" points are of tensor-product type with "
                                  <<n_radial_pts<<" radial points and "
                                  <<n_angular_pts<< " angular points."<<std::endl
                                  <<"But apparently point "<<p
                                  <<" does not fit that scheme: Its radius is "
                                  <<radius <<"but should have "
                                  <<radial_pts[p-n_radial_pts*n_angular_pts]<<".");
                }
            }
        }
      else
        {
          // I don't see any reason to call this function with no points.
          libmesh_error_msg("Calling reinit() with an empty point list is prohibited.\n");
        }

      const std::size_t radial_pts_size = radial_pts.size();
      const std::size_t base_pts_size = pts->size() / radial_pts_size;
      // If we're a tensor product we should have no remainder
      libmesh_assert_equal_to
        (base_pts_size * radial_pts_size, pts->size());


      std::vector<Point> base_pts;
      base_pts.reserve(base_pts_size);
      for (std::size_t p=0, ps=pts->size(); p != ps; p += radial_pts_size)
        {
          Point pt = (*pts)[p];
          pt(Dim-1) = 0;
          base_pts.push_back(pt);
        }

      // init radial shapes
      this->init_radial_shape_functions(inf_elem, &radial_pts);

      // update the base
      this->update_base_elem(inf_elem);

      //FIXME: why not use :
      //this->update_base_elem(inf_elem);
      // -> it looks better to me though...

      // the finite element on the ifem base
      base_fe = FEBase::build(Dim-1, this->fe_type);

      base_fe->calculate_phi = base_fe->calculate_dphi = base_fe->calculate_dphiref = true;
      base_fe->get_xyz();
      base_fe->determine_calculations();

      //TODO: replace the following by base_fe->reinit()!?
      // init base shapes
      base_fe->init_base_shape_functions(base_pts,
                                         base_elem.get());

      // compute the shape functions and map functions of base_fe
      // before using them later in compute_shape_functions.

      if (weights)
        {
          base_fe->_fe_map->compute_map (base_fe->dim, *weights,
                                         base_elem.get(), base_fe->calculate_d2phi);
        }
      else
        {
          std::vector<Real> dummy_weights (base_pts.size(), 1.);
          base_fe->_fe_map->compute_map (base_fe->dim, dummy_weights,
                                         base_elem.get(), base_fe->calculate_d2phi);
        }

      base_fe->compute_shape_functions(base_elem.get(), base_pts);

      this->init_shape_functions (radial_pts, base_pts, inf_elem);

      // finally compute the ifem shapes
      if (weights != nullptr)
        {
          this->compute_shape_functions (inf_elem,base_pts,radial_pts);
        }
      else
        {
          this->compute_shape_functions (inf_elem,base_pts,radial_pts);
        }

    }

}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void
InfFE<Dim, T_radial, T_map>::
init_radial_shape_functions(const Elem * libmesh_dbg_var(inf_elem),
                            const std::vector<Point> * radial_pts)
{
  libmesh_assert(radial_qrule.get() || radial_pts);
  libmesh_assert(inf_elem);

  // Start logging the radial shape function initialization
  LOG_SCOPE("init_radial_shape_functions()", "InfFE");

  // initialize most of the things related to physical approximation
  const Order radial_approx_order = fe_type.radial_order;
  const unsigned int n_radial_approx_shape_functions =
    InfFERadial::n_dofs(radial_approx_order);

  const std::size_t n_radial_qp =
    radial_pts ? radial_pts->size() : radial_qrule->n_points();
  const std::vector<Point> & radial_qp =
    radial_pts ? *radial_pts : radial_qrule->get_points();

  // the radial polynomials (eval)
  mode.resize      (n_radial_approx_shape_functions);
  dmodedv.resize   (n_radial_approx_shape_functions);

  // the (1-v)/2 weight
  som.resize       (n_radial_qp);
  dsomdv.resize    (n_radial_qp);


  for (unsigned int i=0; i<n_radial_approx_shape_functions; ++i)
    {
      mode[i].resize    (n_radial_qp);
      dmodedv[i].resize (n_radial_qp);
    }


  // compute scalar values at radial quadrature points
  for (std::size_t p=0; p<n_radial_qp; ++p)
    {
      som[p] = InfFERadial::decay (Dim, radial_qp[p](0));
      dsomdv[p] = InfFERadial::decay_deriv (radial_qp[p](0));
    }


  // evaluate the mode shapes in radial direction at radial quadrature points
  for (unsigned int i=0; i<n_radial_approx_shape_functions; ++i)
    for (std::size_t p=0; p<n_radial_qp; ++p)
      {
        mode[i][p] = InfFE<Dim,T_radial,T_map>::eval (radial_qp[p](0), radial_approx_order, i);
        dmodedv[i][p] = InfFE<Dim,T_radial,T_map>::eval_deriv (radial_qp[p](0), radial_approx_order, i);
      }

}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::init_shape_functions(const std::vector<Point> & radial_qp,
                                                     const std::vector<Point> & base_qp,
                                                     const Elem * inf_elem)
{
  libmesh_assert(inf_elem);

  // Start logging the radial shape function initialization
  LOG_SCOPE("init_shape_functions()", "InfFE");

  // fast access to some const ints for the radial data
  const unsigned int n_radial_approx_sf  = cast_int<unsigned int>(mode.size());
  const unsigned int n_radial_qp         = cast_int<unsigned int>(som.size());


  // initialize most of the quantities related to mapping

  // The element type and order to use in the base map
  //const Order base_mapping_order = base_elem->default_order();

  // the number of base shape functions used to construct the map
  // (Lagrange shape functions are used for mapping in the base)
  //unsigned int n_base_mapping_shape_functions =
  //  InfFEBase::n_base_mapping_sf(*base_elem,
  //                               base_mapping_order);

  // initialize most of the things related to physical approximation
  unsigned int n_base_approx_shape_functions;
  if (Dim > 1)
    n_base_approx_shape_functions = base_fe->n_shape_functions();
  else
    n_base_approx_shape_functions = 1;


  // update class member field
  _n_total_approx_sf =
     n_radial_approx_sf * n_base_approx_shape_functions;


  // The number of the base quadrature points.
  const unsigned int n_base_qp = cast_int<unsigned int>(base_qp.size());

  // The total number of quadrature points.
  _n_total_qp = n_radial_qp * n_base_qp;


  // initialize the node and shape numbering maps
  {
    // similar for the shapes: the i-th entry stores
    // the associated base/radial shape number
    _radial_shape_index.resize(_n_total_approx_sf);
    _base_shape_index.resize(_n_total_approx_sf);

    const ElemType inf_elem_type = inf_elem->type();

    // fill the shape index map
    for (unsigned int n=0; n<_n_total_approx_sf; ++n)
      {
        compute_shape_indices (this->fe_type,
                               inf_elem_type,
                               n,
                               _base_shape_index[n],
                               _radial_shape_index[n]);
        libmesh_assert_less (_base_shape_index[n], n_base_approx_shape_functions);
        libmesh_assert_less (_radial_shape_index[n], n_radial_approx_sf);
      }
  }

  // resize the base data fields
  //dist.resize(n_base_mapping_shape_functions);

  // resize the total data fields

  // the phase term varies with xi, eta and zeta(v): store it for _all_ qp
  //
  // when computing the phase, we need the base approximations
  // therefore, initialize the phase here, but evaluate it
  // in compute_shape_functions().
  //
  // the weight, though, is only needed at the radial quadrature points, n_radial_qp.
  // but for a uniform interface to the protected data fields
  // the weight data field (which are accessible from the outside) are expanded to _n_total_qp.
  weight.resize      (_n_total_qp);
  weightxr_sq.resize (_n_total_qp);
  dweightdv.resize   (n_radial_qp);
  dweight.resize     (_n_total_qp);
  dweightxr_sq.resize(_n_total_qp);

  dphase.resize      (_n_total_qp);

  // this vector contains the integration weights for the combined quadrature rule
  // if no quadrature rules are given, use only ones.
  _total_qrule_weights.resize(_n_total_qp,1);

  // InfFE's data fields phi, dphi, dphidx, phi_map etc hold the _total_
  // shape and mapping functions, respectively
  {
    JxWxdecay.resize(_n_total_qp);
    xyz.resize(_n_total_qp);
    dxidx_map.resize(_n_total_qp);
    dxidy_map.resize(_n_total_qp);
    dxidz_map.resize(_n_total_qp);
    detadx_map.resize(_n_total_qp);
    detady_map.resize(_n_total_qp);
    detadz_map.resize(_n_total_qp);
    dzetadx_map.resize(_n_total_qp);
    dzetady_map.resize(_n_total_qp);
    dzetadz_map.resize(_n_total_qp);
    dxidx_map_scaled.resize(_n_total_qp);
    dxidy_map_scaled.resize(_n_total_qp);
    dxidz_map_scaled.resize(_n_total_qp);
    detadx_map_scaled.resize(_n_total_qp);
    detady_map_scaled.resize(_n_total_qp);
    detadz_map_scaled.resize(_n_total_qp);
    dzetadx_map_scaled.resize(_n_total_qp);
    dzetady_map_scaled.resize(_n_total_qp);
    dzetadz_map_scaled.resize(_n_total_qp);

    phi.resize     (_n_total_approx_sf);
    dphi.resize    (_n_total_approx_sf);
    dphidx.resize  (_n_total_approx_sf);
    dphidy.resize  (_n_total_approx_sf);
    dphidz.resize  (_n_total_approx_sf);
    dphidxi.resize (_n_total_approx_sf);

    phixr.resize    (_n_total_approx_sf);
    dphixr.resize   (_n_total_approx_sf);
    dphixr_sq.resize(_n_total_approx_sf);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    libmesh_warning("Warning: Second derivatives for Infinite elements"
                    << " are not yet implemented!"
                    << std::endl);

    d2phi.resize     (_n_total_approx_sf);
    d2phidx2.resize  (_n_total_approx_sf);
    d2phidxdy.resize (_n_total_approx_sf);
    d2phidxdz.resize (_n_total_approx_sf);
    d2phidy2.resize  (_n_total_approx_sf);
    d2phidydz.resize (_n_total_approx_sf);
    d2phidz2.resize  (_n_total_approx_sf);
    d2phidxi2.resize (_n_total_approx_sf);

    if (Dim > 1)
      {
        d2phidxideta.resize(_n_total_approx_sf);
        d2phideta2.resize(_n_total_approx_sf);
      }

    if (Dim > 2)
      {
        d2phidetadzeta.resize(_n_total_approx_sf);
        d2phidxidzeta.resize(_n_total_approx_sf);
        d2phidzeta2.resize(_n_total_approx_sf);
      }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

    if (Dim > 1)
      dphideta.resize(_n_total_approx_sf);

    if (Dim == 3)
      dphidzeta.resize(_n_total_approx_sf);

  }

  // collect all the for loops, where inner vectors are
  // resized to the appropriate number of quadrature points
  {
    for (unsigned int i=0; i<_n_total_approx_sf; ++i)
      {
        phi[i].resize         (_n_total_qp);
        dphi[i].resize        (_n_total_qp);
        dphidx[i].resize      (_n_total_qp);
        dphidy[i].resize      (_n_total_qp);
        dphidz[i].resize      (_n_total_qp);
        dphidxi[i].resize     (_n_total_qp);

        phixr[i].resize (_n_total_qp);
        dphixr[i].resize(_n_total_qp);
        dphixr_sq[i].resize(_n_total_qp);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        d2phi[i].resize       (_n_total_qp);
        d2phidx2[i].resize    (_n_total_qp);
        d2phidxdy[i].resize   (_n_total_qp);
        d2phidxdz[i].resize   (_n_total_qp);
        d2phidy2[i].resize    (_n_total_qp);
        d2phidydz[i].resize   (_n_total_qp);
        d2phidy2[i].resize    (_n_total_qp);
        d2phidxi2[i].resize   (_n_total_qp);

        if (Dim > 1)
          {
            d2phidxideta[i].resize   (_n_total_qp);
            d2phideta2[i].resize     (_n_total_qp);
          }
        if (Dim > 2)
          {
            d2phidxidzeta[i].resize  (_n_total_qp);
            d2phidetadzeta[i].resize (_n_total_qp);
            d2phidzeta2[i].resize    (_n_total_qp);
          }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        if (Dim > 1)
          dphideta[i].resize  (_n_total_qp);

        if (Dim == 3)
          dphidzeta[i].resize (_n_total_qp);

      }

  }

  {
    // (a) compute scalar values at _all_ quadrature points  -- for uniform
    //     access from the outside to these fields
    // (b) form a std::vector<Real> which contains the appropriate weights
    //     of the combined quadrature rule!
    libmesh_assert_equal_to (radial_qp.size(), n_radial_qp);

    if (radial_qrule && base_qrule)
      {
        const std::vector<Real> & radial_qw = radial_qrule->get_weights();
        const std::vector<Real> & base_qw = base_qrule->get_weights();
        libmesh_assert_equal_to (radial_qw.size(), n_radial_qp);
        libmesh_assert_equal_to (base_qw.size(), n_base_qp);

        for (unsigned int rp=0; rp<n_radial_qp; ++rp)
          {
            for (unsigned int bp=0; bp<n_base_qp; ++bp)
              {
                weight[bp + rp*n_base_qp] = InfFERadial::D(radial_qp[rp](0));
                weightxr_sq[bp + rp*n_base_qp] = InfFERadial::Dxr_sq(radial_qp[rp](0));
                _total_qrule_weights[bp + rp*n_base_qp] = radial_qw[rp] * base_qw[bp];
              }
            dweightdv[rp] = InfFERadial::D_deriv(radial_qp[rp](0));
          }
      }
    else
      {
        for (unsigned int rp=0; rp<n_radial_qp; ++rp)
          {
            for (unsigned int bp=0; bp<n_base_qp; ++bp)
              {
                weight[bp + rp*n_base_qp] = InfFERadial::D(radial_qp[rp](0));
                weightxr_sq[bp + rp*n_base_qp] = InfFERadial::Dxr_sq(radial_qp[rp](0));
              }
            dweightdv[rp] = InfFERadial::D_deriv(radial_qp[rp](0));
          }
      }
  }
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_shape_functions(const Elem * inf_elem,
                                                        const std::vector<Point> & base_qp,
                                                        const std::vector<Point> & radial_qp
                                                        )
{
  libmesh_assert(inf_elem);
  // at least check whether the base element type is correct.
  // otherwise this version of computing dist would give problems
  libmesh_assert_equal_to (base_elem->type(),
                           InfFEBase::get_elem_type(inf_elem->type()));

  // Start logging the overall computation of shape functions
  LOG_SCOPE("compute_shape_functions()", "InfFE");

  // these vectors are needed later; initialize here already to have access to
  // n_base_qp etc.
  const std::vector<std::vector<Real>> & S  = base_fe->phi;
  const std::vector<std::vector<Real>> & S_map  = (base_fe->get_fe_map()).get_phi_map();

  const unsigned int n_radial_qp = cast_int<unsigned int>(som.size());
  const unsigned int n_base_qp = cast_int<unsigned int>(S_map[0].size());

  libmesh_assert_equal_to (_n_total_qp, n_radial_qp*n_base_qp);
  libmesh_assert_equal_to (_n_total_qp, _total_qrule_weights.size());
#ifdef DEBUG
  if (radial_qrule)
    libmesh_assert_equal_to(n_radial_qp, radial_qrule->n_points());
  if (base_qrule)
    libmesh_assert_equal_to(n_base_qp, base_qrule->n_points());
  libmesh_assert_equal_to(_n_total_qp % n_radial_qp, 0); // "Error in the structure of quadrature pointes!");
#endif


  _n_total_approx_sf = InfFERadial::n_dofs(fe_type.radial_order) *
                                base_fe->n_shape_functions();



  const Point origin = inf_elem->origin();

  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Compute the value of the derivative shape function i at quadrature point p
  switch (dim)
    {
    case 1:
    case 2:
      {
        libmesh_not_implemented();
        break;
      }

    case 3:
      {
        // fast access to the approximation and mapping shapes of base_fe
        const std::vector<std::vector<Real>> & Ss = base_fe->dphidxi;
        const std::vector<std::vector<Real>> & St = base_fe->dphideta;

        const std::vector<Real> & base_dxidx = base_fe->get_dxidx();
        const std::vector<Real> & base_dxidy = base_fe->get_dxidy();
        const std::vector<Real> & base_dxidz = base_fe->get_dxidz();
        const std::vector<Real> & base_detadx = base_fe->get_detadx();
        const std::vector<Real> & base_detady = base_fe->get_detady();
        const std::vector<Real> & base_detadz = base_fe->get_detadz();

        const std::vector<Point> & base_xyz = base_fe->get_xyz();

        libmesh_assert_equal_to (phi.size(), _n_total_approx_sf);
        libmesh_assert_equal_to (dphidxi.size(), _n_total_approx_sf);
        libmesh_assert_equal_to (dphideta.size(), _n_total_approx_sf);
        libmesh_assert_equal_to (dphidzeta.size(), _n_total_approx_sf);

        unsigned int tp=0; // total qp
        for (unsigned int rp=0; rp<n_radial_qp; ++rp)  // over radial qps
          for (unsigned int bp=0; bp<n_base_qp; ++bp)  // over base qps

            { // First compute the map from base element quantities to physical space:
              xyz[tp] = map(inf_elem,Point(base_qp[bp](0),base_qp[bp](1),radial_qp[rp](0)));

              // Compute the shape function derivatives wrt x,y at the
              // quadrature points
              // dxi,eta dx,y,z is dictated by the base elements mapping, but decreases with a/r.
              dxidx_map[tp] = base_dxidx[bp]*som[rp];
              dxidy_map[tp] = base_dxidy[bp]*som[rp];
              dxidz_map[tp] = base_dxidz[bp]*som[rp];

              detadx_map[tp] =  base_detadx[bp]*som[rp];
              detady_map[tp] =  base_detady[bp]*som[rp];
              detadz_map[tp] =  base_detadz[bp]*som[rp];

              // som = a/r
              dxidx_map_scaled[tp] = base_dxidx[bp];
              dxidy_map_scaled[tp] = base_dxidy[bp];
              dxidz_map_scaled[tp] = base_dxidz[bp];

              detadx_map_scaled[tp] =  base_detadx[bp];
              detady_map_scaled[tp] =  base_detady[bp];
              detadz_map_scaled[tp] =  base_detadz[bp];

              const Point r(xyz[tp]-origin);
              const Real a ((base_xyz[bp]-origin).norm());

              // check that 'som' == a/r.
              libmesh_assert_less(std::abs(som[rp] -a/r.norm()) , 1e-7);

              const Point unit_r(r/r.norm());
              const Real dzetadr_map= 2.*a/r.norm_sq();

              // When rescaling, we use r as being normalized to the base_elems coordinate,
              // thus always start with r=1 at the FE/InfFE-boundary.
              const Real dzetadr_mapxr_sq= 2./a;

              Point e_eta(detadx_map[tp],
                          detady_map[tp],
                          detadz_map[tp]);
              Point e_xi(dxidx_map[tp],
                         dxidy_map[tp],
                         dxidz_map[tp]);

              RealGradient normal=e_eta.cross(e_xi);
              normal/=normal.norm();

              // grad_a divided by a/r
              RealGradient grad_a_scaled=unit_r - normal/(normal*unit_r);

              if (!inf_elem->side_ptr(0)->has_affine_map())
                {
                  /**
                   * The full form for 'a' is
                   * a =  (r0*normal)/(normal*unit_r);
                   * where r0 is some point on the base plane(!)
                   * when the base element is not a plane, r0 and normal are functions of space.
                   * Here, some approximation is used:
                   */

                   const unsigned int n_sf = base_fe->n_shape_functions();
                   for (unsigned int i=0; i< n_sf; ++i)
                   {
                       grad_a_scaled += (FE<2,LAGRANGE>::shape_deriv(base_elem.get(),
                                                    base_elem->default_order(),
                                                    i, 0, base_xyz[bp]) * e_xi
                                        +FE<2,LAGRANGE>::shape_deriv(base_elem.get(),
                                                    base_elem->default_order(),
                                                    i, 1, base_xyz[bp]) * e_eta )
                                         *(normal*base_elem->node_ref(i))/(normal*unit_r*a);
                   }

                }

              // dzetadx = dzetadr - 2/r * grad_a
              //         = dzetadr - 2*a/r^2 * grad_a_scaled
              dzetadx_map[tp] =dzetadr_map*(unit_r(0) - grad_a_scaled(0));
              dzetady_map[tp] =dzetadr_map*(unit_r(1) - grad_a_scaled(1));
              dzetadz_map[tp] =dzetadr_map*(unit_r(2) - grad_a_scaled(2));
              dzetadx_map_scaled[tp] =dzetadr_mapxr_sq*(unit_r(0) - grad_a_scaled(0));
              dzetady_map_scaled[tp] =dzetadr_mapxr_sq*(unit_r(1) - grad_a_scaled(1));
              dzetadz_map_scaled[tp] =dzetadr_mapxr_sq*(unit_r(2) - grad_a_scaled(2));

              Real inv_jacxR_pow4 = (dxidx_map_scaled[tp]  *(    detady_map_scaled[tp]*dzetadz_map_scaled[tp]- dzetady_map_scaled[tp]*detadz_map_scaled[tp]) +
                                     detadx_map_scaled[tp] *(dzetady_map_scaled[tp]*     dxidz_map_scaled[tp]-  dxidy_map_scaled[tp]*dzetadz_map_scaled[tp]) +
                                     dzetadx_map_scaled[tp]*(     dxidy_map_scaled[tp]* detadz_map_scaled[tp]- detady_map_scaled[tp]*  dxidz_map_scaled[tp]));

              if (inv_jacxR_pow4 <= 1e-7)
                {
                  libmesh_error_msg("ERROR: negative inverse Jacobian " \
                                    << inv_jacxR_pow4 \
                                    << " at point " \
                                    << xyz[tp] \
                                    << " in element " \
                                    << inf_elem->id());
                }

              JxWxdecay[tp] = _total_qrule_weights[tp]/inv_jacxR_pow4;

              // phase term mu(r)=i*k*(r-a).
              // skip i*k: it is added separately during matrix assembly.
              dphase[tp] = (unit_r - grad_a_scaled*a/r.norm());

              // This is the derivative of the phase term of this infinite element
              dweight[tp](0) = dweightdv[rp] * dzetadx_map[tp];
              dweight[tp](1) = dweightdv[rp] * dzetady_map[tp];
              dweight[tp](2) = dweightdv[rp] * dzetadz_map[tp];

              dweightxr_sq[tp](0) = dweightdv[rp] * dzetadx_map_scaled[tp];
              dweightxr_sq[tp](1) = dweightdv[rp] * dzetady_map_scaled[tp];
              dweightxr_sq[tp](2) = dweightdv[rp] * dzetadz_map_scaled[tp];

              // compute the shape-functions and derivative quantities:
              for (auto i : index_range(phi))
                {

                  // let the index vectors take care of selecting the appropriate base/radial shape
                  unsigned int bi = _base_shape_index  [i];
                  unsigned int ri = _radial_shape_index[i];
                  phi      [i][tp] = S [bi][bp] * mode[ri][rp] * som[rp];
                  dphidxi  [i][tp] = Ss[bi][bp] * mode[ri][rp] * som[rp];
                  dphideta [i][tp] = St[bi][bp] * mode[ri][rp] * som[rp];
                  dphidzeta[i][tp] = S [bi][bp]
                    * (dmodedv[ri][rp] * som[rp] + mode[ri][rp] * dsomdv[rp]);

                  // temporary quantities used below. Not accessible to user-code.
                  const Real dphidxixr = Ss[bi][bp] * mode[ri][rp];
                  const Real dphidetaxr= St[bi][bp] * mode[ri][rp];

                  phixr[i][tp] = S [bi][bp] * mode[ri][rp];

                  // dphi/dx    = (dphi/dxi)*(dxi/dx) + (dphi/deta)*(deta/dx) + (dphi/dzeta)*(dzeta/dx);
                  dphi[i][tp](0) =
                    dphidx[i][tp] = (dphidxi[i][tp]*dxidx_map[tp] +
                                     dphideta[i][tp]*detadx_map[tp] +
                                     dphidzeta[i][tp]*dzetadx_map[tp]);

                  // dphi/dy    = (dphi/dxi)*(dxi/dy) + (dphi/deta)*(deta/dy) + (dphi/dzeta)*(dzeta/dy);
                  dphi[i][tp](1) =
                    dphidy[i][tp] = (dphidxi[i][tp]*dxidy_map[tp] +
                                     dphideta[i][tp]*detady_map[tp] +
                                     dphidzeta[i][tp]*dzetady_map[tp]);

                  // dphi/dz    = (dphi/dxi)*(dxi/dz) + (dphi/deta)*(deta/dz) + (dphi/dzeta)*(dzeta/dz);
                  dphi[i][tp](2) =
                    dphidz[i][tp] = (dphidxi[i][tp]*dxidz_map[tp] +
                                     dphideta[i][tp]*detadz_map[tp] +
                                     dphidzeta[i][tp]*dzetadz_map[tp]);

                  dphixr[i][tp](0)= (dphidxi[i][tp]*dxidx_map_scaled[tp] +
                                     dphideta[i][tp]*detadx_map_scaled[tp] +
                                     dphidzeta[i][tp]*dzetadx_map_scaled[tp]*som[rp]);

                  dphixr[i][tp](1) = (dphidxi[i][tp]*dxidy_map_scaled[tp] +
                                      dphideta[i][tp]*detady_map_scaled[tp] +
                                      dphidzeta[i][tp]*dzetady_map_scaled[tp]*som[rp]);

                  dphixr[i][tp](2) = (dphidxi[i][tp]*dxidz_map_scaled[tp] +
                                      dphideta[i][tp]*detadz_map_scaled[tp] +
                                      dphidzeta[i][tp]*dzetadz_map_scaled[tp]*som[rp]);

                  dphixr_sq[i][tp](0)= (dphidxixr*dxidx_map_scaled[tp] +
                                        dphidetaxr*detadx_map_scaled[tp] +
                                        dphidzeta[i][tp]*dzetadx_map_scaled[tp]);

                  dphixr_sq[i][tp](1) = (dphidxixr*dxidy_map_scaled[tp] +
                                         dphidetaxr*detady_map_scaled[tp] +
                                         dphidzeta[i][tp]*dzetady_map_scaled[tp]);

                  dphixr_sq[i][tp](2) = (dphidxixr*dxidz_map_scaled[tp] +
                                         dphidetaxr*detadz_map_scaled[tp] +
                                         dphidzeta[i][tp]*dzetadz_map_scaled[tp]);

                }

              tp++;
            }

        break;
      }

    default:
      libmesh_error_msg("Unsupported dim = " << dim);
    }
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
bool InfFE<Dim,T_radial,T_map>::shapes_need_reinit() const
{
  return false;
}

} // namespace libMesh


// Explicit instantiations
#include "libmesh/inf_fe_instantiate_1D.h"
#include "libmesh/inf_fe_instantiate_2D.h"
#include "libmesh/inf_fe_instantiate_3D.h"



#endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
