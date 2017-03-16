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
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"

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
  // map type should @e not change), but use an invalid order
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
    base_fe.reset(FEBase::build(Dim-1, fet).release());
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::attach_quadrature_rule (QBase * q)
{
  libmesh_assert(q);
  libmesh_assert(base_fe.get());

  const Order base_int_order   = q->get_order();
  const Order radial_int_order = static_cast<Order>(2 * (static_cast<unsigned int>(fe_type.radial_order.get_order()) + 1) +2);
  const unsigned int qrule_dim = q->get_dim();

  if (Dim != 1)
    {
      // build a Dim-1 quadrature rule of the type that we received
      base_qrule.reset(QBase::build(q->type(), qrule_dim-1, base_int_order).release());
      base_fe->attach_quadrature_rule(base_qrule.get());
    }

  // in radial direction, always use Gauss quadrature
  radial_qrule.reset(new QGauss(1, radial_int_order));

  // currently not used. But maybe helpful to store the QBase *
  // with which we initialized our own quadrature rules
  qrule = q;
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
void InfFE<Dim,T_radial,T_base>::update_base_elem (const Elem * inf_elem)
{
  base_elem.reset(Base::build_elem(inf_elem));
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

  if (pts == libmesh_nullptr)
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

          // initialize the shape functions in the base
          base_fe->init_base_shape_functions(base_fe->qrule->get_points(),
                                             base_elem.get());

          // compute the shape functions and map functions of base_fe
          // before using them later in combine_base_radial.
          base_fe->_fe_map->compute_map (base_fe->dim, base_fe->qrule->get_weights(),
                                         base_elem.get(), base_fe->calculate_d2phi);
          base_fe->compute_shape_functions(base_elem.get(), base_fe->qrule->get_points());

          init_shape_functions_required=true;
        }


      // when either the radial or base part change,
      // we have to init the whole fields
      if (init_shape_functions_required)
        this->init_shape_functions (radial_qrule->get_points(),
                                    base_fe->qrule->get_points(),
                                    inf_elem);

      // computing the distance only works when we have the current
      // base_elem stored.  This happens when fe_type is const,
      // the inf_elem->type remains the same.  Then we have to
      // update the base elem _here_.
      if (update_base_elem_required)
        this->update_base_elem(inf_elem);

      // compute dist (depends on geometry, therefore has to be updated for
      // each and every new element), throw radial and base part together
      this->combine_base_radial (inf_elem);

      this->_fe_map->compute_map (this->dim, _total_qrule_weights, inf_elem, this->calculate_d2phi);

      // Compute the shape functions and the derivatives
      // at all quadrature points.
      this->compute_shape_functions (inf_elem,base_fe->qrule->get_points());
    }

  else // if pts != libmesh_nullptr
    {
      // update the elem_type
      elem_type = inf_elem->type();

      // We'll assume that pts is a tensor product mesh of points.
      // That will handle the pts.size()==1 case that we care about
      // right now, and it will generalize a bit, and it won't break
      // the assumptions elsewhere in InfFE.
      std::vector<Point> radial_pts;
      for (std::size_t p=0; p != pts->size(); ++p)
        {
          Real radius = (*pts)[p](Dim-1);
          if (radial_pts.size() && radial_pts[0](0) == radius)
            break;
          radial_pts.push_back(Point(radius));
        }
      const unsigned int radial_pts_size = radial_pts.size();
      const unsigned int base_pts_size = pts->size() / radial_pts_size;
      // If we're a tensor product we should have no remainder
      libmesh_assert_equal_to
        (base_pts_size * radial_pts_size, pts->size());

      std::vector<Point> base_pts;
      base_pts.reserve(base_pts_size);
      for (std::size_t p=0; p != pts->size(); p += radial_pts_size)
        {
          Point pt = (*pts)[p];
          pt(Dim-1) = 0;
          base_pts.push_back(pt);
        }

      // init radial shapes
      this->init_radial_shape_functions(inf_elem, &radial_pts);

      // update the base
      this->update_base_elem(inf_elem);

      // the finite element on the ifem base
      base_fe.reset(FEBase::build(Dim-1, this->fe_type).release());

      base_fe->calculate_phi = base_fe->calculate_dphi = base_fe->calculate_dphiref = true;
      base_fe->get_xyz();
      base_fe->determine_calculations();

      // init base shapes
      base_fe->init_base_shape_functions(base_pts,
                                         base_elem.get());

      // compute the shape functions and map functions of base_fe
      // before using them later in combine_base_radial.

      if (weights)
        {
          base_fe->_fe_map->compute_map (base_fe->dim, *weights,
                                         base_elem.get(), base_fe->calculate_d2phi);
        }
      else
        {
          std::vector<Real> dummy_weights (pts->size(), 1.);
          base_fe->_fe_map->compute_map (base_fe->dim, dummy_weights,
                                         base_elem.get(), base_fe->calculate_d2phi);
        }

      base_fe->compute_shape_functions(base_elem.get(), *pts);

      this->init_shape_functions (radial_pts, base_pts, inf_elem);

      // combine the base and radial shapes
      this->combine_base_radial (inf_elem);

      // weights
      if (weights != libmesh_nullptr)
        {
          this->_fe_map->compute_map (this->dim, *weights, inf_elem, this->calculate_d2phi);
        }
      else
        {
          std::vector<Real> dummy_weights (pts->size(), 1.);
          this->_fe_map->compute_map (this->dim, dummy_weights, inf_elem, this->calculate_d2phi);
        }

      // finally compute the ifem shapes
      this->compute_shape_functions (inf_elem,*pts);
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

  // initialize most of the things related to mapping

  // The order to use in the radial map (currently independent of the element type)
  const Order radial_mapping_order = Radial::mapping_order();
  const unsigned int n_radial_mapping_shape_functions = Radial::n_dofs(radial_mapping_order);

  // initialize most of the things related to physical approximation
  const Order radial_approx_order = fe_type.radial_order;
  const unsigned int n_radial_approx_shape_functions = Radial::n_dofs(radial_approx_order);

  const unsigned int n_radial_qp =
    radial_pts ? radial_pts->size() : radial_qrule->n_points();
  const std::vector<Point> & radial_qp =
    radial_pts ? *radial_pts : radial_qrule->get_points();

  // resize the radial data fields

  // the radial polynomials (eval)
  mode.resize      (n_radial_approx_shape_functions);
  dmodedv.resize   (n_radial_approx_shape_functions);

  // the (1-v)/2 weight
  som.resize       (n_radial_qp);
  dsomdv.resize    (n_radial_qp);

  // the radial map
  radial_map.resize    (n_radial_mapping_shape_functions);
  dradialdv_map.resize (n_radial_mapping_shape_functions);


  for (unsigned int i=0; i<n_radial_mapping_shape_functions; i++)
    {
      radial_map[i].resize    (n_radial_qp);
      dradialdv_map[i].resize (n_radial_qp);
    }


  for (unsigned int i=0; i<n_radial_approx_shape_functions; i++)
    {
      mode[i].resize    (n_radial_qp);
      dmodedv[i].resize (n_radial_qp);
    }


  // compute scalar values at radial quadrature points
  for (unsigned int p=0; p<n_radial_qp; p++)
    {
      som[p] = Radial::decay (radial_qp[p](0));
      dsomdv[p] = Radial::decay_deriv (radial_qp[p](0));
    }


  // evaluate the mode shapes in radial direction at radial quadrature points
  for (unsigned int i=0; i<n_radial_approx_shape_functions; i++)
    for (unsigned int p=0; p<n_radial_qp; p++)
      {
        mode[i][p] = InfFE<Dim,T_radial,T_map>::eval (radial_qp[p](0), radial_approx_order, i);
        dmodedv[i][p] = InfFE<Dim,T_radial,T_map>::eval_deriv (radial_qp[p](0), radial_approx_order, i);
      }


  // evaluate the mapping functions in radial direction at radial quadrature points
  for (unsigned int i=0; i<n_radial_mapping_shape_functions; i++)
    for (unsigned int p=0; p<n_radial_qp; p++)
      {
        radial_map[i][p] = InfFE<Dim,INFINITE_MAP,T_map>::eval (radial_qp[p](0), radial_mapping_order, i);
        dradialdv_map[i][p] = InfFE<Dim,INFINITE_MAP,T_map>::eval_deriv (radial_qp[p](0), radial_mapping_order, i);
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
  const unsigned int n_radial_mapping_sf = cast_int<unsigned int>(radial_map.size());
  const unsigned int n_radial_approx_sf  = cast_int<unsigned int>(mode.size());
  const unsigned int n_radial_qp         = cast_int<unsigned int>(som.size());


  // initialize most of the things related to mapping

  // The element type and order to use in the base map
  const Order base_mapping_order = base_elem->default_order();
  const ElemType base_mapping_elem_type = base_elem->type();

  // the number of base shape functions used to construct the map
  // (Lagrange shape functions are used for mapping in the base)
  unsigned int n_base_mapping_shape_functions = Base::n_base_mapping_sf(base_mapping_elem_type,
                                                                        base_mapping_order);

  const unsigned int n_total_mapping_shape_functions =
    n_radial_mapping_sf * n_base_mapping_shape_functions;

  // initialize most of the things related to physical approximation
  unsigned int n_base_approx_shape_functions;
  if (Dim > 1)
    n_base_approx_shape_functions = base_fe->n_shape_functions();
  else
    n_base_approx_shape_functions = 1;


  const unsigned int n_total_approx_shape_functions =
    n_radial_approx_sf * n_base_approx_shape_functions;

  // update class member field
  _n_total_approx_sf = n_total_approx_shape_functions;


  // The number of the base quadrature points.
  const unsigned int n_base_qp = base_qp.size();

  // The total number of quadrature points.
  const unsigned int n_total_qp = n_radial_qp * n_base_qp;


  // update class member field
  _n_total_qp = n_total_qp;



  // initialize the node and shape numbering maps
  {
    // these vectors work as follows: the i-th entry stores
    // the associated base/radial node number
    _radial_node_index.resize(n_total_mapping_shape_functions);
    _base_node_index.resize(n_total_mapping_shape_functions);

    // similar for the shapes: the i-th entry stores
    // the associated base/radial shape number
    _radial_shape_index.resize(n_total_approx_shape_functions);
    _base_shape_index.resize(n_total_approx_shape_functions);

    const ElemType inf_elem_type = inf_elem->type();

    // fill the node index map
    for (unsigned int n=0; n<n_total_mapping_shape_functions; n++)
      {
        compute_node_indices (inf_elem_type,
                              n,
                              _base_node_index[n],
                              _radial_node_index[n]);
        libmesh_assert_less (_base_node_index[n], n_base_mapping_shape_functions);
        libmesh_assert_less (_radial_node_index[n], n_radial_mapping_sf);
      }

    // fill the shape index map
    for (unsigned int n=0; n<n_total_approx_shape_functions; n++)
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
  dist.resize(n_base_mapping_shape_functions);

  // resize the total data fields

  // the phase term varies with xi, eta and zeta(v): store it for _all_ qp
  //
  // when computing the phase, we need the base approximations
  // therefore, initialize the phase here, but evaluate it
  // in combine_base_radial().
  //
  // the weight, though, is only needed at the radial quadrature points, n_radial_qp.
  // but for a uniform interface to the protected data fields
  // the weight data field (which are accessible from the outside) are expanded to n_total_qp.
  weight.resize      (n_total_qp);
  dweightdv.resize   (n_total_qp);
  dweight.resize     (n_total_qp);

  dphase.resize      (n_total_qp);
  dphasedxi.resize   (n_total_qp);
  dphasedeta.resize  (n_total_qp);
  dphasedzeta.resize (n_total_qp);

  // this vector contains the integration weights for the combined quadrature rule
  _total_qrule_weights.resize(n_total_qp);

  // InfFE's data fields phi, dphi, dphidx, phi_map etc hold the _total_
  // shape and mapping functions, respectively
  {
    phi.resize     (n_total_approx_shape_functions);
    dphi.resize    (n_total_approx_shape_functions);
    dphidx.resize  (n_total_approx_shape_functions);
    dphidy.resize  (n_total_approx_shape_functions);
    dphidz.resize  (n_total_approx_shape_functions);
    dphidxi.resize (n_total_approx_shape_functions);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    libmesh_do_once(libMesh::err << "Second derivatives for Infinite elements"
                    << " are not yet implemented!"
                    << std::endl);

    d2phi.resize     (n_total_approx_shape_functions);
    d2phidx2.resize  (n_total_approx_shape_functions);
    d2phidxdy.resize (n_total_approx_shape_functions);
    d2phidxdz.resize (n_total_approx_shape_functions);
    d2phidy2.resize  (n_total_approx_shape_functions);
    d2phidydz.resize (n_total_approx_shape_functions);
    d2phidz2.resize  (n_total_approx_shape_functions);
    d2phidxi2.resize (n_total_approx_shape_functions);

    if (Dim > 1)
      {
        d2phidxideta.resize(n_total_approx_shape_functions);
        d2phideta2.resize(n_total_approx_shape_functions);
      }

    if (Dim > 2)
      {
        d2phidetadzeta.resize(n_total_approx_shape_functions);
        d2phidxidzeta.resize(n_total_approx_shape_functions);
        d2phidzeta2.resize(n_total_approx_shape_functions);
      }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

    if (Dim > 1)
      dphideta.resize(n_total_approx_shape_functions);

    if (Dim == 3)
      dphidzeta.resize(n_total_approx_shape_functions);

    std::vector<std::vector<Real> > & phi_map = this->_fe_map->get_phi_map();
    std::vector<std::vector<Real> > & dphidxi_map = this->_fe_map->get_dphidxi_map();

    phi_map.resize(n_total_mapping_shape_functions);
    dphidxi_map.resize(n_total_mapping_shape_functions);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    std::vector<std::vector<Real> > & d2phidxi2_map = this->_fe_map->get_d2phidxi2_map();
    d2phidxi2_map.resize(n_total_mapping_shape_functions);

    if (Dim > 1)
      {
        std::vector<std::vector<Real> > & d2phidxideta_map = this->_fe_map->get_d2phidxideta_map();
        std::vector<std::vector<Real> > & d2phideta2_map = this->_fe_map->get_d2phideta2_map();
        d2phidxideta_map.resize(n_total_mapping_shape_functions);
        d2phideta2_map.resize(n_total_mapping_shape_functions);
      }

    if (Dim == 3)
      {
        std::vector<std::vector<Real> > & d2phidxidzeta_map = this->_fe_map->get_d2phidxidzeta_map();
        std::vector<std::vector<Real> > & d2phidetadzeta_map = this->_fe_map->get_d2phidetadzeta_map();
        std::vector<std::vector<Real> > & d2phidzeta2_map = this->_fe_map->get_d2phidzeta2_map();
        d2phidxidzeta_map.resize(n_total_mapping_shape_functions);
        d2phidetadzeta_map.resize(n_total_mapping_shape_functions);
        d2phidzeta2_map.resize(n_total_mapping_shape_functions);
      }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

    if (Dim > 1)
      {
        std::vector<std::vector<Real> > & dphideta_map = this->_fe_map->get_dphideta_map();
        dphideta_map.resize(n_total_mapping_shape_functions);
      }

    if (Dim == 3)
      {
        std::vector<std::vector<Real> > & dphidzeta_map = this->_fe_map->get_dphidzeta_map();
        dphidzeta_map.resize(n_total_mapping_shape_functions);
      }
  }

  // collect all the for loops, where inner vectors are
  // resized to the appropriate number of quadrature points
  {
    for (unsigned int i=0; i<n_total_approx_shape_functions; i++)
      {
        phi[i].resize         (n_total_qp);
        dphi[i].resize        (n_total_qp);
        dphidx[i].resize      (n_total_qp);
        dphidy[i].resize      (n_total_qp);
        dphidz[i].resize      (n_total_qp);
        dphidxi[i].resize     (n_total_qp);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        d2phi[i].resize       (n_total_qp);
        d2phidx2[i].resize    (n_total_qp);
        d2phidxdy[i].resize   (n_total_qp);
        d2phidxdz[i].resize   (n_total_qp);
        d2phidy2[i].resize    (n_total_qp);
        d2phidydz[i].resize   (n_total_qp);
        d2phidy2[i].resize    (n_total_qp);
        d2phidxi2[i].resize   (n_total_qp);

        if (Dim > 1)
          {
            d2phidxideta[i].resize   (n_total_qp);
            d2phideta2[i].resize     (n_total_qp);
          }
        if (Dim > 2)
          {
            d2phidxidzeta[i].resize  (n_total_qp);
            d2phidetadzeta[i].resize (n_total_qp);
            d2phidzeta2[i].resize    (n_total_qp);
          }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        if (Dim > 1)
          dphideta[i].resize  (n_total_qp);

        if (Dim == 3)
          dphidzeta[i].resize (n_total_qp);

      }

    for (unsigned int i=0; i<n_total_mapping_shape_functions; i++)
      {
        std::vector<std::vector<Real> > & phi_map = this->_fe_map->get_phi_map();
        std::vector<std::vector<Real> > & dphidxi_map = this->_fe_map->get_dphidxi_map();
        phi_map[i].resize         (n_total_qp);
        dphidxi_map[i].resize     (n_total_qp);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        std::vector<std::vector<Real> > & d2phidxi2_map = this->_fe_map->get_d2phidxi2_map();
        d2phidxi2_map[i].resize   (n_total_qp);
        if (Dim > 1)
          {
            std::vector<std::vector<Real> > & d2phidxideta_map = this->_fe_map->get_d2phidxideta_map();
            std::vector<std::vector<Real> > & d2phideta2_map = this->_fe_map->get_d2phideta2_map();
            d2phidxideta_map[i].resize   (n_total_qp);
            d2phideta2_map[i].resize     (n_total_qp);
          }

        if (Dim > 2)
          {
            std::vector<std::vector<Real> > & d2phidxidzeta_map = this->_fe_map->get_d2phidxidzeta_map();
            std::vector<std::vector<Real> > & d2phidetadzeta_map = this->_fe_map->get_d2phidetadzeta_map();
            std::vector<std::vector<Real> > & d2phidzeta2_map = this->_fe_map->get_d2phidzeta2_map();
            d2phidxidzeta_map[i].resize  (n_total_qp);
            d2phidetadzeta_map[i].resize (n_total_qp);
            d2phidzeta2_map[i].resize    (n_total_qp);
          }
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

        if (Dim > 1)
          {
            std::vector<std::vector<Real> > & dphideta_map = this->_fe_map->get_dphideta_map();
            dphideta_map[i].resize  (n_total_qp);
          }

        if (Dim == 3)
          {
            std::vector<std::vector<Real> > & dphidzeta_map = this->_fe_map->get_dphidzeta_map();
            dphidzeta_map[i].resize (n_total_qp);
          }
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

        for (unsigned int rp=0; rp<n_radial_qp; rp++)
          for (unsigned int bp=0; bp<n_base_qp; bp++)
            {
              weight[bp + rp*n_base_qp] = Radial::D(radial_qp[rp](0));
              dweightdv[bp + rp*n_base_qp] = Radial::D_deriv(radial_qp[rp](0));
              _total_qrule_weights[bp + rp*n_base_qp] = radial_qw[rp] * base_qw[bp];
            }
      }
    else
      {
        for (unsigned int rp=0; rp<n_radial_qp; rp++)
          for (unsigned int bp=0; bp<n_base_qp; bp++)
            {
              weight[bp + rp*n_base_qp] = Radial::D(radial_qp[rp](0));
              dweightdv[bp + rp*n_base_qp] = Radial::D_deriv(radial_qp[rp](0));
            }
      }
  }
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::combine_base_radial(const Elem * inf_elem)
{
  libmesh_assert(inf_elem);
  // at least check whether the base element type is correct.
  // otherwise this version of computing dist would give problems
  libmesh_assert_equal_to (base_elem->type(), Base::get_elem_type(inf_elem->type()));

  // Start logging the combination of radial and base parts
  LOG_SCOPE("combine_base_radial()", "InfFE");

  // zero  the phase, since it is to be summed up
  std::fill (dphasedxi.begin(),   dphasedxi.end(),   0.);
  std::fill (dphasedeta.begin(),  dphasedeta.end(),  0.);
  std::fill (dphasedzeta.begin(), dphasedzeta.end(), 0.);


  const unsigned int n_base_mapping_sf = cast_int<unsigned int>(dist.size());
  const Point origin = inf_elem->origin();

  // for each new infinite element, compute the radial distances
  for (unsigned int n=0; n<n_base_mapping_sf; n++)
    dist[n] = Point(base_elem->point(n) - origin).norm();


  switch (Dim)
    {
      // 1D
    case 1:
      {
        libmesh_not_implemented();
        break;
      }

      // 2D
    case 2:
      {
        libmesh_not_implemented();
        break;
      }

      // 3D
    case 3:
      {
        // fast access to the approximation and mapping shapes of base_fe
        const std::vector<std::vector<Real> > & S  = base_fe->phi;
        const std::vector<std::vector<Real> > & Ss = base_fe->dphidxi;
        const std::vector<std::vector<Real> > & St = base_fe->dphideta;
        const std::vector<std::vector<Real> > & S_map  = (base_fe->get_fe_map()).get_phi_map();
        const std::vector<std::vector<Real> > & Ss_map = (base_fe->get_fe_map()).get_dphidxi_map();
        const std::vector<std::vector<Real> > & St_map = (base_fe->get_fe_map()).get_dphideta_map();

        const unsigned int n_radial_qp = som.size();
        if (radial_qrule)
          libmesh_assert_equal_to(n_radial_qp, radial_qrule->n_points());
        const unsigned int n_base_qp = S_map[0].size();
        if (base_qrule)
          libmesh_assert_equal_to(n_base_qp, base_qrule->n_points());

        const unsigned int n_total_mapping_sf =
          cast_int<unsigned int>(radial_map.size()) * n_base_mapping_sf;

        const unsigned int n_total_approx_sf = Radial::n_dofs(fe_type.radial_order) *  base_fe->n_shape_functions();


        // compute the phase term derivatives
        {
          unsigned int tp=0;
          for (unsigned int rp=0; rp<n_radial_qp; rp++)  // over radial qps
            for (unsigned int bp=0; bp<n_base_qp; bp++)  // over base qps
              {
                // sum over all base shapes, to get the average distance
                for (unsigned int i=0; i<n_base_mapping_sf; i++)
                  {
                    dphasedxi[tp]   += Ss_map[i][bp] * dist[i] * radial_map   [1][rp];
                    dphasedeta[tp]  += St_map[i][bp] * dist[i] * radial_map   [1][rp];
                    dphasedzeta[tp] += S_map [i][bp] * dist[i] * dradialdv_map[1][rp];
                  }

                tp++;

              } // loop radial and base qps
        }

        libmesh_assert_equal_to (phi.size(), n_total_approx_sf);
        libmesh_assert_equal_to (dphidxi.size(), n_total_approx_sf);
        libmesh_assert_equal_to (dphideta.size(), n_total_approx_sf);
        libmesh_assert_equal_to (dphidzeta.size(), n_total_approx_sf);

        // compute the overall approximation shape functions,
        // pick the appropriate radial and base shapes through using
        // _base_shape_index and _radial_shape_index
        for (unsigned int rp=0; rp<n_radial_qp; rp++)  // over radial qps
          for (unsigned int bp=0; bp<n_base_qp; bp++)  // over base qps
            for (unsigned int ti=0; ti<n_total_approx_sf; ti++)  // over _all_ approx_sf
              {
                // let the index vectors take care of selecting the appropriate base/radial shape
                const unsigned int bi = _base_shape_index  [ti];
                const unsigned int ri = _radial_shape_index[ti];
                phi      [ti][bp+rp*n_base_qp] = S [bi][bp] * mode[ri][rp] * som[rp];
                dphidxi  [ti][bp+rp*n_base_qp] = Ss[bi][bp] * mode[ri][rp] * som[rp];
                dphideta [ti][bp+rp*n_base_qp] = St[bi][bp] * mode[ri][rp] * som[rp];
                dphidzeta[ti][bp+rp*n_base_qp] = S [bi][bp]
                  * (dmodedv[ri][rp] * som[rp] + mode[ri][rp] * dsomdv[rp]);
              }

        std::vector<std::vector<Real> > & phi_map = this->_fe_map->get_phi_map();
        std::vector<std::vector<Real> > & dphidxi_map = this->_fe_map->get_dphidxi_map();
        std::vector<std::vector<Real> > & dphideta_map = this->_fe_map->get_dphideta_map();
        std::vector<std::vector<Real> > & dphidzeta_map = this->_fe_map->get_dphidzeta_map();

        libmesh_assert_equal_to (phi_map.size(), n_total_mapping_sf);
        libmesh_assert_equal_to (dphidxi_map.size(), n_total_mapping_sf);
        libmesh_assert_equal_to (dphideta_map.size(), n_total_mapping_sf);
        libmesh_assert_equal_to (dphidzeta_map.size(), n_total_mapping_sf);

        // compute the overall mapping functions,
        // pick the appropriate radial and base entries through using
        // _base_node_index and _radial_node_index
        for (unsigned int rp=0; rp<n_radial_qp; rp++) // over radial qps
          for (unsigned int bp=0; bp<n_base_qp; bp++) // over base qps
            for (unsigned int ti=0; ti<n_total_mapping_sf; ti++) // over all mapping shapes
              {
                // let the index vectors take care of selecting the appropriate base/radial mapping shape
                const unsigned int bi = _base_node_index  [ti];
                const unsigned int ri = _radial_node_index[ti];
                phi_map      [ti][bp+rp*n_base_qp] = S_map [bi][bp] * radial_map   [ri][rp];
                dphidxi_map  [ti][bp+rp*n_base_qp] = Ss_map[bi][bp] * radial_map   [ri][rp];
                dphideta_map [ti][bp+rp*n_base_qp] = St_map[bi][bp] * radial_map   [ri][rp];
                dphidzeta_map[ti][bp+rp*n_base_qp] = S_map [bi][bp] * dradialdv_map[ri][rp];
              }

        break;
      }

    default:
      libmesh_error_msg("Unsupported Dim = " << Dim);
    }
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_shape_functions(const Elem *,
                                                        const std::vector<Point> &)
{
  // Start logging the overall computation of shape functions
  LOG_SCOPE("compute_shape_functions()", "InfFE");

  const unsigned int n_total_qp  = _n_total_qp;

  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Compute the value of the derivative shape function i at quadrature point p
  switch (dim)
    {

    case 1:
      {
        libmesh_not_implemented();
        break;
      }

    case 2:
      {
        libmesh_not_implemented();
        break;
      }

    case 3:
      {
        const std::vector<Real> & dxidx_map = this->_fe_map->get_dxidx();
        const std::vector<Real> & dxidy_map = this->_fe_map->get_dxidy();
        const std::vector<Real> & dxidz_map = this->_fe_map->get_dxidz();

        const std::vector<Real> & detadx_map = this->_fe_map->get_detadx();
        const std::vector<Real> & detady_map = this->_fe_map->get_detady();
        const std::vector<Real> & detadz_map = this->_fe_map->get_detadz();

        const std::vector<Real> & dzetadx_map = this->_fe_map->get_dzetadx();
        const std::vector<Real> & dzetady_map = this->_fe_map->get_dzetady();
        const std::vector<Real> & dzetadz_map = this->_fe_map->get_dzetadz();

        // These are _all_ shape functions of this infinite element
        for (std::size_t i=0; i<phi.size(); i++)
          for (unsigned int p=0; p<n_total_qp; p++)
            {
              // dphi/dx    = (dphi/dxi)*(dxi/dx) + (dphi/deta)*(deta/dx) + (dphi/dzeta)*(dzeta/dx);
              dphi[i][p](0) =
                dphidx[i][p] = (dphidxi[i][p]*dxidx_map[p] +
                                dphideta[i][p]*detadx_map[p] +
                                dphidzeta[i][p]*dzetadx_map[p]);

              // dphi/dy    = (dphi/dxi)*(dxi/dy) + (dphi/deta)*(deta/dy) + (dphi/dzeta)*(dzeta/dy);
              dphi[i][p](1) =
                dphidy[i][p] = (dphidxi[i][p]*dxidy_map[p] +
                                dphideta[i][p]*detady_map[p] +
                                dphidzeta[i][p]*dzetady_map[p]);

              // dphi/dz    = (dphi/dxi)*(dxi/dz) + (dphi/deta)*(deta/dz) + (dphi/dzeta)*(dzeta/dz);
              dphi[i][p](2) =
                dphidz[i][p] = (dphidxi[i][p]*dxidz_map[p] +
                                dphideta[i][p]*detadz_map[p] +
                                dphidzeta[i][p]*dzetadz_map[p]);
            }


        // This is the derivative of the phase term of this infinite element
        for (unsigned int p=0; p<n_total_qp; p++)
          {
            // the derivative of the phase term
            dphase[p](0) = (dphasedxi[p]   * dxidx_map[p] +
                            dphasedeta[p]  * detadx_map[p] +
                            dphasedzeta[p] * dzetadx_map[p]);

            dphase[p](1) = (dphasedxi[p]   * dxidy_map[p] +
                            dphasedeta[p]  * detady_map[p] +
                            dphasedzeta[p] * dzetady_map[p]);

            dphase[p](2) = (dphasedxi[p]   * dxidz_map[p] +
                            dphasedeta[p]  * detadz_map[p] +
                            dphasedzeta[p] * dzetadz_map[p]);

            // the derivative of the radial weight - varies only in radial direction,
            // therefore dweightdxi = dweightdeta = 0.
            dweight[p](0) = dweightdv[p] * dzetadx_map[p];

            dweight[p](1) = dweightdv[p] * dzetady_map[p];

            dweight[p](2) = dweightdv[p] * dzetadz_map[p];
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
