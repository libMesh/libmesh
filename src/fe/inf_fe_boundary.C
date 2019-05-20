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
                                        const Real tolerance,
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

  // We need to allow this as well...
  // libmesh_assert_not_equal_to (s, 0);


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

          //FIXME: Do I have to care about the order of my neighbours element?
          //unsigned int side_p_level = elem->p_level();
          //if (elem->neighbor_ptr(s) != nullptr)
          //  side_p_level = std::max(side_p_level, elem->neighbor_ptr(s)->p_level());
          base_qrule->init(side->type(), side->p_level());
        }
      radial_qrule_initialized = true;
    }

  // Initialize the face shape functions
  if (this->get_type() != inf_elem->type() ||
      base_fe->shapes_need_reinit()        ||
      radial_qrule_initialized)
    this->init_face_shape_functions (qrule->get_points(), side.get());


  // compute the face map
  this->_fe_map->compute_face_map(this->dim, _total_qrule_weights, side.get());

  // make a copy of the Jacobian for integration
  const std::vector<Real> JxW_int(this->_fe_map->get_JxW());

  // Find where the integration points are located on the
  // full element.
  std::vector<Point> qp;
  this->inverse_map (inf_elem, this->_fe_map->get_xyz(), qp, tolerance);

  // just to ensure that we are working on the base and not, for numeric reasons,
  // somewhere else...
  if (s==0)
    {
      for (unsigned int p=0; p<qp.size(); p++)
        {
          qp[p](Dim-1)=-1.;
        }
    }

  // compute the shape function and derivative values
  // at the points qp
  this->reinit  (inf_elem, &qp);

  // copy back old data
  this->_fe_map->get_JxW() = JxW_int;

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

  // Initialize the radial shape functions
  this->init_radial_shape_functions(inf_side);

  // Initialize the base shape functions
  if (inf_side->infinite())
    this->update_base_elem(inf_side);
  else
    // in this case, I need the 2D base
    this->update_base_elem(inf_side->parent());

  // Initialize the base quadrature rule
  base_qrule->init(base_elem->type(), inf_side->p_level());

  // base_fe still corresponds to the (dim-1)-dimensional base of the InfFE object,
  // so update the fe_base.
  if (inf_side->infinite())
    {
      libmesh_assert_equal_to (Dim, 3);
      base_fe = FEBase::build(Dim-2, this->fe_type);
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

  // the quadrature weights
  _total_qrule_weights.resize(n_total_qp);

  // now init the shapes for boundary work
  {

    // The element type and order to use in the base map
    const Order    base_mapping_order     ( base_elem->default_order() );
    const ElemType base_mapping_elem_type ( base_elem->type()          );

    // the number of mapping shape functions. For base side it is 1.
    // (Lagrange shape functions are used for mapping in the base)
    const unsigned int n_radial_mapping_sf =
      inf_side->infinite() ? cast_int<unsigned int>(radial_map.size()) : 1;

    const unsigned int n_base_mapping_shape_functions = Base::n_base_mapping_sf(base_mapping_elem_type,
                                                                                base_mapping_order);

    const unsigned int n_total_mapping_shape_functions =
      n_radial_mapping_sf * n_base_mapping_shape_functions;


    // initialize the node and shape numbering maps
    _radial_node_index.resize    (n_total_mapping_shape_functions);
    _base_node_index.resize      (n_total_mapping_shape_functions);

    if (inf_side->infinite())
      {
        const ElemType inf_face_elem_type (inf_side->type());

        // fill the node index map
        for (unsigned int n=0; n<n_total_mapping_shape_functions; n++)
          {
            compute_node_indices (inf_face_elem_type,
                                  n,
                                  _base_node_index[n],
                                  _radial_node_index[n]);

            libmesh_assert_less (_base_node_index[n], n_base_mapping_shape_functions);
            libmesh_assert_less (_radial_node_index[n], n_radial_mapping_sf);
          }
      }
    else
      {
        for (unsigned int n=0; n<n_total_mapping_shape_functions; n++)
          {
            _base_node_index[n] = n;
            _radial_node_index[n] = 0;
          }
      }

    // resize map data fields
    std::vector<std::vector<Real>> & psi_map = this->_fe_map->get_psi();
    std::vector<std::vector<Real>> & dpsidxi_map = this->_fe_map->get_dpsidxi();
    std::vector<std::vector<Real>> & dpsideta_map = this->_fe_map->get_dpsideta();
    psi_map.resize          (n_total_mapping_shape_functions);
    dpsidxi_map.resize      (n_total_mapping_shape_functions);
    dpsideta_map.resize     (n_total_mapping_shape_functions);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    std::vector<std::vector<Real>> & d2psidxi2_map = this->_fe_map->get_d2psidxi2();
    std::vector<std::vector<Real>> & d2psidxideta_map = this->_fe_map->get_d2psidxideta();
    std::vector<std::vector<Real>> & d2psideta2_map = this->_fe_map->get_d2psideta2();
    d2psidxi2_map.resize    (n_total_mapping_shape_functions);
    d2psidxideta_map.resize (n_total_mapping_shape_functions);
    d2psideta2_map.resize   (n_total_mapping_shape_functions);
#endif

    for (unsigned int i=0; i<n_total_mapping_shape_functions; i++)
      {
        psi_map[i].resize         (n_total_qp);
        dpsidxi_map[i].resize     (n_total_qp);
        dpsideta_map[i].resize     (n_total_qp);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        d2psidxi2_map[i].resize   (n_total_qp);
#endif
      }

    // compute shape maps
    if (inf_side->infinite())
      {
        const std::vector<std::vector<Real>> & S_map  = (base_fe->get_fe_map()).get_phi_map();
        const std::vector<std::vector<Real>> & Ss_map = (base_fe->get_fe_map()).get_dphidxi_map();

        for (unsigned int rp=0; rp<n_radial_qp; rp++)  // over radial qps
          for (unsigned int bp=0; bp<n_base_qp; bp++)  // over base qps
            for (unsigned int ti=0; ti<n_total_mapping_shape_functions; ti++)  // over all mapping shapes
              {
                // let the index vectors take care of selecting the appropriate base/radial mapping shape
                const unsigned int bi = _base_node_index  [ti];
                const unsigned int ri = _radial_node_index[ti];
                psi_map          [ti][bp+rp*n_base_qp] = S_map [bi][bp] * radial_map   [ri][rp];
                dpsidxi_map      [ti][bp+rp*n_base_qp] = Ss_map[bi][bp] * radial_map   [ri][rp];
                dpsideta_map     [ti][bp+rp*n_base_qp] = S_map [bi][bp] * dradialdv_map[ri][rp];

                // second derivatives are not implemented for infinite elements
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                d2psidxi2_map    [ti][bp+rp*n_base_qp] = 0.;
                d2psidxideta_map [ti][bp+rp*n_base_qp] = 0.;
                d2psideta2_map   [ti][bp+rp*n_base_qp] = 0.;
#endif
              }

      }
    else
      {
        const std::vector<std::vector<Real>> & S_map  = (base_fe->get_fe_map()).get_phi_map();
        const std::vector<std::vector<Real>> & Ss_map = (base_fe->get_fe_map()).get_dphidxi_map();
        const std::vector<std::vector<Real>> & St_map = (base_fe->get_fe_map()).get_dphideta_map();
        for (unsigned int bp=0; bp<n_base_qp; bp++)  // over base qps
          for (unsigned int ti=0; ti<n_total_mapping_shape_functions; ti++)  // over all mapping shapes
            {
              psi_map      [ti][bp] = S_map[ti][bp] ;
              dpsidxi_map  [ti][bp] = Ss_map[ti][bp] ;
              dpsideta_map [ti][bp] = St_map[ti][bp] ;
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              d2psidxi2_map    [ti][bp] = 0.;
              d2psidxideta_map [ti][bp] = 0.;
              d2psideta2_map   [ti][bp] = 0.;
#endif
            }
      }
  }

  // quadrature rule weights
  {
    const std::vector<Real> & radial_qw = radial_qrule->get_weights();
    const std::vector<Real> & base_qw   = base_qrule->get_weights();

    libmesh_assert_equal_to (radial_qw.size(), n_radial_qp);
    libmesh_assert_equal_to (base_qw.size(), n_base_qp);

    for (unsigned int rp=0; rp<n_radial_qp; rp++)
      for (unsigned int bp=0; bp<n_base_qp; bp++)
        {
          _total_qrule_weights[bp + rp*n_base_qp] = radial_qw[rp] * base_qw[bp];
        }
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

} // namespace libMesh

#endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
