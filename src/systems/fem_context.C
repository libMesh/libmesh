// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/boundary_info.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/time_solver.h"
#include "libmesh/unsteady_solver.h" // For euler_residual

namespace libMesh
{

FEMContext::FEMContext (const System &sys)
  : DiffContext(sys),
    _mesh_sys(NULL),
    _mesh_x_var(0),
    _mesh_y_var(0),
    _mesh_z_var(0),
    side(0), edge(0),
    _boundary_info(sys.get_mesh().boundary_info.get()),
    elem(NULL),
    dim(sys.get_mesh().mesh_dimension()),
    element_qrule(NULL), side_qrule(NULL),
    edge_qrule(NULL)
{
  // We need to know which of our variables has the hardest
  // shape functions to numerically integrate.

  unsigned int nv = sys.n_vars();

  libmesh_assert (nv);
  FEType hardest_fe_type = sys.variable_type(0);

  for (unsigned int i=0; i != nv; ++i)
    {
      FEType fe_type = sys.variable_type(i);

      // FIXME - we don't yet handle mixed finite elements from
      // different families which require different quadrature rules
      // libmesh_assert_equal_to (fe_type.family, hardest_fe_type.family);

      if (fe_type.order > hardest_fe_type.order)
        hardest_fe_type = fe_type;
    }

  // Create an adequate quadrature rule
  element_qrule = hardest_fe_type.default_quadrature_rule
    (dim, sys.extra_quadrature_order).release();
  side_qrule = hardest_fe_type.default_quadrature_rule
    (dim-1, sys.extra_quadrature_order).release();
  if (dim == 3)
    edge_qrule = hardest_fe_type.default_quadrature_rule
      (1, sys.extra_quadrature_order).release();

  // Next, create finite element objects
  // Preserving backward compatibility here for now
  // Should move to the protected/FEAbstract interface
  _element_fe_var.resize(nv);
  _side_fe_var.resize(nv);
  if (dim == 3)
    _edge_fe_var.resize(nv);

  for (unsigned int i=0; i != nv; ++i)
    {
      FEType fe_type = sys.variable_type(i);

      if ( _element_fe[fe_type] == NULL )
        {
          _element_fe[fe_type] = FEAbstract::build(dim, fe_type).release();
          _element_fe[fe_type]->attach_quadrature_rule(element_qrule);
          _side_fe[fe_type] = FEAbstract::build(dim, fe_type).release();
          _side_fe[fe_type]->attach_quadrature_rule(side_qrule);

          if (dim == 3)
            {
              _edge_fe[fe_type] = FEAbstract::build(dim, fe_type).release();
              _edge_fe[fe_type]->attach_quadrature_rule(edge_qrule);
            }
        }
      _element_fe_var[i] = _element_fe[fe_type];
      _side_fe_var[i] = _side_fe[fe_type];
      if (dim == 3)
        _edge_fe_var[i] = _edge_fe[fe_type];

    }
}


FEMContext::~FEMContext()
{
  // We don't want to store AutoPtrs in STL containers, but we don't
  // want to leak memory either
  for (std::map<FEType, FEAbstract *>::iterator i = _element_fe.begin();
       i != _element_fe.end(); ++i)
    delete i->second;
  _element_fe.clear();

  for (std::map<FEType, FEAbstract *>::iterator i = _side_fe.begin();
       i != _side_fe.end(); ++i)
    delete i->second;
  _side_fe.clear();

  for (std::map<FEType, FEAbstract *>::iterator i = _edge_fe.begin();
       i != _edge_fe.end(); ++i)
    delete i->second;
  _edge_fe.clear();

  delete element_qrule;
  element_qrule = NULL;

  delete side_qrule;
  side_qrule = NULL;

  delete edge_qrule;
  side_qrule = NULL;
}



bool FEMContext::has_side_boundary_id(boundary_id_type id) const
{
  return _boundary_info->has_boundary_id(elem, side, id);
}


std::vector<boundary_id_type> FEMContext::side_boundary_ids() const
{
  return _boundary_info->boundary_ids(elem, side);
}



Number FEMContext::interior_value(unsigned int var, unsigned int qp) const
{
  Number u = 0.;

  this->interior_value( var, qp, u );

  return u;
}

template<typename OutputType>
void FEMContext::interior_value(unsigned int var, unsigned int qp,
                                OutputType& u) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > &phi = fe->get_phi();

  // Accumulate solution value
  u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return;
}


template<typename OutputType>
void FEMContext::interior_values (unsigned int var,
                                  const NumericVector<Number> & _system_vector,
                                  std::vector<OutputType>& u_vals) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  const DenseSubVector<Number> &coef = get_localized_subvector(_system_vector, var);

  // Get the finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > &phi = fe->get_phi();

  // Loop over all the q_points on this element
  for (unsigned int qp=0; qp != u_vals.size(); qp++)
    {
      OutputType &u = u_vals[qp];

      // Compute the value at this q_point
      u = 0.;

      for (unsigned int l=0; l != n_dofs; l++)
        u += phi[l][qp] * coef(l);
    }

  return;
}

Gradient FEMContext::interior_gradient(unsigned int var, unsigned int qp) const
{
  Gradient du;

  this->interior_gradient( var, qp, du );

  return du;
}



template<typename OutputType>
void FEMContext::interior_gradient(unsigned int var, unsigned int qp,
                                   OutputType& du) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > &dphi = fe->get_dphi();

  // Accumulate solution derivatives
  du = 0;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return;
}



template<typename OutputType>
void FEMContext::interior_gradients
(unsigned int var,
 const NumericVector<Number> & _system_vector,
 std::vector<OutputType>& du_vals) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  const DenseSubVector<Number> &coef = get_localized_subvector(_system_vector, var);

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > &dphi = fe->get_dphi();

  // Loop over all the q_points in this finite element
  for (unsigned int qp=0; qp != du_vals.size(); qp++)
    {
      OutputType &du = du_vals[qp];

      // Compute the gradient at this q_point
      du = 0;

      for (unsigned int l=0; l != n_dofs; l++)
        du.add_scaled(dphi[l][qp], coef(l));
    }

  return;
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::interior_hessian(unsigned int var, unsigned int qp) const
{
  Tensor d2u;

  this->interior_hessian( var, qp, d2u );

  return d2u;
}

template<typename OutputType>
void FEMContext::interior_hessian(unsigned int var, unsigned int qp,
                                  OutputType& d2u) const
{
  typedef typename TensorTools::MakeReal<
    typename TensorTools::DecrementRank<
    typename TensorTools::DecrementRank<
    OutputType>::type>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &d2phi = fe->get_d2phi();

  // Accumulate solution second derivatives
  d2u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return;
}


template<typename OutputType>
void FEMContext::interior_hessians
(unsigned int var,
 const NumericVector<Number> & _system_vector,
 std::vector<OutputType>& d2u_vals) const
{
  typedef typename TensorTools::MakeReal<
    typename TensorTools::DecrementRank<
    typename TensorTools::DecrementRank<
    OutputType>::type>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  const DenseSubVector<Number> &coef = get_localized_subvector(_system_vector, var);

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &d2phi = fe->get_d2phi();

  // Loop over all the q_points in this finite element
  for (unsigned int qp=0; qp != d2u_vals.size(); qp++)
    {
      OutputType &d2u = d2u_vals[qp];

      // Compute the gradient at this q_point
      d2u = 0;

      for (unsigned int l=0; l != n_dofs; l++)
        d2u.add_scaled(d2phi[l][qp], coef(l));
    }

  return;
}


#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES


template<typename OutputType>
void FEMContext::interior_curl(unsigned int var, unsigned int qp,
                               OutputType& curl_u) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputShape> > &curl_phi = fe->get_curl_phi();

  // Accumulate solution curl
  curl_u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    curl_u.add_scaled(curl_phi[l][qp], coef(l));

  return;
}


template<typename OutputType>
void FEMContext::interior_div(unsigned int var, unsigned int qp,
                              OutputType& div_u) const
{
  typedef typename
    TensorTools::IncrementRank
    <typename TensorTools::MakeReal<OutputType>::type>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputDivergence> > &div_phi = fe->get_div_phi();

  // Accumulate solution curl
  div_u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    div_u += div_phi[l][qp] * coef(l);

  return;
}


Number FEMContext::side_value(unsigned int var, unsigned int qp) const
{
  Number u = 0.;

  this->side_value( var, qp, u );

  return u;
}


template<typename OutputType>
void FEMContext::side_value(unsigned int var, unsigned int qp,
                            OutputType& u) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* the_side_fe = NULL;
  this->get_side_fe<OutputShape>( var, the_side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > &phi = the_side_fe->get_phi();

  // Accumulate solution value
  u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return;
}


template<typename OutputType>
void FEMContext::side_values
(unsigned int var,
 const NumericVector<Number> & _system_vector,
 std::vector<OutputType>& u_vals) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  const DenseSubVector<Number> &coef = get_localized_subvector(_system_vector, var);

  // Get the finite element object
  FEGenericBase<OutputShape>* the_side_fe = NULL;
  this->get_side_fe<OutputShape>( var, the_side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > &phi = the_side_fe->get_phi();

  // Loop over all the q_points on this element
  for (unsigned int qp=0; qp != u_vals.size(); qp++)
    {
      OutputType &u = u_vals[qp];

      // Compute the value at this q_point
      u = 0.;

      for (unsigned int l=0; l != n_dofs; l++)
        u += phi[l][qp] * coef(l);
    }

  return;
}

Gradient FEMContext::side_gradient(unsigned int var, unsigned int qp) const
{
  Gradient du;

  this->side_gradient( var, qp, du );

  return du;
}


template<typename OutputType>
void FEMContext::side_gradient(unsigned int var, unsigned int qp,
                               OutputType& du) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* the_side_fe = NULL;
  this->get_side_fe<OutputShape>( var, the_side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector< typename FEGenericBase<OutputShape>::OutputGradient> > &dphi = the_side_fe->get_dphi();

  // Accumulate solution derivatives
  du = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return;
}



template<typename OutputType>
void FEMContext::side_gradients
(unsigned int var,
 const NumericVector<Number> & _system_vector,
 std::vector<OutputType>& du_vals) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  const DenseSubVector<Number> &coef = get_localized_subvector(_system_vector, var);

  // Get finite element object
  FEGenericBase<OutputShape>* the_side_fe = NULL;
  this->get_side_fe<OutputShape>( var, the_side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > &dphi = the_side_fe->get_dphi();

  // Loop over all the q_points in this finite element
  for (unsigned int qp=0; qp != du_vals.size(); qp++)
    {
      OutputType &du = du_vals[qp];

      du = 0;

      // Compute the gradient at this q_point
      for (unsigned int l=0; l != n_dofs; l++)
        du.add_scaled(dphi[l][qp], coef(l));
    }

  return;
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::side_hessian(unsigned int var, unsigned int qp) const
{
  Tensor d2u;

  this->side_hessian( var, qp, d2u );

  return d2u;
}

template<typename OutputType>
void FEMContext::side_hessian(unsigned int var, unsigned int qp,
                              OutputType& d2u) const
{
  typedef typename TensorTools::MakeReal<
    typename TensorTools::DecrementRank<
    typename TensorTools::DecrementRank<
    OutputType>::type>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* the_side_fe = NULL;
  this->get_side_fe<OutputShape>( var, the_side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &d2phi = the_side_fe->get_d2phi();

  // Accumulate solution second derivatives
  d2u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return;
}


template<typename OutputType>
void FEMContext::side_hessians
(unsigned int var,
 const NumericVector<Number> & _system_vector,
 std::vector<OutputType>& d2u_vals) const
{
  typedef typename TensorTools::MakeReal<
    typename TensorTools::DecrementRank<
    typename TensorTools::DecrementRank<
    OutputType>::type>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  const DenseSubVector<Number> &coef = get_localized_subvector(_system_vector, var);

  // Get finite element object
  FEGenericBase<OutputShape>* the_side_fe = NULL;
  this->get_side_fe<OutputShape>( var, the_side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &d2phi = the_side_fe->get_d2phi();

  // Loop over all the q_points in this finite element
  for (unsigned int qp=0; qp != d2u_vals.size(); qp++)
    {
      OutputType &d2u = d2u_vals[qp];

      // Compute the gradient at this q_point
      d2u = 0;

      for (unsigned int l=0; l != n_dofs; l++)
        d2u.add_scaled(d2phi[l][qp], coef(l));
    }

  return;
}



#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::point_value(unsigned int var, const Point &p) const
{
  Number u = 0.;

  this->point_value( var, p, u );

  return u;
}

template<typename OutputType>
void FEMContext::point_value(unsigned int var, const Point &p,
                             OutputType& u) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Build a FE for calculating u(p)
  AutoPtr<FEGenericBase<OutputShape> > fe_new = this->build_new_fe( fe, p );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<OutputShape> >&  phi = fe_new->get_phi();

  u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][0] * coef(l);

  return;
}



Gradient FEMContext::point_gradient(unsigned int var, const Point &p) const
{
  Gradient grad_u;

  this->point_gradient( var, p, grad_u );

  return grad_u;
}



template<typename OutputType>
void FEMContext::point_gradient(unsigned int var, const Point &p,
                                OutputType& grad_u) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Build a FE for calculating u(p)
  AutoPtr<FEGenericBase<OutputShape> > fe_new = this->build_new_fe( fe, p );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> >&  dphi = fe_new->get_dphi();

  grad_u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    grad_u.add_scaled(dphi[l][0], coef(l));

  return;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

Tensor FEMContext::point_hessian(unsigned int var, const Point &p) const
{
  Tensor hess_u;

  this->point_hessian( var, p, hess_u );

  return hess_u;
}


template<typename OutputType>
void FEMContext::point_hessian(unsigned int var, const Point &p,
                               OutputType& hess_u) const
{
  typedef typename TensorTools::MakeReal<
    typename TensorTools::DecrementRank<
    typename TensorTools::DecrementRank<
    OutputType>::type>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Build a FE for calculating u(p)
  AutoPtr<FEGenericBase<OutputShape> > fe_new = this->build_new_fe( fe, p );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> >&  d2phi = fe_new->get_d2phi();

  hess_u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    hess_u.add_scaled(d2phi[l][0], coef(l));

  return;
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES


template<typename OutputType>
void FEMContext::point_curl(unsigned int var, const Point &p,
                            OutputType& curl_u) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_subsolutions.size(), var);
  libmesh_assert(elem_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Build a FE for calculating u(p)
  AutoPtr<FEGenericBase<OutputShape> > fe_new = this->build_new_fe( fe, p );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputShape> >&  curl_phi = fe_new->get_curl_phi();

  curl_u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    curl_u.add_scaled(curl_phi[l][0], coef(l));

  return;
}



Number FEMContext::fixed_interior_value(unsigned int var, unsigned int qp) const
{
  Number u = 0.;

  this->fixed_interior_value( var, qp, u );

  return u;
}



template<typename OutputType>
void FEMContext::fixed_interior_value(unsigned int var, unsigned int qp,
                                      OutputType& u) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_fixed_subsolutions.size(), var);
  libmesh_assert(elem_fixed_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > &phi = fe->get_phi();

  // Accumulate solution value
  u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return;
}



Gradient FEMContext::fixed_interior_gradient(unsigned int var, unsigned int qp) const
{
  Gradient du;

  this->fixed_interior_gradient( var, qp, du );

  return du;
}


template<typename OutputType>
void FEMContext::FEMContext::fixed_interior_gradient(unsigned int var, unsigned int qp,
                                                     OutputType& du) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_fixed_subsolutions.size(), var);
  libmesh_assert(elem_fixed_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > &dphi = fe->get_dphi();

  // Accumulate solution derivatives
  du = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::fixed_interior_hessian(unsigned int var, unsigned int qp) const
{
  Tensor d2u;

  this->fixed_interior_hessian( var, qp, d2u );

  return d2u;
}


template<typename OutputType>
void FEMContext::fixed_interior_hessian(unsigned int var, unsigned int qp,
                                        OutputType& d2u) const
{
  typedef typename TensorTools::MakeReal<
    typename TensorTools::DecrementRank<
    typename TensorTools::DecrementRank<
    OutputType>::type>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_fixed_subsolutions.size(), var);
  libmesh_assert(elem_fixed_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &d2phi = fe->get_d2phi();

  // Accumulate solution second derivatives
  d2u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::fixed_side_value(unsigned int var, unsigned int qp) const
{
  Number u = 0.;

  this->fixed_side_value( var, qp, u );

  return u;
}


template<typename OutputType>
void FEMContext::fixed_side_value(unsigned int var, unsigned int qp,
                                  OutputType& u) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_fixed_subsolutions.size(), var);
  libmesh_assert(elem_fixed_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* the_side_fe = NULL;
  this->get_side_fe<OutputShape>( var, the_side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > &phi = the_side_fe->get_phi();

  // Accumulate solution value
  u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return;
}



Gradient FEMContext::fixed_side_gradient(unsigned int var, unsigned int qp) const
{
  Gradient du;

  this->fixed_side_gradient( var, qp, du );

  return du;
}


template<typename OutputType>
void FEMContext::fixed_side_gradient(unsigned int var, unsigned int qp,
                                     OutputType& du) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_fixed_subsolutions.size(), var);
  libmesh_assert(elem_fixed_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* the_side_fe = NULL;
  this->get_side_fe<OutputShape>( var, the_side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > &dphi = the_side_fe->get_dphi();

  // Accumulate solution derivatives
  du = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::fixed_side_hessian(unsigned int var, unsigned int qp) const
{
  Tensor d2u;

  this->fixed_side_hessian( var, qp, d2u );

  return d2u;
}

template<typename OutputType>
void FEMContext::fixed_side_hessian(unsigned int var, unsigned int qp,
                                    OutputType& d2u) const
{
  typedef typename TensorTools::MakeReal<
    typename TensorTools::DecrementRank<
    typename TensorTools::DecrementRank<
    OutputType>::type>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_fixed_subsolutions.size(), var);
  libmesh_assert(elem_fixed_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* the_side_fe = NULL;
  this->get_side_fe<OutputShape>( var, the_side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &d2phi = the_side_fe->get_d2phi();

  // Accumulate solution second derivatives
  d2u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::fixed_point_value(unsigned int var, const Point &p) const
{
  Number u = 0.;

  this->fixed_point_value( var, p, u );

  return u;
}

template<typename OutputType>
void FEMContext::fixed_point_value(unsigned int var, const Point &p,
                                   OutputType& u) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_fixed_subsolutions.size(), var);
  libmesh_assert(elem_fixed_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Build a FE for calculating u(p)
  AutoPtr<FEGenericBase<OutputShape> > fe_new = this->build_new_fe( fe, p );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<OutputShape> >&  phi = fe_new->get_phi();

  u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][0] * coef(l);

  return;
}



Gradient FEMContext::fixed_point_gradient(unsigned int var, const Point &p) const
{
  Gradient grad_u;

  this->fixed_point_gradient( var, p, grad_u );

  return grad_u;
}



template<typename OutputType>
void FEMContext::fixed_point_gradient(unsigned int var, const Point &p,
                                      OutputType& grad_u) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_fixed_subsolutions.size(), var);
  libmesh_assert(elem_fixed_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Build a FE for calculating u(p)
  AutoPtr<FEGenericBase<OutputShape> > fe_new = this->build_new_fe( fe, p );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> >&  dphi = fe_new->get_dphi();

  grad_u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    grad_u.add_scaled(dphi[l][0], coef(l));

  return;
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

Tensor FEMContext::fixed_point_hessian(unsigned int var, const Point &p) const
{
  Tensor hess_u;

  this->fixed_point_hessian( var, p, hess_u );

  return hess_u;
}



template<typename OutputType>
void FEMContext::fixed_point_hessian(unsigned int var, const Point &p,
                                     OutputType& hess_u) const
{
  typedef typename TensorTools::MakeReal<
    typename TensorTools::DecrementRank<
    typename TensorTools::DecrementRank<
    OutputType>::type>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  libmesh_assert_greater (dof_indices.size(), var);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices_var[var].size());

  // Get current local coefficients
  libmesh_assert_greater (elem_fixed_subsolutions.size(), var);
  libmesh_assert(elem_fixed_subsolutions[var]);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Build a FE for calculating u(p)
  AutoPtr<FEGenericBase<OutputShape> > fe_new = this->build_new_fe( fe, p );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> >&  d2phi = fe_new->get_d2phi();

  hess_u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    hess_u.add_scaled(d2phi[l][0], coef(l));

  return;
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES





void FEMContext::elem_reinit(Real theta)
{
  // Update the "time" variable of this context object
  this->_update_time_from_system(theta);

  // Handle a moving element if necessary.
  if (_mesh_sys)
    // We assume that the ``default'' state
    // of the mesh is its final, theta=1.0
    // position, so we don't bother with
    // mesh motion in that case.
    if (theta != 1.0)
      {
        // FIXME - ALE is not threadsafe yet!
        libmesh_assert_equal_to (libMesh::n_threads(), 1);

        elem_position_set(theta);
      }
  elem_fe_reinit();
}


void FEMContext::elem_side_reinit(Real theta)
{
  // Update the "time" variable of this context object
  this->_update_time_from_system(theta);

  // Handle a moving element if necessary
  if (_mesh_sys)
    {
      // FIXME - not threadsafe yet!
      elem_position_set(theta);
      side_fe_reinit();
    }
}


void FEMContext::elem_edge_reinit(Real theta)
{
  // Update the "time" variable of this context object
  this->_update_time_from_system(theta);

  // Handle a moving element if necessary
  if (_mesh_sys)
    {
      // FIXME - not threadsafe yet!
      elem_position_set(theta);
      edge_fe_reinit();
    }
}


void FEMContext::nonlocal_reinit(Real theta)
{
  // Update the "time" variable of this context object
  this->_update_time_from_system(theta);

  // We can reuse the Elem FE safely here.
  elem_fe_reinit();
}


void FEMContext::elem_fe_reinit ()
{
  // Initialize all the interior FE objects on elem.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEAbstract *>::iterator local_fe_end = _element_fe.end();
  for (std::map<FEType, FEAbstract *>::iterator i = _element_fe.begin();
       i != local_fe_end; ++i)
    {
      i->second->reinit(elem);
    }
}


void FEMContext::side_fe_reinit ()
{
  // Initialize all the side FE objects on elem/side.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEAbstract *>::iterator local_fe_end = _side_fe.end();
  for (std::map<FEType, FEAbstract *>::iterator i = _side_fe.begin();
       i != local_fe_end; ++i)
    {
      i->second->reinit(elem, side);
    }
}



void FEMContext::edge_fe_reinit ()
{
  libmesh_assert_equal_to (dim, 3);

  // Initialize all the interior FE objects on elem/edge.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEAbstract *>::iterator local_fe_end = _edge_fe.end();
  for (std::map<FEType, FEAbstract *>::iterator i = _edge_fe.begin();
       i != local_fe_end; ++i)
    {
      i->second->edge_reinit(elem, edge);
    }
}



void FEMContext::elem_position_get()
{
  // This is too expensive to call unless we've been asked to move the mesh
  libmesh_assert (_mesh_sys);

  // This will probably break with threading when two contexts are
  // operating on elements which share a node
  libmesh_assert_equal_to (libMesh::n_threads(), 1);

  // If the coordinate data is in our own system, it's already
  // been set up for us
  //  if (_mesh_sys == this->number())
  //    {
  unsigned int n_nodes = elem->n_nodes();
  // For simplicity we demand that mesh coordinates be stored
  // in a format that allows a direct copy
  libmesh_assert(_mesh_x_var == libMesh::invalid_uint ||
                 (_element_fe_var[_mesh_x_var]->get_fe_type().family
                  == LAGRANGE &&
                  _element_fe_var[_mesh_x_var]->get_fe_type().order
                  == elem->default_order()));
  libmesh_assert(_mesh_y_var == libMesh::invalid_uint ||
                 (_element_fe_var[_mesh_y_var]->get_fe_type().family
                  == LAGRANGE &&
                  _element_fe_var[_mesh_y_var]->get_fe_type().order
                  == elem->default_order()));
  libmesh_assert(_mesh_z_var == libMesh::invalid_uint ||
                 (_element_fe_var[_mesh_z_var]->get_fe_type().family
                  == LAGRANGE &&
                  _element_fe_var[_mesh_z_var]->get_fe_type().order
                  == elem->default_order()));

  // Get degree of freedom coefficients from point coordinates
  if (_mesh_x_var != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      (*elem_subsolutions[_mesh_x_var])(i) = elem->point(i)(0);

  if (_mesh_y_var != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      (*elem_subsolutions[_mesh_y_var])(i) = elem->point(i)(1);

  if (_mesh_z_var != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      (*elem_subsolutions[_mesh_z_var])(i) = elem->point(i)(2);
  //    }
  // FIXME - If the coordinate data is not in our own system, someone
  // had better get around to implementing that... - RHS
  //  else
  //    {
  //      libmesh_not_implemented();
  //    }
}



// We can ignore the theta argument in the current use of this
// function, because elem_subsolutions will already have been set to
// the theta value.
//
// To enable loose mesh movement coupling things will need to change.
void FEMContext::_do_elem_position_set(Real)
{
  // This is too expensive to call unless we've been asked to move the mesh
  libmesh_assert (_mesh_sys);

  // This will probably break with threading when two contexts are
  // operating on elements which share a node
  libmesh_assert_equal_to (libMesh::n_threads(), 1);

  // If the coordinate data is in our own system, it's already
  // been set up for us, and we can ignore our input parameter theta
  //  if (_mesh_sys == this->number())
  //    {
  unsigned int n_nodes = elem->n_nodes();
  // For simplicity we demand that mesh coordinates be stored
  // in a format that allows a direct copy
  libmesh_assert(_mesh_x_var == libMesh::invalid_uint ||
                 (_element_fe_var[_mesh_x_var]->get_fe_type().family
                  == LAGRANGE &&
                  elem_subsolutions[_mesh_x_var]->size() == n_nodes));
  libmesh_assert(_mesh_y_var == libMesh::invalid_uint ||
                 (_element_fe_var[_mesh_y_var]->get_fe_type().family
                  == LAGRANGE &&
                  elem_subsolutions[_mesh_y_var]->size() == n_nodes));
  libmesh_assert(_mesh_z_var == libMesh::invalid_uint ||
                 (_element_fe_var[_mesh_z_var]->get_fe_type().family
                  == LAGRANGE &&
                  elem_subsolutions[_mesh_z_var]->size() == n_nodes));

  // Set the new point coordinates
  if (_mesh_x_var != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      const_cast<Elem*>(elem)->point(i)(0) =
        libmesh_real((*elem_subsolutions[_mesh_x_var])(i));

  if (_mesh_y_var != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      const_cast<Elem*>(elem)->point(i)(1) =
        libmesh_real((*elem_subsolutions[_mesh_y_var])(i));

  if (_mesh_z_var != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      const_cast<Elem*>(elem)->point(i)(2) =
        libmesh_real((*elem_subsolutions[_mesh_z_var])(i));
  //    }
  // FIXME - If the coordinate data is not in our own system, someone
  // had better get around to implementing that... - RHS
  //  else
  //    {
  //      libmesh_not_implemented();
  //    }
}





/*
  void FEMContext::reinit(const FEMSystem &sys, Elem *e)
  {
  // Initialize our elem pointer, algebraic objects
  this->pre_fe_reinit(e);

  // Moving the mesh may be necessary
  // Reinitializing the FE objects is definitely necessary
  this->elem_reinit(1.);
  }
*/



void FEMContext::pre_fe_reinit(const System &sys, const Elem *e)
{
  elem = e;

  // Initialize the per-element data for elem.
  sys.get_dof_map().dof_indices (elem, dof_indices);
  const unsigned int n_dofs = cast_int<unsigned int>
    (dof_indices.size());
  const std::size_t n_qoi = sys.qoi.size();

  if (sys.use_fixed_solution)
    elem_fixed_solution.resize(n_dofs);

  // This also resizes elem_solution
  sys.current_local_solution->get(dof_indices, elem_solution.get_values());

  // These resize calls also zero out the residual and jacobian
  elem_residual.resize(n_dofs);
  elem_jacobian.resize(n_dofs, n_dofs);

  elem_qoi_derivative.resize(n_qoi);
  elem_qoi_subderivatives.resize(n_qoi);
  for (std::size_t q=0; q != n_qoi; ++q)
    elem_qoi_derivative[q].resize(n_dofs);

  // Initialize the per-variable data for elem.
  {
    unsigned int sub_dofs = 0;
    for (unsigned int i=0; i != sys.n_vars(); ++i)
      {
        sys.get_dof_map().dof_indices (elem, dof_indices_var[i], i);

        const unsigned int n_dofs_var = cast_int<unsigned int>
          (dof_indices_var[i].size());

        elem_subsolutions[i]->reposition
          (sub_dofs, n_dofs_var);

        if (sys.use_fixed_solution)
          elem_fixed_subsolutions[i]->reposition
            (sub_dofs, n_dofs_var);

        elem_subresiduals[i]->reposition
          (sub_dofs, n_dofs_var);

        for (std::size_t q=0; q != n_qoi; ++q)
          elem_qoi_subderivatives[q][i]->reposition
            (sub_dofs, n_dofs_var);

        for (unsigned int j=0; j != i; ++j)
          {
            const unsigned int n_dofs_var_j =
              cast_int<unsigned int>
              (dof_indices_var[j].size());

            elem_subjacobians[i][j]->reposition
              (sub_dofs, elem_subresiduals[j]->i_off(),
               n_dofs_var, n_dofs_var_j);
            elem_subjacobians[j][i]->reposition
              (elem_subresiduals[j]->i_off(), sub_dofs,
               n_dofs_var_j, n_dofs_var);
          }
        elem_subjacobians[i][i]->reposition
          (sub_dofs, sub_dofs,
           n_dofs_var,
           n_dofs_var);
        sub_dofs += n_dofs_var;
      }
    libmesh_assert_equal_to (sub_dofs, n_dofs);
  }

  // Now do the localization for the user requested vectors
  DiffContext::localized_vectors_iterator localized_vec_it = localized_vectors.begin();
  const DiffContext::localized_vectors_iterator localized_vec_end = localized_vectors.end();

  for(; localized_vec_it != localized_vec_end; ++localized_vec_it)
    {
      const NumericVector<Number>& current_localized_vector = *localized_vec_it->first;
      DenseVector<Number>& target_vector = localized_vec_it->second.first;

      current_localized_vector.get(dof_indices, target_vector.get_values());

      // Initialize the per-variable data for elem.
      unsigned int sub_dofs = 0;
      for (unsigned int i=0; i != sys.n_vars(); ++i)
        {
          const unsigned int n_dofs_var = cast_int<unsigned int>
            (dof_indices_var[i].size());
          sys.get_dof_map().dof_indices (elem, dof_indices_var[i], i);

          localized_vec_it->second.second[i]->reposition
            (sub_dofs, n_dofs_var);

          sub_dofs += n_dofs_var;
        }
      libmesh_assert_equal_to (sub_dofs, n_dofs);
    }
}



void FEMContext::_update_time_from_system(Real theta)
{
  // Update the "time" variable based on the value of theta.  For this
  // to work, we need to know the value of deltat, a pointer to which is now
  // stored by our parent DiffContext class.  Note: get_deltat_value() will
  // assert in debug mode if the requested pointer is NULL.
  const Real deltat = this->get_deltat_value();

  this->time = theta*(this->system_time + deltat) + (1.-theta)*this->system_time;
}



template<typename OutputShape>
AutoPtr<FEGenericBase<OutputShape> > FEMContext::build_new_fe( const FEGenericBase<OutputShape>* fe,
                                                               const Point &p ) const
{
  FEType fe_type = fe->get_fe_type();

  // If we don't have an Elem to evaluate on, then the only functions
  // we can sensibly evaluate are the scalar dofs which are the same
  // everywhere.
  libmesh_assert(elem || fe_type.family == SCALAR);

  unsigned int elem_dim = elem ? elem->dim() : 0;

  FEGenericBase<OutputShape>* fe_new =
    FEGenericBase<OutputShape>::build(dim, fe_type).release();

  // Map the physical co-ordinates to the master co-ordinates using the inverse_map from fe_interface.h
  // Build a vector of point co-ordinates to send to reinit
  Point master_point = elem ?
    FEInterface::inverse_map(elem_dim, fe_type, elem, p) :
    Point(0);

  std::vector<Point> coor(1, master_point);

  // Reinitialize the element and compute the shape function values at coor
  fe_new->reinit (elem, &coor);

  return AutoPtr<FEGenericBase<OutputShape> >(fe_new);
}





// Instantiate member function templates
template void FEMContext::interior_value<Number>(unsigned int, unsigned int, Number&) const;
template void FEMContext::interior_values<Number>(unsigned int, const NumericVector<Number> &,
                                                  std::vector<Number>&) const;
template void FEMContext::interior_value<Gradient>(unsigned int, unsigned int, Gradient&) const;
template void FEMContext::interior_values<Gradient>(unsigned int, const NumericVector<Number> &,
                                                    std::vector<Gradient>&) const;

template void FEMContext::interior_gradient<Gradient>(unsigned int, unsigned int, Gradient&) const;
template void FEMContext::interior_gradients<Gradient>(unsigned int, const NumericVector<Number> &,
                                                       std::vector<Gradient>&) const;
template void FEMContext::interior_gradient<Tensor>(unsigned int, unsigned int, Tensor&) const;
template void FEMContext::interior_gradients<Tensor>(unsigned int, const NumericVector<Number> &,
                                                     std::vector<Tensor>&) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::interior_hessian<Tensor>(unsigned int, unsigned int, Tensor&) const;
template void FEMContext::interior_hessians<Tensor>(unsigned int, const NumericVector<Number> &,
                                                    std::vector<Tensor>&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::interior_hessian<??>(unsigned int, unsigned int, ??&) const;
//template void FEMContext::interior_hessians<??>(unsigned int, const NumericVector<Number> &,
//                                                std::vector<??>&) const;
#endif

template void FEMContext::interior_curl<Gradient>(unsigned int, unsigned int, Gradient&) const;

template void FEMContext::interior_div<Number>(unsigned int, unsigned int, Number&) const;

template void FEMContext::side_value<Number>(unsigned int, unsigned int, Number&) const;
template void FEMContext::side_value<Gradient>(unsigned int, unsigned int, Gradient&) const;
template void FEMContext::side_values<Number>(unsigned int, const NumericVector<Number> &,
                                              std::vector<Number>&) const;
template void FEMContext::side_values<Gradient>(unsigned int, const NumericVector<Number> &,
                                                std::vector<Gradient>&) const;

template void FEMContext::side_gradient<Gradient>(unsigned int, unsigned int, Gradient&) const;
template void FEMContext::side_gradients<Gradient>(unsigned int, const NumericVector<Number> &,
                                                   std::vector<Gradient>&) const;
template void FEMContext::side_gradient<Tensor>(unsigned int, unsigned int, Tensor&) const;
template void FEMContext::side_gradients<Tensor>(unsigned int, const NumericVector<Number> &,
                                                 std::vector<Tensor>&) const;


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::side_hessian<Tensor>(unsigned int, unsigned int, Tensor&) const;
template void FEMContext::side_hessians<Tensor>(unsigned int, const NumericVector<Number> &,
                                                std::vector<Tensor>&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::side_hessian<??>(unsigned int, unsigned int,
//                                           ??&) const;
//template void FEMContext::side_hessians<??>(unsigned int, const NumericVector<Number> &,
//                                            std::vector<??>&) const;
#endif

template void FEMContext::point_value<Number>(unsigned int, const Point&, Number&) const;
template void FEMContext::point_value<Gradient>(unsigned int, const Point&, Gradient&) const;

template void FEMContext::point_gradient<Gradient>(unsigned int, const Point&, Gradient&) const;
template void FEMContext::point_gradient<Tensor>(unsigned int, const Point&, Tensor&) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::point_hessian<Tensor>(unsigned int, const Point&, Tensor&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::point_hessian<??>(unsigned int, const Point&, ??&) const;
#endif

template void FEMContext::point_curl<Gradient>(unsigned int, const Point&, Gradient&) const;

template void FEMContext::fixed_interior_value<Number>(unsigned int, unsigned int, Number&) const;
template void FEMContext::fixed_interior_value<Gradient>(unsigned int, unsigned int, Gradient&) const;

template void FEMContext::fixed_interior_gradient<Gradient>(unsigned int, unsigned int, Gradient&) const;
template void FEMContext::fixed_interior_gradient<Tensor>(unsigned int, unsigned int, Tensor&) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::fixed_interior_hessian<Tensor>(unsigned int, unsigned int, Tensor&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::fixed_interior_hessian<??>(unsigned int, unsigned int, ??&) const;
#endif

template void FEMContext::fixed_side_value<Number>(unsigned int, unsigned int, Number&) const;
template void FEMContext::fixed_side_value<Gradient>(unsigned int, unsigned int, Gradient&) const;

template void FEMContext::fixed_side_gradient<Gradient>(unsigned int, unsigned int, Gradient&) const;
template void FEMContext::fixed_side_gradient<Tensor>(unsigned int, unsigned int, Tensor&) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::fixed_side_hessian<Tensor>(unsigned int, unsigned int, Tensor&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::fixed_side_hessian<??>(unsigned int, unsigned int, ??&) const;
#endif

template void FEMContext::fixed_point_value<Number>(unsigned int, const Point&, Number&) const;
template void FEMContext::fixed_point_value<Gradient>(unsigned int, const Point&, Gradient&) const;

template void FEMContext::fixed_point_gradient<Gradient>(unsigned int, const Point&, Gradient&) const;
template void FEMContext::fixed_point_gradient<Tensor>(unsigned int, const Point&, Tensor&) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::fixed_point_hessian<Tensor>(unsigned int, const Point&, Tensor&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::fixed_point_hessian<??>(unsigned int, const Point&, ??&) const;
#endif

template AutoPtr<FEGenericBase<Real> > FEMContext::build_new_fe( const FEGenericBase<Real>*, const Point & ) const;
template AutoPtr<FEGenericBase<RealGradient> > FEMContext::build_new_fe( const FEGenericBase<RealGradient>*, const Point & ) const;


} // namespace libMesh
