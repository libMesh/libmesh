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
#include "libmesh/diff_system.h"
#include "libmesh/time_solver.h"
#include "libmesh/unsteady_solver.h" // For euler_residual

namespace libMesh
{

FEMContext::FEMContext (const System & sys)
  : DiffContext(sys),
    _mesh_sys(libmesh_nullptr),
    _mesh_x_var(0),
    _mesh_y_var(0),
    _mesh_z_var(0),
    side(0), edge(0),
    _atype(CURRENT),
    _custom_solution(libmesh_nullptr),
    _boundary_info(sys.get_mesh().get_boundary_info()),
    _elem(libmesh_nullptr),
    _dim(sys.get_mesh().mesh_dimension()),
    _elem_dim(0), /* This will be reset in set_elem(). */
    _elem_dims(sys.get_mesh().elem_dimensions()),
    _element_qrule(4,libmesh_nullptr),
    _side_qrule(4,libmesh_nullptr)
{
  // Reserve space for the FEAbstract and QBase objects for each
  // element dimension possiblity (0,1,2,3)
  _element_fe.resize(4);
  _side_fe.resize(4);
  _element_fe_var.resize(4);
  _side_fe_var.resize(4);

  // We need to know which of our variables has the hardest
  // shape functions to numerically integrate.

  unsigned int nv = sys.n_vars();

  libmesh_assert (nv);
  FEType hardest_fe_type = sys.variable_type(0);

  bool have_scalar = false;

  for (unsigned int i=0; i != nv; ++i)
    {
      FEType fe_type = sys.variable_type(i);

      // Make sure we find a non-SCALAR FE family, even in the case
      // where the first variable(s) weren't
      if (hardest_fe_type.family == SCALAR)
        {
          hardest_fe_type.family = fe_type.family;
          hardest_fe_type.order = fe_type.order;
        }

      // FIXME - we don't yet handle mixed finite elements from
      // different families which require different quadrature rules
      // libmesh_assert_equal_to (fe_type.family, hardest_fe_type.family);

      // We need to detect SCALAR's so we can prepare FE objects for
      // them, and so we don't mistake high order scalars as a reason
      // to crank up the quadrature order on other types.
      if( fe_type.family == SCALAR )
        have_scalar = true;
      else if (fe_type.order > hardest_fe_type.order)
        hardest_fe_type = fe_type;
    }

  if(have_scalar)
    // SCALAR FEs have dimension 0 by assumption
    _elem_dims.insert(0);

  for( std::set<unsigned char>::const_iterator dim_it = _elem_dims.begin();
       dim_it != _elem_dims.end(); ++dim_it )
    {
      const unsigned char dim = *dim_it;

      // Create an adequate quadrature rule
      _element_qrule[dim] = hardest_fe_type.default_quadrature_rule
        (dim, sys.extra_quadrature_order).release();
      _side_qrule[dim] = hardest_fe_type.default_quadrature_rule
        (dim-1, sys.extra_quadrature_order).release();
      if (dim == 3)
        _edge_qrule.reset(hardest_fe_type.default_quadrature_rule
                          (1, sys.extra_quadrature_order).release());

      // Next, create finite element objects
      _element_fe_var[dim].resize(nv);
      _side_fe_var[dim].resize(nv);
      if (dim == 3)
        _edge_fe_var.resize(nv);


      for (unsigned int i=0; i != nv; ++i)
        {
          FEType fe_type = sys.variable_type(i);

          if ( _element_fe[dim][fe_type] == libmesh_nullptr )
            {
              _element_fe[dim][fe_type] = FEAbstract::build(dim, fe_type).release();
              _element_fe[dim][fe_type]->attach_quadrature_rule(_element_qrule[dim]);
              _side_fe[dim][fe_type] = FEAbstract::build(dim, fe_type).release();
              _side_fe[dim][fe_type]->attach_quadrature_rule(_side_qrule[dim]);

              if (dim == 3)
                {
                  _edge_fe[fe_type] = FEAbstract::build(dim, fe_type).release();
                  _edge_fe[fe_type]->attach_quadrature_rule(_edge_qrule.get());
                }
            }

          _element_fe_var[dim][i] = _element_fe[dim][fe_type];
          _side_fe_var[dim][i] = _side_fe[dim][fe_type];
          if ((dim) == 3)
            _edge_fe_var[i] = _edge_fe[fe_type];

        }
    }
}

FEMContext::~FEMContext()
{
  // We don't want to store UniquePtrs in STL containers, but we don't
  // want to leak memory either
  for (std::vector<std::map<FEType, FEAbstract *> >::iterator d = _element_fe.begin();
       d != _element_fe.end(); ++d)
    for (std::map<FEType, FEAbstract *>::iterator i = d->begin();
         i != d->end(); ++i)
      delete i->second;

  for (std::vector<std::map<FEType, FEAbstract *> >::iterator d = _side_fe.begin();
       d != _side_fe.end(); ++d)
    for (std::map<FEType, FEAbstract *>::iterator i = d->begin();
         i != d->end(); ++i)
      delete i->second;

  for (std::map<FEType, FEAbstract *>::iterator i = _edge_fe.begin();
       i != _edge_fe.end(); ++i)
    delete i->second;
  _edge_fe.clear();

  for (std::vector<QBase *>::iterator i = _element_qrule.begin();
       i != _element_qrule.end(); ++i)
    delete *i;
  _element_qrule.clear();

  for (std::vector<QBase *>::iterator i = _side_qrule.begin();
       i != _side_qrule.end(); ++i)
    delete *i;
  _side_qrule.clear();
}



bool FEMContext::has_side_boundary_id(boundary_id_type id) const
{
  return _boundary_info.has_boundary_id(&(this->get_elem()), side, id);
}


std::vector<boundary_id_type> FEMContext::side_boundary_ids() const
{
  libmesh_deprecated();
  return _boundary_info.boundary_ids(&(this->get_elem()), side);
}


void FEMContext::side_boundary_ids(std::vector<boundary_id_type> & vec_to_fill) const
{
  _boundary_info.boundary_ids(&(this->get_elem()), side, vec_to_fill);
}



template<typename OutputType,
         typename FEMContext::FENeeded<OutputType>::value_getter fe_getter,
         FEMContext::diff_subsolution_getter subsolution_getter>
void FEMContext::some_value(unsigned int var, unsigned int qp, OutputType & u) const
{
  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  const DenseSubVector<Number> & coef = (this->*subsolution_getter)(var);

  // Get finite element object
  typename FENeeded<OutputType>::value_base * fe = libmesh_nullptr;
  (this->*fe_getter)( var, fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector
                    <typename FENeeded<OutputType>::value_shape> > & phi = fe->get_phi();

  // Accumulate solution value
  u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);
}



template<typename OutputType,
         typename FEMContext::FENeeded<OutputType>::grad_getter fe_getter,
         FEMContext::diff_subsolution_getter subsolution_getter>
void FEMContext::some_gradient(unsigned int var, unsigned int qp, OutputType & du) const
{
  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  const DenseSubVector<Number> & coef = (this->*subsolution_getter)(var);

  // Get finite element object
  typename FENeeded<OutputType>::grad_base * fe = libmesh_nullptr;
  (this->*fe_getter)( var, fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector
                    <typename FENeeded<OutputType>::grad_base::OutputGradient> >
    & dphi = fe->get_dphi();

  // Accumulate solution derivatives
  du = 0;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template<typename OutputType,
         typename FEMContext::FENeeded<OutputType>::hess_getter fe_getter,
         FEMContext::diff_subsolution_getter subsolution_getter>
void FEMContext::some_hessian(unsigned int var, unsigned int qp, OutputType & d2u) const
{
  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  const DenseSubVector<Number> & coef = (this->*subsolution_getter)(var);

  // Get finite element object
  typename FENeeded<OutputType>::hess_base * fe = libmesh_nullptr;
  (this->*fe_getter)( var, fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector
                    <typename FENeeded<OutputType>::hess_base::OutputTensor> >
    & d2phi = fe->get_d2phi();

  // Accumulate solution second derivatives
  d2u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return;
}
#endif



Number FEMContext::interior_value(unsigned int var, unsigned int qp) const
{
  Number u;

  this->interior_value( var, qp, u );

  return u;
}

template<typename OutputType>
void FEMContext::interior_value(unsigned int var, unsigned int qp,
                                OutputType & u) const
{
  this->some_value<OutputType,
                   &FEMContext::get_element_fe<typename TensorTools::MakeReal<OutputType>::type>,
                   &DiffContext::get_elem_solution>(var, qp, u);
}


template<typename OutputType>
void FEMContext::interior_values (unsigned int var,
                                  const NumericVector<Number> & _system_vector,
                                  std::vector<OutputType> & u_vals) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  const DenseSubVector<Number> & coef = get_localized_subvector(_system_vector, var);

  // Get the finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > & phi = fe->get_phi();

  // Loop over all the q_points on this element
  for (std::size_t qp=0; qp != u_vals.size(); qp++)
    {
      OutputType & u = u_vals[qp];

      // Compute the value at this q_point
      u = 0.;

      for (unsigned int l=0; l != n_dofs; l++)
        u += phi[l][qp] * coef(l);
    }

  return;
}

Gradient FEMContext::interior_gradient(unsigned int var,
                                       unsigned int qp) const
{
  Gradient du;

  this->interior_gradient( var, qp, du );

  return du;
}



template<typename OutputType>
void FEMContext::interior_gradient(unsigned int var,
                                   unsigned int qp,
                                   OutputType & du) const
{
  this->some_gradient<OutputType,
                      &FEMContext::get_element_fe<typename TensorTools::MakeReal
                                                  <typename TensorTools::DecrementRank
                                                   <OutputType>::type>::type>,
                      &DiffContext::get_elem_solution>(var, qp, du);
}



template<typename OutputType>
void FEMContext::interior_gradients(unsigned int var,
                                    const NumericVector<Number> & _system_vector,
                                    std::vector<OutputType> & du_vals) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  const DenseSubVector<Number> & coef = get_localized_subvector(_system_vector, var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > & dphi = fe->get_dphi();

  // Loop over all the q_points in this finite element
  for (std::size_t qp=0; qp != du_vals.size(); qp++)
    {
      OutputType & du = du_vals[qp];

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
                                  OutputType & d2u) const
{
  this->some_hessian<OutputType,
                     &FEMContext::get_element_fe
                     <typename TensorTools::MakeReal
                      <typename TensorTools::DecrementRank
                       <typename TensorTools::DecrementRank
                        <OutputType>::type>::type>::type>,
                     &DiffContext::get_elem_solution>(var, qp, d2u);
}


template<typename OutputType>
void FEMContext::interior_hessians(unsigned int var,
                                   const NumericVector<Number> & _system_vector,
                                   std::vector<OutputType> & d2u_vals) const
{
  typedef typename TensorTools::DecrementRank<OutputType>::type Rank1Decrement;
  typedef typename TensorTools::DecrementRank<Rank1Decrement>::type Rank2Decrement;
  typedef typename TensorTools::MakeReal<Rank2Decrement>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  const DenseSubVector<Number> & coef = get_localized_subvector(_system_vector, var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > & d2phi = fe->get_d2phi();

  // Loop over all the q_points in this finite element
  for (std::size_t qp=0; qp != d2u_vals.size(); qp++)
    {
      OutputType & d2u = d2u_vals[qp];

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
                               OutputType & curl_u) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (this->_elem_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputShape> > & curl_phi = fe->get_curl_phi();

  // Accumulate solution curl
  curl_u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    curl_u.add_scaled(curl_phi[l][qp], coef(l));

  return;
}


template<typename OutputType>
void FEMContext::interior_div(unsigned int var, unsigned int qp,
                              OutputType & div_u) const
{
  typedef typename
    TensorTools::IncrementRank
    <typename TensorTools::MakeReal<OutputType>::type>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (this->_elem_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputDivergence> > & div_phi = fe->get_div_phi();

  // Accumulate solution curl
  div_u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    div_u += div_phi[l][qp] * coef(l);

  return;
}


Number FEMContext::side_value(unsigned int var,
                              unsigned int qp) const
{
  Number u = 0.;

  this->side_value( var, qp, u );

  return u;
}


template<typename OutputType>
void FEMContext::side_value(unsigned int var,
                            unsigned int qp,
                            OutputType & u) const
{
  this->some_value<OutputType,
                   &FEMContext::get_side_fe<typename TensorTools::MakeReal<OutputType>::type>,
                   &DiffContext::get_elem_solution>(var, qp, u);
}


template<typename OutputType>
void FEMContext::side_values(unsigned int var,
                             const NumericVector<Number> & _system_vector,
                             std::vector<OutputType> & u_vals) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  const DenseSubVector<Number> & coef = get_localized_subvector(_system_vector, var);

  // Get the finite element object
  FEGenericBase<OutputShape> * the_side_fe = libmesh_nullptr;
  this->get_side_fe<OutputShape>( var, the_side_fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > & phi = the_side_fe->get_phi();

  // Loop over all the q_points on this element
  for (std::size_t qp=0; qp != u_vals.size(); qp++)
    {
      OutputType & u = u_vals[qp];

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
                               OutputType & du) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (this->_elem_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * the_side_fe = libmesh_nullptr;
  this->get_side_fe<OutputShape>( var, the_side_fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector< typename FEGenericBase<OutputShape>::OutputGradient> > & dphi = the_side_fe->get_dphi();

  // Accumulate solution derivatives
  du = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return;
}



template<typename OutputType>
void FEMContext::side_gradients(unsigned int var,
                                const NumericVector<Number> & _system_vector,
                                std::vector<OutputType> & du_vals) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  const DenseSubVector<Number> & coef = get_localized_subvector(_system_vector, var);

  // Get finite element object
  FEGenericBase<OutputShape> * the_side_fe = libmesh_nullptr;
  this->get_side_fe<OutputShape>( var, the_side_fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > & dphi = the_side_fe->get_dphi();

  // Loop over all the q_points in this finite element
  for (std::size_t qp=0; qp != du_vals.size(); qp++)
    {
      OutputType & du = du_vals[qp];

      du = 0;

      // Compute the gradient at this q_point
      for (unsigned int l=0; l != n_dofs; l++)
        du.add_scaled(dphi[l][qp], coef(l));
    }

  return;
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::side_hessian(unsigned int var,
                                unsigned int qp) const
{
  Tensor d2u;

  this->side_hessian( var, qp, d2u );

  return d2u;
}



template<typename OutputType>
void FEMContext::side_hessian(unsigned int var,
                              unsigned int qp,
                              OutputType & d2u) const
{
  this->some_hessian<OutputType,
                     &FEMContext::get_side_fe
                     <typename TensorTools::MakeReal
                      <typename TensorTools::DecrementRank
                       <typename TensorTools::DecrementRank
                        <OutputType>::type>::type>::type>,
                     &DiffContext::get_elem_solution>(var, qp, d2u);
}



template<typename OutputType>
void FEMContext::side_hessians(unsigned int var,
                               const NumericVector<Number> & _system_vector,
                               std::vector<OutputType> & d2u_vals) const
{
  typedef typename TensorTools::DecrementRank<OutputType>::type Rank1Decrement;
  typedef typename TensorTools::DecrementRank<Rank1Decrement>::type Rank2Decrement;
  typedef typename TensorTools::MakeReal<Rank2Decrement>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  const DenseSubVector<Number> & coef = get_localized_subvector(_system_vector, var);

  // Get finite element object
  FEGenericBase<OutputShape> * the_side_fe = libmesh_nullptr;
  this->get_side_fe<OutputShape>( var, the_side_fe, this->get_elem_dim() );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > & d2phi = the_side_fe->get_d2phi();

  // Loop over all the q_points in this finite element
  for (std::size_t qp=0; qp != d2u_vals.size(); qp++)
    {
      OutputType & d2u = d2u_vals[qp];

      // Compute the gradient at this q_point
      d2u = 0;

      for (unsigned int l=0; l != n_dofs; l++)
        d2u.add_scaled(d2phi[l][qp], coef(l));
    }

  return;
}



#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::point_value(unsigned int var, const Point & p) const
{
  Number u = 0.;

  this->point_value( var, p, u );

  return u;
}

template<typename OutputType>
void FEMContext::point_value(unsigned int var,
                             const Point & p,
                             OutputType & u,
                             const Real tolerance) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (this->_elem_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Build a FE for calculating u(p)
  UniquePtr<FEGenericBase<OutputShape> > fe_new =
    this->build_new_fe( fe, p, tolerance );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<OutputShape> > &  phi = fe_new->get_phi();

  u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][0] * coef(l);

  return;
}



Gradient FEMContext::point_gradient(unsigned int var, const Point & p) const
{
  Gradient grad_u;

  this->point_gradient( var, p, grad_u );

  return grad_u;
}



template<typename OutputType>
void FEMContext::point_gradient(unsigned int var,
                                const Point & p,
                                OutputType & grad_u,
                                const Real tolerance) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (this->_elem_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Build a FE for calculating u(p)
  UniquePtr<FEGenericBase<OutputShape> > fe_new =
    this->build_new_fe( fe, p, tolerance );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > &  dphi = fe_new->get_dphi();

  grad_u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    grad_u.add_scaled(dphi[l][0], coef(l));

  return;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

Tensor FEMContext::point_hessian(unsigned int var, const Point & p) const
{
  Tensor hess_u;

  this->point_hessian( var, p, hess_u );

  return hess_u;
}


template<typename OutputType>
void FEMContext::point_hessian(unsigned int var,
                               const Point & p,
                               OutputType & hess_u,
                               const Real tolerance) const
{
  typedef typename TensorTools::DecrementRank<OutputType>::type Rank1Decrement;
  typedef typename TensorTools::DecrementRank<Rank1Decrement>::type Rank2Decrement;
  typedef typename TensorTools::MakeReal<Rank2Decrement>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (this->_elem_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Build a FE for calculating u(p)
  UniquePtr<FEGenericBase<OutputShape> > fe_new =
    this->build_new_fe( fe, p, tolerance );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &  d2phi = fe_new->get_d2phi();

  hess_u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    hess_u.add_scaled(d2phi[l][0], coef(l));

  return;
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES


template<typename OutputType>
void FEMContext::point_curl(unsigned int var,
                            const Point & p,
                            OutputType & curl_u,
                            const Real tolerance) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (this->_elem_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Build a FE for calculating u(p)
  UniquePtr<FEGenericBase<OutputShape> > fe_new =
    this->build_new_fe( fe, p, tolerance );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputShape> > &  curl_phi = fe_new->get_curl_phi();

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
                                      OutputType & u) const
{
  this->some_value<OutputType,
                   &FEMContext::get_element_fe
                   <typename TensorTools::MakeReal<OutputType>::type>,
                   &DiffContext::get_elem_fixed_solution>(var, qp, u);
}



Gradient FEMContext::fixed_interior_gradient(unsigned int var, unsigned int qp) const
{
  Gradient du;

  this->fixed_interior_gradient( var, qp, du );

  return du;
}


template<typename OutputType>
void FEMContext::fixed_interior_gradient(unsigned int var, unsigned int qp,
                                         OutputType & du) const
{
  this->some_gradient
    <OutputType,
     &FEMContext::get_element_fe
     <typename TensorTools::MakeReal
      <typename TensorTools::DecrementRank
       <OutputType>::type>::type>,
     &DiffContext::get_elem_fixed_solution>
    (var, qp, du);
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
                                        OutputType & d2u) const
{
  this->some_hessian<OutputType,
                     &FEMContext::get_element_fe
                     <typename TensorTools::MakeReal
                      <typename TensorTools::DecrementRank
                       <typename TensorTools::DecrementRank
                        <OutputType>::type>::type>::type>,
                     &DiffContext::get_elem_fixed_solution>(var, qp, d2u);
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
                                  OutputType & u) const
{
  this->some_value
    <OutputType,
     &FEMContext::get_side_fe
     <typename TensorTools::MakeReal<OutputType>::type>,
     &DiffContext::get_elem_fixed_solution>
    (var, qp, u);
}



Gradient FEMContext::fixed_side_gradient(unsigned int var, unsigned int qp) const
{
  Gradient du;

  this->fixed_side_gradient( var, qp, du );

  return du;
}


template<typename OutputType>
void FEMContext::fixed_side_gradient(unsigned int var, unsigned int qp,
                                     OutputType & du) const
{
  this->some_gradient<OutputType,
                      &FEMContext::get_side_fe
                      <typename TensorTools::MakeReal
                       <typename TensorTools::DecrementRank
                        <OutputType>::type>::type>,
                      &DiffContext::get_elem_fixed_solution>(var, qp, du);
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
                                    OutputType & d2u) const
{
  this->some_hessian<OutputType,
                     &FEMContext::get_side_fe
                     <typename TensorTools::MakeReal
                      <typename TensorTools::DecrementRank
                       <typename TensorTools::DecrementRank
                        <OutputType>::type>::type>::type>,
                     &DiffContext::get_elem_fixed_solution>(var, qp, d2u);
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::fixed_point_value(unsigned int var, const Point & p) const
{
  Number u = 0.;

  this->fixed_point_value( var, p, u );

  return u;
}

template<typename OutputType>
void FEMContext::fixed_point_value(unsigned int var,
                                   const Point & p,
                                   OutputType & u,
                                   const Real tolerance) const
{
  typedef typename TensorTools::MakeReal<OutputType>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (_elem_fixed_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_fixed_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Build a FE for calculating u(p)
  UniquePtr<FEGenericBase<OutputShape> > fe_new =
    this->build_new_fe( fe, p, tolerance );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<OutputShape> > &  phi = fe_new->get_phi();

  u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][0] * coef(l);

  return;
}



Gradient FEMContext::fixed_point_gradient(unsigned int var, const Point & p) const
{
  Gradient grad_u;

  this->fixed_point_gradient( var, p, grad_u );

  return grad_u;
}



template<typename OutputType>
void FEMContext::fixed_point_gradient(unsigned int var,
                                      const Point & p,
                                      OutputType & grad_u,
                                      const Real tolerance) const
{
  typedef typename TensorTools::MakeReal
    <typename TensorTools::DecrementRank<OutputType>::type>::type
    OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (_elem_fixed_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_fixed_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Build a FE for calculating u(p)
  UniquePtr<FEGenericBase<OutputShape> > fe_new =
    this->build_new_fe( fe, p, tolerance );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > &  dphi = fe_new->get_dphi();

  grad_u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    grad_u.add_scaled(dphi[l][0], coef(l));

  return;
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

Tensor FEMContext::fixed_point_hessian(unsigned int var, const Point & p) const
{
  Tensor hess_u;

  this->fixed_point_hessian( var, p, hess_u );

  return hess_u;
}



template<typename OutputType>
void FEMContext::fixed_point_hessian(unsigned int var,
                                     const Point & p,
                                     OutputType & hess_u,
                                     const Real tolerance) const
{
  typedef typename TensorTools::DecrementRank<OutputType>::type Rank1Decrement;
  typedef typename TensorTools::DecrementRank<Rank1Decrement>::type Rank2Decrement;
  typedef typename TensorTools::MakeReal<Rank2Decrement>::type OutputShape;

  // Get local-to-global dof index lookup
  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices(var).size());

  // Get current local coefficients
  libmesh_assert_greater (_elem_fixed_subsolutions.size(), var);
  const DenseSubVector<Number> & coef = this->get_elem_fixed_solution(var);

  // Get finite element object
  FEGenericBase<OutputShape> * fe = libmesh_nullptr;
  this->get_element_fe<OutputShape>( var, fe, this->get_elem_dim() );

  // Build a FE for calculating u(p)
  UniquePtr<FEGenericBase<OutputShape> > fe_new =
    this->build_new_fe( fe, p, tolerance );

  // Get the values of the shape function derivatives
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &  d2phi = fe_new->get_d2phi();

  hess_u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    hess_u.add_scaled(d2phi[l][0], coef(l));

  return;
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES



template<typename OutputType>
void FEMContext::interior_rate(unsigned int var, unsigned int qp,
                               OutputType & u) const
{
  this->some_value<OutputType,
                   &FEMContext::get_element_fe
                   <typename TensorTools::MakeReal<OutputType>::type>,
                   &DiffContext::get_elem_solution_rate>(var, qp, u);
}



template<typename OutputType>
void FEMContext::side_rate(unsigned int var, unsigned int qp,
                           OutputType & u) const
{
  this->some_value<OutputType,
                   &FEMContext::get_side_fe
                   <typename TensorTools::MakeReal<OutputType>::type>,
                   &DiffContext::get_elem_solution_rate>(var, qp, u);
}

template<typename OutputType>
void FEMContext::interior_accel(unsigned int var, unsigned int qp,
                                OutputType & u) const
{
  this->some_value<OutputType,
                   &FEMContext::get_element_fe
                   <typename TensorTools::MakeReal<OutputType>::type>,
                   &DiffContext::get_elem_solution_accel>(var, qp, u);
}



template<typename OutputType>
void FEMContext::side_accel(unsigned int var, unsigned int qp,
                            OutputType & u) const
{
  this->some_value<OutputType,
                   &FEMContext::get_side_fe
                   <typename TensorTools::MakeReal<OutputType>::type>,
                   &DiffContext::get_elem_solution_accel>(var, qp, u);
}



void FEMContext::elem_reinit(Real theta)
{
  // Update the "time" variable of this context object
  this->_update_time_from_system(theta);

  // Handle a moving element if necessary.
  if (_mesh_sys)
    {
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


void FEMContext::elem_fe_reinit(const std::vector<Point> * const pts)
{
  // Initialize all the interior FE objects on elem.
  // Logging of FE::reinit is done in the FE functions
  // We only reinit the FE objects for the current element
  // dimension
  const unsigned char dim = this->get_elem_dim();

  libmesh_assert( !_element_fe[dim].empty() );

  std::map<FEType, FEAbstract *>::iterator local_fe_end = _element_fe[dim].end();
  for (std::map<FEType, FEAbstract *>::iterator i = _element_fe[dim].begin();
       i != local_fe_end; ++i)
    {
      if(this->has_elem())
        i->second->reinit(&(this->get_elem()), pts);
      else
        // If !this->has_elem(), then we assume we are dealing with a SCALAR variable
        i->second->reinit(libmesh_nullptr);
    }
}


void FEMContext::side_fe_reinit ()
{
  // Initialize all the side FE objects on elem/side.
  // Logging of FE::reinit is done in the FE functions
  // We only reinit the FE objects for the current element
  // dimension
  const unsigned char dim = this->get_elem_dim();

  libmesh_assert( !_side_fe[dim].empty() );

  std::map<FEType, FEAbstract *>::iterator local_fe_end = _side_fe[dim].end();
  for (std::map<FEType, FEAbstract *>::iterator i = _side_fe[dim].begin();
       i != local_fe_end; ++i)
    {
      i->second->reinit(&(this->get_elem()), this->get_side());
    }
}



void FEMContext::edge_fe_reinit ()
{
  libmesh_assert_equal_to (this->get_elem_dim(), 3);

  // Initialize all the interior FE objects on elem/edge.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEAbstract *>::iterator local_fe_end = _edge_fe.end();
  for (std::map<FEType, FEAbstract *>::iterator i = _edge_fe.begin();
       i != local_fe_end; ++i)
    {
      i->second->edge_reinit(&(this->get_elem()), this->get_edge());
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
  unsigned int n_nodes = this->get_elem().n_nodes();

#ifndef NDEBUG
  const unsigned char dim = this->get_elem_dim();

  // For simplicity we demand that mesh coordinates be stored
  // in a format that allows a direct copy
  libmesh_assert(this->get_mesh_x_var() == libMesh::invalid_uint ||
                 (this->get_element_fe(this->get_mesh_x_var(), dim)->get_fe_type().family
                  == LAGRANGE &&
                  this->get_element_fe(this->get_mesh_x_var(), dim)->get_fe_type().order.get_order()
                  == this->get_elem().default_order()));
  libmesh_assert(this->get_mesh_y_var() == libMesh::invalid_uint ||
                 (this->get_element_fe(this->get_mesh_y_var(), dim)->get_fe_type().family
                  == LAGRANGE &&
                  this->get_element_fe(this->get_mesh_y_var(), dim)->get_fe_type().order.get_order()
                  == this->get_elem().default_order()));
  libmesh_assert(this->get_mesh_z_var() == libMesh::invalid_uint ||
                 (this->get_element_fe(this->get_mesh_z_var(), dim)->get_fe_type().family
                  == LAGRANGE &&
                  this->get_element_fe(this->get_mesh_z_var(), dim)->get_fe_type().order.get_order()
                  == this->get_elem().default_order()));
#endif

  // Get degree of freedom coefficients from point coordinates
  if (this->get_mesh_x_var() != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      (this->get_elem_solution(this->get_mesh_x_var()))(i) = this->get_elem().point(i)(0);

  if (this->get_mesh_y_var() != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      (this->get_elem_solution(this->get_mesh_y_var()))(i) = this->get_elem().point(i)(1);

  if (this->get_mesh_z_var() != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      (this->get_elem_solution(this->get_mesh_z_var()))(i) = this->get_elem().point(i)(2);
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
  unsigned int n_nodes = this->get_elem().n_nodes();

#ifndef NDEBUG
  const unsigned char dim = this->get_elem_dim();

  // For simplicity we demand that mesh coordinates be stored
  // in a format that allows a direct copy
  libmesh_assert(this->get_mesh_x_var() == libMesh::invalid_uint ||
                 (this->get_element_fe(this->get_mesh_x_var(), dim)->get_fe_type().family
                  == LAGRANGE &&
                  this->get_elem_solution(this->get_mesh_x_var()).size() == n_nodes));
  libmesh_assert(this->get_mesh_y_var() == libMesh::invalid_uint ||
                 (this->get_element_fe(this->get_mesh_y_var(), dim)->get_fe_type().family
                  == LAGRANGE &&
                  this->get_elem_solution(this->get_mesh_y_var()).size() == n_nodes));
  libmesh_assert(this->get_mesh_z_var() == libMesh::invalid_uint ||
                 (this->get_element_fe(this->get_mesh_z_var(), dim)->get_fe_type().family
                  == LAGRANGE &&
                  this->get_elem_solution(this->get_mesh_z_var()).size() == n_nodes));
#endif

  // Set the new point coordinates
  if (this->get_mesh_x_var() != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      const_cast<Elem &>(this->get_elem()).point(i)(0) =
        libmesh_real(this->get_elem_solution(this->get_mesh_x_var())(i));

  if (this->get_mesh_y_var() != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      const_cast<Elem &>(this->get_elem()).point(i)(1) =
        libmesh_real(this->get_elem_solution(this->get_mesh_y_var())(i));

  if (this->get_mesh_z_var() != libMesh::invalid_uint)
    for (unsigned int i=0; i != n_nodes; ++i)
      const_cast<Elem &>(this->get_elem()).point(i)(2) =
        libmesh_real(this->get_elem_solution(this->get_mesh_z_var())(i));
  //    }
  // FIXME - If the coordinate data is not in our own system, someone
  // had better get around to implementing that... - RHS
  //  else
  //    {
  //      libmesh_not_implemented();
  //    }
}





/*
  void FEMContext::reinit(const FEMSystem & sys, Elem * e)
  {
  // Initialize our elem pointer, algebraic objects
  this->pre_fe_reinit(e);

  // Moving the mesh may be necessary
  // Reinitializing the FE objects is definitely necessary
  this->elem_reinit(1.);
  }
*/



void FEMContext::pre_fe_reinit(const System & sys, const Elem * e)
{
  this->set_elem(e);

  if (algebraic_type() == CURRENT ||
      algebraic_type() == DOFS_ONLY)
    {
      // Initialize the per-element data for elem.
      if(this->has_elem())
        sys.get_dof_map().dof_indices (&(this->get_elem()), this->get_dof_indices());
      else
        // If !this->has_elem(), then we assume we are dealing with a SCALAR variable
        sys.get_dof_map().dof_indices
          (static_cast<Elem*>(libmesh_nullptr), this->get_dof_indices());
    }
#ifdef LIBMESH_ENABLE_AMR
  else if (algebraic_type() == OLD)
    {
      // Initialize the per-element data for elem.
      if(this->has_elem())
        sys.get_dof_map().old_dof_indices (&(this->get_elem()), this->get_dof_indices());
      else
        // If !this->has_elem(), then we assume we are dealing with a SCALAR variable
        sys.get_dof_map().old_dof_indices
          (static_cast<Elem*>(libmesh_nullptr), this->get_dof_indices());
    }
#endif // LIBMESH_ENABLE_AMR

  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices().size());
  const std::size_t n_qoi = sys.qoi.size();

  if (this->algebraic_type() != NONE &&
      this->algebraic_type() != DOFS_ONLY)
    {
      // This also resizes elem_solution
      if (_custom_solution == libmesh_nullptr)
        sys.current_local_solution->get(this->get_dof_indices(), this->get_elem_solution().get_values());
      else
        _custom_solution->get(this->get_dof_indices(), this->get_elem_solution().get_values());

      if (sys.use_fixed_solution)
        this->get_elem_fixed_solution().resize(n_dofs);

      // Only make space for these if we're using DiffSystem
      // This is assuming *only* DiffSystem is using elem_solution_rate/accel
      const DifferentiableSystem * diff_system = dynamic_cast<const DifferentiableSystem *>(&sys);
      if(diff_system)
        {
          // Now, we only need these if the solver is unsteady
          if( !diff_system->get_time_solver().is_steady() )
            {
              this->get_elem_solution_rate().resize(n_dofs);

              // We only need accel space if the TimeSolver is second order
              const UnsteadySolver & time_solver = cast_ref<const UnsteadySolver &>(diff_system->get_time_solver());

              if( time_solver.time_order() >= 2 || !diff_system->get_second_order_vars().empty() )
                this->get_elem_solution_accel().resize(n_dofs);
            }
        }

      if (algebraic_type() != OLD)
        {
          // These resize calls also zero out the residual and jacobian
          this->get_elem_residual().resize(n_dofs);
          this->get_elem_jacobian().resize(n_dofs, n_dofs);

          this->get_qoi_derivatives().resize(n_qoi);
          this->_elem_qoi_subderivatives.resize(n_qoi);
          for (std::size_t q=0; q != n_qoi; ++q)
            (this->get_qoi_derivatives())[q].resize(n_dofs);
        }
    }

  // Initialize the per-variable data for elem.
  {
    unsigned int sub_dofs = 0;
    for (unsigned int i=0; i != sys.n_vars(); ++i)
      {
        if (algebraic_type() == CURRENT ||
            algebraic_type() == DOFS_ONLY)
          {
            if(this->has_elem())
              sys.get_dof_map().dof_indices (&(this->get_elem()), this->get_dof_indices(i), i);
            else
              // If !this->has_elem(), then we assume we are dealing with a SCALAR variable
              sys.get_dof_map().dof_indices
                (static_cast<Elem*>(libmesh_nullptr), this->get_dof_indices(i), i);
          }
#ifdef LIBMESH_ENABLE_AMR
        else if (algebraic_type() == OLD)
          {
            if(this->has_elem())
              sys.get_dof_map().old_dof_indices (&(this->get_elem()), this->get_dof_indices(i), i);
            else
              // If !this->has_elem(), then we assume we are dealing with a SCALAR variable
              sys.get_dof_map().old_dof_indices
                (static_cast<Elem*>(libmesh_nullptr), this->get_dof_indices(i), i);
          }
#endif // LIBMESH_ENABLE_AMR

        if (this->algebraic_type() != NONE &&
            this->algebraic_type() != DOFS_ONLY)
          {
            const unsigned int n_dofs_var = cast_int<unsigned int>
              (this->get_dof_indices(i).size());

            this->get_elem_solution(i).reposition
              (sub_dofs, n_dofs_var);

            // Only make space for these if we're using DiffSystem
            // This is assuming *only* DiffSystem is using elem_solution_rate/accel
            const DifferentiableSystem * diff_system = dynamic_cast<const DifferentiableSystem *>(&sys);
            if(diff_system)
              {
                // Now, we only need these if the solver is unsteady
                if( !diff_system->get_time_solver().is_steady() )
                  {
                    this->get_elem_solution_rate(i).reposition
                      (sub_dofs, n_dofs_var);

                    // We only need accel space if the TimeSolver is second order
                    const UnsteadySolver & time_solver = cast_ref<const UnsteadySolver &>(diff_system->get_time_solver());

                    if( time_solver.time_order() >= 2 || !diff_system->get_second_order_vars().empty() )
                      this->get_elem_solution_accel(i).reposition
                        (sub_dofs, n_dofs_var);
                  }
              }

            if (sys.use_fixed_solution)
              this->get_elem_fixed_solution(i).reposition
                (sub_dofs, n_dofs_var);

            if (algebraic_type() != OLD)
              {
                this->get_elem_residual(i).reposition
                  (sub_dofs, n_dofs_var);

                for (std::size_t q=0; q != n_qoi; ++q)
                  this->get_qoi_derivatives(q,i).reposition
                    (sub_dofs, n_dofs_var);

                for (unsigned int j=0; j != i; ++j)
                  {
                    const unsigned int n_dofs_var_j =
                      cast_int<unsigned int>
                      (this->get_dof_indices(j).size());

                    this->get_elem_jacobian(i,j).reposition
                      (sub_dofs, this->get_elem_residual(j).i_off(),
                       n_dofs_var, n_dofs_var_j);
                    this->get_elem_jacobian(j,i).reposition
                      (this->get_elem_residual(j).i_off(), sub_dofs,
                       n_dofs_var_j, n_dofs_var);
                  }
                this->get_elem_jacobian(i,i).reposition
                  (sub_dofs, sub_dofs,
                   n_dofs_var,
                   n_dofs_var);
              }

            sub_dofs += n_dofs_var;
          }
      }

    if (this->algebraic_type() != NONE &&
        this->algebraic_type() != DOFS_ONLY &&
        this->algebraic_type() != OLD)
      libmesh_assert_equal_to (sub_dofs, n_dofs);
  }

  // Now do the localization for the user requested vectors
  if (this->algebraic_type() != NONE)
    {
      DiffContext::localized_vectors_iterator localized_vec_it = this->_localized_vectors.begin();
      const DiffContext::localized_vectors_iterator localized_vec_end = this->_localized_vectors.end();

      for(; localized_vec_it != localized_vec_end; ++localized_vec_it)
        {
          const NumericVector<Number> & current_localized_vector = *localized_vec_it->first;
          DenseVector<Number> & target_vector = localized_vec_it->second.first;

          current_localized_vector.get(this->get_dof_indices(), target_vector.get_values());

          // Initialize the per-variable data for elem.
          unsigned int sub_dofs = 0;
          for (unsigned int i=0; i != sys.n_vars(); ++i)
            {
              const unsigned int n_dofs_var = cast_int<unsigned int>
                (this->get_dof_indices(i).size());

              // This is redundant with earlier initialization, isn't it? - RHS
              // sys.get_dof_map().dof_indices (&(this->get_elem()), this->get_dof_indices(i), i);

              localized_vec_it->second.second[i]->reposition
                (sub_dofs, n_dofs_var);

              sub_dofs += n_dofs_var;
            }
          libmesh_assert_equal_to (sub_dofs, n_dofs);
        }
    }
}

void FEMContext::set_elem( const Elem * e )
{
  this->_elem = e;

  // If e is NULL, we assume it's SCALAR and set _elem_dim to 0.
  this->_elem_dim = this->_elem ? this->_elem->dim() : 0;
}

void FEMContext::_update_time_from_system(Real theta)
{
  // Update the "time" variable based on the value of theta.  For this
  // to work, we need to know the value of deltat, a pointer to which is now
  // stored by our parent DiffContext class.  Note: get_deltat_value() will
  // assert in debug mode if the requested pointer is NULL.
  const Real deltat = this->get_deltat_value();

  this->set_time(theta*(this->get_system_time() + deltat) + (1.-theta)*this->get_system_time());
}



template<typename OutputShape>
UniquePtr<FEGenericBase<OutputShape> >
FEMContext::build_new_fe( const FEGenericBase<OutputShape>* fe,
                          const Point & p,
                          const Real tolerance) const
{
  FEType fe_type = fe->get_fe_type();

  // If we don't have an Elem to evaluate on, then the only functions
  // we can sensibly evaluate are the scalar dofs which are the same
  // everywhere.
  libmesh_assert(this->has_elem() || fe_type.family == SCALAR);

#ifdef LIBMESH_ENABLE_AMR
  if ((algebraic_type() == OLD) &&
      this->has_elem())
    {
      if (this->get_elem().p_refinement_flag() == Elem::JUST_REFINED)
        fe_type.order = static_cast<Order>(fe_type.order - 1);
      else if (this->get_elem().p_refinement_flag() == Elem::JUST_COARSENED)
        fe_type.order = static_cast<Order>(fe_type.order + 1);
    }
#endif // LIBMESH_ENABLE_AMR

  const unsigned int elem_dim = this->has_elem() ? this->get_elem().dim() : 0;

  FEGenericBase<OutputShape>* fe_new =
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    (this->has_elem() && this->get_elem().infinite()) ?
    FEGenericBase<OutputShape>::build_InfFE(elem_dim, fe_type).release() :
#endif
    FEGenericBase<OutputShape>::build(elem_dim, fe_type).release();

  // Map the physical co-ordinates to the master co-ordinates using the inverse_map from fe_interface.h
  // Build a vector of point co-ordinates to send to reinit
  Point master_point = this->has_elem() ?
    FEInterface::inverse_map (elem_dim,
                              fe_type,
                              &this->get_elem(),
                              p,
                              tolerance) : Point(0);

  std::vector<Point> coor(1, master_point);

  // Reinitialize the element and compute the shape function values at coor
  if(this->has_elem())
    fe_new->reinit (&this->get_elem(), &coor);
  else
    // If !this->has_elem(), then we assume we are dealing with a SCALAR variable
    fe_new->reinit (libmesh_nullptr, &coor);

  return UniquePtr<FEGenericBase<OutputShape> >(fe_new);
}





// Instantiate member function templates
template void FEMContext::interior_value<Number>(unsigned int, unsigned int, Number &) const;
template void FEMContext::interior_values<Number>(unsigned int, const NumericVector<Number> &,
                                                  std::vector<Number> &) const;
template void FEMContext::interior_value<Gradient>(unsigned int, unsigned int, Gradient &) const;
template void FEMContext::interior_values<Gradient>(unsigned int, const NumericVector<Number> &,
                                                    std::vector<Gradient> &) const;

template void FEMContext::interior_gradient<Gradient>(unsigned int, unsigned int, Gradient &) const;
template void FEMContext::interior_gradients<Gradient>(unsigned int, const NumericVector<Number> &,
                                                       std::vector<Gradient> &) const;
template void FEMContext::interior_gradient<Tensor>(unsigned int, unsigned int, Tensor &) const;
template void FEMContext::interior_gradients<Tensor>(unsigned int, const NumericVector<Number> &,
                                                     std::vector<Tensor> &) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::interior_hessian<Tensor>(unsigned int, unsigned int, Tensor &) const;
template void FEMContext::interior_hessians<Tensor>(unsigned int, const NumericVector<Number> &,
                                                    std::vector<Tensor> &) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::interior_hessian<??>(unsigned int, unsigned int, ??&) const;
//template void FEMContext::interior_hessians<??>(unsigned int, const NumericVector<Number> &,
//                                                std::vector<??> &) const;
#endif

template void FEMContext::interior_curl<Gradient>(unsigned int, unsigned int, Gradient &) const;

template void FEMContext::interior_div<Number>(unsigned int, unsigned int, Number &) const;

template void FEMContext::side_value<Number>(unsigned int, unsigned int, Number &) const;
template void FEMContext::side_value<Gradient>(unsigned int, unsigned int, Gradient &) const;
template void FEMContext::side_values<Number>(unsigned int, const NumericVector<Number> &,
                                              std::vector<Number> &) const;
template void FEMContext::side_values<Gradient>(unsigned int, const NumericVector<Number> &,
                                                std::vector<Gradient> &) const;

template void FEMContext::side_gradient<Gradient>(unsigned int, unsigned int, Gradient &) const;
template void FEMContext::side_gradients<Gradient>(unsigned int, const NumericVector<Number> &,
                                                   std::vector<Gradient> &) const;
template void FEMContext::side_gradient<Tensor>(unsigned int, unsigned int, Tensor &) const;
template void FEMContext::side_gradients<Tensor>(unsigned int, const NumericVector<Number> &,
                                                 std::vector<Tensor> &) const;


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::side_hessian<Tensor>(unsigned int, unsigned int, Tensor &) const;
template void FEMContext::side_hessians<Tensor>(unsigned int, const NumericVector<Number> &,
                                                std::vector<Tensor> &) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::side_hessian<??>(unsigned int, unsigned int,
//                                           ??&) const;
//template void FEMContext::side_hessians<??>(unsigned int, const NumericVector<Number> &,
//                                            std::vector<??> &) const;
#endif

template void FEMContext::point_value<Number>(unsigned int, const Point &, Number &, const Real) const;
template void FEMContext::point_value<Gradient>(unsigned int, const Point &, Gradient &, const Real) const;

template void FEMContext::point_gradient<Gradient>(unsigned int, const Point &, Gradient &, const Real) const;
template void FEMContext::point_gradient<Tensor>(unsigned int, const Point &, Tensor &, const Real) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::point_hessian<Tensor>(unsigned int, const Point &, Tensor &, const Real) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::point_hessian<??>(unsigned int, const Point &, ??&) const;
#endif

template void FEMContext::point_curl<Gradient>(unsigned int, const Point &, Gradient &, const Real) const;

template void FEMContext::fixed_interior_value<Number>(unsigned int, unsigned int, Number &) const;
template void FEMContext::fixed_interior_value<Gradient>(unsigned int, unsigned int, Gradient &) const;

template void FEMContext::fixed_interior_gradient<Gradient>(unsigned int, unsigned int, Gradient &) const;
template void FEMContext::fixed_interior_gradient<Tensor>(unsigned int, unsigned int, Tensor &) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::fixed_interior_hessian<Tensor>(unsigned int, unsigned int, Tensor &) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::fixed_interior_hessian<??>(unsigned int, unsigned int, ??&) const;
#endif

template void FEMContext::fixed_side_value<Number>(unsigned int, unsigned int, Number &) const;
template void FEMContext::fixed_side_value<Gradient>(unsigned int, unsigned int, Gradient &) const;

template void FEMContext::fixed_side_gradient<Gradient>(unsigned int, unsigned int, Gradient &) const;
template void FEMContext::fixed_side_gradient<Tensor>(unsigned int, unsigned int, Tensor &) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::fixed_side_hessian<Tensor>(unsigned int, unsigned int, Tensor &) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::fixed_side_hessian<??>(unsigned int, unsigned int, ??&) const;
#endif

template void FEMContext::fixed_point_value<Number>(unsigned int, const Point &, Number &, const Real) const;
template void FEMContext::fixed_point_value<Gradient>(unsigned int, const Point &, Gradient &, const Real) const;

template void FEMContext::fixed_point_gradient<Gradient>(unsigned int, const Point &, Gradient &, const Real) const;
template void FEMContext::fixed_point_gradient<Tensor>(unsigned int, const Point &, Tensor &, const Real) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template void FEMContext::fixed_point_hessian<Tensor>(unsigned int, const Point &, Tensor &, const Real) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::fixed_point_hessian<??>(unsigned int, const Point &, ??&) const;
#endif

template void FEMContext::interior_rate<Number>(unsigned int, unsigned int, Number &) const;
template void FEMContext::interior_rate<Gradient>(unsigned int, unsigned int, Gradient &) const;

template void FEMContext::side_rate<Number>(unsigned int, unsigned int, Number &) const;
template void FEMContext::side_rate<Gradient>(unsigned int, unsigned int, Gradient &) const;

template void FEMContext::interior_accel<Number>(unsigned int, unsigned int, Number &) const;
template void FEMContext::interior_accel<Gradient>(unsigned int, unsigned int, Gradient &) const;

template void FEMContext::side_accel<Number>(unsigned int, unsigned int, Number &) const;
template void FEMContext::side_accel<Gradient>(unsigned int, unsigned int, Gradient &) const;

template UniquePtr<FEGenericBase<Real> >
FEMContext::build_new_fe(const FEGenericBase<Real>*,
                         const Point &,
                         const Real) const;

template UniquePtr<FEGenericBase<RealGradient> >
FEMContext::build_new_fe(const FEGenericBase<RealGradient>*,
                         const Point &,
                         const Real) const;

} // namespace libMesh
