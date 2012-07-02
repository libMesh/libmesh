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



#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "fem_context.h"
#include "libmesh_logging.h"
#include "mesh_base.h"
#include "quadrature.h"
#include "system.h"
#include "time_solver.h"
#include "unsteady_solver.h" // For euler_residual

namespace libMesh
{

FEMContext::FEMContext (const System &sys)
  : DiffContext(sys),
    element_qrule(NULL), side_qrule(NULL),
    edge_qrule(NULL),
    _mesh_sys(sys.get_mesh_system()),
    _mesh_x_var(sys.get_mesh_x_var()),
    _mesh_y_var(sys.get_mesh_y_var()),
    _mesh_z_var(sys.get_mesh_z_var()),
    elem(NULL),
    side(0), edge(0), dim(sys.get_mesh().mesh_dimension())
{
  // We need to know which of our variables has the hardest
  // shape functions to numerically integrate.

  unsigned int n_vars = sys.n_vars();

  libmesh_assert (n_vars);
  FEType hardest_fe_type = sys.variable_type(0);

  for (unsigned int i=0; i != n_vars; ++i)
    {
      FEType fe_type = sys.variable_type(i);

      // FIXME - we don't yet handle mixed finite elements from
      // different families which require different quadrature rules
      // libmesh_assert (fe_type.family == hardest_fe_type.family);

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
  element_fe_var.resize(n_vars);
  side_fe_var.resize(n_vars);
  if (dim == 3)
    edge_fe_var.resize(n_vars);

  _element_fe_var.resize(n_vars);
  _side_fe_var.resize(n_vars);
  if (dim == 3)
    _edge_fe_var.resize(n_vars);

  for (unsigned int i=0; i != n_vars; ++i)
    {
      FEType fe_type = sys.variable_type(i);

      // Preserving backward compatibility here for now
      // Should move to the protected/FEAbstract interface
      if( FEInterface::field_type( fe_type ) == TYPE_SCALAR )
	{
	  if( element_fe[fe_type] == NULL )
	    {
	      element_fe[fe_type] = FEBase::build(dim, fe_type).release();
	      element_fe[fe_type]->attach_quadrature_rule(element_qrule);
	      side_fe[fe_type] = FEBase::build(dim, fe_type).release();
	      side_fe[fe_type]->attach_quadrature_rule(side_qrule);
	      
	      if (dim == 3)
		{
		  edge_fe[fe_type] = FEBase::build(dim, fe_type).release();
		  edge_fe[fe_type]->attach_quadrature_rule(edge_qrule);
		}
	    }
	  element_fe_var[i] = element_fe[fe_type];
	  side_fe_var[i] = side_fe[fe_type];
	  if (dim == 3)
	    edge_fe_var[i] = edge_fe[fe_type];
	}

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
  // Preserving backward compatibility here for now
  // Should move to the protected/FEAbstract interface
  for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
       i != element_fe.end(); ++i)
    delete i->second;
  element_fe.clear();

  for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
       i != side_fe.end(); ++i)
    delete i->second;
  side_fe.clear();

  for (std::map<FEType, FEBase *>::iterator i = edge_fe.begin();
       i != edge_fe.end(); ++i)
    delete i->second;
  edge_fe.clear();


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

  if (edge_qrule)
    delete edge_qrule;
  side_qrule = NULL;
}



Number FEMContext::interior_value(unsigned int var, unsigned int qp) const
{
  Number u = 0.;

  this->interior_value<Real>( var, qp, u );

  return u;
}

template<typename OutputShape> 
void FEMContext::interior_value(unsigned int var, unsigned int qp,
				typename FEGenericBase<OutputShape>::OutputNumber& u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
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



Gradient FEMContext::interior_gradient(unsigned int var, unsigned int qp) const
{
  Gradient du;

  this->interior_gradient<Real>( var, qp, du );

  return du;
}

template<typename OutputShape>
void FEMContext::interior_gradient(unsigned int var, unsigned int qp,
                                   typename FEGenericBase<OutputShape>::OutputNumberGradient& du) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* fe = NULL;
  this->get_element_fe<OutputShape>( var, fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > &dphi = fe->get_dphi();

  // Accumulate solution derivatives
  du = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::interior_hessian(unsigned int var, unsigned int qp) const
{
  Tensor d2u;

  this->interior_hessian<Real>( var, qp, d2u );

  return d2u;
}

template<typename OutputShape>
void FEMContext::interior_hessian(unsigned int var, unsigned int qp,
				  typename FEGenericBase<OutputShape>::OutputNumberTensor& d2u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
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
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::side_value(unsigned int var, unsigned int qp) const
{
  Number u = 0.;

  this->side_value<Real>( var, qp, u );

  return u;
}


template<typename OutputShape> 
void FEMContext::side_value(unsigned int var, unsigned int qp,
			    typename FEGenericBase<OutputShape>::OutputNumber& u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* side_fe = NULL;
  this->get_side_fe<OutputShape>( var, side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > &phi = side_fe->get_phi();

  // Accumulate solution value
  u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return;
}



Gradient FEMContext::side_gradient(unsigned int var, unsigned int qp) const
{
  Gradient du;

  this->side_gradient<Real>( var, qp, du );

  return du;
}


template<typename OutputShape> 
void FEMContext::side_gradient(unsigned int var, unsigned int qp, 
			       typename FEGenericBase<OutputShape>::OutputNumberGradient& du) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* side_fe = NULL;
  this->get_side_fe<OutputShape>( var, side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector< typename FEGenericBase<OutputShape>::OutputGradient> > &dphi = side_fe->get_dphi();

  // Accumulate solution derivatives
  du = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::side_hessian(unsigned int var, unsigned int qp) const
{
  Tensor d2u;

  this->side_hessian<Real>( var, qp, d2u );

  return d2u;
}

template<typename OutputShape>
void FEMContext::side_hessian(unsigned int var, unsigned int qp, 
			      typename FEGenericBase<OutputShape>::OutputNumberTensor& d2u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* side_fe = NULL;
  this->get_side_fe<OutputShape>( var, side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &d2phi = side_fe->get_d2phi();

  // Accumulate solution second derivatives
  d2u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::point_value(unsigned int var, const Point &p) const
{
  Number u = 0.;

  this->point_value<Real>( var, p, u );

  return u;
}

template<typename OutputShape>
void FEMContext::point_value(unsigned int var, const Point &p,
			     typename FEGenericBase<OutputShape>::OutputNumber& u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
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

  this->point_gradient<Real>( var, p, grad_u );

  return grad_u;
}



template<typename OutputShape>
void FEMContext::point_gradient(unsigned int var, const Point &p,
				typename FEGenericBase<OutputShape>::OutputNumberGradient& grad_u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();


  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
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

  this->point_hessian<Real>( var, p, hess_u );

  return hess_u;
}


template<typename OutputShape>
void FEMContext::point_hessian(unsigned int var, const Point &p,
			       typename FEGenericBase<OutputShape>::OutputNumberTensor& hess_u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
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




Number FEMContext::fixed_interior_value(unsigned int var, unsigned int qp) const
{
  Number u = 0.;

  this->fixed_interior_value<Real>( var, qp, u );

  return u;
}



template<typename OutputShape>
void FEMContext::fixed_interior_value(unsigned int var, unsigned int qp,
			              typename FEGenericBase<OutputShape>::OutputNumber& u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
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

  this->fixed_interior_gradient<Real>( var, qp, du );

  return du;
}


template<typename OutputShape>
void FEMContext::FEMContext::fixed_interior_gradient(unsigned int var, unsigned int qp,
						     typename FEGenericBase<OutputShape>::OutputNumberGradient & du) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
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

  this->fixed_interior_hessian<Real>( var, qp, d2u );

  return d2u;
}


template<typename OutputShape>
void FEMContext::fixed_interior_hessian(unsigned int var, unsigned int qp,
					typename FEGenericBase<OutputShape>::OutputNumberTensor& d2u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
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

  this->fixed_side_value<Real>( var, qp, u );

  return u;
}


template<typename OutputShape> 
void FEMContext::fixed_side_value(unsigned int var, unsigned int qp,
			          typename FEGenericBase<OutputShape>::OutputNumber& u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* side_fe = NULL;
  this->get_side_fe<OutputShape>( var, side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<OutputShape> > &phi = side_fe->get_phi();

  // Accumulate solution value
  u = 0.0;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return;
}



Gradient FEMContext::fixed_side_gradient(unsigned int var, unsigned int qp) const
{
  Gradient du;

  this->fixed_side_gradient<Real>( var, qp, du );

  return du;
}


template<typename OutputShape> 
void FEMContext::fixed_side_gradient(unsigned int var, unsigned int qp, 
				     typename FEGenericBase<OutputShape>::OutputNumberGradient& du) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* side_fe = NULL;
  this->get_side_fe<OutputShape>( var, side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputGradient> > &dphi = side_fe->get_dphi();

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

  this->fixed_side_hessian<Real>( var, qp, d2u );

  return d2u;
}

template<typename OutputShape> 
void FEMContext::fixed_side_hessian(unsigned int var, unsigned int qp, 
				    typename FEGenericBase<OutputShape>::OutputNumberTensor& d2u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get finite element object
  FEGenericBase<OutputShape>* side_fe = NULL;
  this->get_side_fe<OutputShape>( var, side_fe );

  // Get shape function values at quadrature point
  const std::vector<std::vector<typename FEGenericBase<OutputShape>::OutputTensor> > &d2phi = side_fe->get_d2phi();

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

  this->fixed_point_value<Real>( var, p, u );

  return u;
}

template<typename OutputShape>
void FEMContext::fixed_point_value(unsigned int var, const Point &p,
			           typename FEGenericBase<OutputShape>::OutputNumber& u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
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

  this->fixed_point_gradient<Real>( var, p, grad_u );

  return grad_u;
}



template<typename OutputShape>
void FEMContext::fixed_point_gradient(unsigned int var, const Point &p,
				      typename FEGenericBase<OutputShape>::OutputNumberGradient& grad_u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
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

  this->fixed_point_hessian<Real>( var, p, hess_u );

  return hess_u;
}



template<typename OutputShape>
void FEMContext::fixed_point_hessian(unsigned int var, const Point &p,
				     typename FEGenericBase<OutputShape>::OutputNumberTensor& hess_u) const
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
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
        libmesh_assert(libMesh::n_threads() == 1);

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


void FEMContext::elem_fe_reinit ()
{
  // Initialize all the interior FE objects on elem.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEBase *>::iterator fe_end = element_fe.end();
  for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
       i != fe_end; ++i)
    {
      i->second->reinit(elem);
    }


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
  std::map<FEType, FEBase *>::iterator fe_end = side_fe.end();
  for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
       i != fe_end; ++i)
    {
      i->second->reinit(elem, side);
    }

  std::map<FEType, FEAbstract *>::iterator local_fe_end = _side_fe.end();
  for (std::map<FEType, FEAbstract *>::iterator i = _side_fe.begin();
       i != local_fe_end; ++i)
    {
      i->second->reinit(elem, side);
    }
}



void FEMContext::edge_fe_reinit ()
{
  libmesh_assert(dim == 3);

  // Initialize all the interior FE objects on elem/edge.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEBase *>::iterator fe_end = edge_fe.end();
  for (std::map<FEType, FEBase *>::iterator i = edge_fe.begin();
       i != fe_end; ++i)
    {
      i->second->edge_reinit(elem, edge);
    }

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
  libmesh_assert(libMesh::n_threads() == 1);

  // If the coordinate data is in our own system, it's already
  // been set up for us
//  if (_mesh_sys == this->number())
//    {
      unsigned int n_nodes = elem->n_nodes();
      // For simplicity we demand that mesh coordinates be stored
      // in a format that allows a direct copy
      libmesh_assert(_mesh_x_var == libMesh::invalid_uint ||
		     (element_fe_var[_mesh_x_var]->get_fe_type().family
                      == LAGRANGE &&
                      element_fe_var[_mesh_x_var]->get_fe_type().order
                      == elem->default_order()));
      libmesh_assert(_mesh_y_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_y_var]->get_fe_type().family
                      == LAGRANGE &&
                      element_fe_var[_mesh_y_var]->get_fe_type().order
                      == elem->default_order()));
      libmesh_assert(_mesh_z_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_z_var]->get_fe_type().family
                      == LAGRANGE &&
                      element_fe_var[_mesh_z_var]->get_fe_type().order
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
  libmesh_assert(libMesh::n_threads() == 1);

  // If the coordinate data is in our own system, it's already
  // been set up for us, and we can ignore our input parameter theta
//  if (_mesh_sys == this->number())
//    {
      unsigned int n_nodes = elem->n_nodes();
      // For simplicity we demand that mesh coordinates be stored
      // in a format that allows a direct copy
      libmesh_assert(_mesh_x_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_x_var]->get_fe_type().family
                      == LAGRANGE &&
                      elem_subsolutions[_mesh_x_var]->size() == n_nodes));
      libmesh_assert(_mesh_y_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_y_var]->get_fe_type().family
                      == LAGRANGE &&
                      elem_subsolutions[_mesh_y_var]->size() == n_nodes));
      libmesh_assert(_mesh_z_var == libMesh::invalid_uint ||
                     (element_fe_var[_mesh_z_var]->get_fe_type().family
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
  unsigned int n_dofs = dof_indices.size();
  unsigned int n_qoi = sys.qoi.size();

  elem_solution.resize(n_dofs);
  if (sys.use_fixed_solution)
    elem_fixed_solution.resize(n_dofs);

  for (unsigned int i=0; i != n_dofs; ++i)
    elem_solution(i) = sys.current_solution(dof_indices[i]);

  // These resize calls also zero out the residual and jacobian
  elem_residual.resize(n_dofs);
  elem_jacobian.resize(n_dofs, n_dofs);

  elem_qoi_derivative.resize(n_qoi);
  elem_qoi_subderivatives.resize(n_qoi);
  for (unsigned int q=0; q != n_qoi; ++q)
    elem_qoi_derivative[q].resize(n_dofs);

  // Initialize the per-variable data for elem.
  unsigned int sub_dofs = 0;
  for (unsigned int i=0; i != sys.n_vars(); ++i)
    {
      sys.get_dof_map().dof_indices (elem, dof_indices_var[i], i);

      elem_subsolutions[i]->reposition
        (sub_dofs, dof_indices_var[i].size());

      if (sys.use_fixed_solution)
        elem_fixed_subsolutions[i]->reposition
          (sub_dofs, dof_indices_var[i].size());

      elem_subresiduals[i]->reposition
        (sub_dofs, dof_indices_var[i].size());

      for (unsigned int q=0; q != n_qoi; ++q)
        elem_qoi_subderivatives[q][i]->reposition
          (sub_dofs, dof_indices_var[i].size());

      for (unsigned int j=0; j != i; ++j)
        {
          elem_subjacobians[i][j]->reposition
            (sub_dofs, elem_subresiduals[j]->i_off(),
             dof_indices_var[i].size(),
             dof_indices_var[j].size());
          elem_subjacobians[j][i]->reposition
            (elem_subresiduals[j]->i_off(), sub_dofs,
             dof_indices_var[j].size(),
             dof_indices_var[i].size());
        }
      elem_subjacobians[i][i]->reposition
        (sub_dofs, sub_dofs,
         dof_indices_var[i].size(),
         dof_indices_var[i].size());
      sub_dofs += dof_indices_var[i].size();
    }
  libmesh_assert(sub_dofs == n_dofs);
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
  AutoPtr<FEGenericBase<OutputShape> > fe_new(FEGenericBase<OutputShape>::build(elem->dim(), fe_type));

  // Map the physical co-ordinates to the master co-ordinates using the inverse_map from fe_interface.h
  // Build a vector of point co-ordinates to send to reinit
  std::vector<Point> coor(1, FEInterface::inverse_map(dim, fe_type, elem, p));

  // Reinitialize the element and compute the shape function values at coor
  fe_new->reinit (elem, &coor);

  return fe_new;
}





// Instantiate member function templates
template void FEMContext::interior_value<Real>(unsigned int, unsigned int, Number&) const;
template void FEMContext::interior_value<RealGradient>(unsigned int, unsigned int, Gradient&) const;

template void FEMContext::interior_gradient<Real>(unsigned int, unsigned int, FEGenericBase<Real>::OutputNumberGradient&) const;
template void FEMContext::interior_gradient<RealGradient>(unsigned int, unsigned int, FEGenericBase<RealGradient>::OutputNumberGradient&) const;

template void FEMContext::interior_hessian<Real>(unsigned int, unsigned int, FEGenericBase<Real>::OutputNumberTensor&) const;

//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::interior_hessian<RealGradient>(unsigned int, unsigned int, FEGenericBase<RealGradient>::OutputNumberTensor&) const;

template void FEMContext::side_value<Real>(unsigned int, unsigned int, Number&) const;
template void FEMContext::side_value<RealGradient>(unsigned int, unsigned int, Gradient&) const;

template void FEMContext::side_gradient<Real>(unsigned int, unsigned int, FEGenericBase<Real>::OutputNumberGradient&) const;
template void FEMContext::side_gradient<RealGradient>(unsigned int, unsigned int, FEGenericBase<RealGradient>::OutputNumberGradient&) const;

template void FEMContext::side_hessian<Real>(unsigned int, unsigned int, FEGenericBase<Real>::OutputNumberTensor&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::side_hessian<RealGradient>(unsigned int, unsigned int, FEGenericBase<RealGradient>::OutputNumberTensor&) const;

template void FEMContext::point_value<Real>(unsigned int, const Point&, Number&) const;
template void FEMContext::point_value<RealGradient>(unsigned int, const Point&, Gradient&) const;

template void FEMContext::point_gradient<Real>(unsigned int, const Point&, FEGenericBase<Real>::OutputNumberGradient&) const;
template void FEMContext::point_gradient<RealGradient>(unsigned int, const Point&, FEGenericBase<RealGradient>::OutputNumberGradient&) const;

template void FEMContext::point_hessian<Real>(unsigned int, const Point&, FEGenericBase<Real>::OutputNumberTensor&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::point_hessian<RealGradient>(unsigned int, const Point&, FEGenericBase<RealGradient>::OutputNumberTensor&) const;

template void FEMContext::fixed_interior_value<Real>(unsigned int, unsigned int, Number&) const;
template void FEMContext::fixed_interior_value<RealGradient>(unsigned int, unsigned int, Gradient&) const;

template void FEMContext::fixed_interior_gradient<Real>(unsigned int, unsigned int, FEGenericBase<Real>::OutputNumberGradient&) const;
template void FEMContext::fixed_interior_gradient<RealGradient>(unsigned int, unsigned int, FEGenericBase<RealGradient>::OutputNumberGradient&) const;

template void FEMContext::fixed_interior_hessian<Real>(unsigned int, unsigned int, FEGenericBase<Real>::OutputNumberTensor&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::fixed_interior_hessian<RealGradient>(unsigned int, unsigned int, FEGenericBase<RealGradient>::OutputNumberTensor&) const;

template void FEMContext::fixed_side_value<Real>(unsigned int, unsigned int, Number&) const;
template void FEMContext::fixed_side_value<RealGradient>(unsigned int, unsigned int, Gradient&) const;

template void FEMContext::fixed_side_gradient<Real>(unsigned int, unsigned int, FEGenericBase<Real>::OutputNumberGradient&) const;
template void FEMContext::fixed_side_gradient<RealGradient>(unsigned int, unsigned int, FEGenericBase<RealGradient>::OutputNumberGradient&) const;

template void FEMContext::fixed_side_hessian<Real>(unsigned int, unsigned int, FEGenericBase<Real>::OutputNumberTensor&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::fixed_side_hessian<RealGradient>(unsigned int, unsigned int, FEGenericBase<RealGradient>::OutputRealTensor&) const;

template void FEMContext::fixed_point_value<Real>(unsigned int, const Point&, Number&) const;
template void FEMContext::fixed_point_value<RealGradient>(unsigned int, const Point&, Gradient&) const;

template void FEMContext::fixed_point_gradient<Real>(unsigned int, const Point&, FEGenericBase<Real>::OutputNumberGradient&) const;
template void FEMContext::fixed_point_gradient<RealGradient>(unsigned int, const Point&, FEGenericBase<RealGradient>::OutputNumberGradient&) const;

template void FEMContext::fixed_point_hessian<Real>(unsigned int, const Point&, FEGenericBase<Real>::OutputNumberTensor&) const;
//FIXME: Not everything is implemented yet for second derivatives of RealGradients
//template void FEMContext::fixed_point_hessian<RealGradient>(unsigned int, const Point&, FEGenericBase<RealGradient>::OutputRealTensor&) const;

template AutoPtr<FEGenericBase<Real> > FEMContext::build_new_fe( const FEGenericBase<Real>*, const Point & ) const;
template AutoPtr<FEGenericBase<RealGradient> > FEMContext::build_new_fe( const FEGenericBase<RealGradient>*, const Point & ) const;


} // namespace libMesh
