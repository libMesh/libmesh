// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
    edge_qrule(NULL), neighbor_qrule(NULL),
    _mesh_sys(sys.get_mesh_system()),
    _mesh_x_var(sys.get_mesh_x_var()),
    _mesh_y_var(sys.get_mesh_y_var()),
    _mesh_z_var(sys.get_mesh_z_var()),
    elem(NULL), neigh(NULL),  
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
  neighbor_qrule = hardest_fe_type.default_quadrature_rule
    (dim-1, sys.extra_quadrature_order).release();
  if (dim == 3)
    edge_qrule = hardest_fe_type.default_quadrature_rule
      (1, sys.extra_quadrature_order).release();

  // Next, create finite element objects
  element_fe_var.resize(n_vars);
  side_fe_var.resize(n_vars);
  neighbor_fe_var.resize(n_vars);
  if (dim == 3)
    edge_fe_var.resize(n_vars);

  for (unsigned int i=0; i != n_vars; ++i)
    {
      FEType fe_type = sys.variable_type(i);
      if (element_fe[fe_type] == NULL)
        {
          element_fe[fe_type] = FEBase::build(dim, fe_type).release();
          element_fe[fe_type]->attach_quadrature_rule(element_qrule);
          side_fe[fe_type] = FEBase::build(dim, fe_type).release();
          side_fe[fe_type]->attach_quadrature_rule(side_qrule);
          neighbor_fe[fe_type] = FEBase::build(dim, fe_type).release();
          neighbor_fe[fe_type]->attach_quadrature_rule(neighbor_qrule);

          if (dim == 3)
            {
              edge_fe[fe_type] = FEBase::build(dim, fe_type).release();
              edge_fe[fe_type]->attach_quadrature_rule(edge_qrule);
            }
        }
      element_fe_var[i] = element_fe[fe_type];
      side_fe_var[i] = side_fe[fe_type];
      neighbor_fe_var[i] = neighbor_fe[fe_type];
      if (dim == 3)
        edge_fe_var[i] = edge_fe[fe_type];
    }
}


FEMContext::~FEMContext()
{
  // We don't want to store AutoPtrs in STL containers, but we don't
  // want to leak memory either
  for (std::map<FEType, FEBase *>::iterator i = element_fe.begin();
       i != element_fe.end(); ++i)
    delete i->second;
  element_fe.clear();

  for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
       i != side_fe.end(); ++i)
    delete i->second;
  side_fe.clear();

  for (std::map<FEType, FEBase *>::iterator i = neighbor_fe.begin();
       i != neighbor_fe.end(); ++i)
    delete i->second;
  neighbor_fe.clear();

  for (std::map<FEType, FEBase *>::iterator i = edge_fe.begin();
       i != edge_fe.end(); ++i)
    delete i->second;
  edge_fe.clear();

  delete element_qrule;
  element_qrule = NULL;

  delete side_qrule;
  side_qrule = NULL;

  delete neighbor_qrule;
  neighbor_qrule = NULL;

  if (edge_qrule)
    delete edge_qrule;
  side_qrule = NULL;
}



Number FEMContext::interior_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    element_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return u;
}



Gradient FEMContext::interior_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    element_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return du;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::interior_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealTensor> > &d2phi =
    element_fe_var[var]->get_d2phi();

  // Accumulate solution second derivatives
  Tensor d2u;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return d2u;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::side_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    side_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return u;
}



Gradient FEMContext::side_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    side_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return du;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::side_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealTensor> > &d2phi =
    side_fe_var[var]->get_d2phi();

  // Accumulate solution second derivatives
  Tensor d2u;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return d2u;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::neighbor_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (neigh_dof_indices.size() > var);
  const unsigned int n_dofs = neigh_dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (neigh_subsolutions.size() > var);
  libmesh_assert (neigh_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *neigh_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    neighbor_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u = 0.;

  FEType fe_type = neighbor_fe_var[var]->get_fe_type();
  // Point p_master = FEInterface::inverse_map(dim, fe_type, elem, p);

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return u;
}



Gradient FEMContext::neighbor_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (neigh_dof_indices.size() > var);
  const unsigned int n_dofs = neigh_dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (neigh_subsolutions.size() > var);
  libmesh_assert (neigh_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *neigh_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    neighbor_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return du;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::neighbor_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (neigh_dof_indices.size() > var);
  const unsigned int n_dofs = neigh_dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (neigh_subsolutions.size() > var);
  libmesh_assert (neigh_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *neigh_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealTensor> > &d2phi =
    neighbor_fe_var[var]->get_d2phi();

  // Accumulate solution second derivatives
  Tensor d2u;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return d2u;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::point_value(unsigned int var, const Point &p)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  
  // Get current local coefficients
  libmesh_assert (elem_subsolutions.size() > var);
  libmesh_assert (elem_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_subsolutions[var];

  Number u = 0.;

  FEType fe_type = element_fe_var[var]->get_fe_type();
  Point p_master = FEInterface::inverse_map(dim, fe_type, elem, p);

  for (unsigned int l=0; l != n_dofs; l++)
    u += FEInterface::shape(dim, fe_type, elem, l, p_master)
         * coef(l);

  return u;
}



Number FEMContext::fixed_interior_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    element_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return u;
}



Gradient FEMContext::fixed_interior_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    element_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return du;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::fixed_interior_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealTensor> > &d2phi =
    element_fe_var[var]->get_d2phi();

  // Accumulate solution second derivatives
  Tensor d2u;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return d2u;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::fixed_side_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    side_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u += phi[l][qp] * coef(l);

  return u;
}



Gradient FEMContext::fixed_side_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    side_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du;

  for (unsigned int l=0; l != n_dofs; l++)
    du.add_scaled(dphi[l][qp], coef(l));

  return du;
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor FEMContext::fixed_side_hessian(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealTensor> > &d2phi =
    side_fe_var[var]->get_d2phi();

  // Accumulate solution second derivatives
  Tensor d2u;

  for (unsigned int l=0; l != n_dofs; l++)
    d2u.add_scaled(d2phi[l][qp], coef(l));

  return d2u;
}
#endif // ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES



Number FEMContext::fixed_point_value(unsigned int var, const Point &p)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  
  // Get current local coefficients
  libmesh_assert (elem_fixed_subsolutions.size() > var);
  libmesh_assert (elem_fixed_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_fixed_subsolutions[var];

  Number u = 0.;

  FEType fe_type = element_fe_var[var]->get_fe_type();
  Point p_master = FEInterface::inverse_map(dim, fe_type, elem, p);

  for (unsigned int l=0; l != n_dofs; l++)
    u += FEInterface::shape(dim, fe_type, elem, l, p_master)
         * coef(l);

  return u;
}



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
}


void FEMContext::side_fe_reinit ()
{
  // Initialize all the interior FE objects on elem/side.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEBase *>::iterator fe_end = side_fe.end();
  for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
       i != fe_end; ++i)
    {
      i->second->reinit(elem, side);
    }

  // Set pointer to neighbor Elem
  neigh = elem->neighbor(side);

  // Initialize all the interior FE objects on neighbor/side.
  // Logging of FE::reinit is done in the FE functions
  if (compute_neighbor_values)
    {
      // We will simulatneously iterate through side and neighbor FETypes.
      std::map<FEType, FEBase *>::iterator neigh_fe_end = neighbor_fe.end();
      std::map<FEType, FEBase *>::iterator i_side = side_fe.begin();
      // For loop takes care of iterating through neighbor
      for (std::map<FEType, FEBase *>::iterator i_neigh = neighbor_fe.begin();
           i_neigh != neigh_fe_end; ++i_neigh)
        {
          // The quadrature points on the side
          std::vector<Point > qface_point = i_side->second->get_xyz();
          // The points on the neighbor that we will reinit
          std::vector<Point> qface_neighbor_point;
          FEInterface::inverse_map (elem->dim(), i_side->first, neigh,
            qface_point, qface_neighbor_point);
          i_neigh->second->reinit(neigh, &qface_neighbor_point);
          // Iterate to next side FEType
          ++i_side;
        }
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
          elem->point(i)(0) =
            libmesh_real((*elem_subsolutions[_mesh_x_var])(i));

      if (_mesh_y_var != libMesh::invalid_uint)
        for (unsigned int i=0; i != n_nodes; ++i)
          elem->point(i)(1) =
            libmesh_real((*elem_subsolutions[_mesh_y_var])(i));

      if (_mesh_z_var != libMesh::invalid_uint)
        for (unsigned int i=0; i != n_nodes; ++i)
          elem->point(i)(2) =
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



void FEMContext::pre_fe_reinit(const System &sys, Elem *e)
{
  elem = e;

  // Initialize the per-element data for elem.
  sys.get_dof_map().dof_indices (elem, dof_indices);
  if (compute_neighbor_values)
    sys.get_dof_map().dof_indices (neigh, neigh_dof_indices);
  unsigned int n_dofs = dof_indices.size();
  unsigned int neigh_n_dofs;
  if (compute_neighbor_values)
    neigh_n_dofs = neigh_dof_indices.size();
  unsigned int n_qoi = sys.qoi.size();

  elem_solution.resize(n_dofs);
  if (compute_neighbor_values)
    neigh_solution.resize(n_dofs);
  if (sys.use_fixed_solution)
    elem_fixed_solution.resize(n_dofs);

  for (unsigned int i=0; i != n_dofs; ++i)
    elem_solution(i) = sys.current_solution(dof_indices[i]);
  if (compute_neighbor_values)
    for (unsigned int i=0; i != neigh_n_dofs; ++i)
      neigh_solution(i) = sys.current_solution(neigh_dof_indices[i]);

  // These resize calls also zero out the residual and jacobian
  elem_residual.resize(n_dofs);
  elem_jacobian.resize(n_dofs, n_dofs);
  if (compute_neighbor_values)
    {
      elem_residual.resize(n_dofs);
      elem_jacobian.resize(n_dofs, n_dofs);
    }

  elem_qoi_derivative.resize(n_qoi);
  elem_qoi_subderivatives.resize(n_qoi);
  for (unsigned int q=0; q != n_qoi; ++q)
    elem_qoi_derivative[q].resize(n_dofs);

  // Initialize the per-variable data for elem.
  unsigned int sub_dofs = 0;
  unsigned int neigh_sub_dofs = 0;
  for (unsigned int i=0; i != sys.n_vars(); ++i)
    {
      sys.get_dof_map().dof_indices (elem, dof_indices_var[i], i);
      if (compute_neighbor_values)
        sys.get_dof_map().dof_indices (neigh, neigh_dof_indices_var[i], i);

      elem_subsolutions[i]->reposition
        (sub_dofs, dof_indices_var[i].size());
      if (compute_neighbor_values)
        neigh_subsolutions[i]->reposition
          (neigh_sub_dofs, neigh_dof_indices_var[i].size());

      if (sys.use_fixed_solution)
        elem_fixed_subsolutions[i]->reposition
          (sub_dofs, dof_indices_var[i].size());

      elem_subresiduals[i]->reposition
        (sub_dofs, dof_indices_var[i].size());
      if (compute_neighbor_values)
        neigh_subresiduals[i]->reposition
          (neigh_sub_dofs, neigh_dof_indices_var[i].size());

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
          if (compute_neighbor_values)
            {
              neigh_subjacobians[i][j]->reposition
                (neigh_sub_dofs, neigh_subresiduals[j]->i_off(),
                 neigh_dof_indices_var[i].size(),
                 neigh_dof_indices_var[j].size());
              neigh_subjacobians[j][i]->reposition
                (neigh_subresiduals[j]->i_off(), neigh_sub_dofs,
                 neigh_dof_indices_var[j].size(),
                 neigh_dof_indices_var[i].size());
            }
        }
      elem_subjacobians[i][i]->reposition
        (sub_dofs, sub_dofs,
         dof_indices_var[i].size(),
         dof_indices_var[i].size());
      if (compute_neighbor_values)
        neigh_subjacobians[i][i]->reposition
          (neigh_sub_dofs, neigh_sub_dofs,
           neigh_dof_indices_var[i].size(),
           neigh_dof_indices_var[i].size());
      sub_dofs += dof_indices_var[i].size();
      neigh_sub_dofs += dof_indices_var[i].size();
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


} // namespace libMesh
