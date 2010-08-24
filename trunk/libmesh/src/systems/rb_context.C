// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "rb_context.h"
#include "rb_system.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "equation_systems.h"
#include "mesh_base.h"
#include "quadrature.h"
#include "dof_map.h"

namespace libMesh
{



RBContext::RBContext (const RBSystem& sys) :
  dof_indices_var(sys.n_vars()),
  element_qrule(NULL), side_qrule(NULL),
  elem(NULL),
  side(0),
  dim(sys.get_mesh().mesh_dimension())
{
  /* ---- Borrowed from DiffContext ctor ---- */
  unsigned int n_vars = sys.n_vars();
  libmesh_assert (n_vars);

  elem_subsolutions.reserve(n_vars);
  elem_subvectors.reserve(n_vars);
  elem_submatrices.resize(n_vars);

  for (unsigned int i=0; i != n_vars; ++i)
    {
      elem_subsolutions.push_back(new DenseSubVector<Number>(elem_solution));
      elem_subvectors.push_back(new DenseSubVector<Number>(elem_vector));
      elem_submatrices[i].reserve(n_vars);

      for (unsigned int j=0; j != n_vars; ++j)
        elem_submatrices[i].push_back
          (new DenseSubMatrix<Number>(elem_matrix));
    }

  /* ---- Borrowed from FEMContext ctor ---- */
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

  // Next, create finite element objects
  element_fe_var.resize(n_vars);
  side_fe_var.resize(n_vars);
  for (unsigned int i=0; i != n_vars; ++i)
    {
      FEType fe_type = sys.variable_type(i);
      if (element_fe[fe_type] == NULL)
        {
          element_fe[fe_type] = FEBase::build(dim, fe_type).release();
          element_fe[fe_type]->attach_quadrature_rule(element_qrule);
          side_fe[fe_type] = FEBase::build(dim, fe_type).release();
          side_fe[fe_type]->attach_quadrature_rule(side_qrule);
        }
      element_fe_var[i] = element_fe[fe_type];
      side_fe_var[i] = side_fe[fe_type];
    }
}

RBContext::~RBContext ()
{
  /* ---- Borrowed from FEMContext dtor ---- */
  for (unsigned int i=0; i != elem_subsolutions.size(); ++i)
    {
      delete elem_subsolutions[i];
      delete elem_subvectors[i];

      for (unsigned int j=0; j != elem_submatrices[i].size(); ++j)
        delete elem_submatrices[i][j];
    }

  /* ---- Borrowed from FEMContext dtor ---- */
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

  delete element_qrule;
  element_qrule = NULL;

  delete side_qrule;
  side_qrule = NULL;
}

/* ---- All code below is borrowed from FEMContext ---- */

Number RBContext::interior_value(unsigned int var, unsigned int qp)
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

Gradient RBContext::interior_gradient(unsigned int var, unsigned int qp)
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

Number RBContext::side_value(unsigned int var, unsigned int qp)
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



Gradient RBContext::side_gradient(unsigned int var, unsigned int qp)
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

Number RBContext::point_value(unsigned int var, Point &p)
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

void RBContext::elem_fe_reinit ()
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


void RBContext::elem_side_fe_reinit ()
{
  // Initialize all the interior FE objects on elem/side.
  // Logging of FE::reinit is done in the FE functions
  std::map<FEType, FEBase *>::iterator fe_end = side_fe.end();
  for (std::map<FEType, FEBase *>::iterator i = side_fe.begin();
       i != fe_end; ++i)
    {
      i->second->reinit(elem, side);
    }
}

void RBContext::reinit(RBSystem &sys, Elem *e)
{
  elem = e;
  // Initialize the per-element data for elem.
  sys.get_dof_map().dof_indices (elem, dof_indices);
  unsigned int n_dofs = dof_indices.size();

  elem_solution.resize(n_dofs);

  for (unsigned int i=0; i != n_dofs; ++i)
    elem_solution(i) = sys.current_solution(dof_indices[i]);

  // These resize calls also zero out the vector and matrix
  elem_vector.resize(n_dofs);
  elem_matrix.resize(n_dofs, n_dofs);

  // Initialize the per-variable data for elem.
  unsigned int sub_dofs = 0;
  for (unsigned int i=0; i != sys.n_vars(); ++i)
    {
      sys.get_dof_map().dof_indices (elem, dof_indices_var[i], i);

      elem_subsolutions[i]->reposition
        (sub_dofs, dof_indices_var[i].size());

      elem_subvectors[i]->reposition
        (sub_dofs, dof_indices_var[i].size());

      for (unsigned int j=0; j != i; ++j)
        {
          elem_submatrices[i][j]->reposition
            (sub_dofs, elem_subvectors[j]->i_off(),
             dof_indices_var[i].size(),
             dof_indices_var[j].size());
          elem_submatrices[j][i]->reposition
            (elem_subvectors[j]->i_off(), sub_dofs,
             dof_indices_var[j].size(),
             dof_indices_var[i].size());
        }
      elem_submatrices[i][i]->reposition
        (sub_dofs, sub_dofs,
         dof_indices_var[i].size(),
         dof_indices_var[i].size());
      sub_dofs += dof_indices_var[i].size();
    }
  libmesh_assert(sub_dofs == n_dofs);

  // Reinitializing the FE objects is definitely necessary
  this->elem_fe_reinit();
}

} // namespace libMesh
