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

// Configuration data
#include "libmesh_config.h"

// Depends on the QNTransientRBSystem, which requires SLEPC
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "qn_transient_rb_context.h"
#include "qn_transient_rb_system.h"
#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "equation_systems.h"
#include "mesh.h"
#include "quadrature.h"
#include "dof_map.h"

namespace libMesh
{

QNTransientRBContext::QNTransientRBContext (const QNTransientRBSystem& sys) :
  RBContext(sys)
{
  unsigned int n_vars = sys.n_vars();
  elem_old_subsolutions.reserve(n_vars);

  for (unsigned int i=0; i != n_vars; ++i)
  {
    elem_old_subsolutions.push_back(new DenseSubVector<Number>(elem_old_solution));
  }
}

QNTransientRBContext::~QNTransientRBContext ()
{
  for (unsigned int i=0; i != elem_old_subsolutions.size(); ++i)
  {
    delete elem_old_subsolutions[i];
  }
}

Number QNTransientRBContext::old_interior_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_old_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    element_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u_old = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u_old += phi[l][qp] * coef(l);

  return u_old;
}

Gradient QNTransientRBContext::old_interior_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_old_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    element_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient old_du;

  for (unsigned int l=0; l != n_dofs; l++)
    old_du.add_scaled(dphi[l][qp], coef(l));

  return old_du;
}

Number QNTransientRBContext::old_side_value(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_old_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<Real> > &phi =
    side_fe_var[var]->get_phi();

  // Accumulate solution value
  Number u_old = 0.;

  for (unsigned int l=0; l != n_dofs; l++)
    u_old += phi[l][qp] * coef(l);

  return u_old;
}

Gradient QNTransientRBContext::old_side_gradient(unsigned int var, unsigned int qp)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();

  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_old_subsolutions[var];

  // Get shape function values at quadrature point
  const std::vector<std::vector<RealGradient> > &dphi =
    side_fe_var[var]->get_dphi();

  // Accumulate solution derivatives
  Gradient du_old;

  for (unsigned int l=0; l != n_dofs; l++)
    du_old.add_scaled(dphi[l][qp], coef(l));

  return du_old;
}

Number QNTransientRBContext::old_point_value(unsigned int var, Point &p)
{
  // Get local-to-global dof index lookup
  libmesh_assert (dof_indices.size() > var);
  const unsigned int n_dofs = dof_indices_var[var].size();


  // Get current local coefficients
  libmesh_assert (elem_old_subsolutions.size() > var);
  libmesh_assert (elem_old_subsolutions[var] != NULL);
  DenseSubVector<Number> &coef = *elem_old_subsolutions[var];

  Number u_old = 0.;

  FEType fe_type = element_fe_var[var]->get_fe_type();
  Point p_master = FEInterface::inverse_map(dim, fe_type, elem, p);

  for (unsigned int l=0; l != n_dofs; l++)
    u_old += FEInterface::shape(dim, fe_type, elem, l, p_master)
         * coef(l);

  return u_old;
}

void QNTransientRBContext::reinit(RBSystem &sys, Elem *e)
{
  // sys must be a QNTransientRBSystem
  QNTransientRBSystem& qn_system = libmesh_cast_ref<QNTransientRBSystem&>(sys);

  // Copy most of the code from RBContext::reinit,
  // except replace call to current_solution
  // with current_newton_iterate, and add
  // old_solution data

  elem = e;
  // Initialize the per-element data for elem.
  qn_system.get_dof_map().dof_indices (elem, dof_indices);
  unsigned int n_dofs = dof_indices.size();

  elem_old_solution.resize(n_dofs);
  elem_solution.resize(n_dofs);

  for (unsigned int i=0; i != n_dofs; ++i)
  {
    elem_old_solution(i) = qn_system.old_solution(dof_indices[i]);
    elem_solution(i) = (*qn_system.current_newton_iterate)(dof_indices[i]);
  }

  // These resize calls also zero out the vector and matrix
  elem_vector.resize(n_dofs);
  elem_matrix.resize(n_dofs, n_dofs);

  // Initialize the per-variable data for elem.
  unsigned int sub_dofs = 0;
  for (unsigned int i=0; i != qn_system.n_vars(); ++i)
    {
      qn_system.get_dof_map().dof_indices (elem, dof_indices_var[i], i);

      elem_old_subsolutions[i]->reposition
        (sub_dofs, dof_indices_var[i].size());

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


#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
