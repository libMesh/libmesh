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



#include "libmesh/boundary_info.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/dg_fem_context.h"

namespace libMesh
{

DGFEMContext::DGEMContext (const System &sys)
  : FEMContext(sys),
    _neighbor(NULL),
    _neighbor_side_qrule(NULL),
    _dg_terms_active(false)
{
  // Create neighbor_side_qrule to match the type of side_qrule
  // We on't actually use this quadrature rule, since we reinit
  // _neighbor_side_fe objects based on inverse mapping the
  // quadrature points on elem's side
//  _neighbor_side_qrule = QBase::build(side_qrule.type(),
//                                      side_qrule.get_dim(),
//                                      side_qrule.get_order()).release();

  unsigned int nv = sys.n_vars();
  libmesh_assert (nv);
  
  _neighbor_side_fe_var.resize(nv);
  for (unsigned int i=0; i != nv; ++i)
    {
      FEType fe_type = sys.variable_type(i);

      if ( _neighbor_side_fe[fe_type] == NULL )
	{
	  _neighbor_side_fe[fe_type] = FEAbstract::build(dim, fe_type).release();
//	  _neighbor_side_fe[fe_type]->attach_quadrature_rule(_neighbor_side_qrule);
	}
      _neighbor_side_fe_var[i] = _neighbor_side_fe[fe_type];
    }
}

DGFEMContext::~DGFEMContext()
{
  for (std::map<FEType, FEBase *>::iterator i = _neighbor_side_fe.begin();
       i != _neighbor_side_fe.end(); ++i)
    delete i->second;
  _neighbor_side_fe.clear();
  
  if(_neighbor_side_qrule)
  {
    delete _neighbor_side_qrule;
  }
  _neighbor_side_qrule = NULL;
}

void DGEMContext::pre_fe_reinit(const System& sys, const Elem *e)
{
  FEMContext::pre_fe_reinit(sys, e);

  // By default we assume that the DG terms are inactive
  // They are only active if neighbor_side_fe_reinit is called
  _dg_terms_active = false;
}

void DGEMContext::neighbor_side_fe_reinit ()
{
  // Call this *after* side_fe_reinit

  // Initialize all the neighbor side FE objects based on inverse mapping
  // the quadrature points on the current side
  std::vector<Point> qface_neighbor_points;
  std::map<FEType, FEAbstract *>::iterator local_fe_end = _neighbor_side_fe.end();
  for (std::map<FEType, FEAbstract *>::iterator i = _neighbor_side_fe.begin();
       i != local_fe_end; ++i)
    {
      FEInterface::inverse_map (dim,
                                i->first,
                                get_neighbor(),
                                get_side_qrule().get_xyz(),
                                qface_neighbor_points);

      i->second->reinit(get_neighbor(), &qface_neighbor_points);
    }

  // Also, initialize data required for DG assembly. In FEMContext the analogue of this is
  // in pre_fe_reinit. In the DG context, it makes sense to resize data structures here,
  // since we may not always need to assemble DG terms. That is, this should enable us to use
  // DGFEMContext in a continuous Galerkin context without a significant performance hit.

  // Initialize the per-element data for elem.
  sys.get_dof_map().dof_indices (get_neighbor(), _neighbor_dof_indices);
  
  const unsigned int n_dofs = dof_indices.size();
  const unsigned int n_neighbor_dofs = libmesh_cast_int<unsigned int>
    (_neighbor_dof_indices.size());

  // These resize calls also zero out the residual and jacobian
  _neighbor_residual.resize(n_neighbor_dofs);
  _elem_neighbor_jacobian.resize(n_dofs, n_neighbor_dofs);
  _neighbor_elem_jacobian.resize(n_neighbor_dofs, n_dofs);
  _neighbor_neighbor_jacobian.resize(n_neighbor_dofs, n_neighbor_dofs);
  
  // Initialize the per-variable data for elem.
  {
    unsigned int sub_dofs = 0;
    for (unsigned int i=0; i != sys.n_vars(); ++i)
      {
        sys.get_dof_map().dof_indices (get_neighbor(), _neighbor_dof_indices_var[i], i);

        const unsigned int n_dofs_var = libmesh_cast_int<unsigned int>
          (_neighbor_dof_indices_var[i].size());

        _neighbor_subresiduals[i]->reposition
          (sub_dofs, n_dofs_var);

        for (unsigned int j=0; j != i; ++j)
          {
            const unsigned int n_dofs_var_j =
	      libmesh_cast_int<unsigned int>
                (dof_indices_var[j].size());

            _elem_neighbor_subjacobians[i][j]->reposition
              (sub_dofs, _neighbor_subresiduals[j]->i_off(),
               n_dofs_var, n_dofs_var_j);
            _elem_neighbor_subjacobians[j][i]->reposition
              (_neighbor_subresiduals[j]->i_off(), sub_dofs,
               n_dofs_var_j, n_dofs_var);

            _neighbor_elem_subjacobians[i][j]->reposition
              (sub_dofs, _neighbor_subresiduals[j]->i_off(),
               n_dofs_var, n_dofs_var_j);
            _neighbor_elem_subjacobians[j][i]->reposition
              (_neighbor_subresiduals[j]->i_off(), sub_dofs,
               n_dofs_var_j, n_dofs_var);

            _neighbor_neighbor_subjacobians[i][j]->reposition
              (sub_dofs, _neighbor_subresiduals[j]->i_off(),
               n_dofs_var, n_dofs_var_j);
            _neighbor_neighbor_subjacobians[j][i]->reposition
              (_neighbor_subresiduals[j]->i_off(), sub_dofs,
               n_dofs_var_j, n_dofs_var);
          }
        _elem_neighbor_subjacobians[i][i]->reposition
          (sub_dofs, sub_dofs,
           n_dofs_var,
           n_dofs_var);
        _neighbor_elem_subjacobians[i][i]->reposition
          (sub_dofs, sub_dofs,
           n_dofs_var,
           n_dofs_var);
        _neighbor_neighbor_subjacobians[i][i]->reposition
          (sub_dofs, sub_dofs,
           n_dofs_var,
           n_dofs_var);
        sub_dofs += n_dofs_var;
      }
    libmesh_assert_equal_to (sub_dofs, n_dofs);
  }
  
  // Set boolean flag to indicate that the DG terms are active on this element
  _dg_terms_active = true;
}

}

