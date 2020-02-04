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

#include "libmesh/dg_fem_context.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"

namespace libMesh
{

DGFEMContext::DGFEMContext (const System & sys)
  : FEMContext(sys),
    _neighbor(nullptr),
    _neighbor_dof_indices_var(sys.n_vars()),
    _dg_terms_active(false)
{
  unsigned int nv = sys.n_vars();
  libmesh_assert (nv);

  _neighbor_subresiduals.reserve(nv);
  _elem_elem_subjacobians.resize(nv);
  _elem_neighbor_subjacobians.resize(nv);
  _neighbor_elem_subjacobians.resize(nv);
  _neighbor_neighbor_subjacobians.resize(nv);

  for (unsigned int i=0; i != nv; ++i)
    {
      _neighbor_subresiduals.emplace_back(libmesh_make_unique<DenseSubVector<Number>>(_neighbor_residual));
      _elem_elem_subjacobians[i].reserve(nv);
      _elem_neighbor_subjacobians[i].reserve(nv);
      _neighbor_elem_subjacobians[i].reserve(nv);
      _neighbor_neighbor_subjacobians[i].reserve(nv);

      for (unsigned int j=0; j != nv; ++j)
        {
          _elem_elem_subjacobians[i].emplace_back(libmesh_make_unique<DenseSubMatrix<Number>>(_elem_elem_jacobian));
          _elem_neighbor_subjacobians[i].emplace_back(libmesh_make_unique<DenseSubMatrix<Number>>(_elem_neighbor_jacobian));
          _neighbor_elem_subjacobians[i].emplace_back(libmesh_make_unique<DenseSubMatrix<Number>>(_neighbor_elem_jacobian));
          _neighbor_neighbor_subjacobians[i].emplace_back(libmesh_make_unique<DenseSubMatrix<Number>>(_neighbor_neighbor_jacobian));
        }
    }

  _neighbor_side_fe_var.resize(nv);
  for (unsigned int i=0; i != nv; ++i)
    {
      FEType fe_type = sys.variable_type(i);

      if (_neighbor_side_fe[fe_type] == nullptr)
        _neighbor_side_fe[fe_type] = FEAbstract<>::build(this->_dim, fe_type);

      _neighbor_side_fe_var[i] = _neighbor_side_fe[fe_type].get();
    }
}

DGFEMContext::~DGFEMContext()
{
}

void DGFEMContext::side_fe_reinit()
{
  FEMContext::side_fe_reinit();

  // By default we assume that the DG terms are inactive
  // They are only active if neighbor_side_fe_reinit is called
  _dg_terms_active = false;
}

void DGFEMContext::neighbor_side_fe_reinit ()
{
  // Call this *after* side_fe_reinit

  // Initialize all the neighbor side FE objects based on inverse mapping
  // the quadrature points on the current side
  std::vector<Point> qface_side_points;
  std::vector<Point> qface_neighbor_points;
  for (auto & pr : _neighbor_side_fe)
    {
      FEType neighbor_side_fe_type = pr.first;
      FEAbstract<> * side_fe = _side_fe[this->get_dim()][neighbor_side_fe_type].get();
      qface_side_points = side_fe->get_xyz();

      FEMap::inverse_map (this->get_dim(), &get_neighbor(),
                          qface_side_points, qface_neighbor_points);

      pr.second->reinit(&get_neighbor(), &qface_neighbor_points);
    }

  // Set boolean flag to indicate that the DG terms are active on this element
  _dg_terms_active = true;

  // Also, initialize data required for DG assembly on the current side,
  // analogously to FEMContext::pre_fe_reinit

  // Initialize the per-element data for elem.
  get_system().get_dof_map().dof_indices (&get_neighbor(), _neighbor_dof_indices);

  const unsigned int n_dofs = cast_int<unsigned int>
    (this->get_dof_indices().size());
  const unsigned int n_neighbor_dofs = cast_int<unsigned int>
    (_neighbor_dof_indices.size());

  // These resize calls also zero out the residual and jacobian
  _neighbor_residual.resize(n_neighbor_dofs);
  _elem_elem_jacobian.resize(n_dofs, n_dofs);
  _elem_neighbor_jacobian.resize(n_dofs, n_neighbor_dofs);
  _neighbor_elem_jacobian.resize(n_neighbor_dofs, n_dofs);
  _neighbor_neighbor_jacobian.resize(n_neighbor_dofs, n_neighbor_dofs);

  // Initialize the per-variable data for elem.
  {
    unsigned int sub_dofs = 0;
    for (auto i : IntRange<unsigned int>(0, get_system().n_vars()))
      {
        get_system().get_dof_map().dof_indices (&get_neighbor(), _neighbor_dof_indices_var[i], i);

        const unsigned int n_dofs_var = cast_int<unsigned int>
          (_neighbor_dof_indices_var[i].size());

        _neighbor_subresiduals[i]->reposition
          (sub_dofs, n_dofs_var);

        for (unsigned int j=0; j != i; ++j)
          {
            const unsigned int n_dofs_var_j =
              cast_int<unsigned int>
              (this->get_dof_indices(j).size());

            _elem_elem_subjacobians[i][j]->reposition
              (sub_dofs, _neighbor_subresiduals[j]->i_off(),
               n_dofs_var, n_dofs_var_j);
            _elem_elem_subjacobians[j][i]->reposition
              (_neighbor_subresiduals[j]->i_off(), sub_dofs,
               n_dofs_var_j, n_dofs_var);

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
        _elem_elem_subjacobians[i][i]->reposition
          (sub_dofs, sub_dofs,
           n_dofs_var,
           n_dofs_var);
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

}

}
