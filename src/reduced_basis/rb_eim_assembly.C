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

// rbOOmit includes
#include "libmesh/rb_eim_assembly.h"
#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_evaluation.h"

// libMesh includes
#include "libmesh/fem_context.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/int_range.h"

namespace libMesh
{

RBEIMAssembly::RBEIMAssembly(RBEIMConstruction & rb_eim_con_in,
                             unsigned int basis_function_index_in)
  : _rb_eim_con(rb_eim_con_in),
    _basis_function_index(basis_function_index_in),
    _ghosted_basis_function(NumericVector<Number>::build(rb_eim_con_in.get_explicit_system().comm()))
{
  // localize the vector that stores the basis function for this assembly object,
  // i.e. the vector that is used in evaluate_basis_function_at_quad_pts
#ifdef LIBMESH_ENABLE_GHOSTED
  _ghosted_basis_function->init (_rb_eim_con.get_explicit_system().n_dofs(),
                                 _rb_eim_con.get_explicit_system().n_local_dofs(),
                                 _rb_eim_con.get_explicit_system().get_dof_map().get_send_list(),
                                 false,
                                 GHOSTED);
  _rb_eim_con.get_rb_evaluation().get_basis_function(_basis_function_index).
    localize(*_ghosted_basis_function,
             _rb_eim_con.get_explicit_system().get_dof_map().get_send_list());
#else
  _ghosted_basis_function->init (_rb_eim_con.get_explicit_system().n_dofs(), false, SERIAL);
  _rb_eim_con.get_rb_evaluation().get_basis_function(_basis_function_index).
    localize(*_ghosted_basis_function);
#endif

  initialize_fe();
}

RBEIMAssembly::~RBEIMAssembly()
{
}

void RBEIMAssembly::evaluate_basis_function(unsigned int var,
                                            const Elem & element,
                                            const std::vector<Point> & points,
                                            std::vector<Number> & values)
{
  LOG_SCOPE("evaluate_basis_function", "RBEIMAssembly");

  const std::vector<std::vector<Real>> & phi = get_fe().get_phi();

  // The FE object caches data, hence we recompute as little as
  // possible on the call to reinit.
  get_fe().reinit (&element, &points);

  std::vector<dof_id_type> dof_indices_var;

  DofMap & dof_map = get_rb_eim_construction().get_explicit_system().get_dof_map();
  dof_map.dof_indices (&element, dof_indices_var, var);

  libmesh_assert(dof_indices_var.size() == phi.size());

  unsigned int n_points = points.size();
  values.resize(n_points);

  for (unsigned int pt_index=0; pt_index<n_points; pt_index++)
    {
      values[pt_index] = 0.;
      for (auto i : index_range(dof_indices_var))
        values[pt_index] += (*_ghosted_basis_function)(dof_indices_var[i]) * phi[i][pt_index];
    }
}

RBEIMConstruction & RBEIMAssembly::get_rb_eim_construction()
{
  return _rb_eim_con;
}

NumericVector<Number> & RBEIMAssembly::get_ghosted_basis_function()
{
  return *_ghosted_basis_function;
}

FEBase & RBEIMAssembly::get_fe()
{
  return *_fe;
}

void RBEIMAssembly::initialize_fe()
{
  DofMap & dof_map = get_rb_eim_construction().get_explicit_system().get_dof_map();

  const unsigned int dim =
    get_rb_eim_construction().get_mesh().mesh_dimension();

  FEType fe_type = dof_map.variable_type(0);
  _fe = FEBase::build(dim, fe_type);

  // Pre-request the shape function for efficieny's sake
  _fe->get_phi();
}

}
