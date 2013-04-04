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

// libMesh includes
#include "libmesh/fem_context.h"
#include "libmesh/quadrature.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"

namespace libMesh
{

RBEIMAssembly::RBEIMAssembly(RBEIMConstruction& rb_eim_con_in,
                             unsigned int basis_function_index_in)
: _rb_eim_con(rb_eim_con_in),
  _basis_function_index(basis_function_index_in),
  _ghosted_basis_function(NumericVector<Number>::build(libMesh::default_solver_package(), rb_eim_con_in.system().communicator()))
{
  // localize the vector that stores the basis function for this assembly object,
  // i.e. the vector that is used in evaluate_basis_function_at_quad_pts
#ifdef LIBMESH_ENABLE_GHOSTED
  _ghosted_basis_function->init (_rb_eim_con.n_dofs(), _rb_eim_con.n_local_dofs(),
                                 _rb_eim_con.get_dof_map().get_send_list(), false, GHOSTED);
  _rb_eim_con.get_rb_evaluation().get_basis_function(_basis_function_index).
      localize(*_ghosted_basis_function, _rb_eim_con.get_dof_map().get_send_list());
#else
  _ghosted_basis_function->init (_rb_eim_con.n_dofs(), false, SERIAL);
  _rb_eim_con.get_rb_evaluation().get_basis_function(_basis_function_index).
      localize(*_ghosted_basis_function);
#endif
}

void RBEIMAssembly::evaluate_basis_function(unsigned int var,
                                            const Elem& element,
                                            const std::vector<Point>& qpoints,
                                            std::vector<Number>& values)
{
  DofMap& dof_map = _rb_eim_con.get_dof_map();

  // Get local coordinates to feed these into compute_data().
  // Note that the fe_type can safely be used from the 0-variable,
  // since the inverse mapping is the same for all FEFamilies
  std::vector<Point> mapped_qpoints;
  FEInterface::inverse_map (_rb_eim_con.get_mesh().mesh_dimension(),
   		            dof_map.variable_type(0),
                            &element,
                            qpoints,
                            mapped_qpoints);

  const FEType& fe_type = dof_map.variable_type(var);

  std::vector<dof_id_type> dof_indices_var;
  dof_map.dof_indices (&element, dof_indices_var, var);

  values.resize(mapped_qpoints.size());

  for(unsigned int qp=0; qp<mapped_qpoints.size(); qp++)
  {
    FEComputeData data (_rb_eim_con.get_equation_systems(), mapped_qpoints[qp]);
    FEInterface::compute_data (_rb_eim_con.get_mesh().mesh_dimension(), fe_type, &element, data);

    values[qp] = 0.;
    for (unsigned int i=0; i<dof_indices_var.size(); i++)
      values[qp] += (*_ghosted_basis_function)(dof_indices_var[i]) * data.shape[i];
  }

}

RBEIMConstruction& RBEIMAssembly::get_rb_eim_construction()
{
  return _rb_eim_con;
}

}
