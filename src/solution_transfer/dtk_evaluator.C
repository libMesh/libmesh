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

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_TRILINOS_HAVE_DTK

#include "libmesh/dtk_evaluator.h"

#include "libmesh/dof_map.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/system.h"

namespace libMesh
{

DTKEvaluator::DTKEvaluator(System & in_sys,
                           std::string var_name):
  sys(in_sys),
  current_local_solution(*sys.current_local_solution),
  es(in_sys.get_equation_systems()),
  mesh(sys.get_mesh()),
  dim(mesh.mesh_dimension()),
  dof_map(sys.get_dof_map()),
  var_num(sys.variable_number(var_name)),
  fe_type(dof_map.variable_type(var_num))
{}

DTKEvaluator::FieldContainerType
DTKEvaluator::evaluate(const Teuchos::ArrayRCP<int> & elements,
                       const Teuchos::ArrayRCP<double> & coords)
{
  unsigned int num_values = elements.size();

  Teuchos::ArrayRCP<Number> values(num_values);
  DataTransferKit::FieldContainer<Number> evaluations(values, 1);

  for(unsigned int i=0; i<num_values; i++)
    {
      Elem * elem = mesh.elem_ptr(elements[i]);

      Point p;

      for(unsigned int j=0; j<dim; j++)
        p(j) = coords[(j*num_values)+i];

      const Point mapped_point(FEInterface::inverse_map(dim, dof_map.variable_type(0), elem, p));

      FEComputeData data (es, mapped_point);
      FEInterface::compute_data (dim, fe_type, elem, data);

      std::vector<dof_id_type> dof_indices;
      dof_map.dof_indices(elem, dof_indices, var_num);

      Number value = 0;

      for (std::size_t j=0; j<dof_indices.size(); j++)
        value += current_local_solution(dof_indices[j]) * data.shape[j];

      values[i] = value;
    }

  return evaluations;
}

} // namespace libMesh

#endif
