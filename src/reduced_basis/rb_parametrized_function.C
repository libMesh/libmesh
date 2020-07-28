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

// libmesh includes
#include "libmesh/rb_parametrized_function.h"
#include "libmesh/int_range.h"
#include "libmesh/point.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/utility.h"

namespace libMesh
{

RBParametrizedFunction::RBParametrizedFunction()
:
requires_xyz_perturbations(false),
fd_delta(1.e-6)
{}

RBParametrizedFunction::~RBParametrizedFunction() = default;

Number
RBParametrizedFunction::evaluate_comp(const RBParameters & mu,
                                      unsigned int comp,
                                      const Point & xyz,
                                      subdomain_id_type subdomain_id,
                                      const std::vector<Point> & xyz_perturb)
{
  std::vector<Number> values = evaluate(mu, xyz, subdomain_id, xyz_perturb);

  libmesh_error_msg_if(comp >= values.size(), "Error: Invalid value of comp");

  return values[comp];
}

void RBParametrizedFunction::vectorized_evaluate(const RBParameters & mu,
                                                 const std::unordered_map<dof_id_type, std::vector<Point>> & all_xyz,
                                                 const std::unordered_map<dof_id_type, subdomain_id_type> & sbd_ids,
                                                 const std::unordered_map<dof_id_type, std::vector<std::vector<Point>> > & all_xyz_perturb,
                                                 std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & output)
{
  LOG_SCOPE("vectorized_evaluate()", "RBParametrizedFunction");

  output.clear();

  auto n_comp = get_n_components();

  // Dummy vector to be used when xyz perturbations are not required
  std::vector<Point> empty_perturbs;

  for (const auto & xyz_pair : all_xyz)
    {
      dof_id_type elem_id = xyz_pair.first;
      const std::vector<Point> & xyz_vec = xyz_pair.second;

      subdomain_id_type subdomain_id = libmesh_map_find(sbd_ids, elem_id);

      // The amount of data to be stored for each component
      auto n_qp = xyz_vec.size();

      // Get a reference to the array[n_comp][n_qp] of data for this
      // Elem, and properly resize it.
      auto & array = output[elem_id];
      array.resize(n_comp);
      for (auto comp : index_range(array))
        array[comp].resize(n_qp);

      for (unsigned int qp : index_range(xyz_vec))
        {
          std::vector<Number> evaluated_values_at_qp;
          if (requires_xyz_perturbations)
            {
              const auto & qps_and_perturbs =
                libmesh_map_find(all_xyz_perturb, elem_id);

              libmesh_error_msg_if(qp >= qps_and_perturbs.size(), "Error: Invalid qp");

              evaluated_values_at_qp = evaluate(mu, xyz_vec[qp], subdomain_id, qps_and_perturbs[qp]);
            }
          else
            {
              evaluated_values_at_qp = evaluate(mu, xyz_vec[qp], subdomain_id, empty_perturbs);
            }
          // Copy evaluated_values_at_qp into array at position "qp"
          for (auto comp : make_range(n_comp))
            array[comp][qp] = evaluated_values_at_qp[comp];
        }
    }
}

void RBParametrizedFunction::preevaluate_parametrized_function(const RBParameters & mu,
                                                               const std::unordered_map<dof_id_type, std::vector<Point>> & all_xyz,
                                                               const std::unordered_map<dof_id_type, subdomain_id_type> & sbd_ids,
                                                               const std::unordered_map<dof_id_type, std::vector<std::vector<Point>> > & all_xyz_perturb)
{
  vectorized_evaluate(mu, all_xyz, sbd_ids, all_xyz_perturb, preevaluated_values);
}

Number RBParametrizedFunction::lookup_preevaluated_value(unsigned int comp,
                                                         dof_id_type elem_id,
                                                         unsigned int qp) const
{
  const auto & values =
    libmesh_map_find(preevaluated_values, elem_id);

  libmesh_error_msg_if(comp >= values.size(), "Error: invalid comp");
  libmesh_error_msg_if(qp >= values[comp].size(), "Error: invalid qp");

  return values[comp][qp];
}

}
