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
#include "libmesh/rb_parameters.h"

namespace libMesh
{

RBParametrizedFunction::RBParametrizedFunction()
:
requires_xyz_perturbations(false),
is_lookup_table(false),
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

void RBParametrizedFunction::vectorized_evaluate(const std::vector<RBParameters> & mus,
                                                 const std::vector<Point> & all_xyz,
                                                 const std::vector<subdomain_id_type> & sbd_ids,
                                                 const std::vector<std::vector<Point>> & all_xyz_perturb,
                                                 std::vector<std::vector<std::vector<Number>>> & output)
{
  LOG_SCOPE("vectorized_evaluate()", "RBParametrizedFunction");

  output.clear();
  unsigned int n_points = all_xyz.size();

  libmesh_error_msg_if(sbd_ids.size() != n_points, "Error: invalid vector sizes");
  libmesh_error_msg_if(requires_xyz_perturbations && (all_xyz_perturb.size() != n_points), "Error: invalid vector sizes");

  // Dummy vector to be used when xyz perturbations are not required
  std::vector<Point> empty_perturbs;

  output.resize(mus.size());
  for ( unsigned int mu_index : index_range(mus))
    {
      output[mu_index].resize(n_points);
      for (unsigned int point_index=0; point_index<n_points; point_index++)
        {
          if (requires_xyz_perturbations)
            {
              output[mu_index][point_index] = evaluate(mus[mu_index],
                                                       all_xyz[point_index],
                                                       sbd_ids[point_index],
                                                       all_xyz_perturb[point_index]);
            }
          else
            {
              output[mu_index][point_index] = evaluate(mus[mu_index],
                                                       all_xyz[point_index],
                                                       sbd_ids[point_index],
                                                       empty_perturbs);
            }
        }
    }
}

void RBParametrizedFunction::preevaluate_parametrized_function_on_mesh(const RBParameters & mu,
                                                                       const std::unordered_map<dof_id_type, std::vector<Point>> & all_xyz,
                                                                       const std::unordered_map<dof_id_type, subdomain_id_type> & sbd_ids,
                                                                       const std::unordered_map<dof_id_type, std::vector<std::vector<Point>> > & all_xyz_perturb)
{
  mesh_to_preevaluated_values_map.clear();

  unsigned int n_points = 0;
  for (const auto & xyz_pair : all_xyz)
  {
    const std::vector<Point> & xyz_vec = xyz_pair.second;
    n_points += xyz_vec.size();
  }

  std::vector<Point> all_xyz_vec(n_points);
  std::vector<subdomain_id_type> sbd_ids_vec(n_points);
  std::vector<std::vector<Point>> all_xyz_perturb_vec(n_points);

  // Empty vector to be used when xyz perturbations are not required
  std::vector<Point> empty_perturbs;

  unsigned int counter = 0;
  for (const auto & xyz_pair : all_xyz)
    {
      dof_id_type elem_id = xyz_pair.first;
      const std::vector<Point> & xyz_vec = xyz_pair.second;

      subdomain_id_type subdomain_id = libmesh_map_find(sbd_ids, elem_id);

      // The amount of data to be stored for each component
      auto n_qp = xyz_vec.size();
      mesh_to_preevaluated_values_map[elem_id].resize(n_qp);

      for (unsigned int qp : index_range(xyz_vec))
        {
          mesh_to_preevaluated_values_map[elem_id][qp] = counter;

          all_xyz_vec[counter] = xyz_vec[qp];
          sbd_ids_vec[counter] = subdomain_id;

          if (requires_xyz_perturbations)
            {
              const auto & qps_and_perturbs =
                libmesh_map_find(all_xyz_perturb, elem_id);
              libmesh_error_msg_if(qp >= qps_and_perturbs.size(), "Error: Invalid qp");

              all_xyz_perturb_vec[counter] = qps_and_perturbs[qp];
            }
          else
            {
              all_xyz_perturb_vec[counter] = empty_perturbs;
            }

          counter++;
        }
    }

  std::vector<RBParameters> mus {mu};
  vectorized_evaluate(mus,
                      all_xyz_vec,
                      sbd_ids_vec,
                      all_xyz_perturb_vec,
                      preevaluated_values);
}

Number RBParametrizedFunction::lookup_preevaluated_value_on_mesh(unsigned int comp,
                                                                 dof_id_type elem_id,
                                                                 unsigned int qp) const
{
  const std::vector<unsigned int> & indices_at_qps =
    libmesh_map_find(mesh_to_preevaluated_values_map, elem_id);

  libmesh_error_msg_if(qp >= indices_at_qps.size(), "Error: invalid qp");

  unsigned int index = indices_at_qps[qp];
  libmesh_error_msg_if(preevaluated_values.size() != 1, "Error: we expect only one parameter index");
  libmesh_error_msg_if(index >= preevaluated_values[0].size(), "Error: invalid index");

  return preevaluated_values[0][index][comp];
}

void RBParametrizedFunction::initialize_lookup_table()
{
  // No-op by default, override in subclasses as needed
}

}
