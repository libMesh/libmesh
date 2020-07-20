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

#include "libmesh/rb_parametrized_function.h"
#include "libmesh/int_range.h"
#include "libmesh/point.h"

namespace libMesh
{

RBParametrizedFunction::RBParametrizedFunction()
:
requires_xyz_perturbations(false),
fd_delta(1.e-6)
{}

RBParametrizedFunction::~RBParametrizedFunction() {}

void RBParametrizedFunction::vectorized_evaluate(const RBParameters & mu,
                                                 const std::unordered_map<dof_id_type, std::vector<Point>> & all_xyz,
                                                 const std::unordered_map<dof_id_type, subdomain_id_type> & sbd_ids,
                                                 const std::unordered_map<dof_id_type, std::vector<std::vector<Point>> > & all_xyz_perturb,
                                                 std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & output)
{
  output.clear();

  for (auto it : all_xyz)
    {
      dof_id_type elem_id = it.first;
      const std::vector<Point> & xyz_vec = it.second;

      auto sbd_it = sbd_ids.find(elem_id);
      if (sbd_it == sbd_ids.end())
        libmesh_error_msg("Error: elem_id not found");
      subdomain_id_type subdomain_id = sbd_it->second;

      std::vector<std::vector<Number>> values(get_n_components());
      for (unsigned int comp=0; comp<get_n_components(); comp++)
        {
          values[comp].resize(xyz_vec.size());
          for (unsigned int qp : index_range(xyz_vec))
            {
               if (requires_xyz_perturbations)
                 {
                   auto xyz_perturb_it = all_xyz_perturb.find(elem_id);
                   if (xyz_perturb_it == all_xyz_perturb.end())
                     {
                       libmesh_error_msg("Error: elem_id not found");
                     }
                   const std::vector<std::vector<Point>> & qps_and_perturbs = xyz_perturb_it->second;

                   if (qp >= qps_and_perturbs.size())
                     libmesh_error_msg("Error: Invalid qp");

                   values[comp][qp] = evaluate(mu, comp, xyz_vec[qp], subdomain_id, qps_and_perturbs[qp]);
                 }
               else
                 {
                   std::vector<Point> empty_perturbs;
                   values[comp][qp] = evaluate(mu, comp, xyz_vec[qp], subdomain_id, empty_perturbs);
                 }
            }
        }
      output[elem_id] = values;
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
  const auto elem_it = preevaluated_values.find(elem_id);
  if (elem_it == preevaluated_values.end())
    libmesh_error_msg("Error: elem_id not found");

  const std::vector<std::vector<Number>> & values = elem_it->second;

  if (comp >= values.size())
    libmesh_error_msg("Error: invalid comp");

  if (qp >= values[comp].size())
    libmesh_error_msg("Error: invalid qp");

  return values[comp][qp];
}

}
