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
#include "libmesh/system.h"
#include "libmesh/elem.h"
#include "libmesh/fem_context.h"
#include "libmesh/quadrature.h"

namespace libMesh
{

void VectorizedEvalInput::clear()
{
  all_xyz.clear();
  elem_ids.clear();
  qps.clear();
  sbd_ids.clear();
  all_xyz_perturb.clear();
  phi_i_qp.clear();
  side_indices.clear();
  boundary_ids.clear();
  node_ids.clear();
}

RBParametrizedFunction::RBParametrizedFunction()
:
requires_xyz_perturbations(false),
requires_all_elem_qp_data(false),
requires_all_elem_center_data(false),
is_lookup_table(false),
fd_delta(1.e-6),
_is_nodal_boundary(false)
{}

RBParametrizedFunction::~RBParametrizedFunction() = default;

Number
RBParametrizedFunction::evaluate_comp(const RBParameters & mu,
                                      unsigned int comp,
                                      const Point & xyz,
                                      dof_id_type elem_id,
                                      unsigned int qp,
                                      subdomain_id_type subdomain_id,
                                      const std::vector<Point> & xyz_perturb,
                                      const std::vector<Real> & phi_i_qp)
{
  std::vector<Number> values = evaluate(mu, xyz, elem_id, qp, subdomain_id, xyz_perturb, phi_i_qp);

  libmesh_error_msg_if(comp >= values.size(), "Error: Invalid value of comp");

  return values[comp];
}

Number
RBParametrizedFunction::side_evaluate_comp(const RBParameters & mu,
                                           unsigned int comp,
                                           const Point & xyz,
                                           dof_id_type elem_id,
                                           unsigned int side_index,
                                           unsigned int qp,
                                           subdomain_id_type subdomain_id,
                                           boundary_id_type boundary_id,
                                           const std::vector<Point> & xyz_perturb,
                                           const std::vector<Real> & phi_i_qp)
{
  std::vector<Number> values = side_evaluate(mu, xyz, elem_id, side_index, qp, subdomain_id, boundary_id, xyz_perturb, phi_i_qp);

  libmesh_error_msg_if(comp >= values.size(), "Error: Invalid value of comp");

  return values[comp];
}

Number
RBParametrizedFunction::node_evaluate_comp(const RBParameters & mu,
                                           unsigned int comp,
                                           const Point & xyz,
                                           dof_id_type node_id,
                                           boundary_id_type boundary_id)
{
  std::vector<Number> values = node_evaluate(mu, xyz, node_id, boundary_id);

  libmesh_error_msg_if(comp >= values.size(), "Error: Invalid value of comp");

  return values[comp];
}

std::vector<Number>
RBParametrizedFunction::evaluate(const RBParameters & /*mu*/,
                                 const Point & /*xyz*/,
                                 dof_id_type /*elem_id*/,
                                 unsigned int /*qp*/,
                                 subdomain_id_type /*subdomain_id*/,
                                 const std::vector<Point> & /*xyz_perturb*/,
                                 const std::vector<Real> & /*phi_i_qp*/)
{
  // This method should be overridden in subclasses, so we just give a not implemented error message here
  libmesh_not_implemented();

  return std::vector<Number>();
}

std::vector<Number>
RBParametrizedFunction::side_evaluate(const RBParameters & /*mu*/,
                                      const Point & /*xyz*/,
                                      dof_id_type /*elem_id*/,
                                      unsigned int /*side_index*/,
                                      unsigned int /*qp*/,
                                      subdomain_id_type /*subdomain_id*/,
                                      boundary_id_type /*boundary_id*/,
                                      const std::vector<Point> & /*xyz_perturb*/,
                                      const std::vector<Real> & /*phi_i_qp*/)
{
  // This method should be overridden in subclasses, so we just give a not implemented error message here
  libmesh_not_implemented();

  return std::vector<Number>();
}

std::vector<Number>
RBParametrizedFunction::node_evaluate(const RBParameters & /*mu*/,
                                      const Point & /*xyz*/,
                                      dof_id_type /*node_id*/,
                                      boundary_id_type /*boundary_id*/)
{
  // This method should be overridden in subclasses, so we just give a not implemented error message here
  libmesh_not_implemented();

  return std::vector<Number>();
}

void RBParametrizedFunction::vectorized_evaluate(const std::vector<RBParameters> & mus,
                                                 const VectorizedEvalInput & v,
                                                 std::vector<std::vector<std::vector<Number>>> & output)
{
  LOG_SCOPE("vectorized_evaluate()", "RBParametrizedFunction");

  output.clear();
  unsigned int n_points = v.all_xyz.size();

  libmesh_error_msg_if(v.sbd_ids.size() != n_points, "Error: invalid vector sizes");
  libmesh_error_msg_if(requires_xyz_perturbations && (v.all_xyz_perturb.size() != n_points), "Error: invalid vector sizes");

  // Dummy vector to be used when xyz perturbations are not required
  std::vector<Point> empty_perturbs;

  // The number of components returned by this RBParametrizedFunction
  auto n_components = this->get_n_components();

  // We first loop over all mus and all n_points, calling evaluate()
  // for each and storing the results.  It is easier to first
  // pre-compute all the values before filling output, since, in the
  // case of multi-sample RBParameters, the ordering of the loops is a
  // bit complicated otherwise.
  std::vector<std::vector<std::vector<Number>>> all_evals(mus.size());
  for (auto mu_index : index_range(mus))
    {
      // Allocate enough space to store all points for the current mu
      all_evals[mu_index].resize(n_points);
      for (auto point_index : index_range(all_evals[mu_index]))
        {
          // Evaluate all samples for the current mu at the current interpolation point
          all_evals[mu_index][point_index] =
            this->evaluate(mus[mu_index],
                           v.all_xyz[point_index],
                           v.elem_ids[point_index],
                           v.qps[point_index],
                           v.sbd_ids[point_index],
                           requires_xyz_perturbations ? v.all_xyz_perturb[point_index] : empty_perturbs,
                           v.phi_i_qp[point_index]);

          // The vector returned by evaluate() should contain:
          // n_components * mus[mu_index].n_samples()
          // entries. That is, for multi-sample RBParameters objects,
          // the vector will be packed with entries as follows:
          // [sample0_component0, sample0_component1, ..., sample0_componentN,
          //  sample1_component0, sample1_component1, ..., sample1_componentN,
          //  ...
          //  sampleM_component0, sampleM_component1, ..., sampleM_componentN]
          auto n_samples = mus[mu_index].n_samples();
          auto received_data = all_evals[mu_index][point_index].size();
          libmesh_error_msg_if(received_data != n_components * n_samples,
                               "Recieved " << received_data <<
                               " evaluated values but expected to receive " << n_components * n_samples);
        }
    }

  // TODO: move this code for computing the total number of samples
  // represented by a std::vector of RBParameters objects to a helper
  // function.
  unsigned int output_size = 0;
  for (const auto & mu : mus)
    output_size += mu.n_samples();

  output.resize(output_size);

  // We use traditional for-loops here (rather than range-based) so that we can declare and
  // increment multiple loop counters all within the local scope of the for-loop.
  for (auto [mu_index, output_index] = std::make_tuple(0u, 0u); mu_index < mus.size(); ++mu_index)
    {
      auto n_samples = mus[mu_index].n_samples();
      for (auto mu_sample_idx = 0u; mu_sample_idx < n_samples; ++mu_sample_idx, ++output_index)
        {
          output[output_index].resize(n_points);
          for (auto point_index : make_range(n_points))
            {
              output[output_index][point_index].resize(n_components);

              for (auto comp : make_range(n_components))
                output[output_index][point_index][comp] = all_evals[mu_index][point_index][n_components*mu_sample_idx + comp];
            }
        }
    }
}

void RBParametrizedFunction::side_vectorized_evaluate(const std::vector<RBParameters> & mus,
                                                      const VectorizedEvalInput & v,
                                                      std::vector<std::vector<std::vector<Number>>> & output)
{
  LOG_SCOPE("side_vectorized_evaluate()", "RBParametrizedFunction");

  output.clear();
  unsigned int n_points = v.all_xyz.size();

  libmesh_error_msg_if(v.sbd_ids.size() != n_points, "Error: invalid vector sizes");
  libmesh_error_msg_if(requires_xyz_perturbations && (v.all_xyz_perturb.size() != n_points), "Error: invalid vector sizes");

  // Dummy vector to be used when xyz perturbations are not required
  std::vector<Point> empty_perturbs;

  // The number of components returned by this RBParametrizedFunction
  auto n_components = this->get_n_components();

  // We first loop over all mus and all n_points, calling side_evaluate()
  // for each and storing the results.  It is easier to first
  // pre-compute all the values before filling output, since, in the
  // case of multi-sample RBParameters, the ordering of the loops is a
  // bit complicated otherwise.
  std::vector<std::vector<std::vector<Number>>> all_evals(mus.size());
  for (auto mu_index : index_range(mus))
    {
      // Allocate enough space to store all points for the current mu
      all_evals[mu_index].resize(n_points);
      for (auto point_index : index_range(all_evals[mu_index]))
        {
          // Evaluate all samples for the current mu at the current interpolation point
          all_evals[mu_index][point_index] =
            this->side_evaluate(mus[mu_index],
                                v.all_xyz[point_index],
                                v.elem_ids[point_index],
                                v.side_indices[point_index],
                                v.qps[point_index],
                                v.sbd_ids[point_index],
                                v.boundary_ids[point_index],
                                requires_xyz_perturbations ? v.all_xyz_perturb[point_index] : empty_perturbs,
                                v.phi_i_qp[point_index]);

          // The vector returned by side_evaluate() should contain:
          // n_components * mus[mu_index].n_samples()
          // entries. That is, for multi-sample RBParameters objects,
          // the vector will be packed with entries as follows:
          // [sample0_component0, sample0_component1, ..., sample0_componentN,
          //  sample1_component0, sample1_component1, ..., sample1_componentN,
          //  ...
          //  sampleM_component0, sampleM_component1, ..., sampleM_componentN]
          auto n_samples = mus[mu_index].n_samples();
          auto received_data = all_evals[mu_index][point_index].size();
          libmesh_error_msg_if(received_data != n_components * n_samples,
                               "Recieved " << received_data <<
                               " evaluated values but expected to receive " << n_components * n_samples);
        }
    }

  // TODO: move this code for computing the total number of samples
  // represented by a std::vector of RBParameters objects to a helper
  // function.
  unsigned int output_size = 0;
  for (const auto & mu : mus)
    output_size += mu.n_samples();

  output.resize(output_size);

  // We use traditional for-loops here (rather than range-based) so that we can declare and
  // increment multiple loop counters all within the local scope of the for-loop.
  for (auto [mu_index, output_index] = std::make_tuple(0u, 0u); mu_index < mus.size(); ++mu_index)
    {
      auto n_samples = mus[mu_index].n_samples();
      for (auto mu_sample_idx = 0u; mu_sample_idx < n_samples; ++mu_sample_idx, ++output_index)
        {
          output[output_index].resize(n_points);
          for (auto point_index : make_range(n_points))
            {
              output[output_index][point_index].resize(n_components);

              for (auto comp : make_range(n_components))
                output[output_index][point_index][comp] = all_evals[mu_index][point_index][n_components*mu_sample_idx + comp];
            }
        }
    }
}

void RBParametrizedFunction::node_vectorized_evaluate(const std::vector<RBParameters> & mus,
                                                      const VectorizedEvalInput & v,
                                                      std::vector<std::vector<std::vector<Number>>> & output)
{
  LOG_SCOPE("node_vectorized_evaluate()", "RBParametrizedFunction");

  output.clear();
  unsigned int n_points = v.all_xyz.size();

  // The number of components returned by this RBParametrizedFunction
  auto n_components = this->get_n_components();

  // We first loop over all mus and all n_points, calling node_evaluate()
  // for each and storing the results.  It is easier to first
  // pre-compute all the values before filling output, since, in the
  // case of multi-sample RBParameters, the ordering of the loops is a
  // bit complicated otherwise.
  std::vector<std::vector<std::vector<Number>>> all_evals(mus.size());
  for (auto mu_index : index_range(mus))
    {
      // Allocate enough space to store all points for the current mu
      all_evals[mu_index].resize(n_points);
      for (auto point_index : index_range(all_evals[mu_index]))
        {
          // Evaluate all samples for the current mu at the current interpolation point
          all_evals[mu_index][point_index] =
            this->node_evaluate(mus[mu_index],
                                v.all_xyz[point_index],
                                v.node_ids[point_index],
                                v.boundary_ids[point_index]);

          // The vector returned by node_evaluate() should contain:
          // n_components * mus[mu_index].n_samples()
          // entries. That is, for multi-sample RBParameters objects,
          // the vector will be packed with entries as follows:
          // [sample0_component0, sample0_component1, ..., sample0_componentN,
          //  sample1_component0, sample1_component1, ..., sample1_componentN,
          //  ...
          //  sampleM_component0, sampleM_component1, ..., sampleM_componentN]
          auto n_samples = mus[mu_index].n_samples();
          auto received_data = all_evals[mu_index][point_index].size();
          libmesh_error_msg_if(received_data != n_components * n_samples,
                               "Recieved " << received_data <<
                               " evaluated values but expected to receive " << n_components * n_samples);
        }
    }

  // TODO: move this code for computing the total number of samples
  // represented by a std::vector of RBParameters objects to a helper
  // function.
  unsigned int output_size = 0;
  for (const auto & mu : mus)
    output_size += mu.n_samples();

  output.resize(output_size);

  // We use traditional for-loops here (rather than range-based) so that we can declare and
  // increment multiple loop counters all within the local scope of the for-loop.
  for (auto [mu_index, output_index] = std::make_tuple(0u, 0u); mu_index < mus.size(); ++mu_index)
    {
      auto n_samples = mus[mu_index].n_samples();
      for (auto mu_sample_idx = 0u; mu_sample_idx < n_samples; ++mu_sample_idx, ++output_index)
        {
          output[output_index].resize(n_points);
          for (auto point_index : make_range(n_points))
            {
              output[output_index][point_index].resize(n_components);

              for (auto comp : make_range(n_components))
                output[output_index][point_index][comp] = all_evals[mu_index][point_index][n_components*mu_sample_idx + comp];
            }
        }
    }
}

void RBParametrizedFunction::preevaluate_parametrized_function_on_mesh(const RBParameters & mu,
                                                                       const std::unordered_map<dof_id_type, std::vector<Point>> & all_xyz,
                                                                       const std::unordered_map<dof_id_type, subdomain_id_type> & sbd_ids,
                                                                       const std::unordered_map<dof_id_type, std::vector<std::vector<Point>> > & all_xyz_perturb,
                                                                       const System & sys)
{
  mesh_to_preevaluated_values_map.clear();

  unsigned int n_elems = all_xyz.size();
  unsigned int n_points = 0;
  for (const auto & xyz_pair : all_xyz)
  {
    const std::vector<Point> & xyz_vec = xyz_pair.second;
    n_points += xyz_vec.size();
  }

  VectorizedEvalInput v;
  v.all_xyz.resize(n_points);
  v.elem_ids.resize(n_points);
  v.qps.resize(n_points);
  v.sbd_ids.resize(n_points);
  v.all_xyz_perturb.resize(n_points);
  v.phi_i_qp.resize(n_points);
  v.elem_types.resize(n_points);

  if (requires_all_elem_qp_data)
    {
      v.elem_id_to_local_index.clear();
      v.JxW_all_qp.resize(n_elems);
      v.phi_i_all_qp.resize(n_elems);
    }
  if (requires_all_elem_center_data)
    {
      v.dxyzdxi_elem_center.resize(n_elems);
      v.dxyzdeta_elem_center.resize(n_elems);
      // At the time of writing, the quadrature order is only used in conjunction with element
      // center data so we should not compute it elsewhere for now.
      v.qrule_orders.resize(n_elems);
    }

  // Empty vector to be used when xyz perturbations are not required
  std::vector<Point> empty_perturbs;

  // In order to compute phi_i_qp, we initialize a FEMContext
  FEMContext con(sys);
  for (auto dim : con.elem_dimensions())
    {
      auto fe = con.get_element_fe(/*var=*/0, dim);
      fe->get_JxW();
      fe->get_phi();
      fe->get_dxyzdxi();
      fe->get_dxyzdeta();
    }

  unsigned int counter = 0;
  unsigned int elem_counter = 0;
  for (const auto & [elem_id, xyz_vec] : all_xyz)
    {
      subdomain_id_type subdomain_id = libmesh_map_find(sbd_ids, elem_id);

      // The amount of data to be stored for each component
      auto n_qp = xyz_vec.size();
      mesh_to_preevaluated_values_map[elem_id].resize(n_qp);

      // Also initialize phi in order to compute phi_i_qp
      const Elem & elem_ref = sys.get_mesh().elem_ref(elem_id);
      con.pre_fe_reinit(sys, &elem_ref);

      auto elem_fe = con.get_element_fe(/*var=*/0, elem_ref.dim());
      const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();
      const std::vector<Real> & JxW = elem_fe->get_JxW();
      const auto & dxyzdxi = elem_fe->get_dxyzdxi();
      const auto & dxyzdeta = elem_fe->get_dxyzdeta();

      elem_fe->reinit(&elem_ref);

      for (auto qp : index_range(xyz_vec))
        {
          mesh_to_preevaluated_values_map[elem_id][qp] = counter;

          v.all_xyz[counter] = xyz_vec[qp];
          v.elem_ids[counter] = elem_id;
          v.qps[counter] = qp;
          v.sbd_ids[counter] = subdomain_id;
          v.elem_types[counter] = elem_ref.type();

          v.phi_i_qp[counter].resize(phi.size());
          for(auto i : index_range(phi))
            v.phi_i_qp[counter][i] = phi[i][qp];

          if (requires_xyz_perturbations)
            {
              const auto & qps_and_perturbs =
                libmesh_map_find(all_xyz_perturb, elem_id);
              libmesh_error_msg_if(qp >= qps_and_perturbs.size(), "Error: Invalid qp");

              v.all_xyz_perturb[counter] = qps_and_perturbs[qp];
            }
          else
            {
              v.all_xyz_perturb[counter] = empty_perturbs;
            }

          if (requires_all_elem_qp_data)
            {
              if (v.elem_id_to_local_index.count(elem_id) == 0)
                {
                  // In this case we store data for all qps on this element
                  // at each point.
                  v.JxW_all_qp[elem_counter].resize(JxW.size());
                  for(auto i : index_range(JxW))
                    v.JxW_all_qp[elem_counter][i] = JxW[i];

                  v.phi_i_all_qp[elem_counter].resize(phi.size());
                  for(auto i : index_range(phi))
                  {
                    v.phi_i_all_qp[elem_counter][i].resize(phi[i].size());
                    for(auto j : index_range(phi[i]))
                      v.phi_i_all_qp[elem_counter][i][j] = phi[i][j];
                  }
                  v.elem_id_to_local_index[elem_id] = elem_counter;
                }
            }

          counter++;
        }

      // Here we presume that if requires_all_elem_center_data is set to true
      // then requires_all_elem_qp_data is also set to true so that elem_id_to_local_index entry
      // for this elem has already been set.
      if (requires_all_elem_center_data)
        {
          // Get data derivatives at vertex average
          std::vector<Point> nodes = { elem_ref.reference_elem()->vertex_average() };
          elem_fe->reinit (&elem_ref, &nodes);
          // Set qrule_order here to prevent calling getter multiple times.
          Order qrule_order = con.get_element_qrule().get_order();

          // We add qrule_order in this loop as it is used in conjunction with elem center
          // quantities for now.
          v.qrule_orders[elem_counter] = qrule_order;
          Point dxyzdxi_pt, dxyzdeta_pt;
          if (con.get_elem_dim()>0)
            dxyzdxi_pt = dxyzdxi[0];
          if (con.get_elem_dim()>1)
            dxyzdeta_pt = dxyzdeta[0];
          // Here we do an implicit conversion from RealGradient which is a VectorValue<Real>
          // which in turn is a TypeVector<T> to a Point which is a TypeVector<Real>.
          // They are essentially the same thing. This helps us limiting the number of includes
          // in serialization and deserialization as RealGradient is a typedef and we cannot
          // forward declare typedefs. As a result we leverage the fact that point.h is already
          // included in most places we need RealGradient.
          v.dxyzdxi_elem_center[elem_counter] = dxyzdxi_pt;
          v.dxyzdeta_elem_center[elem_counter] = dxyzdeta_pt;
        }
      elem_counter++;
    }

  std::vector<RBParameters> mus {mu};
  vectorized_evaluate(mus, v, preevaluated_values);

  preevaluate_parametrized_function_cleanup();
}

void RBParametrizedFunction::preevaluate_parametrized_function_on_mesh_sides(const RBParameters & mu,
                                                                             const std::map<std::pair<dof_id_type,unsigned int>, std::vector<Point>> & side_all_xyz,
                                                                             const std::map<std::pair<dof_id_type,unsigned int>, subdomain_id_type> & sbd_ids,
                                                                             const std::map<std::pair<dof_id_type,unsigned int>, boundary_id_type> & side_boundary_ids,
                                                                             const std::map<std::pair<dof_id_type,unsigned int>, unsigned int> & side_types,
                                                                             const std::map<std::pair<dof_id_type,unsigned int>, std::vector<std::vector<Point>> > & side_all_xyz_perturb,
                                                                             const System & sys)
{
  mesh_to_preevaluated_side_values_map.clear();

  unsigned int n_points = 0;
  for (const auto & xyz_pair : side_all_xyz)
  {
    const std::vector<Point> & xyz_vec = xyz_pair.second;
    n_points += xyz_vec.size();
  }

  VectorizedEvalInput v;
  v.all_xyz.resize(n_points);
  v.elem_ids.resize(n_points);
  v.side_indices.resize(n_points);
  v.qps.resize(n_points);
  v.sbd_ids.resize(n_points);
  v.boundary_ids.resize(n_points);
  v.all_xyz_perturb.resize(n_points);
  v.phi_i_qp.resize(n_points);

  // Empty vector to be used when xyz perturbations are not required
  std::vector<Point> empty_perturbs;

  // In order to compute phi_i_qp, we initialize a FEMContext
  FEMContext con(sys);
  for (auto dim : con.elem_dimensions())
    {
      auto fe = con.get_element_fe(/*var=*/0, dim);
      fe->get_phi();

      auto side_fe = con.get_side_fe(/*var=*/0, dim);
      side_fe->get_phi();
    }

  unsigned int counter = 0;
  for (const auto & xyz_pair : side_all_xyz)
    {
      auto elem_side_pair = xyz_pair.first;
      dof_id_type elem_id = elem_side_pair.first;
      unsigned int side_index = elem_side_pair.second;

      const std::vector<Point> & xyz_vec = xyz_pair.second;

      subdomain_id_type subdomain_id = libmesh_map_find(sbd_ids, elem_side_pair);
      boundary_id_type boundary_id = libmesh_map_find(side_boundary_ids, elem_side_pair);
      unsigned int side_type = libmesh_map_find(side_types, elem_side_pair);

      // The amount of data to be stored for each component
      auto n_qp = xyz_vec.size();
      mesh_to_preevaluated_side_values_map[elem_side_pair].resize(n_qp);

      const Elem & elem_ref = sys.get_mesh().elem_ref(elem_id);
      con.pre_fe_reinit(sys, &elem_ref);

      // side_type == 0 --> standard side
      // side_type == 1 --> shellface
      if (side_type == 0)
        {
          std::unique_ptr<const Elem> elem_side;
          elem_ref.build_side_ptr(elem_side, side_index);

          auto side_fe = con.get_side_fe(/*var=*/0, elem_ref.dim());
          side_fe->reinit(&elem_ref, side_index);

          const std::vector<std::vector<Real>> & phi = side_fe->get_phi();
          for (auto qp : index_range(xyz_vec))
            {
              mesh_to_preevaluated_side_values_map[elem_side_pair][qp] = counter;

              v.all_xyz[counter] = xyz_vec[qp];
              v.elem_ids[counter] = elem_side_pair.first;
              v.side_indices[counter] = elem_side_pair.second;
              v.qps[counter] = qp;
              v.sbd_ids[counter] = subdomain_id;
              v.boundary_ids[counter] = boundary_id;

              v.phi_i_qp[counter].resize(phi.size());
              for(auto i : index_range(phi))
                v.phi_i_qp[counter][i] = phi[i][qp];

              if (requires_xyz_perturbations)
                {
                  const auto & qps_and_perturbs =
                    libmesh_map_find(side_all_xyz_perturb, elem_side_pair);
                  libmesh_error_msg_if(qp >= qps_and_perturbs.size(), "Error: Invalid qp");

                  v.all_xyz_perturb[counter] = qps_and_perturbs[qp];
                }
              else
                {
                  v.all_xyz_perturb[counter] = empty_perturbs;
                }
              counter++;
            }
        }
      else if (side_type == 1)
        {
          auto elem_fe = con.get_element_fe(/*var=*/0, elem_ref.dim());
          const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();

          elem_fe->reinit(&elem_ref);

          for (auto qp : index_range(xyz_vec))
            {
              mesh_to_preevaluated_side_values_map[elem_side_pair][qp] = counter;

              v.all_xyz[counter] = xyz_vec[qp];
              v.elem_ids[counter] = elem_side_pair.first;
              v.side_indices[counter] = elem_side_pair.second;
              v.qps[counter] = qp;
              v.sbd_ids[counter] = subdomain_id;
              v.boundary_ids[counter] = boundary_id;

              v.phi_i_qp[counter].resize(phi.size());
              for(auto i : index_range(phi))
                v.phi_i_qp[counter][i] = phi[i][qp];

              if (requires_xyz_perturbations)
                {
                  const auto & qps_and_perturbs =
                    libmesh_map_find(side_all_xyz_perturb, elem_side_pair);
                  libmesh_error_msg_if(qp >= qps_and_perturbs.size(), "Error: Invalid qp");

                  v.all_xyz_perturb[counter] = qps_and_perturbs[qp];
                }
              else
                {
                  v.all_xyz_perturb[counter] = empty_perturbs;
                }
              counter++;
            }
        }
      else
        libmesh_error_msg ("Unrecognized side_type: " << side_type);
    }

  std::vector<RBParameters> mus {mu};
  side_vectorized_evaluate(mus, v, preevaluated_values);

  preevaluate_parametrized_function_cleanup();
}

void RBParametrizedFunction::preevaluate_parametrized_function_on_mesh_nodes(const RBParameters & mu,
                                                                             const std::unordered_map<dof_id_type, Point> & all_xyz,
                                                                             const std::unordered_map<dof_id_type, boundary_id_type> & node_boundary_ids,
                                                                             const System & /*sys*/)
{
  mesh_to_preevaluated_node_values_map.clear();

  unsigned int n_points = all_xyz.size();

  VectorizedEvalInput v;
  v.all_xyz.resize(n_points);
  v.node_ids.resize(n_points);
  v.boundary_ids.resize(n_points);

  // Empty vector to be used when xyz perturbations are not required
  std::vector<Point> empty_perturbs;

  unsigned int counter = 0;
  for (const auto & [node_id, p] : all_xyz)
    {
      boundary_id_type boundary_id = libmesh_map_find(node_boundary_ids, node_id);

      mesh_to_preevaluated_node_values_map[node_id] = counter;

      v.all_xyz[counter] = p;
      v.node_ids[counter] = node_id;
      v.boundary_ids[counter] = boundary_id;

      counter++;
    }

  std::vector<RBParameters> mus {mu};
  node_vectorized_evaluate(mus, v, preevaluated_values);

  preevaluate_parametrized_function_cleanup();
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

Number RBParametrizedFunction::lookup_preevaluated_side_value_on_mesh(unsigned int comp,
                                                                      dof_id_type elem_id,
                                                                      unsigned int side_index,
                                                                      unsigned int qp) const
{
  const std::vector<unsigned int> & indices_at_qps =
    libmesh_map_find(mesh_to_preevaluated_side_values_map, std::make_pair(elem_id,side_index));

  libmesh_error_msg_if(qp >= indices_at_qps.size(), "Error: invalid qp");

  unsigned int index = indices_at_qps[qp];
  libmesh_error_msg_if(preevaluated_values.size() != 1, "Error: we expect only one parameter index");
  libmesh_error_msg_if(index >= preevaluated_values[0].size(), "Error: invalid index");

  return preevaluated_values[0][index][comp];
}

Number RBParametrizedFunction::lookup_preevaluated_node_value_on_mesh(unsigned int comp,
                                                                      dof_id_type node_id) const
{
  unsigned int index =
    libmesh_map_find(mesh_to_preevaluated_node_values_map, node_id);

  libmesh_error_msg_if(preevaluated_values.size() != 1, "Error: we expect only one parameter index");
  libmesh_error_msg_if(index >= preevaluated_values[0].size(), "Error: invalid index");

  return preevaluated_values[0][index][comp];
}

void RBParametrizedFunction::initialize_lookup_table()
{
  // No-op by default, override in subclasses as needed
}

Number RBParametrizedFunction::get_parameter_independent_data(const std::string & property_name,
                                                              subdomain_id_type sbd_id) const
{
  return libmesh_map_find(libmesh_map_find(_parameter_independent_data, property_name), sbd_id);
}

void RBParametrizedFunction::get_spatial_indices(std::vector<std::vector<unsigned int>> & /*spatial_indices*/,
                                                 const VectorizedEvalInput & /*v*/)
{
  // No-op by default
}

void RBParametrizedFunction::initialize_spatial_indices(const std::vector<std::vector<unsigned int>> & /*spatial_indices*/,
                                                        const VectorizedEvalInput & /*v*/)
{
  // No-op by default
}

void RBParametrizedFunction::preevaluate_parametrized_function_cleanup()
{
  // No-op by default
}

const std::set<boundary_id_type> & RBParametrizedFunction::get_parametrized_function_boundary_ids() const
{
  return _parametrized_function_boundary_ids;
}

void RBParametrizedFunction::set_parametrized_function_boundary_ids(const std::set<boundary_id_type> & boundary_ids, bool is_nodal_boundary)
{
  _parametrized_function_boundary_ids = boundary_ids;
  _is_nodal_boundary = is_nodal_boundary;
}

bool RBParametrizedFunction::on_mesh_sides() const
{
  return !get_parametrized_function_boundary_ids().empty() && !_is_nodal_boundary;
}

bool RBParametrizedFunction::on_mesh_nodes() const
{
  return !get_parametrized_function_boundary_ids().empty() && _is_nodal_boundary;
}

}
