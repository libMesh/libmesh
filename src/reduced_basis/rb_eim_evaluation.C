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
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_eim_theta.h"
#include "libmesh/rb_parametrized_function.h"
#include "libmesh/rb_evaluation.h"
#include "libmesh/utility.h" // Utility::mkdir

// libMesh includes
#include "libmesh/xdr_cxx.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"
#include "timpi/parallel_implementation.h"

// C++ includes
#include <sstream>
#include <fstream>
#include <numeric> // std::accumulate
#include <iterator> // std::advance

namespace libMesh
{

RBEIMEvaluation::RBEIMEvaluation(const Parallel::Communicator & comm)
:
ParallelObject(comm),
eim_error_indicator_normalization(RBEIMEvaluation::MAX_RHS),
limit_eim_error_indicator_to_one(true),
_rb_eim_solves_N(0),
_preserve_rb_eim_solutions(false),
_is_eim_error_indicator_active(false)
{
}

RBEIMEvaluation::~RBEIMEvaluation() = default;

void RBEIMEvaluation::clear()
{
  _vec_eval_input.clear();

  _interpolation_points_comp.clear();
  _interpolation_points_spatial_indices.clear();


  _interpolation_matrix.resize(0,0);

  // Delete any RBTheta objects that were created
  _rb_eim_theta_objects.clear();
}

void RBEIMEvaluation::resize_data_structures(const unsigned int Nmax)
{
  // Resize the data structures relevant to the EIM system
  _vec_eval_input.clear();
  _vec_eval_input.all_xyz.clear();

  _interpolation_points_comp.clear();
  _interpolation_points_spatial_indices.clear();

  _interpolation_matrix.resize(Nmax,Nmax);
}

void RBEIMEvaluation::set_parametrized_function(std::unique_ptr<RBParametrizedFunction> pf)
{
  _parametrized_function = std::move(pf);
}

RBParametrizedFunction & RBEIMEvaluation::get_parametrized_function()
{
  libmesh_error_msg_if(!_parametrized_function, "Parametrized function not initialized yet");

  return *_parametrized_function;
}

const RBParametrizedFunction & RBEIMEvaluation::get_parametrized_function() const
{
  libmesh_error_msg_if(!_parametrized_function, "Parametrized function not initialized yet");

  return *_parametrized_function;
}

DenseVector<Number> RBEIMEvaluation::rb_eim_solve(DenseVector<Number> & EIM_rhs)
{
  LOG_SCOPE("rb_eim_solve()", "RBEIMEvaluation");

  libmesh_error_msg_if(EIM_rhs.size() > get_n_basis_functions(),
                       "Error: N cannot be larger than the number of basis functions in rb_solve");

  libmesh_error_msg_if(EIM_rhs.size()==0, "Error: N must be greater than 0 in rb_solve");

  const unsigned int N = EIM_rhs.size();
  DenseVector<Number> rb_eim_solution(N);
  DenseMatrix<Number> interpolation_matrix_N;
  _interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  interpolation_matrix_N.lu_solve(EIM_rhs, rb_eim_solution);

  return rb_eim_solution;
}

void RBEIMEvaluation::rb_eim_solves(const std::vector<RBParameters> & mus,
                                    unsigned int N)
{
  if (_preserve_rb_eim_solutions)
    {
      // In this case we preserve _rb_eim_solutions and hence we
      // just return immediately so that we skip updating
      // _rb_eim_solutions below. This is relevant in cases where
      // we set up _rb_eim_solutions elsewhere and we don't want
      // to override it.
      return;
    }

  libmesh_error_msg_if(N > get_n_basis_functions(),
    "Error: N cannot be larger than the number of basis functions in rb_eim_solves");
  libmesh_error_msg_if(N==0, "Error: N must be greater than 0 in rb_eim_solves");

  // If mus and N are the same as before, then we return early
  if ((_rb_eim_solves_mus == mus) && (_rb_eim_solves_N == N))
    return;

  LOG_SCOPE("rb_eim_solves()", "RBEIMEvaluation");

  _rb_eim_solves_mus = mus;
  _rb_eim_solves_N = N;

  if (get_parametrized_function().is_lookup_table)
    {
      _rb_eim_solutions.resize(mus.size());
      for (auto mu_index : index_range(mus))
        {
          Real lookup_table_param =
            mus[mu_index].get_value(get_parametrized_function().lookup_table_param_name);

          // Cast lookup_table_param to an unsigned integer so that we can use
          // it as an index into the EIM rhs values obtained from the lookup table.
          unsigned int lookup_table_index =
            cast_int<unsigned int>(std::round(lookup_table_param));

          DenseVector<Number> values;
          _eim_solutions_for_training_set[lookup_table_index].get_principal_subvector(N, values);
          _rb_eim_solutions[mu_index] = values;
        }

      return;
    }

  // output all comps indexing is as follows:
  //   mu index --> interpolation point index --> component index --> value.
  std::vector<std::vector<std::vector<Number>>> output_all_comps;
  if (get_parametrized_function().on_mesh_sides())
    get_parametrized_function().side_vectorized_evaluate(mus, _vec_eval_input, output_all_comps);
  else if (get_parametrized_function().on_mesh_nodes())
    get_parametrized_function().node_vectorized_evaluate(mus, _vec_eval_input, output_all_comps);
  else
    get_parametrized_function().vectorized_evaluate(mus, _vec_eval_input, output_all_comps);

  // Previously we did one RB-EIM solve per input mu, but now we do
  // one RB-EIM solve per input mu, per sample. In order for this to
  // work, we require that all the input mu objects have the same
  // number of samples.
  auto n_samples_0 = mus[0].n_samples();
  for (const auto & mu : mus)
    libmesh_error_msg_if(mu.n_samples() != n_samples_0, "All RBParameters objects must have same n_samples()");

  // After we verified that all mus have the same number of samples,
  // the total number of RB-EIM solves is simply the number of mus
  // times the number of samples.
  unsigned int num_rb_eim_solves = mus.size() * n_samples_0;

  // A special case is when we are passed a single RBParameters object
  // with no parameters stored on it. In this case, we effectively
  // have Theta(mu) == const, and therefore we still need to do at
  // least one RB-EIM solve. In this case, we require that there is
  // only one entry in "mus" for simplicity.
  if (num_rb_eim_solves == 0 && mus[0].n_parameters() == 0)
    {
      libmesh_error_msg_if(mus.size() != 1, "Must pass in only a single RBParameters object when solving with no parameters.");
      num_rb_eim_solves = 1;
    }

  std::vector<std::vector<Number>> evaluated_values_at_interp_points(num_rb_eim_solves);

  std::vector<Number> evaluated_values_at_err_indicator_point;
  if (_is_eim_error_indicator_active)
    evaluated_values_at_err_indicator_point.resize(num_rb_eim_solves);

  // In this loop, counter goes from 0 to num_rb_eim_solves.  The
  // purpose of this loop is to strip out the "columns" of the
  // output_all_comps array into rows.
  {
    unsigned int counter = 0;
    for (auto mu_index : index_range(mus))
      for (auto sample_index : make_range(mus[mu_index].n_samples()))
      {
        // Ignore compiler warnings about unused loop index
        libmesh_ignore(sample_index);

        evaluated_values_at_interp_points[counter].resize(N);

        for (unsigned int interp_pt_index=0; interp_pt_index<N; interp_pt_index++)
          {
            unsigned int comp = _interpolation_points_comp[interp_pt_index];

            // This line of code previously used "mu_index", now we use
            // "counter" handle the multi-sample RBParameters case.
            evaluated_values_at_interp_points[counter][interp_pt_index] =
              output_all_comps[counter][interp_pt_index][comp];
          }

        if (_is_eim_error_indicator_active)
          {
            unsigned int comp = _interpolation_points_comp[N];

            evaluated_values_at_err_indicator_point[counter] =
              output_all_comps[counter][N][comp];
          }

        counter++;
      }

    // Throw an error if we didn't do the required number of solves for
    // some reason
    libmesh_error_msg_if(counter != num_rb_eim_solves,
                        "We should have done " << num_rb_eim_solves <<
                        " solves, instead we did " << counter);
  }

  DenseMatrix<Number> interpolation_matrix_N;
  _interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  // The number of RB EIM solutions is equal to the size of the
  // "evaluated_values_at_interp_points" vector which we determined
  // earlier.
  _rb_eim_solutions.resize(num_rb_eim_solves);
  if (_is_eim_error_indicator_active)
    _rb_eim_error_indicators.resize(num_rb_eim_solves);

  {
    unsigned int counter = 0;
    for (auto mu_index : index_range(mus))
      for (auto sample_index : make_range(mus[mu_index].n_samples()))
      {
        // Ignore compiler warnings about unused loop index
        libmesh_ignore(sample_index);

        DenseVector<Number> EIM_rhs = evaluated_values_at_interp_points[counter];
        interpolation_matrix_N.lu_solve(EIM_rhs, _rb_eim_solutions[counter]);

        // If we're using the EIM error indicator, then we compute it via the approach
        // proposed in Proposition 3.3 of "An empirical interpolation method: application
        // to efficient reduced-basis discretization of partial differential equations",
        // Barrault et al.
        if (_is_eim_error_indicator_active)
          {
            Number error_indicator_rhs = evaluated_values_at_err_indicator_point[counter];
            _rb_eim_error_indicators[counter] =
              get_eim_error_indicator(
                error_indicator_rhs, _rb_eim_solutions[counter], EIM_rhs);
          }

        counter++;
      }
  }
}

void RBEIMEvaluation::initialize_interpolation_points_spatial_indices()
{
  _interpolation_points_spatial_indices.clear();

  get_parametrized_function().get_spatial_indices(_interpolation_points_spatial_indices,
                                                  _vec_eval_input);
}

void RBEIMEvaluation::initialize_param_fn_spatial_indices()
{
  get_parametrized_function().initialize_spatial_indices(_interpolation_points_spatial_indices,
                                                         _vec_eval_input);
}

unsigned int RBEIMEvaluation::get_n_basis_functions() const
{
  if (get_parametrized_function().on_mesh_sides())
    return _local_side_eim_basis_functions.size();
  else if (get_parametrized_function().on_mesh_nodes())
    return _local_node_eim_basis_functions.size();
  else
    return _local_eim_basis_functions.size();
}

unsigned int RBEIMEvaluation::get_n_interpolation_points() const
{
  return _vec_eval_input.all_xyz.size();
}

unsigned int RBEIMEvaluation::get_n_elems() const
{
  return _vec_eval_input.elem_id_to_local_index.size();
}

void RBEIMEvaluation::set_n_basis_functions(unsigned int n_bfs)
{
  if (get_parametrized_function().on_mesh_sides())
    _local_side_eim_basis_functions.resize(n_bfs);
  else if (get_parametrized_function().on_mesh_nodes())
    _local_node_eim_basis_functions.resize(n_bfs);
  else
    _local_eim_basis_functions.resize(n_bfs);
}

void RBEIMEvaluation::decrement_vector(QpDataMap & v,
                                       const DenseVector<Number> & coeffs)
{
  LOG_SCOPE("decrement_vector()", "RBEIMEvaluation");

  libmesh_error_msg_if(get_n_basis_functions() != coeffs.size(),
                       "Error: Number of coefficients should match number of basis functions");

  for (auto & [elem_id, v_comp_and_qp] : v)
    {
      for (const auto & comp : index_range(v_comp_and_qp))
        for (unsigned int qp : index_range(v_comp_and_qp[comp]))
          for (unsigned int i : index_range(_local_eim_basis_functions))
            {
              // Check that entry (elem_id,comp,qp) exists in _local_eim_basis_functions so that
              // we get a clear error message if there is any missing data
              const auto & basis_comp_and_qp = libmesh_map_find(_local_eim_basis_functions[i], elem_id);

              libmesh_error_msg_if(comp >= basis_comp_and_qp.size(), "Error: Invalid comp");
              libmesh_error_msg_if(qp >= basis_comp_and_qp[comp].size(), "Error: Invalid qp");

              v_comp_and_qp[comp][qp] -= coeffs(i) * basis_comp_and_qp[comp][qp];
            }
    }
}

void RBEIMEvaluation::side_decrement_vector(SideQpDataMap & v,
                                            const DenseVector<Number> & coeffs)
{
  LOG_SCOPE("side_decrement_vector()", "RBEIMEvaluation");

  libmesh_error_msg_if(get_n_basis_functions() != coeffs.size(),
                       "Error: Number of coefficients should match number of basis functions");

  for (auto & [elem_and_side, v_comp_and_qp] : v)
    {
      for (const auto & comp : index_range(v_comp_and_qp))
        for (unsigned int qp : index_range(v_comp_and_qp[comp]))
          for (unsigned int i : index_range(_local_side_eim_basis_functions))
            {
              // Check that entry (elem_and_side,comp,qp) exists in _local_side_eim_basis_functions so that
              // we get a clear error message if there is any missing data
              const auto & basis_comp_and_qp = libmesh_map_find(_local_side_eim_basis_functions[i], elem_and_side);

              libmesh_error_msg_if(comp >= basis_comp_and_qp.size(), "Error: Invalid comp");
              libmesh_error_msg_if(qp >= basis_comp_and_qp[comp].size(), "Error: Invalid qp");

              v_comp_and_qp[comp][qp] -= coeffs(i) * basis_comp_and_qp[comp][qp];
            }
    }
}

void RBEIMEvaluation::node_decrement_vector(NodeDataMap & v,
                                            const DenseVector<Number> & coeffs)
{
  LOG_SCOPE("node_decrement_vector()", "RBEIMEvaluation");

  libmesh_error_msg_if(get_n_basis_functions() != coeffs.size(),
                       "Error: Number of coefficients should match number of basis functions");

  for (auto & [node_id, v_comps] : v)
    {
      for (const auto & comp : index_range(v_comps))
        for (unsigned int i : index_range(_local_node_eim_basis_functions))
          {
            // Check that entry (node_id,comp) exists in _local_node_eim_basis_functions so that
            // we get a clear error message if there is any missing data
            const auto & basis_comp = libmesh_map_find(_local_node_eim_basis_functions[i], node_id);

            libmesh_error_msg_if(comp >= basis_comp.size(), "Error: Invalid comp");

            v_comps[comp] -= coeffs(i) * basis_comp[comp];
          }
    }
}

void RBEIMEvaluation::initialize_eim_theta_objects()
{
  // Initialize the rb_theta objects that access the solution from this rb_eim_evaluation
  _rb_eim_theta_objects.clear();
  for (auto i : make_range(get_n_basis_functions()))
    _rb_eim_theta_objects.emplace_back(build_eim_theta(i));
}

std::vector<std::unique_ptr<RBTheta>> & RBEIMEvaluation::get_eim_theta_objects()
{
  return _rb_eim_theta_objects;
}

std::unique_ptr<RBTheta> RBEIMEvaluation::build_eim_theta(unsigned int index)
{
  return std::make_unique<RBEIMTheta>(*this, index);
}

void RBEIMEvaluation::get_parametrized_function_values_at_qps(
  const QpDataMap & pf,
  dof_id_type elem_id,
  unsigned int comp,
  std::vector<Number> & values)
{
  LOG_SCOPE("get_parametrized_function_values_at_qps()", "RBEIMConstruction");

  values.clear();

  if (const auto it = pf.find(elem_id);
      it != pf.end())
  {
    const auto & comps_and_qps_on_elem = it->second;
    libmesh_error_msg_if(comp >= comps_and_qps_on_elem.size(),
                         "Invalid comp index: " << comp);

    values = comps_and_qps_on_elem[comp];
  }
}

void RBEIMEvaluation::get_parametrized_function_side_values_at_qps(
  const SideQpDataMap & pf,
  dof_id_type elem_id,
  unsigned int side_index,
  unsigned int comp,
  std::vector<Number> & values)
{
  LOG_SCOPE("get_parametrized_function_side_values_at_qps()", "RBEIMConstruction");

  values.clear();

  if (const auto it = pf.find(std::make_pair(elem_id, side_index));
      it != pf.end())
  {
    const auto & comps_and_qps_on_elem = it->second;
    libmesh_error_msg_if(comp >= comps_and_qps_on_elem.size(),
                         "Invalid comp index: " << comp);

    values = comps_and_qps_on_elem[comp];
  }
}

Number RBEIMEvaluation::get_parametrized_function_node_local_value(
  const NodeDataMap & pf,
  dof_id_type node_id,
  unsigned int comp)
{
  LOG_SCOPE("get_parametrized_function_node_local_value()", "RBEIMConstruction");

  if (const auto it = pf.find(node_id);
      it != pf.end())
    {
      const std::vector<Number> & vec = it->second;
      libmesh_error_msg_if (comp >= vec.size(), "Error: Invalid comp index");
      return vec[comp];
    }
  else
    return 0.;
}

Number RBEIMEvaluation::get_parametrized_function_value(
  const Parallel::Communicator & comm,
  const QpDataMap & pf,
  dof_id_type elem_id,
  unsigned int comp,
  unsigned int qp)
{
  std::vector<Number> values;
  get_parametrized_function_values_at_qps(pf, elem_id, comp, values);

  // In parallel, values should only be non-empty on one processor
  Number value = 0.;
  if (!values.empty())
  {
    libmesh_error_msg_if(qp >= values.size(), "Error: Invalid qp index");

    value = values[qp];
  }
  comm.sum(value);

  return value;
}

Number RBEIMEvaluation::get_parametrized_function_side_value(
  const Parallel::Communicator & comm,
  const SideQpDataMap & pf,
  dof_id_type elem_id,
  unsigned int side_index,
  unsigned int comp,
  unsigned int qp)
{
  std::vector<Number> values;
  get_parametrized_function_side_values_at_qps(pf, elem_id, side_index, comp, values);

  // In parallel, values should only be non-empty on one processor
  Number value = 0.;
  if (!values.empty())
  {
    libmesh_error_msg_if(qp >= values.size(), "Error: Invalid qp index");

    value = values[qp];
  }
  comm.sum(value);

  return value;
}

Number RBEIMEvaluation::get_parametrized_function_node_value(
  const Parallel::Communicator & comm,
  const NodeDataMap & pf,
  dof_id_type node_id,
  unsigned int comp)
{
  LOG_SCOPE("get_parametrized_function_node_value()", "RBEIMConstruction");

  Number value = get_parametrized_function_node_local_value(pf, node_id, comp);
  comm.sum(value);

  return value;
}

void RBEIMEvaluation::get_eim_basis_function_values_at_qps(unsigned int basis_function_index,
                                                           dof_id_type elem_id,
                                                           unsigned int comp,
                                                           std::vector<Number> & values) const
{
  libmesh_error_msg_if(basis_function_index >= _local_eim_basis_functions.size(),
                       "Invalid basis function index: " << basis_function_index);

  get_parametrized_function_values_at_qps(
    _local_eim_basis_functions[basis_function_index],
    elem_id,
    comp,
    values);
}

void RBEIMEvaluation::get_eim_basis_function_side_values_at_qps(unsigned int basis_function_index,
                                                                dof_id_type elem_id,
                                                                unsigned int side_index,
                                                                unsigned int comp,
                                                                std::vector<Number> & values) const
{
  libmesh_error_msg_if(basis_function_index >= _local_side_eim_basis_functions.size(),
                       "Invalid basis function index: " << basis_function_index);

  get_parametrized_function_side_values_at_qps(
    _local_side_eim_basis_functions[basis_function_index],
    elem_id,
    side_index,
    comp,
    values);
}

Number RBEIMEvaluation::get_eim_basis_function_node_local_value(unsigned int basis_function_index,
                                                                dof_id_type node_id,
                                                                unsigned int comp) const
{
  libmesh_error_msg_if(basis_function_index >= _local_node_eim_basis_functions.size(),
                       "Invalid basis function index: " << basis_function_index);

  return get_parametrized_function_node_local_value(
    _local_node_eim_basis_functions[basis_function_index],
    node_id,
    comp);
}

Number RBEIMEvaluation::get_eim_basis_function_node_value(unsigned int basis_function_index,
                                                          dof_id_type node_id,
                                                          unsigned int comp) const
{
  libmesh_error_msg_if(basis_function_index >= _local_node_eim_basis_functions.size(),
                       "Invalid basis function index: " << basis_function_index);

  return get_parametrized_function_node_value(
    comm(),
    _local_node_eim_basis_functions[basis_function_index],
    node_id,
    comp);
}

Number RBEIMEvaluation::get_eim_basis_function_value(unsigned int basis_function_index,
                                                     dof_id_type elem_id,
                                                     unsigned int comp,
                                                     unsigned int qp) const
{
  libmesh_error_msg_if(basis_function_index >= _local_eim_basis_functions.size(),
                       "Invalid basis function index: " << basis_function_index);

  return get_parametrized_function_value(
    comm(),
    _local_eim_basis_functions[basis_function_index],
    elem_id,
    comp,
    qp);
}

Number RBEIMEvaluation::get_eim_basis_function_side_value(unsigned int basis_function_index,
                                                          dof_id_type elem_id,
                                                          unsigned int side_index,
                                                          unsigned int comp,
                                                          unsigned int qp) const
{
  libmesh_error_msg_if(basis_function_index >= _local_side_eim_basis_functions.size(),
                       "Invalid side basis function index: " << basis_function_index);

  return get_parametrized_function_side_value(
    comm(),
    _local_side_eim_basis_functions[basis_function_index],
    elem_id,
    side_index,
    comp,
    qp);
}

const RBEIMEvaluation::QpDataMap &
RBEIMEvaluation::get_basis_function(unsigned int i) const
{
  return _local_eim_basis_functions[i];
}

const RBEIMEvaluation::SideQpDataMap &
RBEIMEvaluation::get_side_basis_function(unsigned int i) const
{
  return _local_side_eim_basis_functions[i];
}

const RBEIMEvaluation::NodeDataMap &
RBEIMEvaluation::get_node_basis_function(unsigned int i) const
{
  return _local_node_eim_basis_functions[i];
}

void RBEIMEvaluation::set_rb_eim_solutions(const std::vector<DenseVector<Number>> & rb_eim_solutions)
{
  _rb_eim_solutions = rb_eim_solutions;
}

const std::vector<DenseVector<Number>> & RBEIMEvaluation::get_rb_eim_solutions() const
{
  return _rb_eim_solutions;
}

std::vector<Number> RBEIMEvaluation::get_rb_eim_solutions_entries(unsigned int index) const
{
  LOG_SCOPE("get_rb_eim_solutions_entries()", "RBEIMEvaluation");

  std::vector<Number> rb_eim_solutions_entries(_rb_eim_solutions.size());
  for (unsigned int mu_index : index_range(_rb_eim_solutions))
    {
      libmesh_error_msg_if(index >= _rb_eim_solutions[mu_index].size(),
                           "Error: Requested solution index " << index <<
                           ", but only have " << _rb_eim_solutions[mu_index].size() << " entries.");
      rb_eim_solutions_entries[mu_index] = _rb_eim_solutions[mu_index](index);
    }

  return rb_eim_solutions_entries;
}

const std::vector<DenseVector<Number>> & RBEIMEvaluation::get_eim_solutions_for_training_set() const
{
  return _eim_solutions_for_training_set;
}

std::vector<DenseVector<Number>> & RBEIMEvaluation::get_eim_solutions_for_training_set()
{
  return _eim_solutions_for_training_set;
}

const std::vector<std::pair<Real,Real>> & RBEIMEvaluation::get_rb_eim_error_indicators() const
{
  return _rb_eim_error_indicators;
}

void RBEIMEvaluation::add_interpolation_points_xyz(Point p)
{
  _vec_eval_input.all_xyz.emplace_back(p);
}

void RBEIMEvaluation::add_interpolation_points_comp(unsigned int comp)
{
  _interpolation_points_comp.emplace_back(comp);
}

void RBEIMEvaluation::add_interpolation_points_subdomain_id(subdomain_id_type sbd_id)
{
  _vec_eval_input.sbd_ids.emplace_back(sbd_id);
}

void RBEIMEvaluation::add_interpolation_points_boundary_id(boundary_id_type b_id)
{
  _vec_eval_input.boundary_ids.emplace_back(b_id);
}

void RBEIMEvaluation::add_interpolation_points_xyz_perturbations(const std::vector<Point> & perturbs)
{
  _vec_eval_input.all_xyz_perturb.emplace_back(perturbs);
}

void RBEIMEvaluation::add_interpolation_points_elem_id(dof_id_type elem_id)
{
  _vec_eval_input.elem_ids.emplace_back(elem_id);
}

void RBEIMEvaluation::add_interpolation_points_side_index(unsigned int side_index)
{
  _vec_eval_input.side_indices.emplace_back(side_index);
}

void RBEIMEvaluation::add_interpolation_points_node_id(dof_id_type node_id)
{
  _vec_eval_input.node_ids.emplace_back(node_id);
}

void RBEIMEvaluation::add_interpolation_points_qp(unsigned int qp)
{
  _vec_eval_input.qps.emplace_back(qp);
}

void RBEIMEvaluation::add_interpolation_points_elem_type(ElemType elem_type)
{
  _vec_eval_input.elem_types.emplace_back(elem_type);
}

void RBEIMEvaluation::add_interpolation_points_JxW_all_qp(const std::vector<Real> & JxW_all_qp)
{
  _vec_eval_input.JxW_all_qp.emplace_back(JxW_all_qp);
}

void RBEIMEvaluation::add_interpolation_points_phi_i_all_qp(const std::vector<std::vector<Real>> & phi_i_all_qp)
{
  _vec_eval_input.phi_i_all_qp.emplace_back(phi_i_all_qp);
}

void RBEIMEvaluation::add_interpolation_points_phi_i_qp(const std::vector<Real> & phi_i_qp)
{
  _vec_eval_input.phi_i_qp.emplace_back(phi_i_qp);
}

void RBEIMEvaluation::add_elem_center_dxyzdxi(const Point & dxyzdxi)
{
  _vec_eval_input.dxyzdxi_elem_center.emplace_back(dxyzdxi);
}

void RBEIMEvaluation::add_elem_center_dxyzdeta(const Point & dxyzdeta)
{
  _vec_eval_input.dxyzdeta_elem_center.emplace_back(dxyzdeta);
}

void RBEIMEvaluation::add_interpolation_points_qrule_order(Order qrule_order)
{
  _vec_eval_input.qrule_orders.emplace_back(qrule_order);
}

void RBEIMEvaluation::add_interpolation_points_spatial_indices(const std::vector<unsigned int> & spatial_indices)
{
  _interpolation_points_spatial_indices.emplace_back(spatial_indices);
}

void RBEIMEvaluation::add_elem_id_local_index_map_entry(const dof_id_type & elem_id, const unsigned int local_index)
{
  libmesh_error_msg_if(_vec_eval_input.elem_id_to_local_index.count(elem_id) == 1, "Entry already added, duplicate detected.");

  _vec_eval_input.elem_id_to_local_index[elem_id] = local_index;
}

Point RBEIMEvaluation::get_interpolation_points_xyz(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.all_xyz.size(), "Error: Invalid index");

  return _vec_eval_input.all_xyz[index];
}

unsigned int RBEIMEvaluation::get_interpolation_points_comp(unsigned int index) const
{
  libmesh_error_msg_if(index >= _interpolation_points_comp.size(), "Error: Invalid index");

  return _interpolation_points_comp[index];
}

subdomain_id_type RBEIMEvaluation::get_interpolation_points_subdomain_id(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.sbd_ids.size(), "Error: Invalid index");

  return _vec_eval_input.sbd_ids[index];
}

boundary_id_type RBEIMEvaluation::get_interpolation_points_boundary_id(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.boundary_ids.size(), "Error: Invalid index");

  return _vec_eval_input.boundary_ids[index];
}

const std::vector<Point> & RBEIMEvaluation::get_interpolation_points_xyz_perturbations(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.all_xyz_perturb.size(), "Error: Invalid index");

  return _vec_eval_input.all_xyz_perturb[index];
}

dof_id_type RBEIMEvaluation::get_interpolation_points_elem_id(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.elem_ids.size(), "Error: Invalid index");

  return _vec_eval_input.elem_ids[index];
}

unsigned int RBEIMEvaluation::get_interpolation_points_side_index(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.side_indices.size(), "Error: Invalid index");

  return _vec_eval_input.side_indices[index];
}

dof_id_type RBEIMEvaluation::get_interpolation_points_node_id(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.node_ids.size(), "Error: Invalid index");

  return _vec_eval_input.node_ids[index];
}

unsigned int RBEIMEvaluation::get_interpolation_points_qp(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.qps.size(), "Error: Invalid index");

  return _vec_eval_input.qps[index];
}

ElemType RBEIMEvaluation::get_interpolation_points_elem_type(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.elem_types.size(), "Error: Invalid index");

  return _vec_eval_input.elem_types[index];
}

const std::vector<Real> & RBEIMEvaluation::get_interpolation_points_JxW_all_qp(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.JxW_all_qp.size(), "Error: Invalid index");

  return _vec_eval_input.JxW_all_qp[index];
}

const std::map<dof_id_type, unsigned int> & RBEIMEvaluation::get_elem_id_to_local_index_map() const
{
  return _vec_eval_input.elem_id_to_local_index;
}

const std::vector<std::vector<Real>> & RBEIMEvaluation::get_interpolation_points_phi_i_all_qp(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.phi_i_all_qp.size(), "Error: Invalid index");

  return _vec_eval_input.phi_i_all_qp[index];
}

const std::vector<Real> & RBEIMEvaluation::get_interpolation_points_phi_i_qp(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.phi_i_qp.size(), "Error: Invalid index");

  return _vec_eval_input.phi_i_qp[index];
}

const Point & RBEIMEvaluation::get_elem_center_dxyzdxi(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.dxyzdxi_elem_center.size(), "Error: Invalid index");

  return _vec_eval_input.dxyzdxi_elem_center[index];
}

const Point & RBEIMEvaluation::get_elem_center_dxyzdeta(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.dxyzdeta_elem_center.size(), "Error: Invalid index");

  return _vec_eval_input.dxyzdeta_elem_center[index];
}

Order RBEIMEvaluation::get_interpolation_points_qrule_order(unsigned int index) const
{
  libmesh_error_msg_if(index >= _vec_eval_input.qrule_orders.size(), "Error: Invalid index");

  return _vec_eval_input.qrule_orders[index];
}

const std::vector<unsigned int> & RBEIMEvaluation::get_interpolation_points_spatial_indices(unsigned int index) const
{
  libmesh_error_msg_if(index >= _interpolation_points_spatial_indices.size(), "Error: Invalid index");

  return _interpolation_points_spatial_indices[index];
}

unsigned int RBEIMEvaluation::get_n_interpolation_points_spatial_indices() const
{
  return _interpolation_points_spatial_indices.size();
}

void RBEIMEvaluation::set_interpolation_matrix_entry(unsigned int i, unsigned int j, Number value)
{
  libmesh_error_msg_if((i >= _interpolation_matrix.m()) || (j >= _interpolation_matrix.n()),
                       "Error: Invalid matrix indices");

  _interpolation_matrix(i,j) = value;
}

const DenseMatrix<Number> & RBEIMEvaluation::get_interpolation_matrix() const
{
  return _interpolation_matrix;
}

void RBEIMEvaluation::set_preserve_rb_eim_solutions(bool preserve_rb_eim_solutions)
{
  _preserve_rb_eim_solutions = preserve_rb_eim_solutions;
}

bool RBEIMEvaluation::get_preserve_rb_eim_solutions() const
{
  return _preserve_rb_eim_solutions;
}

void RBEIMEvaluation::add_basis_function(
  const QpDataMap & bf)
{
  _local_eim_basis_functions.emplace_back(bf);
}

void RBEIMEvaluation::add_interpolation_data(
  Point p,
  unsigned int comp,
  dof_id_type elem_id,
  subdomain_id_type subdomain_id,
  unsigned int qp,
  const std::vector<Point> & perturbs,
  const std::vector<Real> & phi_i_qp,
  ElemType elem_type,
  const std::vector<Real> & JxW_all_qp,
  const std::vector<std::vector<Real>> & phi_i_all_qp,
  Order qrule_order,
  const Point & dxyzdxi_elem_center,
  const Point & dxyzdeta_elem_center)
{
  _vec_eval_input.all_xyz.emplace_back(p);
  _interpolation_points_comp.emplace_back(comp);
  _vec_eval_input.elem_ids.emplace_back(elem_id);
  _vec_eval_input.sbd_ids.emplace_back(subdomain_id);
  _vec_eval_input.qps.emplace_back(qp);
  _vec_eval_input.all_xyz_perturb.emplace_back(perturbs);
  _vec_eval_input.phi_i_qp.emplace_back(phi_i_qp);
  _vec_eval_input.elem_types.emplace_back(elem_type);

  // The following quantities are indexed by elem id to prevent duplicated data
  // If an entry is already present then we should not need to add that data again as it is element
  // based so it should be identical between 2 points if they refer to the same element id.
  if (_vec_eval_input.elem_id_to_local_index.count(elem_id) == 0)
  {
    unsigned int local_index = _vec_eval_input.JxW_all_qp.size();
    // Maybe check that getting local index from a a different source is still ok.
    _vec_eval_input.elem_id_to_local_index[elem_id] = local_index;
    _vec_eval_input.JxW_all_qp.emplace_back(JxW_all_qp);
    _vec_eval_input.phi_i_all_qp.emplace_back(phi_i_all_qp);
    _vec_eval_input.qrule_orders.emplace_back(qrule_order);
    _vec_eval_input.dxyzdxi_elem_center.emplace_back(dxyzdxi_elem_center);
    _vec_eval_input.dxyzdeta_elem_center.emplace_back(dxyzdeta_elem_center);
  }
}

void RBEIMEvaluation::add_side_basis_function(
  const SideQpDataMap & side_bf)
{
  _local_side_eim_basis_functions.emplace_back(side_bf);
}

void RBEIMEvaluation::add_side_interpolation_data(
  Point p,
  unsigned int comp,
  dof_id_type elem_id,
  unsigned int side_index,
  subdomain_id_type subdomain_id,
  boundary_id_type boundary_id,
  unsigned int qp,
  const std::vector<Point> & perturbs,
  const std::vector<Real> & phi_i_qp)
{
  _vec_eval_input.all_xyz.emplace_back(p);
  _interpolation_points_comp.emplace_back(comp);
  _vec_eval_input.elem_ids.emplace_back(elem_id);
  _vec_eval_input.side_indices.emplace_back(side_index);
  _vec_eval_input.sbd_ids.emplace_back(subdomain_id);
  _vec_eval_input.boundary_ids.emplace_back(boundary_id);
  _vec_eval_input.qps.emplace_back(qp);
  _vec_eval_input.all_xyz_perturb.emplace_back(perturbs);
  _vec_eval_input.phi_i_qp.emplace_back(phi_i_qp);

  // Add dummy values for the other properties, which are unused in the
  // node case.
  _vec_eval_input.elem_types.emplace_back(INVALID_ELEM);
  _vec_eval_input.JxW_all_qp.emplace_back();
  _vec_eval_input.phi_i_all_qp.emplace_back();
  _vec_eval_input.qrule_orders.emplace_back(INVALID_ORDER);
  _vec_eval_input.dxyzdxi_elem_center.emplace_back();
  _vec_eval_input.dxyzdeta_elem_center.emplace_back();
}

void RBEIMEvaluation::add_node_basis_function(
  const NodeDataMap & node_bf)
{
  _local_node_eim_basis_functions.emplace_back(node_bf);
}

void RBEIMEvaluation::add_node_interpolation_data(
  Point p,
  unsigned int comp,
  dof_id_type node_id,
  boundary_id_type boundary_id)
{
  _vec_eval_input.all_xyz.emplace_back(p);
  _interpolation_points_comp.emplace_back(comp);
  _vec_eval_input.node_ids.emplace_back(node_id);
  _vec_eval_input.boundary_ids.emplace_back(boundary_id);

  // Add dummy values for the other properties, which are unused in the
  // node case.
  _vec_eval_input.elem_ids.emplace_back(0);
  _vec_eval_input.side_indices.emplace_back(0);
  _vec_eval_input.sbd_ids.emplace_back(0);
  _vec_eval_input.qps.emplace_back(0);
  _vec_eval_input.all_xyz_perturb.emplace_back();
  _vec_eval_input.phi_i_qp.emplace_back();
  _vec_eval_input.elem_types.emplace_back(INVALID_ELEM);
  _vec_eval_input.JxW_all_qp.emplace_back();
  _vec_eval_input.phi_i_all_qp.emplace_back();
  _vec_eval_input.qrule_orders.emplace_back(INVALID_ORDER);
  _vec_eval_input.dxyzdxi_elem_center.emplace_back();
  _vec_eval_input.dxyzdeta_elem_center.emplace_back();
}

void RBEIMEvaluation::
write_out_basis_functions(const std::string & directory_name,
                          bool write_binary_basis_functions)
{
  LOG_SCOPE("write_out_basis_functions()", "RBEIMEvaluation");

  if (get_parametrized_function().on_mesh_sides())
    write_out_side_basis_functions(directory_name, write_binary_basis_functions);
  else if (get_parametrized_function().on_mesh_nodes())
    write_out_node_basis_functions(directory_name, write_binary_basis_functions);
  else
    write_out_interior_basis_functions(directory_name, write_binary_basis_functions);
}

void RBEIMEvaluation::
write_out_interior_basis_functions(const std::string & directory_name,
                                   bool write_binary_basis_functions)
{
  LOG_SCOPE("write_out_interior_basis_functions()", "RBEIMEvaluation");

  // Quick return if there is no work to do. Note: make sure all procs
  // agree there is no work to do.
  bool is_empty = _local_eim_basis_functions.empty();
  this->comm().verify(is_empty);

  if (is_empty)
    return;

  // Gather basis function data from other procs, storing it in
  // _local_eim_basis_functions, so that we can then print everything
  // from processor 0.
  this->gather_bfs();

  // Write values from processor 0 only.
  if (this->processor_id() == 0)
    {
      std::vector<unsigned int> n_qp_per_elem;
      auto interior_basis_function_sizes =
        get_interior_basis_function_sizes(n_qp_per_elem);

      // Make a directory to store all the data files
      Utility::mkdir(directory_name.c_str());

      // Create filename
      std::ostringstream file_name;
      const std::string basis_function_suffix = (write_binary_basis_functions ? ".xdr" : ".dat");
      file_name << directory_name << "/" << "bf_data" << basis_function_suffix;

      // Create XDR writer object
      Xdr xdr(file_name.str(), write_binary_basis_functions ? ENCODE : WRITE);

      // Write number of basis functions to file. Note: the
      // Xdr::data() function takes non-const references, so you can't
      // pass e.g. vec.size() to that interface.
      auto n_bf = libmesh_map_find(interior_basis_function_sizes,"n_bf");
      xdr.data(n_bf, "# Number of basis functions");

      // We assume that each basis function has data for the same
      // number of elements as basis function 0, which is equal to the
      // size of the map.
      auto n_elem = libmesh_map_find(interior_basis_function_sizes,"n_elem");
      xdr.data(n_elem, "# Number of elements");

      // We assume that each element has the same number of variables,
      // and we get the number of vars from the first element of the
      // first basis function.
      auto n_vars = libmesh_map_find(interior_basis_function_sizes,"n_vars");
      xdr.data(n_vars, "# Number of variables");

      // We assume that the list of elements for each basis function
      // is the same as basis function 0. We also assume that all vars
      // have the same number of qps.
      xdr.data(n_qp_per_elem, "# Number of QPs per Elem");

      // The total amount of qp data for each var is the sum of the
      // entries in the "n_qp_per_elem" array.
      auto n_qp_data = libmesh_map_find(interior_basis_function_sizes,"n_qp_data");

      // Now we construct a vector for each basis function, for each
      // variable which is ordered according to:
      // [ [qp vals for Elem 0], [qp vals for Elem 1], ... [qp vals for Elem N] ]
      // and write it to file.
      for (auto bf_index : index_range(_local_eim_basis_functions))
        {
          auto qp_data = get_interior_basis_function_as_vec_helper(n_vars, n_qp_data, bf_index);

          // Write all the var values for this bf
          for (auto var : index_range(qp_data))
            xdr.data_stream(qp_data[var].data(), qp_data[var].size(), /*line_break=*/qp_data[var].size());
        }
    }
}

std::map<std::string,std::size_t> RBEIMEvaluation::
get_interior_basis_function_sizes(std::vector<unsigned int> & n_qp_per_elem)
{
  std::map<std::string,std::size_t> interior_basis_function_sizes;

  // Write number of basis functions to file. Note: the
  // Xdr::data() function takes non-const references, so you can't
  // pass e.g. vec.size() to that interface.
  auto n_bf = _local_eim_basis_functions.size();
  interior_basis_function_sizes["n_bf"] = n_bf;

  // We assume that each basis function has data for the same
  // number of elements as basis function 0, which is equal to the
  // size of the map.
  auto n_elem = _local_eim_basis_functions[0].size();
  interior_basis_function_sizes["n_elem"] = n_elem;

  // We assume that each element has the same number of variables,
  // and we get the number of vars from the first element of the
  // first basis function.
  auto n_vars = _local_eim_basis_functions[0].begin()->second.size();
  interior_basis_function_sizes["n_vars"] = n_vars;

  // We assume that the list of elements for each basis function
  // is the same as basis function 0. We also assume that all vars
  // have the same number of qps.
  n_qp_per_elem.clear();
  n_qp_per_elem.reserve(n_elem);
  dof_id_type expected_elem_id = 0;
  for (const auto & [actual_elem_id, array] : _local_eim_basis_functions[0])
    {
      // Note: Currently we require that the Elems are numbered
      // contiguously from [0..n_elem).  This allows us to avoid
      // writing the Elem ids to the Xdr file, but if we need to
      // generalize this assumption later, we can.
      libmesh_error_msg_if(actual_elem_id != expected_elem_id++,
                            "RBEIMEvaluation currently assumes a contiguous Elem numbering starting from 0.");

      // array[n_vars][n_qp] per Elem. We get the number of QPs
      // for variable 0, assuming they are all the same.
      n_qp_per_elem.push_back(array[0].size());
    }

  // The total amount of qp data for each var is the sum of the
  // entries in the "n_qp_per_elem" array.
  auto n_qp_data =
    std::accumulate(n_qp_per_elem.begin(),
                    n_qp_per_elem.end(),
                    0u);
  interior_basis_function_sizes["n_qp_data"] = n_qp_data;

  return interior_basis_function_sizes;
}

std::vector<std::vector<Number>> RBEIMEvaluation::
get_interior_basis_function_as_vec_helper(
  unsigned int n_vars,
  unsigned int n_qp_data,
  unsigned int bf_index)
{
  LOG_SCOPE("get_interior_basis_function_as_vec_helper()", "RBEIMEvaluation");

  std::vector<std::vector<Number>> qp_data(n_vars);

  // Reserve enough capacity in qp_data in order to do the insertions below
  // without further memory allocation.
  for (auto var : index_range(qp_data))
    qp_data[var].reserve(n_qp_data);

  // Now we construct a vector for each basis function, for each
  // variable which is ordered according to:
  // [ [qp vals for Elem 0], [qp vals for Elem 1], ... [qp vals for Elem N] ]
  // and write it to file.

  libmesh_error_msg_if(bf_index >= _local_eim_basis_functions.size(), "bf_index not valid");
  for (const auto & pr : _local_eim_basis_functions[bf_index])
    {
      // array[n_vars][n_qp] per Elem
      const auto & array = pr.second;
      for (auto var : index_range(array))
        {
          // Insert all qp values for this var
          qp_data[var].insert(/*insert at*/qp_data[var].end(),
                              /*data start*/array[var].begin(),
                              /*data end*/array[var].end());
        }
    }

  return qp_data;
}

std::vector<std::vector<std::vector<Number>>> RBEIMEvaluation::
get_interior_basis_functions_as_vecs()
{
  LOG_SCOPE("get_interior_basis_function_as_vec()", "RBEIMEvaluation");

  std::vector<std::vector<std::vector<Number>>> interior_basis_functions;

  // Quick return if there is no work to do. Note: make sure all procs
  // agree there is no work to do.
  bool is_empty = _local_eim_basis_functions.empty();
  this->comm().verify(is_empty);

  if (is_empty)
    return interior_basis_functions;

  // Gather basis function data from other procs, storing it in
  // _local_eim_basis_functions, so that we can then print everything
  // from processor 0.
  this->gather_bfs();

  if (this->processor_id() == 0)
    {
      std::vector<unsigned int> n_qp_per_elem;
      std::map<std::string,std::size_t> interior_basis_function_sizes =
        get_interior_basis_function_sizes(n_qp_per_elem);

      for (auto bf_index : index_range(_local_eim_basis_functions))
        interior_basis_functions.emplace_back(
          get_interior_basis_function_as_vec_helper(
            libmesh_map_find(interior_basis_function_sizes,"n_vars"),
            libmesh_map_find(interior_basis_function_sizes,"n_qp_data"),
            bf_index));
    }

  return interior_basis_functions;
}

void RBEIMEvaluation::
write_out_side_basis_functions(const std::string & directory_name,
                               bool write_binary_basis_functions)
{
  LOG_SCOPE("write_out_side_basis_functions()", "RBEIMEvaluation");

  // Quick return if there is no work to do. Note: make sure all procs
  // agree there is no work to do.
  bool is_empty = _local_side_eim_basis_functions.empty();
  this->comm().verify(is_empty);

  if (is_empty)
    return;

  // Gather basis function data from other procs, storing it in
  // _local_side_eim_basis_functions, so that we can then print everything
  // from processor 0.
  this->side_gather_bfs();

  // Write values from processor 0 only.
  if (this->processor_id() == 0)
    {
      // Make a directory to store all the data files
      Utility::mkdir(directory_name.c_str());

      // Create filename
      std::ostringstream file_name;
      const std::string basis_function_suffix = (write_binary_basis_functions ? ".xdr" : ".dat");
      file_name << directory_name << "/" << "bf_data" << basis_function_suffix;

      // Create XDR writer object
      Xdr xdr(file_name.str(), write_binary_basis_functions ? ENCODE : WRITE);

      // Write number of basis functions to file. Note: the
      // Xdr::data() function takes non-const references, so you can't
      // pass e.g. vec.size() to that interface.
      auto n_bf = _local_side_eim_basis_functions.size();
      xdr.data(n_bf, "# Number of basis functions");

      // We assume that each basis function has data for the same
      // number of (elem,side) pairs as basis function 0, which is equal to the
      // size of the map.
      auto n_elem = _local_side_eim_basis_functions[0].size();
      xdr.data(n_elem, "# Number of (elem,side) pairs");

      // We assume that each element has the same number of variables,
      // and we get the number of vars from the first element of the
      // first basis function.
      auto n_vars = _local_side_eim_basis_functions[0].begin()->second.size();
      xdr.data(n_vars, "# Number of variables");

      // We write out the following arrays:
      // - element IDs
      // - side indices
      // - n_qp_per_elem_side
      std::vector<unsigned int> n_qp_per_elem_side;
      std::vector<unsigned int> elem_ids;
      std::vector<unsigned int> side_indices;
      elem_ids.reserve(n_elem);
      side_indices.reserve(n_elem);
      n_qp_per_elem_side.reserve(n_elem);
      for (const auto & [elem_side_pair, array] : _local_side_eim_basis_functions[0])
        {
          elem_ids.push_back(elem_side_pair.first);
          side_indices.push_back(elem_side_pair.second);

          // array[n_vars][n_qp] per Elem. We get the number of QPs
          // for variable 0, assuming they are all the same.
          n_qp_per_elem_side.push_back(array[0].size());
        }
      xdr.data(elem_ids, "# Elem IDs");
      xdr.data(side_indices, "# Side indices");
      xdr.data(n_qp_per_elem_side, "# Number of QPs per Elem");

      // The total amount of qp data for each var is the sum of the
      // entries in the "n_qp_per_elem" array.
      auto n_qp_data =
        std::accumulate(n_qp_per_elem_side.begin(),
                        n_qp_per_elem_side.end(),
                        0u);

      // Reserve space to store contiguous vectors of qp data for each var
      std::vector<std::vector<Number>> qp_data(n_vars);
      for (auto var : index_range(qp_data))
        qp_data[var].reserve(n_qp_data);

      // Now we construct a vector for each basis function, for each
      // variable which is ordered according to:
      // [ [qp vals for Elem 0], [qp vals for Elem 1], ... [qp vals for Elem N] ]
      // and write it to file.
      for (auto bf : index_range(_local_side_eim_basis_functions))
        {
          // Clear any data from previous bf
          for (auto var : index_range(qp_data))
            qp_data[var].clear();

          for (const auto & pr : _local_side_eim_basis_functions[bf])
            {
              // array[n_vars][n_qp] per Elem
              const auto & array = pr.second;
              for (auto var : index_range(array))
                {
                  // Insert all qp values for this var
                  qp_data[var].insert(/*insert at*/qp_data[var].end(),
                                      /*data start*/array[var].begin(),
                                      /*data end*/array[var].end());
                }
            }

          // Write all the var values for this bf
          for (auto var : index_range(qp_data))
            xdr.data_stream(qp_data[var].data(), qp_data[var].size(), /*line_break=*/qp_data[var].size());
        }
    }
}

void RBEIMEvaluation::
write_out_node_basis_functions(const std::string & directory_name,
                               bool write_binary_basis_functions)
{
  LOG_SCOPE("write_out_node_basis_functions()", "RBEIMEvaluation");

  // Quick return if there is no work to do. Note: make sure all procs
  // agree there is no work to do.
  bool is_empty = _local_node_eim_basis_functions.empty();
  this->comm().verify(is_empty);

  if (is_empty)
    return;

  // Gather basis function data from other procs, storing it in
  // _local_node_eim_basis_functions, so that we can then print everything
  // from processor 0.
  this->node_gather_bfs();

  // Write values from processor 0 only.
  if (this->processor_id() == 0)
    {
      // Make a directory to store all the data files
      Utility::mkdir(directory_name.c_str());

      // Create filename
      std::ostringstream file_name;
      const std::string basis_function_suffix = (write_binary_basis_functions ? ".xdr" : ".dat");
      file_name << directory_name << "/" << "bf_data" << basis_function_suffix;

      // Create XDR writer object
      Xdr xdr(file_name.str(), write_binary_basis_functions ? ENCODE : WRITE);

      // Write number of basis functions to file. Note: the
      // Xdr::data() function takes non-const references, so you can't
      // pass e.g. vec.size() to that interface.
      auto n_bf = _local_node_eim_basis_functions.size();
      xdr.data(n_bf, "# Number of basis functions");

      // We assume that each basis function has data for the same
      // number of elements as basis function 0, which is equal to the
      // size of the map.
      auto n_node = _local_node_eim_basis_functions[0].size();
      xdr.data(n_node, "# Number of nodes");

      // We assume that each element has the same number of variables,
      // and we get the number of vars from the first element of the
      // first basis function.
      auto n_vars = _local_node_eim_basis_functions[0].begin()->second.size();
      xdr.data(n_vars, "# Number of variables");

      // We write out the following arrays:
      // - node IDs
      std::vector<unsigned int> node_ids;
      node_ids.reserve(n_node);
      for (const auto & pr : _local_node_eim_basis_functions[0])
        {
          node_ids.push_back(pr.first);
        }
      xdr.data(node_ids, "# Node IDs");

      // Now we construct a vector for each basis function, for each
      // variable which is ordered according to:
      // [ [val for Node 0], [val for Node 1], ... [val for Node N] ]
      // and write it to file.

      std::vector<std::vector<Number>> var_data(n_vars);
      for (unsigned int var=0; var<n_vars; var++)
        var_data[var].resize(n_node);

      for (auto bf : index_range(_local_node_eim_basis_functions))
        {
          unsigned int node_counter = 0;
          for (const auto & pr : _local_node_eim_basis_functions[bf])
            {
              // array[n_vars] per Node
              const auto & array = pr.second;
              for (auto var : index_range(array))
                {
                  // Based on the error check above, we know that node_id is numbered
                  // contiguously from [0..nodes], so we can use it as the vector
                  // index here.
                  var_data[var][node_counter] = array[var];
                }

              node_counter++;
            }

          // Write all the var values for this bf
          for (auto var : index_range(var_data))
            xdr.data_stream(var_data[var].data(), var_data[var].size(), /*line_break=*/var_data[var].size());
        }
    }
}

void RBEIMEvaluation::
read_in_basis_functions(const System & sys,
                        const std::string & directory_name,
                        bool read_binary_basis_functions)
{
  LOG_SCOPE("read_in_basis_functions()", "RBEIMEvaluation");

  // Return early without reading in anything if there are no basis functions
  if (get_n_basis_functions() == 0)
    return;

  if (get_parametrized_function().on_mesh_sides())
    read_in_side_basis_functions(sys, directory_name, read_binary_basis_functions);
  else if (get_parametrized_function().on_mesh_nodes())
    read_in_node_basis_functions(sys, directory_name, read_binary_basis_functions);
  else
    read_in_interior_basis_functions(sys, directory_name, read_binary_basis_functions);
}

void RBEIMEvaluation::
read_in_interior_basis_functions(const System & sys,
                                 const std::string & directory_name,
                                 bool read_binary_basis_functions)
{
  LOG_SCOPE("read_in_interior_basis_functions()", "RBEIMEvaluation");

  // Read values on processor 0 only.
  if (sys.comm().rank() == 0)
    {
      // Create filename
      std::ostringstream file_name;
      const std::string basis_function_suffix = (read_binary_basis_functions ? ".xdr" : ".dat");
      file_name << directory_name << "/" << "bf_data" << basis_function_suffix;

      // Create XDR reader object
      Xdr xdr(file_name.str(), read_binary_basis_functions ? DECODE : READ);

      // Read in the number of basis functions. The comment parameter
      // is ignored when reading.
      std::size_t n_bf;
      xdr.data(n_bf);

      // Read in the number of elements
      std::size_t n_elem;
      xdr.data(n_elem);

      // Read in the number of variables.
      std::size_t n_vars;
      xdr.data(n_vars);

      // Read in vector containing the number of QPs per elem. We can
      // create this vector with the required size or let it be read
      // from the file and sized for us.
      std::vector<unsigned int> n_qp_per_elem(n_elem);
      xdr.data(n_qp_per_elem);

      // The total amount of qp data for each var is the sum of the
      // entries in the "n_qp_per_elem" array.
      auto n_qp_data =
        std::accumulate(n_qp_per_elem.begin(),
                        n_qp_per_elem.end(),
                        0u);

      // Allocate space to store all required basis functions,
      // clearing any data that may have been there previously.
      //
      // TODO: Do we need to also write out/read in Elem ids?
      // Or can we assume they will always be contiguously
      // numbered (at least on proc 0)?
      _local_eim_basis_functions.clear();
      _local_eim_basis_functions.resize(n_bf);
      for (auto i : index_range(_local_eim_basis_functions))
        for (std::size_t elem_id=0; elem_id<n_elem; ++elem_id)
          {
            auto & array = _local_eim_basis_functions[i][elem_id];
            array.resize(n_vars);
          }

      // Allocate temporary storage for one var's worth of qp data.
      std::vector<Number> qp_data;

      // Read in data for each basis function
      for (auto i : index_range(_local_eim_basis_functions))
        {
          // Reference to the data map for the current basis function.
          auto & bf_map = _local_eim_basis_functions[i];

          for (std::size_t var=0; var<n_vars; ++var)
            {
              qp_data.clear();
              qp_data.resize(n_qp_data);

              // Read data using data_stream() since that is
              // (currently) how we write it out. The "line_break"
              // parameter of data_stream() is ignored while reading.
              xdr.data_stream(qp_data.data(), qp_data.size());

              // Iterate over the qp_data vector, filling in the
              // "small" vectors for each Elem.
              auto cursor = qp_data.begin();
              for (std::size_t elem_id=0; elem_id<n_elem; ++elem_id)
                {
                  // Get reference to the [n_vars][n_qp] array for
                  // this Elem. We assign() into the vector of
                  // quadrature point values, which allocates space if
                  // it doesn't already exist.
                  auto & array = bf_map[elem_id];
                  array[var].assign(cursor, cursor + n_qp_per_elem[elem_id]);
                  std::advance(cursor, n_qp_per_elem[elem_id]);
                }
            } // end for (var)
        } // end for (i)
    } // end if processor 0

  // Distribute the basis function information to the processors that require it
  this->distribute_bfs(sys);
}

void RBEIMEvaluation::
read_in_side_basis_functions(const System & sys,
                             const std::string & directory_name,
                             bool read_binary_basis_functions)
{
  LOG_SCOPE("read_in_basis_functions()", "RBEIMEvaluation");

  // Read values on processor 0 only.
  if (sys.comm().rank() == 0)
    {
      // Create filename
      std::ostringstream file_name;
      const std::string basis_function_suffix = (read_binary_basis_functions ? ".xdr" : ".dat");
      file_name << directory_name << "/" << "bf_data" << basis_function_suffix;

      // Create XDR reader object
      Xdr xdr(file_name.str(), read_binary_basis_functions ? DECODE : READ);

      // Read in the number of basis functions. The comment parameter
      // is ignored when reading.
      std::size_t n_bf;
      xdr.data(n_bf);

      // Read in the number of elements
      std::size_t n_elem_side;
      xdr.data(n_elem_side);

      // Read in the number of variables.
      std::size_t n_vars;
      xdr.data(n_vars);

      std::vector<unsigned int> elem_ids(n_elem_side);
      xdr.data(elem_ids);
      std::vector<unsigned int> side_indices(n_elem_side);
      xdr.data(side_indices);

      // Read in vector containing the number of QPs per elem. We can
      // create this vector with the required size or let it be read
      // from the file and sized for us.
      std::vector<unsigned int> n_qp_per_elem_side(n_elem_side);
      xdr.data(n_qp_per_elem_side);

      // The total amount of qp data for each var is the sum of the
      // entries in the "n_qp_per_elem" array.
      auto n_qp_data =
        std::accumulate(n_qp_per_elem_side.begin(),
                        n_qp_per_elem_side.end(),
                        0u);

      // Allocate space to store all required basis functions,
      // clearing any data that may have been there previously.
      _local_side_eim_basis_functions.clear();
      _local_side_eim_basis_functions.resize(n_bf);
      for (auto i : index_range(_local_side_eim_basis_functions))
        for (std::size_t elem_side_idx=0; elem_side_idx<n_elem_side; ++elem_side_idx)
          {
            unsigned int elem_id = elem_ids[elem_side_idx];
            unsigned int side_index = side_indices[elem_side_idx];
            auto elem_side_pair = std::make_pair(elem_id, side_index);

            auto & array = _local_side_eim_basis_functions[i][elem_side_pair];
            array.resize(n_vars);
          }

      // Allocate temporary storage for one var's worth of qp data.
      std::vector<Number> qp_data;

      // Read in data for each basis function
      for (auto i : index_range(_local_side_eim_basis_functions))
        {
          // Reference to the data map for the current basis function.
          auto & bf_map = _local_side_eim_basis_functions[i];

          for (std::size_t var=0; var<n_vars; ++var)
            {
              qp_data.clear();
              qp_data.resize(n_qp_data);

              // Read data using data_stream() since that is
              // (currently) how we write it out. The "line_break"
              // parameter of data_stream() is ignored while reading.
              xdr.data_stream(qp_data.data(), qp_data.size());

              // Iterate over the qp_data vector, filling in the
              // "small" vectors for each Elem.
              auto cursor = qp_data.begin();
              for (std::size_t elem_side_idx=0; elem_side_idx<n_elem_side; ++elem_side_idx)
                {
                  unsigned int elem_id = elem_ids[elem_side_idx];
                  unsigned int side_index = side_indices[elem_side_idx];
                  auto elem_side_pair = std::make_pair(elem_id, side_index);

                  // Get reference to the [n_vars][n_qp] array for
                  // this Elem. We assign() into the vector of
                  // quadrature point values, which allocates space if
                  // it doesn't already exist.
                  auto & array = bf_map[elem_side_pair];
                  array[var].assign(cursor, cursor + n_qp_per_elem_side[elem_side_idx]);
                  std::advance(cursor, n_qp_per_elem_side[elem_side_idx]);
                }
            } // end for (var)
        } // end for (i)
    } // end if processor 0

  // Distribute the basis function information to the processors that require it
  this->side_distribute_bfs(sys);
}

void RBEIMEvaluation::
read_in_node_basis_functions(const System & sys,
                             const std::string & directory_name,
                             bool read_binary_basis_functions)
{
  LOG_SCOPE("read_in_node_basis_functions()", "RBEIMEvaluation");

  // Read values on processor 0 only.
  if (sys.comm().rank() == 0)
    {
      // Create filename
      std::ostringstream file_name;
      const std::string basis_function_suffix = (read_binary_basis_functions ? ".xdr" : ".dat");
      file_name << directory_name << "/" << "bf_data" << basis_function_suffix;

      // Create XDR reader object
      Xdr xdr(file_name.str(), read_binary_basis_functions ? DECODE : READ);

      // Read in the number of basis functions. The comment parameter
      // is ignored when reading.
      std::size_t n_bf;
      xdr.data(n_bf);

      // Read in the number of nodes
      std::size_t n_node;
      xdr.data(n_node);

      // Read in the number of variables.
      std::size_t n_vars;
      xdr.data(n_vars);

      std::vector<unsigned int> node_ids(n_node);
      xdr.data(node_ids);

      // Allocate space to store all required basis functions,
      // clearing any data that may have been there previously.
      //
      // TODO: Do we need to also write out/read in Node ids?
      // Or can we assume they will always be contiguously
      // numbered (at least on proc 0)?
      _local_node_eim_basis_functions.clear();
      _local_node_eim_basis_functions.resize(n_bf);
      for (auto i : index_range(_local_node_eim_basis_functions))
        for (auto node_id : node_ids)
          {
            auto & array = _local_node_eim_basis_functions[i][node_id];
            array.resize(n_vars);
          }

      // Read data into node_value from xdr
      std::vector<Number> node_value;

      // Read in data for each basis function
      for (auto i : index_range(_local_node_eim_basis_functions))
        {
          // Reference to the data map for the current basis function.
          auto & bf_map = _local_node_eim_basis_functions[i];

          for (std::size_t var=0; var<n_vars; ++var)
            {
              node_value.clear();
              node_value.resize(n_node);

              // Read data using data_stream() since that is
              // (currently) how we write it out. The "line_break"
              // parameter of data_stream() is ignored while reading.
              xdr.data_stream(node_value.data(), node_value.size());

              for (unsigned int node_counter=0; node_counter<n_node; node_counter++)
                {
                  auto & array = bf_map[node_ids[node_counter]];
                  array[var] = node_value[node_counter];
                }
            } // end for (var)
        } // end for (i)
    } // end if processor 0

  // Distribute the basis function information to the processors that require it
  this->node_distribute_bfs(sys);
}

void RBEIMEvaluation::print_local_eim_basis_functions() const
{
  for (auto bf : index_range(_local_eim_basis_functions))
    {
      libMesh::out << "Interior basis function " << bf << std::endl;
      for (const auto & [elem_id, array] : _local_eim_basis_functions[bf])
        {
          libMesh::out << "Elem " << elem_id << std::endl;
          for (auto var : index_range(array))
            {
              libMesh::out << "Variable " << var << std::endl;
              for (auto qp : index_range(array[var]))
                libMesh::out << array[var][qp] << " ";
              libMesh::out << std::endl;
            }
        }
    }

  for (auto bf : index_range(_local_side_eim_basis_functions))
    {
      libMesh::out << "Side basis function " << bf << std::endl;
      for (const auto & [pr, array] : _local_side_eim_basis_functions[bf])
        {
          const auto & elem_id = pr.first;
          const auto & side_index = pr.second;
          libMesh::out << "Elem " << elem_id << ", Side " << side_index << std::endl;
          for (auto var : index_range(array))
            {
              libMesh::out << "Variable " << var << std::endl;
              for (auto qp : index_range(array[var]))
                libMesh::out << array[var][qp] << " ";
              libMesh::out << std::endl;
            }
        }
    }

  for (auto bf : index_range(_local_node_eim_basis_functions))
    {
      libMesh::out << "Node basis function " << bf << std::endl;
      for (const auto & [node_id, array] : _local_node_eim_basis_functions[bf])
        {
          libMesh::out << "Node " << node_id << std::endl;
          for (auto var : index_range(array))
            {
              libMesh::out << "Variable " << var << ": " << array[var] << std::endl;
            }
          libMesh::out << std::endl;
        }
    }
}

void RBEIMEvaluation::gather_bfs()
{
  // We need to gather _local_eim_basis_functions data from other
  // procs for printing.
  //
  // Ideally, this could be accomplished by simply calling:
  // this->comm().gather(/*root_id=*/0, _local_eim_basis_functions);
  //
  // but the data structure seems to be too complicated for this to
  // work automatically. (I get some error about the function called
  // being "private within this context".) Therefore, we have to
  // gather the information manually.

  // So we can avoid calling this many times below
  auto n_procs = this->n_processors();

  // In serial there's nothing to gather
  if (n_procs == 1)
    return;

  // Current assumption is that the number of basis functions stored on
  // each processor is the same, the only thing that differs is the number
  // of elements, so make sure that is the case now.
  auto n_bf = _local_eim_basis_functions.size();
  this->comm().verify(n_bf);

  // This function should never be called if there are no basis
  // functions, so if it was, something went wrong.
  libmesh_error_msg_if(!n_bf, "RBEIMEvaluation::gather_bfs() should not be called with 0 basis functions.");

  // The number of variables should be the same on all processors
  // and we can get this from _local_eim_basis_functions. However,
  // it may be that some processors have no local elements, so on
  // those processors we cannot look up the size from
  // _local_eim_basis_functions. As a result we use comm().max(n_vars)
  // to make sure all processors agree on the final value.
  std::size_t n_vars =
    _local_eim_basis_functions[0].empty() ? 0 : _local_eim_basis_functions[0].begin()->second.size();
  this->comm().max(n_vars);

  // Gather list of Elem ids stored on each processor to proc 0.  We
  // use basis function 0 as an example and assume all the basis
  // functions are distributed similarly.
  std::vector<dof_id_type> elem_ids;
  elem_ids.reserve(_local_eim_basis_functions[0].size());
  for (const auto & pr : _local_eim_basis_functions[0])
    elem_ids.push_back(pr.first);
  this->comm().gather(/*root_id=*/0, elem_ids);

  // Store the number of qps per Elem on this processor. Again, use
  // basis function 0 (and variable 0) to get this information, then
  // apply it to all basis functions.
  std::vector<unsigned int> n_qp_per_elem;
  n_qp_per_elem.reserve(_local_eim_basis_functions[0].size());
  for (const auto & pr : _local_eim_basis_functions[0])
    {
      // array[n_vars][n_qp] per Elem. We get the number of QPs
      // for variable 0, assuming they are all the same.
      const auto & array = pr.second;
      n_qp_per_elem.push_back(array[0].size());
    }

  // Before gathering, compute the total amount of local qp data for
  // each var, which is the sum of the entries in the "n_qp_per_elem" array.
  // This will be used to reserve space in a vector below.
  auto n_local_qp_data =
    std::accumulate(n_qp_per_elem.begin(),
                    n_qp_per_elem.end(),
                    0u);

  // Gather the number of qps per Elem for each processor onto processor 0.
  this->comm().gather(/*root_id=*/0, n_qp_per_elem);

  // Sanity check: On processor 0, this checks that we have gathered the same number
  // of elem ids and qp counts.
  libmesh_error_msg_if(elem_ids.size() != n_qp_per_elem.size(),
                       "Must gather same number of Elem ids as qps per Elem.");

  // Reserve space to store contiguous vectors of qp data for each var
  std::vector<std::vector<Number>> gathered_qp_data(n_vars);
  for (auto var : index_range(gathered_qp_data))
    gathered_qp_data[var].reserve(n_local_qp_data);

  // Now we construct a vector for each basis function, for each
  // variable, which is ordered according to:
  // [ [qp vals for Elem 0], [qp vals for Elem 1], ... [qp vals for Elem N] ]
  // and gather it to processor 0.
  for (auto bf : index_range(_local_eim_basis_functions))
    {
      // Clear any data from previous bf
      for (auto var : index_range(gathered_qp_data))
        gathered_qp_data[var].clear();

      for (const auto & pr : _local_eim_basis_functions[bf])
        {
          // array[n_vars][n_qp] per Elem
          const auto & array = pr.second;
          for (auto var : index_range(array))
            {
              // Insert all qp values for this var
              gathered_qp_data[var].insert(/*insert at*/gathered_qp_data[var].end(),
                                           /*data start*/array[var].begin(),
                                           /*data end*/array[var].end());
            }
        }

      // Reference to the data map for the current basis function.
      auto & bf_map = _local_eim_basis_functions[bf];

      for (auto var : index_range(gathered_qp_data))
        {
          // For each var, gather gathered_qp_data[var] onto processor
          // 0. There apparently is not a gather overload for
          // vector-of-vectors...
          this->comm().gather(/*root_id=*/0, gathered_qp_data[var]);

          // On processor 0, iterate over the gathered_qp_data[var]
          // vector we just gathered, filling in the "small" vectors
          // for each Elem. Note: here we ignore the fact that we
          // already have the data on processor 0 and just overwrite
          // it, this makes the indexing logic a bit simpler.
          if (this->processor_id() == 0)
            {
              auto cursor = gathered_qp_data[var].begin();
              for (auto i : index_range(elem_ids))
                {
                  auto elem_id = elem_ids[i];
                  auto n_qp_this_elem = n_qp_per_elem[i];

                  // Get reference to the [n_vars][n_qp] array for
                  // this Elem. We assign() into the vector of
                  // quadrature point values, which allocates space if
                  // it doesn't already exist.
                  auto & array = bf_map[elem_id];

                  // Possibly allocate space if this is data for a new
                  // element we haven't seen before.
                  if (array.empty())
                    array.resize(n_vars);

                  array[var].assign(cursor, cursor + n_qp_this_elem);
                  std::advance(cursor, n_qp_this_elem);
                }
            }
        }
    } // end loop over basis functions
}

void RBEIMEvaluation::side_gather_bfs()
{
  // We need to gather _local_side_eim_basis_functions data from other
  // procs for printing.
  //
  // Ideally, this could be accomplished by simply calling:
  // this->comm().gather(/*root_id=*/0, _local_side_eim_basis_functions);
  //
  // but the data structure seems to be too complicated for this to
  // work automatically. (I get some error about the function called
  // being "private within this context".) Therefore, we have to
  // gather the information manually.

  // So we can avoid calling this many times below
  auto n_procs = this->n_processors();

  // In serial there's nothing to gather
  if (n_procs == 1)
    return;

  // Current assumption is that the number of basis functions stored on
  // each processor is the same, the only thing that differs is the number
  // of elements, so make sure that is the case now.
  auto n_bf = _local_side_eim_basis_functions.size();
  this->comm().verify(n_bf);

  // This function should never be called if there are no basis
  // functions, so if it was, something went wrong.
  libmesh_error_msg_if(!n_bf, "SideRBEIMEvaluation::gather_bfs() should not be called with 0 basis functions.");

  // The number of variables should be the same on all processors
  // and we can get this from _local_side_eim_basis_functions. However,
  // it may be that some processors have no local elements, so on
  // those processors we cannot look up the size from
  // _local_side_eim_basis_functions. As a result we use comm().max(n_vars)
  // to make sure all processors agree on the final value.
  std::size_t n_vars =
    _local_side_eim_basis_functions[0].empty() ? 0 : _local_side_eim_basis_functions[0].begin()->second.size();
  this->comm().max(n_vars);

  // Gather list of (elem,side) pairs stored on each processor to proc 0.  We
  // use basis function 0 as an example and assume all the basis
  // functions are distributed similarly.
  std::vector<std::pair<dof_id_type,unsigned int>> elem_side_pairs;
  elem_side_pairs.reserve(_local_side_eim_basis_functions[0].size());
  for (const auto & pr : _local_side_eim_basis_functions[0])
    elem_side_pairs.push_back(pr.first);
  this->comm().gather(/*root_id=*/0, elem_side_pairs);

  // Store the number of qps per Elem on this processor. Again, use
  // basis function 0 (and variable 0) to get this information, then
  // apply it to all basis functions.
  std::vector<unsigned int> n_qp_per_elem_side;
  n_qp_per_elem_side.reserve(_local_side_eim_basis_functions[0].size());
  for (const auto & pr : _local_side_eim_basis_functions[0])
    {
      // array[n_vars][n_qp] per (elem,side). We get the number of QPs
      // for variable 0, assuming they are all the same.
      const auto & array = pr.second;
      n_qp_per_elem_side.push_back(array[0].size());
    }

  // Before gathering, compute the total amount of local qp data for
  // each var, which is the sum of the entries in the "n_qp_per_elem_side" array.
  // This will be used to reserve space in a vector below.
  auto n_local_qp_data =
    std::accumulate(n_qp_per_elem_side.begin(),
                    n_qp_per_elem_side.end(),
                    0u);

  // Gather the number of qps per Elem for each processor onto processor 0.
  this->comm().gather(/*root_id=*/0, n_qp_per_elem_side);

  // Sanity check: On processor 0, this checks that we have gathered the same number
  // of (elem,side) pairs and qp counts.
  libmesh_error_msg_if(elem_side_pairs.size() != n_qp_per_elem_side.size(),
                       "Must gather same number of Elem ids as qps per Elem.");

  // Reserve space to store contiguous vectors of qp data for each var
  std::vector<std::vector<Number>> gathered_qp_data(n_vars);
  for (auto var : index_range(gathered_qp_data))
    gathered_qp_data[var].reserve(n_local_qp_data);

  // Now we construct a vector for each basis function, for each
  // variable, which is ordered according to:
  // [ [qp vals for Elem 0], [qp vals for Elem 1], ... [qp vals for Elem N] ]
  // and gather it to processor 0.
  for (auto bf : index_range(_local_side_eim_basis_functions))
    {
      // Clear any data from previous bf
      for (auto var : index_range(gathered_qp_data))
        gathered_qp_data[var].clear();

      for (const auto & pr : _local_side_eim_basis_functions[bf])
        {
          // array[n_vars][n_qp] per (elem,side) pair
          const auto & array = pr.second;
          for (auto var : index_range(array))
            {
              // Insert all qp values for this var
              gathered_qp_data[var].insert(/*insert at*/gathered_qp_data[var].end(),
                                           /*data start*/array[var].begin(),
                                           /*data end*/array[var].end());
            }
        }

      // Reference to the data map for the current basis function.
      auto & bf_map = _local_side_eim_basis_functions[bf];

      for (auto var : index_range(gathered_qp_data))
        {
          // For each var, gather gathered_qp_data[var] onto processor
          // 0. There apparently is not a gather overload for
          // vector-of-vectors...
          this->comm().gather(/*root_id=*/0, gathered_qp_data[var]);

          // On processor 0, iterate over the gathered_qp_data[var]
          // vector we just gathered, filling in the "small" vectors
          // for each Elem. Note: here we ignore the fact that we
          // already have the data on processor 0 and just overwrite
          // it, this makes the indexing logic a bit simpler.
          if (this->processor_id() == 0)
            {
              auto cursor = gathered_qp_data[var].begin();
              for (auto i : index_range(elem_side_pairs))
                {
                  auto elem_side_pair = elem_side_pairs[i];
                  auto n_qp_this_elem_side = n_qp_per_elem_side[i];

                  // Get reference to the [n_vars][n_qp] array for
                  // this Elem. We assign() into the vector of
                  // quadrature point values, which allocates space if
                  // it doesn't already exist.
                  auto & array = bf_map[elem_side_pair];

                  // Possibly allocate space if this is data for a new
                  // element we haven't seen before.
                  if (array.empty())
                    array.resize(n_vars);

                  array[var].assign(cursor, cursor + n_qp_this_elem_side);
                  std::advance(cursor, n_qp_this_elem_side);
                }
            }
        }
    } // end loop over basis functions
}

void RBEIMEvaluation::node_gather_bfs()
{
  // We need to gather _local_node_eim_basis_functions data from other
  // procs for printing.
  //
  // Ideally, this could be accomplished by simply calling:
  // this->comm().gather(/*root_id=*/0, _local_node_eim_basis_functions);
  //
  // but the data structure seems to be too complicated for this to
  // work automatically. (I get some error about the function called
  // being "private within this context".) Therefore, we have to
  // gather the information manually.

  // So we can avoid calling this many times below
  auto n_procs = this->n_processors();

  // In serial there's nothing to gather
  if (n_procs == 1)
    return;

  // Current assumption is that the number of basis functions stored on
  // each processor is the same, the only thing that differs is the number
  // of elements, so make sure that is the case now.
  auto n_bf = _local_node_eim_basis_functions.size();
  this->comm().verify(n_bf);

  // This function should never be called if there are no basis
  // functions, so if it was, something went wrong.
  libmesh_error_msg_if(!n_bf, "RBEIMEvaluation::gather_bfs() should not be called with 0 basis functions.");

  // The number of variables should be the same on all processors
  // and we can get this from _local_eim_basis_functions. However,
  // it may be that some processors have no local elements, so on
  // those processors we cannot look up the size from
  // _local_eim_basis_functions. As a result we use comm().max(n_vars)
  // to make sure all processors agree on the final value.
  std::size_t n_vars =
    _local_node_eim_basis_functions[0].empty() ? 0 : _local_node_eim_basis_functions[0].begin()->second.size();
  this->comm().max(n_vars);

  // Gather list of Node ids stored on each processor to proc 0.  We
  // use basis function 0 as an example and assume all the basis
  // functions are distributed similarly.
  unsigned int n_local_nodes = _local_node_eim_basis_functions[0].size();
  std::vector<dof_id_type> node_ids;
  node_ids.reserve(n_local_nodes);
  for (const auto & pr : _local_node_eim_basis_functions[0])
    node_ids.push_back(pr.first);
  this->comm().gather(/*root_id=*/0, node_ids);

  // Now we construct a vector for each basis function, for each
  // variable, which is ordered according to:
  // [ [val for Node 0], [val for Node 1], ... [val for Node N] ]
  // and gather it to processor 0.
  std::vector<std::vector<Number>> gathered_node_data(n_vars);
  for (auto bf : index_range(_local_node_eim_basis_functions))
    {
      // Clear any data from previous bf
      for (auto var : index_range(gathered_node_data))
        {
          gathered_node_data[var].clear();
          gathered_node_data[var].resize(n_local_nodes);
        }

      unsigned int local_node_idx = 0;
      for (const auto & pr : _local_node_eim_basis_functions[bf])
        {
          // array[n_vars] per Node
          const auto & array = pr.second;
          for (auto var : index_range(array))
            {
              gathered_node_data[var][local_node_idx] = array[var];
            }

          local_node_idx++;
        }

      // Reference to the data map for the current basis function.
      auto & bf_map = _local_node_eim_basis_functions[bf];

      for (auto var : index_range(gathered_node_data))
        {
          // For each var, gather gathered_qp_data[var] onto processor
          // 0. There apparently is not a gather overload for
          // vector-of-vectors...
          this->comm().gather(/*root_id=*/0, gathered_node_data[var]);

          // On processor 0, iterate over the gathered_qp_data[var]
          // vector we just gathered, filling in the "small" vectors
          // for each Elem. Note: here we ignore the fact that we
          // already have the data on processor 0 and just overwrite
          // it, this makes the indexing logic a bit simpler.
          if (this->processor_id() == 0)
            {
              auto cursor = gathered_node_data[var].begin();
              for (auto i : index_range(node_ids))
                {
                  auto node_id = node_ids[i];

                  // Get reference to the [n_vars] array for
                  // this Node. We assign() into the vector of
                  // node values, which allocates space if
                  // it doesn't already exist.
                  auto & array = bf_map[node_id];

                  // Possibly allocate space if this is data for a new
                  // node we haven't seen before.
                  if (array.empty())
                    array.resize(n_vars);

                  // There is only one value per variable per node, so
                  // we set the value by de-referencing cursor, and
                  // then advance the cursor by 1.
                  array[var] = *cursor;
                  std::advance(cursor, 1);
                }
            }
        }
    } // end loop over basis functions
}



void RBEIMEvaluation::distribute_bfs(const System & sys)
{
  // So we can avoid calling these many times below
  auto n_procs = sys.comm().size();
  auto rank = sys.comm().rank();

  // In serial there's nothing to distribute
  if (n_procs == 1)
    return;

  // Broadcast the number of basis functions from proc 0. After
  // distributing, all procs should have the same number of basis
  // functions.
  auto n_bf = _local_eim_basis_functions.size();
  sys.comm().broadcast(n_bf);

  // Allocate enough space to store n_bf basis functions on non-zero ranks
  if (rank != 0)
    _local_eim_basis_functions.resize(n_bf);

  // Broadcast the number of variables from proc 0. After
  // distributing, all procs should have the same number of variables.
  auto n_vars = _local_eim_basis_functions[0].begin()->second.size();
  sys.comm().broadcast(n_vars);

  // Construct lists of elem ids owned by different processors
  const MeshBase & mesh = sys.get_mesh();

  std::vector<dof_id_type> gathered_local_elem_ids;
  gathered_local_elem_ids.reserve(mesh.n_elem());
  for (const auto & elem : mesh.active_local_element_ptr_range())
    gathered_local_elem_ids.push_back(elem->id());

  // I _think_ the local elem ids are likely to already be sorted in
  // ascending order, since that is how they are stored on the Mesh,
  // but we can always just guarantee this to be on the safe side as
  // well.
  std::sort(gathered_local_elem_ids.begin(), gathered_local_elem_ids.end());

  // Gather the number of local elems from all procs to proc 0
  auto n_local_elems = gathered_local_elem_ids.size();
  std::vector<std::size_t> gathered_n_local_elems = {n_local_elems};
  sys.comm().gather(/*root_id=*/0, gathered_n_local_elems);

  // Gather the elem ids owned by each processor onto processor 0.
  sys.comm().gather(/*root_id=*/0, gathered_local_elem_ids);

  // Construct vectors of "start" and "one-past-the-end" indices into
  // the gathered_local_elem_ids vector for each proc. Only valid on
  // processor 0.
  std::vector<std::size_t> start_elem_ids_index, end_elem_ids_index;

  if (rank == 0)
    {
      start_elem_ids_index.resize(n_procs);
      start_elem_ids_index[0] = 0;
      for (processor_id_type p=1; p<n_procs; ++p)
        start_elem_ids_index[p] = start_elem_ids_index[p-1] + gathered_n_local_elems[p-1];

      end_elem_ids_index.resize(n_procs);
      end_elem_ids_index[n_procs - 1] = gathered_local_elem_ids.size();
      for (processor_id_type p=0; p<n_procs - 1; ++p)
        end_elem_ids_index[p] = start_elem_ids_index[p+1];
    }

  // On processor 0, using basis function 0 and variable 0, prepare a
  // vector with the number of qps per Elem.  Then scatter this vector
  // out to the processors that require it. The order of this vector
  // matches the gathered_local_elem_ids ordering. The counts will be
  // gathered_n_local_elems, since there will be one qp count per Elem.
  std::vector<unsigned int> n_qp_per_elem_data;

  // On rank 0, the "counts" vector holds the number of floating point values that
  // are to be scattered to each proc. It is only required on proc 0.
  std::vector<int> counts;

  if (rank == 0)
    {
      n_qp_per_elem_data.reserve(gathered_local_elem_ids.size());
      counts.resize(n_procs);

      auto & bf_map = _local_eim_basis_functions[0];

      for (processor_id_type p=0; p<n_procs; ++p)
        {
          for (auto e : make_range(start_elem_ids_index[p], end_elem_ids_index[p]))
            {
              auto elem_id = gathered_local_elem_ids[e];

              // Get reference to array[n_vars][n_qp] for current Elem.
              // Throws an error if the required elem_id is not found.
              const auto & array = libmesh_map_find(bf_map, elem_id);

              auto n_qps = array[0].size();

              // We use var==0 to set the number of qps for all vars
              n_qp_per_elem_data.push_back(n_qps);

              // Accumulate the count for this proc
              counts[p] += n_qps;
            } // end for (e)
        } // end for proc_id
    } // if (rank == 0)

  // Now scatter the n_qp_per_elem_data to all procs (must call the
  // scatter on all procs, it is a collective).
  {
    std::vector<unsigned int> recv;
    std::vector<int> tmp(gathered_n_local_elems.begin(), gathered_n_local_elems.end());
    sys.comm().scatter(n_qp_per_elem_data, tmp, recv, /*root_id=*/0);

    // Now swap n_qp_per_elem_data and recv. All processors now have a
    // vector of length n_local_elems containing the number of
    // quadarature points per Elem.
    n_qp_per_elem_data.swap(recv);
  }

  // For each basis function and each variable, build a vector
  // of qp data in the Elem ordering given by the
  // gathered_local_elem_ids, then call
  //
  // sys.comm().scatter(data, counts, recv, /*root_id=*/0);
  std::vector<std::vector<Number>> qp_data(n_vars);
  if (rank == 0)
    {
      // The total amount of qp data is given by summing the entries
      // of the "counts" vector.
      auto n_qp_data =
        std::accumulate(counts.begin(), counts.end(), 0u);

      // On processor 0, reserve enough space to hold all the qp
      // data for a single basis function for each var.
      for (auto var : index_range(qp_data))
        qp_data[var].reserve(n_qp_data);
    }

  // The recv_qp_data vector will be used on the receiving end of all
  // the scatters below.
  std::vector<Number> recv_qp_data;

  // Loop from 0..n_bf on _all_ procs, since the scatters inside this
  // loop are collective.
  for (auto bf : make_range(n_bf))
    {
      // Prepare data for scattering (only on proc 0)
      if (rank == 0)
        {
          // Reference to the data map for the current basis function.
          auto & bf_map = _local_eim_basis_functions[bf];

          // Clear any data from previous bf
          for (auto var : index_range(qp_data))
            qp_data[var].clear();

          for (processor_id_type p=0; p<n_procs; ++p)
            {
              for (auto e : make_range(start_elem_ids_index[p], end_elem_ids_index[p]))
                {
                  auto elem_id = gathered_local_elem_ids[e];

                  // Get reference to array[n_vars][n_qp] for current Elem.
                  // Throws an error if the required elem_id is not found.
                  const auto & array = libmesh_map_find(bf_map, elem_id);

                  for (auto var : index_range(array))
                    {
                      // Insert all qp values for this var
                      qp_data[var].insert(/*insert at*/qp_data[var].end(),
                                          /*data start*/array[var].begin(),
                                          /*data end*/array[var].end());
                    } // end for (var)
                } // end for (e)
            } // end for proc_id
        } // end if rank==0

      // Perform the scatters (all procs)
      for (auto var : make_range(n_vars))
        {
          // Do the scatter for the current var
          sys.comm().scatter(qp_data[var], counts, recv_qp_data, /*root_id=*/0);

          if (rank != 0)
            {
              // Store the scattered data we received in _local_eim_basis_functions[bf]
              auto & bf_map = _local_eim_basis_functions[bf];
              auto cursor = recv_qp_data.begin();

              for (auto i : index_range(gathered_local_elem_ids))
                {
                  auto elem_id = gathered_local_elem_ids[i];
                  auto n_qp_this_elem = n_qp_per_elem_data[i];
                  auto & array = bf_map[elem_id];

                  // Create space to store the data if it doesn't already exist.
                  if (array.empty())
                    array.resize(n_vars);

                  array[var].assign(cursor, cursor + n_qp_this_elem);
                  std::advance(cursor, n_qp_this_elem);
                }
            } // if (rank != 0)
        } // end for (var)
    } // end for (bf)

  // Now that the scattering is done, delete non-local Elem
  // information from processor 0's _local_eim_basis_functions data
  // structure.
  if (rank == 0)
    {
      for (processor_id_type p=1; p<n_procs; ++p)
        {
          for (auto e : make_range(start_elem_ids_index[p], end_elem_ids_index[p]))
            {
              auto elem_id = gathered_local_elem_ids[e];

              // Delete this Elem's information from every basis function.
              for (auto & bf_map : _local_eim_basis_functions)
                bf_map.erase(elem_id);
            } // end for (e)
        } // end for proc_id
    } // if (rank == 0)
}

void RBEIMEvaluation::side_distribute_bfs(const System & sys)
{
  // So we can avoid calling these many times below
  auto n_procs = sys.comm().size();
  auto rank = sys.comm().rank();

  // In serial there's nothing to distribute
  if (n_procs == 1)
    return;

  // Broadcast the number of basis functions from proc 0. After
  // distributing, all procs should have the same number of basis
  // functions.
  auto n_bf = _local_side_eim_basis_functions.size();
  sys.comm().broadcast(n_bf);

  // Allocate enough space to store n_bf basis functions on non-zero ranks
  if (rank != 0)
    _local_side_eim_basis_functions.resize(n_bf);

  // Broadcast the number of variables from proc 0. After
  // distributing, all procs should have the same number of variables.
  auto n_vars = _local_side_eim_basis_functions[0].begin()->second.size();
  sys.comm().broadcast(n_vars);

  const std::set<boundary_id_type> & parametrized_function_boundary_ids =
    get_parametrized_function().get_parametrized_function_boundary_ids();

  // Construct lists of elem ids owned by different processors
  const MeshBase & mesh = sys.get_mesh();

  // BoundaryInfo and related data structures
  const auto & binfo = mesh.get_boundary_info();
  std::vector<boundary_id_type> side_boundary_ids;

  std::vector<dof_id_type> gathered_local_elem_ids;
  std::vector<dof_id_type> gathered_local_side_indices;
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      for (unsigned int side = 0; side != elem->n_sides(); ++side)
        {
          // skip non-boundary elements
          if (!elem->neighbor_ptr(side))
            {
              binfo.boundary_ids(elem, side, side_boundary_ids);

              bool has_side_boundary_id = false;
              for (boundary_id_type side_boundary_id : side_boundary_ids)
                if (parametrized_function_boundary_ids.count(side_boundary_id))
                  {
                    has_side_boundary_id = true;
                    break;
                  }

              if (has_side_boundary_id)
                {
                  gathered_local_elem_ids.push_back(elem->id());
                  gathered_local_side_indices.push_back(side);
                }
            }
        }

        // In the case of 2D elements, we also check the shellfaces
        if (elem->dim() == 2)
          for (unsigned int shellface_index=0; shellface_index<2; shellface_index++)
            {
              binfo.shellface_boundary_ids(elem, shellface_index, side_boundary_ids);

              bool has_side_boundary_id = false;
              for (boundary_id_type side_boundary_id : side_boundary_ids)
                if (parametrized_function_boundary_ids.count(side_boundary_id))
                  {
                    has_side_boundary_id = true;
                    break;
                  }

              if (has_side_boundary_id)
                {
                  // We use shellface_index as the side_index since shellface boundary conditions
                  // are stored separately from side boundary conditions in BoundaryInfo.
                  gathered_local_elem_ids.push_back(elem->id());
                  gathered_local_side_indices.push_back(shellface_index);
                }
            }
    }

  // Gather the number of local elems from all procs to proc 0
  auto n_local_elems = gathered_local_elem_ids.size();
  std::vector<std::size_t> gathered_n_local_elems = {n_local_elems};
  sys.comm().gather(/*root_id=*/0, gathered_n_local_elems);

  // Gather the (elem,side) owned by each processor onto processor 0.
  sys.comm().gather(/*root_id=*/0, gathered_local_elem_ids);
  sys.comm().gather(/*root_id=*/0, gathered_local_side_indices);

  // Construct vectors of "start" and "one-past-the-end" indices into
  // the gathered_local_elem_ids vector for each proc. Only valid on
  // processor 0.
  std::vector<std::size_t> start_elem_ids_index, end_elem_ids_index;

  if (rank == 0)
    {
      start_elem_ids_index.resize(n_procs);
      start_elem_ids_index[0] = 0;
      for (processor_id_type p=1; p<n_procs; ++p)
        start_elem_ids_index[p] = start_elem_ids_index[p-1] + gathered_n_local_elems[p-1];

      end_elem_ids_index.resize(n_procs);
      end_elem_ids_index[n_procs - 1] = gathered_local_elem_ids.size();
      for (processor_id_type p=0; p<n_procs - 1; ++p)
        end_elem_ids_index[p] = start_elem_ids_index[p+1];
    }

  // On processor 0, using basis function 0 and variable 0, prepare a
  // vector with the number of qps per Elem.  Then scatter this vector
  // out to the processors that require it. The order of this vector
  // matches the gathered_local_elem_ids ordering. The counts will be
  // gathered_n_local_elems, since there will be one qp count per Elem.
  std::vector<unsigned int> n_qp_per_elem_data;

  // On rank 0, the "counts" vector holds the number of floating point values that
  // are to be scattered to each proc. It is only required on proc 0.
  std::vector<int> counts;

  if (rank == 0)
    {
      n_qp_per_elem_data.reserve(gathered_local_elem_ids.size());
      counts.resize(n_procs);

      auto & bf_map = _local_side_eim_basis_functions[0];

      for (processor_id_type p=0; p<n_procs; ++p)
        {
          for (auto e : make_range(start_elem_ids_index[p], end_elem_ids_index[p]))
            {
              auto elem_id = gathered_local_elem_ids[e];
              auto side_index = gathered_local_side_indices[e];

              // Get reference to array[n_vars][n_qp] for current Elem.
              // Throws an error if the required elem_id is not found.
              const auto & array = libmesh_map_find(bf_map, std::make_pair(elem_id, side_index));

              auto n_qps = array[0].size();

              // We use var==0 to set the number of qps for all vars
              n_qp_per_elem_data.push_back(n_qps);

              // Accumulate the count for this proc
              counts[p] += n_qps;
            } // end for (e)
        } // end for proc_id
    } // if (rank == 0)

  // Now scatter the n_qp_per_elem_data to all procs (must call the
  // scatter on all procs, it is a collective).
  {
    std::vector<unsigned int> recv;
    std::vector<int> tmp(gathered_n_local_elems.begin(), gathered_n_local_elems.end());
    sys.comm().scatter(n_qp_per_elem_data, tmp, recv, /*root_id=*/0);

    // Now swap n_qp_per_elem_data and recv. All processors now have a
    // vector of length n_local_elems containing the number of
    // quadarature points per Elem.
    n_qp_per_elem_data.swap(recv);
  }

  // For each basis function and each variable, build a vector
  // of qp data in the Elem ordering given by the
  // gathered_local_elem_ids, then call
  //
  // sys.comm().scatter(data, counts, recv, /*root_id=*/0);
  std::vector<std::vector<Number>> qp_data(n_vars);
  if (rank == 0)
    {
      // The total amount of qp data is given by summing the entries
      // of the "counts" vector.
      auto n_qp_data =
        std::accumulate(counts.begin(), counts.end(), 0u);

      // On processor 0, reserve enough space to hold all the qp
      // data for a single basis function for each var.
      for (auto var : index_range(qp_data))
        qp_data[var].reserve(n_qp_data);
    }

  // The recv_qp_data vector will be used on the receiving end of all
  // the scatters below.
  std::vector<Number> recv_qp_data;

  // Loop from 0..n_bf on _all_ procs, since the scatters inside this
  // loop are collective.
  for (auto bf : make_range(n_bf))
    {
      // Prepare data for scattering (only on proc 0)
      if (rank == 0)
        {
          // Reference to the data map for the current basis function.
          auto & bf_map = _local_side_eim_basis_functions[bf];

          // Clear any data from previous bf
          for (auto var : index_range(qp_data))
            qp_data[var].clear();

          for (processor_id_type p=0; p<n_procs; ++p)
            {
              for (auto e : make_range(start_elem_ids_index[p], end_elem_ids_index[p]))
                {
                  auto elem_id = gathered_local_elem_ids[e];
                  auto side_index = gathered_local_side_indices[e];

                  // Get reference to array[n_vars][n_qp] for current Elem.
                  // Throws an error if the required (elem,side) is not found.
                  const auto & array = libmesh_map_find(bf_map, std::make_pair(elem_id, side_index));

                  for (auto var : index_range(array))
                    {
                      // Insert all qp values for this var
                      qp_data[var].insert(/*insert at*/qp_data[var].end(),
                                          /*data start*/array[var].begin(),
                                          /*data end*/array[var].end());
                    } // end for (var)
                } // end for (e)
            } // end for proc_id
        } // end if rank==0

      // Perform the scatters (all procs)
      for (auto var : make_range(n_vars))
        {
          // Do the scatter for the current var
          sys.comm().scatter(qp_data[var], counts, recv_qp_data, /*root_id=*/0);

          if (rank != 0)
            {
              // Store the scattered data we received in _local_side_eim_basis_functions[bf]
              auto & bf_map = _local_side_eim_basis_functions[bf];
              auto cursor = recv_qp_data.begin();

              for (auto i : index_range(gathered_local_elem_ids))
                {
                  auto elem_id = gathered_local_elem_ids[i];
                  auto side_index = gathered_local_side_indices[i];
                  auto n_qp_this_elem = n_qp_per_elem_data[i];
                  auto & array = bf_map[std::make_pair(elem_id, side_index)];

                  // Create space to store the data if it doesn't already exist.
                  if (array.empty())
                    array.resize(n_vars);

                  array[var].assign(cursor, cursor + n_qp_this_elem);
                  std::advance(cursor, n_qp_this_elem);
                }
            } // if (rank != 0)
        } // end for (var)
    } // end for (bf)

  // Now that the scattering is done, delete non-local Elem
  // information from processor 0's _local_side_eim_basis_functions data
  // structure.
  if (rank == 0)
    {
      for (processor_id_type p=1; p<n_procs; ++p)
        {
          for (auto e : make_range(start_elem_ids_index[p], end_elem_ids_index[p]))
            {
              auto elem_id = gathered_local_elem_ids[e];
              auto side_index = gathered_local_side_indices[e];

              // Delete this Elem's information from every basis function.
              for (auto & bf_map : _local_side_eim_basis_functions)
                bf_map.erase(std::make_pair(elem_id, side_index));
            } // end for (e)
        } // end for proc_id
    } // if (rank == 0)
}

void RBEIMEvaluation::node_distribute_bfs(const System & sys)
{
  // So we can avoid calling these many times below
  auto n_procs = sys.comm().size();
  auto rank = sys.comm().rank();

  // In serial there's nothing to distribute
  if (n_procs == 1)
    return;

  // Broadcast the number of basis functions from proc 0. After
  // distributing, all procs should have the same number of basis
  // functions.
  auto n_bf = _local_node_eim_basis_functions.size();
  sys.comm().broadcast(n_bf);

  // Allocate enough space to store n_bf basis functions on non-zero ranks
  if (rank != 0)
    _local_node_eim_basis_functions.resize(n_bf);

  // Broadcast the number of variables from proc 0. After
  // distributing, all procs should have the same number of variables.
  auto n_vars = _local_node_eim_basis_functions[0].begin()->second.size();
  sys.comm().broadcast(n_vars);

  // Construct lists of elem ids owned by different processors
  const MeshBase & mesh = sys.get_mesh();

  std::vector<dof_id_type> gathered_local_node_ids;
  {
    const std::set<boundary_id_type> & parametrized_function_boundary_ids =
      get_parametrized_function().get_parametrized_function_boundary_ids();

    const auto & binfo = mesh.get_boundary_info();

    // Make a set with all the nodes that have nodesets. Use
    // a set so that we don't have any duplicate entries. We
    // deal with duplicate entries below by getting all boundary
    // IDs on each node.
    std::set<dof_id_type> nodes_with_nodesets;
    for (const auto & t : binfo.build_node_list())
      nodes_with_nodesets.insert(std::get<0>(t));

    // To be filled in by BoundaryInfo calls in loop below
    std::vector<boundary_id_type> node_boundary_ids;

    for(dof_id_type node_id : nodes_with_nodesets)
      {
        const Node * node = mesh.node_ptr(node_id);

        if (node->processor_id() != mesh.comm().rank())
          continue;

        binfo.boundary_ids(node, node_boundary_ids);

        bool has_node_boundary_id = false;
        for(boundary_id_type node_boundary_id : node_boundary_ids)
          if(parametrized_function_boundary_ids.count(node_boundary_id))
            {
              has_node_boundary_id = true;
              break;
            }

        if(has_node_boundary_id)
          {
            gathered_local_node_ids.push_back(node_id);
          }
      }
  }

  // I _think_ the local node ids are likely to already be sorted in
  // ascending order, since that is how they are stored on the Mesh,
  // but we can always just guarantee this to be on the safe side as
  // well.
  std::sort(gathered_local_node_ids.begin(), gathered_local_node_ids.end());

  // Gather the number of local nodes from all procs to proc 0
  auto n_local_nodes = gathered_local_node_ids.size();
  std::vector<std::size_t> gathered_n_local_nodes = {n_local_nodes};
  sys.comm().gather(/*root_id=*/0, gathered_n_local_nodes);

  // Gather the node ids owned by each processor onto processor 0.
  sys.comm().gather(/*root_id=*/0, gathered_local_node_ids);

  // Construct vectors of "start" and "one-past-the-end" indices into
  // the gathered_local_node_ids vector for each proc. Only valid on
  // processor 0.
  std::vector<std::size_t> start_node_ids_index, end_node_ids_index;

  if (rank == 0)
    {
      start_node_ids_index.resize(n_procs);
      start_node_ids_index[0] = 0;
      for (processor_id_type p=1; p<n_procs; ++p)
        start_node_ids_index[p] = start_node_ids_index[p-1] + gathered_n_local_nodes[p-1];

      end_node_ids_index.resize(n_procs);
      end_node_ids_index[n_procs - 1] = gathered_local_node_ids.size();
      for (processor_id_type p=0; p<n_procs - 1; ++p)
        end_node_ids_index[p] = start_node_ids_index[p+1];
    }

  // On processor 0, using basis function 0 and variable 0, prepare a
  // vector with the nodes.  Then scatter this vector
  // out to the processors that require it. The order of this vector
  // matches the gathered_local_node_ids ordering. The counts will be
  // gathered_n_local_nodes.

  // On rank 0, the "counts" vector holds the number of floating point values that
  // are to be scattered to each proc. It is only required on proc 0.
  std::vector<int> counts;

  if (rank == 0)
    {
      counts.resize(n_procs);

      for (processor_id_type p=0; p<n_procs; ++p)
        {
          auto node_ids_range = (end_node_ids_index[p] - start_node_ids_index[p]);

          // Accumulate the count for this proc
          counts[p] += node_ids_range;
        } // end for proc_id
    } // if (rank == 0)

  // The recv_node_data vector will be used on the receiving end of all
  // the scatters below.
  std::vector<Number> recv_node_data;

  // For each basis function and each variable, build a vector
  // data in the Node ordering given by the
  // gathered_local_node_ids, then call
  //
  // sys.comm().scatter(data, counts, recv, /*root_id=*/0);
  std::vector<std::vector<Number>> node_data(n_vars);

  // We also reserve space in node_data, since we will push_back into it below.
  int count_sum = std::accumulate(counts.begin(), counts.end(), 0);
  for (auto var : index_range(node_data))
    node_data[var].reserve(count_sum);

  // Loop from 0..n_bf on _all_ procs, since the scatters inside this
  // loop are collective.
  for (auto bf : make_range(n_bf))
    {
      // Prepare data for scattering (only on proc 0)
      if (rank == 0)
        {
          // Reference to the data map for the current basis function.
          auto & bf_map = _local_node_eim_basis_functions[bf];

          // Clear any data from previous bf (this does not change the capacity
          // that was reserved above).
          for (auto var : index_range(node_data))
              node_data[var].clear();

          for (processor_id_type p=0; p<n_procs; ++p)
            {
              for (auto n : make_range(start_node_ids_index[p], end_node_ids_index[p]))
                {
                  auto node_id = gathered_local_node_ids[n];

                  // Get reference to array[n_vars] for current Node.
                  // Throws an error if the required node_id is not found.
                  const auto & array = libmesh_map_find(bf_map, node_id);

                  for (auto var : index_range(array))
                    node_data[var].push_back(array[var]);
                } // end for (n)
            } // end for proc_id
        } // end if rank==0

      // Perform the scatters (all procs)
      for (auto var : make_range(n_vars))
        {
          // Do the scatter for the current var
          sys.comm().scatter(node_data[var], counts, recv_node_data, /*root_id=*/0);

          if (rank != 0)
            {
              // Store the scattered data we received in _local_eim_basis_functions[bf]
              auto & bf_map = _local_node_eim_basis_functions[bf];
              auto cursor = recv_node_data.begin();

              for (auto i : index_range(gathered_local_node_ids))
                {
                  auto node_id = gathered_local_node_ids[i];
                  auto & array = bf_map[node_id];

                  // Create space to store the data if it doesn't already exist.
                  if (array.empty())
                    array.resize(n_vars);

                  // There is only one value per variable per node, so
                  // we set the value by de-referencing cursor, and
                  // then advance the cursor by 1.
                  array[var] = *cursor;
                  std::advance(cursor, 1);
                }
            } // if (rank != 0)
        } // end for (var)
    } // end for (bf)

  // Now that the scattering is done, delete non-local Elem
  // information from processor 0's _local_eim_basis_functions data
  // structure.
  if (rank == 0)
    {
      for (processor_id_type p=1; p<n_procs; ++p)
        {
          for (auto n : make_range(start_node_ids_index[p], end_node_ids_index[p]))
            {
              auto node_id = gathered_local_node_ids[n];

              // Delete this Node's information from every basis function.
              for (auto & bf_map : _local_node_eim_basis_functions)
                bf_map.erase(node_id);
            } // end for (n)
        } // end for proc_id
    } // if (rank == 0)
}

void RBEIMEvaluation::project_qp_data_vector_onto_system(System & /*sys*/,
                                                         const std::vector<Number> & /*bf_data*/,
                                                         const EIMVarGroupPlottingInfo & /*eim_vargroup*/,
                                                         const std::map<std::string,std::string> & /*extra_options*/)
{
  // No-op by default, implement in subclasses if needed
}

const std::vector<EIMVarGroupPlottingInfo> & RBEIMEvaluation::get_eim_vars_to_project_and_write() const
{
  return _eim_vars_to_project_and_write;
}

const std::set<unsigned int> & RBEIMEvaluation::scale_components_in_enrichment() const
{
  return _scale_components_in_enrichment;
}

bool RBEIMEvaluation::use_eim_error_indicator() const
{
  // Return false by default, but we override this in subclasses
  // for cases where we want to use the error indicator.
  return false;
}

void RBEIMEvaluation::set_eim_error_indicator_active(bool is_active)
{
  // We skip setting _is_eim_error_indicator_active in the case that
  // we have no parameters, since we do not use the EIM error indicator
  // in that case. We also check if the number of interpolation points
  // is larger than the number of EIM basis functions, since that is
  // also always the case when the error indicator is active.
  if ((get_n_params() > 0) && (get_n_interpolation_points() > get_n_basis_functions()))
    _is_eim_error_indicator_active = (is_active && use_eim_error_indicator());
}

std::pair<Real,Real> RBEIMEvaluation::get_eim_error_indicator(
  Number error_indicator_rhs,
  const DenseVector<Number> & eim_solution,
  const DenseVector<Number> & eim_rhs)
{
  DenseVector<Number> coeffs;
  _error_indicator_interpolation_row.get_principal_subvector(eim_solution.size(), coeffs);

  Number EIM_val_at_error_indicator_pt = coeffs.dot(eim_solution);
  Real error_indicator_val =
    std::real(error_indicator_rhs - EIM_val_at_error_indicator_pt);

  Real normalization = 0.;
  if (eim_error_indicator_normalization == RESIDUAL_SUM)
  {
    // This normalization is based on the sum of terms from the "EIM residual" calculation
    // used in the calculation of error_indicator_val. This ensures that the error indicator
    // will always be less than or equal to one, which is a useful property for an error
    // indicator.
    normalization =
      std::abs(error_indicator_rhs) + std::abs(EIM_val_at_error_indicator_pt);
  }
  else if (eim_error_indicator_normalization == RESIDUAL_RHS)
  {
    // Normalize with respect to the right-hand side from the "EIM residual" calculation.
    normalization = std::abs(error_indicator_rhs);
  }
  else if (eim_error_indicator_normalization == MAX_RHS)
  {
    // Normalize the error indicator based on the max-norm of the EIM RHS vector.
    // This approach handles the case where different EIM variables have different
    // magnitudes well, i.e. if error_indicator_val is based on a
    // "small magnitude" variable, then by normalizing based on the entire
    // RHS vector (which will typically include values from multiple different
    // EIM variables) we will effectively scale down the error indicator
    // corresponding to small variables, which is typically what we want.
    normalization = std::max(eim_rhs.linfty_norm(), std::abs(error_indicator_rhs));
  }
  else
  {
    libmesh_error_msg("unsupported eim_error_indicator_normalization");
  }

  // We avoid NaNs by setting normalization to 1 in the case that it is exactly 0.
  // But we return the "original normalization" as well (as opposed to the modified
  // normalization) since that can be useful information since it can indicate that
  // the EIM approximation was identically zero, for example.
  Real orig_normalization = normalization;
  if (normalization == 0.)
    normalization = 1.;

  // Return the relative error indicator, and the normalization that we used. By returning
  // the normalization, we can subsequently recover the absolute error indicator if
  // desired.
  //
  // We also optionally clamp the relative error indicator to 1.0 since this typically
  // indicators 100% error and hence it may be the maximum value that we want to see.
  Real rel_error_indicator = std::abs(error_indicator_val) / normalization;
  if (limit_eim_error_indicator_to_one && (rel_error_indicator > 1.0))
    rel_error_indicator = 1.0;

  return std::make_pair(rel_error_indicator, orig_normalization);
}

const VectorizedEvalInput & RBEIMEvaluation::get_vec_eval_input() const
{
  return _vec_eval_input;
}

const DenseVector<Number> & RBEIMEvaluation::get_error_indicator_interpolation_row() const
{
  return _error_indicator_interpolation_row;
}

void RBEIMEvaluation::set_error_indicator_interpolation_row(const DenseVector<Number> & extra_point_row)
{
  _error_indicator_interpolation_row = extra_point_row;
}

} // namespace libMesh
