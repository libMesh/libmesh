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

// C++ includes
#include <fstream>
#include <sstream>

// LibMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/equation_systems.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fem_context.h"
#include "libmesh/elem.h"
#include "libmesh/int_range.h"

// rbOOmit includes
#include "libmesh/rb_construction_base.h"
#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_parameters.h"
#include "libmesh/rb_parametrized_function.h"

// C++ include
#include <limits>
#include <memory>
#include <random>

namespace libMesh
{

namespace
{
// We make an anonymous namespace for local helper functions. The helper functions below
// are used to do some basic operators for QpDataMap and SideQpDataMap.

// Recall that we have:
// typedef std::map<dof_id_type, std::vector<std::vector<Number>>> QpDataMap;
// typedef std::map<std::pair<dof_id_type,unsigned int>, std::vector<std::vector<Number>>> SideQpDataMap;

// Implement u <- u + k*v
template <typename DataMap>
void add(DataMap & u, const Number k, const DataMap & v)
{
  for (auto & [key, vec_vec_u] : u)
    {
      const std::vector<std::vector<Number>> & vec_vec_v = libmesh_map_find(v, key);

      libmesh_error_msg_if (vec_vec_u.size() != vec_vec_v.size(), "Size mismatch");

      for (auto i : index_range(vec_vec_u))
        {
          libmesh_error_msg_if (vec_vec_u[i].size() != vec_vec_v[i].size(), "Size mismatch");
          for (auto j : index_range(vec_vec_u[i]))
            {
              vec_vec_u[i][j] += k*vec_vec_v[i][j];
            }
        }
    }
}

void add_node_data_map(RBEIMConstruction::NodeDataMap & u, const Number k, const RBEIMConstruction::NodeDataMap & v)
{
  for (auto & [key, vec_u] : u)
    {
      const std::vector<Number> & vec_v = libmesh_map_find(v, key);

      libmesh_error_msg_if (vec_u.size() != vec_v.size(), "Size mismatch");

      for (auto i : index_range(vec_u))
        {
          vec_u[i] += k*vec_v[i];
        }
    }
}

}

RBEIMConstruction::RBEIMConstruction (EquationSystems & es,
                                      const std::string & name_in,
                                      const unsigned int number_in)
  : RBConstructionBase(es, name_in, number_in),
    best_fit_type_flag(PROJECTION_BEST_FIT),
    _Nmax(0),
    _set_Nmax_from_n_snapshots(false),
    _Nmax_from_n_snapshots_increment(0),
    _rel_training_tolerance(1.e-4),
    _abs_training_tolerance(1.e-12),
    _max_abs_value_in_training_set(0.),
    _max_abs_value_in_training_set_index(0)
{
  // The training set should be the same on all processors in the
  // case of EIM training.
  serial_training_set = true;
}

RBEIMConstruction::~RBEIMConstruction () = default;

void RBEIMConstruction::clear()
{
  RBConstructionBase::clear();

  _rb_eim_assembly_objects.clear();

  _local_parametrized_functions_for_training.clear();
  _local_quad_point_locations.clear();
  _local_quad_point_JxW.clear();
  _local_quad_point_subdomain_ids.clear();

  _local_side_parametrized_functions_for_training.clear();
  _local_side_quad_point_locations.clear();
  _local_side_quad_point_JxW.clear();
  _local_side_quad_point_subdomain_ids.clear();
  _local_side_quad_point_boundary_ids.clear();
  _local_side_quad_point_side_types.clear();

  _local_node_parametrized_functions_for_training.clear();
  _local_node_locations.clear();
  _local_node_boundary_ids.clear();

  _eim_projection_matrix.resize(0,0);
}

void RBEIMConstruction::set_rb_eim_evaluation(RBEIMEvaluation & rb_eim_eval_in)
{
  _rb_eim_eval = &rb_eim_eval_in;
}

RBEIMEvaluation & RBEIMConstruction::get_rb_eim_evaluation()
{
  libmesh_error_msg_if(!_rb_eim_eval, "Error: RBEIMEvaluation object hasn't been initialized yet");
  return *_rb_eim_eval;
}

const RBEIMEvaluation & RBEIMConstruction::get_rb_eim_evaluation() const
{
  libmesh_error_msg_if(!_rb_eim_eval, "Error: RBEIMEvaluation object hasn't been initialized yet");
  return *_rb_eim_eval;
}

void RBEIMConstruction::set_best_fit_type_flag (const std::string & best_fit_type_string)
{
  if (best_fit_type_string == "projection")
    {
      best_fit_type_flag = PROJECTION_BEST_FIT;
    }
  else if (best_fit_type_string == "eim")
    {
      best_fit_type_flag = EIM_BEST_FIT;
    }
  else if (best_fit_type_string == "pod")
    {
      best_fit_type_flag = POD_BEST_FIT;
    }
  else
    libmesh_error_msg("Error: invalid best_fit_type in input file");
}

void RBEIMConstruction::print_info()
{
  // Print out info that describes the current setup
  libMesh::out << std::endl << "RBEIMConstruction parameters:" << std::endl;
  libMesh::out << "system name: " << this->name() << std::endl;
  libMesh::out << "Nmax: " << get_Nmax() << std::endl;
  if (_set_Nmax_from_n_snapshots)
  {
    libMesh::out << "Overruling Nmax based on number of snapshots, with increment set to "
      << _Nmax_from_n_snapshots_increment
      << std::endl;
  }
  libMesh::out << "Greedy relative error tolerance: " << get_rel_training_tolerance() << std::endl;
  libMesh::out << "Greedy absolute error tolerance: " << get_abs_training_tolerance() << std::endl;
  libMesh::out << "Number of parameters: " << get_n_params() << std::endl;
  for (const auto & pr : get_parameters())
    if (!is_discrete_parameter(pr.first))
      {
        libMesh::out <<   "Parameter " << pr.first
                     << ": Min = " << get_parameter_min(pr.first)
                     << ", Max = " << get_parameter_max(pr.first) << std::endl;
      }

  print_discrete_parameter_values();
  libMesh::out << "n_training_samples: " << get_n_training_samples() << std::endl;
  libMesh::out << "quiet mode? " << is_quiet() << std::endl;

  if (best_fit_type_flag == PROJECTION_BEST_FIT)
    {
      libMesh::out << "EIM best fit type: projection" << std::endl;
    }
  else
    if (best_fit_type_flag == EIM_BEST_FIT)
      {
        libMesh::out << "EIM best fit type: eim" << std::endl;
      }
  libMesh::out << std::endl;
}

void RBEIMConstruction::initialize_eim_construction()
{
  initialize_parametrized_functions_in_training_set();
}

void RBEIMConstruction::process_parameters_file (const std::string & parameters_filename)
{
  // First read in data from input_filename
  GetPot infile(parameters_filename);

  std::string best_fit_type_string = infile("best_fit_type","projection");
  set_best_fit_type_flag(best_fit_type_string);

  const unsigned int n_training_samples = infile("n_training_samples",0);
  const bool deterministic_training = infile("deterministic_training",false);
  unsigned int training_parameters_random_seed_in =
    static_cast<unsigned int>(-1);
  training_parameters_random_seed_in = infile("training_parameters_random_seed",
                                              training_parameters_random_seed_in);
  const bool quiet_mode_in = infile("quiet_mode", quiet_mode);
  const unsigned int Nmax_in = infile("Nmax", _Nmax);
  const Real rel_training_tolerance_in = infile("rel_training_tolerance",
                                                _rel_training_tolerance);
  const Real abs_training_tolerance_in = infile("abs_training_tolerance",
                                                _abs_training_tolerance);

  // Read in the parameters from the input file too
  unsigned int n_continuous_parameters = infile.vector_variable_size("parameter_names");
  RBParameters mu_min_in;
  RBParameters mu_max_in;
  for (unsigned int i=0; i<n_continuous_parameters; i++)
    {
      // Read in the parameter names
      std::string param_name = infile("parameter_names", "NONE", i);

      {
        Real min_val = infile(param_name, 0., 0);
        mu_min_in.set_value(param_name, min_val);
      }

      {
        Real max_val = infile(param_name, 0., 1);
        mu_max_in.set_value(param_name, max_val);
      }
    }

  std::map<std::string, std::vector<Real>> discrete_parameter_values_in;

  unsigned int n_discrete_parameters = infile.vector_variable_size("discrete_parameter_names");
  for (unsigned int i=0; i<n_discrete_parameters; i++)
    {
      std::string param_name = infile("discrete_parameter_names", "NONE", i);

      unsigned int n_vals_for_param = infile.vector_variable_size(param_name);
      std::vector<Real> vals_for_param(n_vals_for_param);
      for (auto j : make_range(vals_for_param.size()))
        vals_for_param[j] = infile(param_name, 0., j);

      discrete_parameter_values_in[param_name] = vals_for_param;
    }

  std::map<std::string,bool> log_scaling_in;
  // For now, just set all entries to false.
  // TODO: Implement a decent way to specify log-scaling true/false
  // in the input text file
  for (const auto & pr : mu_min_in)
    log_scaling_in[pr.first] = false;

  // Set the parameters that have been read in
  set_rb_construction_parameters(n_training_samples,
                                 deterministic_training,
                                 training_parameters_random_seed_in,
                                 quiet_mode_in,
                                 Nmax_in,
                                 rel_training_tolerance_in,
                                 abs_training_tolerance_in,
                                 mu_min_in,
                                 mu_max_in,
                                 discrete_parameter_values_in,
                                 log_scaling_in);
}

void RBEIMConstruction::set_rb_construction_parameters(unsigned int n_training_samples_in,
                                                       bool deterministic_training_in,
                                                       int training_parameters_random_seed_in,
                                                       bool quiet_mode_in,
                                                       unsigned int Nmax_in,
                                                       Real rel_training_tolerance_in,
                                                       Real abs_training_tolerance_in,
                                                       const RBParameters & mu_min_in,
                                                       const RBParameters & mu_max_in,
                                                       const std::map<std::string, std::vector<Real>> & discrete_parameter_values_in,
                                                       const std::map<std::string,bool> & log_scaling_in,
                                                       std::map<std::string, std::vector<RBParameter>> * training_sample_list)
{
  // Read in training_parameters_random_seed value.  This is used to
  // seed the RNG when picking the training parameters.  By default the
  // value is -1, which means use std::time to seed the RNG.
  set_training_random_seed(training_parameters_random_seed_in);

  // Set quiet mode
  set_quiet_mode(quiet_mode_in);

  // Initialize RB parameters
  set_Nmax(Nmax_in);

  set_rel_training_tolerance(rel_training_tolerance_in);
  set_abs_training_tolerance(abs_training_tolerance_in);

  if (get_rb_eim_evaluation().get_parametrized_function().is_lookup_table)
    {
      const std::string & lookup_table_param_name =
        get_rb_eim_evaluation().get_parametrized_function().lookup_table_param_name;

      libmesh_error_msg_if(!discrete_parameter_values_in.count(lookup_table_param_name),
        "Lookup table parameter should be discrete");

      // Make an editable copy of discrete_parameters_values_in.
      std::map<std::string, std::vector<Real>> discrete_parameter_values_final(
          discrete_parameter_values_in);

      std::vector<Real> & lookup_table_param_values =
        libmesh_map_find(discrete_parameter_values_final, lookup_table_param_name);

      // Overwrite the discrete values for lookup_table_param to make sure that
      // it is: 0, 1, 2, ..., size-1.
      std::iota(lookup_table_param_values.begin(), lookup_table_param_values.end(), 0);

      // Also, overwrite n_training_samples_in to make sure it matches
      // lookup_table_size so that we will get full coverage of the
      // lookup table in our training set.
      n_training_samples_in = lookup_table_param_values.size();

      // Initialize the parameter ranges and the parameters themselves
      initialize_parameters(mu_min_in, mu_max_in, discrete_parameter_values_final);
    }
  else
    {
      // Initialize the parameter ranges and the parameters themselves
      initialize_parameters(mu_min_in, mu_max_in, discrete_parameter_values_in);
    }

  bool updated_deterministic_training = deterministic_training_in;
  if (training_sample_list && (this->get_parameters_min().n_parameters() > 3))
    {
      // In this case we force deterministic_training to be false because
      // a) deterministic training samples are not currrently supported with
      //    more than 3 parameters, and
      // b) we will overwrite the training samples anyway in the call to
      //    load_training_set() below, so we do not want to generate an
      //    error due to deterministic training sample generation when
      //    the samples will be overwritten anyway.
      updated_deterministic_training = false;
    }

  initialize_training_parameters(this->get_parameters_min(),
                                 this->get_parameters_max(),
                                 n_training_samples_in,
                                 log_scaling_in,
                                 updated_deterministic_training);

  if (training_sample_list)
    {
      // Note that we must call initialize_training_parameters() before
      // load_training_set() in order to initialize the parameter vectors.
      load_training_set(*training_sample_list);
    }


  if (get_rb_eim_evaluation().get_parametrized_function().is_lookup_table)
    {
      // Also, now that we've initialized the training set, overwrite the training
      // samples to ensure that we have full coverage of the lookup tbale.
      const std::string & lookup_table_param_name =
        get_rb_eim_evaluation().get_parametrized_function().lookup_table_param_name;

      // Fill the lookup_table_training_samples with sequential single-entry vectors,
      // i.e. {{0.0}, {1.0}, {2.0}, ...}
      Real val = 0.0;
      std::vector<RBParameter> lookup_table_training_samples(n_training_samples_in, {val});
      for (auto & vec : lookup_table_training_samples)
        {
          vec[0] = val;
          val += 1.0;   // Could use val++, but better to be explicit for doubles.
        }

      set_training_parameter_values(lookup_table_param_name, lookup_table_training_samples);
    }
}

Real RBEIMConstruction::train_eim_approximation()
{
  if (_normalize_solution_snapshots)
    apply_normalization_to_solution_snapshots();

  if(best_fit_type_flag == POD_BEST_FIT)
    {
      train_eim_approximation_with_POD();
      return 0.;
    }
  else
    {
      return train_eim_approximation_with_greedy();
    }
}

Real RBEIMConstruction::train_eim_approximation_with_greedy()
{
  LOG_SCOPE("train_eim_approximation_with_greedy()", "RBEIMConstruction");

  RBEIMEvaluation & rbe = get_rb_eim_evaluation();

  // We need space for one extra interpolation point if we're using the
  // EIM error indicator.
  unsigned int max_matrix_size = rbe.use_eim_error_indicator() ? get_Nmax()+1 : get_Nmax();
  _eim_projection_matrix.resize(max_matrix_size,max_matrix_size);

  rbe.initialize_parameters(*this);
  rbe.resize_data_structures(max_matrix_size);

  // If we are continuing from a previous training run,
  // we might already be at the max number of basis functions.
  // If so, we can just return.
  libmesh_error_msg_if(rbe.get_n_basis_functions() > 0,
                       "Error: We currently only support EIM training starting from an empty basis");

  libMesh::out << std::endl << "---- Performing Greedy EIM basis enrichment ----" << std::endl;

  // Initialize greedy_error so that we do not incorrectly set is_zero_bf=true on
  // the first iteration.
  Real greedy_error = -1.;
  std::vector<RBParameters> greedy_param_list;

  // Initialize the current training index to the index that corresponds
  // to the largest (in terms of infinity norm) function in the training set.
  // We do this to ensure that the first EIM basis function is not zero.
  unsigned int current_training_index = _max_abs_value_in_training_set_index;
  set_params_from_training_set(current_training_index);

  // We use this boolean to indicate if we will run one more iteration
  // before exiting the loop below. We use this when computing the EIM
  // error indicator, which requires one extra EIM iteration.
  bool exit_on_next_iteration = false;

  // We also initialize a boolean to keep track of whether we have
  // reached "n_samples" EIM basis functions, since we need to
  // handle the EIM error indicator in a special way in this case.
  bool bfs_equals_n_samples = false;

  while (true)
    {
      if (rbe.get_n_basis_functions() >= get_n_training_samples())
        {
          libMesh::out << "Number of basis functions (" << rbe.get_n_basis_functions()
            << ") equals number of training samples." << std::endl;

          bfs_equals_n_samples = true;

          // If exit_on_next_iteration==true then we don't exit yet, since
          // we still need to add data for the error indicator before exiting.
          if (!exit_on_next_iteration)
            break;
        }

      libMesh::out << "Greedily selected parameter vector:" << std::endl;
      print_parameters();
      greedy_param_list.emplace_back(get_parameters());

      libMesh::out << "Enriching the EIM approximation" << std::endl;
      libmesh_try
        {
          bool is_zero_bf = bfs_equals_n_samples || (greedy_error == 0.);

          // If is_zero_bf==true then we add an "extra point" because we
          // cannot add a usual EIM interpolation point in that case since
          // the full EIM space is already covered. This is necessary when we
          // want to add an extra point for error indicator purposes in the
          // is_zero_bf==true case, for example.
          std::unique_ptr<EimPointData> eim_point_data;
          if (is_zero_bf)
              eim_point_data = std::make_unique<EimPointData>(get_random_point_from_training_sample());

          // If exit_on_next_iteration==true then we do not add a basis function in
          // that case since in that case we only need to add data for the EIM error
          // indicator.
          enrich_eim_approximation(current_training_index,
                                   /*add_basis_function*/ !exit_on_next_iteration,
                                   eim_point_data.get());
          update_eim_matrices(/*set_error_indicator*/ exit_on_next_iteration);

          libMesh::out << std::endl << "---- Basis dimension: "
                      << rbe.get_n_basis_functions() << " ----" << std::endl;

          if (get_rb_eim_evaluation().get_parametrized_function().is_lookup_table &&
              best_fit_type_flag == EIM_BEST_FIT)
            {
              // If this is a lookup table and we're using "EIM best fit" then we
              // need to update the eim_solutions after each EIM enrichment so that
              // we can call rb_eim_eval.rb_eim_solve() from within compute_max_eim_error().
              store_eim_solutions_for_training_set();
            }

          libMesh::out << "Computing EIM error on training set" << std::endl;
          std::tie(greedy_error, current_training_index) = compute_max_eim_error();
          set_params_from_training_set(current_training_index);

          libMesh::out << "Maximum EIM error is " << greedy_error << std::endl << std::endl;
        }
#ifdef LIBMESH_ENABLE_EXCEPTIONS
      catch (const std::exception & e)
        {
          // If we hit an exception when performing the enrichment for the error indicator, then
          // we just continue and skip the error indicator. Otherwise we rethrow the exception.
          if (exit_on_next_iteration)
            {
              std::cout << "Exception occurred when enriching basis for error indicator hence we skip the error indicator in this case" << std::endl;
              break;
            }
          else
            throw;
        }
#endif

      if (exit_on_next_iteration)
        {
          libMesh::out << "Extra EIM iteration for error indicator is complete, hence exiting EIM training now" << std::endl;
          break;
        }

      // Convergence and/or termination tests
      {
        bool exit_condition_satisfied = false;

        if (rbe.get_n_basis_functions() >= this->get_Nmax())
          {
            libMesh::out << "Maximum number of basis functions reached: Nmax = "
                          << get_Nmax() << std::endl;
            exit_condition_satisfied = true;
          }

        // We consider the relative tolerance as relative to the maximum value in the training
        // set, since we assume that this maximum value provides a relevant scaling.
        if (!exit_condition_satisfied)
          if (greedy_error < (get_rel_training_tolerance() * get_max_abs_value_in_training_set()))
            {
              libMesh::out << "Relative error tolerance reached." << std::endl;
              exit_condition_satisfied = true;
            }

        if (!exit_condition_satisfied)
          if (greedy_error < get_abs_training_tolerance())
            {
              libMesh::out << "Absolute error tolerance reached." << std::endl;
              exit_condition_satisfied = true;
            }

        bool has_parameters = (get_parameters().n_parameters() > 0);
        if (!exit_condition_satisfied)
        {
          bool do_exit = false;
          // In the check for repeated parameters we have to make sure this isn't a case
          // with no parameters, since in that case we would always report repeated
          // parameters.
          for (auto & param : greedy_param_list)
            if (param == get_parameters() && has_parameters)
              {
                libMesh::out << "Exiting greedy because the same parameters were selected twice"
                             << std::endl;
                do_exit = true;
                break;
              }

          if (do_exit)
            exit_condition_satisfied = true;
        }

        if (exit_condition_satisfied)
          {
            // If we're using the EIM error indicator then we need to run
            // one extra EIM iteration since we use the extra EIM point
            // to obtain our error indicator. If we're not using the EIM
            // error indicator, then we just exit now.
            if (get_rb_eim_evaluation().use_eim_error_indicator() && has_parameters)
              {
                exit_on_next_iteration = true;
                libMesh::out << "EIM error indicator is active, hence we will run one extra EIM iteration before exiting"
                             << std::endl;
              }
            else
              break;
          }
      }
    } // end while(true)

  if (rbe.get_parametrized_function().is_lookup_table &&
      best_fit_type_flag != EIM_BEST_FIT)
    {
      // We only enter here if best_fit_type_flag != EIM_BEST_FIT because we
      // already called this above in the EIM_BEST_FIT case.
      store_eim_solutions_for_training_set();
    }

  return greedy_error;
}

void RBEIMConstruction::apply_normalization_to_solution_snapshots()
{
  LOG_SCOPE("apply_normalization_to_solution_snapshots()", "RBEIMConstruction");

  libMesh::out << "Normalizing solution snapshots" << std::endl;

  bool apply_comp_scaling = !get_rb_eim_evaluation().scale_components_in_enrichment().empty();
  unsigned int n_snapshots = get_n_training_samples();
  RBEIMEvaluation & rbe = get_rb_eim_evaluation();

  for (unsigned int i=0; i<n_snapshots; i++)
    {
      if (rbe.get_parametrized_function().on_mesh_sides())
        {
          Real norm_val = std::sqrt(std::real(side_inner_product(
            _local_side_parametrized_functions_for_training[i],
            _local_side_parametrized_functions_for_training[i],
            apply_comp_scaling)));

          if (norm_val > 0.)
            scale_parametrized_function(_local_side_parametrized_functions_for_training[i], 1./norm_val);
        }
      else if (rbe.get_parametrized_function().on_mesh_nodes())
        {
          Real norm_val = std::sqrt(std::real(node_inner_product(
            _local_node_parametrized_functions_for_training[i],
            _local_node_parametrized_functions_for_training[i],
            apply_comp_scaling)));

          if (norm_val > 0.)
            scale_node_parametrized_function(_local_node_parametrized_functions_for_training[i], 1./norm_val);
        }
      else
        {
          Real norm_val = std::sqrt(std::real(inner_product(
            _local_parametrized_functions_for_training[i],
            _local_parametrized_functions_for_training[i],
            apply_comp_scaling)));

          if (norm_val > 0.)
            scale_parametrized_function(_local_parametrized_functions_for_training[i], 1./norm_val);
        }
    }
}

Real RBEIMConstruction::train_eim_approximation_with_POD()
{
  LOG_SCOPE("train_eim_approximation_with_POD()", "RBEIMConstruction");

  RBEIMEvaluation & rbe = get_rb_eim_evaluation();

  unsigned int n_snapshots = get_n_training_samples();

  // If _set_Nmax_from_n_snapshots=true, then we overrule Nmax.
  if (_set_Nmax_from_n_snapshots)
  {
    int updated_Nmax = (static_cast<int>(n_snapshots) + _Nmax_from_n_snapshots_increment);

    // We only overrule _Nmax if updated_Nmax is positive, since if Nmax=0 then we'll skip
    // training here entirely, which is typically not what we want.
    if (updated_Nmax > 0)
      _Nmax = static_cast<unsigned int>(updated_Nmax);
  }

  // _eim_projection_matrix is not used in the POD case, but we resize it here in any case
  // to be consistent with what we do in train_eim_approximation_with_greedy().
  // We need space for one extra interpolation point if we're using the
  // EIM error indicator.
  unsigned int max_matrix_size = rbe.use_eim_error_indicator() ? get_Nmax()+1 : get_Nmax();
  _eim_projection_matrix.resize(max_matrix_size,max_matrix_size);

  rbe.initialize_parameters(*this);
  rbe.resize_data_structures(get_Nmax());

  libmesh_error_msg_if(rbe.get_n_basis_functions() > 0,
                       "Error: We currently only support EIM training starting from an empty basis");

  libMesh::out << std::endl << "---- Performing POD EIM basis enrichment ----" << std::endl;

  // Set up the POD "correlation matrix". This enables us to compute the POD via the
  // "method of snapshots", in which we compute a low rank representation of the
  // n_snapshots x n_snapshots matrix.
  DenseMatrix<Number> correlation_matrix(n_snapshots,n_snapshots);

  std::cout << "Start computing correlation matrix" << std::endl;

  bool apply_comp_scaling = !get_rb_eim_evaluation().scale_components_in_enrichment().empty();
  for (unsigned int i=0; i<n_snapshots; i++)
    {
      for (unsigned int j=0; j<=i; j++)
        {
          Number inner_prod = 0.;
          if (rbe.get_parametrized_function().on_mesh_sides())
            {
              inner_prod = side_inner_product(
                _local_side_parametrized_functions_for_training[i],
                _local_side_parametrized_functions_for_training[j],
                apply_comp_scaling);
            }
          else if (rbe.get_parametrized_function().on_mesh_nodes())
            {
              inner_prod = node_inner_product(
                _local_node_parametrized_functions_for_training[i],
                _local_node_parametrized_functions_for_training[j],
                apply_comp_scaling);
            }
          else
            {
              inner_prod = inner_product(
                _local_parametrized_functions_for_training[i],
                _local_parametrized_functions_for_training[j],
                apply_comp_scaling);
            }


          correlation_matrix(i,j) = inner_prod;
          if(i != j)
            {
              correlation_matrix(j,i) = libmesh_conj(inner_prod);
            }
        }

      // Print out every 10th row so that we can see the progress
      if ( (i+1) % 10 == 0)
        std::cout << "Finished row " << (i+1) << " of " << n_snapshots << std::endl;
    }
  std::cout << "Finished computing correlation matrix" << std::endl;

  // Compute SVD of correlation matrix.
  // Let Y = U S V^T, then the SVD below corresponds
  // to Y^T Y = V S U^T U S V^T = V S^2 V^T.
  // The POD basis we use is then given by U, which
  // we can compute via U = Y V S^{-1}, which is what
  // we compute below.
  //
  // Note that the formulation remains the same in the
  // case that we use a weighted inner product (as
  // in the case that we used apply_comp_scaling=true
  // when computing the correlation matrix), see (1.28)
  // from the lecture notes from Volkwein on POD for more
  // details.
  DenseVector<Real> sigma( n_snapshots );
  DenseMatrix<Number> U( n_snapshots, n_snapshots );
  DenseMatrix<Number> VT( n_snapshots, n_snapshots );
  correlation_matrix.svd(sigma, U, VT );

  // We use this boolean to indicate if we will run one more iteration
  // before exiting the loop below. We use this when computing the EIM
  // error indicator, which requires one extra EIM iteration.
  bool exit_on_next_iteration = false;

  // Add dominant vectors from the POD as basis functions.
  unsigned int j = 0;
  Real rel_err = 0.;

  // We also initialize a boolean to keep track of whether we have
  // reached "n_snapshots" EIM basis functions, since we need to
  // handle the EIM error indicator in a special way in this case.
  bool j_equals_n_snapshots = false;
  while (true)
    {
      bool exit_condition_satisfied = false;

      if ((j == 0) && (sigma(0) == 0.))
      {
        libMesh::out << "Terminating EIM POD with empty basis because first singular value is zero" << std::endl;
        exit_condition_satisfied = true;
        rel_err = 0.;
      }
      else if (j >= n_snapshots)
        {
          libMesh::out << "Number of basis functions equals number of training samples." << std::endl;
          exit_condition_satisfied = true;
          j_equals_n_snapshots = true;

          // In this case we set the rel. error to be zero since we've filled up the
          // entire space. We cannot use the formula below for rel_err since
          // sigma(n_snapshots) is not defined.
          rel_err = 0.;
        }
      else
        {
          // The "energy" error in the POD approximation is determined by the first omitted
          // singular value, i.e. sigma(j). We normalize by sigma(0), which gives the total
          // "energy", in order to obtain a relative error.
          rel_err = std::sqrt(sigma(j)) / std::sqrt(sigma(0));
        }

      if (exit_on_next_iteration)
        {
          libMesh::out << "Extra EIM iteration for error indicator is complete, POD error norm for extra iteration: " << rel_err << std::endl;
          break;
        }

      libMesh::out << "Number of basis functions: " << j
                   << ", POD error norm: " << rel_err << std::endl;

      if (!exit_condition_satisfied)
        if (j >= get_Nmax())
          {
            libMesh::out << "Maximum number of basis functions (" << j << ") reached." << std::endl;
            exit_condition_satisfied = true;
          }

      if (!exit_condition_satisfied)
        if (rel_err < get_rel_training_tolerance())
          {
            libMesh::out << "Training tolerance reached." << std::endl;
            exit_condition_satisfied = true;
          }

      if (exit_condition_satisfied)
        {
          // If we're using the EIM error indicator then we need to run
          // one extra EIM iteration since we use the extra EIM point
          // to obtain our error indicator. If we're not using the EIM
          // error indicator, then we just exit now.
          bool has_parameters = (get_parameters().n_parameters() > 0);
          if (get_rb_eim_evaluation().use_eim_error_indicator() && has_parameters)
            {
              exit_on_next_iteration = true;
              libMesh::out << "EIM error indicator is active, hence we will run one extra EIM iteration before exiting"
                            << std::endl;
            }
          else
            break;
        }

      bool is_zero_bf = j_equals_n_snapshots || (rel_err == 0.);
      if (rbe.get_parametrized_function().on_mesh_sides())
        {
          // Make a "zero clone" by copying to get the same data layout, and then scaling by zero
          SideQpDataMap v = _local_side_parametrized_functions_for_training[0];

          if (!is_zero_bf)
            {
              scale_parametrized_function(v, 0.);

              for ( unsigned int i=0; i<n_snapshots; ++i )
                add(v, U.el(i, j), _local_side_parametrized_functions_for_training[i] );

              Real norm_v = std::sqrt(sigma(j));
              scale_parametrized_function(v, 1./norm_v);
            }

          libmesh_try
            {
              // If is_zero_bf==true then we add an "extra point" because we cannot
              // add a usual EIM interpolation point in that case since the full EIM
              // space is already covered. This is necessary when we want to add an
              // extra point for error indicator purposes in the is_zero_bf==true
              // case, for example.
              std::unique_ptr<EimPointData> eim_point_data;
              if (is_zero_bf)
                  eim_point_data = std::make_unique<EimPointData>(get_random_point(v));

              // If exit_on_next_iteration==true then we do not add a basis function in
              // that case since in that case we only need to add data for the EIM error
              // indicator.
              bool is_linearly_dependent = enrich_eim_approximation_on_sides(v,
                                                                             /*add_basis_function*/ !exit_on_next_iteration,
                                                                             eim_point_data.get());

              if (is_linearly_dependent && !is_zero_bf)
                {
                  // In this case we detected that v is actually linearly dependent and that is_zero_bf
                  // was previously not correct --- it should have been true. We typically
                  // catch this earlier (e.g. by checking rel_err) but in some cases we do not catch
                  // this until we call the enrichment method. In this situation we update is_zero_bf
                  // to true and call the enrichment again.
                  is_zero_bf = true;
                  eim_point_data = std::make_unique<EimPointData>(get_random_point(v));

                  enrich_eim_approximation_on_sides(v,
                                                    /*add_basis_function*/ !exit_on_next_iteration,
                                                    eim_point_data.get());
                }

              update_eim_matrices(/*set_error_indicator*/ exit_on_next_iteration);
            }
#ifdef LIBMESH_ENABLE_EXCEPTIONS
          catch (const std::exception & e)
            {
              // If we hit an exception when performing the enrichment for the error indicator, then
              // we just continue and skip the error indicator. Otherwise we rethrow the exception.
              if (exit_on_next_iteration)
                {
                  std::cout << "Exception occurred when enriching basis for error indicator hence we skip the error indicator in this case" << std::endl;
                  break;
                }
              else
                throw;
            }
#endif
        }
      else if (rbe.get_parametrized_function().on_mesh_nodes())
        {
          // Make a "zero clone" by copying to get the same data layout, and then scaling by zero
          NodeDataMap v = _local_node_parametrized_functions_for_training[0];

          if (!is_zero_bf)
            {
              scale_node_parametrized_function(v, 0.);

              for ( unsigned int i=0; i<n_snapshots; ++i )
                add_node_data_map(v, U.el(i, j), _local_node_parametrized_functions_for_training[i] );

              Real norm_v = std::sqrt(sigma(j));
              scale_node_parametrized_function(v, 1./norm_v);
            }

          libmesh_try
            {
              // If is_zero_bf==true then we add an "extra point" because we cannot
              // add a usual EIM interpolation point in that case since the full EIM
              // space is already covered. This is necessary when we want to add an
              // extra point for error indicator purposes in the is_zero_bf==true
              // case, for example.
              std::unique_ptr<EimPointData> eim_point_data;
              if (is_zero_bf)
                  eim_point_data = std::make_unique<EimPointData>(get_random_point(v));

              // If exit_on_next_iteration==true then we do not add a basis function in
              // that case since in that case we only need to add data for the EIM error
              // indicator.
              bool is_linearly_dependent = enrich_eim_approximation_on_nodes(v,
                                                                             /*add_basis_function*/ !exit_on_next_iteration,
                                                                             eim_point_data.get());

              if (is_linearly_dependent && !is_zero_bf)
                {
                  // In this case we detected that v is actually linearly dependent and that is_zero_bf
                  // was previously not correct --- it should have been true. We typically
                  // catch this earlier (e.g. by checking rel_err) but in some cases we do not catch
                  // this until we call the enrichment method. In this situation we update is_zero_bf
                  // to true and call the enrichment again.
                  is_zero_bf = true;
                  eim_point_data = std::make_unique<EimPointData>(get_random_point(v));

                  enrich_eim_approximation_on_nodes(v,
                                                    /*add_basis_function*/ !exit_on_next_iteration,
                                                    eim_point_data.get());
                }

              update_eim_matrices(/*set_error_indicator*/ exit_on_next_iteration);
            }
#ifdef LIBMESH_ENABLE_EXCEPTIONS
          catch (const std::exception & e)
            {
              // If we hit an exception when performing the enrichment for the error indicator, then
              // we just continue and skip the error indicator. Otherwise we rethrow the exception.
              if (exit_on_next_iteration)
                {
                  std::cout << "Exception occurred when enriching basis for error indicator hence we skip the error indicator in this case" << std::endl;
                  break;
                }
              else
                throw;
            }
#endif
        }
      else
        {
          // Make a "zero clone" by copying to get the same data layout, and then scaling by zero
          QpDataMap v = _local_parametrized_functions_for_training[0];

          if (!is_zero_bf)
            {
              scale_parametrized_function(v, 0.);

              for ( unsigned int i=0; i<n_snapshots; ++i )
                add(v, U.el(i, j), _local_parametrized_functions_for_training[i] );

              Real norm_v = std::sqrt(sigma(j));
              scale_parametrized_function(v, 1./norm_v);
            }

          libmesh_try
            {
              // If is_zero_bf==true then we add an "extra point" because we cannot
              // add a usual EIM interpolation point in that case since the full EIM
              // space is already covered. This is necessary when we want to add an
              // extra point for error indicator purposes in the is_zero_bf==true
              // case, for example.
              std::unique_ptr<EimPointData> eim_point_data;
              if (is_zero_bf)
                  eim_point_data = std::make_unique<EimPointData>(get_random_point(v));

              // If exit_on_next_iteration==true then we do not add a basis function in
              // that case since in that case we only need to add data for the EIM error
              // indicator.
              bool is_linearly_dependent = enrich_eim_approximation_on_interiors(v,
                                                                                 /*add_basis_function*/ !exit_on_next_iteration,
                                                                                 eim_point_data.get());

              if (is_linearly_dependent && !is_zero_bf)
                {
                  // In this case we detected that v is actually linearly dependent and that is_zero_bf
                  // was previously not correct --- it should have been true. We typically
                  // catch this earlier (e.g. by checking rel_err) but in some cases we do not catch
                  // this until we call the enrichment method. In this situation we update is_zero_bf
                  // to true and call the enrichment again.
                  is_zero_bf = true;
                  eim_point_data = std::make_unique<EimPointData>(get_random_point(v));

                  enrich_eim_approximation_on_interiors(v,
                                                        /*add_basis_function*/ !exit_on_next_iteration,
                                                        eim_point_data.get());
                }

              update_eim_matrices(/*set_error_indicator*/ exit_on_next_iteration);
            }
#ifdef LIBMESH_ENABLE_EXCEPTIONS
          catch (const std::exception & e)
            {
              // If we hit an exception when performing the enrichment for the error indicator, then
              // we just continue and skip the error indicator. Otherwise we rethrow the exception.
              if (exit_on_next_iteration)
                {
                  std::cout << "Exception occurred when enriching basis for error indicator hence we skip the error indicator in this case" << std::endl;
                  break;
                }
              else
                throw;
            }
#endif
        }

      if (is_zero_bf)
        {
          // In this case we exit here instead of increment j and continuing because
          // if we've encountered a zero EIM basis function then we must not have
          // any more valid data to add.
          std::cout << "Zero basis function encountered, hence exiting." << std::endl;
          break;
        }

      j++;
    }
  libMesh::out << std::endl;

  return rel_err;
}

void RBEIMConstruction::initialize_eim_assembly_objects()
{
  _rb_eim_assembly_objects.clear();
  for (auto i : make_range(get_rb_eim_evaluation().get_n_basis_functions()))
    _rb_eim_assembly_objects.push_back(build_eim_assembly(i));
}

std::vector<std::unique_ptr<ElemAssembly>> & RBEIMConstruction::get_eim_assembly_objects()
{
  return _rb_eim_assembly_objects;
}

void RBEIMConstruction::init_context(FEMContext & c)
{
  // Pre-request FE data for all element dimensions present in the
  // mesh.  Note: we currently pre-request FE data for all variables
  // in the current system but in some cases that may be overkill, for
  // example if only variable 0 is used.
  const System & sys = c.get_system();
  const MeshBase & mesh = sys.get_mesh();

  for (unsigned int dim=1; dim<=3; ++dim)
    if (mesh.elem_dimensions().count(dim))
      for (auto var : make_range(sys.n_vars()))
      {
        auto fe = c.get_element_fe(var, dim);
        fe->get_JxW();
        fe->get_xyz();
        fe->get_phi();

        auto side_fe = c.get_side_fe(var, dim);
        side_fe->get_JxW();
        side_fe->get_xyz();
        side_fe->get_phi();
      }
}

void RBEIMConstruction::set_rel_training_tolerance(Real new_training_tolerance)
{
  _rel_training_tolerance = new_training_tolerance;
}

Real RBEIMConstruction::get_rel_training_tolerance()
{
  return _rel_training_tolerance;
}

void RBEIMConstruction::set_abs_training_tolerance(Real new_training_tolerance)
{
  _abs_training_tolerance = new_training_tolerance;
}

Real RBEIMConstruction::get_abs_training_tolerance()
{
  return _abs_training_tolerance;
}

unsigned int RBEIMConstruction::get_Nmax() const
{
  return _Nmax;
}

void RBEIMConstruction::set_Nmax(unsigned int Nmax)
{
  _Nmax = Nmax;
}

void RBEIMConstruction::enable_set_Nmax_from_n_snapshots(int increment)
{
  _set_Nmax_from_n_snapshots = true;
  _Nmax_from_n_snapshots_increment = increment;
}

void RBEIMConstruction::disable_set_Nmax_from_n_snapshots()
{
  _set_Nmax_from_n_snapshots = false;
  _Nmax_from_n_snapshots_increment = 0;
}

Real RBEIMConstruction::get_max_abs_value_in_training_set() const
{
  return _max_abs_value_in_training_set;
}

void RBEIMConstruction::store_eim_solutions_for_training_set()
{
  LOG_SCOPE("store_eim_solutions_for_training_set()", "RBEIMConstruction");

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  std::vector<DenseVector<Number>> & eim_solutions = get_rb_eim_evaluation().get_eim_solutions_for_training_set();
  eim_solutions.clear();
  eim_solutions.resize(get_n_training_samples());

  unsigned int RB_size = get_rb_eim_evaluation().get_n_basis_functions();

  for (auto i : make_range(get_n_training_samples()))
    {
      if (eim_eval.get_parametrized_function().on_mesh_sides())
        {
          const auto & local_side_pf = _local_side_parametrized_functions_for_training[i];

          if (RB_size > 0)
            {
              // Get the right-hand side vector for the EIM approximation
              // by sampling the parametrized function (stored in solution)
              // at the interpolation points.
              DenseVector<Number> EIM_rhs(RB_size);
              for (unsigned int j=0; j<RB_size; j++)
                {
                  EIM_rhs(j) =
                    RBEIMEvaluation::get_parametrized_function_side_value(comm(),
                                                                          local_side_pf,
                                                                          eim_eval.get_interpolation_points_elem_id(j),
                                                                          eim_eval.get_interpolation_points_side_index(j),
                                                                          eim_eval.get_interpolation_points_comp(j),
                                                                          eim_eval.get_interpolation_points_qp(j));
                }
              eim_solutions[i] = eim_eval.rb_eim_solve(EIM_rhs);
            }
        }
      else if (eim_eval.get_parametrized_function().on_mesh_nodes())
        {
          const auto & local_node_pf = _local_node_parametrized_functions_for_training[i];

          if (RB_size > 0)
            {
              // Get the right-hand side vector for the EIM approximation
              // by sampling the parametrized function (stored in solution)
              // at the interpolation points.
              DenseVector<Number> EIM_rhs(RB_size);
              for (unsigned int j=0; j<RB_size; j++)
                {
                  EIM_rhs(j) =
                    RBEIMEvaluation::get_parametrized_function_node_value(comm(),
                                                                          local_node_pf,
                                                                          eim_eval.get_interpolation_points_node_id(j),
                                                                          eim_eval.get_interpolation_points_comp(j));
                }
              eim_solutions[i] = eim_eval.rb_eim_solve(EIM_rhs);
            }
        }
      else
        {
          const auto & local_pf = _local_parametrized_functions_for_training[i];

          if (RB_size > 0)
            {
              // Get the right-hand side vector for the EIM approximation
              // by sampling the parametrized function (stored in solution)
              // at the interpolation points.
              DenseVector<Number> EIM_rhs(RB_size);
              for (unsigned int j=0; j<RB_size; j++)
                {
                  EIM_rhs(j) =
                    RBEIMEvaluation::get_parametrized_function_value(comm(),
                                                                    local_pf,
                                                                    eim_eval.get_interpolation_points_elem_id(j),
                                                                    eim_eval.get_interpolation_points_comp(j),
                                                                    eim_eval.get_interpolation_points_qp(j));
                }
              eim_solutions[i] = eim_eval.rb_eim_solve(EIM_rhs);
            }
        }
    }
}

const RBEIMEvaluation::QpDataMap & RBEIMConstruction::get_parametrized_function_from_training_set(unsigned int training_index) const
{
  libmesh_error_msg_if(training_index >= _local_parametrized_functions_for_training.size(),
                       "Invalid index: " << training_index);
  return _local_parametrized_functions_for_training[training_index];
}

const RBEIMEvaluation::SideQpDataMap & RBEIMConstruction::get_side_parametrized_function_from_training_set(unsigned int training_index) const
{
  libmesh_error_msg_if(training_index >= _local_side_parametrized_functions_for_training.size(),
                       "Invalid index: " << training_index);
  return _local_side_parametrized_functions_for_training[training_index];
}

const RBEIMEvaluation::NodeDataMap & RBEIMConstruction::get_node_parametrized_function_from_training_set(unsigned int training_index) const
{
  libmesh_error_msg_if(training_index >= _local_node_parametrized_functions_for_training.size(),
                       "Invalid index: " << training_index);
  return _local_node_parametrized_functions_for_training[training_index];
}

const std::unordered_map<dof_id_type, std::vector<Real> > & RBEIMConstruction::get_local_quad_point_JxW()
{
  return _local_quad_point_JxW;
}

const std::map<std::pair<dof_id_type,unsigned int>, std::vector<Real> > & RBEIMConstruction::get_local_side_quad_point_JxW()
{
  return _local_side_quad_point_JxW;
}

unsigned int RBEIMConstruction::get_n_parametrized_functions_for_training() const
{
  if (get_rb_eim_evaluation().get_parametrized_function().on_mesh_sides())
    return _local_side_parametrized_functions_for_training.size();
  else if (get_rb_eim_evaluation().get_parametrized_function().on_mesh_nodes())
    return _local_node_parametrized_functions_for_training.size();
  else
    return _local_parametrized_functions_for_training.size();
}

void RBEIMConstruction::reinit_eim_projection_matrix()
{
  RBEIMEvaluation & rbe = get_rb_eim_evaluation();

  // We need space for one extra interpolation point if we're using the
  // EIM error indicator.
  unsigned int max_matrix_size = rbe.use_eim_error_indicator() ? get_Nmax()+1 : get_Nmax();
  _eim_projection_matrix.resize(max_matrix_size,max_matrix_size);
}

std::pair<Real,unsigned int> RBEIMConstruction::compute_max_eim_error()
{
  LOG_SCOPE("compute_max_eim_error()", "RBEIMConstruction");

  if (get_n_params() == 0)
    {
      // Just return 0 if we have no parameters.
      return std::make_pair(0.,0);
    }

  // keep track of the maximum error
  unsigned int max_err_index = 0;
  Real max_err = 0.;

  libmesh_error_msg_if(get_n_training_samples() != get_local_n_training_samples(),
                       "Error: Training samples should be the same on all procs");

  const unsigned int RB_size = get_rb_eim_evaluation().get_n_basis_functions();

  if(best_fit_type_flag == PROJECTION_BEST_FIT)
    {
      for (auto training_index : make_range(get_n_training_samples()))
        {
          if (get_rb_eim_evaluation().get_parametrized_function().on_mesh_sides())
            {
              // Make a copy of the pre-computed solution for the specified training sample
              // since we will modify it below to compute the best fit error.
              SideQpDataMap solution_copy = _local_side_parametrized_functions_for_training[training_index];

              // Perform an L2 projection in order to find the best approximation to
              // the parametrized function from the current EIM space.
              DenseVector<Number> best_fit_rhs(RB_size);
              for (unsigned int i=0; i<RB_size; i++)
                {
                  best_fit_rhs(i) = side_inner_product(solution_copy,
                                                       get_rb_eim_evaluation().get_side_basis_function(i),
                                                       /*apply_comp_scaling*/ false);
                }

              // Now compute the best fit by an LU solve
              DenseMatrix<Number> RB_inner_product_matrix_N(RB_size);
              _eim_projection_matrix.get_principal_submatrix(RB_size, RB_inner_product_matrix_N);

              DenseVector<Number> best_fit_coeffs;
              RB_inner_product_matrix_N.lu_solve(best_fit_rhs, best_fit_coeffs);

              get_rb_eim_evaluation().side_decrement_vector(solution_copy, best_fit_coeffs);
              Real best_fit_error = get_max_abs_value(solution_copy);

              if (best_fit_error > max_err)
                {
                  max_err_index = training_index;
                  max_err = best_fit_error;
                }
            }
          else if (get_rb_eim_evaluation().get_parametrized_function().on_mesh_nodes())
            {
              // Make a copy of the pre-computed solution for the specified training sample
              // since we will modify it below to compute the best fit error.
              NodeDataMap solution_copy = _local_node_parametrized_functions_for_training[training_index];

              // Perform an L2 projection in order to find the best approximation to
              // the parametrized function from the current EIM space.
              DenseVector<Number> best_fit_rhs(RB_size);
              for (unsigned int i=0; i<RB_size; i++)
                {
                  best_fit_rhs(i) = node_inner_product(solution_copy,
                                                       get_rb_eim_evaluation().get_node_basis_function(i),
                                                       /*apply_comp_scaling*/ false);
                }

              // Now compute the best fit by an LU solve
              DenseMatrix<Number> RB_inner_product_matrix_N(RB_size);
              _eim_projection_matrix.get_principal_submatrix(RB_size, RB_inner_product_matrix_N);

              DenseVector<Number> best_fit_coeffs;
              RB_inner_product_matrix_N.lu_solve(best_fit_rhs, best_fit_coeffs);

              get_rb_eim_evaluation().node_decrement_vector(solution_copy, best_fit_coeffs);
              Real best_fit_error = get_node_max_abs_value(solution_copy);

              if (best_fit_error > max_err)
                {
                  max_err_index = training_index;
                  max_err = best_fit_error;
                }
            }
          else
            {
              // Make a copy of the pre-computed solution for the specified training sample
              // since we will modify it below to compute the best fit error.
              QpDataMap solution_copy = _local_parametrized_functions_for_training[training_index];

              // Perform an L2 projection in order to find the best approximation to
              // the parametrized function from the current EIM space.
              DenseVector<Number> best_fit_rhs(RB_size);
              for (unsigned int i=0; i<RB_size; i++)
                {
                  best_fit_rhs(i) = inner_product(solution_copy,
                                                  get_rb_eim_evaluation().get_basis_function(i),
                                                  /*apply_comp_scaling*/ false);
                }

              // Now compute the best fit by an LU solve
              DenseMatrix<Number> RB_inner_product_matrix_N(RB_size);
              _eim_projection_matrix.get_principal_submatrix(RB_size, RB_inner_product_matrix_N);

              DenseVector<Number> best_fit_coeffs;
              RB_inner_product_matrix_N.lu_solve(best_fit_rhs, best_fit_coeffs);

              get_rb_eim_evaluation().decrement_vector(solution_copy, best_fit_coeffs);
              Real best_fit_error = get_max_abs_value(solution_copy);

              if (best_fit_error > max_err)
                {
                  max_err_index = training_index;
                  max_err = best_fit_error;
                }
            }
        }
    }
  else if(best_fit_type_flag == EIM_BEST_FIT)
    {
      // Perform EIM solve in order to find the approximation to solution
      // (rb_eim_solve provides the EIM basis function coefficients used below)

      std::vector<RBParameters> training_parameters_copy(get_n_training_samples());
      for (auto training_index : make_range(get_n_training_samples()))
        {
          training_parameters_copy[training_index] = get_params_from_training_set(training_index);
        }

      get_rb_eim_evaluation().rb_eim_solves(training_parameters_copy, RB_size);
      const std::vector<DenseVector<Number>> & rb_eim_solutions = get_rb_eim_evaluation().get_rb_eim_solutions();

      for (auto training_index : make_range(get_n_training_samples()))
        {
          const DenseVector<Number> & best_fit_coeffs = rb_eim_solutions[training_index];

          if (get_rb_eim_evaluation().get_parametrized_function().on_mesh_sides())
            {
              SideQpDataMap solution_copy = _local_side_parametrized_functions_for_training[training_index];
              get_rb_eim_evaluation().side_decrement_vector(solution_copy, best_fit_coeffs);
              Real best_fit_error = get_max_abs_value(solution_copy);

              if (best_fit_error > max_err)
                {
                  max_err_index = training_index;
                  max_err = best_fit_error;
                }
            }
          else if (get_rb_eim_evaluation().get_parametrized_function().on_mesh_nodes())
            {
              NodeDataMap solution_copy = _local_node_parametrized_functions_for_training[training_index];
              get_rb_eim_evaluation().node_decrement_vector(solution_copy, best_fit_coeffs);
              Real best_fit_error = get_node_max_abs_value(solution_copy);

              if (best_fit_error > max_err)
                {
                  max_err_index = training_index;
                  max_err = best_fit_error;
                }
            }
          else
            {
              QpDataMap solution_copy = _local_parametrized_functions_for_training[training_index];
              get_rb_eim_evaluation().decrement_vector(solution_copy, best_fit_coeffs);
              Real best_fit_error = get_max_abs_value(solution_copy);

              if (best_fit_error > max_err)
                {
                  max_err_index = training_index;
                  max_err = best_fit_error;
                }
            }
        }
    }
  else
    {
      libmesh_error_msg("EIM best fit type not recognized");
    }

  return std::make_pair(max_err,max_err_index);
}

void RBEIMConstruction::initialize_parametrized_functions_in_training_set()
{
  LOG_SCOPE("initialize_parametrized_functions_in_training_set()", "RBEIMConstruction");

  libmesh_error_msg_if(!serial_training_set,
                       "Error: We must have serial_training_set==true in "
                       "RBEIMConstruction::initialize_parametrized_functions_in_training_set");

  libMesh::out << "Initializing parametrized functions in training set..." << std::endl;

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  if (eim_eval.get_parametrized_function().is_lookup_table)
    eim_eval.get_parametrized_function().initialize_lookup_table();

  // Store the locations of all quadrature points
  initialize_qp_data();

  // Keep track of the largest value in our parametrized functions
  // in the training set. We can use this value for normalization
  // purposes, for example.
  _max_abs_value_in_training_set = 0.;

  unsigned int n_comps = eim_eval.get_parametrized_function().get_n_components();

  // Keep track of the maximum value per component. This will allow
  // us to scale the components to all have a similar magnitude,
  // which is helpful during the error assessment for the basis
  // enrichment to ensure that components with smaller magnitude
  // are not ignored.
  std::vector<Real> max_abs_value_per_component_in_training_set(n_comps);

  if (eim_eval.get_parametrized_function().on_mesh_sides())
    {
      _local_side_parametrized_functions_for_training.resize( get_n_training_samples() );
      for (auto i : make_range(get_n_training_samples()))
        {
          libMesh::out << "Initializing parametrized function for training sample "
            << (i+1) << " of " << get_n_training_samples() << std::endl;

          set_params_from_training_set(i);

          eim_eval.get_parametrized_function().preevaluate_parametrized_function_on_mesh_sides(get_parameters(),
                                                                                               _local_side_quad_point_locations,
                                                                                               _local_side_quad_point_subdomain_ids,
                                                                                               _local_side_quad_point_boundary_ids,
                                                                                               _local_side_quad_point_side_types,
                                                                                               _local_side_quad_point_locations_perturbations,
                                                                                               *this);

          for (const auto & [elem_side_pair, xyz_vector] : _local_side_quad_point_locations)
          {
            std::vector<std::vector<Number>> comps_and_qps(n_comps);
            for (unsigned int comp=0; comp<n_comps; comp++)
              {
                comps_and_qps[comp].resize(xyz_vector.size());
                for (unsigned int qp : index_range(xyz_vector))
                  {
                    Number value =
                      eim_eval.get_parametrized_function().lookup_preevaluated_side_value_on_mesh(comp,
                                                                                                  elem_side_pair.first,
                                                                                                  elem_side_pair.second,
                                                                                                  qp);
                    comps_and_qps[comp][qp] = value;

                    Real abs_value = std::abs(value);
                    if (abs_value > _max_abs_value_in_training_set)
                      {
                        _max_abs_value_in_training_set = abs_value;
                        _max_abs_value_in_training_set_index = i;
                      }

                    if (abs_value > max_abs_value_per_component_in_training_set[comp])
                      max_abs_value_per_component_in_training_set[comp] = abs_value;
                  }
              }

            _local_side_parametrized_functions_for_training[i][elem_side_pair] = comps_and_qps;
          }
        }

      libMesh::out << "Parametrized functions in training set initialized" << std::endl;

      unsigned int max_id = 0;
      comm().maxloc(_max_abs_value_in_training_set, max_id);
      comm().broadcast(_max_abs_value_in_training_set_index, max_id);
      libMesh::out << "Maximum absolute value in the training set: "
        << _max_abs_value_in_training_set << std::endl << std::endl;

      // Calculate the maximum value for each component in the training set
      // across all components
      comm().max(max_abs_value_per_component_in_training_set);

      // We store the maximum value across all components divided by the maximum value for this component
      // so that when we scale using these factors all components should have a magnitude on the same
      // order as the maximum component.
      _component_scaling_in_training_set.resize(n_comps);
      for (unsigned int i : make_range(n_comps))
        {
          if ((eim_eval.scale_components_in_enrichment().count(i) == 0) ||
               max_abs_value_per_component_in_training_set[i] == 0.)
            _component_scaling_in_training_set[i] = 1.;
          else
            _component_scaling_in_training_set[i] = _max_abs_value_in_training_set / max_abs_value_per_component_in_training_set[i];
        }
    }
  else if (eim_eval.get_parametrized_function().on_mesh_nodes())
    {
      _local_node_parametrized_functions_for_training.resize( get_n_training_samples() );
      for (auto i : make_range(get_n_training_samples()))
        {
          libMesh::out << "Initializing parametrized function for training sample "
            << (i+1) << " of " << get_n_training_samples() << std::endl;

          set_params_from_training_set(i);

          eim_eval.get_parametrized_function().preevaluate_parametrized_function_on_mesh_nodes(get_parameters(),
                                                                                               _local_node_locations,
                                                                                               _local_node_boundary_ids,
                                                                                               *this);

          for (const auto & pr : _local_node_locations)
          {
            const auto & node_id = pr.first;

            std::vector<Number> comps(n_comps);
            for (unsigned int comp=0; comp<n_comps; comp++)
              {
                Number value =
                  eim_eval.get_parametrized_function().lookup_preevaluated_node_value_on_mesh(comp,
                                                                                              node_id);
                comps[comp] = value;

                Real abs_value = std::abs(value);
                if (abs_value > _max_abs_value_in_training_set)
                  {
                    _max_abs_value_in_training_set = abs_value;
                    _max_abs_value_in_training_set_index = i;
                  }

                if (abs_value > max_abs_value_per_component_in_training_set[comp])
                  max_abs_value_per_component_in_training_set[comp] = abs_value;
              }

            _local_node_parametrized_functions_for_training[i][node_id] = comps;
          }
        }

      libMesh::out << "Parametrized functions in training set initialized" << std::endl;

      unsigned int max_id = 0;
      comm().maxloc(_max_abs_value_in_training_set, max_id);
      comm().broadcast(_max_abs_value_in_training_set_index, max_id);
      libMesh::out << "Maximum absolute value in the training set: "
        << _max_abs_value_in_training_set << std::endl << std::endl;

      // Calculate the maximum value for each component in the training set
      // across all components
      comm().max(max_abs_value_per_component_in_training_set);

      // We store the maximum value across all components divided by the maximum value for this component
      // so that when we scale using these factors all components should have a magnitude on the same
      // order as the maximum component.
      _component_scaling_in_training_set.resize(n_comps);
      for (unsigned int i : make_range(n_comps))
        {
          if ((eim_eval.scale_components_in_enrichment().count(i) == 0) ||
               max_abs_value_per_component_in_training_set[i] == 0.)
            _component_scaling_in_training_set[i] = 1.;
          else
            _component_scaling_in_training_set[i] = _max_abs_value_in_training_set / max_abs_value_per_component_in_training_set[i];
        }
    }
  else
    {
      _local_parametrized_functions_for_training.resize( get_n_training_samples() );
      for (auto i : make_range(get_n_training_samples()))
        {
          libMesh::out << "Initializing parametrized function for training sample "
            << (i+1) << " of " << get_n_training_samples() << std::endl;

          set_params_from_training_set(i);

          eim_eval.get_parametrized_function().preevaluate_parametrized_function_on_mesh(get_parameters(),
                                                                                        _local_quad_point_locations,
                                                                                        _local_quad_point_subdomain_ids,
                                                                                        _local_quad_point_locations_perturbations,
                                                                                        *this);

          for (const auto & [elem_id, xyz_vector] : _local_quad_point_locations)
          {
            std::vector<std::vector<Number>> comps_and_qps(n_comps);
            for (unsigned int comp=0; comp<n_comps; comp++)
              {
                comps_and_qps[comp].resize(xyz_vector.size());
                for (unsigned int qp : index_range(xyz_vector))
                  {
                    Number value =
                      eim_eval.get_parametrized_function().lookup_preevaluated_value_on_mesh(comp, elem_id, qp);
                    comps_and_qps[comp][qp] = value;

                    Real abs_value = std::abs(value);
                    if (abs_value > _max_abs_value_in_training_set)
                      {
                        _max_abs_value_in_training_set = abs_value;
                        _max_abs_value_in_training_set_index = i;
                      }

                    if (abs_value > max_abs_value_per_component_in_training_set[comp])
                      max_abs_value_per_component_in_training_set[comp] = abs_value;
                  }
              }

            _local_parametrized_functions_for_training[i][elem_id] = comps_and_qps;
          }
        }

      libMesh::out << "Parametrized functions in training set initialized" << std::endl;

      unsigned int max_id = 0;
      comm().maxloc(_max_abs_value_in_training_set, max_id);
      comm().broadcast(_max_abs_value_in_training_set_index, max_id);
      libMesh::out << "Maximum absolute value in the training set: "
        << _max_abs_value_in_training_set << std::endl << std::endl;

      // Calculate the maximum value for each component in the training set
      // across all components
      comm().max(max_abs_value_per_component_in_training_set);

      // We store the maximum value across all components divided by the maximum value for this component
      // so that when we scale using these factors all components should have a magnitude on the same
      // order as the maximum component.
      _component_scaling_in_training_set.resize(n_comps);
      for (unsigned int i : make_range(n_comps))
        {
          if ((eim_eval.scale_components_in_enrichment().count(i) == 0) ||
               max_abs_value_per_component_in_training_set[i] == 0.)
            _component_scaling_in_training_set[i] = 1.;
          else
            _component_scaling_in_training_set[i] = _max_abs_value_in_training_set / max_abs_value_per_component_in_training_set[i];
        }
    }
    // This function does nothing if rb_property_map from RBParametrizedFunction
    // is empty which would result in an empty rb_property_map in VectorizedEvalInput
    // stored in RBEIMEvaluation.
    eim_eval.initialize_rb_property_map();
}

void RBEIMConstruction::initialize_qp_data()
{
  LOG_SCOPE("initialize_qp_data()", "RBEIMConstruction");

  if (!get_rb_eim_evaluation().get_parametrized_function().requires_xyz_perturbations)
    {
      libMesh::out << "Initializing quadrature point locations" << std::endl;
    }
  else
    {
      libMesh::out << "Initializing quadrature point and perturbation locations" << std::endl;
    }

  // Compute truth representation via L2 projection
  const MeshBase & mesh = this->get_mesh();

  FEMContext context(*this);
  init_context(context);

  if (get_rb_eim_evaluation().get_parametrized_function().on_mesh_sides())
    {
      const std::set<boundary_id_type> & parametrized_function_boundary_ids =
        get_rb_eim_evaluation().get_parametrized_function().get_parametrized_function_boundary_ids();
      libmesh_error_msg_if (parametrized_function_boundary_ids.empty(),
                            "Need to have non-empty boundary IDs to initialize side data");

      _local_side_quad_point_locations.clear();
      _local_side_quad_point_subdomain_ids.clear();
      _local_side_quad_point_boundary_ids.clear();
      _local_side_quad_point_JxW.clear();
      _local_side_quad_point_side_types.clear();

      _local_side_quad_point_locations_perturbations.clear();

      // BoundaryInfo and related data structures
      const auto & binfo = mesh.get_boundary_info();
      std::vector<boundary_id_type> side_boundary_ids;

      for (const auto & elem : mesh.active_local_element_ptr_range())
        {
          dof_id_type elem_id = elem->id();

          context.pre_fe_reinit(*this, elem);

          // elem_fe is used for shellface data
          auto elem_fe = context.get_element_fe(/*var=*/0, elem->dim());
          const std::vector<Real> & JxW = elem_fe->get_JxW();
          const std::vector<Point> & xyz = elem_fe->get_xyz();

          // side_fe is used for element side data
          auto side_fe = context.get_side_fe(/*var=*/0, elem->dim());
          const std::vector<Real> & JxW_side = side_fe->get_JxW();
          const std::vector< Point > & xyz_side = side_fe->get_xyz();

          for (context.side = 0;
               context.side != context.get_elem().n_sides();
               ++context.side)
            {
              // skip non-boundary elements
              if(!context.get_elem().neighbor_ptr(context.side))
                {
                  binfo.boundary_ids(elem, context.side, side_boundary_ids);

                  bool has_side_boundary_id = false;
                  boundary_id_type matching_boundary_id = BoundaryInfo::invalid_id;
                  for (boundary_id_type side_boundary_id : side_boundary_ids)
                    if(parametrized_function_boundary_ids.count(side_boundary_id))
                      {
                        has_side_boundary_id = true;
                        matching_boundary_id = side_boundary_id;
                        break;
                      }

                  if(has_side_boundary_id)
                  {
                    context.get_side_fe(/*var=*/0, elem->dim())->reinit(elem, context.side);

                    auto elem_side_pair = std::make_pair(elem_id, context.side);

                    _local_side_quad_point_locations[elem_side_pair] = xyz_side;
                    _local_side_quad_point_JxW[elem_side_pair] = JxW_side;
                    _local_side_quad_point_subdomain_ids[elem_side_pair] = elem->subdomain_id();
                    _local_side_quad_point_boundary_ids[elem_side_pair] = matching_boundary_id;

                    // This is a standard side (not a shellface) so set side type to 0
                    _local_side_quad_point_side_types[elem_side_pair] = 0;

                    if (get_rb_eim_evaluation().get_parametrized_function().requires_xyz_perturbations)
                      {
                        Real fd_delta = get_rb_eim_evaluation().get_parametrized_function().fd_delta;

                        std::vector<std::vector<Point>> xyz_perturb_vec_at_qps;

                        for (const Point & xyz_qp : xyz_side)
                          {
                            std::vector<Point> xyz_perturb_vec;
                            if (elem->dim() == 3)
                              {
                                // In this case we have a 3D element, and hence the side is 2D.
                                //
                                // We use the following approach to perturb xyz:
                                //  1) inverse map xyz to the reference element
                                //  2) perturb on the reference element in the (xi,eta) "directions"
                                //  3) map the perturbed points back to the physical element
                                // This approach is necessary to ensure that the perturbed points
                                // are still in the element's side.

                                std::unique_ptr<const Elem> elem_side;
                                elem->build_side_ptr(elem_side, context.side);

                                Point xi_eta =
                                  FEMap::inverse_map(elem_side->dim(),
                                                     elem_side.get(),
                                                     xyz_qp,
                                                     /*Newton iteration tolerance*/ TOLERANCE,
                                                     /*secure*/ true);

                                // Inverse map should map back to a 2D reference domain
                                libmesh_assert(std::abs(xi_eta(2)) < TOLERANCE);

                                Point xi_eta_perturb = xi_eta;

                                xi_eta_perturb(0) += fd_delta;
                                Point xyz_perturb_0 =
                                  FEMap::map(elem_side->dim(),
                                             elem_side.get(),
                                             xi_eta_perturb);
                                xi_eta_perturb(0) -= fd_delta;

                                xi_eta_perturb(1) += fd_delta;
                                Point xyz_perturb_1 =
                                  FEMap::map(elem_side->dim(),
                                             elem_side.get(),
                                             xi_eta_perturb);
                                xi_eta_perturb(1) -= fd_delta;

                                // Finally, we rescale xyz_perturb_0 and xyz_perturb_1 so that
                                // (xyz_perturb - xyz_qp).norm() == fd_delta, since this is
                                // required in order to compute finite differences correctly.
                                Point unit_0 = (xyz_perturb_0-xyz_qp).unit();
                                Point unit_1 = (xyz_perturb_1-xyz_qp).unit();

                                xyz_perturb_vec.emplace_back(xyz_qp + fd_delta*unit_0);
                                xyz_perturb_vec.emplace_back(xyz_qp + fd_delta*unit_1);
                              }
                            else
                              {
                                // We current do nothing for sides of dim=2 or dim=1 elements
                                // since we have no need for this capability so far.
                                // Support for these cases could be added if it is needed.
                              }

                            xyz_perturb_vec_at_qps.emplace_back(xyz_perturb_vec);
                          }

                        _local_side_quad_point_locations_perturbations[elem_side_pair] = xyz_perturb_vec_at_qps;
                      }
                  }
                }
            }

            // In the case of 2D elements, we also check the shellfaces
            if (elem->dim() == 2)
              for (unsigned int shellface_index=0; shellface_index<2; shellface_index++)
                {
                  binfo.shellface_boundary_ids(elem, shellface_index, side_boundary_ids);

                  bool has_side_boundary_id = false;
                  boundary_id_type matching_boundary_id = BoundaryInfo::invalid_id;
                  for (boundary_id_type side_boundary_id : side_boundary_ids)
                    if(parametrized_function_boundary_ids.count(side_boundary_id))
                      {
                        has_side_boundary_id = true;
                        matching_boundary_id = side_boundary_id;
                        break;
                      }

                  if(has_side_boundary_id)
                  {
                    context.elem_fe_reinit();

                    // We use shellface_index as the side_index since shellface boundary conditions
                    // are stored separately from side boundary conditions in BoundaryInfo.
                    auto elem_side_pair = std::make_pair(elem_id, shellface_index);

                    _local_side_quad_point_locations[elem_side_pair] = xyz;
                    _local_side_quad_point_JxW[elem_side_pair] = JxW;
                    _local_side_quad_point_subdomain_ids[elem_side_pair] = elem->subdomain_id();
                    _local_side_quad_point_boundary_ids[elem_side_pair] = matching_boundary_id;

                    // This is a shellface (not a standard side) so set side type to 1
                    _local_side_quad_point_side_types[elem_side_pair] = 1;

                    if (get_rb_eim_evaluation().get_parametrized_function().requires_xyz_perturbations)
                      {
                        Real fd_delta = get_rb_eim_evaluation().get_parametrized_function().fd_delta;

                        std::vector<std::vector<Point>> xyz_perturb_vec_at_qps;

                        for (const Point & xyz_qp : xyz)
                          {
                            std::vector<Point> xyz_perturb_vec;
                            // Here we follow the same approach as above for getting xyz_perturb_vec,
                            // except that we are using the element itself instead of its side.
                            {
                              Point xi_eta =
                                FEMap::inverse_map(elem->dim(),
                                                   elem,
                                                   xyz_qp,
                                                   /*Newton iteration tolerance*/ TOLERANCE,
                                                   /*secure*/ true);

                              // Inverse map should map back to a 2D reference domain
                              libmesh_assert(std::abs(xi_eta(2)) < TOLERANCE);

                              Point xi_eta_perturb = xi_eta;

                              xi_eta_perturb(0) += fd_delta;
                              Point xyz_perturb_0 =
                                FEMap::map(elem->dim(),
                                            elem,
                                            xi_eta_perturb);
                              xi_eta_perturb(0) -= fd_delta;

                              xi_eta_perturb(1) += fd_delta;
                              Point xyz_perturb_1 =
                                FEMap::map(elem->dim(),
                                            elem,
                                            xi_eta_perturb);
                              xi_eta_perturb(1) -= fd_delta;

                              // Finally, we rescale xyz_perturb_0 and xyz_perturb_1 so that
                              // (xyz_perturb - xyz_qp).norm() == fd_delta, since this is
                              // required in order to compute finite differences correctly.
                              Point unit_0 = (xyz_perturb_0-xyz_qp).unit();
                              Point unit_1 = (xyz_perturb_1-xyz_qp).unit();

                              xyz_perturb_vec.emplace_back(xyz_qp + fd_delta*unit_0);
                              xyz_perturb_vec.emplace_back(xyz_qp + fd_delta*unit_1);
                            }

                            xyz_perturb_vec_at_qps.emplace_back(xyz_perturb_vec);
                          }

                        _local_side_quad_point_locations_perturbations[elem_side_pair] = xyz_perturb_vec_at_qps;
                      }
                  }
                }
        }
    }
  else if (get_rb_eim_evaluation().get_parametrized_function().on_mesh_nodes())
    {
      const std::set<boundary_id_type> & parametrized_function_boundary_ids =
        get_rb_eim_evaluation().get_parametrized_function().get_parametrized_function_boundary_ids();
      libmesh_error_msg_if (parametrized_function_boundary_ids.empty(),
                            "Need to have non-empty boundary IDs to initialize node data");

      _local_node_locations.clear();
      _local_node_boundary_ids.clear();

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

      for (dof_id_type node_id : nodes_with_nodesets)
        {
          const Node * node = mesh.node_ptr(node_id);

          if (node->processor_id() != mesh.comm().rank())
            continue;

          binfo.boundary_ids(node, node_boundary_ids);

          bool has_node_boundary_id = false;
          boundary_id_type matching_boundary_id = BoundaryInfo::invalid_id;
          for (boundary_id_type node_boundary_id : node_boundary_ids)
            if(parametrized_function_boundary_ids.count(node_boundary_id))
              {
                has_node_boundary_id = true;
                matching_boundary_id = node_boundary_id;
                break;
              }

          if(has_node_boundary_id)
            {
              _local_node_locations[node_id] = *node;
              _local_node_boundary_ids[node_id] = matching_boundary_id;
            }
        }
    }
  else
    {
      _local_quad_point_locations.clear();
      _local_quad_point_subdomain_ids.clear();
      _local_quad_point_JxW.clear();

      _local_quad_point_locations_perturbations.clear();

      for (const auto & elem : mesh.active_local_element_ptr_range())
        {
          auto elem_fe = context.get_element_fe(/*var=*/0, elem->dim());
          const std::vector<Real> & JxW = elem_fe->get_JxW();
          const std::vector<Point> & xyz = elem_fe->get_xyz();

          dof_id_type elem_id = elem->id();

          context.pre_fe_reinit(*this, elem);
          context.elem_fe_reinit();

          _local_quad_point_locations[elem_id] = xyz;
          _local_quad_point_JxW[elem_id] = JxW;
          _local_quad_point_subdomain_ids[elem_id] = elem->subdomain_id();

          if (get_rb_eim_evaluation().get_parametrized_function().requires_xyz_perturbations)
            {
              Real fd_delta = get_rb_eim_evaluation().get_parametrized_function().fd_delta;

              std::vector<std::vector<Point>> xyz_perturb_vec_at_qps;

              for (const Point & xyz_qp : xyz)
                {
                  std::vector<Point> xyz_perturb_vec;
                  if (elem->dim() == 3)
                    {
                      Point xyz_perturb = xyz_qp;

                      xyz_perturb(0) += fd_delta;
                      xyz_perturb_vec.emplace_back(xyz_perturb);
                      xyz_perturb(0) -= fd_delta;

                      xyz_perturb(1) += fd_delta;
                      xyz_perturb_vec.emplace_back(xyz_perturb);
                      xyz_perturb(1) -= fd_delta;

                      xyz_perturb(2) += fd_delta;
                      xyz_perturb_vec.emplace_back(xyz_perturb);
                      xyz_perturb(2) -= fd_delta;
                    }
                  else if (elem->dim() == 2)
                    {
                      // In this case we assume that we have a 2D element
                      // embedded in 3D space. In this case we have to use
                      // the following approach to perturb xyz:
                      //  1) inverse map xyz to the reference element
                      //  2) perturb on the reference element in the (xi,eta) "directions"
                      //  3) map the perturbed points back to the physical element
                      // This approach is necessary to ensure that the perturbed points
                      // are still in the element.

                      Point xi_eta =
                        FEMap::inverse_map(elem->dim(),
                                          elem,
                                          xyz_qp,
                                          /*Newton iteration tolerance*/ TOLERANCE,
                                          /*secure*/ true);

                      // Inverse map should map back to a 2D reference domain
                      libmesh_assert(std::abs(xi_eta(2)) < TOLERANCE);

                      Point xi_eta_perturb = xi_eta;

                      xi_eta_perturb(0) += fd_delta;
                      Point xyz_perturb_0 =
                        FEMap::map(elem->dim(),
                                  elem,
                                  xi_eta_perturb);
                      xi_eta_perturb(0) -= fd_delta;

                      xi_eta_perturb(1) += fd_delta;
                      Point xyz_perturb_1 =
                        FEMap::map(elem->dim(),
                                  elem,
                                  xi_eta_perturb);
                      xi_eta_perturb(1) -= fd_delta;

                      // Finally, we rescale xyz_perturb_0 and xyz_perturb_1 so that
                      // (xyz_perturb - xyz_qp).norm() == fd_delta, since this is
                      // required in order to compute finite differences correctly.
                      Point unit_0 = (xyz_perturb_0-xyz_qp).unit();
                      Point unit_1 = (xyz_perturb_1-xyz_qp).unit();

                      xyz_perturb_vec.emplace_back(xyz_qp + fd_delta*unit_0);
                      xyz_perturb_vec.emplace_back(xyz_qp + fd_delta*unit_1);
                    }
                  else
                    {
                      // We current do nothing in the dim=1 case since
                      // we have no need for this capability so far.
                      // Support for this case could be added if it is
                      // needed.
                    }

                  xyz_perturb_vec_at_qps.emplace_back(xyz_perturb_vec);
                }

              _local_quad_point_locations_perturbations[elem_id] = xyz_perturb_vec_at_qps;
            }
        }
    }
}

Number
RBEIMConstruction::inner_product(const QpDataMap & v, const QpDataMap & w, bool apply_comp_scaling)
{
  LOG_SCOPE("inner_product()", "RBEIMConstruction");

  Number val = 0.;

  for (const auto & [elem_id, v_comp_and_qp] : v)
    {
      const auto & w_comp_and_qp = libmesh_map_find(w, elem_id);
      const auto & JxW = libmesh_map_find(_local_quad_point_JxW, elem_id);

      for (const auto & comp : index_range(v_comp_and_qp))
        {
          const std::vector<Number> & v_qp = v_comp_and_qp[comp];
          const std::vector<Number> & w_qp = w_comp_and_qp[comp];

          Real comp_scaling = 1.;
          if (apply_comp_scaling)
            {
              // We square the component scaling here because it occurs twice in
              // the inner product calculation below.
              comp_scaling = std::pow(_component_scaling_in_training_set[comp], 2.);
            }

          for (unsigned int qp : index_range(JxW))
            val += JxW[qp] * comp_scaling * v_qp[qp] * libmesh_conj(w_qp[qp]);
        }
    }

  comm().sum(val);
  return val;
}

Number
RBEIMConstruction::side_inner_product(const SideQpDataMap & v, const SideQpDataMap & w, bool apply_comp_scaling)
{
  LOG_SCOPE("side_inner_product()", "RBEIMConstruction");

  Number val = 0.;

  for (const auto & [elem_and_side, v_comp_and_qp] : v)
    {
      const auto & w_comp_and_qp = libmesh_map_find(w, elem_and_side);
      const auto & JxW = libmesh_map_find(_local_side_quad_point_JxW, elem_and_side);

      for (const auto & comp : index_range(v_comp_and_qp))
        {
          const std::vector<Number> & v_qp = v_comp_and_qp[comp];
          const std::vector<Number> & w_qp = w_comp_and_qp[comp];

          Real comp_scaling = 1.;
          if (apply_comp_scaling)
            {
              // We square the component scaling here because it occurs twice in
              // the inner product calculation below.
              comp_scaling = std::pow(_component_scaling_in_training_set[comp], 2.);
            }

          for (unsigned int qp : index_range(JxW))
            val += JxW[qp] * comp_scaling * v_qp[qp] * libmesh_conj(w_qp[qp]);
        }
    }

  comm().sum(val);
  return val;
}

Number
RBEIMConstruction::node_inner_product(const NodeDataMap & v, const NodeDataMap & w, bool apply_comp_scaling)
{
  LOG_SCOPE("node_inner_product()", "RBEIMConstruction");

  Number val = 0.;

  for (const auto & [node_id, v_comps] : v)
    {
      const auto & w_comps = libmesh_map_find(w, node_id);

      for (const auto & comp : index_range(v_comps))
        {
          // There is no quadrature rule on nodes, so we just multiply the values directly.
          // Hence we effectively work with the Euclidean inner product in this case.

          Real comp_scaling = 1.;
          if (apply_comp_scaling)
            {
              // We square the component scaling here because it occurs twice in
              // the inner product calculation below.
              comp_scaling = std::pow(_component_scaling_in_training_set[comp], 2.);
            }

          val += comp_scaling * v_comps[comp] * libmesh_conj(w_comps[comp]);
        }
    }

  comm().sum(val);
  return val;
}

Real RBEIMConstruction::get_node_max_abs_value(const NodeDataMap & v) const
{
  Real max_value = 0.;

  for (const auto & pr : v)
    {
      const auto & values = pr.second;
      for (const auto & comp : index_range(values))
        {
          const auto & value = values[comp];

          Real comp_scaling = 1.;
          if (get_rb_eim_evaluation().scale_components_in_enrichment().count(comp))
            {
              // Make sure that _component_scaling_in_training_set is initialized
              libmesh_error_msg_if(comp >= _component_scaling_in_training_set.size(),
                                  "Invalid vector index");
              comp_scaling = _component_scaling_in_training_set[comp];
            }

          max_value = std::max(max_value, std::abs(value * comp_scaling));
        }
    }

  comm().max(max_value);
  return max_value;
}

void RBEIMConstruction::enrich_eim_approximation(unsigned int training_index,
                                                 bool add_basis_function,
                                                 EimPointData * eim_point_data)
{
  LOG_SCOPE("enrich_eim_approximation()", "RBEIMConstruction");

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  set_params_from_training_set(training_index);

  if (eim_eval.get_parametrized_function().on_mesh_sides())
    enrich_eim_approximation_on_sides(_local_side_parametrized_functions_for_training[training_index],
                                      add_basis_function,
                                      eim_point_data);
  else if (eim_eval.get_parametrized_function().on_mesh_nodes())
    enrich_eim_approximation_on_nodes(_local_node_parametrized_functions_for_training[training_index],
                                      add_basis_function,
                                      eim_point_data);
  else
    {
      enrich_eim_approximation_on_interiors(_local_parametrized_functions_for_training[training_index],
                                            add_basis_function,
                                            eim_point_data);
    }
}

bool RBEIMConstruction::enrich_eim_approximation_on_sides(const SideQpDataMap & side_pf,
                                                          bool add_basis_function,
                                                          EimPointData * eim_point_data)
{
  // Make a copy of the input parametrized function, since we will modify this below
  // to give us a new basis function.
  SideQpDataMap local_pf = side_pf;

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  // If we have at least one basis function, then we need to use
  // rb_eim_solve() to find the EIM interpolation error. Otherwise,
  // just use solution as is.
  if (!eim_point_data && (eim_eval.get_n_basis_functions() > 0))
    {
      // Get the right-hand side vector for the EIM approximation
      // by sampling the parametrized function (stored in solution)
      // at the interpolation points.
      unsigned int RB_size = eim_eval.get_n_basis_functions();
      DenseVector<Number> EIM_rhs(RB_size);
      for (unsigned int i=0; i<RB_size; i++)
        {
          EIM_rhs(i) =
            RBEIMEvaluation::get_parametrized_function_side_value(comm(),
                                                                  local_pf,
                                                                  eim_eval.get_interpolation_points_elem_id(i),
                                                                  eim_eval.get_interpolation_points_side_index(i),
                                                                  eim_eval.get_interpolation_points_comp(i),
                                                                  eim_eval.get_interpolation_points_qp(i));
        }

      eim_eval.set_parameters( get_parameters() );
      DenseVector<Number> rb_eim_solution = eim_eval.rb_eim_solve(EIM_rhs);

      // Load the "EIM residual" into solution by subtracting
      // the EIM approximation
      eim_eval.side_decrement_vector(local_pf, rb_eim_solution);
    }

  // Find the quadrature point at which local_pf (which now stores
  // the "EIM residual") has maximum absolute value
  Number optimal_value = 0.;
  Point optimal_point;
  unsigned int optimal_comp = 0;
  dof_id_type optimal_elem_id = DofObject::invalid_id;
  unsigned int optimal_side_index = 0;
  subdomain_id_type optimal_subdomain_id = 0;
  boundary_id_type optimal_boundary_id = 0;
  unsigned int optimal_qp = 0;
  std::vector<Point> optimal_point_perturbs;
  std::vector<Real> optimal_point_phi_i_qp;

  // Initialize largest_abs_value to be negative so that it definitely gets updated.
  Real largest_abs_value = -1.;

  // In order to compute phi_i_qp, we initialize a FEMContext
  FEMContext con(*this);
  init_context(con);

  for (const auto & [elem_and_side, comp_and_qp] : local_pf)
    {
      dof_id_type elem_id = elem_and_side.first;
      unsigned int side_index = elem_and_side.second;

      const Elem & elem_ref = get_mesh().elem_ref(elem_id);
      con.pre_fe_reinit(*this, &elem_ref);

      unsigned int side_type = libmesh_map_find(_local_side_quad_point_side_types, elem_and_side);

      std::vector<std::vector<Real>> phi;
      // side_type == 0 --> standard side
      // side_type == 1 --> shellface
      if (side_type == 0)
        {
          // TODO: We only want the "dofs on side" entries
          // from phi_side. Could do this by initing an FE object
          // on the side itself, rather than using get_side_fe().
          auto side_fe = con.get_side_fe(/*var=*/ 0);
          side_fe->reinit(&elem_ref, side_index);

          phi = side_fe->get_phi();
        }
      else if (side_type == 1)
        {
          con.elem_fe_reinit();

          auto elem_fe = con.get_element_fe(/*var=*/0, elem_ref.dim());
          phi = elem_fe->get_phi();
        }
      else
        libmesh_error_msg ("Unrecognized side_type: " << side_type);

      for (const auto & comp : index_range(comp_and_qp))
        {
          const std::vector<Number> & qp_values = comp_and_qp[comp];

          for (auto qp : index_range(qp_values))
            {
              Number value = qp_values[qp];
              Real abs_value = std::abs(value);

              bool update_optimal_point = false;
              if (!eim_point_data)
                update_optimal_point = (abs_value > largest_abs_value);
              else
                update_optimal_point = (elem_id == eim_point_data->elem_id) &&
                                       (side_index == eim_point_data->side_index) &&
                                       (comp == eim_point_data->comp_index) &&
                                       (qp == eim_point_data->qp_index);

              if (update_optimal_point)
                {
                  largest_abs_value = abs_value;
                  optimal_value = value;
                  optimal_comp = comp;
                  optimal_elem_id = elem_id;
                  optimal_side_index = side_index;
                  optimal_qp = qp;

                  optimal_point_phi_i_qp.resize(phi.size());
                  for (auto i : index_range(phi))
                    optimal_point_phi_i_qp[i] = phi[i][qp];

                  const auto & point_list =
                    libmesh_map_find(_local_side_quad_point_locations, elem_and_side);

                  libmesh_error_msg_if(qp >= point_list.size(), "Error: Invalid qp");

                  optimal_point = point_list[qp];

                  optimal_subdomain_id = libmesh_map_find(_local_side_quad_point_subdomain_ids, elem_and_side);
                  optimal_boundary_id = libmesh_map_find(_local_side_quad_point_boundary_ids, elem_and_side);

                  if (get_rb_eim_evaluation().get_parametrized_function().requires_xyz_perturbations)
                    {
                      const auto & perturb_list =
                        libmesh_map_find(_local_side_quad_point_locations_perturbations, elem_and_side);

                      libmesh_error_msg_if(qp >= perturb_list.size(), "Error: Invalid qp");

                      optimal_point_perturbs = perturb_list[qp];
                    }
                }
            }
        }
    }

  // Find out which processor has the largest of the abs values
  // and broadcast from that processor.
  unsigned int proc_ID_index;
  this->comm().maxloc(largest_abs_value, proc_ID_index);

  this->comm().broadcast(optimal_value, proc_ID_index);
  this->comm().broadcast(optimal_point, proc_ID_index);
  this->comm().broadcast(optimal_comp, proc_ID_index);
  this->comm().broadcast(optimal_elem_id, proc_ID_index);
  this->comm().broadcast(optimal_side_index, proc_ID_index);
  this->comm().broadcast(optimal_subdomain_id, proc_ID_index);
  this->comm().broadcast(optimal_boundary_id, proc_ID_index);
  this->comm().broadcast(optimal_qp, proc_ID_index);
  this->comm().broadcast(optimal_point_perturbs, proc_ID_index);
  this->comm().broadcast(optimal_point_phi_i_qp, proc_ID_index);

  libmesh_error_msg_if(optimal_elem_id == DofObject::invalid_id, "Error: Invalid element ID");

  if (add_basis_function)
    {
      if (optimal_value == 0.)
        {
          libMesh::out << "Encountered linearly dependent data in EIM enrichment, hence skip adding new basis function" << std::endl;
          return true;
        }

      // Scale local_pf so that its largest value is 1.0
      scale_parametrized_function(local_pf, 1./optimal_value);

      // Add local_pf as the new basis function and store data
      // associated with the interpolation point.
      eim_eval.add_side_basis_function(local_pf);
    }

  eim_eval.add_side_interpolation_data(optimal_point,
                                       optimal_comp,
                                       optimal_elem_id,
                                       optimal_side_index,
                                       optimal_subdomain_id,
                                       optimal_boundary_id,
                                       optimal_qp,
                                       optimal_point_perturbs,
                                       optimal_point_phi_i_qp);

  // In this case we did not encounter a linearly dependent basis function, so return false
  return false;
}

bool RBEIMConstruction::enrich_eim_approximation_on_nodes(const NodeDataMap & node_pf,
                                                          bool add_basis_function,
                                                          EimPointData * eim_point_data)
{
  // Make a copy of the input parametrized function, since we will modify this below
  // to give us a new basis function.
  NodeDataMap local_pf = node_pf;

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  // If we have at least one basis function, then we need to use
  // rb_eim_solve() to find the EIM interpolation error. Otherwise,
  // just use solution as is.
  if (!eim_point_data && (eim_eval.get_n_basis_functions() > 0))
    {
      // Get the right-hand side vector for the EIM approximation
      // by sampling the parametrized function (stored in solution)
      // at the interpolation points.
      unsigned int RB_size = eim_eval.get_n_basis_functions();
      DenseVector<Number> EIM_rhs(RB_size);
      for (unsigned int i=0; i<RB_size; i++)
        {
          EIM_rhs(i) =
            RBEIMEvaluation::get_parametrized_function_node_value(comm(),
                                                                  local_pf,
                                                                  eim_eval.get_interpolation_points_node_id(i),
                                                                  eim_eval.get_interpolation_points_comp(i));
        }

      eim_eval.set_parameters( get_parameters() );
      DenseVector<Number> rb_eim_solution = eim_eval.rb_eim_solve(EIM_rhs);

      // Load the "EIM residual" into solution by subtracting
      // the EIM approximation
      eim_eval.node_decrement_vector(local_pf, rb_eim_solution);
    }

  // Find the quadrature point at which local_pf (which now stores
  // the "EIM residual") has maximum absolute value
  Number optimal_value = 0.;
  Point optimal_point;
  unsigned int optimal_comp = 0;
  dof_id_type optimal_node_id = DofObject::invalid_id;
  boundary_id_type optimal_boundary_id = 0;

  // Initialize largest_abs_value to be negative so that it definitely gets updated.
  Real largest_abs_value = -1.;

  for (const auto & [node_id, values] : local_pf)
    {
      for (unsigned int comp : index_range(values))
        {
          Number value = values[comp];
          Real abs_value = std::abs(value);

          bool update_optimal_point = false;
          if (!eim_point_data)
            update_optimal_point = (abs_value > largest_abs_value);
          else
            update_optimal_point = (node_id == eim_point_data->node_id) &&
                                   (comp == eim_point_data->comp_index);

          if (update_optimal_point)
            {
              largest_abs_value = abs_value;
              optimal_value = value;
              optimal_comp = comp;
              optimal_node_id = node_id;

              optimal_point = libmesh_map_find(_local_node_locations, node_id);

              optimal_boundary_id = libmesh_map_find(_local_node_boundary_ids, node_id);
            }
        }
    }

  // Find out which processor has the largest of the abs values
  // and broadcast from that processor.
  unsigned int proc_ID_index;
  this->comm().maxloc(largest_abs_value, proc_ID_index);

  this->comm().broadcast(optimal_value, proc_ID_index);
  this->comm().broadcast(optimal_point, proc_ID_index);
  this->comm().broadcast(optimal_comp, proc_ID_index);
  this->comm().broadcast(optimal_node_id, proc_ID_index);
  this->comm().broadcast(optimal_boundary_id, proc_ID_index);

  libmesh_error_msg_if(optimal_node_id == DofObject::invalid_id, "Error: Invalid node ID");

  if (add_basis_function)
    {
      if (optimal_value == 0.)
        {
          libMesh::out << "Encountered linearly dependent data in EIM enrichment, hence skip adding new basis function" << std::endl;
          return true;
        }

      // Scale local_pf so that its largest value is 1.0
      scale_node_parametrized_function(local_pf, 1./optimal_value);

      // Add local_pf as the new basis function and store data
      // associated with the interpolation point.
      eim_eval.add_node_basis_function(local_pf);
    }

  eim_eval.add_node_interpolation_data(optimal_point,
                                       optimal_comp,
                                       optimal_node_id,
                                       optimal_boundary_id);

  // In this case we did not encounter a linearly dependent basis function, so return false
  return false;
}

bool RBEIMConstruction::enrich_eim_approximation_on_interiors(const QpDataMap & interior_pf,
                                                              bool add_basis_function,
                                                              EimPointData * eim_point_data)
{
  // Make a copy of the input parametrized function, since we will modify this below
  // to give us a new basis function.
  QpDataMap local_pf = interior_pf;

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  // If we have at least one basis function, then we need to use
  // rb_eim_solve() to find the EIM interpolation error. Otherwise,
  // just use solution as is.
  if (!eim_point_data && (eim_eval.get_n_basis_functions() > 0))
    {
      // Get the right-hand side vector for the EIM approximation
      // by sampling the parametrized function (stored in solution)
      // at the interpolation points.
      unsigned int RB_size = eim_eval.get_n_basis_functions();
      DenseVector<Number> EIM_rhs(RB_size);
      for (unsigned int i=0; i<RB_size; i++)
        {
          EIM_rhs(i) =
            RBEIMEvaluation::get_parametrized_function_value(comm(),
                                                            local_pf,
                                                            eim_eval.get_interpolation_points_elem_id(i),
                                                            eim_eval.get_interpolation_points_comp(i),
                                                            eim_eval.get_interpolation_points_qp(i));
        }

      eim_eval.set_parameters( get_parameters() );
      DenseVector<Number> rb_eim_solution = eim_eval.rb_eim_solve(EIM_rhs);

      // Load the "EIM residual" into solution by subtracting
      // the EIM approximation
      eim_eval.decrement_vector(local_pf, rb_eim_solution);
    }

  // Find the quadrature point at which local_pf (which now stores
  // the "EIM residual") has maximum absolute value
  Number optimal_value = 0.;
  Point optimal_point;
  unsigned int optimal_comp = 0;
  dof_id_type optimal_elem_id = DofObject::invalid_id;
  subdomain_id_type optimal_subdomain_id = 0;
  unsigned int optimal_qp = 0;
  std::vector<Point> optimal_point_perturbs;
  std::vector<Real> optimal_point_phi_i_qp;
  ElemType optimal_elem_type = INVALID_ELEM;
  std::vector<Real> optimal_JxW_all_qp;
  std::vector<std::vector<Real>> optimal_phi_i_all_qp;
  Order optimal_qrule_order = INVALID_ORDER;
  Point optimal_dxyzdxi_elem_center;
  Point optimal_dxyzdeta_elem_center;

  // Initialize largest_abs_value to be negative so that it definitely gets updated.
  Real largest_abs_value = -1.;

  // In order to compute phi_i_qp, we initialize a FEMContext
  FEMContext con(*this);
  for (auto dim : con.elem_dimensions())
    {
      auto fe = con.get_element_fe(/*var=*/0, dim);
      fe->get_phi();
      fe->get_JxW();
      fe->get_dxyzdxi();
      fe->get_dxyzdeta();
    }

  for (const auto & [elem_id, comp_and_qp] : local_pf)
    {
      // Also initialize phi in order to compute phi_i_qp
      const Elem & elem_ref = get_mesh().elem_ref(elem_id);
      con.pre_fe_reinit(*this, &elem_ref);

      auto elem_fe = con.get_element_fe(/*var=*/0, elem_ref.dim());
      const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();
      const auto & JxW = elem_fe->get_JxW();
      const auto & dxyzdxi = elem_fe->get_dxyzdxi();
      const auto & dxyzdeta = elem_fe->get_dxyzdeta();

      elem_fe->reinit(&elem_ref);

      for (const auto & comp : index_range(comp_and_qp))
        {
          const std::vector<Number> & qp_values = comp_and_qp[comp];

          for (auto qp : index_range(qp_values))
            {
              Number value = qp_values[qp];
              Real abs_value = std::abs(value);

              if (get_rb_eim_evaluation().scale_components_in_enrichment().count(comp))
                abs_value *= _component_scaling_in_training_set[comp];

              bool update_optimal_point = false;
              if (!eim_point_data)
                update_optimal_point = (abs_value > largest_abs_value);
              else
                update_optimal_point = (elem_id == eim_point_data->elem_id) &&
                                       (comp == eim_point_data->comp_index) &&
                                       (qp == eim_point_data->qp_index);

              if (update_optimal_point)
                {
                  largest_abs_value = abs_value;
                  optimal_value = value;
                  optimal_comp = comp;
                  optimal_elem_id = elem_id;
                  optimal_qp = qp;
                  optimal_elem_type = elem_ref.type();

                  optimal_point_phi_i_qp.resize(phi.size());
                  for (auto i : index_range(phi))
                    optimal_point_phi_i_qp[i] = phi[i][qp];

                  const auto & point_list =
                    libmesh_map_find(_local_quad_point_locations, elem_id);

                  libmesh_error_msg_if(qp >= point_list.size(), "Error: Invalid qp");

                  optimal_point = point_list[qp];

                  optimal_subdomain_id = libmesh_map_find(_local_quad_point_subdomain_ids, elem_id);

                  if (get_rb_eim_evaluation().get_parametrized_function().requires_xyz_perturbations)
                    {
                      const auto & perturb_list =
                        libmesh_map_find(_local_quad_point_locations_perturbations, elem_id);

                      libmesh_error_msg_if(qp >= perturb_list.size(), "Error: Invalid qp");

                      optimal_point_perturbs = perturb_list[qp];
                    }

                  if (get_rb_eim_evaluation().get_parametrized_function().requires_all_elem_qp_data)
                    {
                      optimal_JxW_all_qp = JxW;
                      optimal_phi_i_all_qp = phi;
                    }

                  if (get_rb_eim_evaluation().get_parametrized_function().requires_all_elem_center_data)
                    {
                      optimal_qrule_order = con.get_element_qrule().get_order();
                      // Get data derivatives at vertex average
                      std::vector<Point> nodes = { elem_ref.reference_elem()->vertex_average() };
                      elem_fe->reinit (&elem_ref, &nodes);

                      Point dxyzdxi_pt, dxyzdeta_pt;
                      if (con.get_elem_dim()>0)
                        dxyzdxi_pt = dxyzdxi[0];
                      if (con.get_elem_dim()>1)
                        dxyzdeta_pt = dxyzdeta[0];

                      optimal_dxyzdxi_elem_center = dxyzdxi_pt;
                      optimal_dxyzdeta_elem_center = dxyzdeta_pt;

                      elem_fe->reinit(&elem_ref);
                    }
                }
            }
        }
    }

  // Find out which processor has the largest of the abs values
  // and broadcast from that processor.
  unsigned int proc_ID_index;
  this->comm().maxloc(largest_abs_value, proc_ID_index);

  this->comm().broadcast(optimal_value, proc_ID_index);
  this->comm().broadcast(optimal_point, proc_ID_index);
  this->comm().broadcast(optimal_comp, proc_ID_index);
  this->comm().broadcast(optimal_elem_id, proc_ID_index);
  this->comm().broadcast(optimal_subdomain_id, proc_ID_index);
  this->comm().broadcast(optimal_qp, proc_ID_index);
  this->comm().broadcast(optimal_point_perturbs, proc_ID_index);
  this->comm().broadcast(optimal_point_phi_i_qp, proc_ID_index);
  this->comm().broadcast(optimal_JxW_all_qp, proc_ID_index);
  this->comm().broadcast(optimal_phi_i_all_qp, proc_ID_index);
  this->comm().broadcast(optimal_dxyzdxi_elem_center, proc_ID_index);
  this->comm().broadcast(optimal_dxyzdeta_elem_center, proc_ID_index);

  // Cast optimal_elem_type to an int in order to broadcast it
  {
    int optimal_elem_type_int = static_cast<int>(optimal_elem_type);
    this->comm().broadcast(optimal_elem_type_int, proc_ID_index);
    optimal_elem_type = static_cast<ElemType>(optimal_elem_type_int);
  }

  // Cast optimal_qrule_order to an int in order to broadcast it
  {
    int optimal_qrule_order_int = static_cast<int>(optimal_qrule_order);
    this->comm().broadcast(optimal_qrule_order_int, proc_ID_index);
    optimal_qrule_order = static_cast<Order>(optimal_qrule_order_int);
  }

  libmesh_error_msg_if(optimal_elem_id == DofObject::invalid_id, "Error: Invalid element ID");

  if (add_basis_function)
    {
      if (optimal_value == 0.)
        {
          libMesh::out << "Encountered linearly dependent data in EIM enrichment, hence skip adding new basis function" << std::endl;
          return true;
        }

      // Scale local_pf so that its largest value is 1.0
      scale_parametrized_function(local_pf, 1./optimal_value);

      // Add local_pf as the new basis function and store data
      // associated with the interpolation point.
      eim_eval.add_basis_function(local_pf);
    }

  eim_eval.add_interpolation_data(optimal_point,
                                  optimal_comp,
                                  optimal_elem_id,
                                  optimal_subdomain_id,
                                  optimal_qp,
                                  optimal_point_perturbs,
                                  optimal_point_phi_i_qp,
                                  optimal_elem_type,
                                  optimal_JxW_all_qp,
                                  optimal_phi_i_all_qp,
                                  optimal_qrule_order,
                                  optimal_dxyzdxi_elem_center,
                                  optimal_dxyzdeta_elem_center);

  // In this case we did not encounter a linearly dependent basis function, so return false
  return false;
}

void RBEIMConstruction::update_eim_matrices(bool set_eim_error_indicator)
{
  LOG_SCOPE("update_eim_matrices()", "RBEIMConstruction");

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();
  unsigned int RB_size = eim_eval.get_n_basis_functions();

  libmesh_assert_msg(RB_size >= 1, "Must have at least 1 basis function.");

  if (set_eim_error_indicator)
    {
      // Here we have RB_size EIM basis functions, and RB_size+1 interpolation points,
      // since we should have added one extra interpolation point for the EIM error
      // indicator. As a result, we use RB_size as the index to access the (RB_size+1)^th
      // interpolation point in the calls to eim_eval.get_interpolation_points_*.
      DenseVector<Number> extra_point_row(RB_size);

      if (eim_eval.get_parametrized_function().on_mesh_sides())
        {
          // update the EIM interpolation matrix
          for (unsigned int j=0; j<RB_size; j++)
            {
              // Evaluate the basis functions at the new interpolation point in order
              // to update the interpolation matrix
              Number value =
                eim_eval.get_eim_basis_function_side_value(j,
                                                           eim_eval.get_interpolation_points_elem_id(RB_size),
                                                           eim_eval.get_interpolation_points_side_index(RB_size),
                                                           eim_eval.get_interpolation_points_comp(RB_size),
                                                           eim_eval.get_interpolation_points_qp(RB_size));
              extra_point_row(j) = value;
            }
        }
      else if (eim_eval.get_parametrized_function().on_mesh_nodes())
        {
          // update the EIM interpolation matrix
          for (unsigned int j=0; j<RB_size; j++)
            {
              // Evaluate the basis functions at the new interpolation point in order
              // to update the interpolation matrix
              Number value =
                eim_eval.get_eim_basis_function_node_value(j,
                                                           eim_eval.get_interpolation_points_node_id(RB_size),
                                                           eim_eval.get_interpolation_points_comp(RB_size));
              extra_point_row(j) = value;
            }
        }
      else
        {
          // update the EIM interpolation matrix
          for (unsigned int j=0; j<RB_size; j++)
            {
              // Evaluate the basis functions at the new interpolation point in order
              // to update the interpolation matrix
              Number value =
                eim_eval.get_eim_basis_function_value(j,
                                                      eim_eval.get_interpolation_points_elem_id(RB_size),
                                                      eim_eval.get_interpolation_points_comp(RB_size),
                                                      eim_eval.get_interpolation_points_qp(RB_size));
              extra_point_row(j) = value;
            }
        }

      eim_eval.set_error_indicator_interpolation_row(extra_point_row);
      return;
    }

  if (eim_eval.get_parametrized_function().on_mesh_sides())
    {
      // update the matrix that is used to evaluate L2 projections
      // into the EIM approximation space
      for (unsigned int i=(RB_size-1); i<RB_size; i++)
        {
          for (unsigned int j=0; j<RB_size; j++)
            {
              Number value = side_inner_product(eim_eval.get_side_basis_function(j),
                                                eim_eval.get_side_basis_function(i),
                                                /*apply_comp_scaling*/ false);

              _eim_projection_matrix(i,j) = value;
              if (i!=j)
                {
                  // The inner product matrix is assumed to be hermitian
                  _eim_projection_matrix(j,i) = libmesh_conj(value);
                }
            }
        }

      // update the EIM interpolation matrix
      for (unsigned int j=0; j<RB_size; j++)
        {
          // Evaluate the basis functions at the new interpolation point in order
          // to update the interpolation matrix
          Number value =
            eim_eval.get_eim_basis_function_side_value(j,
                                                       eim_eval.get_interpolation_points_elem_id(RB_size-1),
                                                       eim_eval.get_interpolation_points_side_index(RB_size-1),
                                                       eim_eval.get_interpolation_points_comp(RB_size-1),
                                                       eim_eval.get_interpolation_points_qp(RB_size-1));
          eim_eval.set_interpolation_matrix_entry(RB_size-1, j, value);
        }
    }
  else if (eim_eval.get_parametrized_function().on_mesh_nodes())
    {
      // update the matrix that is used to evaluate L2 projections
      // into the EIM approximation space
      for (unsigned int i=(RB_size-1); i<RB_size; i++)
        {
          for (unsigned int j=0; j<RB_size; j++)
            {
              Number value = node_inner_product(eim_eval.get_node_basis_function(j),
                                                eim_eval.get_node_basis_function(i),
                                                /*apply_comp_scaling*/ false);

              _eim_projection_matrix(i,j) = value;
              if (i!=j)
                {
                  // The inner product matrix is assumed to be hermitian
                  _eim_projection_matrix(j,i) = libmesh_conj(value);
                }
            }
        }

      // update the EIM interpolation matrix
      for (unsigned int j=0; j<RB_size; j++)
        {
          // Evaluate the basis functions at the new interpolation point in order
          // to update the interpolation matrix
          Number value =
            eim_eval.get_eim_basis_function_node_value(j,
                                                       eim_eval.get_interpolation_points_node_id(RB_size-1),
                                                       eim_eval.get_interpolation_points_comp(RB_size-1));
          eim_eval.set_interpolation_matrix_entry(RB_size-1, j, value);
        }
    }
  else
    {
      // update the matrix that is used to evaluate L2 projections
      // into the EIM approximation space
      for (unsigned int i=(RB_size-1); i<RB_size; i++)
        {
          for (unsigned int j=0; j<RB_size; j++)
            {
              Number value = inner_product(eim_eval.get_basis_function(j),
                                           eim_eval.get_basis_function(i),
                                           /*apply_comp_scaling*/ false);

              _eim_projection_matrix(i,j) = value;
              if (i!=j)
                {
                  // The inner product matrix is assumed to be hermitian
                  _eim_projection_matrix(j,i) = libmesh_conj(value);
                }
            }
        }

      // update the EIM interpolation matrix
      for (unsigned int j=0; j<RB_size; j++)
        {
          // Evaluate the basis functions at the new interpolation point in order
          // to update the interpolation matrix
          Number value =
            eim_eval.get_eim_basis_function_value(j,
                                                  eim_eval.get_interpolation_points_elem_id(RB_size-1),
                                                  eim_eval.get_interpolation_points_comp(RB_size-1),
                                                  eim_eval.get_interpolation_points_qp(RB_size-1));
          eim_eval.set_interpolation_matrix_entry(RB_size-1, j, value);
        }
    }
}

void RBEIMConstruction::scale_node_parametrized_function(NodeDataMap & local_pf,
                                                         Number scaling_factor)
{
  for (auto & pr : local_pf)
    {
      auto & values = pr.second;
      for ( auto & value : values)
        value *= scaling_factor;
    }
}

unsigned int RBEIMConstruction::get_random_int_0_to_n(unsigned int n)
{
  // std::random_device seed;
  // std::mt19937 gen{seed()};
  // We do not use a random seed here, since we generally prefer our results
  // to reproducible, rather than fully random. If desired we could provide an
  // option to use the random seed approach (commented out above).
  std::default_random_engine gen;
  std::uniform_int_distribution<> dist{0, static_cast<int>(n)};
  return dist(gen);
}

EimPointData RBEIMConstruction::get_random_point(const QpDataMap & v)
{
  EimPointData eim_point_data;

  // If we have more than one process, then we need to do a parallel union
  // of v to make sure that we have data from all processors. Our approach
  // here is to set v_ptr to either v or global_v, depending on whether we
  // are in parallel or serial. The purpose of this approach is to avoid
  // making a copy of v in the case that this is a serial job.
  QpDataMap const * v_ptr = nullptr;
  QpDataMap global_v;
  if (comm().size() > 1)
  {
    global_v = v;

    // We only use global_v on proc 0, so we set the second argument of
    // set_union() to zero here to indicate that we only need the result
    // on proc 0.
    comm().set_union(global_v, 0);
    v_ptr = &global_v;
  }
  else
  {
    v_ptr = &v;
  }

  bool error_finding_new_element = false;
  if (comm().rank() == 0)
    {
      const VectorizedEvalInput & vec_eval_input = get_rb_eim_evaluation().get_vec_eval_input();

      {
        std::set<dof_id_type> previous_elem_ids(vec_eval_input.elem_ids.begin(), vec_eval_input.elem_ids.end());

        // We ensure that we select a point that has not been selected previously
        // by setting up new_elem_ids to contain only elements that are not in
        // previous_elem_ids, and then selecting the elem_id at random from new_elem_ids.
        // We give an error if there are no elements in new_elem_ids. This is potentially
        // an overzealous assertion since we could pick an element that has already
        // been selected as long as we pick a (comp_index, qp_index) that has not already
        // been selected for that element.
        //
        // However, in general we do not expect all elements to be selected in the EIM
        // training, so it is reasonable to use the simple assertion below. Moreover, by
        // ensuring that we choose a new element we should typically ensure that the
        // randomly selected point has some separation from the previous EIM points, which
        // is typically desirable if we want EIM evaluations that are independent from
        // the EIM points (e.g. for EIM error indicator purposes).
        std::set<dof_id_type> new_elem_ids;
        for (const auto & v_pair : *v_ptr)
          if (previous_elem_ids.count(v_pair.first) == 0)
            new_elem_ids.insert(v_pair.first);

        // If new_elem_ids is empty then we set error_finding_new_element to true.
        // We then broadcast the value of error_finding_new_element to all processors
        // below in order to ensure that all processors agree on whether or not
        // there was an error.
        error_finding_new_element = (new_elem_ids.empty());

        if (!error_finding_new_element)
          {
            unsigned int random_elem_idx = get_random_int_0_to_n(new_elem_ids.size()-1);

            auto item = new_elem_ids.begin();
            std::advance(item, random_elem_idx);
            eim_point_data.elem_id = *item;
          }
      }

      if (!error_finding_new_element)
        {
          {
            const auto & vars_and_qps = libmesh_map_find(*v_ptr,eim_point_data.elem_id);
            eim_point_data.comp_index = get_random_int_0_to_n(vars_and_qps.size()-1);
          }

          {
            const auto & qps = libmesh_map_find(*v_ptr,eim_point_data.elem_id)[eim_point_data.comp_index];
            eim_point_data.qp_index = get_random_int_0_to_n(qps.size()-1);
          }
        }
    }

  comm().broadcast(error_finding_new_element);
  libmesh_error_msg_if(error_finding_new_element, "Could not find new element in get_random_point()");

  // Broadcast the values computed above from rank 0
  comm().broadcast(eim_point_data.elem_id);
  comm().broadcast(eim_point_data.comp_index);
  comm().broadcast(eim_point_data.qp_index);

  return eim_point_data;
}

EimPointData RBEIMConstruction::get_random_point(const SideQpDataMap & v)
{
  EimPointData eim_point_data;

  // If we have more than one process, then we need to do a parallel union
  // of v to make sure that we have data from all processors. Our approach
  // here is to set v_ptr to either v or global_v, depending on whether we
  // are in parallel or serial. The purpose of this approach is to avoid
  // making a copy of v in the case that this is a serial job.
  SideQpDataMap const * v_ptr = nullptr;
  SideQpDataMap global_v;
  if (comm().size() > 1)
  {
    global_v = v;

    // We only use global_v on proc 0, so we set the second argument of
    // set_union() to zero here to indicate that we only need the result
    // on proc 0.
    comm().set_union(global_v, 0);
    v_ptr = &global_v;
  }
  else
  {
    v_ptr = &v;
  }

  bool error_finding_new_element_and_side = false;
  if (comm().rank() == 0)
    {
      const VectorizedEvalInput & vec_eval_input = get_rb_eim_evaluation().get_vec_eval_input();

      std::pair<dof_id_type,unsigned int> elem_and_side;
      {
        std::set<std::pair<dof_id_type,unsigned int>> previous_elem_and_side_ids;
        for (const auto idx : index_range(vec_eval_input.elem_ids))
          {
            previous_elem_and_side_ids.insert(
              std::make_pair(vec_eval_input.elem_ids[idx],
                             vec_eval_input.side_indices[idx]));
          }

        // See discussion above in the QpDataMap case for the justification
        // of how we set up new_elem_and_side_ids below.
        std::set<std::pair<dof_id_type,unsigned int>> new_elem_and_side_ids;
        for (const auto & v_pair : *v_ptr)
          if (previous_elem_and_side_ids.count(v_pair.first) == 0)
            new_elem_and_side_ids.insert(v_pair.first);

        // If new_elem_and_side_ids is empty then we set error_finding_new_element_and_side
        // to true. We then broadcast the value of error_finding_new_element_and_side to all
        // processors below in order to ensure that all processors agree on whether
        // or not there was an error.
        error_finding_new_element_and_side = (new_elem_and_side_ids.empty());

        if (!error_finding_new_element_and_side)
          {
            unsigned int random_elem_and_side_idx = get_random_int_0_to_n(new_elem_and_side_ids.size()-1);

            auto item = new_elem_and_side_ids.begin();
            std::advance(item, random_elem_and_side_idx);
            elem_and_side = *item;
            eim_point_data.elem_id = elem_and_side.first;
            eim_point_data.side_index = elem_and_side.second;
          }
      }

      if (!error_finding_new_element_and_side)
        {
          {
            const auto & vars_and_qps = libmesh_map_find(*v_ptr,elem_and_side);
            eim_point_data.comp_index = get_random_int_0_to_n(vars_and_qps.size()-1);
          }

          {
            const auto & qps = libmesh_map_find(*v_ptr,elem_and_side)[eim_point_data.comp_index];
            eim_point_data.qp_index = get_random_int_0_to_n(qps.size()-1);
          }
        }
    }

  comm().broadcast(error_finding_new_element_and_side);
  libmesh_error_msg_if(error_finding_new_element_and_side, "Could not find new (element,side) in get_random_point()");

  // Broadcast the values computed above from rank 0
  comm().broadcast(eim_point_data.elem_id);
  comm().broadcast(eim_point_data.side_index);
  comm().broadcast(eim_point_data.comp_index);
  comm().broadcast(eim_point_data.qp_index);

  return eim_point_data;
}

EimPointData RBEIMConstruction::get_random_point(const NodeDataMap & v)
{
  EimPointData eim_point_data;

  // If we have more than one process, then we need to do a parallel union
  // of v to make sure that we have data from all processors. Our approach
  // here is to set v_ptr to either v or global_v, depending on whether we
  // are in parallel or serial. The purpose of this approach is to avoid
  // making a copy of v in the case that this is a serial job.
  NodeDataMap const * v_ptr = nullptr;
  NodeDataMap global_v;
  if (comm().size() > 1)
  {
    global_v = v;

    // We only use global_v on proc 0, so we set the second argument of
    // set_union() to zero here to indicate that we only need the result
    // on proc 0.
    comm().set_union(global_v, 0);
    v_ptr = &global_v;
  }
  else
  {
    v_ptr = &v;
  }

  bool error_finding_new_node = false;
  if (comm().rank() == 0)
    {
      const VectorizedEvalInput & vec_eval_input = get_rb_eim_evaluation().get_vec_eval_input();

      {
        std::set<dof_id_type> previous_node_ids(vec_eval_input.node_ids.begin(), vec_eval_input.node_ids.end());

        // See discussion above in the QpDataMap case for the justification
        // of how we set up new_node_ids below.
        std::set<dof_id_type> new_node_ids;
        for (const auto & v_pair : *v_ptr)
          if (previous_node_ids.count(v_pair.first) == 0)
            new_node_ids.insert(v_pair.first);

        // If new_node_ids is empty then we set error_finding_new_node
        // to true. We then broadcast the value of error_finding_new_node to all
        // processors below in order to ensure that all processors agree on whether
        // or not there was an error.
        error_finding_new_node = (new_node_ids.empty());

        if (!error_finding_new_node)
          {
            unsigned int random_node_idx = get_random_int_0_to_n(new_node_ids.size()-1);

            auto item = new_node_ids.begin();
            std::advance(item, random_node_idx);
            eim_point_data.node_id = *item;
          }
      }

      if (!error_finding_new_node)
        {
          const auto & vars = libmesh_map_find(*v_ptr,eim_point_data.node_id);
          eim_point_data.comp_index = get_random_int_0_to_n(vars.size()-1);
        }
    }

  comm().broadcast(error_finding_new_node);
  libmesh_error_msg_if(error_finding_new_node, "Could not find new node in get_random_point()");

  // Broadcast the values computed above from rank 0
  comm().broadcast(eim_point_data.node_id);
  comm().broadcast(eim_point_data.comp_index);

  return eim_point_data;
}

EimPointData RBEIMConstruction::get_random_point_from_training_sample()
{
  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  if (eim_eval.get_parametrized_function().on_mesh_sides())
    return get_random_point(_local_side_parametrized_functions_for_training[0]);
  else if (eim_eval.get_parametrized_function().on_mesh_nodes())
    return get_random_point(_local_node_parametrized_functions_for_training[0]);
  else
    return get_random_point(_local_parametrized_functions_for_training[0]);
}

} // namespace libMesh
