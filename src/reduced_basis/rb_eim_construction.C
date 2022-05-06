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
#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_parametrized_function.h"

// C++ include
#include <limits>
#include <memory>

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

// Implement u <- k*u
template <typename DataMap>
void scale(DataMap & u, const Number k)
{
  for (auto & it : u)
    {
      std::vector<std::vector<Number>> & outer_vec = it.second;
      for (auto & inner_vec : outer_vec)
        for (auto & value : inner_vec)
          {
            value *= k;
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

void scale_node_data_map(RBEIMConstruction::NodeDataMap & u, const Number k)
{
  for (auto & it : u)
    {
      std::vector<Number> & vec = it.second;
      for (auto & value : vec)
        {
          value *= k;
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
                                                       unsigned int training_parameters_random_seed_in,
                                                       bool quiet_mode_in,
                                                       unsigned int Nmax_in,
                                                       Real rel_training_tolerance_in,
                                                       Real abs_training_tolerance_in,
                                                       RBParameters mu_min_in,
                                                       RBParameters mu_max_in,
                                                       std::map<std::string, std::vector<Real>> discrete_parameter_values_in,
                                                       std::map<std::string,bool> log_scaling_in,
                                                       std::map<std::string, std::vector<Number>> * training_sample_list)
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

      std::vector<Real> & lookup_table_param_values =
        libmesh_map_find(discrete_parameter_values_in, lookup_table_param_name);

      // Overwrite the discrete values for lookup_table_param to make sure that
      // it is: 0, 1, 2, ..., size-1.
      std::iota(lookup_table_param_values.begin(), lookup_table_param_values.end(), 0);

      // Also, overwrite n_training_samples_in to make sure it matches
      // lookup_table_size so that we will get full coverage of the
      // lookup table in our training set.
      n_training_samples_in = lookup_table_param_values.size();
    }

  // Initialize the parameter ranges and the parameters themselves
  initialize_parameters(mu_min_in, mu_max_in, discrete_parameter_values_in);

  initialize_training_parameters(this->get_parameters_min(),
                                 this->get_parameters_max(),
                                 n_training_samples_in,
                                 log_scaling_in,
                                 deterministic_training_in);

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

      std::vector<Number> lookup_table_training_samples(n_training_samples_in);
      std::iota(lookup_table_training_samples.begin(), lookup_table_training_samples.end(), 0);

      set_training_parameter_values(lookup_table_param_name, lookup_table_training_samples);
    }
}

Real RBEIMConstruction::train_eim_approximation()
{
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

  _eim_projection_matrix.resize(get_Nmax(),get_Nmax());

  RBEIMEvaluation & rbe = get_rb_eim_evaluation();
  rbe.initialize_parameters(*this);
  rbe.resize_data_structures(get_Nmax());

  // If we are continuing from a previous training run,
  // we might already be at the max number of basis functions.
  // If so, we can just return.
  libmesh_error_msg_if(rbe.get_n_basis_functions() > 0,
                       "Error: We currently only support EIM training starting from an empty basis");

  libMesh::out << std::endl << "---- Performing Greedy EIM basis enrichment ----" << std::endl;
  Real greedy_error = 0.;
  std::vector<RBParameters> greedy_param_list;

  // Initialize the current training index to the index that corresponds
  // to the largest (in terms of infinity norm) function in the training set.
  // We do this to ensure that the first EIM basis function is not zero.
  unsigned int current_training_index = _max_abs_value_in_training_set_index;
  set_params_from_training_set(current_training_index);
  while (true)
    {
      libMesh::out << "Greedily selected parameter vector:" << std::endl;
      print_parameters();
      greedy_param_list.emplace_back(get_parameters());

      libMesh::out << "Enriching the EIM approximation" << std::endl;
      enrich_eim_approximation(current_training_index);
      update_eim_matrices();

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

      // Convergence and/or termination tests
      {
        if (rbe.get_n_basis_functions() >= this->get_Nmax())
          {
            libMesh::out << "Maximum number of basis functions reached: Nmax = "
                          << get_Nmax() << std::endl;
            break;
          }

        // We consider the relative tolerance as relative to the maximum value in the training
        // set, since we assume that this maximum value provides a relevant scaling.
        if (greedy_error < (get_rel_training_tolerance() * get_max_abs_value_in_training_set()))
          {
            libMesh::out << "Relative error tolerance reached." << std::endl;
            break;
          }

        if (greedy_error < get_abs_training_tolerance())
          {
            libMesh::out << "Absolute error tolerance reached." << std::endl;
            break;
          }

        if (rbe.get_n_basis_functions() >= this->get_Nmax())
          {
            libMesh::out << "Maximum number of basis functions reached: Nmax = "
                         << get_Nmax() << std::endl;
            break;
          }

        {
          bool do_exit = false;
          for (auto & param : greedy_param_list)
            if (param == get_parameters())
              {
                libMesh::out << "Exiting greedy because the same parameters were selected twice"
                             << std::endl;
                do_exit = true;
                break;
              }

          if (do_exit)
            break; // out of while
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

Real RBEIMConstruction::train_eim_approximation_with_POD()
{
  LOG_SCOPE("train_eim_approximation_with_POD()", "RBEIMConstruction");

  // _eim_projection_matrix is not used in the POD case, but we resize it here in any case
  // to be consistent with what we do in train_eim_approximation_with_greedy().
  _eim_projection_matrix.resize(get_Nmax(),get_Nmax());

  RBEIMEvaluation & rbe = get_rb_eim_evaluation();
  rbe.initialize_parameters(*this);
  rbe.resize_data_structures(get_Nmax());

  libmesh_error_msg_if(rbe.get_n_basis_functions() > 0,
                       "Error: We currently only support EIM training starting from an empty basis");

  libMesh::out << std::endl << "---- Performing POD EIM basis enrichment ----" << std::endl;

  // Set up the POD "correlation matrix"
  unsigned int n_snapshots = get_n_training_samples();
  DenseMatrix<Number> correlation_matrix(n_snapshots,n_snapshots);

  std::cout << "Start computing correlation matrix" << std::endl;
  for (unsigned int i=0; i<n_snapshots; i++)
    {
      for (unsigned int j=0; j<=i; j++)
        {
          Number inner_prod = 0.;
          if (rbe.get_parametrized_function().on_mesh_sides())
            {
              inner_prod = side_inner_product(
                _local_side_parametrized_functions_for_training[i],
                _local_side_parametrized_functions_for_training[j]);
            }
          else if (rbe.get_parametrized_function().on_mesh_nodes())
            {
              inner_prod = node_inner_product(
                _local_node_parametrized_functions_for_training[i],
                _local_node_parametrized_functions_for_training[j]);
            }
          else
            {
              inner_prod = inner_product(
                _local_parametrized_functions_for_training[i],
                _local_parametrized_functions_for_training[j]);
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

  // compute SVD of correlation matrix
  DenseVector<Real> sigma( n_snapshots );
  DenseMatrix<Number> U( n_snapshots, n_snapshots );
  DenseMatrix<Number> VT( n_snapshots, n_snapshots );
  correlation_matrix.svd(sigma, U, VT );

  libmesh_error_msg_if(sigma(0) == 0., "Zero singular value encountered in POD construction");

  // Add dominant vectors from the POD as basis functions.
  unsigned int j = 0;
  Real rel_err = 0.;
  while (true)
    {
      if (j >= get_Nmax() || j >= n_snapshots)
        {
          libMesh::out << "Maximum number of basis functions (" << j << ") reached." << std::endl;
          break;
        }

      // The "energy" error in the POD approximation is determined by the first omitted
      // singular value, i.e. sigma(j). We normalize by sigma(0), which gives the total
      // "energy", in order to obtain a relative error.
      rel_err = std::sqrt(sigma(j)) / std::sqrt(sigma(0));

      libMesh::out << "Number of basis functions: " << j
                   << ", POD error norm: " << rel_err << std::endl;

      if (rel_err < get_rel_training_tolerance())
        {
          libMesh::out << "Training tolerance reached." << std::endl;
          break;
        }

      if (rbe.get_parametrized_function().on_mesh_sides())
        {
          // Make a "zero clone" by copying to get the same data layout, and then scaling by zero
          SideQpDataMap v = _local_side_parametrized_functions_for_training[j];
          scale(v, 0.);

          for ( unsigned int i=0; i<n_snapshots; ++i )
            add(v, U.el(i, j), _local_side_parametrized_functions_for_training[i] );

          Real norm_v = std::sqrt(sigma(j));
          scale(v, 1./norm_v);

          enrich_eim_approximation_on_sides(v);
          update_eim_matrices();
        }
      else if (rbe.get_parametrized_function().on_mesh_nodes())
        {
          // Make a "zero clone" by copying to get the same data layout, and then scaling by zero
          NodeDataMap v = _local_node_parametrized_functions_for_training[j];
          scale_node_data_map(v, 0.);

          for ( unsigned int i=0; i<n_snapshots; ++i )
            add_node_data_map(v, U.el(i, j), _local_node_parametrized_functions_for_training[i] );

          Real norm_v = std::sqrt(sigma(j));
          scale_node_data_map(v, 1./norm_v);

          enrich_eim_approximation_on_nodes(v);
          update_eim_matrices();
        }
      else
        {
          // Make a "zero clone" by copying to get the same data layout, and then scaling by zero
          QpDataMap v = _local_parametrized_functions_for_training[j];
          scale(v, 0.);

          for ( unsigned int i=0; i<n_snapshots; ++i )
            add(v, U.el(i, j), _local_parametrized_functions_for_training[i] );

          Real norm_v = std::sqrt(sigma(j));
          scale(v, 1./norm_v);

          // We leave v_obs_vals empty for now, but we can support this later
          // by accumulating v_obs_vals in the same way that we accumulate v.
          std::vector<std::vector<Number>> v_obs_vals;
          enrich_eim_approximation_on_interiors(v, v_obs_vals);
          update_eim_matrices();
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

        auto side_fe = c.get_side_fe(var, dim);
        side_fe->get_JxW();
        side_fe->get_xyz();
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
  _eim_projection_matrix.resize(get_Nmax(),get_Nmax());
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
                  best_fit_rhs(i) = side_inner_product(solution_copy, get_rb_eim_evaluation().get_side_basis_function(i));
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
                  best_fit_rhs(i) = node_inner_product(solution_copy, get_rb_eim_evaluation().get_node_basis_function(i));
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
                  best_fit_rhs(i) = inner_product(solution_copy, get_rb_eim_evaluation().get_basis_function(i));
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
      for(unsigned int i : make_range(n_comps))
        {
          if (max_abs_value_per_component_in_training_set[i] == 0.)
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
      for(unsigned int i : make_range(n_comps))
        {
          if (max_abs_value_per_component_in_training_set[i] == 0.)
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
      for(unsigned int i : make_range(n_comps))
        {
          if (max_abs_value_per_component_in_training_set[i] == 0.)
            _component_scaling_in_training_set[i] = 1.;
          else
            _component_scaling_in_training_set[i] = _max_abs_value_in_training_set / max_abs_value_per_component_in_training_set[i];
        }

      _parametrized_functions_for_training_obs_values.resize( get_n_training_samples() );

      // Finally, we also evaluate the parametrized functions for training at the "observation points"
      if (eim_eval.get_n_observation_points() > 0)
        {
          std::vector<dof_id_type> observation_points_elem_ids;
          std::vector<subdomain_id_type> observation_points_sbd_ids;
          initialize_observation_points_data(observation_points_elem_ids, observation_points_sbd_ids);

          for (auto i : make_range(get_n_training_samples()))
            {
              libMesh::out << "Initializing observation values for training sample "
                << (i+1) << " of " << get_n_training_samples() << std::endl;

              set_params_from_training_set(i);

              _parametrized_functions_for_training_obs_values[i] =
                eim_eval.get_parametrized_function().evaluate_at_observation_points(get_parameters(),
                                                                                    eim_eval.get_observation_points(),
                                                                                    observation_points_elem_ids,
                                                                                    observation_points_sbd_ids,
                                                                                    *this);

              libmesh_error_msg_if(_parametrized_functions_for_training_obs_values[i].size() != eim_eval.get_n_observation_points(),
                                  "Number of observation values should match number of observation points");
            }
        }
    }
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
                  for(boundary_id_type side_boundary_id : side_boundary_ids)
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
                  for(boundary_id_type side_boundary_id : side_boundary_ids)
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

      for(dof_id_type node_id : nodes_with_nodesets)
        {
          const Node * node = mesh.node_ptr(node_id);

          if (node->processor_id() != mesh.comm().rank())
            continue;

          binfo.boundary_ids(node, node_boundary_ids);

          bool has_node_boundary_id = false;
          boundary_id_type matching_boundary_id = BoundaryInfo::invalid_id;
          for(boundary_id_type node_boundary_id : node_boundary_ids)
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

void RBEIMConstruction::initialize_observation_points_data(
  std::vector<dof_id_type> & observation_points_elem_ids,
  std::vector<subdomain_id_type> & observation_points_sbd_ids)
{
  LOG_SCOPE("initialize_observation_points_data()", "RBEIMConstruction");

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  if (eim_eval.get_n_observation_points() == 0)
    return;

  libmesh_error_msg_if (eim_eval.get_parametrized_function().on_mesh_sides(),
                        "initialize_observation_points_data() not currently supported for data on element sides");

  libMesh::out << "Initializing observation point locations" << std::endl;

  const MeshBase & mesh = this->get_mesh();

  const std::vector<Point> & observation_points = eim_eval.get_observation_points();

  std::unique_ptr<PointLocatorBase> point_locator = mesh.sub_point_locator();

  observation_points_elem_ids.resize(observation_points.size());
  observation_points_sbd_ids.resize(observation_points.size());

  for (unsigned int obs_pt_index : index_range(observation_points))
    {
      const Point & p = observation_points[obs_pt_index];
      const Elem * elem = (*point_locator)(p);

      libmesh_error_msg_if (!elem, "No element containing observation found");

      observation_points_elem_ids[obs_pt_index] = elem->id();
      observation_points_sbd_ids[obs_pt_index] = elem->subdomain_id();
    }
}

Number
RBEIMConstruction::inner_product(const QpDataMap & v, const QpDataMap & w)
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

          for (unsigned int qp : index_range(JxW))
            val += JxW[qp] * v_qp[qp] * libmesh_conj(w_qp[qp]);
        }
    }

  comm().sum(val);
  return val;
}

Number
RBEIMConstruction::side_inner_product(const SideQpDataMap & v, const SideQpDataMap & w)
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

          for (unsigned int qp : index_range(JxW))
            val += JxW[qp] * v_qp[qp] * libmesh_conj(w_qp[qp]);
        }
    }

  comm().sum(val);
  return val;
}

Number
RBEIMConstruction::node_inner_product(const NodeDataMap & v, const NodeDataMap & w)
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
          val += v_comps[comp] * libmesh_conj(w_comps[comp]);
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

          // If scale_components_in_enrichment() returns true then we
          // apply a scaling to give an approximately uniform scaling
          // for all components.
          Real comp_scaling = 1.;
          if (get_rb_eim_evaluation().scale_components_in_enrichment())
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

void RBEIMConstruction::enrich_eim_approximation(unsigned int training_index)
{
  LOG_SCOPE("enrich_eim_approximation()", "RBEIMConstruction");

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  set_params_from_training_set(training_index);

  if (eim_eval.get_parametrized_function().on_mesh_sides())
    enrich_eim_approximation_on_sides(_local_side_parametrized_functions_for_training[training_index]);
  else if (eim_eval.get_parametrized_function().on_mesh_nodes())
    enrich_eim_approximation_on_nodes(_local_node_parametrized_functions_for_training[training_index]);
  else
    {
      bool has_obs_vals = (eim_eval.get_n_observation_points() > 0);

      std::vector<std::vector<Number>> new_bf_obs_vals;
      if (has_obs_vals)
        new_bf_obs_vals = _parametrized_functions_for_training_obs_values[training_index];

      enrich_eim_approximation_on_interiors(_local_parametrized_functions_for_training[training_index],
                                            new_bf_obs_vals);
    }
}

void RBEIMConstruction::enrich_eim_approximation_on_sides(const SideQpDataMap & side_pf)
{
  // Make a copy of the input parametrized function, since we will modify this below
  // to give us a new basis function.
  SideQpDataMap local_pf = side_pf;

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  // If we have at least one basis function, then we need to use
  // rb_eim_solve() to find the EIM interpolation error. Otherwise,
  // just use solution as is.
  if (eim_eval.get_n_basis_functions() > 0)
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

              if (abs_value > largest_abs_value)
                {
                  largest_abs_value = abs_value;
                  optimal_value = value;
                  optimal_comp = comp;
                  optimal_elem_id = elem_id;
                  optimal_side_index = side_index;
                  optimal_qp = qp;

                  optimal_point_phi_i_qp.resize(phi.size());
                  for(auto i : index_range(phi))
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
  this->comm().broadcast(optimal_subdomain_id, proc_ID_index);
  this->comm().broadcast(optimal_boundary_id, proc_ID_index);
  this->comm().broadcast(optimal_qp, proc_ID_index);
  this->comm().broadcast(optimal_point_perturbs, proc_ID_index);
  this->comm().broadcast(optimal_point_phi_i_qp, proc_ID_index);

  libmesh_error_msg_if(optimal_elem_id == DofObject::invalid_id, "Error: Invalid element ID");

  libmesh_error_msg_if(optimal_value == 0., "New EIM basis function should not be zero");

  // Scale local_pf so that its largest value is 1.0
  scale_parametrized_function(local_pf, 1./optimal_value);

  // Add local_pf as the new basis function and store data
  // associated with the interpolation point.
  eim_eval.add_side_basis_function_and_interpolation_data(local_pf,
                                                          optimal_point,
                                                          optimal_comp,
                                                          optimal_elem_id,
                                                          optimal_side_index,
                                                          optimal_subdomain_id,
                                                          optimal_boundary_id,
                                                          optimal_qp,
                                                          optimal_point_perturbs,
                                                          optimal_point_phi_i_qp);
}

void RBEIMConstruction::enrich_eim_approximation_on_nodes(const NodeDataMap & node_pf)
{
  // Make a copy of the input parametrized function, since we will modify this below
  // to give us a new basis function.
  NodeDataMap local_pf = node_pf;

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  // If we have at least one basis function, then we need to use
  // rb_eim_solve() to find the EIM interpolation error. Otherwise,
  // just use solution as is.
  if (eim_eval.get_n_basis_functions() > 0)
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

          if (abs_value > largest_abs_value)
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

  libmesh_error_msg_if(optimal_value == 0., "New EIM basis function should not be zero");

  // Scale local_pf so that its largest value is 1.0
  scale_node_parametrized_function(local_pf, 1./optimal_value);

  // Add local_pf as the new basis function and store data
  // associated with the interpolation point.
  eim_eval.add_node_basis_function_and_interpolation_data(local_pf,
                                                          optimal_point,
                                                          optimal_comp,
                                                          optimal_node_id,
                                                          optimal_boundary_id);
}

void RBEIMConstruction::enrich_eim_approximation_on_interiors(const QpDataMap & interior_pf,
                                                              const std::vector<std::vector<Number>> & interior_pf_obs_values)
{
  // Make a copy of the input parametrized function, since we will modify this below
  // to give us a new basis function.
  QpDataMap local_pf = interior_pf;

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();

  bool has_obs_vals = (eim_eval.get_n_observation_points() > 0);

  std::vector<std::vector<Number>> new_bf_obs_vals;
  if (has_obs_vals)
    {
      new_bf_obs_vals = interior_pf_obs_values;
    }

  // If we have at least one basis function, then we need to use
  // rb_eim_solve() to find the EIM interpolation error. Otherwise,
  // just use solution as is.
  if (eim_eval.get_n_basis_functions() > 0)
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

      if(has_obs_vals)
        {
          for (unsigned int i=0; i<RB_size; i++)
            for (unsigned int j=0; j<eim_eval.get_n_observation_points(); j++)
              for (unsigned int k=0; k<new_bf_obs_vals[j].size(); k++)
                new_bf_obs_vals[j][k] -= rb_eim_solution(i) * eim_eval.get_observation_values(i,j)[k];
        }
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

  // Initialize largest_abs_value to be negative so that it definitely gets updated.
  Real largest_abs_value = -1.;

  // In order to compute phi_i_qp, we initialize a FEMContext
  FEMContext con(*this);
  for (auto dim : con.elem_dimensions())
    {
      auto fe = con.get_element_fe(/*var=*/0, dim);
      fe->get_phi();
    }

  for (const auto & [elem_id, comp_and_qp] : local_pf)
    {
      // Also initialize phi in order to compute phi_i_qp
      const Elem & elem_ref = get_mesh().elem_ref(elem_id);
      con.pre_fe_reinit(*this, &elem_ref);

      auto elem_fe = con.get_element_fe(/*var=*/0, elem_ref.dim());
      const std::vector<std::vector<Real>> & phi = elem_fe->get_phi();

      elem_fe->reinit(&elem_ref);

      for (const auto & comp : index_range(comp_and_qp))
        {
          const std::vector<Number> & qp_values = comp_and_qp[comp];

          for (auto qp : index_range(qp_values))
            {
              Number value = qp_values[qp];
              Real abs_value = std::abs(value);

              if (abs_value > largest_abs_value)
                {
                  largest_abs_value = abs_value;
                  optimal_value = value;
                  optimal_comp = comp;
                  optimal_elem_id = elem_id;
                  optimal_qp = qp;

                  optimal_point_phi_i_qp.resize(phi.size());
                  for(auto i : index_range(phi))
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

  libmesh_error_msg_if(optimal_elem_id == DofObject::invalid_id, "Error: Invalid element ID");

  libmesh_error_msg_if(optimal_value == 0., "New EIM basis function should not be zero");

  // Scale local_pf so that its largest value is 1.0
  scale_parametrized_function(local_pf, 1./optimal_value);

  // Add local_pf as the new basis function and store data
  // associated with the interpolation point.
  eim_eval.add_basis_function_and_interpolation_data(local_pf,
                                                    optimal_point,
                                                    optimal_comp,
                                                    optimal_elem_id,
                                                    optimal_subdomain_id,
                                                    optimal_qp,
                                                    optimal_point_perturbs,
                                                    optimal_point_phi_i_qp);

  if (has_obs_vals)
    {
      // Apply the scame scaling to new_bf_obs_vals as we did to
      // the new basis function itself
      for (unsigned int i=0; i<new_bf_obs_vals.size(); i++)
        for (unsigned int j=0; j<new_bf_obs_vals[i].size(); j++)
          new_bf_obs_vals[i][j] *= 1./optimal_value;

      eim_eval.add_observation_values_for_basis_function(new_bf_obs_vals);
    }
}

void RBEIMConstruction::update_eim_matrices()
{
  LOG_SCOPE("update_eim_matrices()", "RBEIMConstruction");

  RBEIMEvaluation & eim_eval = get_rb_eim_evaluation();
  unsigned int RB_size = eim_eval.get_n_basis_functions();

  libmesh_assert_msg(RB_size >= 1, "Must have at least 1 basis function.");

  if (eim_eval.get_parametrized_function().on_mesh_sides())
    {
      // update the matrix that is used to evaluate L2 projections
      // into the EIM approximation space
      for (unsigned int i=(RB_size-1); i<RB_size; i++)
        {
          for (unsigned int j=0; j<RB_size; j++)
            {
              Number value = side_inner_product(eim_eval.get_side_basis_function(j),
                                                eim_eval.get_side_basis_function(i));

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
                                                eim_eval.get_node_basis_function(i));

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
                                           eim_eval.get_basis_function(i));

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

} // namespace libMesh
