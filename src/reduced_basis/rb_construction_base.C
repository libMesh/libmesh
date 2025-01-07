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
#include <algorithm>
#include <cstddef>
#include <ctime>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
#include <iterator>
#include <memory>
#include <numeric>

// rbOOmit includes
#include "libmesh/rb_construction_base.h"
#include "libmesh/rb_parameters.h"

// libMesh includes
#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/parallel.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/int_range.h"
#include "libmesh/utility.h"

// Nvidia C++ whining about destroying incomplete unique_ptr<T> Base::foo types
#include "libmesh/dof_map.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "timpi/communicator.h"

// Anonymous namespace
namespace
{

/*
 * Helper function to divide a vector of samples across processors.
 * Returns a pair with num_local_samples and first_local_index.
 */
std::pair<unsigned int, unsigned int> calculate_n_local_samples_and_index(
    const libMesh::Parallel::Communicator &communicator,
    const unsigned int n_global_samples,
    const bool serial)
{
  unsigned int n_local_samples = n_global_samples;
  unsigned int first_local_index = 0;

  if (serial || communicator.size() == 1)
    return {n_local_samples, first_local_index};

  // Calculate the number of training parameters local to this processor
  unsigned int quotient  = n_global_samples/communicator.size();
  unsigned int remainder = n_global_samples%communicator.size();
  if (communicator.rank() < remainder)
    {
      n_local_samples = (quotient + 1);
      first_local_index = communicator.rank()*(quotient+1);
    }
  else
    {
      n_local_samples = quotient;
      first_local_index = communicator.rank()*quotient + remainder;
    }
    return {n_local_samples, first_local_index};
}
}  // end anonymous namespace

namespace libMesh
{

// ------------------------------------------------------------
// RBConstructionBase implementation


template <class Base>
RBConstructionBase<Base>::RBConstructionBase (EquationSystems & es,
                                              const std::string & name_in,
                                              const unsigned int number_in)
  : Base(es, name_in, number_in),
    quiet_mode(true),
    serial_training_set(false),
    _normalize_solution_snapshots(false),
    _training_parameters_initialized(false),
    _first_local_index(0),
    _n_local_training_samples(0),
    _n_global_training_samples(0),
    _training_parameters_random_seed(-1) // by default, use std::time to seed RNG
{
}

template <class Base>
RBConstructionBase<Base>::~RBConstructionBase () = default;

template <class Base>
void RBConstructionBase<Base>::clear ()
{
  // clear the parent data
  Base::clear();
  RBParametrized::clear();
  _training_parameters.clear();
}

template <class Base>
void RBConstructionBase<Base>::init_data ()
{
  Base::init_data();

  // Initialize the inner product storage vector, which is useful for
  // storing intermediate results when evaluating inner products
  inner_product_storage_vector = NumericVector<Number>::build(this->comm());
  inner_product_storage_vector->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
}

template <class Base>
void RBConstructionBase<Base>::get_global_max_error_pair(const Parallel::Communicator & communicator,
                                                         std::pair<numeric_index_type, Real> & error_pair)
{
  // Set error_pair.second to the maximum global value and also
  // find which processor contains the maximum value
  unsigned int proc_ID_index;
  communicator.maxloc(error_pair.second, proc_ID_index);

  // Then broadcast error_pair.first from proc_ID_index
  communicator.broadcast(error_pair.first, proc_ID_index);
}

template <class Base>
void RBConstructionBase<Base>::set_normalize_solution_snapshots(bool value)
{
  _normalize_solution_snapshots = value;
}

template <class Base>
numeric_index_type RBConstructionBase<Base>::get_n_training_samples() const
{
  libmesh_error_msg_if(!_training_parameters_initialized,
                       "Error: training parameters must first be initialized.");

  // First we check if there are no parameters here, and in that case we
  // return 1 since a single training sample is sufficient to generate an
  // RB approximation if there are no parameters. Note that in parallel,
  // and when we don't have a serial training set, set return comm().size()
  // so that each processor is assigned a single (empty) training sample.
  if (_training_parameters.empty())
    {
      if (serial_training_set)
        return 1;
      else
        return this->comm().size();
    }

  return _n_global_training_samples;
}

template <class Base>
numeric_index_type RBConstructionBase<Base>::get_local_n_training_samples() const
{
  libmesh_error_msg_if(!_training_parameters_initialized,
                       "Error: training parameters must first be initialized.");

  // First we check if there are no parameters here, and in that case we
  // return 1 for both serial and parallel training sets. This is consistent
  // with get_n_training_samples(), and avoids accessing
  // training_parameters.begin() when training_parameters is empty.
  if (_training_parameters.empty())
    return 1;

  return _n_local_training_samples;
}

template <class Base>
numeric_index_type RBConstructionBase<Base>::get_first_local_training_index() const
{
  libmesh_error_msg_if(!_training_parameters_initialized,
                       "Error: training parameters must first be initialized.");

  // First we check if there are no parameters here, and in that case we
  // return 0 for a serial training set and comm().rank() for a parallel
  // training set. This is consistent with get_n_training_samples(), and
  // avoids accessing training_parameters.begin() when training_parameters
  // is empty.
  if (_training_parameters.empty())
    {
      if (serial_training_set)
        return 0;
      else
        return this->comm().rank();
    }

  return _first_local_index;
}

template <class Base>
numeric_index_type RBConstructionBase<Base>::get_last_local_training_index() const
{
  libmesh_error_msg_if(!_training_parameters_initialized,
                       "Error: training parameters must first be initialized.");

  if (_training_parameters.empty())
    return 0;

  return _first_local_index + _n_local_training_samples;
}

template <class Base>
void RBConstructionBase<Base>::set_params_from_training_set(unsigned int global_index)
{
  set_parameters(get_params_from_training_set(global_index));
}

template <class Base>
RBParameters RBConstructionBase<Base>::get_params_from_training_set(unsigned int global_index)
{
  libmesh_error_msg_if(!_training_parameters_initialized,
                       "Error: training parameters must first be initialized.");

  // If the _training_parameters are empty, return an empty RBParameters.
  // Otherwise, create a new RBParameters object from the single sample requested.
  RBParameters params;
  if (!_training_parameters.empty())
    {
      libmesh_error_msg_if((global_index < this->get_first_local_training_index()) ||
                           (global_index >= this->get_last_local_training_index()),
                           "Error: index "
                           << global_index
                           << " must be within range: "
                           << this->get_first_local_training_index()
                           << " - "
                           << this->get_last_local_training_index());

      const numeric_index_type local_index = global_index - get_first_local_training_index();
      for (const auto & [param_name, sample_vector] : _training_parameters)
        params.set_value(param_name, sample_vector[local_index]);

      // Copy all extra values into the new RBParameters.
      // We assume that the samples may be indexed differently for extra parameters,
      // so we don't just copy the local_index value.
      const auto & mine = get_parameters();
      for (const auto & [key, extra_sample_vector] :
           as_range(mine.extra_begin(), mine.extra_end()))
        {
          for (const auto idx : index_range(extra_sample_vector))
            params.set_extra_value(key, idx, extra_sample_vector[idx]);
        }
    }

  return params;
}

template <class Base>
void RBConstructionBase<Base>::set_params_from_training_set_and_broadcast(unsigned int global_index)
{
  libmesh_error_msg_if(!_training_parameters_initialized,
                       "Error: training parameters must first be initialized.");

  processor_id_type root_id = 0;
  if ((this->get_first_local_training_index() <= global_index) &&
      (global_index < this->get_last_local_training_index()))
    {
      // Set parameters on only one processor
      set_params_from_training_set(global_index);

      // set root_id, only non-zero on one processor
      root_id = this->processor_id();
    }

  // broadcast
  this->comm().max(root_id);
  broadcast_parameters(root_id);
}

template <class Base>
void RBConstructionBase<Base>::initialize_training_parameters(const RBParameters & mu_min,
                                                              const RBParameters & mu_max,
                                                              const unsigned int n_global_training_samples,
                                                              const std::map<std::string,bool> & log_param_scale,
                                                              const bool deterministic)
{
  if (!is_quiet())
    {
      // Print out some info about the training set initialization
      libMesh::out << "Initializing training parameters with "
                  << (deterministic ? "deterministic " : "random " )
                  << "training set..." << std::endl;

      for (const auto & pr : log_param_scale)
        libMesh::out << "Parameter "
                     << pr.first
                     << ": log scaling = "
                     << pr.second
                     << std::endl;

      libMesh::out << std::endl;
    }

  if (deterministic)
    {
      const auto [first_local_index, last_local_index] =
        generate_training_parameters_deterministic(this->comm(),
                                                   log_param_scale,
                                                   _training_parameters,
                                                   n_global_training_samples,
                                                   mu_min,
                                                   mu_max,
                                                   serial_training_set);
      _first_local_index = first_local_index;
      _n_local_training_samples = last_local_index-first_local_index;
    }
  else
    {
      // Generate random training samples for all parameters
      const auto [first_local_index, last_local_index] =
        generate_training_parameters_random(this->comm(),
                                            log_param_scale,
                                            _training_parameters,
                                            n_global_training_samples,
                                            mu_min,
                                            mu_max,
                                            this->_training_parameters_random_seed,
                                            serial_training_set);
      _first_local_index = first_local_index;
      _n_local_training_samples = last_local_index-first_local_index;
    }
  _n_global_training_samples = _n_local_training_samples;

  if (!serial_training_set)
    this->comm().sum(_n_global_training_samples);

  // For each parameter that only allows discrete values, we "snap" to the nearest
  // allowable discrete value
  if (get_n_discrete_params() > 0)
    {
      for (auto & [param_name, sample_vector] : _training_parameters)
        {
          if (is_discrete_parameter(param_name))
            {
              const std::vector<Real> & discrete_values =
                libmesh_map_find(get_discrete_parameter_values(), param_name);

              for (const auto sample_idx : index_range(sample_vector))
                {
                  // Round all values to the closest discrete value.
                  std::vector<Real> discretized_vector(sample_vector[sample_idx].size());
                  std::transform(sample_vector[sample_idx].cbegin(),
                                 sample_vector[sample_idx].cend(),
                                 discretized_vector.begin(),
                                 [&discrete_values](const Real & val) {
                                   return get_closest_value(val, discrete_values);
                                 });
                  sample_vector[sample_idx] = discretized_vector;
                }
            }
        }
    }

  _training_parameters_initialized = true;
}

template <class Base>
void RBConstructionBase<Base>::load_training_set(const std::map<std::string, std::vector<RBParameter>> & new_training_set)
{
  // Make sure we're running this on all processors at the same time
  libmesh_parallel_only(this->comm());

  // First, make sure that an initial training set has already been generated
  libmesh_error_msg_if(!_training_parameters_initialized,
                       "Error: load_training_set cannot be used to initialize parameters");

  // Make sure that the training set has the correct number of parameters
  const unsigned int n_params = get_n_params();
  libmesh_error_msg_if(new_training_set.size() > n_params,
                       "Error: new_training_set should not have more than get_n_params() parameters.");

  // Check that (new_training_set.size() == get_n_params()) is the same on all processes so that
  // we go into the same branch of the "if" statement below on all processes.
  const bool size_matches = (new_training_set.size() == n_params);
  this->comm().verify(size_matches);

  if (size_matches)
    {
      // If new_training_set stores values for all parameters, then we overwrite
      // _training_parameters with new_training_set.

      // Get the number of local and global training parameters
      _first_local_index = 0;
      _n_local_training_samples =
        cast_int<numeric_index_type>(new_training_set.begin()->second.size());
      _n_global_training_samples = _n_local_training_samples;

      if (!serial_training_set)
        {
          this->comm().sum(_n_global_training_samples);

          // Set the first/last indices.
          std::vector<numeric_index_type> local_sizes (this->n_processors(), 0);
          local_sizes[this->processor_id()] = _n_local_training_samples;
          this->comm().sum(local_sizes);

          // first_local_index is the sum of local_sizes
          // for all processor ids less than ours
          for (auto p : make_range(this->processor_id()))
            _first_local_index += local_sizes[p];
        }

      // Ensure that the parameters are the same.
      for (const auto & pr : _training_parameters)
        libmesh_error_msg_if(!new_training_set.count(pr.first),
                             "Parameters must be identical in order to overwrite dataset.");

      // Copy the values from the new_training_set to the internal training_parameters.
      _training_parameters = new_training_set;
    }
  else
    {
      // If new_training_set stores values for a subset of the parameters, then we keep the
      // length of training_parameters unchanged and overwrite the entries of the specified
      // parameters from new_training_set. Note that we repeatedly loop over new_training_set
      // to fill up the entire length of the sample_vector.
      for (auto & [param_name, sample_vector]: _training_parameters)
        {
          if (new_training_set.count(param_name))
            {
              for (const auto i : make_range(get_local_n_training_samples()))
                {
                  const unsigned int num_new_samples = libmesh_map_find(new_training_set,param_name).size();
                  libmesh_error_msg_if (num_new_samples==0, "new_training_set set should not be empty");

                  const unsigned int new_training_set_index = i % num_new_samples;
                  sample_vector[i] = libmesh_map_find(new_training_set,param_name)[new_training_set_index];
                }
            }
        }
    }
}

template <class Base>
void RBConstructionBase<Base>::set_training_parameter_values(
  const std::string & param_name, const std::vector<RBParameter> & values)
{
  libmesh_error_msg_if(!_training_parameters_initialized,
    "Training parameters must be initialized before calling set_training_parameter_values");
  libmesh_error_msg_if(values.size() != get_local_n_training_samples(),
    "Inconsistent sizes");

  // Copy the new data, overwriting the old data.
  auto & training_vector = libmesh_map_find(_training_parameters, param_name);
  training_vector = values;
}


template <class Base>
std::pair<std::size_t, std::size_t>
RBConstructionBase<Base>::generate_training_parameters_random(const Parallel::Communicator & communicator,
                                                              const std::map<std::string, bool> & log_param_scale,
                                                              std::map<std::string, std::vector<RBParameter>> & local_training_parameters_in,
                                                              const unsigned int n_global_training_samples_in,
                                                              const RBParameters & min_parameters,
                                                              const RBParameters & max_parameters,
                                                              const int training_parameters_random_seed,
                                                              const bool serial_training_set)
{
  const unsigned int num_params = min_parameters.n_parameters();
  libmesh_error_msg_if(num_params!=max_parameters.n_parameters(),
    "Number of parameters must be identical for min/max.");

  // Clear training_parameters_in
  local_training_parameters_in.clear();

  if (num_params == 0)
    return {0,0};

  if (training_parameters_random_seed < 0)
    {
      if (!serial_training_set)
        {
          // seed the random number generator with the system time
          // and the processor ID so that the seed is different
          // on different processors
          std::srand( static_cast<unsigned>( std::time(0)*(1+communicator.rank()) ));
        }
      else
        {
          // seed the random number generator with the system time
          // only so that the seed is the same on all processors
          //
          // Note that we broadcast the time on processor 0 to make
          // sure all processors agree.
          unsigned int current_time = static_cast<unsigned>( std::time(0) );
          communicator.broadcast(current_time, 0);
          std::srand(current_time);
        }
    }
  else
    {
      if (!serial_training_set)
        {
          // seed the random number generator with the provided value
          // and the processor ID so that the seed is different
          // on different processors
          std::srand( static_cast<unsigned>( training_parameters_random_seed*(1+communicator.rank()) ));
        }
      else
        {
          // seed the random number generator with the provided value
          // so that the seed is the same on all processors
          std::srand( static_cast<unsigned>( training_parameters_random_seed ));
        }
    }

  // TODO - we don't support vector-data here yet. This would only apply in the case where
  //        min or max are vector-valued, and all the generated points need to stay within those ranges.
  //        But typically we expect that if we're calling this function, we only have 1 min and 1 max,
  //        so the generated values are single-valued as well. The .get_value() calls will throw an error
  //        if this is not the case.

  // initialize training_parameters_in
  const auto & [n_local_training_samples, first_local_index] =
      calculate_n_local_samples_and_index(communicator, n_global_training_samples_in,
                                          serial_training_set);
  for (const auto & pr : min_parameters)
    local_training_parameters_in[pr.first] = std::vector<RBParameter>(n_local_training_samples);

  // finally, set the values
  for (auto & [param_name, sample_vector] : local_training_parameters_in)
    {
      for (auto i : make_range(n_local_training_samples))
        {
          Real random_number = static_cast<Real>(std::rand()) / RAND_MAX; // in range [0,1]

          // Generate log10 scaled training parameters
          if (libmesh_map_find(log_param_scale, param_name))
            {
              Real log_min   = std::log10(min_parameters.get_value(param_name));
              Real log_range = std::log10(max_parameters.get_value(param_name) / min_parameters.get_value(param_name));

              sample_vector[i] = {std::pow(Real(10.), log_min + random_number*log_range )};
            }
          // Generate linearly scaled training parameters
          else
            {
              sample_vector[i] = {
                  random_number * (max_parameters.get_value(param_name) -
                                   min_parameters.get_value(param_name)) +
                  min_parameters.get_value(param_name)};
            }
        }
    }
  return {first_local_index, first_local_index+n_local_training_samples};
}

template <class Base>
std::pair<std::size_t, std::size_t>
RBConstructionBase<Base>::generate_training_parameters_deterministic(const Parallel::Communicator & communicator,
                                                                     const std::map<std::string, bool> & log_param_scale,
                                                                     std::map<std::string, std::vector<RBParameter>> & local_training_parameters_in,
                                                                     const unsigned int n_global_training_samples_in,
                                                                     const RBParameters & min_parameters,
                                                                     const RBParameters & max_parameters,
                                                                     const bool serial_training_set)
{
  libmesh_assert_equal_to ( min_parameters.n_parameters(), max_parameters.n_parameters() );
  const unsigned int num_params = min_parameters.n_parameters();

  if (num_params == 0)
    return {0,0};

  if (num_params > 3)
    libmesh_not_implemented_msg("ERROR: Deterministic training sample generation "
                                "not implemented for more than three parameters.");

  // TODO - we don't support vector-data here yet. This would only apply in the case where
  //        min or max are vector-valued, and all the generated points need to stay within those ranges.
  //        But typically we expect that if we're calling this function, we only have 1 min and 1 max,
  //        so the generated values are single-valued as well. The .get_value() calls will throw an error
  //        if this is not the case.

  // Reinitialize training_parameters_in (but don't remove existing keys!)
  const auto &[n_local_training_samples, first_local_index] =
      calculate_n_local_samples_and_index(communicator, n_global_training_samples_in,
                                         serial_training_set);
  const auto last_local_index = first_local_index + n_local_training_samples;
  for (const auto & pr : min_parameters)
    local_training_parameters_in[pr.first] = std::vector<RBParameter>(n_local_training_samples);

  // n_training_samples_per_param has 3 entries, but entries after "num_params"
  // are unused so we just set their value to 1. We need to set it to 1 (rather
  // than 0) so that we don't skip the inner part of the triply-nested loop over
  // n_training_samples_per_param below.
  std::vector<unsigned int> n_training_samples_per_param(3);
  for (unsigned int param=0; param<3; param++)
    {
      if (param < num_params)
        {
          n_training_samples_per_param[param] =
            static_cast<unsigned int>( std::round(std::pow(static_cast<Real>(n_global_training_samples_in), 1./num_params)) );
        }
      else
        {
          n_training_samples_per_param[param] = 1;
        }
    }

  {
    // The current implementation assumes that we have the same number of
    // samples in each parameter, so we check that n_training_samples_in
    // is consistent with this assumption.
    unsigned int total_samples_check = 1;
    for (unsigned int n_samples : n_training_samples_per_param)
      {
        total_samples_check *= n_samples;
      }

    libmesh_error_msg_if(total_samples_check != n_global_training_samples_in,
                         "Error: Number of training samples = "
                         << n_global_training_samples_in
                         << " does not enable a uniform grid of samples with "
                         << num_params << " parameters. Try "
                         << total_samples_check << " samples instead?");
  }

  // First we make a list of training samples associated with each parameter,
  // then we take a tensor product to obtain the final set of training samples.
  std::vector<std::vector<Real>> training_samples_per_param(num_params);
  {
    unsigned int i = 0;
    for (const auto & pr : min_parameters)
      {
        const std::string & param_name = pr.first;
        const bool use_log_scaling = libmesh_map_find(log_param_scale, param_name);
        Real min_param = min_parameters.get_value(param_name);
        Real max_param = max_parameters.get_value(param_name);

        training_samples_per_param[i].resize(n_training_samples_per_param[i]);

        for (unsigned int j=0; j<n_training_samples_per_param[i]; j++)
          {
            // Generate log10 scaled training parameters
            if (use_log_scaling)
              {
                Real epsilon = 1.e-6; // Prevent rounding errors triggering asserts
                Real log_min   = std::log10(min_param + epsilon);
                Real log_range = std::log10( (max_param-epsilon) / (min_param+epsilon) );
                Real step_size = log_range /
                  std::max((unsigned int)1,(n_training_samples_per_param[i]-1));

                if (j<(n_training_samples_per_param[i]-1))
                  {
                    training_samples_per_param[i][j] = std::pow(10., log_min + j*step_size );
                  }
                else
                  {
                    // due to rounding error, the last parameter can be slightly
                    // bigger than max_parameters, hence snap back to the max
                    training_samples_per_param[i][j] = max_param;
                  }
              }
            else
              {
                // Generate linearly scaled training parameters
                Real step_size = (max_param - min_param) /
                  std::max((unsigned int)1,(n_training_samples_per_param[i]-1));
                training_samples_per_param[i][j] = j*step_size + min_param;
              }

          }
        i++;
      }
  }

  // Now load into training_samples_in
  {
    std::vector<unsigned int> indices(3);
    unsigned int index_count = 0;
    for (indices[0]=0; indices[0]<n_training_samples_per_param[0]; indices[0]++)
      {
        for (indices[1]=0; indices[1]<n_training_samples_per_param[1]; indices[1]++)
          {
            for (indices[2]=0; indices[2]<n_training_samples_per_param[2]; indices[2]++)
              {
                unsigned int param_count = 0;
                for (const auto & pr : min_parameters)
                  {
                    std::vector<RBParameter> & training_vector =
                      libmesh_map_find(local_training_parameters_in, pr.first);
                    if (first_local_index <= index_count && index_count < last_local_index)
                      training_vector[index_count - first_local_index] =
                        {training_samples_per_param[param_count][indices[param_count]]};

                    param_count++;
                  }
                index_count++;
              }
          }
      }
  }
  return {first_local_index, first_local_index+n_local_training_samples};
}


template <class Base>
void RBConstructionBase<Base>::broadcast_parameters(const unsigned int proc_id)
{
  libmesh_assert_less (proc_id, this->n_processors());

  // create a copy of the current parameters
  RBParameters current_parameters = get_parameters();
  libmesh_error_msg_if(current_parameters.n_samples()!=1,
      "Only single-sample RBParameter objects can be broadcast.");

  // Serialize the current_parameters to current_parameters_vector in order to broadcast.
  // We handle multiple samples and vector values.
  // However, the vector values are assumed to remain the same size across samples.
  const std::size_t nparams = current_parameters.n_parameters();
  const std::size_t nsamples = current_parameters.n_samples();

  // First we get the sizes of all the parameter value vectors.
  std::vector<std::size_t> param_value_sizes;
  param_value_sizes.reserve(nparams);
  for (const auto & pr : current_parameters)
    param_value_sizes.push_back(pr.second[0].size());

  // Broadcast the sizes vector and reserve memory.
  this->comm().broadcast(param_value_sizes, proc_id);
  std::size_t buffsize = std::accumulate(param_value_sizes.cbegin(), param_value_sizes.cend(), 0ul);
  std::vector<Real> serialized_parameters;
  serialized_parameters.reserve(buffsize);

  // Then we serialize the parameters/sample/value vectors into a single vector.
  for (const auto & pr : current_parameters)
    {
      for (const auto sample_idx : make_range(nsamples))
        serialized_parameters.insert(serialized_parameters.end(),
                                     pr.second[sample_idx].cbegin(),
                                     pr.second[sample_idx].cend());
    }

  // Do the broadcasts.
  this->comm().broadcast(serialized_parameters, proc_id);

  // Deserialize into the copy of the RBParameters object.
  std::size_t param_idx = 0;
  auto val_idx = serialized_parameters.cbegin();
  for (const auto & pr : current_parameters)
    {
      const std::size_t param_value_size = param_value_sizes[param_idx];
      for (const auto sample_idx: make_range(nsamples))
        {
          auto end_val_idx = std::next(val_idx,param_value_size);
          RBParameter sample_val(val_idx, end_val_idx);
          current_parameters.set_value(pr.first, sample_idx, sample_val);
          val_idx = end_val_idx;
        }
      ++param_idx;
    }

  // Overwrite the parameters globally.
  set_parameters(current_parameters);
}

template <class Base>
void RBConstructionBase<Base>::set_training_random_seed(int seed)
{
  this->_training_parameters_random_seed = seed;
}

// Template specializations

// EigenSystem is only defined if we have SLEPc
#if defined(LIBMESH_HAVE_SLEPC)
template class LIBMESH_EXPORT RBConstructionBase<CondensedEigenSystem>;
#endif

template class LIBMESH_EXPORT RBConstructionBase<LinearImplicitSystem>;
template class LIBMESH_EXPORT RBConstructionBase<System>;

} // namespace libMesh
