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
#include <ctime>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>

// rbOOmit includes
#include "libmesh/rb_construction_base.h"

// libMesh includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/parallel.h"
#include "libmesh/petsc_linear_solver.h"
// Includes for template instantiation
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/linear_implicit_system.h"

namespace libMesh
{

// ------------------------------------------------------------
// RBConstructionBase implementation


template <class Base>
RBConstructionBase<Base>::RBConstructionBase (EquationSystems& es,
                                              const std::string& name_in,
                                              const unsigned int number_in)
  : Base(es, name_in, number_in),
    serial_training_set(false),
    alternative_solver("unchanged"),
    training_parameters_initialized(false),
    training_parameters_random_seed(-1), // by default, use std::time to seed RNG
    _deterministic_training_parameter_name("NONE"),
    _deterministic_training_parameter_repeats(1)
{
  training_parameters.clear();
}

template <class Base>
RBConstructionBase<Base>::~RBConstructionBase ()
{
  this->clear();
}

template <class Base>
void RBConstructionBase<Base>::clear ()
{
  // clear the parent data
  Base::clear();
  RBParametrized::clear();

  std::map< std::string, NumericVector<Number>* >::iterator it           = training_parameters.begin();
  std::map< std::string, NumericVector<Number>* >::const_iterator it_end = training_parameters.end();

  for( ; it != it_end; ++it)
  {
    NumericVector<Number>* training_vector = it->second;
    if(training_vector)
    {
      delete training_vector;
      training_vector = NULL;
    }
  }
  training_parameters.clear();
}

template <class Base>
void RBConstructionBase<Base>::init_data ()
{
  Base::init_data();

  // Initialize the inner product storage vector, which is useful for
  // storing intermediate results when evaluating inner products
  inner_product_storage_vector = NumericVector<Number>::build(libMesh::default_solver_package(), this->communicator());
  inner_product_storage_vector->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
}

template <class Base>
void RBConstructionBase<Base>::get_global_max_error_pair(const Parallel::Communicator &communicator,
							 std::pair<unsigned int, Real>& error_pair)
{
  // Set error_pair.second to the maximum global value and also
  // find which processor contains the maximum value
  unsigned int proc_ID_index;
  communicator.maxloc(error_pair.second, proc_ID_index);

  // Then broadcast error_pair.first from proc_ID_index
  communicator.broadcast(error_pair.first, proc_ID_index);
}

template <class Base>
numeric_index_type RBConstructionBase<Base>::get_n_training_samples() const
{
  libmesh_assert(training_parameters_initialized);

  if(training_parameters.empty())
    return 0;

  return training_parameters.begin()->second->size();
}

template <class Base>
numeric_index_type RBConstructionBase<Base>::get_local_n_training_samples() const
{
  libmesh_assert(training_parameters_initialized);
  return training_parameters.begin()->second->local_size();
}

template <class Base>
numeric_index_type RBConstructionBase<Base>::get_first_local_training_index() const
{
  libmesh_assert(training_parameters_initialized);
  return training_parameters.begin()->second->first_local_index();
}

template <class Base>
numeric_index_type RBConstructionBase<Base>::get_last_local_training_index() const
{
  libmesh_assert(training_parameters_initialized);
  return training_parameters.begin()->second->last_local_index();
}

template <class Base>
void RBConstructionBase<Base>::set_params_from_training_set(unsigned int index)
{
  set_parameters(get_params_from_training_set(index));
}

template <class Base>
RBParameters RBConstructionBase<Base>::get_params_from_training_set(unsigned int index)
{
  libmesh_assert(training_parameters_initialized);

  libmesh_assert( (this->get_first_local_training_index() <= index) &&
                  (index < this->get_last_local_training_index()) );

  RBParameters params;
  std::map< std::string, NumericVector<Number>* >::const_iterator it     = training_parameters.begin();
  std::map< std::string, NumericVector<Number>* >::const_iterator it_end = training_parameters.end();
  for( ; it != it_end; ++it)
  {
    std::string param_name = it->first;
    Real param_value = libmesh_real( ( *(it->second) )(index) );

    params.set_value(param_name, param_value);
  }

  return params;
}

template <class Base>
void RBConstructionBase<Base>::set_params_from_training_set_and_broadcast(unsigned int index)
{
  libmesh_assert(training_parameters_initialized);

  unsigned int root_id = 0;
  if( (this->get_first_local_training_index() <= index) &&
      (index < this->get_last_local_training_index()) )
  {
    // Set parameters on only one processor
    set_params_from_training_set(index);

    // set root_id, only non-zero on one processor
    root_id = this->processor_id();
  }

  // broadcast
  this->communicator().max(root_id);
  broadcast_parameters(root_id);
}

template <class Base>
void RBConstructionBase<Base>::initialize_training_parameters(const RBParameters& mu_min,
                                                              const RBParameters& mu_max,
                                                              unsigned int n_training_samples,
                                                              std::map<std::string,bool> log_param_scale,
                                                              bool deterministic)
{
  // Print out some info about the training set initialization
  libMesh::out << "Initializing training parameters with "
               << (deterministic ? "deterministic " : "random " )
               << "training set..." << std::endl;

  std::map<std::string,bool>::iterator it           = log_param_scale.begin();
  std::map<std::string,bool>::const_iterator it_end = log_param_scale.end();
  for(; it != it_end; ++it)
  {
    libMesh::out << "Parameter " << it->first
                 << ": log scaling = " << it->second << std::endl;
  }
  libMesh::out << std::endl;

  if(deterministic)
  {
    generate_training_parameters_deterministic(this->communicator(),
					       log_param_scale,
                                               training_parameters,
                                               n_training_samples,
                                               mu_min,
                                               mu_max,
                                               serial_training_set);
  }
  else
  {
    if(get_deterministic_training_parameter_name() == "NONE")
    {
      // Generate random training samples for all parameters
      generate_training_parameters_random(this->communicator(),
					  log_param_scale,
                                          training_parameters,
                                          n_training_samples,
                                          mu_min,
                                          mu_max,
                                          this->training_parameters_random_seed,
                                          serial_training_set);
    }
    else
    {
      // Here we generate a "partially random" training set.
      // Generate deterministic training samples for specified parameter, random for the rest.
      generate_training_parameters_partially_random(this->communicator(),
						    get_deterministic_training_parameter_name(),
                                                    get_deterministic_training_parameter_repeats(),
                                                    log_param_scale,
                                                    training_parameters,
                                                    n_training_samples,
                                                    mu_min,
                                                    mu_max,
                                                    this->training_parameters_random_seed,
                                                    serial_training_set);
    }

  }

  training_parameters_initialized = true;
}

template <class Base>
void RBConstructionBase<Base>::load_training_set(std::map< std::string, std::vector<Number> >& new_training_set)
{
  // First, make sure that an initial training set has already been
  // generated
  if(!training_parameters_initialized)
  {
    libMesh::out << "Error: load_training_set cannot be used to initialize parameters"
                 << std::endl;
    libmesh_error();
  }

  // Make sure that the training set has the correct number of parameters
  if(new_training_set.size() != get_n_params())
  {
    libMesh::out << "Error: Incorrect number of parameters in load_training_set."
                 << std::endl;
    libmesh_error();
  }

  // Clear the training set
  std::map< std::string, NumericVector<Number>* >::iterator it           = training_parameters.begin();
  std::map< std::string, NumericVector<Number>* >::const_iterator it_end = training_parameters.end();
  for( ; it != it_end; ++it)
  {
    NumericVector<Number>* training_vector = it->second;
    if(training_vector)
    {
      delete training_vector;
      training_vector = NULL;
    }
  }

  // Get the number of local and global training parameters
  numeric_index_type n_local_training_samples  = new_training_set.begin()->second.size();
  numeric_index_type n_global_training_samples = n_local_training_samples;
  this->communicator().sum(n_global_training_samples);

  it = training_parameters.begin();
  for( ; it != it_end; ++it)
  {
    it->second = NumericVector<Number>::build(libMesh::default_solver_package(), this->communicator()).release();
    it->second->init(n_global_training_samples, n_local_training_samples, false, libMeshEnums::PARALLEL);
  }

  it = training_parameters.begin();
  for( ; it != it_end; ++it)
  {
    std::string param_name = it->first;
    NumericVector<Number>* training_vector = it->second;

    numeric_index_type first_index = training_vector->first_local_index();
    for(numeric_index_type i=0; i<n_local_training_samples; i++)
    {
      numeric_index_type index = first_index + i;
      training_vector->set(index, new_training_set[param_name][i]);
    }
  }
}


template <class Base>
void RBConstructionBase<Base>::generate_training_parameters_random(const Parallel::Communicator &communicator,
								   std::map<std::string, bool> log_param_scale,
                                                                   std::map< std::string, NumericVector<Number>* >& training_parameters_in,
                                                                   unsigned int n_training_samples_in,
                                                                   const RBParameters& min_parameters,
                                                                   const RBParameters& max_parameters,
                                                                   int training_parameters_random_seed,
                                                                   bool serial_training_set)
{
  libmesh_assert_equal_to ( min_parameters.n_parameters(), max_parameters.n_parameters() );
  const unsigned int num_params = min_parameters.n_parameters();

  // Clear training_parameters_in
  {
    std::map< std::string, NumericVector<Number>* >::iterator it           = training_parameters_in.begin();
    std::map< std::string, NumericVector<Number>* >::const_iterator it_end = training_parameters_in.end();

    for( ; it != it_end; ++it)
    {
      NumericVector<Number>* training_vector = it->second;
      if(training_vector)
      {
        delete training_vector;
        training_vector = NULL;
      }
    }
    training_parameters_in.clear();
  }

  if (num_params == 0)
    return;

  if (training_parameters_random_seed < 0)
  {
    if(!serial_training_set)
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
      std::srand( static_cast<unsigned>( std::time(0) ));
    }
  }
  else
  {
    if(!serial_training_set)
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

  // initialize training_parameters_in
  {
    RBParameters::const_iterator it     = min_parameters.begin();
    RBParameters::const_iterator it_end = min_parameters.end();
    for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      training_parameters_in[param_name] = NumericVector<Number>::build(libMesh::default_solver_package(), communicator).release();

      if(!serial_training_set)
      {
        // Calculate the number of training parameters local to this processor
        unsigned int n_local_training_samples;
        unsigned int quotient  = n_training_samples_in/communicator.size();
        unsigned int remainder = n_training_samples_in%communicator.size();
        if(communicator.rank() < remainder)
          n_local_training_samples = (quotient + 1);
        else
          n_local_training_samples = quotient;

        training_parameters_in[param_name]->init(n_training_samples_in, n_local_training_samples, false, libMeshEnums::PARALLEL);
      }
      else
      {
        training_parameters_in[param_name]->init(n_training_samples_in, false, libMeshEnums::SERIAL);
      }
    }
  }

  // finally, set the values
  {
    std::map< std::string, NumericVector<Number>* >::iterator it           = training_parameters_in.begin();
    std::map< std::string, NumericVector<Number>* >::const_iterator it_end = training_parameters_in.end();

    for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      NumericVector<Number>* training_vector = it->second;

      numeric_index_type first_index = training_vector->first_local_index();
      for(numeric_index_type i=0; i<training_vector->local_size(); i++)
      {
        numeric_index_type index = first_index + i;
        Real random_number = ((double)std::rand())/RAND_MAX; // in range [0,1]

        // Generate log10 scaled training parameters
        if(log_param_scale[param_name])
        {
          Real log_min   = log10(min_parameters.get_value(param_name));
          Real log_range = log10(max_parameters.get_value(param_name) / min_parameters.get_value(param_name));

          training_vector->set(index, pow(10., log_min + random_number*log_range ) );
        }
        // Generate linearly scaled training parameters
        else
        {
          training_vector->set(index, random_number*(max_parameters.get_value(param_name) - min_parameters.get_value(param_name))
                                        + min_parameters.get_value(param_name));
        }
      }
    }
  }
}

template <class Base>
void RBConstructionBase<Base>::generate_training_parameters_partially_random(const Parallel::Communicator &communicator,
									     const std::string& deterministic_training_parameter_name,
                                                                             const unsigned int deterministic_training_parameter_repeats,
                                                                             std::map<std::string, bool> log_param_scale,
                                                                             std::map< std::string, NumericVector<Number>* >& training_parameters_in,
                                                                             unsigned int n_deterministic_training_samples_in,
                                                                             const RBParameters& min_parameters,
                                                                             const RBParameters& max_parameters,
                                                                             int training_parameters_random_seed,
                                                                             bool serial_training_set)
{
  libmesh_assert_equal_to ( min_parameters.n_parameters(), max_parameters.n_parameters() );
  const unsigned int num_params = min_parameters.n_parameters();

  // Clear training_parameters_in
  {
    std::map< std::string, NumericVector<Number>* >::iterator it           = training_parameters_in.begin();
    std::map< std::string, NumericVector<Number>* >::const_iterator it_end = training_parameters_in.end();

    for( ; it != it_end; ++it)
    {
      NumericVector<Number>* training_vector = it->second;
      if(training_vector)
      {
        delete training_vector;
        training_vector = NULL;
      }
    }
    training_parameters_in.clear();
  }

  if (num_params == 0)
    return;

  if (training_parameters_random_seed < 0)
  {
    if(!serial_training_set)
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
      std::srand( static_cast<unsigned>( std::time(0) ));
    }
  }
  else
  {
    if(!serial_training_set)
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

  // initialize training_parameters_in
  const unsigned int n_training_samples_in = n_deterministic_training_samples_in * deterministic_training_parameter_repeats;
  {
    RBParameters::const_iterator it     = min_parameters.begin();
    RBParameters::const_iterator it_end = min_parameters.end();
    for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      training_parameters_in[param_name] = NumericVector<Number>::build(libMesh::default_solver_package(), communicator).release();

      if(!serial_training_set)
      {
        // Calculate the number of training parameters local to this processor
        unsigned int n_local_training_samples;
        unsigned int quotient  = n_training_samples_in/communicator.size();
        unsigned int remainder = n_training_samples_in%communicator.size();
        if(communicator.rank() < remainder)
          n_local_training_samples = (quotient + 1);
        else
          n_local_training_samples = quotient;

        training_parameters_in[param_name]->init(n_training_samples_in, n_local_training_samples, false, libMeshEnums::PARALLEL);
      }
      else
      {
        training_parameters_in[param_name]->init(n_training_samples_in, false, libMeshEnums::SERIAL);
      }
    }
  }

  // finally, set the values
  bool found_deterministic_parameter = false;
  {
    std::map< std::string, NumericVector<Number>* >::iterator it           = training_parameters_in.begin();
    std::map< std::string, NumericVector<Number>* >::const_iterator it_end = training_parameters_in.end();

    for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      NumericVector<Number>* training_vector = it->second;

      if(param_name == deterministic_training_parameter_name)
      {
        found_deterministic_parameter = true;

        // Copy code from deterministic training
        bool use_log_scaling = log_param_scale[param_name];
        Real min_param = min_parameters.get_value(param_name);
        Real max_param = max_parameters.get_value(param_name);

        numeric_index_type first_index = training_vector->first_local_index();
        for(numeric_index_type i=0; i<training_vector->local_size(); i++)
        {
          numeric_index_type index = first_index+i;
          if(use_log_scaling)
          {
            Real epsilon = 1.e-6; // Prevent rounding errors triggering asserts
            Real log_min   = log10(min_param + epsilon);
            Real log_range = log10( (max_param-epsilon) / (min_param+epsilon) );
            Real step_size = log_range /
              std::max((unsigned int)1,(n_deterministic_training_samples_in-1));

            if(index<(n_training_samples_in-1))
            {
              training_vector->set(index, pow(10., log_min + (index % n_deterministic_training_samples_in)*step_size ));
            }
            else
            {
              // due to rounding error, the last parameter can be slightly
              // bigger than max_parameters, hence snap back to the max
              training_vector->set(index, max_param);
            }
          }
          else
          {
            // Generate linearly scaled training parameters
            Real step_size = (max_param - min_param) /
              std::max((unsigned int)1,(n_deterministic_training_samples_in-1));
            training_vector->set(index, (index % n_deterministic_training_samples_in)*step_size + min_param);
          }
        }

      }
      else // Otherwise, generate random parameters
      {
        numeric_index_type first_index = training_vector->first_local_index();
        for(numeric_index_type i=0; i<training_vector->local_size(); i++)
        {
          numeric_index_type index = first_index + i;
          Real random_number = ((double)std::rand())/RAND_MAX; // in range [0,1]

          // Generate log10 scaled training parameters
          if(log_param_scale[param_name])
          {
            Real log_min   = log10(min_parameters.get_value(param_name));
            Real log_range = log10(max_parameters.get_value(param_name) / min_parameters.get_value(param_name));

            training_vector->set(index, pow(10., log_min + random_number*log_range ) );
          }
          // Generate linearly scaled training parameters
          else
          {
            training_vector->set(index, random_number*(max_parameters.get_value(param_name) - min_parameters.get_value(param_name))
                                          + min_parameters.get_value(param_name));
          }
        }
      }

    }
  }

  if(!found_deterministic_parameter)
  {
    libMesh::out << "Error: The deterministic parameter " << deterministic_training_parameter_name
                 << " was not found! Check the name and try again.";
    libmesh_error();
  }
}

template <class Base>
void RBConstructionBase<Base>::generate_training_parameters_deterministic(const Parallel::Communicator &communicator,
									  std::map<std::string, bool> log_param_scale,
                                                                          std::map< std::string, NumericVector<Number>* >& training_parameters_in,
                                                                          unsigned int n_training_samples_in,
                                                                          const RBParameters& min_parameters,
                                                                          const RBParameters& max_parameters,
                                                                          bool serial_training_set)
{
  libmesh_assert_equal_to ( min_parameters.n_parameters(), max_parameters.n_parameters() );
  const unsigned int num_params = min_parameters.n_parameters();

  if (num_params == 0)
    return;

  if(num_params > 2)
  {
    libMesh::out << "ERROR: Deterministic training sample generation "
                 << " not implemented for more than two parameters." << std::endl;
    libmesh_not_implemented();
  }

  // Clear training_parameters_in
  {
    std::map< std::string, NumericVector<Number>* >::iterator it           = training_parameters_in.begin();
    std::map< std::string, NumericVector<Number>* >::const_iterator it_end = training_parameters_in.end();

    for( ; it != it_end; ++it)
    {
      NumericVector<Number>* training_vector = it->second;
      if(training_vector)
      {
        delete training_vector;
        training_vector = NULL;
      }
    }
  }

  // Initialize training_parameters_in
  {
    RBParameters::const_iterator it     = min_parameters.begin();
    RBParameters::const_iterator it_end = min_parameters.end();
    for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      training_parameters_in[param_name] = NumericVector<Number>::build(libMesh::default_solver_package(), communicator).release();

      if(!serial_training_set)
      {
        // Calculate the number of training parameters local to this processor
        unsigned int n_local_training_samples;
        unsigned int quotient  = n_training_samples_in/communicator.size();
        unsigned int remainder = n_training_samples_in%communicator.size();
        if(communicator.rank() < remainder)
          n_local_training_samples = (quotient + 1);
        else
          n_local_training_samples = quotient;

        training_parameters_in[param_name]->init(n_training_samples_in, n_local_training_samples, false, libMeshEnums::PARALLEL);
      }
      else
      {
        training_parameters_in[param_name]->init(n_training_samples_in, false, libMeshEnums::SERIAL);
      }
    }
  }

  if(num_params == 1)
  {
    NumericVector<Number>* training_vector = training_parameters_in.begin()->second;
    bool use_log_scaling = log_param_scale.begin()->second;
    Real min_param = min_parameters.begin()->second;
    Real max_param = max_parameters.begin()->second;

    numeric_index_type first_index = training_vector->first_local_index();
    for(numeric_index_type i=0; i<training_vector->local_size(); i++)
    {
      numeric_index_type index = first_index+i;
      if(use_log_scaling)
      {
        Real epsilon = 1.e-6; // Prevent rounding errors triggering asserts
        Real log_min   = log10(min_param + epsilon);
        Real log_range = log10( (max_param-epsilon) / (min_param+epsilon) );
        Real step_size = log_range /
          std::max((unsigned int)1,(n_training_samples_in-1));

        if(index<(n_training_samples_in-1))
        {
          training_vector->set(index, pow(10., log_min + index*step_size ));
        }
        else
        {
          // due to rounding error, the last parameter can be slightly
          // bigger than max_parameters, hence snap back to the max
          training_vector->set(index, max_param);
        }
      }
      else
      {
        // Generate linearly scaled training parameters
        Real step_size = (max_param - min_param) /
          std::max((unsigned int)1,(n_training_samples_in-1));
        training_vector->set(index, index*step_size + min_param);
      }
    }
  }


  // This is for two parameters
  if(num_params == 2)
  {
    // First make sure n_training_samples_in is a square number
    unsigned int n_training_parameters_per_var = static_cast<unsigned int>( std::sqrt(static_cast<Real>(n_training_samples_in)) );
    if( (n_training_parameters_per_var*n_training_parameters_per_var) != n_training_samples_in)
    {
      libMesh::out << "Error: Number of training parameters = " << n_training_samples_in << "." << std::endl
                   << "Deterministic training set generation with two parameters requires " << std::endl
                   << "the number of training parameters to be a perfect square." << std::endl;
      libmesh_error();
    }

    // make a matrix to store all the parameters, put them in vector form afterwards
    std::vector< std::vector<Real> > training_parameters_matrix(num_params);

    RBParameters::const_iterator it     = min_parameters.begin();
    RBParameters::const_iterator it_end = min_parameters.end();
    unsigned int i = 0;
    for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      Real min_param         = it->second;
      bool use_log_scaling = log_param_scale[param_name];
      Real max_param = max_parameters.get_value(param_name);

      training_parameters_matrix[i].resize(n_training_parameters_per_var);

      for(unsigned int j=0; j<n_training_parameters_per_var; j++)
      {
          // Generate log10 scaled training parameters
          if(use_log_scaling)
          {
            Real epsilon = 1.e-6; // Prevent rounding errors triggering asserts
            Real log_min   = log10(min_param + epsilon);
            Real log_range = log10( (max_param-epsilon) / (min_param+epsilon) );
            Real step_size = log_range /
              std::max((unsigned int)1,(n_training_parameters_per_var-1));

            if(j<(n_training_parameters_per_var-1))
            {
              training_parameters_matrix[i][j] = pow(10., log_min + j*step_size );
            }
            else
            {
              // due to rounding error, the last parameter can be slightly
              // bigger than max_parameters, hence snap back to the max
              training_parameters_matrix[i][j] = max_param;
            }
          }
          else
          {
            // Generate linearly scaled training parameters
            Real step_size = (max_param - min_param) /
              std::max((unsigned int)1,(n_training_parameters_per_var-1));
            training_parameters_matrix[i][j] = j*step_size + min_param;
          }

      }
      i++;
    }

    // now load into training_samples_in:
    std::map<std::string, NumericVector<Number>*>::iterator new_it = training_parameters_in.begin();

    NumericVector<Number>* training_vector_0 = new_it->second;
    new_it++;
    NumericVector<Number>* training_vector_1 = new_it->second;

    for(unsigned int index1=0; index1<n_training_parameters_per_var; index1++)
    {
      for(unsigned int index2=0; index2<n_training_parameters_per_var; index2++)
      {
        unsigned int index = index1*n_training_parameters_per_var + index2;

        if( (training_vector_0->first_local_index() <= index) &&
            (index < training_vector_0->last_local_index()) )
        {
          training_vector_0->set(index, training_parameters_matrix[0][index1]);
          training_vector_1->set(index, training_parameters_matrix[1][index2]);
        }
      }
    }

//     libMesh::out << "n_training_samples = " << n_training_samples_in << std::endl;
//     for(unsigned int index=0; index<n_training_samples_in; index++)
//     {
//         libMesh::out << "training parameters for index="<<index<<":"<<std::endl;
//         for(unsigned int param=0; param<num_params; param++)
//         {
//           libMesh::out << " " << (*training_parameters_in[param])(index);
//         }
//         libMesh::out << std::endl << std::endl;
//     }

  }
}


template <class Base>
std::pair<std::string,std::string>
RBConstructionBase<Base>::set_alternative_solver
  (AutoPtr<LinearSolver<Number> >&
#ifdef LIBMESH_HAVE_PETSC
    ls
#endif
  )
{
  // It seems that setting it this generic way has no effect...
  // PreconditionerType orig_pc = this->linear_solver->preconditioner_type();
  // this->linear_solver->set_preconditioner_type(AMG_PRECOND);
  // so we do it the "hard" way.
  std::string orig_petsc_pc_type_string, orig_petsc_ksp_type_string;

#ifdef LIBMESH_HAVE_PETSC
  // ... but we can set it the "hard" way
  PetscLinearSolver<Number>* petsc_linear_solver =
    libmesh_cast_ptr<PetscLinearSolver<Number>*>(ls.get());

  // Note: #define PCType char*, and PCGetType just sets a pointer.  We'll use
  // the string below to make a real copy, and set the PC back to its original
  // type at the end of the function.
#if PETSC_VERSION_LESS_THAN(3,0,0) || !PETSC_VERSION_RELEASE
  // Pre-3.0 and petsc-dev (as of October 2012) use non-const versions
  PCType orig_petsc_pc_type;
  KSPType orig_petsc_ksp_type;
#else
  const PCType orig_petsc_pc_type;
  const KSPType orig_petsc_ksp_type;
#endif
  int ierr = 0;

  if (petsc_linear_solver)
    {
      PC pc = petsc_linear_solver->pc();
      ierr = PCGetType(pc, &orig_petsc_pc_type); CHKERRABORT(libMesh::COMM_WORLD,ierr);

      KSP ksp = petsc_linear_solver->ksp();
      ierr = KSPGetType(ksp, &orig_petsc_ksp_type); CHKERRABORT(libMesh::COMM_WORLD,ierr);

      // libMesh::out << "orig_petsc_pc_type (before)=" << orig_petsc_pc_type << std::endl;
      // Make actual copies of the original PC and KSP types
      orig_petsc_pc_type_string = orig_petsc_pc_type;
      orig_petsc_ksp_type_string = orig_petsc_ksp_type;

#ifdef LIBMESH_HAVE_PETSC_HYPRE
      // Set solver/PC combo specified in input file...
      if (this->alternative_solver == "amg")
	{
	  // Set HYPRE and boomeramg PC types
	  ierr = PCSetType(pc, PCHYPRE); CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  ierr = PCHYPRESetType(pc, "boomeramg"); CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}
#endif // LIBMESH_HAVE_PETSC_HYPRE
      if (this->alternative_solver == "mumps")
	{
	  // We'll use MUMPS... TODO: configure libmesh to detect
	  // when MUMPS is available via PETSc.

	  // No initial guesses can be specified with KSPPREONLY.  We
	  // can leave the solver as gmres or whatever and it should
	  // converge in 1 iteration.  Otherwise, to use KSPPREONLY,
	  // you may need to do:
	  // KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	  // ierr = KSPSetType(ksp, KSPPREONLY); CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  // Need to call the equivalent for the command line options:
	  // -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps
	  ierr = PCSetType(pc, PCLU); CHKERRABORT(libMesh::COMM_WORLD,ierr);
#if !(PETSC_VERSION_LESS_THAN(3,0,0))
	  ierr = PCFactorSetMatSolverPackage(pc,"mumps"); CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif
	}
    }
  else
    {
      // Otherwise, the cast failed and we are not using PETSc...
      libMesh::out << "You are not using PETSc, so don't know how to set AMG PC." << std::endl;
      libMesh::out << "Returning empty string!" << std::endl;
    }
#endif // LIBMESH_HAVE_PETSC

  return std::make_pair(orig_petsc_pc_type_string, orig_petsc_ksp_type_string);
}




template <class Base>
void RBConstructionBase<Base>::reset_alternative_solver(
#ifdef LIBMESH_HAVE_PETSC
  AutoPtr<LinearSolver<Number> >& ls,
  const std::pair<std::string,std::string>& orig
#else
  AutoPtr<LinearSolver<Number> >&,
  const std::pair<std::string,std::string>&
#endif
)
{
#ifdef LIBMESH_HAVE_PETSC

  // If we never switched, we don't need to do anything...
  if (this->alternative_solver != "unchanged")
    {
      // this->linear_solver->set_preconditioner_type(orig_pc);
      // Set PC back to its previous type
      PetscLinearSolver<Number>* petsc_linear_solver =
	libmesh_cast_ptr<PetscLinearSolver<Number>*>(ls.get());

      int ierr = 0;
      PC pc;
      KSP ksp;

      if (petsc_linear_solver)
	{
	  pc = petsc_linear_solver->pc();
	  ierr = PCSetType(pc, orig.first.c_str()); CHKERRABORT(libMesh::COMM_WORLD,ierr);

	  ksp = petsc_linear_solver->ksp();
	  ierr = KSPSetType(ksp, orig.second.c_str()); CHKERRABORT(libMesh::COMM_WORLD,ierr);
	}
    }

#endif
}


template <class Base>
void RBConstructionBase<Base>::broadcast_parameters(unsigned int proc_id)
{
  libmesh_assert_less (proc_id, this->n_processors());

  // create a copy of the current parameters
  RBParameters current_parameters = get_parameters();

  // copy current_parameters to current_parameters_vector in order to broadcast
  std::vector<Real> current_parameters_vector;

  RBParameters::const_iterator it           = current_parameters.begin();
  RBParameters::const_iterator it_end = current_parameters.end();

  for( ; it != it_end; ++it)
  {
    current_parameters_vector.push_back(it->second);
  }

  // do the broadcast
  this->communicator().broadcast(current_parameters_vector, proc_id);

  // update the copy of the RBParameters object
  it = current_parameters.begin();
  unsigned int count = 0;
  for( ; it != it_end; ++it)
  {
    std::string param_name = it->first;
    current_parameters.set_value(param_name, current_parameters_vector[count]);
    count++;
  }

  // set the parameters globally
  set_parameters(current_parameters);
}

template <class Base>
void RBConstructionBase<Base>::set_training_random_seed(unsigned int seed)
{
  this->training_parameters_random_seed = seed;
}

template <class Base>
void RBConstructionBase<Base>::set_deterministic_training_parameter_name(const std::string name_in)
{
  this->_deterministic_training_parameter_name = name_in;
}

template <class Base>
const std::string& RBConstructionBase<Base>::get_deterministic_training_parameter_name() const
{
  return this->_deterministic_training_parameter_name;
}

template <class Base>
void RBConstructionBase<Base>::set_deterministic_training_parameter_repeats(unsigned int repeats)
{
  this->_deterministic_training_parameter_repeats = repeats;
}

template <class Base>
unsigned int RBConstructionBase<Base>::get_deterministic_training_parameter_repeats() const
{
  return this->_deterministic_training_parameter_repeats;
}

// Template specializations

// EigenSystem is only defined if we have SLEPc
#if defined(LIBMESH_HAVE_SLEPC)
template class RBConstructionBase<CondensedEigenSystem>;
#endif

template class RBConstructionBase<LinearImplicitSystem>;

} // namespace libMesh
