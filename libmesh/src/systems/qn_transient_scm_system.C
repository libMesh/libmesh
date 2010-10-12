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

// Configuration data
#include "libmesh_config.h"

// This class derives from RBSCMSystem, which requires SLEPC
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "qn_transient_scm_system.h"
#include "qn_transient_rb_system.h"

#include "numeric_vector.h"
#include "libmesh_logging.h"
#include "dof_map.h"
#include "getpot.h"
#include "equation_systems.h"
#include "parallel.h"

namespace libMesh
{

QNTransientSCMSystem::QNTransientSCMSystem (EquationSystems& es,
                                            const std::string& name,
                                            const unsigned int number)
  : Parent(es, name, number),
    n_bfs(0),
    _k(0),
    _K(0),
    dt(0.05),
    n_time_samples(1)
{
  training_RB_coeffs.clear();
}

QNTransientSCMSystem::~QNTransientSCMSystem ()
{
  this->clear();
}

void QNTransientSCMSystem::clear()
{
  for(unsigned int i=0; i<training_RB_coeffs.size(); i++)
  {
    if(training_RB_coeffs[i])
    {
      delete training_RB_coeffs[i];
      training_RB_coeffs[i] = NULL;
    }
  }
  training_RB_coeffs.resize(0);

  Parent::clear();
}

void QNTransientSCMSystem::init_data()
{
  // Set parameters before calling Parent::init
  GetPot infile(parameters_filename);

  const unsigned int K_in     = infile("K", _K);
  const Real dt_in            = infile("dt", dt);

  set_K(K_in);
  set_dt(dt_in);

  n_time_samples = infile("n_SCM_time_samples", n_time_samples);

  Parent::init_data();

  libMesh::out << std::endl << "QNTransientSCMSystem parameters:" << std::endl;
  libMesh::out << "K: " << get_K() << std::endl;
  libMesh::out << "dt: " << get_dt() << std::endl;
  libMesh::out << "Number of training samples in time for each parameter: " << n_time_samples << std::endl;
  libMesh::out << std::endl;
}


void QNTransientSCMSystem::resize_to_new_n_bfs()
{
  EquationSystems& es = this->get_equation_systems();
  QNTransientRBSystem& RB_system = es.get_system<QNTransientRBSystem>(RB_system_name);

  // Set n_bfs to be the number of basis functions in the
  // associated RB system
  this->set_n_basis_functions(RB_system.get_n_basis_functions());

  // Clear and resize training_RB_coeffs
  for(unsigned int i=0; i<training_RB_coeffs.size(); i++)
  {
    if(training_RB_coeffs[i])
    {
      delete training_RB_coeffs[i];
      training_RB_coeffs[i] = NULL;
    }
  }
  training_RB_coeffs.resize(get_n_basis_functions());
  for(unsigned int i=0; i<get_n_basis_functions(); i++)
  {
    training_RB_coeffs[i] = NumericVector<Number>::build().release();
    training_RB_coeffs[i]->init(get_n_training_samples(), get_local_n_training_samples(), false, libMeshEnums::PARALLEL);
  }

  // Need to resize the vectors storing SCM data such as B_min/B_max
  this->resize_SCM_vectors();
}

void QNTransientSCMSystem::resize_SCM_vectors()
{
  for(unsigned int i=0; i<C_J_RB_coeffs.size(); i++)
    C_J_RB_coeffs[i].clear();
  C_J_RB_coeffs.clear();

  Parent::resize_SCM_vectors();
}

void QNTransientSCMSystem::perform_SCM_greedy()
{
  // Solve the RB system at each different mu in the SCM training set
  libMesh::out << "Performing RB solves at each parameter in SCM training set..." << std::endl << std::endl;
  EquationSystems& es = this->get_equation_systems();
  QNTransientRBSystem& RB_system = es.get_system<QNTransientRBSystem>(RB_system_name);
  const unsigned int RB_size = RB_system.get_n_basis_functions();

  // Now resize the system based on the number of basis functions in the
  // associated RBSystem
  resize_to_new_n_bfs();

  std::vector<Real> previous_mu(get_n_params()-1);
  for(unsigned int i=this->get_first_local_training_index();
                   i<this->get_last_local_training_index(); i++)
  {
    // Construct a parameter vector that does not include time
    // so we can pass it to system
    std::vector<Real> eigen_mu = this->get_training_parameter(i);
    std::vector<Real> mu(this->get_n_params()-1);
    for(unsigned int j=0; j<mu.size(); j++)
      mu[j] = eigen_mu[j];

    // Create a bool that indicates whether previous_mu == mu
    bool same_mu;
    if(i == this->get_first_local_training_index())
      same_mu = false;
    else
    {
      same_mu = true;
      for(unsigned int j=0; j<mu.size(); j++)
        same_mu = same_mu && (previous_mu[j] == mu[j]);
    }

    // Don't need to re-solve unless mu has changed!
    if( !same_mu )
    {
      RB_system.set_current_parameters(mu);
      RB_system.RB_solve(RB_size);
    }

    unsigned int time_level = round( eigen_mu[this->get_n_params()-1] );

    std::vector<Number> RB_u_euler_theta(RB_size);
    for(unsigned int j=0; j<RB_size; j++)
    {
      // The training parameters should be generated such that
      // we never get time_level == 0
      libmesh_assert(time_level > 0);

      RB_u_euler_theta[j] = RB_system.get_euler_theta() *RB_system.RB_temporal_solution_data[time_level]  (j)
                + (1.-RB_system.get_euler_theta())*RB_system.RB_temporal_solution_data[time_level-1](j);
    }
    this->set_training_RB_coeffs(i, RB_u_euler_theta);

    previous_mu = RB_system.get_current_parameters();
  }

  Parent::perform_SCM_greedy();
}

void QNTransientSCMSystem::attach_theta_c(theta_q_fptr theta_c_in)
{
  theta_c = theta_c_in;
}

Number QNTransientSCMSystem::eval_theta_c()
{
  libmesh_assert(theta_c != NULL);

  return theta_c(current_parameters);
}

Number QNTransientSCMSystem::eval_theta_q_a(unsigned int q)
{
  if(q < get_n_basis_functions())
  {
    libmesh_assert(theta_c != NULL);
    return eval_theta_c()*current_RB_coeffs[q];
  }
  else
  {
    return Parent::eval_theta_q_a(q-get_n_basis_functions());
  }
}

std::vector<Number> QNTransientSCMSystem::get_training_RB_coeffs(unsigned int index)
{
  libmesh_assert(training_parameters_initialized);

  libmesh_assert( (this->get_first_local_training_index() <= index) &&
                  (index < this->get_last_local_training_index()) );

  std::vector<Number> coeffs(this->get_n_basis_functions());
  for(unsigned int i=0; i<coeffs.size(); i++)
    coeffs[i] = (*training_RB_coeffs[i])(index);

  return coeffs;
}

void QNTransientSCMSystem::set_training_RB_coeffs(unsigned int index, std::vector<Number> coeffs)
{
  libmesh_assert(training_parameters_initialized);

  libmesh_assert( (this->get_first_local_training_index() <= index) &&
                  (index < this->get_last_local_training_index()) );

  for(unsigned int i=0; i<coeffs.size(); i++)
  {
    training_RB_coeffs[i]->set(index, coeffs[i]);
  }
}

void QNTransientSCMSystem::add_scaled_symm_Aq(unsigned int q_a, Number scalar)
{
  START_LOG("add_scaled_symm_Aq()", "QNTransientSCMSystem");

  if(q_a < get_n_basis_functions())
  {
    // Load the trilinear form operators from the RBSystem
    EquationSystems& es = this->get_equation_systems();
    QNTransientRBSystem& rb_system = es.get_system<QNTransientRBSystem>(RB_system_name);

    rb_system.add_scaled_Cn(scalar, q_a, matrix_A, true);
  }
  else
  {
    Parent::add_scaled_symm_Aq(q_a-get_n_basis_functions(), scalar);
  }

  STOP_LOG("add_scaled_symm_Aq()", "QNTransientSCMSystem");
}

void QNTransientSCMSystem::load_training_parameter_locally(unsigned int index)
{
  Parent::load_training_parameter_locally(index);

  // Also, load RB_coeffs
  set_current_RB_coeffs( get_training_RB_coeffs(index) );
}

void QNTransientSCMSystem::load_training_parameter_globally(unsigned int index)
{
  Parent::load_training_parameter_globally(index);

  // Also, load RB_coeffs
  unsigned int root_id = 0;
  std::vector<Number> new_RB_coeffs(get_n_basis_functions());
  if( (this->get_first_local_training_index() <= index) &&
      (index < this->get_last_local_training_index()) )
  {
    new_RB_coeffs = get_training_RB_coeffs(index);
    root_id = libMesh::processor_id(); // Only non-zero on one processor
  }
  Parallel::sum(root_id);
  Parallel::broadcast(new_RB_coeffs, root_id);

  // Store in C_J_RB_coeffs
  C_J_RB_coeffs.push_back( new_RB_coeffs );

  // and set current_RB_coeffs
  set_current_RB_coeffs( new_RB_coeffs );
}

void QNTransientSCMSystem::save_current_parameters()
{
  Parent::save_current_parameters();

  // Also, save RB_coeffs
  saved_RB_coeffs = current_RB_coeffs;
}

void QNTransientSCMSystem::reload_current_parameters()
{
  Parent::reload_current_parameters();

  // Also, reload RB_coeffs
  set_current_RB_coeffs( saved_RB_coeffs );
}

void QNTransientSCMSystem::set_current_parameters_from_C_J(unsigned int C_J_index)
{
  Parent::set_current_parameters_from_C_J(C_J_index);

  // Also, set RB_coeffs
  current_RB_coeffs = C_J_RB_coeffs[C_J_index];
}

void QNTransientSCMSystem::initialize_training_parameters(const std::vector<Real>& mu_min_vector,
                                                          const std::vector<Real>& mu_max_vector,
                                                          const unsigned int n_param_samples,
                                                          const std::vector<bool> log_param_scale,
                                                          const bool deterministic)
{
  if( (_K % n_time_samples) != 0)
  {
    libMesh::out << "ERROR: K must be divisible by n_SCM_time_samples" << std::endl;
    libmesh_error();
  }

  this->mu_min_vector = mu_min_vector;
  this->mu_max_vector = mu_max_vector;

  // Since time is a parameter here, we need at least one other parameter
  libmesh_assert(get_n_params() > 1);

  // Initialize a temporary training set (distributed across processors)
  // that doesn't include time
  std::vector< NumericVector<Number>* > dist_temp_training_parameters;

  // Strip off the time parameter from log_param_scale, mu_min/max_vector
  unsigned int temp_n_params = log_param_scale.size()-1;
  std::vector<bool> temp_log_param_scale(temp_n_params);
  std::vector<Real> temp_mu_min_vector(temp_n_params);
  std::vector<Real> temp_mu_max_vector(temp_n_params);
  for(unsigned int i=0; i<temp_n_params; i++)
  {
    temp_log_param_scale[i] = log_param_scale[i];
    temp_mu_min_vector[i]   = mu_min_vector[i];
    temp_mu_max_vector[i]   = mu_max_vector[i];
  }

  if(deterministic)
  {
    generate_training_parameters_deterministic(temp_log_param_scale,
                                               dist_temp_training_parameters,
                                               n_param_samples,
                                               temp_mu_min_vector,
                                               temp_mu_max_vector);
  }
  else
  {
    generate_training_parameters_random(temp_log_param_scale,
                                        dist_temp_training_parameters,
                                        n_param_samples,
                                        temp_mu_min_vector,
                                        temp_mu_max_vector,
					this->training_parameters_random_seed);
  }

  // dist_temp_training_parameters is now a distributed training set, so
  // first localize and then delete the distributed vectors
  std::vector< std::vector<Number> >
    temp_training_parameters( dist_temp_training_parameters.size() );
  for(unsigned int i=0; i<dist_temp_training_parameters.size(); i++)
  {
    dist_temp_training_parameters[i]->localize( temp_training_parameters[i] );
    delete dist_temp_training_parameters[i];
    dist_temp_training_parameters[i] = NULL;
  }
  dist_temp_training_parameters.clear();

  // Now add time to the parameter training set...
  unsigned int n_training_samples_in = n_time_samples*n_param_samples;

  // Calculate the number of training parameters local to this processor
  unsigned int n_local_training_samples;
  unsigned int quotient  = n_training_samples_in/libMesh::n_processors();
  unsigned int remainder = n_training_samples_in%libMesh::n_processors();
  if(libMesh::processor_id() < remainder)
    n_local_training_samples = (quotient + 1);
  else
    n_local_training_samples = quotient;

  // Initialize num_params NumericVectors
  training_parameters.resize(get_n_params());

  for(unsigned int i=0; i<get_n_params(); i++)
  {
    training_parameters[i] = NumericVector<Number>::build().release();
    training_parameters[i]->init(n_training_samples_in, n_local_training_samples, false, libMeshEnums::PARALLEL);
  }

  int time_step_size = round ( _K / n_time_samples );

  unsigned int first_local_index = training_parameters[0]->first_local_index();
  unsigned int last_local_index  = training_parameters[0]->last_local_index();
  unsigned int count = 0;
  for(unsigned int i=0; i<temp_training_parameters[0].size(); i++)
  {
    std::vector<Real> params_i( get_n_params()-1 );
    for(unsigned int j=0; j<(get_n_params()-1); j++)
      params_i[j] = libmesh_real(temp_training_parameters[j][i]);
    params_i.push_back(0.);

    for(unsigned int j=1; j<=n_time_samples; j++)
    {
      if( (first_local_index <= count) && (count < last_local_index) )
      {
        for(unsigned int index=0; index<(get_n_params()-1); index++)
        {
          training_parameters[index]->set(count, params_i[index]);
        }

        // Now overwrite the time entry with the appropriate value
        training_parameters[get_n_params()-1]->set(count, j*time_step_size);
      }
      count++;
    }
  }

  training_parameters_initialized = true;
}

void QNTransientSCMSystem::load_training_set(std::vector< std::vector<Number> >& new_training_set)
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
  // (excluding the time "parameter")
  if(new_training_set.size() != (get_n_params()-1))
  {
    libMesh::out << "Error: Incorrect number of parameters in load_training_set."
                 << std::endl;
    libmesh_error();
  }

  // Make sure n_time_samples is valid
  if( (_K % n_time_samples) != 0)
  {
    libMesh::out << "ERROR: K must be divisible by n_SCM_time_samples" << std::endl;
    libmesh_error();
  }

  // Initialize a temporary training set (distributed across processors)
  // that doesn't include time. This is where we load the training set into.
  std::vector< NumericVector<Number>* > dist_temp_training_parameters(get_n_params()-1);
  unsigned int n_global_training_samples;

  {
    // Get the number of local and global training parameters
    unsigned int n_local_training_samples  = new_training_set[0].size();
    n_global_training_samples = n_local_training_samples;
    Parallel::sum(n_global_training_samples);

    for(unsigned int i=0; i<dist_temp_training_parameters.size(); i++)
    {
      dist_temp_training_parameters[i] = NumericVector<Number>::build().release();
      dist_temp_training_parameters[i]->init(n_global_training_samples, n_local_training_samples, false, libMeshEnums::PARALLEL);
    }

    for(unsigned int j=0; j<dist_temp_training_parameters.size(); j++)
    {
      unsigned int first_index = dist_temp_training_parameters[j]->first_local_index();
      for(unsigned int i=0; i<n_local_training_samples; i++)
      {
        unsigned int index = first_index + i;
        dist_temp_training_parameters[j]->set(index, new_training_set[j][i]);
      }
    }
  }

  // THE CODE BELOW IS MOSTLY COPIED FROM QNTransientSCMSystem::initialize_training_parameters
  // dist_temp_training_parameters is now a distributed training set, so
  // first localize and then delete the distributed vectors
  std::vector< std::vector<Number> >
    temp_training_parameters( dist_temp_training_parameters.size() );
  for(unsigned int i=0; i<dist_temp_training_parameters.size(); i++)
  {
    dist_temp_training_parameters[i]->localize( temp_training_parameters[i] );
    delete dist_temp_training_parameters[i];
    dist_temp_training_parameters[i] = NULL;
  }
  dist_temp_training_parameters.clear();

  // Now add time to the parameter training set...
  unsigned int n_training_samples_in = n_time_samples*n_global_training_samples;

  // Calculate the number of training parameters local to this processor
  unsigned int n_local_training_samples;
  unsigned int quotient  = n_training_samples_in/libMesh::n_processors();
  unsigned int remainder = n_training_samples_in%libMesh::n_processors();
  if(libMesh::processor_id() < remainder)
    n_local_training_samples = (quotient + 1);
  else
    n_local_training_samples = quotient;

  // Clear the training set
  for(unsigned int i=0; i<training_parameters.size(); i++)
  {
    if(training_parameters[i])
    {
      delete training_parameters[i];
      training_parameters[i] = NULL;
    }
  }

  // Reinitialize the training set
  training_parameters.resize(get_n_params());
  for(unsigned int i=0; i<get_n_params(); i++)
  {
    training_parameters[i] = NumericVector<Number>::build().release();
    training_parameters[i]->init(n_training_samples_in, n_local_training_samples, false, libMeshEnums::PARALLEL);
  }

  int time_step_size = round ( _K / n_time_samples );

  unsigned int first_local_index = training_parameters[0]->first_local_index();
  unsigned int last_local_index  = training_parameters[0]->last_local_index();
  unsigned int count = 0;
  for(unsigned int i=0; i<temp_training_parameters[0].size(); i++)
  {
    std::vector<Real> params_i( get_n_params()-1 );
    for(unsigned int j=0; j<(get_n_params()-1); j++)
      params_i[j] = libmesh_real(temp_training_parameters[j][i]);
    params_i.push_back(0.);

    for(unsigned int j=1; j<=n_time_samples; j++)
    {
      if( (first_local_index <= count) && (count < last_local_index) )
      {
        for(unsigned int index=0; index<(get_n_params()-1); index++)
        {
          training_parameters[index]->set(count, params_i[index]);
        }

        // Now overwrite the time entry with the appropriate value
        training_parameters[get_n_params()-1]->set(count, j*time_step_size);
      }
      count++;
    }
  }

//  // Print out the training set
//  libMesh::out << "n_training_samples = " << n_training_samples_in << std::endl;
//  for(unsigned int index=0; index<n_training_samples_in; index++)
//  {
//    libMesh::out << "training parameters for index="<<index<<":"<<std::endl;
//    for(unsigned int param=0; param<get_n_params(); param++)
//    {
//      libMesh::out << " " << (*training_parameters[param])(index);
//    }
//    libMesh::out << std::endl << std::endl;
//  }

  // Finally, clear and reinitialize vectors to store the RB coefficients
  for(unsigned int i=0; i<training_RB_coeffs.size(); i++)
  {
    if(training_RB_coeffs[i])
    {
      delete training_RB_coeffs[i];
      training_RB_coeffs[i] = NULL;
    }
  }
  training_RB_coeffs.resize(get_n_basis_functions());
  for(unsigned int i=0; i<get_n_basis_functions(); i++)
  {
    training_RB_coeffs[i] = NumericVector<Number>::build().release();
    training_RB_coeffs[i]->init(n_training_samples_in, n_local_training_samples, false, libMeshEnums::PARALLEL);
  }

  training_parameters_initialized = true;
}

void QNTransientSCMSystem::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "QNTransientSCMSystem");

  Parent::write_offline_data_to_files(directory_name);

  const unsigned int precision_level = 14;

  if(libMesh::processor_id() == 0)
  {
    // Write out C_J_RB_coeffs
    std::ofstream C_J_RB_coeffs_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/C_J_RB_coeffs.dat";
      C_J_RB_coeffs_out.open(file_name.str().c_str());
    }
    if ( !C_J_RB_coeffs_out.good() )
    {
      libMesh::err << "Error opening C_J_RB_coeffs.dat" << std::endl;
      libmesh_error();
    }
    C_J_RB_coeffs_out.precision(precision_level);
    for(unsigned int i=0; i<C_J_RB_coeffs.size(); i++)
    {
      for(unsigned int j=0; j<get_n_basis_functions(); j++)
      {
        C_J_RB_coeffs_out << std::scientific << C_J_RB_coeffs[i][j] << " ";
      }
    }
    C_J_RB_coeffs_out << std::endl;
    C_J_RB_coeffs_out.close();
  }

  STOP_LOG("write_offline_data_to_files()", "QNTransientSCMSystem");
}

void QNTransientSCMSystem::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "QNTransientSCMSystem");

  // Need to make sure we know how many basis functions are in the associated system
  resize_to_new_n_bfs();
  if(get_n_basis_functions() == 0)
  {
    libMesh::out << "Error: The RBSystem associated to this QNTransientSCMSystem needs to be initialized first.";
    libmesh_error();
  }

  Parent::read_offline_data_from_files(directory_name);

  // Read in C_J_RB_coeffs
  std::ifstream C_J_RB_coeffs_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/C_J_RB_coeffs.dat";
    C_J_RB_coeffs_in.open(file_name.str().c_str());
  }
  if ( !C_J_RB_coeffs_in.good() )
  {
    libMesh::err << "Error opening C_J_RB_coeffs.dat" << std::endl;
    libmesh_error();
  }
  C_J_RB_coeffs.resize( C_J_stability_vector.size() );
  for(unsigned int i=0; i<C_J_RB_coeffs.size(); i++)
  {
    C_J_RB_coeffs[i].resize(get_n_basis_functions());
    for(unsigned int j=0; j<get_n_basis_functions(); j++)
    {
      C_J_RB_coeffs_in >> C_J_RB_coeffs[i][j];
    }
  }
  C_J_RB_coeffs_in.close();

  STOP_LOG("read_offline_data_from_files()", "QNTransientSCMSystem");
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
