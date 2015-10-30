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
#include "rb_construction_base.h"

// libMesh includes
#include "libmesh_logging.h"
#include "numeric_vector.h"
#include "equation_systems.h"
#include "parallel.h"
#include "petsc_linear_solver.h"
// Includes for template instantiation
#include "condensed_eigen_system.h"
#include "linear_implicit_system.h"

// C++ includes
#include <ctime>

namespace libMesh
{

// ------------------------------------------------------------
// RBConstructionBase implementation


template <class Base>
RBConstructionBase<Base>::RBConstructionBase (EquationSystems& es,
                                              const std::string& name,
                                              const unsigned int number)
  : Base(es, name, number),
    rb_theta_expansion(NULL),
    training_parameters_initialized(false),
    training_parameters_random_seed(-1), // by default, use std::time to seed RNG
    serial_training_set(false),
    alternative_solver("unchanged")
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

  for(unsigned int i=0; i<training_parameters.size(); i++)
  {
    if(training_parameters[i])
    {
      delete training_parameters[i];
      training_parameters[i] = NULL;
    }
  }
  training_parameters.resize(0);
}

template <class Base>
void RBConstructionBase<Base>::set_parameter_range(std::vector<Real> mu_min_in,
                                                   std::vector<Real> mu_max_in)
{
  libmesh_assert( mu_min_in.size() == get_n_params() &&
                  mu_max_in.size() == get_n_params() );

  mu_min_vector = mu_min_in;
  mu_max_vector = mu_max_in;
}

template <class Base>
Real RBConstructionBase<Base>::get_parameter_min(unsigned int i) const
{
  libmesh_assert(i < get_n_params());

  return mu_min_vector[i];
}

template <class Base>
Real RBConstructionBase<Base>::get_parameter_max(unsigned int i) const
{
  libmesh_assert(i < get_n_params());

  return mu_max_vector[i];
}

template <class Base>
void RBConstructionBase<Base>::set_current_parameters(const std::vector<Real>& params)
{
  libmesh_assert( valid_params(params) );

  RBParametrizedObject::set_current_parameters(params);
}

template <class Base>
bool RBConstructionBase<Base>::valid_params(const std::vector<Real>& params)
{
  bool valid = ( params.size() == get_n_params() );

  if(!valid)
  {
   return false;
  }
  else
  {
    for(unsigned int i=0; i<params.size(); i++)
    {
      valid = valid && ( (mu_min_vector[i] <= params[i]) &&
                         (params[i] <= mu_max_vector[i]) );
    }
  }

  return valid;
}

template <class Base>
void RBConstructionBase<Base>::init_data ()
{
  Base::init_data();

  // Initialize the inner product storage vector, which is useful for
  // storing intermediate results when evaluating inner products
  inner_product_storage_vector = NumericVector<Number>::build();
  inner_product_storage_vector->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
}

template <class Base>
void RBConstructionBase<Base>::get_global_max_error_pair(std::pair<unsigned int, Real>& error_pair)
{
  // Set error_pair.second to the maximum global value and also
  // find which processor contains the maximum value
  unsigned int proc_ID_index;
  Parallel::maxloc(error_pair.second, proc_ID_index);

  // Then broadcast error_pair.first from proc_ID_index
  Parallel::broadcast(error_pair.first, proc_ID_index);
}

template <class Base>
std::vector<Real> RBConstructionBase<Base>::get_training_parameter(unsigned int index) const
{
  libmesh_assert(training_parameters_initialized);

  libmesh_assert( (this->get_first_local_training_index() <= index) &&
                  (index < this->get_last_local_training_index()) );

  std::vector<Real> parameter(this->get_n_params());
  for(unsigned int i=0; i<parameter.size(); i++)
    parameter[i] = libmesh_real((*training_parameters[i])(index));

  return parameter;
}

template <class Base>
void RBConstructionBase<Base>::load_training_parameter_locally(unsigned int index)
{
  set_current_parameters( get_training_parameter(index) );
}

template <class Base>
void RBConstructionBase<Base>::load_training_parameter_globally(unsigned int index)
{
  libmesh_assert(training_parameters_initialized);

  unsigned int root_id = 0;
  std::vector<Real> new_param(get_n_params());
  if( (this->get_first_local_training_index() <= index) &&
      (index < this->get_last_local_training_index()) )
  {
    new_param = this->get_training_parameter(index);
    root_id = libMesh::processor_id(); // Only non-zero on one processor
  }
  Parallel::sum(root_id);
  Parallel::broadcast(new_param, root_id);

  set_current_parameters(new_param);
}

template <class Base>
void RBConstructionBase<Base>::initialize_training_parameters(const std::vector<Real>& mu_min_vector,
                                                              const std::vector<Real>& mu_max_vector,
                                                              const unsigned int n_training_samples,
                                                              const std::vector<bool> log_param_scale,
                                                              const bool deterministic)
{
  set_parameter_range(mu_min_vector, mu_max_vector);

  // Print out some info about the training set initialization
  libMesh::out << "Initializing training parameters with "
               << (deterministic ? "deterministic " : "random " )
               << "training set..." << std::endl;
  for(unsigned int i=0; i<get_n_params(); i++)
  {
    libMesh::out << "Parameter " << i
                 << ": log scaling = " << log_param_scale[i] << std::endl;
  }
  libMesh::out << std::endl;

  if(deterministic)
  {
    generate_training_parameters_deterministic(log_param_scale,
                                               training_parameters,
                                               n_training_samples,
                                               mu_min_vector,
                                               mu_max_vector,
                                               serial_training_set);
  }
  else
  {
    generate_training_parameters_random(log_param_scale,
                                        training_parameters,
                                        n_training_samples,
                                        mu_min_vector,
                                        mu_max_vector,
					this->training_parameters_random_seed,
					serial_training_set);
  }

  training_parameters_initialized = true;
}

template <class Base>
void RBConstructionBase<Base>::load_training_set(std::vector< std::vector<Number> >& new_training_set)
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
  for(unsigned int i=0; i<training_parameters.size(); i++)
  {
    if(training_parameters[i])
    {
      delete training_parameters[i];
      training_parameters[i] = NULL;
    }
  }

  // Get the number of local and global training parameters
  unsigned int n_local_training_samples  = new_training_set[0].size();
  unsigned int n_global_training_samples = n_local_training_samples;
  Parallel::sum(n_global_training_samples);

  for(unsigned int i=0; i<get_n_params(); i++)
  {
    training_parameters[i] = NumericVector<Number>::build().release();
    training_parameters[i]->init(n_global_training_samples, n_local_training_samples, false, libMeshEnums::PARALLEL);
  }

  for(unsigned int j=0; j<get_n_params(); j++)
  {
    unsigned int first_index = training_parameters[j]->first_local_index();
    for(unsigned int i=0; i<n_local_training_samples; i++)
    {
      unsigned int index = first_index + i;
      training_parameters[j]->set(index, new_training_set[j][i]);
    }
  }
}


template <class Base>
void RBConstructionBase<Base>::generate_training_parameters_random(const std::vector<bool> log_param_scale,
                                                                   std::vector< NumericVector<Number>* >& training_parameters_in,
                                                                   const unsigned int n_training_samples_in,
                                                                   const std::vector<Real>& min_parameters,
                                                                   const std::vector<Real>& max_parameters,
						                   int training_parameters_random_seed,
						                   bool serial_training_set)
{
  libmesh_assert( min_parameters.size() == max_parameters.size() );
  const unsigned int num_params = min_parameters.size();

  // Clear training_parameters_in
  for(unsigned int i=0; i<training_parameters_in.size(); i++)
  {
    if(training_parameters_in[i])
    {
      delete training_parameters_in[i];
      training_parameters_in[i] = NULL;
    }
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
        std::srand( static_cast<unsigned>( std::time(0)*(1+libMesh::processor_id()) ));
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
        std::srand( static_cast<unsigned>( training_parameters_random_seed*(1+libMesh::processor_id()) ));
      }
      else
      {
        // seed the random number generator with the provided value
        // so that the seed is the same on all processors
        std::srand( static_cast<unsigned>( training_parameters_random_seed ));
      }
    }


  // Initialize num_params NumericVectors
  training_parameters_in.resize(num_params);

  for(unsigned int i=0; i<num_params; i++)
  {
    training_parameters_in[i] = NumericVector<Number>::build().release();
    if(!serial_training_set)
    {
      // Calculate the number of training parameters local to this processor
      unsigned int n_local_training_samples;
      unsigned int quotient  = n_training_samples_in/libMesh::n_processors();
      unsigned int remainder = n_training_samples_in%libMesh::n_processors();
      if(libMesh::processor_id() < remainder)
        n_local_training_samples = (quotient + 1);
      else
        n_local_training_samples = quotient;

      training_parameters_in[i]->init(n_training_samples_in, n_local_training_samples, false, libMeshEnums::PARALLEL);
    }
    else
    {
      training_parameters_in[i]->init(n_training_samples_in, false, libMeshEnums::SERIAL);
    }
  }

  for(unsigned int j=0; j<num_params; j++)
  {
    unsigned int first_index = training_parameters_in[j]->first_local_index();
    for(unsigned int i=0; i<training_parameters_in[j]->local_size(); i++)
    {
      unsigned int index = first_index + i;
      Real random_number = ((double)std::rand())/RAND_MAX; // in range [0,1]

      // Generate log10 scaled training parameters
      if(log_param_scale[j])
      {
        Real log_min   = log10(min_parameters[j]);
        Real log_range = log10(max_parameters[j] / min_parameters[j]);

        training_parameters_in[j]->set(index, pow(10., log_min + random_number*log_range ) );
      }
      // Generate linearly scaled training parameters
      else
      {
        training_parameters_in[j]->set(index, random_number*(max_parameters[j] - min_parameters[j])
                                            + min_parameters[j]);
      }
    }
  }
}

template <class Base>
void RBConstructionBase<Base>::generate_training_parameters_deterministic(const std::vector<bool> log_param_scale,
                                                                          std::vector< NumericVector<Number>* >& training_parameters_in,
                                                                          const unsigned int n_training_samples_in,
                                                                          const std::vector<Real>& min_parameters,
                                                                          const std::vector<Real>& max_parameters,
                                                                          bool serial_training_set)
{
  libmesh_assert( min_parameters.size() == max_parameters.size() );
  const unsigned int num_params = min_parameters.size();

  if (num_params == 0)
    return;

  if(num_params > 2)
  {
    libMesh::out << "ERROR: Deterministic training sample generation "
                 << " not implemented for more than two parameters." << std::endl;
    libmesh_not_implemented();
  }

  // Clear training_parameters_in
  for(unsigned int i=0; i<training_parameters_in.size(); i++)
  {
    if(training_parameters_in[i])
    {
      delete training_parameters_in[i];
      training_parameters_in[i] = NULL;
    }
  }

  // Initialize num_params NumericVectors
  training_parameters_in.resize(num_params);

  for(unsigned int i=0; i<num_params; i++)
  {
    training_parameters_in[i] = NumericVector<Number>::build().release();
    if(!serial_training_set)
    {
      // Calculate the number of training parameters local to this processor
      unsigned int n_local_training_samples;
      unsigned int quotient  = n_training_samples_in/libMesh::n_processors();
      unsigned int remainder = n_training_samples_in%libMesh::n_processors();
      if(libMesh::processor_id() < remainder)
        n_local_training_samples = (quotient + 1);
      else
        n_local_training_samples = quotient;

      training_parameters_in[i]->init(n_training_samples_in, n_local_training_samples, false, libMeshEnums::PARALLEL);
    }
    else
    {
      training_parameters_in[i]->init(n_training_samples_in, false, libMeshEnums::SERIAL);
    }
  }


  if(num_params == 1)
  {
    unsigned int first_index = training_parameters_in[0]->first_local_index();
    for(unsigned int i=0; i<training_parameters_in[0]->local_size(); i++)
    {
      unsigned int index = first_index+i;
      if(log_param_scale[0])
      {
        Real epsilon = 1.e-6; // Prevent rounding errors triggering asserts
        Real log_min   = log10(min_parameters[0] + epsilon);
        Real log_range = log10( (max_parameters[0]-epsilon) / (min_parameters[0]+epsilon) );
        Real step_size = log_range /
          std::max((unsigned int)1,(n_training_samples_in-1));

        if(index<(n_training_samples_in-1))
        {
          training_parameters_in[0]->set(index, pow(10., log_min + index*step_size ));
        }
        else
        {
          // due to rounding error, the last parameter can be slightly
          // bigger than max_parameters, hence snap back to the max
          training_parameters_in[0]->set(index, max_parameters[0]);
        }
      }
      else
      {
        // Generate linearly scaled training parameters
        Real step_size = (max_parameters[0] - min_parameters[0]) /
          std::max((unsigned int)1,(n_training_samples_in-1));
        training_parameters_in[0]->set(index, index*step_size + min_parameters[0]);
      }
    }
  }


  // This is for two parameters
  if(num_params == 2)
  {
    // First make sure n_training_samples_in is a square number
    unsigned int n_training_parameters_per_var = static_cast<unsigned int>( std::sqrt(n_training_samples_in) );
    if( (n_training_parameters_per_var*n_training_parameters_per_var) != n_training_samples_in)
    {
      libMesh::out << "Error: Number of training parameters = " << n_training_samples_in << "." << std::endl
                   << "Deterministic training set generation with two parameters requires " << std::endl
                   << "the number of training parameters to be a perfect square." << std::endl;
      libmesh_error();
    }

    // make a matrix to store all the parameters, put them in vector form afterwards
    std::vector< std::vector<Real> > training_parameters_matrix(num_params);

    for(unsigned int i=0; i<num_params; i++)
    {
      training_parameters_matrix[i].resize(n_training_parameters_per_var);

      for(unsigned int j=0; j<n_training_parameters_per_var; j++)
      {
          // Generate log10 scaled training parameters
          if(log_param_scale[i])
          {
            Real epsilon = 1.e-6; // Prevent rounding errors triggering asserts
            Real log_min   = log10(min_parameters[i] + epsilon);
            Real log_range = log10( (max_parameters[i]-epsilon) / (min_parameters[i]+epsilon) );
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
              training_parameters_matrix[i][j] = max_parameters[i];
            }
          }
          else
          {
            // Generate linearly scaled training parameters
            Real step_size = (max_parameters[i] - min_parameters[i]) /
              std::max((unsigned int)1,(n_training_parameters_per_var-1));
            training_parameters_matrix[i][j] = j*step_size + min_parameters[i];
          }

      }
    }

    // now load into training_samples_in:
    for(unsigned int index1=0; index1<n_training_parameters_per_var; index1++)
    {
      for(unsigned int index2=0; index2<n_training_parameters_per_var; index2++)
      {
        unsigned int index = index1*n_training_parameters_per_var + index2;

        if( (training_parameters_in[0]->first_local_index() <= index) &&
            (index < training_parameters_in[0]->last_local_index()) )
        {
          training_parameters_in[0]->set(index, training_parameters_matrix[0][index1]);
          training_parameters_in[1]->set(index, training_parameters_matrix[1][index2]);
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
RBConstructionBase<Base>::set_alternative_solver(AutoPtr<LinearSolver<Number> >& ls)
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
#if PETSC_VERSION_LESS_THAN(3,0,0)
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
void RBConstructionBase<Base>::reset_alternative_solver(AutoPtr<LinearSolver<Number> >& ls,
					                const std::pair<std::string,std::string>& orig)
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



// Template specializations

// EigenSystem is only defined if we have SLEPc
#if defined(LIBMESH_HAVE_SLEPC)
template class RBConstructionBase<CondensedEigenSystem>;
#endif

template class RBConstructionBase<LinearImplicitSystem>;

} // namespace libMesh
