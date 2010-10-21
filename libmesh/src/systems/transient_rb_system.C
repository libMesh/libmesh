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

// LibMesh includes

#include "numeric_vector.h"
#include "sparse_matrix.h"
#include "dof_map.h"
#include "libmesh_logging.h"
#include "transient_rb_system.h"
#include "linear_solver.h"
#include "o_string_stream.h"
#include "equation_systems.h"
#include "gmv_io.h"
#include "getpot.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "xdr_cxx.h"
#include "timestamp.h"
#include "parallel.h"

#include <fstream>
#include <sstream>

// For checking for the existence of files
#include <sys/stat.h>

// Need SLEPc to get the POD eigenvalues
#if defined(LIBMESH_HAVE_SLEPC)
// LAPACK include (via SLEPc)
#include <petscsys.h>
#include <slepcblaslapack.h>
#endif // LIBMESH_HAVE_SLEPC

namespace libMesh
{

TransientRBSystem::TransientRBSystem (EquationSystems& es,
		    const std::string& name,
		    const unsigned int number)
  : Parent(es, name, number),
    L2_matrix(SparseMatrix<Number>::build()),
    nonzero_initialization(false),
    init_filename(""),
    dt(0.05),
    euler_theta(1.),
    _k(0),
    _K(0),
    POD_tol(-1.),
    max_truth_solves(-1),
    L2_assembly(NULL)
{
  old_RB_solution.resize(0);

  temporal_data.resize(0);

  initial_L2_error_all_N.resize(0);

  Fq_Mq_representor_norms.clear();
  Mq_Mq_representor_norms.clear();
  Aq_Mq_representor_norms.clear();
  
  // Clear the theta and assembly vectors so that we can push_back
  theta_q_m_vector.clear();
  M_q_intrr_assembly_vector.clear();
  M_q_bndry_assembly_vector.clear();
}



TransientRBSystem::~TransientRBSystem ()
{
  this->clear();
}


void TransientRBSystem::clear()
{
  Parent::clear();
  
  // clear the mass matrices
  for(unsigned int q=0; q<get_Q_m(); q++)
  {
    if(M_q_vector[q])
    {
      delete M_q_vector[q];
      M_q_vector[q] = NULL;
    }
  }

  // clear the temporal_data
  for(unsigned int i=0; i<temporal_data.size(); i++)
  {
    if(temporal_data[i])
    {
      temporal_data[i]->clear();
      delete temporal_data[i];
      temporal_data[i] = NULL;
    }
  }
  temporal_data.resize(0);

  // Clear the basis-function-dependent data using a call to
  // the non-virtual helper function
  clear_basis_helper();
}

void TransientRBSystem::clear_basis_function_dependent_data()
{
  // Call non-virtual helper function that performs the clear
  clear_basis_helper();

  Parent::clear_basis_function_dependent_data();
}

void TransientRBSystem::clear_basis_helper()
{
  // Finally delete the representors
  for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
  {
    for(unsigned int i=0; i<M_q_representor[q_m].size(); i++)
    {
      if(M_q_representor[q_m][i])
      {
        delete M_q_representor[q_m][i];
        M_q_representor[q_m][i] = NULL;
      }
    }
  }
}

void TransientRBSystem::init_data ()
{
  // Indicate that we need to compute the RB
  // inner product matrix in this case
  compute_RB_inner_product = true;

  // Read in data from parameters_filename
  GetPot infile(parameters_filename);

  nonzero_initialization = infile("nonzero_initialization",nonzero_initialization);
  init_filename = infile("init_filename",init_filename);

  unsigned int K_in                    = infile("K", _K);
  const Real dt_in                     = infile("dt", dt);
  const Real euler_theta_in            = infile("euler_theta", euler_theta);
  const Real POD_tol_in                = infile("POD_tol", POD_tol);
  const int max_truth_solves_in        = infile("max_truth_solves", max_truth_solves);
  const unsigned int delta_N_in        = infile("delta_N", delta_N);

  set_K(K_in);
  set_dt(dt_in);
  set_euler_theta(euler_theta_in);
  set_POD_tol(POD_tol_in);
  set_max_truth_solves(max_truth_solves_in);
  set_delta_N(delta_N_in);

  // We read this in here also (the member variable is initialized in RBSystem, but
  // we need to know its value here as well in order to initialize properly)
  initialize_calN_dependent_data = infile("initialize_calN_dependent_data",
                                          initialize_calN_dependent_data);
  if(initialize_calN_dependent_data)
  {
    // Call the Parent's initialization routine.
    Parent::init_data();
  }
  else
  {
    // If we do not initialize calN-dependent data, then we
    // must skip initialization of Parent (=TransientSystem<RBSystem>)
    // and just initialize RBSystem.
    RBSystem::init_data();
  }

  libMesh::out << std::endl << "TransientRBSystem parameters:" << std::endl;
  libMesh::out << "Q_m: " << get_Q_m() << std::endl;
  libMesh::out << "Number of time-steps: " << get_K() << std::endl;
  libMesh::out << "dt: " << get_dt() << std::endl;
  libMesh::out << "euler_theta (time discretization parameter): " << get_euler_theta() << std::endl;
  if(get_POD_tol() > 0.)
    libMesh::out << "POD_tol: " << get_POD_tol() << std::endl;
  if(max_truth_solves > 0)
    libMesh::out << "Maximum number of truth solves: " << max_truth_solves << std::endl;
  libMesh::out << "delta_N (number of basis functions to add each POD-Greedy step): " << get_delta_N() << std::endl;
  if(nonzero_initialization)
  {
    libMesh::out << "Reading initial condition from " << init_filename << std::endl;
  }
  else
  {
    libMesh::out << "Using zero initial condition" << std::endl;
  }
  libMesh::out << std::endl;
}

void TransientRBSystem::initialize_RB_system(bool online_mode)
{
  // Resize vectors for storing calN-dependent data but only
  // initialize if initialize_calN_dependent_data == true
  unsigned int n_time_levels = get_K()+1;
  temporal_data.resize(n_time_levels);
  
  // Resize vectors for storing calN-dependent data but only
  // initialize if initialize_calN_dependent_data == true
  M_q_vector.resize(get_Q_m());
  
  if(initialize_calN_dependent_data)
  {
    // Only initialize the mass matrices if we
    // are not in low-memory mode
    if(!low_memory_mode)
    {
      DofMap& dof_map = this->get_dof_map();

      dof_map.attach_matrix(*L2_matrix);
      L2_matrix->init();
      L2_matrix->zero();
      
      for(unsigned int q=0; q<get_Q_m(); q++)
      {
        // Initialize the memory for the matrices
        M_q_vector[q] = SparseMatrix<Number>::build().release();
        dof_map.attach_matrix(*M_q_vector[q]);
        M_q_vector[q]->init();
        M_q_vector[q]->zero();
      }
    }

    for(unsigned int i=0; i<n_time_levels; i++)
    {
      temporal_data[i] = (NumericVector<Number>::build().release());
      temporal_data[i]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    }
  }

  Parent::initialize_RB_system(online_mode);
  
  // Resize M_q_representor, has to be after Parent::init_data call
  // since Nmax is read in in RBSystem
  M_q_representor.resize(get_Q_m());
  for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
  {
    M_q_representor[q_m].resize(Nmax);
  }


  // Now allocate the N (i.e. RB) dependent data structures

  RB_L2_matrix.resize(Nmax,Nmax);

  // Allocate dense matrices for RB solves
  RB_M_q_vector.resize(get_Q_m());
  for(unsigned int q=0; q<get_Q_m(); q++)
  {
    // Initialize the memory for the RB matrices
    RB_M_q_vector[q].resize(Nmax,Nmax);
  }

  // Initialize vectors for the norms of the representors
  Fq_Mq_representor_norms.resize(get_Q_f());
  for(unsigned int i=0; i<get_Q_f(); i++)
  {
    Fq_Mq_representor_norms[i].resize(get_Q_m());
    for(unsigned int j=0; j<get_Q_m(); j++)
    {
      Fq_Mq_representor_norms[i][j].resize(Nmax);
    }
  }

  unsigned int Q_m_hat = get_Q_m()*(get_Q_m()+1)/2;
  Mq_Mq_representor_norms.resize(Q_m_hat);
  for(unsigned int i=0; i<Q_m_hat; i++)
  {
    Mq_Mq_representor_norms[i].resize(Nmax);
    for(unsigned int j=0; j<Nmax; j++)
    {
      Mq_Mq_representor_norms[i][j].resize(Nmax);
    }
  }

  Aq_Mq_representor_norms.resize(get_Q_a());
  for(unsigned int i=0; i<get_Q_a(); i++)
  {
    Aq_Mq_representor_norms[i].resize(get_Q_m());
    for(unsigned int j=0; j<get_Q_m(); j++)
    {
      Aq_Mq_representor_norms[i][j].resize(Nmax);
      for(unsigned int k=0; k<Nmax; k++)
      {
        Aq_Mq_representor_norms[i][j][k].resize(Nmax);
      }
    }
  }

  // Initialize the initial condition storage
  RB_initial_condition_all_N.resize(Nmax);
  for(unsigned int i=0; i<RB_initial_condition_all_N.size(); i++)
  {
    // The i^th row holds a vector of lenght i+1
    RB_initial_condition_all_N[i].resize(i+1);
  }

  RB_ic_proj_rhs_all_N.resize(Nmax);
  initial_L2_error_all_N.resize(Nmax);

  // Resize the RB output vectors
  RB_outputs_all_k.resize(get_n_outputs());
  RB_output_error_bounds_all_k.resize(get_n_outputs());
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    RB_outputs_all_k[n].resize(_K+1);
    RB_output_error_bounds_all_k[n].resize(_K+1);
  }

  // and the truth output vectors
  truth_outputs_all_k.resize(get_n_outputs());
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    truth_outputs_all_k[n].resize(_K+1);
  }
  
  // Now update RB_ic_proj_rhs_all_N if necessary.
  // This allows us to compute the L2 projection
  // of the initial condition into the RB space
  // so that we can continue to enrich a given RB
  // space.
  if(!online_mode)
    if(get_n_basis_functions() > 0)
    {
      // Load the initial condition into the solution vector
      initialize_truth();

      AutoPtr< NumericVector<Number> > temp1 = NumericVector<Number>::build();
      temp1->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

      // First compute the right-hand side vector for the L2 projection
      if(!low_memory_mode)
      {
        L2_matrix->vector_mult(*temp1, *solution);
      }
      else
      {
        assemble_L2_matrix(matrix);
        matrix->vector_mult(*temp1, *solution);
      }

      for(unsigned int i=0; i<get_n_basis_functions(); i++)
      {
        RB_ic_proj_rhs_all_N(i) = temp1->dot(*basis_functions[i]);
      }
    }
}

SparseMatrix<Number>* TransientRBSystem::get_M_q(unsigned int q)
{
  if(low_memory_mode)
  {
    libMesh::err << "Error: The affine matrices are not store in low-memory mode." << std::endl;
    libmesh_error();
  }

  if(q >= get_Q_m())
  {
    libMesh::err << "Error: We must have q < Q_m in get_M_q."
                 << std::endl;
    libmesh_error();
  }

  return M_q_vector[q];
}

Number TransientRBSystem::eval_theta_q_m(unsigned int q)
{
  if(q >= get_Q_m())
  {
    libMesh::err << "Error: We must have q < Q_m in eval_theta_q_m."
                 << std::endl;
    libmesh_error();
  }

  libmesh_assert(theta_q_m_vector[q] != NULL);

  return theta_q_m_vector[q](theta_data);
}

void TransientRBSystem::assemble_L2_matrix(SparseMatrix<Number>* input_matrix)
{
  input_matrix->zero();
  add_scaled_matrix_and_vector(1., L2_assembly, NULL, input_matrix, NULL);
}

void TransientRBSystem::assemble_mass_matrix(SparseMatrix<Number>* input_matrix)
{
  input_matrix->zero();
  add_scaled_mass_matrix(1., input_matrix);
}

void TransientRBSystem::add_scaled_mass_matrix(Number scalar, SparseMatrix<Number>* input_matrix)
{
  if(!low_memory_mode)
  {
    for(unsigned int q=0; q<get_Q_m(); q++)
      input_matrix->add(scalar * eval_theta_q_m(q), *get_M_q(q));
  }
  else
  {
    for(unsigned int q=0; q<get_Q_m(); q++)
      add_scaled_matrix_and_vector(scalar * eval_theta_q_m(q),
                                   M_q_intrr_assembly_vector[q],
                                   M_q_bndry_assembly_vector[q],
                                   input_matrix,
                                   NULL);
  }
}

void TransientRBSystem::mass_matrix_scaled_matvec(Number scalar,
                                                  NumericVector<Number>& dest,
                                                  NumericVector<Number>& arg)
{
  START_LOG("mass_matrix_scaled_matvec()", "TransientRBSystem");
  
  dest.zero();

  AutoPtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build();
  temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
        
  for(unsigned int q=0; q<get_Q_m(); q++)
  {
    if(!low_memory_mode)
    {
      get_M_q(q)->vector_mult(*temp_vec, arg);
    }
    else
    {
      assemble_scaled_matvec(1.,
                             M_q_intrr_assembly_vector[q],
                             M_q_bndry_assembly_vector[q],
                             *temp_vec,
                             arg);
    }
    dest.add(scalar * eval_theta_q_m(q), *temp_vec);
  }

  STOP_LOG("mass_matrix_scaled_matvec()", "TransientRBSystem");
}

void TransientRBSystem::truth_assembly()
{
  START_LOG("truth_assembly()", "TransientRBSystem");

  this->matrix->zero();
  this->rhs->zero();

  if(!low_memory_mode)
  {
    // We should have already assembled the matrices
    // and vectors in the affine expansion, so
    // just use them

    add_scaled_mass_matrix(1./get_dt(), matrix);
    mass_matrix_scaled_matvec(1./get_dt(), *rhs, *current_local_solution);
    
    AutoPtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build();
    temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    {
      matrix->add(get_euler_theta()*eval_theta_q_a(q_a), *get_A_q(q_a));

      get_A_q(q_a)->vector_mult(*temp_vec, *current_local_solution);
      temp_vec->scale( -(1.-get_euler_theta())*eval_theta_q_a(q_a) );
      rhs->add(*temp_vec);
    }

    for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
    {
      *temp_vec = *get_F_q(q_f);
      temp_vec->scale( eval_theta_q_f(q_f) );
      rhs->add(*temp_vec);
    }
    zero_dirichlet_dofs_on_rhs();

    if(constrained_problem)
      matrix->add(1., *constraint_matrix);

  }
  else
  {
    // In low memory mode we do not store the matrices
    // from the affine expansion, so need to assemble

    // For efficiency (i.e. to avoid doing Q_a+Q_f loops
    // over the mesh) we do not use add_scaled_matrix_and_vector
    // here

    const MeshBase& mesh = this->get_mesh();
    
    std::vector<FEMContext*> Mq_context(get_Q_m());
    for(unsigned int q_m=0; q_m<Mq_context.size(); q_m++)
    {
      Mq_context[q_m] = this->build_context().release();
      this->init_context(*Mq_context[q_m]);
    }

    std::vector<FEMContext*> Aq_context(get_Q_a());
    for(unsigned int q_a=0; q_a<Aq_context.size(); q_a++)
    {
      Aq_context[q_a] = this->build_context().release();
      this->init_context(*Aq_context[q_a]);
    }

    std::vector<FEMContext*> Fq_context(get_Q_f());
    for(unsigned int q_f=0; q_f<Fq_context.size(); q_f++)
    {
      Fq_context[q_f] = this->build_context().release();
      this->init_context(*Fq_context[q_f]);
    }

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
      for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
      {
        Mq_context[q_m]->pre_fe_reinit(*this, *el);
        Mq_context[q_m]->elem_fe_reinit();
        if( M_q_intrr_assembly_vector[q_m] != NULL )
        {
          this->M_q_intrr_assembly_vector[q_m](*Mq_context[q_m], *this);
          // Now overwrite the local matrix with a matrix multiplication
          Mq_context[q_m]->elem_jacobian.vector_mult(Mq_context[q_m]->elem_residual, Mq_context[q_m]->elem_solution);
        }
      }

      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        Aq_context[q_a]->pre_fe_reinit(*this, *el);
        Aq_context[q_a]->elem_fe_reinit();
        if( A_q_intrr_assembly_vector[q_a] != NULL )
        {
          this->A_q_intrr_assembly_vector[q_a](*Aq_context[q_a], *this);
          // Now overwrite the local matrix with a matrix multiplication
          Aq_context[q_a]->elem_jacobian.vector_mult(Aq_context[q_a]->elem_residual, Aq_context[q_a]->elem_solution);
        }
      }

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      {
        Fq_context[q_f]->pre_fe_reinit(*this, *el);
        Fq_context[q_f]->elem_fe_reinit();
        if(F_q_intrr_assembly_vector[q_f] != NULL)
          this->F_q_intrr_assembly_vector[q_f](*Fq_context[q_f], *this);
      }

      for (Aq_context[0]->side = 0;
            Aq_context[0]->side != Aq_context[0]->elem->n_sides();
            ++Aq_context[0]->side)
      {
        // May not need to apply fluxes on non-boundary elements
        if( (Aq_context[0]->elem->neighbor(Aq_context[0]->side) != NULL) && !impose_internal_fluxes )
          continue;

        for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
        {
          Mq_context[q_m]->side = Mq_context[0]->side;

          if( M_q_bndry_assembly_vector[q_m] != NULL )
          {
            Mq_context[q_m]->side_fe_reinit();
            this->M_q_bndry_assembly_vector[q_m](*Mq_context[q_m], *this);
          }
        }

        for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
        {
          Aq_context[q_a]->side = Aq_context[0]->side;

          if( A_q_bndry_assembly_vector[q_a] != NULL )
          {
            Aq_context[q_a]->side_fe_reinit();
            this->A_q_bndry_assembly_vector[q_a](*Aq_context[q_a], *this);
          }
        }

        // Impose boundary terms, e.g. Neuman BCs
        for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
        {
          // Update the side information for all contexts
          Fq_context[q_f]->side = Aq_context[0]->side;

          if( F_q_bndry_assembly_vector[q_f] != NULL )
          {
            Fq_context[q_f]->side_fe_reinit();
            this->F_q_bndry_assembly_vector[q_f](*Fq_context[q_f], *this);
          }
        }
      }

      // Constrain the dofs to impose hanging node or periodic constraints
      for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
      {
        this->get_dof_map().constrain_element_matrix_and_vector
          (Mq_context[q_m]->elem_jacobian, Mq_context[q_m]->elem_residual, Mq_context[q_m]->dof_indices);
      }

      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        this->get_dof_map().constrain_element_matrix_and_vector
          (Aq_context[q_a]->elem_jacobian, Aq_context[q_a]->elem_residual, Aq_context[q_a]->dof_indices);
      }

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      {
        this->get_dof_map().constrain_element_matrix_and_vector
          (Fq_context[q_f]->elem_jacobian, Fq_context[q_f]->elem_residual, Fq_context[q_f]->dof_indices);
      }
        
      // Apply Dirichlet boundary conditions, we assume zero Dirichlet BCs
      // Note that this cannot be inside the side-loop since non-boundary
      // elements may contain boundary dofs
      std::set<unsigned int>::const_iterator iter;
      for(unsigned int n=0; n<Aq_context[0]->dof_indices.size(); n++)
      {
	iter = global_dirichlet_dofs_set.find( Aq_context[0]->dof_indices[n] );
	if(iter != global_dirichlet_dofs_set.end())
	{
	  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
	  {
	    Aq_context[q_a]->elem_jacobian.condense
	      (n,n,0.,Aq_context[q_a]->elem_residual);
	  }

	  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
	  {
	    Fq_context[q_f]->elem_jacobian.condense
	      (n,n,0.,Fq_context[q_f]->elem_residual);
	  }

	  for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
	  {
	    Mq_context[q_m]->elem_jacobian.condense
	      (n,n,0.,Mq_context[q_m]->elem_residual);
	  }
	}
      }

      // Finally, add local matrices/vectors to the global system
      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        Aq_context[q_a]->elem_jacobian *= get_euler_theta()*eval_theta_q_a(q_a);
        this->matrix->add_matrix (Aq_context[q_a]->elem_jacobian,
                                  Aq_context[q_a]->dof_indices);
        Aq_context[q_a]->elem_residual *= -(1.-get_euler_theta())*eval_theta_q_a(q_a);
        this->rhs->add_vector    (Aq_context[q_a]->elem_residual,
                                  Aq_context[q_a]->dof_indices);
      }

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      {
        Fq_context[q_f]->elem_residual *= eval_theta_q_f(q_f);
        this->rhs->add_vector (Fq_context[q_f]->elem_residual,
                               Fq_context[q_f]->dof_indices);
      }

      for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
      {
        Mq_context[q_m]->elem_jacobian *= 1./get_dt()*eval_theta_q_m(q_m);
        this->matrix->add_matrix (Mq_context[q_m]->elem_jacobian,
                                  Mq_context[q_m]->dof_indices);
        Mq_context[q_m]->elem_residual *= 1./get_dt()*eval_theta_q_m(q_m);
        this->rhs->add_vector    (Mq_context[q_m]->elem_residual,
                                  Mq_context[q_m]->dof_indices);
      }
    }

    if(constrained_problem)
      add_scaled_matrix_and_vector(1., constraint_assembly, NULL, matrix, NULL);

    // Delete all the ptrs to FEMContexts!
    for(unsigned int q_a=0; q_a<Aq_context.size(); q_a++)
    {
      delete Aq_context[q_a];
      Aq_context[q_a] = NULL;
    }
    Aq_context.clear();

    for(unsigned int q_f=0; q_f<Fq_context.size(); q_f++)
    {
      delete Fq_context[q_f];
      Fq_context[q_f] = NULL;
    }
    Fq_context.clear();

    for(unsigned int q_m=0; q_m<Mq_context.size(); q_m++)
    {
      delete Mq_context[q_m];
      Mq_context[q_m] = NULL;
    }
    Mq_context.clear();
  }

  this->matrix->close();
  this->rhs->close();

  STOP_LOG("truth_assembly()", "TransientRBSystem");
}

void TransientRBSystem::attach_L2_assembly(affine_assembly_fptr L2_assembly_in)
{
  L2_assembly = L2_assembly_in;
}

void TransientRBSystem::attach_M_q(theta_q_fptr theta_q_m,
                                   affine_assembly_fptr M_q_intrr_assembly,
                                   affine_assembly_fptr M_q_bndry_assembly)
{
  theta_q_m_vector.push_back(theta_q_m);
  M_q_intrr_assembly_vector.push_back(M_q_intrr_assembly);
  M_q_bndry_assembly_vector.push_back(M_q_bndry_assembly);
}

void TransientRBSystem::assemble_Mq_matrix(unsigned int q, SparseMatrix<Number>* input_matrix)
{
  if(q >= get_Q_m())
  {
    libMesh::err << "Error: We must have q < Q_m in assemble_Mq_matrix."
                 << std::endl;
    libmesh_error();
  }

  input_matrix->zero();
  add_scaled_matrix_and_vector(1., M_q_intrr_assembly_vector[q], M_q_bndry_assembly_vector[q], input_matrix, NULL);
}

void TransientRBSystem::assemble_all_affine_operators()
{
  Parent::assemble_all_affine_operators();

  for(unsigned int q=0; q<get_Q_m(); q++)
    assemble_Mq_matrix(q, get_M_q(q));
}

void TransientRBSystem::assemble_misc_matrices()
{
  assemble_L2_matrix(L2_matrix.get());

  Parent::assemble_misc_matrices();
}

Real TransientRBSystem::truth_solve(int write_interval)
{
  START_LOG("truth_solve()", "TransientRBSystem");

  const MeshBase& mesh = get_mesh();

//   // NumericVector for computing true L2 error
//   AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build();
//   temp->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

  // Apply initial condition again.
  initialize_truth();
  set_time_level(0);

  // Now compute the truth outputs
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    truth_outputs_all_k[n][_k] = 0.;
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      truth_outputs_all_k[n][_k] += eval_theta_q_l(n,q_l)*get_output_vector(n,q_l)->dot(*solution);
    }
  }

  // Load initial projection error into temporal_data dense matrix
  set_error_temporal_data();

  if(reuse_preconditioner)
  {
    // For the first solve, make sure we generate a new preconditioner
    linear_solver->same_preconditioner = false;
  }

  for(unsigned int time_level=1; time_level<=_K; time_level++)
  {
    set_time_level(time_level); // update the member variable _k

    *old_local_solution = *current_local_solution;

    // We assume that the truth assembly has been attached to the system
    truth_assembly();
    solve();

    // Make sure we didn't max out the number of iterations
    if( (this->n_linear_iterations() >=
        this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
        (this->final_linear_residual() >
        this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
    {
      libMesh::out << "Warning: Linear solver may not have converged! Final linear residual = "
                   << this->final_linear_residual() << ", number of iterations = "
                   << this->n_linear_iterations() << std::endl << std::endl;
//       libmesh_error();
    }

    if(reuse_preconditioner)
      {
	linear_solver->same_preconditioner = true;
      }

    // Now compute the truth outputs
    for(unsigned int n=0; n<get_n_outputs(); n++)
    {
      truth_outputs_all_k[n][_k] = 0.;
      for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      {
        truth_outputs_all_k[n][_k] += eval_theta_q_l(n,q_l)*get_output_vector(n,q_l)->dot(*solution);
      }
    }

    // load projection error into column _k of temporal_data matrix
    set_error_temporal_data();

    if ( (write_interval > 0) && (_k%write_interval == 0) )
      {
        libMesh::out << std::endl << "Truth solve, plotting time step " << _k << std::endl;

        OStringStream file_name;

        // We write the file name in the gmv auto-read format.
        file_name << "truth.gmv.";
        OSSRealzeroright(file_name,3,0, _k);

        GMVIO(mesh).write_equation_systems (file_name.str(),
                                            this->get_equation_systems());
      }
  }

  if(reuse_preconditioner)
    {
      linear_solver->same_preconditioner = false;
    }

  // Get the L2 norm of the truth solution at time-level _K
  // Useful for normalizing our true error data
  if(!low_memory_mode)
  {
    L2_matrix->vector_mult(*inner_product_storage_vector, *solution);
  }
  else
  {
    assemble_L2_matrix(matrix);
    matrix->vector_mult(*inner_product_storage_vector, *solution);
  }

  Real final_truth_L2_norm = libmesh_real(std::sqrt(inner_product_storage_vector->dot(*solution)));


  STOP_LOG("truth_solve()", "TransientRBSystem");

  return final_truth_L2_norm;
}

bool TransientRBSystem::greedy_termination_test(Real training_greedy_error, int count)
{
  if ( (get_max_truth_solves()>0) && (count >= get_max_truth_solves()) )
    {
      libMesh::out << "Maximum number of truth solves reached: max = "
                   << count << std::endl;
      return true;
    }

  return Parent::greedy_termination_test(training_greedy_error, count);
}

Number TransientRBSystem::set_error_temporal_data()
{
  START_LOG("set_error_temporal_data()", "TransientRBSystem");

  // first compute the projection of solution onto the current
  // RB space

  if(basis_functions.size() == 0)
  {
    // If the basis is empty, then the error is the solution itself
    temporal_data[_k]->zero();
    temporal_data[_k]->add(1., *solution);
  }
  else
  {
    unsigned int RB_size = get_n_basis_functions();

    AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build();
    temp->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

    // First compute the right-hand side vector for the projection
    if(!low_memory_mode)
    {
      inner_product_matrix->vector_mult(*temp, *solution);
    }
    else
    {
      assemble_inner_product_matrix(matrix);
      matrix->vector_mult(*temp, *solution);
    }

    zero_dirichlet_dofs_on_vector(*temp);

    // Do not assume that RB_stiffness matrix is diagonal,
    // diagonality degrades as N increases

    // Get an appropriately sized copy of RB_inner_product_matrix
    DenseMatrix<Number> RB_inner_product_matrix_N(RB_size,RB_size);
    for(unsigned int i=0; i<RB_size; i++)
      for(unsigned int j=0; j<RB_size; j++)
      {
        RB_inner_product_matrix_N(i,j) = RB_inner_product_matrix(i,j);
      }

    // Compute the projection RHS
    DenseVector<Number> RB_proj_rhs(RB_size);
    for(unsigned int i=0; i<RB_size; i++)
    {
      RB_proj_rhs(i) = temp->dot(get_bf(i));
    }

    DenseVector<Number> RB_proj(RB_size);

    // Now solve the linear system
    RB_inner_product_matrix_N.lu_solve(RB_proj_rhs, RB_proj);

    // Load the RB projection into temp
    temp->zero();
    for(unsigned int i=0; i<RB_size; i++)
    {
      temp->add(RB_proj(i), get_bf(i));
    }

    temp->add(-1., *solution);

    // Now temp holds the projection error, store in temporal_data
    *(temporal_data[_k]) = *temp;
  }

  STOP_LOG("set_error_temporal_data()", "TransientRBSystem");

  // return the square of the X norm of the truth solution
  if(!low_memory_mode)
  {
    inner_product_matrix->vector_mult(*inner_product_storage_vector,*solution);
  }
  else
  {
    assemble_inner_product_matrix(matrix);
    matrix->vector_mult(*inner_product_storage_vector,*solution);
  }

  return solution->dot(*inner_product_storage_vector);
}

void TransientRBSystem::initialize_truth ()
{
  if (nonzero_initialization)
  {
    // Use System::read_serialized_data to read the initial condition
    // into this->solution
    Xdr IC_data(init_filename, READ);
    read_serialized_data(IC_data, false);
  }
  else
  {
    // Otherwise zero out the solution as a default
    this->solution->zero();
  }
  this->solution->close();
  this->update();
}

void TransientRBSystem::add_IC_to_RB_space()
{
  START_LOG("add_IC_to_RB_space()", "TransientRBSystem");

  if (get_n_basis_functions() > 0)
  {
    libMesh::out << "Error: Should not call TransientRBSystem::add_IC_to_RB_space() "
                 << "on a system that already contains basis functions." << std::endl;
    libmesh_error();
  }
  if (!nonzero_initialization)
  {
    libMesh::out << "Error: Should not call TransientRBSystem::add_IC_to_RB_space() "
                 << "when nonzero_initialization==false." << std::endl;
    libmesh_error();
  }

  initialize_truth();

  // load the new basis function into the basis_functions vector.
  basis_functions.push_back( NumericVector<Number>::build().release() );
  NumericVector<Number>& current_bf = *(basis_functions[basis_functions.size()-1]);
  current_bf.init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
  current_bf = *solution;

  // We can just set the norm to 1.
  if(!low_memory_mode)
  {
    inner_product_matrix->vector_mult(*inner_product_storage_vector,*solution);
  }
  else
  {
    assemble_inner_product_matrix(matrix);
    matrix->vector_mult(*inner_product_storage_vector,*solution);
  }
  Real current_bf_norm = libmesh_real(std::sqrt( current_bf.dot(*inner_product_storage_vector) ));
  current_bf.scale(1./current_bf_norm);

  unsigned int saved_delta_N = get_delta_N();
  set_delta_N(1);
  update_system();
  set_delta_N(saved_delta_N);

  STOP_LOG("add_IC_to_RB_space()", "TransientRBSystem");
}

void TransientRBSystem::enrich_RB_space()
{
// Need SLEPc to get the POD eigenvalues
#if defined(LIBMESH_HAVE_SLEPC)
  START_LOG("enrich_RB_space()", "TransientRBSystem");

  // With the "method of snapshots", the size of
  // the eigenproblem is determined by the number
  // of time-steps (rather than the number of spatial dofs).
  int eigen_size = temporal_data.size();
  int LDA = eigen_size; // The leading order of correlation_matrix
  std::vector<Number> correlation_matrix(LDA*eigen_size);

  // set values of the correlation matrix
  if(low_memory_mode)
    assemble_inner_product_matrix(matrix);

  for(int i=0; i<eigen_size; i++)
  {
    if(!low_memory_mode)
    {
      inner_product_matrix->vector_mult(*inner_product_storage_vector, *temporal_data[i]);
    }
    else
    {
      matrix->vector_mult(*inner_product_storage_vector, *temporal_data[i]);
    }

    for(int j=i; j<eigen_size; j++)
    {
       // Scale the inner products by the number of time-steps to normalize the
       // POD energy norm appropriately
      Number inner_prod = (temporal_data[j]->dot(*inner_product_storage_vector)) / (Real)(_K+1);

      // Fill upper triangular part of correlation_matrix
      correlation_matrix[j*eigen_size+i] = inner_prod;
    }
  }

  // Call the LAPACK eigensolver
  char JOBZ = 'V'; // Compute eigenvectors and eigenvalues
  char RANGE = 'A'; // Compute all eigenvalues
  char UPLO = 'U'; // Upper triangular symmetric matrix
  
  Real VL = 0.; // Not used when RANGE = A
  Real VU = 0.; // Not used when RANGE = A
  
  int IL = 0; // Not used when RANGE = A
  int IU = 0; // Not used when RANGE = A
  
  double ABSTOL = 1.e-14; // Absolute tolerance for eigensolver
  
  int M = 0; // (output) The total number of evals found
  
  std::vector<Real> W(eigen_size); // (output) the eigenvalues

  int LDZ = eigen_size; // The leading order of Z
  std::vector<Number> Z(LDZ*eigen_size); // (output) the eigenvectors

  std::vector<int> ISUPPZ(2*eigen_size); // Indicates which evecs in Z are nonzero
  
  // Work array, sized according to lapack documentation
  int LWORK = 26*eigen_size;
  std::vector<Number> WORK(LWORK);
  
  // Work array, sized according to lapack documentation
  int LIWORK = 10*eigen_size;
  std::vector<int> IWORK(LIWORK);
  
  int INFO = 0;
  
#ifdef LIBMESH_USE_REAL_NUMBERS // Use real numbers
  // Call the eigensolver for symmetric eigenvalue problems.
  // NOTE: evals in W are in ascending order
  LAPACKsyevr_(&JOBZ, &RANGE, &UPLO, &eigen_size, &correlation_matrix[0],
               &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W[0], &Z[0], &LDZ,
               &ISUPPZ[0], &WORK[0], &LWORK, &IWORK[0], &LIWORK, &INFO );
#endif

#ifdef LIBMESH_USE_COMPLEX_NUMBERS // Use complex numbers
  // Need some extra data in the complex case
  
  // Work array, sized according to lapack documentation
  int LRWORK = 24*eigen_size;
  std::vector<Real> RWORK(LRWORK);

  // Call the eigensolver for symmetric eigenvalue problems.
  // NOTE: evals in W are in ascending order
  LAPACKsyevr_(&JOBZ, &RANGE, &UPLO, &eigen_size, &correlation_matrix[0],
               &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W[0], &Z[0], &LDZ,
               &ISUPPZ[0], &WORK[0], &LWORK, &RWORK[0], &LRWORK, &IWORK[0],
               &LIWORK, &INFO );
#endif

  if (INFO != 0)
  {
    libMesh::out << "Error in LAPACK syev eigensolver routine, INFO = " << INFO << std::endl;
    libmesh_error();
  }

  // eval and evec now hold the sorted eigenvalues/eigenvectors
  libMesh::out << std::endl << "POD Eigenvalues:" << std::endl;
  for(unsigned int i=0; i<=2; i++)
  {
    libMesh::out << "eigenvalue " << i << " = " << W[eigen_size-1-i] << std::endl;
  }
  libMesh::out << "eigenvalue K = " << W[0] << std::endl;
  libMesh::out << std::endl;

  // Now load the new basis functions
  unsigned int count = 0;
  while (true)
    {
      // load the new basis function into the basis_functions vector.
      basis_functions.push_back( NumericVector<Number>::build().release() );
      NumericVector<Number>& current_bf = *(basis_functions[basis_functions.size()-1]);
      current_bf.init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
      current_bf.zero();

      // Perform the matrix multiplication of temporal data with
      // the next POD eigenvector
      for(int j=0; j<eigen_size; j++)
	{
	  int Z_row = j;
	  int Z_col = (eigen_size-1-count);
	  int Z_index = Z_col*eigen_size + Z_row;
	  current_bf.add(Z[Z_index], *temporal_data[j]);
	}

      // We just set the norm to 1.
      if(!low_memory_mode)
	{
	  inner_product_matrix->vector_mult(*inner_product_storage_vector,current_bf);
	}
      else
	{
	  matrix->vector_mult(*inner_product_storage_vector,current_bf);
	}
      Real current_bf_norm = std::abs( std::sqrt( current_bf.dot(*inner_product_storage_vector) ) );
      current_bf.scale(1./current_bf_norm);

      // Increment count here since we use the incremented counter
      // in the if clauses below
      count++;

      // If positive POD_tol, we use it to determine the number of basis functions
      // to add, and then break the loop when POD_tol is satisfied, or after Nmax
      // basis functions have been added. Else we break the loop after delta_N
      // (or Nmax) new basis functions.
      if (POD_tol > 0.)
	{
	  set_delta_N(1);

	  // We need to define the updated RB system matrices before the RB solve
	  update_system();
	  Real error_bound = RB_solve(get_n_basis_functions());

	  if ( (error_bound <= POD_tol) || (get_n_basis_functions()==get_Nmax()) )
	    {
	      set_delta_N(0);
	      break;
	    }
	}
      else
	{
	  if (count == get_delta_N())
	  {
	    break;
	  }
	  else
	  if (get_n_basis_functions()==get_Nmax())
	  {
	    set_delta_N(count);
	    break;
	  }
	}
    }

  STOP_LOG("enrich_RB_space()", "TransientRBSystem");
#else
  libmesh_not_implemented();
#endif
}


void TransientRBSystem::update_system()
{
  // If delta_N is set to zero, there is nothing to update
  if(get_delta_N() == 0)
    return;

  Parent::update_system();

  libMesh::out << "Updating RB initial conditions" << std::endl;
  update_RB_initial_condition_all_N();
}

Real TransientRBSystem::RB_solve(unsigned int N)
{
  START_LOG("RB_solve()", "TransientRBSystem");

  if(N > get_n_basis_functions())
  {
    libMesh::err << "ERROR: N cannot be larger than the number "
                 << "of basis functions in RB_solve" << std::endl;
    libmesh_error();
  }
  if(N==0)
  {
    libMesh::err << "ERROR: N must be greater than 0 in RB_solve" << std::endl;
    libmesh_error();
  }

  // First assemble the mass matrix
  DenseMatrix<Number> RB_mass_matrix_N(N,N);
  RB_mass_matrix_N.zero();
  DenseMatrix<Number> RB_M_q_m;
  for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
  {
    RB_M_q_vector[q_m].get_principal_submatrix(N, RB_M_q_m);
    RB_mass_matrix_N.add(eval_theta_q_m(q_m), RB_M_q_m);
  }

  DenseMatrix<Number> RB_LHS_matrix(N,N);
  RB_LHS_matrix.zero();

  DenseMatrix<Number> RB_RHS_matrix(N,N);
  RB_RHS_matrix.zero();

  RB_LHS_matrix.add(1./dt, RB_mass_matrix_N);
  RB_RHS_matrix.add(1./dt, RB_mass_matrix_N);

  DenseMatrix<Number> RB_A_q_a;
  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    RB_A_q_vector[q_a].get_principal_submatrix(N, RB_A_q_a);

    RB_LHS_matrix.add(       euler_theta*eval_theta_q_a(q_a), RB_A_q_a);
    RB_RHS_matrix.add( -(1.-euler_theta)*eval_theta_q_a(q_a), RB_A_q_a);
  }

  // Set system time level to 0
  error_bound_all_k.resize(_K+1);
  set_time_level(0); // Sets the member variable _k to zero

  // Resize/clear the solution vector
  RB_solution.resize(N);

  // Load the initial condition into RB_solution
  RB_solution = RB_initial_condition_all_N[N-1];

  // Resize/clear the old solution vector
  old_RB_solution.resize(N);

  // Initialize the RB rhs
  DenseVector<Number> RB_rhs(N);
  RB_rhs.zero();

  // Initialize the vectors storing solution data
  RB_temporal_solution_data.resize(_K+1);
  for(unsigned int time_level=0; time_level<=_K; time_level++)
  {
    RB_temporal_solution_data[time_level].resize(N);
  }
  // and load the _k=0 data
  RB_temporal_solution_data[_k] = RB_solution;

  Real error_bound_sum = pow( initial_L2_error_all_N[N-1], 2.);

  // Set error bound at _k=0
  error_bound_all_k[_k] = std::sqrt(error_bound_sum);

  // Compute the outputs and associated error bounds at _k=0
  DenseVector<Number> RB_output_vector_N;
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    RB_outputs_all_k[n][_k] = 0.;
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
      RB_outputs_all_k[n][_k] += eval_theta_q_l(n,q_l)*RB_output_vector_N.dot(RB_solution);
    }

    RB_output_error_bounds_all_k[n][_k] = error_bound_all_k[_k] * eval_output_dual_norm(n);
  }

  Real alpha_LB = get_SCM_lower_bound();
  
  // Precompute time-invariant parts of the dual norm of the residual.
  cache_online_residual_terms(N);

  for(unsigned int time_level=1; time_level<=_K; time_level++)
  {
    set_time_level(time_level); // This updates the member variable _k
    old_RB_solution = RB_solution;

    // Compute RB_rhs, as RB_LHS_matrix x old_RB_solution
    RB_RHS_matrix.vector_mult(RB_rhs, old_RB_solution);

    // Add forcing terms
    DenseVector<Number> RB_F_q_f;
    for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
    {
      RB_F_q_vector[q_f].get_principal_subvector(N, RB_F_q_f);
      RB_rhs.add(eval_theta_q_f(q_f), RB_F_q_f);
    }

    RB_LHS_matrix.lu_solve(RB_rhs, RB_solution);

    // Save RB_solution for current time level
    RB_temporal_solution_data[_k] = RB_solution;

    // Evaluate the dual norm of the residual for RB_solution_vector
//    Real epsilon_N = uncached_compute_residual_dual_norm(N);
    Real epsilon_N = compute_residual_dual_norm(N);

    error_bound_sum += residual_scaling_numer(alpha_LB) * pow(epsilon_N, 2.);

    // store error bound at time-level _k
    error_bound_all_k[_k] = std::sqrt(error_bound_sum/residual_scaling_denom(alpha_LB));

    // Now compute the outputs and associated errors
    DenseVector<Number> RB_output_vector_N;
    for(unsigned int n=0; n<get_n_outputs(); n++)
    {
      RB_outputs_all_k[n][_k] = 0.;
      for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      {
        RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
        RB_outputs_all_k[n][_k] += eval_theta_q_l(n,q_l)*RB_output_vector_N.dot(RB_solution);
      }

      RB_output_error_bounds_all_k[n][_k] = error_bound_all_k[_k] * eval_output_dual_norm(n);
    }
  }

  // Now compute the L2 norm of the RB solution at time-level _K
  // to normalize the error bound
  // We reuse RB_rhs here
  DenseMatrix<Number> RB_L2_matrix_N;
  RB_L2_matrix.get_principal_submatrix(N,RB_L2_matrix_N);
  RB_L2_matrix_N.vector_mult(RB_rhs, RB_solution);
  Real final_RB_L2_norm = libmesh_real(std::sqrt(RB_solution.dot(RB_rhs)));

  STOP_LOG("RB_solve()", "TransientRBSystem");

  return ( return_rel_error_bound ? error_bound_all_k[_K]/final_RB_L2_norm : error_bound_all_k[_K] );
}

Real TransientRBSystem::residual_scaling_numer(Real)
{
  return get_dt();
}

Real TransientRBSystem::residual_scaling_denom(Real alpha_LB)
{
  return alpha_LB;
}

void TransientRBSystem::assemble_matrix_for_output_dual_solves()
{
  // By default we use the L2 matrix for transient problems
  
  if(!low_memory_mode)
  {
    matrix->zero();
    matrix->add(1., *L2_matrix);
  }
  else
  {
    assemble_L2_matrix(matrix);
  }
}

void TransientRBSystem::load_RB_solution()
{
  START_LOG("load_RB_solution()", "TransientRBSystem");

  solution->zero();

  DenseVector<Number> RB_solution_vector_k = RB_temporal_solution_data[_k];

  if(RB_solution_vector_k.size() > basis_functions.size())
  {
    libMesh::err << "ERROR: System contains " << basis_functions.size() << " basis functions."
                 << " RB_solution vector constains " << RB_solution.size() << " entries."
                 << " RB_solution in TransientRBSystem::load_RB_solution is too long!" << std::endl;
    libmesh_error();
  }

  for(unsigned int i=0; i<RB_solution_vector_k.size(); i++)
  {
    solution->add(RB_solution_vector_k(i), *basis_functions[i]);
  }

  update();

  STOP_LOG("load_RB_solution()", "TransientRBSystem");
}

void TransientRBSystem::update_RB_system_matrices()
{
  START_LOG("update_RB_system_matrices()", "TransientRBSystem");

  RBSystem::update_RB_system_matrices();

  unsigned int RB_size = get_n_basis_functions();

  AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build();
  temp->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
  
  for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
  {
    for(unsigned int j=0; j<RB_size; j++)
    {
      Number value = 0.;
      
      // Compute reduced L2 matrix
      temp->zero();
      if(!low_memory_mode)
      {
        L2_matrix->vector_mult(*temp, *basis_functions[j]);
      }
      else
      {
        assemble_L2_matrix(matrix);
        matrix->vector_mult(*temp, *basis_functions[j]);
      }
      value = basis_functions[i]->dot(*temp);
      RB_L2_matrix(i,j) = value;
      if(i!=j)
      {
        // The L2 matrix is assumed
        // to be symmetric
        RB_L2_matrix(j,i) = value;
      }

      for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
      {
        // Compute reduced M_q matrix
        temp->zero();
        if(!low_memory_mode)
        {
          get_M_q(q_m)->vector_mult(*temp, *basis_functions[j]);
        }
        else
        {
          assemble_Mq_matrix(q_m,matrix);
          matrix->vector_mult(*temp, *basis_functions[j]);
        }

        value = (*basis_functions[i]).dot(*temp);
        RB_M_q_vector[q_m](i,j) = value;

        if(i!=j)
        {
          // Each mass matrix term is assumed
          // to be symmetric
          RB_M_q_vector[q_m](j,i) = value;
        }
      }
      
    }
  }

  STOP_LOG("update_RB_system_matrices()", "TransientRBSystem");
}




void TransientRBSystem::update_residual_terms(bool compute_inner_products)
{
  START_LOG("update_residual_terms()", "TransientRBSystem");

  RBSystem::update_residual_terms(compute_inner_products);

  unsigned int RB_size = get_n_basis_functions();

  if(!low_memory_mode)
  {
    matrix->zero();
    matrix->add(1., *inner_product_matrix);
    if(constrained_problem)
      matrix->add(1., *constraint_matrix);
  }
  else
  {
    assemble_inner_product_matrix(matrix);
    if(constrained_problem)
      add_scaled_matrix_and_vector(1., constraint_assembly, NULL, matrix, NULL);
  }

  if(reuse_preconditioner)
  {
    // For the first solve, make sure we generate a new preconditioner
    // Actually, this is not necessary as we can use the preconditioner
    // generated by RBSystem::update_residual_terms()
    linear_solver->same_preconditioner = true;
  }

  for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
  {
    for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
    {
      // Initialize the vectors when we need them
      M_q_representor[q_m][i] = (NumericVector<Number>::build().release());
      M_q_representor[q_m][i]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

      rhs->zero();
      if(!low_memory_mode)
      {
        M_q_vector[q_m]->vector_mult(*rhs, *basis_functions[i]);
      }
      else
      {
        assemble_scaled_matvec(1.,
                               M_q_intrr_assembly_vector[q_m],
                               M_q_bndry_assembly_vector[q_m],
                               *rhs,
                               *basis_functions[i]);
      }
      zero_dirichlet_dofs_on_rhs();

      solution->zero();

      if (!quiet)
        libMesh::out << "Starting solve i="
                     << i << " in TransientRBSystem::update_residual_terms() at "
                     << Utility::get_timestamp() << std::endl;
      solve();

      if (!quiet)
        {
          libMesh::out << "Finished solve i="
                       << i << " in TransientRBSystem::update_residual_terms() at "
                       << Utility::get_timestamp() << std::endl;

          libMesh::out << this->n_linear_iterations()
                       << " iterations, final residual "
                       << this->final_linear_residual() << std::endl;
        }

      // Make sure we didn't max out the number of iterations
      if( (this->n_linear_iterations() >=
          this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
          (this->final_linear_residual() >
          this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
      {
        libMesh::out << "Warning: Linear solver may not have converged! Final linear residual = "
                     << this->final_linear_residual() << ", number of iterations = "
                     << this->n_linear_iterations() << std::endl << std::endl;
  //       libmesh_error();
      }

      *M_q_representor[q_m][i] = *solution;

      if(reuse_preconditioner)
      {
        linear_solver->same_preconditioner = true;
      }
    }
  }

  if(reuse_preconditioner)
  {
    linear_solver->same_preconditioner = false;
  }


  // Now compute and store the inner products if requested
  if (compute_inner_products)
    {
      if(low_memory_mode && constrained_problem)
	assemble_inner_product_matrix(matrix);

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
	{
	  if(!low_memory_mode)
	    {
	      inner_product_matrix->vector_mult(*inner_product_storage_vector, *F_q_representor[q_f]);
	    }
	  else
	    {
	      matrix->vector_mult(*inner_product_storage_vector, *F_q_representor[q_f]);
	    }

	  for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
	    {
	      for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
		{
		  Fq_Mq_representor_norms[q_f][q_m][i] =
		    M_q_representor[q_m][i]->dot(*inner_product_storage_vector);
		} // end for q_m
	    } // end for i
	} // end for q_f

      unsigned int q=0;
      for(unsigned int q_m1=0; q_m1<get_Q_m(); q_m1++)
	{
	  for(unsigned int q_m2=q_m1; q_m2<get_Q_m(); q_m2++)
	    {
	      for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
		{
		  for(unsigned int j=0; j<RB_size; j++)
		    {
		      if(!low_memory_mode)
			{
			  inner_product_matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m2][j]);
			}
		      else
			{
			  matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m2][j]);
			}
		      Mq_Mq_representor_norms[q][i][j] = M_q_representor[q_m1][i]->dot(*inner_product_storage_vector);

		      if(i != j)
			{
			  if(!low_memory_mode)
			    {
			      inner_product_matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m2][i]);
			    }
			  else
			    {
			      matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m2][i]);
			    }
			  Mq_Mq_representor_norms[q][j][i] = M_q_representor[q_m1][j]->dot(*inner_product_storage_vector);
			}
		    } // end for j
		} // end for i
	      q++;
	    } // end for q_m2
	} // end for q_m1


      for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
	{
	  for(unsigned int j=0; j<RB_size; j++)
	    {
	      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
		{
		  for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
		    {
		      if(!low_memory_mode)
			{
			  inner_product_matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m][j]);
			}
		      else
			{
			  matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m][j]);
			}
		      Aq_Mq_representor_norms[q_a][q_m][i][j] =
			A_q_representor[q_a][i]->dot(*inner_product_storage_vector);

		      if(i != j)
			{
			  if(!low_memory_mode)
			    {
			      inner_product_matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m][i]);
			    }
			  else
			    {
			      matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m][i]);
			    }
			  Aq_Mq_representor_norms[q_a][q_m][j][i] =
			    A_q_representor[q_a][j]->dot(*inner_product_storage_vector);
			}
		    } // end for q_m
		} // end for q_a
	    } // end for j
	} // end for i
    } // end if (compute_inner_products)

  STOP_LOG("update_residual_terms()", "TransientRBSystem");
}


void TransientRBSystem::update_RB_initial_condition_all_N()
{
  START_LOG("update_RB_initial_condition_all_N()", "TransientRBSystem");

  // Load the initial condition into the solution vector
  initialize_truth();

  AutoPtr< NumericVector<Number> > temp1 = NumericVector<Number>::build();
  temp1->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

  AutoPtr< NumericVector<Number> > temp2 = NumericVector<Number>::build();
  temp2->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);


  unsigned int RB_size = basis_functions.size();

  // First compute the right-hand side vector for the L2 projection
  if(!low_memory_mode)
  {
    L2_matrix->vector_mult(*temp1, *solution);
  }
  else
  {
    assemble_L2_matrix(matrix);
    matrix->vector_mult(*temp1, *solution);
  }

  for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
  {
    RB_ic_proj_rhs_all_N(i) = temp1->dot(*basis_functions[i]);
  }


  // Now compute the projection for each N
  DenseMatrix<Number> RB_L2_matrix_N;
  DenseVector<Number> RB_rhs_N;
  for(unsigned int N=(RB_size-delta_N); N<RB_size; N++)
  {
    // We have to index here by N+1 since the loop index is zero-based.
    RB_L2_matrix.get_principal_submatrix(N+1, RB_L2_matrix_N);
    
    RB_ic_proj_rhs_all_N.get_principal_subvector(N+1, RB_rhs_N);

    DenseVector<Number> RB_ic_N(N+1);

    // Now solve the linear system
    RB_L2_matrix_N.lu_solve(RB_rhs_N, RB_ic_N);

    // Load RB_ic_N into RB_initial_condition_all_N
    RB_initial_condition_all_N[N] = RB_ic_N;

    // Compute the L2 error for the RB initial condition
    // This part is dependent on the truth space.

    // load the RB solution into temp1
    temp1->zero();
    for(unsigned int i=0; i<N+1; i++)
    {
      temp1->add(RB_ic_N(i), *basis_functions[i]);
    }

    // subtract truth initial condition from RB_ic_N
    temp1->add(-1., *solution);

    // Compute L2 norm error, i.e. sqrt(M(solution,solution))
    temp2->zero();
    if(!low_memory_mode)
    {
      L2_matrix->vector_mult(*temp2, *temp1);
    }
    else
    {
      matrix->vector_mult(*temp2, *temp1);
    }
    initial_L2_error_all_N[N] = libmesh_real(std::sqrt(temp2->dot(*temp1)));
  }

  STOP_LOG("update_RB_initial_condition_all_N()", "TransientRBSystem");
}

void TransientRBSystem::cache_online_residual_terms(const unsigned int N)
{
  START_LOG("cache_online_residual_terms()", "TransientRBSystem");

  cached_Fq_term = 0.;
  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<get_Q_f(); q_f1++)
  {
    Number cached_theta_q_f1 = eval_theta_q_f(q_f1);
    for(unsigned int q_f2=q_f1; q_f2<get_Q_f(); q_f2++)
    {
      Real delta = (q_f1==q_f2) ? 1. : 2.;
      cached_Fq_term += delta*cached_theta_q_f1*eval_theta_q_f(q_f2) * Fq_representor_norms[q];

      q++;
    }
  }

  cached_Fq_Aq_vector.resize(N);
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    Number cached_theta_q_f = eval_theta_q_f(q_f);
    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    {
      Number cached_theta_q_a = eval_theta_q_a(q_a);
      for(unsigned int i=0; i<N; i++)
      {
        cached_Fq_Aq_vector(i) += 2.*cached_theta_q_f*cached_theta_q_a*
                                  Fq_Aq_representor_norms[q_f][q_a][i];
      }
    }
  }

  cached_Aq_Aq_matrix.resize(N,N);
  q=0;
  for(unsigned int q_a1=0; q_a1<get_Q_a(); q_a1++)
  {
    Number cached_theta_q_a1 = eval_theta_q_a(q_a1);
    for(unsigned int q_a2=q_a1; q_a2<get_Q_a(); q_a2++)
    {
      Number cached_theta_q_a2 = eval_theta_q_a(q_a2);
      Real delta = (q_a1==q_a2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          cached_Aq_Aq_matrix(i,j) += delta*
                                      cached_theta_q_a1*cached_theta_q_a2*
                                      Aq_Aq_representor_norms[q][i][j];
        }
      }
      q++;
    }
  }

  cached_Fq_Mq_vector.resize(N);
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    Number cached_theta_q_f = eval_theta_q_f(q_f);
    for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
    {
      Number cached_theta_q_m = eval_theta_q_m(q_m);
      for(unsigned int i=0; i<N; i++)
      {
        cached_Fq_Mq_vector(i) += 2.*cached_theta_q_f * cached_theta_q_m * Fq_Mq_representor_norms[q_f][q_m][i];
      }
    }
  }

  cached_Aq_Mq_matrix.resize(N,N);
  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    Number cached_theta_q_a = eval_theta_q_a(q_a);
    
    for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
    {
      Number cached_theta_q_m = eval_theta_q_m(q_m);
      
      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          cached_Aq_Mq_matrix(i,j) += 2.*cached_theta_q_a*cached_theta_q_m*Aq_Mq_representor_norms[q_a][q_m][i][j];
        }
      }
    }
  }
  
  cached_Mq_Mq_matrix.resize(N,N);
  q=0;
  for(unsigned int q_m1=0; q_m1<get_Q_m(); q_m1++)
  {
    Number cached_theta_q_m1 = eval_theta_q_m(q_m1);
    for(unsigned int q_m2=q_m1; q_m2<get_Q_m(); q_m2++)
    {
      Number cached_theta_q_m2 = eval_theta_q_m(q_m2);
      Real delta = (q_m1==q_m2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          cached_Mq_Mq_matrix(i,j) += delta*
                                      cached_theta_q_m1*cached_theta_q_m2*
                                      Mq_Mq_representor_norms[q][i][j];
        }
      }
      q++;
    }
  }
  
  STOP_LOG("cache_online_residual_terms()", "TransientRBSystem");
}

Real TransientRBSystem::compute_residual_dual_norm(const unsigned int N)
{
  START_LOG("compute_residual_dual_norm()", "TransientRBSystem");

  // This assembly assumes we have already called cache_online_residual_terms
  // and that the RB_solve parameter is constant in time

  const Real dt = get_dt();

  DenseVector<Number> RB_u_euler_theta(N);
  DenseVector<Number> mass_coeffs(N);
  for(unsigned int i=0; i<N; i++)
  {
    RB_u_euler_theta(i)  = euler_theta*RB_solution(i) + (1.-euler_theta)*old_RB_solution(i);
    mass_coeffs(i) = -(RB_solution(i) - old_RB_solution(i))/dt;
  }

  Number residual_norm_sq = cached_Fq_term;
  
  residual_norm_sq += RB_u_euler_theta.dot(cached_Fq_Aq_vector);
  residual_norm_sq += mass_coeffs.dot(cached_Fq_Mq_vector);

  for(unsigned int i=0; i<N; i++)
    for(unsigned int j=0; j<N; j++)
    {
      residual_norm_sq += RB_u_euler_theta(i)*RB_u_euler_theta(j)*cached_Aq_Aq_matrix(i,j);
      residual_norm_sq += mass_coeffs(i)*mass_coeffs(j)*cached_Mq_Mq_matrix(i,j);
      residual_norm_sq += RB_u_euler_theta(i)*mass_coeffs(j)*cached_Aq_Mq_matrix(i,j);
    }


  if(libmesh_real(residual_norm_sq) < 0)
  {
    libMesh::out << "Warning: Square of residual norm is negative "
                 << "in TransientRBSystem::compute_residual_dual_norm()" << std::endl;

    // Sometimes this is negative due to rounding error,
    // but error is on the order of 1.e-10, so shouldn't
    // affect result
//    libmesh_error();
     residual_norm_sq = std::abs(residual_norm_sq);
  }

  STOP_LOG("compute_residual_dual_norm()", "TransientRBSystem");

  return libmesh_real(std::sqrt( residual_norm_sq ));
}


Real TransientRBSystem::uncached_compute_residual_dual_norm(const unsigned int N)
{
  START_LOG("uncached_compute_residual_dual_norm()", "TransientRBSystem");

//   // This is the "slow" way of computing the residual, but it is useful
//   // for validating the "fast" way.
//   // Note that this only works in serial since otherwise each processor will
//   // have a different parameter value during the Greedy training.
//
//   // Assemble the right-hand side to find the Reisz representor
//   // of the residual in the X norm
//   AutoPtr< NumericVector<Number> > RB_sol = NumericVector<Number>::build();
//   RB_sol->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
//   RB_sol->zero();
//
//   AutoPtr< NumericVector<Number> > ghosted_temp = NumericVector<Number>::build();
//   ghosted_temp->init (this->n_dofs(), this->n_local_dofs(),
//                       this->get_dof_map().get_send_list(), false,
//                       GHOSTED);
//
//   AutoPtr< NumericVector<Number> > parallel_temp = NumericVector<Number>::build();
//   parallel_temp->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
//
//   // Store current_local_solution, since we don't want to corrupt it
//   *ghosted_temp = *current_local_solution;
//
//   for(unsigned int i=0; i<N; i++)
//   {
//     RB_sol->add(RB_solution(i), *basis_functions[i]);
//     parallel_temp->add(old_RB_solution(i), *basis_functions[i]);
//   }
//
//   // Load parallel_temp into current_local_solution in order to do assembly
//   const std::vector<unsigned int>& send_list = this->get_dof_map().get_send_list();
//   parallel_temp->localize (*current_local_solution, send_list);
//
//   // Load the system_matrix
//   this->truth_assembly();
//
//   // Restore current_local_solution
//   *current_local_solution = *ghosted_temp;
//
//   matrix->vector_mult(*parallel_temp, *RB_sol);
//   rhs->add(-1., *parallel_temp);
//   rhs->close();
//
//   zero_dirichlet_dofs_on_rhs();
//
//   // Then solve the system to get the Reisz representor
//   matrix->zero();
//   if(!low_memory_mode)
//   {
//     matrix->add(1., *inner_product_matrix);
//     if(constrained_problem)
//       matrix->add(1., *constraint_matrix);
//   }
//   else
//   {
//     assemble_inner_product_matrix(matrix);
//     if(constrained_problem)
//       add_scaled_matrix_and_vector(1., constraint_assembly, NULL, matrix, NULL);
//   }
//
//
//   solution->zero();
//   solve();
//   // Make sure we didn't max out the number of iterations
//   if( (this->n_linear_iterations() >=
//        this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
//       (this->final_linear_residual() >
//        this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
//   {
//     libMesh::out << "Warning: Linear solver may not have converged! Final linear residual = "
//                  << this->final_linear_residual() << ", number of iterations = "
//                  << this->n_linear_iterations() << std::endl << std::endl;
// //     libmesh_error();
//   }
//
//  if(!low_memory_mode)
//  {
//    inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
//  }
//  else
//  {
//    assemble_inner_product_matrix(matrix);
//    matrix->vector_mult(*inner_product_storage_vector, *solution);
//  }
//
//  Real slow_residual_norm_sq = solution->dot(*inner_product_storage_vector);





  // Use the stored representor inner product values
  // to evaluate the residual norm

  const Real dt = get_dt();

  std::vector<Number> RB_u_euler_theta(N);
  std::vector<Number> mass_coeffs(N);
  for(unsigned int i=0; i<N; i++)
  {
    RB_u_euler_theta[i]  = euler_theta*RB_solution(i) + (1.-euler_theta)*old_RB_solution(i);
    mass_coeffs[i] = -(RB_solution(i) - old_RB_solution(i))/dt;
  }

  Number residual_norm_sq = 0.;

  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<get_Q_f(); q_f1++)
  {
    Number cached_theta_q_f1 = eval_theta_q_f(q_f1);
    for(unsigned int q_f2=q_f1; q_f2<get_Q_f(); q_f2++)
    {
      Real delta = (q_f1==q_f2) ? 1. : 2.;
      residual_norm_sq += delta*cached_theta_q_f1*eval_theta_q_f(q_f2) * Fq_representor_norms[q];

      q++;
    }
  }

  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    Number cached_theta_q_f = eval_theta_q_f(q_f);
    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    {
      Number cached_theta_q_a = eval_theta_q_a(q_a);
      for(unsigned int i=0; i<N; i++)
      {
        residual_norm_sq += 2.*RB_u_euler_theta[i]*cached_theta_q_f*cached_theta_q_a*
                               Fq_Aq_representor_norms[q_f][q_a][i];
      }
    }
  }

  q=0;
  for(unsigned int q_a1=0; q_a1<get_Q_a(); q_a1++)
  {
    Number cached_theta_q_a1 = eval_theta_q_a(q_a1);
    for(unsigned int q_a2=q_a1; q_a2<get_Q_a(); q_a2++)
    {
      Number cached_theta_q_a2 = eval_theta_q_a(q_a2);
      Real delta = (q_a1==q_a2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          residual_norm_sq += delta*RB_u_euler_theta[i]*RB_u_euler_theta[j]*
                              cached_theta_q_a1*cached_theta_q_a2*
                              Aq_Aq_representor_norms[q][i][j];
        }
      }
      q++;
    }
  }

  // Now add the terms due to the time-derivative
  q=0;
  for(unsigned int q_m1=0; q_m1<get_Q_m(); q_m1++)
  {
    Number cached_theta_q_m1 = eval_theta_q_m(q_m1);
    for(unsigned int q_m2=q_m1; q_m2<get_Q_m(); q_m2++)
    {
      Number cached_theta_q_m2 = eval_theta_q_m(q_m2);
      Real delta = (q_m1==q_m2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          residual_norm_sq += delta*mass_coeffs[i]*mass_coeffs[j]*
                              cached_theta_q_m1*cached_theta_q_m2*
                              Mq_Mq_representor_norms[q][i][j];
        }
      }
      q++;
    }
  }

  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    Number cached_theta_q_f = eval_theta_q_f(q_f);
    for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
    {
      Number cached_theta_q_m = eval_theta_q_m(q_m);
      for(unsigned int i=0; i<N; i++)
      {
        residual_norm_sq += 2.*mass_coeffs[i]*cached_theta_q_f * cached_theta_q_m * Fq_Mq_representor_norms[q_f][q_m][i];
      }
    }
  }

  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    Number cached_theta_q_a = eval_theta_q_a(q_a);
    
    for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
    {
      Number cached_theta_q_m = eval_theta_q_m(q_m);
      
      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          residual_norm_sq += 2.*RB_u_euler_theta[i]*mass_coeffs[j]*
                                         cached_theta_q_a*cached_theta_q_m*
                                         Aq_Mq_representor_norms[q_a][q_m][i][j];
        }
      }
    }
  }
  
  if(libmesh_real(residual_norm_sq) < 0)
  {
    libMesh::out << "Warning: Square of residual norm is negative "
                 << "in TransientRBSystem::compute_residual_dual_norm()" << std::endl;

    // Sometimes this is negative due to rounding error,
    // but error is on the order of 1.e-10, so shouldn't
    // affect result
//    libmesh_error();
     residual_norm_sq = std::abs(residual_norm_sq);
  }

//   libMesh::out << "slow residual_sq = " << slow_residual_norm_sq
//                << ", fast residual_sq = " << residual_norm_sq << std::endl;

  STOP_LOG("uncached_compute_residual_dual_norm()", "TransientRBSystem");

  return libmesh_real(std::sqrt( residual_norm_sq ));
}


void TransientRBSystem::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "TransientRBSystem");

  Parent::write_offline_data_to_files(directory_name);

  const unsigned int n_bfs = get_n_basis_functions();
  libmesh_assert( n_bfs <= Nmax );

  const unsigned int precision_level = 14;

  if(libMesh::processor_id() == 0)
  {
    // Write out the L2 matrix
    std::ofstream RB_L2_matrix_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_L2_matrix.dat";
      RB_L2_matrix_out.open(file_name.str().c_str());
    }
    RB_L2_matrix_out.precision(precision_level);
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<n_bfs; j++)
      {
        RB_L2_matrix_out << std::scientific << RB_L2_matrix(i,j) << " ";
      }
    }
    RB_L2_matrix_out.close();
          
    // Write out the M_q matrices
    for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_M_";
      OSSRealzeroright(file_name,3,0,q_m);
      file_name << ".dat";
      std::ofstream RB_M_q_m_out(file_name.str().c_str());

      if ( !RB_M_q_m_out.good() )
      {
        libMesh::err << "Error opening RB_M_" << q_m << ".dat" << std::endl;
        libmesh_error();
      }

      RB_M_q_m_out.precision(precision_level);
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int j=0; j<n_bfs; j++)
        {
          RB_M_q_m_out << std::scientific << RB_M_q_vector[q_m](i,j) << " ";
        }
      }
      RB_M_q_m_out.close();
    }

    // Write out the initial condition data
    // and the initial L2 error for all N
    std::ofstream initial_conditions_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/initial_conditions.dat";
      initial_conditions_out.open(file_name.str().c_str());
    }
    std::ofstream initial_L2_error_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/initial_L2_error.dat";
      initial_L2_error_out.open(file_name.str().c_str());
    }
    if (!initial_conditions_out.good() || !initial_L2_error_out.good())
    {
      libMesh::err << "Error opening initial conditions output files" << std::endl;
      libmesh_error();
    }

    initial_conditions_out.precision(precision_level);
    initial_L2_error_out.precision(precision_level);

    for(unsigned int i=0; i<n_bfs; i++)
    {
      initial_L2_error_out << initial_L2_error_all_N[i] << " ";
      for(unsigned int j=0; j<=i; j++)
      {
        initial_conditions_out << std::scientific << RB_initial_condition_all_N[i](j) << " ";
      }
    }
    initial_conditions_out.close();
    initial_L2_error_out.close();

    // Next write out the Fq_Mq representor norm data
    std::ofstream RB_Fq_Mq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Fq_Mq_norms.dat";
      RB_Fq_Mq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Fq_Mq_norms_out.good() )
    {
      libMesh::err << "Error opening Fq_Mq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Fq_Mq_norms_out.precision(precision_level);
    for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
    {
      for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
      {
        for(unsigned int i=0; i<n_bfs; i++)
        {
          RB_Fq_Mq_norms_out << std::scientific << Fq_Mq_representor_norms[q_f][q_m][i] << " ";
        }
      }
    }
    RB_Fq_Mq_norms_out.close();

    // Next write out the Mq_Mq representor norm data
    std::ofstream RB_Mq_Mq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Mq_Mq_norms.dat";
      RB_Mq_Mq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Mq_Mq_norms_out.good() )
    {
      libMesh::err << "Error opening RB_Mq_Mq_norms_out.dat" << std::endl;
      libmesh_error();
    }
    RB_Mq_Mq_norms_out.precision(precision_level);
    unsigned int Q_m_hat = get_Q_m()*(get_Q_m()+1)/2;
    for(unsigned int q=0; q<Q_m_hat; q++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int j=0; j<n_bfs; j++)
        {
          RB_Mq_Mq_norms_out << std::scientific << Mq_Mq_representor_norms[q][i][j] << " ";
        }
      }
    }
    RB_Mq_Mq_norms_out.close();

    // Next write out the Aq_Mq representor norm data
    std::ofstream RB_Aq_Mq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Aq_Mq_norms.dat";
      RB_Aq_Mq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Aq_Mq_norms_out.good() )
    {
      libMesh::err << "Error opening Aq_Mq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Aq_Mq_norms_out.precision(precision_level);
    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    {
      for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
      {
        for(unsigned int i=0; i<n_bfs; i++)
        {
          for(unsigned int j=0; j<n_bfs; j++)
          {
            RB_Aq_Mq_norms_out << std::scientific << Aq_Mq_representor_norms[q_a][q_m][i][j] << " ";
          }
        }
      }
    }
    RB_Aq_Mq_norms_out.close();
  }

  // Write out the residual representors to file if requested
  if (store_representors)
    {
      // Write out the M_q_representors.  These are useful to have when restarting,
      // so you don't have to recompute them all over again.  There should be
      // this->get_n_basis_functions() of these.
      libMesh::out << "Writing out the M_q_representors..." << std::endl;

      std::ostringstream file_name;
      const std::string residual_representor_suffix = (write_binary_residual_representors ? ".xdr" : ".dat");

      // Residual representors written out to their own separate directory
      std::string residual_representors_dir = "residual_representors";

      const unsigned int istop  = this->get_n_basis_functions();
      const unsigned int istart = istop-delta_N;

      for (unsigned int q=0; q<M_q_representor.size(); ++q)
        for (unsigned int i=istart; i<istop; ++i)
	{
	  libMesh::out << "Writing out M_q_representor[" << q << "][" << i << "]..." << std::endl;
	  libmesh_assert(M_q_representor[q][i] != NULL);

	  file_name.str(""); // reset filename
	  file_name << residual_representors_dir << "/M_q_representor" << i << residual_representor_suffix;

	  {
	    // No need to copy!
	    //*solution = *(M_q_representor[q][i]);
	    M_q_representor[q][i]->swap(*solution);

	    Xdr mr_data(file_name.str(),
			write_binary_residual_representors ? ENCODE : WRITE);

	    write_serialized_data(mr_data, false);

	    // Synchronize before moving on
	    Parallel::barrier();

	    // Swap back.
	    M_q_representor[q][i]->swap(*solution);

	    // TODO: bzip the resulting file?  See $LIBMESH_DIR/src/mesh/unstructured_mesh.C
	    // for the system call, be sure to do it only on one processor, etc.
	  }
	}
    } // end if store_representors

  STOP_LOG("write_offline_data_to_files()", "TransientRBSystem");
}




void TransientRBSystem::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "TransientRBSystem");

  Parent::read_offline_data_from_files(directory_name);

  // First, find out how many basis functions we had when Greedy terminated
  // This was set in RBSystem::read_offline_data_from_files
  unsigned int n_bfs = this->get_n_basis_functions();
  
  std::ifstream RB_L2_matrix_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/RB_L2_matrix.dat";
    RB_L2_matrix_in.open(file_name.str().c_str());
  }
  for(unsigned int i=0; i<n_bfs; i++)
  {
    for(unsigned int j=0; j<n_bfs; j++)
    {
      Number value;
      RB_L2_matrix_in >> value;
      RB_L2_matrix(i,j) = value;
    }
  }
  RB_L2_matrix_in.close();

  // Read in the M_q matrices
  for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
  {
    OStringStream file_name;
    file_name << directory_name << "/RB_M_";
    OSSRealzeroright(file_name,3,0,q_m);
    file_name << ".dat";
    std::ifstream RB_M_q_m_in(file_name.str().c_str());

    if ( !RB_M_q_m_in.good() )
    {
      libMesh::err << "Error opening RB_M_" << q_m << ".dat" << std::endl;
      libmesh_error();
    }

    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<n_bfs; j++)
      {
        Number value;
        RB_M_q_m_in >> value;
        RB_M_q_vector[q_m](i,j) = value;
      }
    }
    RB_M_q_m_in.close();
  }


  // Read in the initial condition data
  // and the initial L2 error for all N
  std::ifstream initial_conditions_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/initial_conditions.dat";
    initial_conditions_in.open(file_name.str().c_str());
  }
  std::ifstream initial_L2_error_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/initial_L2_error.dat";
    initial_L2_error_in.open(file_name.str().c_str());
  }
  if (!initial_conditions_in.good() || !initial_L2_error_in.good())
  {
    libMesh::err << "Error opening initial conditions output files" << std::endl;
    libmesh_error();
  }

  for(unsigned int i=0; i<n_bfs; i++)
  {
    initial_L2_error_in >> initial_L2_error_all_N[i];
    for(unsigned int j=0; j<=i; j++)
    {
      initial_conditions_in >> RB_initial_condition_all_N[i](j);
    }
  }
  initial_conditions_in.close();
  initial_L2_error_in.close();

  // Next read in the Fq_Mq representor norm data
  std::ifstream RB_Fq_Mq_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Fq_Mq_norms.dat";
    RB_Fq_Mq_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Fq_Mq_norms_in.good() )
  {
    libMesh::err << "Error opening Fq_Mq_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        RB_Fq_Mq_norms_in >> Fq_Mq_representor_norms[q_f][q_m][i];
      }
    }
  }
  RB_Fq_Mq_norms_in.close();

  // Next read in the Mq_Mq representor norm data
  std::ifstream RB_Mq_Mq_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Mq_Mq_norms.dat";
    RB_Mq_Mq_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Mq_Mq_norms_in.good() )
  {
    libMesh::err << "Error opening RB_Mq_Mq_norms_in.dat" << std::endl;
    libmesh_error();
  }
  unsigned int Q_m_hat = get_Q_m()*(get_Q_m()+1)/2;
  for(unsigned int q=0; q<Q_m_hat; q++)
  {
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<n_bfs; j++)
      {
        RB_Mq_Mq_norms_in >> Mq_Mq_representor_norms[q][i][j];
      }
    }
  }
  RB_Mq_Mq_norms_in.close();

  // Next read in the Aq_Mq representor norm data
  std::ifstream RB_Aq_Mq_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Aq_Mq_norms.dat";
    RB_Aq_Mq_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Aq_Mq_norms_in.good() )
  {
    libMesh::err << "Error opening Aq_Mq_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int j=0; j<n_bfs; j++)
        {
          RB_Aq_Mq_norms_in >> Aq_Mq_representor_norms[q_a][q_m][i][j];
        }
      }
    }
  }
  RB_Aq_Mq_norms_in.close();

  // Read in the representors if requested
  if (store_representors)
    {
      const std::string residual_representors_dir = "residual_representors";
      const std::string residual_representor_suffix =
	(read_binary_residual_representors ? ".xdr" : ".dat");

      std::ostringstream file_name;
      struct stat stat_info;
      

      libMesh::out << "Reading in the M_q_representors..." << std::endl;

      // Read in the A_q representors.  The class makes room for [Q_m][Nmax] of these.  We are going to
      // read in [Q_m][this->get_n_basis_functions()].  FIXME:
      // should we be worried about leaks in the locations where we're about to fill entries?
      for (unsigned int i=0; i<M_q_representor.size(); ++i)
	for (unsigned int j=0; j<M_q_representor[i].size(); ++j)
	  {
	    if (M_q_representor[i][j] != NULL)
	      {
		libMesh::out << "Error, must delete existing M_q_representor before reading in from file."
			     << std::endl;
		libmesh_error();
	      }
	  }

      // Now ready to read them in from file!
      for (unsigned int i=0; i<M_q_representor.size(); ++i)
	for (unsigned int j=0; j<this->get_n_basis_functions(); ++j)
	  {
	    file_name.str(""); // reset filename
	    file_name << residual_representors_dir
		      << "/M_q_representor" << i << "_" << j << residual_representor_suffix;

	    // On processor zero check to be sure the file exists
	    if (libMesh::processor_id() == 0)
	    {
	      int stat_result = stat(file_name.str().c_str(), &stat_info);

	      if (stat_result != 0)
		{
		  libMesh::out << "File does not exist: " << file_name.str() << std::endl;
		  libmesh_error();
		}
	    }

	    Xdr aqr_data(file_name.str(),
			 read_binary_residual_representors ? DECODE : READ);

	    read_serialized_data(aqr_data, false);

	    M_q_representor[i][j] = NumericVector<Number>::build().release();
	    M_q_representor[i][j]->init (this->n_dofs(), this->n_local_dofs(),
					 false, libMeshEnums::PARALLEL);

	    // No need to copy, just swap
	    //*M_q_representor[i][j] = *solution;
	    M_q_representor[i][j]->swap(*solution);
	  }
    } // end if (store_representors)

  STOP_LOG("read_offline_data_from_files()", "TransientRBSystem");
}

} // namespace libMesh
