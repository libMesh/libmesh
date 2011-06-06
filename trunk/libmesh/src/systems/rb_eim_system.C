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
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dof_map.h"
#include "libmesh_logging.h"
#include "equation_systems.h"
#include "parallel.h"
#include "parallel_algebra.h"
#include "fe.h"
#include "quadrature.h"
#include "utility.h"
#include "fe_interface.h"
#include "fe_compute_data.h"
#include "getpot.h"
#include <fstream>
#include <sstream>
#include "o_string_stream.h"
#include "gmv_io.h"

#include "fem_context.h"
#include "rb_eim_system.h"
#include "rb_eim_evaluation.h"

namespace libMesh
{

RBEIMSystem::RBEIMSystem (EquationSystems& es,
		          const std::string& name,
		          const unsigned int number)
  : Parent(es, name, number),
    best_fit_type_flag(PROJECTION_BEST_FIT),
    mesh_function(NULL),
    performing_extra_greedy_step(false),
    current_bf_index(0)
{
  use_empty_RB_solve_in_greedy = false;
}

RBEIMSystem::~RBEIMSystem ()
{
  delete mesh_function;
  mesh_function = NULL;
}

std::string RBEIMSystem::system_type () const
{
  return "RBEIMSystem";
}

void RBEIMSystem::process_parameters_file (const std::string& parameters_filename)
{
  // Indicate that we need to compute the RB
  // inner product matrix in this case
  compute_RB_inner_product = true;
  
  // Indicate that we need the training set
  // for the Greedy to be the same on all
  // processors
  serial_training_set = true;
  
  Parent::process_parameters_file(parameters_filename);
  
  if(n_vars() != get_n_parametrized_functions() )
  {
    libMesh::out << "Error: The number of parametrized_functions must match n_vars in RBEIMSystem."
                 << std::endl;
    libmesh_error();
  }
  
  GetPot infile(parameters_filename);
  
  std::string best_fit_type_string = infile("best_fit_type","projection");
  
  if(best_fit_type_string == "projection")
  {
    best_fit_type_flag = PROJECTION_BEST_FIT;
  }
  else
  if(best_fit_type_string == "eim")
  {
    best_fit_type_flag = EIM_BEST_FIT;
  }
  else
  {
    libMesh::out << "Error: invalid best_fit_type in input file" << std::endl;
    libmesh_error();
  }
  
  libMesh::out << std::endl << "RBEIMSystem parameters:" << std::endl;
  libMesh::out << "best fit type: " << best_fit_type_string << std::endl;
  libMesh::out << "number of parametrized functions: " << get_n_parametrized_functions() << std::endl;
  libMesh::out << std::endl;
}

void RBEIMSystem::initialize_RB_system(RBEvaluation* rb_evaluation_in)
{
  Parent::initialize_RB_system(rb_evaluation_in);

  // initialize a serial vector that we will use for MeshFunction evaluations
  serialized_vector = NumericVector<Number>::build();
  serialized_vector->init (this->n_dofs(), false, SERIAL);

  // Initialize the MeshFunction for interpolating the
  // solution vector at quadrature points
  std::vector<unsigned int> vars(n_vars());
  Utility::iota(vars.begin(), vars.end(), 0); // By default use all variables
  mesh_function = new MeshFunction(get_equation_systems(), *serialized_vector, get_dof_map(), vars);
  mesh_function->init();
    
  // initialize the vector that stores the _current_ basis function,
  // i.e. the vector that is used in evaluate_basis_function
  current_ghosted_bf = NumericVector<Number>::build();
#ifdef LIBMESH_ENABLE_GHOSTED
  current_ghosted_bf->init (this->n_dofs(), this->n_local_dofs(),
                            get_dof_map().get_send_list(), false, GHOSTED);
#else
  current_ghosted_bf->init (this->n_dofs(), false, SERIAL);
#endif

  // Load up the inner product matrix
  // We only need one matrix in this class, so we
  // can set matrix to inner_product_matrix here
  if(!low_memory_mode)
  {
    matrix->zero();
    matrix->add(1., *inner_product_matrix);
  }
  else
  {
    assemble_inner_product_matrix(matrix);
  }

}

Number RBEIMSystem::evaluate_parametrized_function(unsigned int index, const Point& p)
{
  if(index >= get_n_parametrized_functions())
  {
    libMesh::err << "Error: We must have index < get_n_parametrized_functions() in evaluate_parametrized_function."
                 << std::endl;
    libmesh_error();
  }

  return parametrized_functions[index]->evaluate(get_current_parameters(), p);
}

unsigned int RBEIMSystem::get_n_affine_functions() const
{
  return n_vars() * rb_eval->get_n_basis_functions();
}

std::vector<Number> RBEIMSystem::evaluate_basis_function(unsigned int bf_index,
                                                         Elem& element,
                                                         const std::vector<Point>& qpoints)
{
  START_LOG("evaluate_current_basis_function()", "RBEIMSystem");
  
  // Load up basis function bf_index (does nothing if bf_index is already loaded)
  set_current_basis_function(bf_index);

  // Get local coordinates to feed these into compute_data().  
  // Note that the fe_type can safely be used from the 0-variable,
  // since the inverse mapping is the same for all FEFamilies
  std::vector<Point> mapped_qpoints;
  FEInterface::inverse_map (get_mesh().mesh_dimension(), 
   		            get_dof_map().variable_type(0),
                            &element, 
                            qpoints,
                            mapped_qpoints);

  const unsigned int current_variable_number = current_bf_index % n_vars();
  const FEType& fe_type = get_dof_map().variable_type(current_variable_number);
  
  std::vector<unsigned int> dof_indices_var;
  get_dof_map().dof_indices (&element, dof_indices_var, current_variable_number);

  std::vector<Number> values(dof_indices_var.size());

  for(unsigned int qp=0; qp<mapped_qpoints.size(); qp++)
  {
    FEComputeData data (get_equation_systems(), mapped_qpoints[qp]);
    FEInterface::compute_data (get_mesh().mesh_dimension(), fe_type, &element, data);

    values[qp] = 0.;
    for (unsigned int i=0; i<dof_indices_var.size(); i++)
      values[qp] += (*current_ghosted_bf)(dof_indices_var[i]) * data.shape[i];
  }

  STOP_LOG("evaluate_current_basis_function()", "RBEIMSystem");

  return values;
}

void RBEIMSystem::set_current_basis_function(unsigned int basis_function_index_in)
{
  START_LOG("set_current_basis_function()", "RBEIMSystem");

  if(basis_function_index_in > get_n_affine_functions())
  {
    libMesh::out << "Error: index cannot be larger than the number of affine functions in evaluate_affine_function"
                 << std::endl;
    libmesh_error();
  }
  
  if(basis_function_index_in != current_bf_index)
  {
    // Set member variable current_bf_index
    current_bf_index = basis_function_index_in;
    
    // First determine the basis function index implied by function_index
    unsigned int basis_function_id = current_bf_index/n_vars();
  
    // and create a ghosted version of the appropriate basis function
    rb_eval->get_basis_function(basis_function_id).localize
      (*current_ghosted_bf, this->get_dof_map().get_send_list());
  }
  
  STOP_LOG("set_current_basis_function()", "RBEIMSystem");
}

void RBEIMSystem::enrich_RB_space()
{
  START_LOG("enrich_RB_space()", "RBEIMSystem");

  RBEIMEvaluation* eim_eval = libmesh_cast_ptr<RBEIMEvaluation*>(rb_eval);

  // If we have at least one basis function we need to use
  // RB_solve, otherwise just use new_bf as is
  if(rb_eval->get_n_basis_functions() > 0)
  {
    // get the right-hand side vector for the EIM approximation
    // by sampling the parametrized function (stored in solution)
    // at the interpolation points
    unsigned int RB_size = rb_eval->get_n_basis_functions();
    DenseVector<Number> EIM_rhs(RB_size);
    for(unsigned int i=0; i<RB_size; i++)
    {
      EIM_rhs(i) = (*mesh_function)(eim_eval->interpolation_points[i]);
    }

    eim_eval->set_current_parameters( get_current_parameters() );
    eim_eval->RB_solve(EIM_rhs);

    // Load the "EIM residual" into solution by subtracting
    // the EIM approximation
    for(unsigned int i=0; i<rb_eval->get_n_basis_functions(); i++)
    {
      solution->add(-eim_eval->RB_solution(i), rb_eval->get_basis_function(i));
    }
  }

  // need to update since context uses current_local_solution
  update();

  // Find the quadrature point at which solution (which now stores
  // the "EIM residual") has maximum absolute value
  // by looping over the mesh
  Point optimal_point;
  Number optimal_value = 0.;
  unsigned int optimal_var;
  
  // Compute truth representation via projection
  const MeshBase& mesh = this->get_mesh();

  AutoPtr<FEMContext> c = this->build_context();
  FEMContext &context  = libmesh_cast_ref<FEMContext&>(*c);

  this->init_context(context);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    context.pre_fe_reinit(*this, *el);
    context.elem_fe_reinit();
      
    for(unsigned int var=0; var<n_vars(); var++)
    {
      unsigned int n_qpoints = context.element_qrule->n_points();

      for(unsigned int qp=0; qp<n_qpoints; qp++)
      {
        Number value = context.interior_value(var, qp);
        
        if( std::abs(value) > std::abs(optimal_value) )
        {
          optimal_value = value;
          optimal_point = context.element_fe_var[var]->get_xyz()[qp];
          optimal_var = var;
        }

      }
    }
  }
  
  Real global_abs_value = std::abs(optimal_value);
  unsigned int proc_ID_index;
  Parallel::maxloc(global_abs_value, proc_ID_index);
  
  // Broadcast the optimal point from proc_ID_index
  Parallel::broadcast(optimal_point, proc_ID_index);

  // Also broadcast the corresponding optimal_var and optimal_value
  Parallel::broadcast(optimal_var, proc_ID_index);
  Parallel::broadcast(optimal_value, proc_ID_index);

  // Scale the solution
  solution->scale(1./optimal_value);

  // Store optimal point in interpolation_points
  if(!performing_extra_greedy_step)
  {
    eim_eval->interpolation_points.push_back(optimal_point);
    eim_eval->interpolation_points_var.push_back(optimal_var);

    NumericVector<Number>* new_bf = NumericVector<Number>::build().release();
    new_bf->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    *new_bf = *solution;
    rb_eval->basis_functions.push_back( new_bf );
  }
  else
  {
    eim_eval->extra_interpolation_point = optimal_point;
    eim_eval->extra_interpolation_point_var = optimal_var;
  }
  
  STOP_LOG("enrich_RB_space()", "RBEIMSystem");
}

Real RBEIMSystem::compute_best_fit_error()
{
  START_LOG("compute_best_fit_error()", "RBEIMSystem");
  
  const unsigned int RB_size = rb_eval->get_n_basis_functions();
  
  // load the parametrized function into the solution vector
  truth_solve(-1);

  switch(best_fit_type_flag)
  {
    case(PROJECTION_BEST_FIT):
    {
      // compute the rhs by performing inner products
      DenseVector<Number> best_fit_rhs(RB_size);
      for(unsigned int i=0; i<RB_size; i++)
      {
        if(!low_memory_mode)
        {
          inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
        }
        else // In low memory mode we loaded the inner-product matrix into matrix during initialization
        {
          matrix->vector_mult(*inner_product_storage_vector, *solution);
        }
        best_fit_rhs(i) = inner_product_storage_vector->dot(rb_eval->get_basis_function(i));
      }

      // Now compute the best fit by an LU solve
      rb_eval->RB_solution.resize(RB_size);
      DenseMatrix<Number> RB_inner_product_matrix_N(RB_size);
      rb_eval->RB_inner_product_matrix.get_principal_submatrix(RB_size, RB_inner_product_matrix_N);

      RB_inner_product_matrix_N.lu_solve(best_fit_rhs, rb_eval->RB_solution);
      break;
    }
    case(EIM_BEST_FIT):
    {
      // Turn off error estimation here, we use the linfty norm instead
      rb_eval->evaluate_RB_error_bound = false;
      rb_eval->set_current_parameters( get_current_parameters() );
      rb_eval->RB_solve(RB_size);
      rb_eval->evaluate_RB_error_bound = true;
      break;
    }
    default:
    {
      libMesh::out << "Should not reach here" << std::endl;
      libmesh_error();
    }
  }
  
  // load the error into solution
  for(unsigned int i=0; i<rb_eval->get_n_basis_functions(); i++)
  {
    solution->add(-rb_eval->RB_solution(i), rb_eval->get_basis_function(i));
  }

  Real best_fit_error = solution->linfty_norm();
  
  STOP_LOG("compute_best_fit_error()", "RBEIMSystem");
  
  return best_fit_error;
}

Real RBEIMSystem::truth_solve(int plot_solution)
{
  START_LOG("truth_solve()", "RBEIMSystem");
        
//  matrix should have been set to inner_product_matrix during initialization
//  if(!low_memory_mode)
//  {
//    matrix->zero();
//    matrix->add(1., *inner_product_matrix);
//  }
//  else
//  {
//    assemble_inner_product_matrix(matrix);
//  }

  // Compute truth representation via projection
  const MeshBase& mesh = this->get_mesh();

  AutoPtr<FEMContext> c = this->build_context();
  FEMContext &context  = libmesh_cast_ref<FEMContext&>(*c);

  this->init_context(context);
  
  rhs->zero();

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    context.pre_fe_reinit(*this, *el);
    context.elem_fe_reinit();
      
    for(unsigned int var=0; var<n_vars(); var++)
    {
      const std::vector<Real> &JxW =
        context.element_fe_var[var]->get_JxW();

      const std::vector<std::vector<Real> >& phi =
        context.element_fe_var[var]->get_phi();

      const std::vector<Point> &xyz =
        context.element_fe_var[var]->get_xyz();

      unsigned int n_qpoints = context.element_qrule->n_points();
      unsigned int n_var_dofs = context.dof_indices_var[var].size();

      DenseSubVector<Number>& subresidual_var = *context.elem_subresiduals[var];

      for(unsigned int qp=0; qp<n_qpoints; qp++)
        for(unsigned int i=0; i != n_var_dofs; i++)
          subresidual_var(i) += JxW[qp] * evaluate_parametrized_function(var, xyz[qp]) * phi[i][qp];
    }

    // Apply constraints, e.g. periodic constraints
    this->get_dof_map().constrain_element_vector(context.elem_residual, context.dof_indices);

    // Add element vector to global vector
    rhs->add_vector(context.elem_residual, context.dof_indices);
  }
  
  // Solve to find the best fit, then solution stores the truth representation
  // of the function to be approximated
  solve();
  update(); // put the solution into current_local_solution as well in case we want to plot it
  solution->localize(*serialized_vector);
  
  if(reuse_preconditioner)
  {
    // After we've done a solve we can now reuse the preconditioner
    // because the matrix is not changing
    linear_solver->reuse_preconditioner(true);
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
//     libmesh_error();
  }

  if(plot_solution > 0)
  {
    const MeshBase& mesh = get_mesh();
    GMVIO(mesh).write_equation_systems ("truth.gmv",
                                        this->get_equation_systems());
  }
  
  STOP_LOG("truth_solve()", "RBEIMSystem");
  
  return 0.;
}

void RBEIMSystem::init_context(FEMContext &c)
{
  // default implementation of init_context
  // for compute_best_fit
  for(unsigned int var=0; var<n_vars(); var++)
  {
    c.element_fe_var[var]->get_JxW();
    c.element_fe_var[var]->get_phi();
    c.element_fe_var[var]->get_xyz();
  }
}

AutoPtr<RBEvaluation> RBEIMSystem::build_rb_evaluation()
{
  return AutoPtr<RBEvaluation>(new RBEIMEvaluation(*this));
}

void RBEIMSystem::update_RB_system_matrices()
{
  START_LOG("update_RB_system_matrices()", "RBEIMSystem");
  
  Parent::update_RB_system_matrices();

  unsigned int RB_size = rb_eval->get_n_basis_functions();
  
  RBEIMEvaluation* eim_eval = libmesh_cast_ptr<RBEIMEvaluation*>(rb_eval);
  
  // update the EIM interpolation matrix
  for(unsigned int j=0; j<RB_size; j++)
  {
    // Sample the basis functions at the
    // new interpolation point
    rb_eval->get_basis_function(j).localize(*serialized_vector);

    if(!performing_extra_greedy_step)
    {
      eim_eval->interpolation_matrix(RB_size-1,j) =
        (*mesh_function)(eim_eval->interpolation_points[RB_size-1]);
    }
    else
    {
      eim_eval->extra_interpolation_matrix_row(j) =
        (*mesh_function)(eim_eval->extra_interpolation_point);
    }
  }

  STOP_LOG("update_RB_system_matrices()", "RBEIMSystem");
}

void RBEIMSystem::update_system()
{
  libMesh::out << "Updating RB matrices" << std::endl;
  update_RB_system_matrices();
}

bool RBEIMSystem::greedy_termination_test(Real training_greedy_error, int)
{
  if(performing_extra_greedy_step)
  {
    libMesh::out << "Extra Greedy iteration finished." << std::endl;
    performing_extra_greedy_step = false;
    return true;
  }

  performing_extra_greedy_step = false;

  if(training_greedy_error < get_training_tolerance())
  {
    libMesh::out << "Specified error tolerance reached." << std::endl
                 << "Perform one more Greedy iteration for error bounds." << std::endl;
    performing_extra_greedy_step = true;
    return false;
  }

  if(rb_eval->get_n_basis_functions() >= this->get_Nmax())
  {
    libMesh::out << "Maximum number of basis functions reached: Nmax = "
              << get_Nmax() << "." << std::endl
              << "Perform one more Greedy iteration for error bounds." << std::endl;
    performing_extra_greedy_step = true;
    return false;
  }
  
  return false;
}

} // namespace libMesh

