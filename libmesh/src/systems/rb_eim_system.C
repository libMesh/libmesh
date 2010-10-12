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

namespace libMesh
{

RBEIMSystem::RBEIMSystem (EquationSystems& es,
		          const std::string& name,
		          const unsigned int number)
  : Parent(es, name, number),
    best_fit_type_flag(PROJECTION_BEST_FIT),
    mesh_function(NULL),
    performing_extra_greedy_step(false),
    current_variable_number(0),
    eval_error_estimate(false)
{}

RBEIMSystem::~RBEIMSystem ()
{
  delete mesh_function;
  mesh_function = NULL;
}

std::string RBEIMSystem::system_type () const
{
  return "RBEIMSystem";
}

void RBEIMSystem::init_data ()
{
  // Indicate that we need to compute the RB
  // inner product matrix in this case
  compute_RB_inner_product = true;
  
  // Indicate that we need the training set
  // for the Greedy to be the same on all
  // processors
  serial_training_set = true;

  Parent::init_data();
  
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

void RBEIMSystem::initialize_RB_system(bool online_mode)
{
  Parent::initialize_RB_system(online_mode);
  
  // Resize the data structures relevant to the EIM system
  interpolation_points.clear();
  interpolation_points_var.clear();
  interpolation_matrix.resize(Nmax,Nmax);
  
  // Resize the "extra" row due to the "extra Greedy step"
  extra_interpolation_matrix_row.resize(Nmax);
  
  if(initialize_calN_dependent_data)
  {
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
    // i.e. the basis funciton that is interpolated by evaluate_current_affine_function
    current_ghosted_bf = NumericVector<Number>::build();
#ifdef LIBMESH_ENABLE_GHOSTED
    current_ghosted_bf->init (this->n_dofs(), this->n_local_dofs(),
                              get_dof_map().get_send_list(), false, GHOSTED);
#else
    current_ghosted_bf->init (this->n_dofs(), false, SERIAL);
#endif

    if(!online_mode)
    {
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

  }
}

Number RBEIMSystem::evaluate_parametrized_function(unsigned int var, const Point& p)
{
  if(var >= get_n_parametrized_functions())
  {
    libMesh::err << "Error: We must have var < get_n_parametrized_functions() in evaluate_parametrized_function."
                 << std::endl;
    libmesh_error();
  }

  return parametrized_functions[var](p, *this);
}

unsigned int RBEIMSystem::get_n_affine_functions() const
{
  return n_vars() * get_n_basis_functions();
}

std::vector<Number> RBEIMSystem::evaluate_current_affine_function(Elem& element,
                                                                  const std::vector<Point>& qpoints)
{
  START_LOG("evaluate_current_affine_function()", "RBEIMSystem");

  // Get local coordinates to feed these into compute_data().  
  // Note that the fe_type can safely be used from the 0-variable,
  // since the inverse mapping is the same for all FEFamilies
  std::vector<Point> mapped_qpoints;
  FEInterface::inverse_map (get_mesh().mesh_dimension(), 
   		            get_dof_map().variable_type(0),
                            &element, 
                            qpoints,
                            mapped_qpoints);

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

  STOP_LOG("evaluate_current_affine_function()", "RBEIMSystem");

  return values;
}

void RBEIMSystem::cache_ghosted_basis_function(unsigned int function_index)
{
  START_LOG("cache_ghosted_basis_function()", "RBEIMSystem");

  if(function_index > get_n_affine_functions())
  {
    libMesh::out << "Error: index cannot be larger than the number of affine functions in evaluate_affine_function"
                 << std::endl;
    libmesh_error();
  }
        
  // First determine the basis function index implied by function_index
  unsigned int bf_index  = function_index/n_vars();
  
  // and create a ghosted version of the appropriate basis function
  basis_functions[bf_index]->localize
    (*current_ghosted_bf, this->get_dof_map().get_send_list());
    
  // Finally, store the index of the variable number that we will be considering
  current_variable_number = function_index % n_vars();
  
  STOP_LOG("cache_ghosted_basis_function()", "RBEIMSystem");
}

Real RBEIMSystem::RB_solve(unsigned int N)
{
  START_LOG("RB_solve()", "RBEIMSystem");

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

  // Get the rhs by sampling parametrized_function
  // at the first N interpolation_points
  DenseVector<Number> EIM_rhs(N);
  for(unsigned int i=0; i<N; i++)
  {
    EIM_rhs(i) = evaluate_parametrized_function(interpolation_points_var[i], interpolation_points[i]);
  }
  
  

  DenseMatrix<Number> interpolation_matrix_N;
  interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);
  
  interpolation_matrix_N.lu_solve(EIM_rhs, RB_solution);

  Real error_estimate = -1.;
  if(eval_error_estimate)
  {
    // Compute the a posteriori error bound
    // First, sample the parametrized function at x_{N+1}
    Number g_at_next_x;
    if(N == get_n_basis_functions())
      g_at_next_x = evaluate_parametrized_function(extra_interpolation_point_var, extra_interpolation_point);
    else
      g_at_next_x = evaluate_parametrized_function(interpolation_points_var[N], interpolation_points[N]);

    // Next, evaluate the EIM approximation at x_{N+1}
    Number EIM_approx_at_next_x = 0.;
    for(unsigned int j=0; j<N; j++)
      if(N == get_n_basis_functions())
        EIM_approx_at_next_x += RB_solution(j) * extra_interpolation_matrix_row(j);
      else
        EIM_approx_at_next_x += RB_solution(j) * interpolation_matrix(N,j);
      
    error_estimate = std::abs(g_at_next_x - EIM_approx_at_next_x);
  }
  
  STOP_LOG("RB_solve()", "RBEIMSystem");
  
  return error_estimate;
}

void RBEIMSystem::RB_solve(DenseVector<Number>& EIM_rhs)
{
  START_LOG("RB_solve()", "RBEIMSystem");
  
  if(EIM_rhs.size() > get_n_basis_functions())
  {
    libMesh::err << "ERROR: N cannot be larger than the number "
                 << "of basis functions in RB_solve" << std::endl;
    libmesh_error();
  }
  if(EIM_rhs.size()==0)
  {
    libMesh::err << "ERROR: N must be greater than 0 in RB_solve" << std::endl;
    libmesh_error();
  }
  
  const unsigned int N = EIM_rhs.size();
  DenseMatrix<Number> interpolation_matrix_N;
  interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);
  
  interpolation_matrix_N.lu_solve(EIM_rhs, RB_solution);
  
  STOP_LOG("RB_solve()", "RBEIMSystem");
}

void RBEIMSystem::enrich_RB_space()
{
  START_LOG("enrich_RB_space()", "RBEIMSystem");
  
  // If we have at least one basis function we need to use
  // RB_solve, otherwise just use new_bf as is
  if(get_n_basis_functions() > 0)
  {
    // get the right-hand side vector for the EIM approximation
    // by sampling the parametrized function (stored in solution)
    // at the interpolation points
    unsigned int RB_size = get_n_basis_functions();
    DenseVector<Number> EIM_rhs(RB_size);
    for(unsigned int i=0; i<RB_size; i++)
    {
      EIM_rhs(i) = (*mesh_function)(interpolation_points[i]);
    }
  
    RB_solve(EIM_rhs);

    // Load the "EIM residual" into solution by subtracting
    // the EIM approximation
    for(unsigned int i=0; i<get_n_basis_functions(); i++)
    {
      solution->add(-RB_solution(i), *basis_functions[i]);
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
  // First need to store optimal_point as a vector
  std::vector<Real> global_optimal_point(3);
  global_optimal_point[0] = optimal_point(0);
  global_optimal_point[1] = optimal_point(1);
  global_optimal_point[2] = optimal_point(2);
  Parallel::broadcast(global_optimal_point, proc_ID_index);
  optimal_point(0) = global_optimal_point[0];
  optimal_point(1) = global_optimal_point[1];
  optimal_point(2) = global_optimal_point[2];
  // Also broadcast the corresponding optimal_var and optimal_value
  Parallel::broadcast(optimal_var, proc_ID_index);
  Parallel::broadcast(optimal_value, proc_ID_index);

  // Scale the solution
  solution->scale(1./optimal_value);

  // Store optimal point in interpolation_points
  if(!performing_extra_greedy_step)
  {
    interpolation_points.push_back(optimal_point);
    interpolation_points_var.push_back(optimal_var);

    NumericVector<Number>* new_bf = NumericVector<Number>::build().release();
    new_bf->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    *new_bf = *solution;
    basis_functions.push_back( new_bf );
  }
  else
  {
    extra_interpolation_point = optimal_point;
    extra_interpolation_point_var = optimal_var;
  }
  
  STOP_LOG("enrich_RB_space()", "RBEIMSystem");
}

Real RBEIMSystem::compute_best_fit_error()
{
  START_LOG("compute_best_fit_error()", "RBEIMSystem");
  
  const unsigned int RB_size = get_n_basis_functions();
  
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
        best_fit_rhs(i) = inner_product_storage_vector->dot(*basis_functions[i]);
      }

      // Now compute the best fit by an LU solve
      RB_solution.resize(RB_size);
      DenseMatrix<Number> RB_inner_product_matrix_N(RB_size);
      RB_inner_product_matrix.get_principal_submatrix(RB_size, RB_inner_product_matrix_N);

      RB_inner_product_matrix_N.lu_solve(best_fit_rhs, RB_solution);
      break;
    }
    case(EIM_BEST_FIT):
    {
      // Turn off error estimation here, we use the linfty norm instead
      eval_error_estimate = false;
      RB_solve(RB_size);
      eval_error_estimate = true;
      break;
    }
    default:
    {
      libMesh::out << "Should not reach here" << std::endl;
      libmesh_error();
    }
  }
  
  // load the error into solution
  for(unsigned int i=0; i<get_n_basis_functions(); i++)
  {
    solution->add(-RB_solution(i), *basis_functions[i]);
  }

  Real best_fit_error = solution->linfty_norm();
  
  STOP_LOG("compute_best_fit_error()", "RBEIMSystem");
  
  return best_fit_error;
}

Real RBEIMSystem::truth_solve(int plot_solution)
{
  START_LOG("truth_solve()", "RBEIMSystem");
        
  if(!initialize_calN_dependent_data)
  {
    libMesh::err << "Error: We must initialize the calN dependent "
                 << "data structures in order to load the truth solution."
                 << std::endl;
    libmesh_error();
  }
  
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
    linear_solver->same_preconditioner = true;
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

void RBEIMSystem::clear_basis_function_dependent_data()
{
  Parent::clear_basis_function_dependent_data();
  
  interpolation_points.clear();
  interpolation_points_var.clear();
}

void RBEIMSystem::update_RB_system_matrices()
{
  START_LOG("update_RB_system_matrices()", "RBEIMSystem");
  
  Parent::update_RB_system_matrices();

  unsigned int RB_size = get_n_basis_functions();
  
  // update the EIM interpolation matrix
  for(unsigned int j=0; j<RB_size; j++)
  {
    // Sample the basis functions at the
    // new interpolation point
    basis_functions[j]->localize(*serialized_vector);

    if(!performing_extra_greedy_step)
    {
      interpolation_matrix(RB_size-1,j) = (*mesh_function)(interpolation_points[RB_size-1]);
    }
    else
    {
      extra_interpolation_matrix_row(j) = (*mesh_function)(extra_interpolation_point);
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

  if(get_n_basis_functions() >= this->get_Nmax())
  {
    libMesh::out << "Maximum number of basis functions reached: Nmax = "
              << get_Nmax() << "." << std::endl
              << "Perform one more Greedy iteration for error bounds." << std::endl;
    performing_extra_greedy_step = true;
    return false;
  }
  
  return false;
}

void RBEIMSystem::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "RBEIMSystem");

  Parent::write_offline_data_to_files(directory_name);

  const unsigned int n_bfs = get_n_basis_functions();
  libmesh_assert( n_bfs <= Nmax );

  const unsigned int precision_level = 14;

  if(libMesh::processor_id() == 0)
  {
    // Next write out the interpolation_matrix
    std::ofstream interpolation_matrix_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/interpolation_matrix.dat";
      interpolation_matrix_out.open(file_name.str().c_str());
    }
    if ( !interpolation_matrix_out.good() )
    {
      libMesh::err << "Error opening interpolation_matrix.dat" << std::endl;
      libmesh_error();
    }
    interpolation_matrix_out.precision(precision_level);
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<=i; j++)
      {
        interpolation_matrix_out << std::scientific
          << interpolation_matrix(i,j) << " ";
      }
    }
    
    // Also, write out the "extra" row
    std::ofstream extra_interpolation_matrix_row_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/extra_interpolation_matrix_row.dat";
      extra_interpolation_matrix_row_out.open(file_name.str().c_str());
    }
    if ( !extra_interpolation_matrix_row_out.good() )
    {
      libMesh::err << "Error opening extra_interpolation_matrix_row.dat" << std::endl;
      libmesh_error();
    }
    extra_interpolation_matrix_row_out.precision(precision_level);
    for(unsigned int j=0; j<n_bfs; j++)
      extra_interpolation_matrix_row_out << std::scientific
          << extra_interpolation_matrix_row(j) << " ";
    extra_interpolation_matrix_row_out.close();
    
    // Next write out interpolation_points
    std::ofstream interpolation_points_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/interpolation_points.dat";
      interpolation_points_out.open(file_name.str().c_str());
    }
    if ( !interpolation_points_out.good() )
    {
      libMesh::err << "Error opening interpolation_points.dat" << std::endl;
      libmesh_error();
    }
    interpolation_points_out.precision(precision_level);
    for(unsigned int i=0; i<n_bfs; i++)
      interpolation_points_out << std::scientific
          << interpolation_points[i](0) << " "
          << interpolation_points[i](1) << " "
          << interpolation_points[i](2) << " ";

    // Also, write out the "extra" interpolation point
    std::ofstream extra_interpolation_point_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/extra_interpolation_point.dat";
      extra_interpolation_point_out.open(file_name.str().c_str());
    }
    if ( !extra_interpolation_point_out.good() )
    {
      libMesh::err << "Error opening extra_interpolation_point.dat" << std::endl;
      libmesh_error();
    }
    extra_interpolation_point_out.precision(precision_level);
    extra_interpolation_point_out << std::scientific
          << extra_interpolation_point(0) << " "
          << extra_interpolation_point(1) << " "
          << extra_interpolation_point(2) << " ";
    extra_interpolation_point_out.close();
    
    // Next write out interpolation_points_var
    std::ofstream interpolation_points_var_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/interpolation_points_var.dat";
      interpolation_points_var_out.open(file_name.str().c_str());
    }
    if ( !interpolation_points_var_out.good() )
    {
      libMesh::err << "Error opening interpolation_points_var.dat" << std::endl;
      libmesh_error();
    }
    interpolation_points_var_out.precision(precision_level);
    for(unsigned int i=0; i<n_bfs; i++)
      interpolation_points_var_out << std::scientific
          << interpolation_points_var[i] << " ";

    // Also, write out the "extra" interpolation variable
    std::ofstream extra_interpolation_point_var_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/extra_interpolation_point_var.dat";
      extra_interpolation_point_var_out.open(file_name.str().c_str());
    }
    if ( !extra_interpolation_point_var_out.good() )
    {
      libMesh::err << "Error opening extra_interpolation_point_var.dat" << std::endl;
      libmesh_error();
    }
    extra_interpolation_point_var_out.precision(precision_level);
    extra_interpolation_point_var_out << std::scientific
          << extra_interpolation_point_var << " ";
    extra_interpolation_point_var_out.close();
  }

  STOP_LOG("write_offline_data_to_files()", "RBEIMSystem");
}

void RBEIMSystem::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "RBEIMSystem");

  Parent::read_offline_data_from_files(directory_name);
  
  // First, find out how many basis functions we had when Greedy terminated
  // This was set in RBSystem::read_offline_data_from_files
  unsigned int n_bfs = this->get_n_basis_functions();
  
  // Read in the interpolation matrix
  std::ifstream interpolation_matrix_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/interpolation_matrix.dat";
    interpolation_matrix_in.open(file_name.str().c_str());
  }
  if ( !interpolation_matrix_in.good() )
  {
    libMesh::err << "Error opening interpolation_matrix.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<n_bfs; i++)
  {
    for(unsigned int j=0; j<=i; j++)
    {
      Number value;
      interpolation_matrix_in >> value;
      interpolation_matrix(i,j) = value;
    }
  }

  // Also, read in the "extra" row
  std::ifstream extra_interpolation_matrix_row_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/extra_interpolation_matrix_row.dat";
    extra_interpolation_matrix_row_in.open(file_name.str().c_str());
  }
  if ( !extra_interpolation_matrix_row_in.good() )
  {
    libMesh::err << "Error opening extra_interpolation_matrix_row.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int j=0; j<n_bfs; j++)
  {
    Number value;
    extra_interpolation_matrix_row_in >> value;
    extra_interpolation_matrix_row(j) = value;
  }
  extra_interpolation_matrix_row_in.close();

  // Next read in interpolation_points
  std::ifstream interpolation_points_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/interpolation_points.dat";
    interpolation_points_in.open(file_name.str().c_str());
  }
  if ( !interpolation_points_in.good() )
  {
    libMesh::err << "Error opening interpolation_points.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<n_bfs; i++)
  {
    Real x_val, y_val, z_val;
    interpolation_points_in >> x_val;
    interpolation_points_in >> y_val;
    interpolation_points_in >> z_val;
    Point p(x_val, y_val, z_val);
    interpolation_points.push_back(p);
  }
  interpolation_points_in.close();
  
  // Also, read in the extra interpolation point
  std::ifstream extra_interpolation_point_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/extra_interpolation_point.dat";
    extra_interpolation_point_in.open(file_name.str().c_str());
  }
  if ( !extra_interpolation_point_in.good() )
  {
    libMesh::err << "Error opening extra_interpolation_point.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<n_bfs; i++)
  {
    Real x_val, y_val, z_val;
    extra_interpolation_point_in >> x_val;
    extra_interpolation_point_in >> y_val;
    extra_interpolation_point_in >> z_val;
    Point p(x_val, y_val, z_val);
    extra_interpolation_point = p;
  }
  extra_interpolation_point_in.close();
  

  // Next read in interpolation_points_var
  std::ifstream interpolation_points_var_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/interpolation_points_var.dat";
    interpolation_points_var_in.open(file_name.str().c_str());
  }
  if ( !interpolation_points_var_in.good() )
  {
    libMesh::err << "Error opening interpolation_points_var.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<=n_bfs; i++)
  {
    unsigned int var;
    interpolation_points_var_in >> var;
    interpolation_points_var.push_back(var);
  }
  interpolation_points_var_in.close();
  
  // Also, read in extra_interpolation_point_var
  std::ifstream extra_interpolation_point_var_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/extra_interpolation_point_var.dat";
    extra_interpolation_point_var_in.open(file_name.str().c_str());
  }
  if ( !extra_interpolation_point_var_in.good() )
  {
    libMesh::err << "Error opening extra_interpolation_point_var.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<=n_bfs; i++)
  {
    unsigned int var;
    extra_interpolation_point_var_in >> var;
    extra_interpolation_point_var = var;
  }
  extra_interpolation_point_var_in.close();
  
  STOP_LOG("read_offline_data_from_files()", "RBEIMSystem");
}

} // namespace libMesh

