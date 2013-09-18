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
#include <fstream>
#include <sstream>
#include "libmesh/exodusII_io.h"
#include "libmesh/fem_context.h"

#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_evaluation.h"

namespace libMesh
{

RBEIMConstruction::RBEIMConstruction (EquationSystems& es,
		                      const std::string& name_in,
		                      const unsigned int number_in)
  : Parent(es, name_in, number_in),
    best_fit_type_flag(PROJECTION_BEST_FIT),
    _parametrized_functions_in_training_set_initialized(false),
    _mesh_function(NULL),
    _performing_extra_greedy_step(false)
{
  // We cannot do rb_solve with an empty
  // "rb space" with EIM
  use_empty_rb_solve_in_greedy = false;

  // Indicate that we need to compute the RB
  // inner product matrix in this case
  compute_RB_inner_product = true;

  // Indicate that we need the training set
  // for the Greedy to be the same on all
  // processors
  serial_training_set = true;

  // attach empty RBAssemblyExpansion object
  set_rb_assembly_expansion(_empty_rb_assembly_expansion);
}

RBEIMConstruction::~RBEIMConstruction ()
{
  this->clear();
}

void RBEIMConstruction::clear()
{
  Parent::clear();

  // clear the mesh function
  delete _mesh_function;
  _mesh_function = NULL;

  // clear the eim assembly vector
  for(unsigned int i=0; i<_rb_eim_assembly_objects.size(); i++)
  {
    delete _rb_eim_assembly_objects[i];
  }
  _rb_eim_assembly_objects.clear();

  // clear the parametrized functions from the training set
  for(unsigned int i=0; i<_parametrized_functions_in_training_set.size(); i++)
  {
    if (_parametrized_functions_in_training_set[i])
    {
      _parametrized_functions_in_training_set[i]->clear();
      delete _parametrized_functions_in_training_set[i];
      _parametrized_functions_in_training_set[i] = NULL;
    }
  }
  _parametrized_functions_in_training_set_initialized = false;
}

void RBEIMConstruction::process_parameters_file (const std::string& parameters_filename)
{
  Parent::process_parameters_file(parameters_filename);

  GetPot infile(parameters_filename);

  std::string best_fit_type_string = infile("best_fit_type","projection");
  set_best_fit_type_flag(best_fit_type_string);
}

void RBEIMConstruction::set_best_fit_type_flag (const std::string& best_fit_type_string)
{
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
}

void RBEIMConstruction::print_info()
{
  Parent::print_info();

  // Print out setup info
  libMesh::out << std::endl << "RBEIMConstruction parameters:" << std::endl;
  if(best_fit_type_flag == PROJECTION_BEST_FIT)
  {
    libMesh::out << "best fit type: projection" << std::endl;
  }
  else
  if(best_fit_type_flag == EIM_BEST_FIT)
  {
    libMesh::out << "best fit type: eim" << std::endl;
  }
  libMesh::out << std::endl;
}

void RBEIMConstruction::init_data()
{
  // set the coupling matrix to be diagonal
  _coupling_matrix.resize(n_vars());
  for(unsigned int var1=0; var1<n_vars(); var1++)
    for(unsigned int var2=0; var2<n_vars(); var2++)
    {
      unsigned char value = (var1==var2) ? 1 : 0;
      _coupling_matrix(var1,var2) = value;
    }

  this->get_dof_map()._dof_coupling = &_coupling_matrix;

  Parent::init_data();
}

void RBEIMConstruction::initialize_rb_construction()
{
  Parent::initialize_rb_construction();

  // initialize a serial vector that we will use for MeshFunction evaluations
  _ghosted_meshfunction_vector = NumericVector<Number>::build(this->comm());
  _ghosted_meshfunction_vector->init (this->n_dofs(), this->n_local_dofs(),
                                      this->get_dof_map().get_send_list(), false,
                                      GHOSTED);

  // Initialize the MeshFunction for interpolating the
  // solution vector at quadrature points
  std::vector<unsigned int> vars;
  get_all_variable_numbers(vars);
  _mesh_function = new MeshFunction(get_equation_systems(),
                                    *_ghosted_meshfunction_vector,
                                    get_dof_map(),
                                    vars);
  _mesh_function->init();

  // Load up the inner product matrix
  // We only need one matrix in this class, so we
  // can set matrix to inner_product_matrix here
  {
    matrix->zero();
    matrix->close();
    matrix->add(1., *inner_product_matrix);
  }

}

Real RBEIMConstruction::train_reduced_basis(const std::string& directory_name,
                                            const bool resize_rb_eval_data)
{
  // precompute all the parametrized functions that we'll use in the greedy
  initialize_parametrized_functions_in_training_set();

  return Parent::train_reduced_basis(directory_name, resize_rb_eval_data);
}

Number RBEIMConstruction::evaluate_mesh_function(unsigned int var_number,
                                                 Point p)
{
  DenseVector<Number> values;
  (*_mesh_function)(p,
                    /*time*/ 0.,
                    values);

  // We evaluated the mesh function, but it will only return a valid set of values on one processor
  // (values will be empty on all other processors) so we need to broadcast those valid values
  // to all processors.
  Number value = 0;
  unsigned int root_id=0;
  if(values.size() != 0)
  {
    root_id = this->processor_id();
    value = values(var_number);
  }

  // root_id may be non-zero on more than one processor due to ghost elements
  // so use this->comm().max to get just one proc id
  this->comm().max(root_id);

  // Then broadcast the result
  this->comm().broadcast(value, root_id);

  return value;
}

void RBEIMConstruction::initialize_eim_assembly_objects()
{
  _rb_eim_assembly_objects.clear();
  for(unsigned int i=0; i<get_rb_evaluation().get_n_basis_functions(); i++)
  {
    _rb_eim_assembly_objects.push_back( build_eim_assembly(i).release() );
  }
}

std::vector<ElemAssembly*> RBEIMConstruction::get_eim_assembly_objects()
{
  return _rb_eim_assembly_objects;
}

void RBEIMConstruction::enrich_RB_space()
{
  START_LOG("enrich_RB_space()", "RBEIMConstruction");

  // put solution in _ghosted_meshfunction_vector so we can access it from the mesh function
  // this allows us to compute EIM_rhs appropriately
  solution->localize(*_ghosted_meshfunction_vector, this->get_dof_map().get_send_list());

  RBEIMEvaluation& eim_eval = libmesh_cast_ref<RBEIMEvaluation&>(get_rb_evaluation());

  // If we have at least one basis function we need to use
  // rb_solve to find the EIM interpolation error, otherwise just use solution as is
  if(get_rb_evaluation().get_n_basis_functions() > 0)
  {
    // get the right-hand side vector for the EIM approximation
    // by sampling the parametrized function (stored in solution)
    // at the interpolation points
    unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();
    DenseVector<Number> EIM_rhs(RB_size);
    for(unsigned int i=0; i<RB_size; i++)
    {
      EIM_rhs(i) = evaluate_mesh_function( eim_eval.interpolation_points_var[i],
                                           eim_eval.interpolation_points[i] );
    }

    eim_eval.set_parameters( get_parameters() );
    eim_eval.rb_solve(EIM_rhs);

    // Load the "EIM residual" into solution by subtracting
    // the EIM approximation
    for(unsigned int i=0; i<get_rb_evaluation().get_n_basis_functions(); i++)
    {
      solution->add(-eim_eval.RB_solution(i), get_rb_evaluation().get_basis_function(i));
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
  subdomain_id_type optimal_subdomain;

  // Compute truth representation via projection
  const MeshBase& mesh = this->get_mesh();

  AutoPtr<DGFEMContext> c = this->build_context();
  DGFEMContext &context  = libmesh_cast_ref<DGFEMContext&>(*c);

  this->init_context(context);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    context.pre_fe_reinit(*this, *el);
    context.elem_fe_reinit();

    for(unsigned int var=0; var<n_vars(); var++)
    {
      unsigned int n_qpoints = context.get_element_qrule().n_points();

      for(unsigned int qp=0; qp<n_qpoints; qp++)
      {
        Number value = context.interior_value(var, qp);

        if( std::abs(value) > std::abs(optimal_value) )
        {
          FEBase* elem_fe = NULL;
          context.get_element_fe( var, elem_fe );
          optimal_value = value;
          optimal_point = elem_fe->get_xyz()[qp];
          optimal_var = var;
          optimal_subdomain = (*el)->subdomain_id();
        }

      }
    }
  }

  Real global_abs_value = std::abs(optimal_value);
  unsigned int proc_ID_index;
  this->comm().maxloc(global_abs_value, proc_ID_index);

  // Broadcast the optimal point from proc_ID_index
  this->comm().broadcast(optimal_point, proc_ID_index);

  // Also broadcast the corresponding optimal_var, optimal_value, and optimal_subdomain
  this->comm().broadcast(optimal_var, proc_ID_index);
  this->comm().broadcast(optimal_value, proc_ID_index);
  this->comm().broadcast(optimal_subdomain, proc_ID_index);

  // Scale the solution
  solution->scale(1./optimal_value);

  // Store optimal point in interpolation_points
  if(!_performing_extra_greedy_step)
  {
    eim_eval.interpolation_points.push_back(optimal_point);
    eim_eval.interpolation_points_var.push_back(optimal_var);
    eim_eval.interpolation_points_subdomain.push_back(optimal_subdomain);

    NumericVector<Number>* new_bf = NumericVector<Number>::build(this->comm()).release();
    new_bf->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    *new_bf = *solution;
    get_rb_evaluation().basis_functions.push_back( new_bf );
  }
  else
  {
    eim_eval.extra_interpolation_point = optimal_point;
    eim_eval.extra_interpolation_point_var = optimal_var;
    eim_eval.extra_interpolation_point_subdomain = optimal_subdomain;
  }

  STOP_LOG("enrich_RB_space()", "RBEIMConstruction");
}

void RBEIMConstruction::initialize_parametrized_functions_in_training_set()
{
  if(!serial_training_set)
  {
    libMesh::err << "Error: We must have serial_training_set==true in "
                 << "RBEIMConstruction::initialize_parametrized_functions_in_training_set"
                 << std::endl;
    libmesh_error();
  }

  libMesh::out << "Initializing parametrized functions in training set..." << std::endl;
  // initialize rb_eval's parameters
  get_rb_evaluation().initialize_parameters(*this);

  _parametrized_functions_in_training_set.resize( get_n_training_samples() );
  for(unsigned int i=0; i<get_n_training_samples(); i++)
  {
    set_params_from_training_set(i);
    truth_solve(-1);

    _parametrized_functions_in_training_set[i] = solution->clone().release();

    libMesh::out << "Completed solve for training sample " << (i+1) << " of " << get_n_training_samples() << std::endl;
  }

  _parametrized_functions_in_training_set_initialized = true;

  libMesh::out << "Parametrized functions in training set initialized" << std::endl << std::endl;
}


Real RBEIMConstruction::compute_best_fit_error()
{
  START_LOG("compute_best_fit_error()", "RBEIMConstruction");

  const unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  // load up the parametrized function for the current parameters
  truth_solve(-1);

  switch(best_fit_type_flag)
  {
    // Perform an L2 projection in order to find an approximation to solution (from truth_solve above)
    case(PROJECTION_BEST_FIT):
    {
      // compute the rhs by performing inner products
      DenseVector<Number> best_fit_rhs(RB_size);
      for(unsigned int i=0; i<RB_size; i++)
      {
        inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
        
        best_fit_rhs(i) = inner_product_storage_vector->dot(get_rb_evaluation().get_basis_function(i));
      }

      // Now compute the best fit by an LU solve
      get_rb_evaluation().RB_solution.resize(RB_size);
      DenseMatrix<Number> RB_inner_product_matrix_N(RB_size);
      get_rb_evaluation().RB_inner_product_matrix.get_principal_submatrix(RB_size, RB_inner_product_matrix_N);

      RB_inner_product_matrix_N.lu_solve(best_fit_rhs, get_rb_evaluation().RB_solution);
      break;
    }
    // Perform EIM solve in order to find the approximation to solution
    // (rb_solve provides the EIM basis function coefficients used below)
    case(EIM_BEST_FIT):
    {
      // Turn off error estimation for this rb_solve, we use the linfty norm instead
      get_rb_evaluation().evaluate_RB_error_bound = false;
      get_rb_evaluation().set_parameters( get_parameters() );
      get_rb_evaluation().rb_solve(RB_size);
      get_rb_evaluation().evaluate_RB_error_bound = true;
      break;
    }
    default:
    {
      libMesh::out << "Should not reach here" << std::endl;
      libmesh_error();
    }
  }

  // load the error into solution
  for(unsigned int i=0; i<get_rb_evaluation().get_n_basis_functions(); i++)
  {
    solution->add(-get_rb_evaluation().RB_solution(i), get_rb_evaluation().get_basis_function(i));
  }

  Real best_fit_error = solution->linfty_norm();

  STOP_LOG("compute_best_fit_error()", "RBEIMConstruction");

  return best_fit_error;
}

Real RBEIMConstruction::truth_solve(int plot_solution)
{
  START_LOG("truth_solve()", "RBEIMConstruction");

//  matrix should have been set to inner_product_matrix during initialization
//  {
//    matrix->zero();
//    matrix->add(1., *inner_product_matrix);
//  }

  int training_parameters_found_index = -1;
  if( _parametrized_functions_in_training_set_initialized )
  {
    // Check if parameters are in the training set. If so, we can just load the
    // solution from _parametrized_functions_in_training_set

    for(unsigned int i=0; i<get_n_training_samples(); i++)
    {
      if(get_parameters() == get_params_from_training_set(i))
      {
        training_parameters_found_index = i;
        break;
      }
    }
  }

  // If the parameters are in the training set, just copy the solution vector
  if(training_parameters_found_index >= 0)
  {
    *solution = *_parametrized_functions_in_training_set[training_parameters_found_index];
    update(); // put the solution into current_local_solution as well
  }
  // Otherwise, we have to compute the projection
  else
  {
    RBEIMEvaluation& eim_eval = libmesh_cast_ref<RBEIMEvaluation&>(get_rb_evaluation());
    eim_eval.set_parameters( get_parameters() );

    // Compute truth representation via projection
    const MeshBase& mesh = this->get_mesh();

    AutoPtr<DGFEMContext> c = this->build_context();
    DGFEMContext &context  = libmesh_cast_ref<DGFEMContext&>(*c);

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
        FEBase* elem_fe = NULL;
        context.get_element_fe( var, elem_fe );
        const std::vector<Real> &JxW = elem_fe->get_JxW();

        const std::vector<std::vector<Real> >& phi = elem_fe->get_phi();

        const std::vector<Point> &xyz = elem_fe->get_xyz();

        unsigned int n_qpoints = context.get_element_qrule().n_points();
        unsigned int n_var_dofs = libmesh_cast_int<unsigned int>
	  (context.get_dof_indices( var ).size());

        DenseSubVector<Number>& subresidual_var = context.get_elem_residual( var );

        for(unsigned int qp=0; qp<n_qpoints; qp++)
          for(unsigned int i=0; i != n_var_dofs; i++)
            subresidual_var(i) += JxW[qp] * eim_eval.evaluate_parametrized_function(var, xyz[qp], (*el)->subdomain_id()) * phi[i][qp];
      }

      // Apply constraints, e.g. periodic constraints
      this->get_dof_map().constrain_element_vector(context.get_elem_residual(), context.get_dof_indices() );

      // Add element vector to global vector
      rhs->add_vector(context.get_elem_residual(), context.get_dof_indices() );
    }

    // Solve to find the best fit, then solution stores the truth representation
    // of the function to be approximated
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
  //     libmesh_error();
    }

    if(reuse_preconditioner)
    {
      // After we've done a solve we can now reuse the preconditioner
      // because the matrix is not changing
      linear_solver->reuse_preconditioner(true);
    }
  }

  if(plot_solution > 0)
  {
#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(get_mesh()).write_equation_systems ("truth.e",
                                                    this->get_equation_systems());
#endif
  }

  STOP_LOG("truth_solve()", "RBEIMConstruction");

  return 0.;
}

void RBEIMConstruction::init_context(FEMContext &c)
{
  // default implementation of init_context
  // for compute_best_fit
  for(unsigned int var=0; var<n_vars(); var++)
  {
    FEBase* elem_fe = NULL;
    c.get_element_fe( var, elem_fe );
    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_xyz();
  }
}

void RBEIMConstruction::update_RB_system_matrices()
{
  START_LOG("update_RB_system_matrices()", "RBEIMConstruction");

  Parent::update_RB_system_matrices();

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  RBEIMEvaluation& eim_eval = libmesh_cast_ref<RBEIMEvaluation&>(get_rb_evaluation());

  // update the EIM interpolation matrix
  for(unsigned int j=0; j<RB_size; j++)
  {
    // Sample the basis functions at the
    // new interpolation point
    get_rb_evaluation().get_basis_function(j).localize(*_ghosted_meshfunction_vector, this->get_dof_map().get_send_list());

    if(!_performing_extra_greedy_step)
    {
      eim_eval.interpolation_matrix(RB_size-1,j) =
        evaluate_mesh_function( eim_eval.interpolation_points_var[RB_size-1],
                                eim_eval.interpolation_points[RB_size-1] );
    }
    else
    {
      eim_eval.extra_interpolation_matrix_row(j) =
        evaluate_mesh_function( eim_eval.extra_interpolation_point_var,
                                eim_eval.extra_interpolation_point );
    }
  }

  STOP_LOG("update_RB_system_matrices()", "RBEIMConstruction");
}

Real RBEIMConstruction::get_RB_error_bound()
{
  return compute_best_fit_error();
}

void RBEIMConstruction::update_system()
{
  libMesh::out << "Updating RB matrices" << std::endl;
  update_RB_system_matrices();
}

bool RBEIMConstruction::greedy_termination_test(Real training_greedy_error, int)
{
  if(_performing_extra_greedy_step)
  {
    libMesh::out << "Extra Greedy iteration finished." << std::endl;
    _performing_extra_greedy_step = false;
    return true;
  }

  _performing_extra_greedy_step = false;

  if(training_greedy_error < get_training_tolerance())
  {
    libMesh::out << "Specified error tolerance reached." << std::endl
                 << "Perform one more Greedy iteration for error bounds." << std::endl;
    _performing_extra_greedy_step = true;
    return false;
  }

  if(get_rb_evaluation().get_n_basis_functions() >= this->get_Nmax())
  {
    libMesh::out << "Maximum number of basis functions reached: Nmax = "
              << get_Nmax() << "." << std::endl
              << "Perform one more Greedy iteration for error bounds." << std::endl;
    _performing_extra_greedy_step = true;
    return false;
  }

  return false;
}

} // namespace libMesh
