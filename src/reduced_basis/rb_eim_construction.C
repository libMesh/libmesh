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

#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_evaluation.h"

namespace libMesh
{

RBEIMConstruction::RBEIMConstruction (EquationSystems & es,
                                      const std::string & name_in,
                                      const unsigned int number_in)
  : Parent(es, name_in, number_in),
    best_fit_type_flag(PROJECTION_BEST_FIT),
    _parametrized_functions_in_training_set_initialized(false),
    _point_locator_tol(TOLERANCE)
{
  _explicit_system_name = name_in + "_explicit_sys";

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

  // We set implicit_neighbor_dofs = false. This is important when we use
  // DISCONTINUOUS basis functions, since by default libMesh sets
  // implicit_neighbor_dofs = true for "all discontinuous" systems, which
  // results in a larger sparsity: dofs on neighboring elements are added
  // to the sparsity pattern since this is typically required for DG or FV
  // discretizations. Since we're only doing L2 projects here, we do not
  // need extra dofs in the sparsity pattern, so we set implicit neighbor
  // dofs to false.
  get_dof_map().set_implicit_neighbor_dofs(false);
}

RBEIMConstruction::~RBEIMConstruction ()
{
  this->clear();
}

void RBEIMConstruction::clear()
{
  Parent::clear();

  // clear the mesh function
  _mesh_function.reset();

  // clear the eim assembly vector
  for (std::size_t i=0; i<_rb_eim_assembly_objects.size(); i++)
    delete _rb_eim_assembly_objects[i];
  _rb_eim_assembly_objects.clear();

  // clear the parametrized functions from the training set
  for (std::size_t i=0; i<_parametrized_functions_in_training_set.size(); i++)
    {
      if (_parametrized_functions_in_training_set[i])
        {
          _parametrized_functions_in_training_set[i]->clear();
          delete _parametrized_functions_in_training_set[i];
          _parametrized_functions_in_training_set[i] = libmesh_nullptr;
        }
    }
  _parametrized_functions_in_training_set_initialized = false;

  for (std::size_t i=0; i<_matrix_times_bfs.size(); i++)
    {
      if(_matrix_times_bfs[i])
        {
          _matrix_times_bfs[i]->clear();
          delete _matrix_times_bfs[i];
          _matrix_times_bfs[i] = libmesh_nullptr;
        }
    }
}

void RBEIMConstruction::process_parameters_file (const std::string & parameters_filename)
{
  Parent::process_parameters_file(parameters_filename);

  GetPot infile(parameters_filename);

  std::string best_fit_type_string = infile("best_fit_type","projection");
  set_best_fit_type_flag(best_fit_type_string);
}

void RBEIMConstruction::set_best_fit_type_flag (const std::string & best_fit_type_string)
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
      libmesh_error_msg("Error: invalid best_fit_type in input file");
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
  // Add the ExplicitSystem that we use to store the EIM basis functions
  get_equation_systems().add_system<ExplicitSystem>(_explicit_system_name);

  init_implicit_system();
  init_explicit_system();

  Parent::init_data();
}

void RBEIMConstruction::initialize_rb_construction(bool skip_matrix_assembly,
                                                   bool skip_vector_assembly)
{
  Parent::initialize_rb_construction(skip_matrix_assembly, skip_vector_assembly);

  // initialize a serial vector that we will use for MeshFunction evaluations
  _ghosted_meshfunction_vector = NumericVector<Number>::build(this->comm());
  _ghosted_meshfunction_vector->init (get_explicit_system().n_dofs(), get_explicit_system().n_local_dofs(),
                                      get_explicit_system().get_dof_map().get_send_list(), false,
                                      GHOSTED);

  // Initialize the MeshFunction for interpolating the
  // solution vector at quadrature points
  std::vector<unsigned int> vars;
  get_explicit_system().get_all_variable_numbers(vars);
  _mesh_function.reset(new MeshFunction(get_equation_systems(),
                                        *_ghosted_meshfunction_vector,
                                        get_explicit_system().get_dof_map(),
                                        vars));
  _mesh_function->init();

  // inner_product_solver performs solves with the same matrix every time
  // hence we can set reuse_preconditioner(true).
  inner_product_solver->reuse_preconditioner(true);

  init_dof_map_between_systems();
}

Real RBEIMConstruction::train_reduced_basis(const bool resize_rb_eval_data)
{
  // precompute all the parametrized functions that we'll use in the greedy
  initialize_parametrized_functions_in_training_set();

  return Parent::train_reduced_basis(resize_rb_eval_data);
}

Number RBEIMConstruction::evaluate_mesh_function(unsigned int var_number,
                                                 Point p)
{
  _mesh_function->set_point_locator_tolerance( get_point_locator_tol() );

  DenseVector<Number> values;
  (*_mesh_function)(p,
                    /*time*/ 0.,
                    values);

  // We evaluated the mesh function, but it will only return a valid set of values on one processor
  // (values will be empty on all other processors) so we need to broadcast those valid values
  // to all processors.
  Number value = 0;
  unsigned int root_id=0;
  unsigned int check_for_valid_value = 0;
  if(values.size() != 0)
    {
      root_id = this->processor_id();
      value = values(var_number);
      check_for_valid_value = 1;
    }

  // If this sum is zero, then we didn't enter the if block above on any processor. In that
  // case we should throw an error.
  this->comm().sum(check_for_valid_value);
  if(check_for_valid_value == 0)
    {
      libmesh_error_msg("MeshFunction evaluation failed on all processors");
    }

  // root_id may be non-zero on more than one processor due to ghost elements
  // so use this->comm().max to get just one proc id
  this->comm().max(root_id);

  // Then broadcast the result
  this->comm().broadcast(value, root_id);

  return value;
}

void RBEIMConstruction::set_point_locator_tol(Real point_locator_tol)
{
  _point_locator_tol = point_locator_tol;
}

Real RBEIMConstruction::get_point_locator_tol() const
{
  return _point_locator_tol;
}

void RBEIMConstruction::initialize_eim_assembly_objects()
{
  _rb_eim_assembly_objects.clear();
  for(unsigned int i=0; i<get_rb_evaluation().get_n_basis_functions(); i++)
    {
      _rb_eim_assembly_objects.push_back( build_eim_assembly(i).release() );
    }
}

ExplicitSystem& RBEIMConstruction::get_explicit_system()
{
  return get_equation_systems().get_system<ExplicitSystem>(_explicit_system_name);
}

void RBEIMConstruction::load_basis_function(unsigned int i)
{
  LOG_SCOPE("load_basis_function()", "RBEIMConstruction");

  libmesh_assert_less (i, get_rb_evaluation().get_n_basis_functions());

  *get_explicit_system().solution = get_rb_evaluation().get_basis_function(i);

  get_explicit_system().update();
}

void RBEIMConstruction::load_rb_solution()
{
  LOG_SCOPE("load_rb_solution()", "RBEIMConstruction");

  solution->zero();

  if(get_rb_evaluation().RB_solution.size() > get_rb_evaluation().get_n_basis_functions())
    libmesh_error_msg("ERROR: System contains " << get_rb_evaluation().get_n_basis_functions() << " basis functions." \
                      << " RB_solution vector constains " << get_rb_evaluation().RB_solution.size() << " entries." \
                      << " RB_solution in RBConstruction::load_rb_solution is too long!");

  for (unsigned int i=0; i<get_rb_evaluation().RB_solution.size(); i++)
    get_explicit_system().solution->add(get_rb_evaluation().RB_solution(i),
                                        get_rb_evaluation().get_basis_function(i));

  get_explicit_system().update();
}

std::vector<ElemAssembly *> RBEIMConstruction::get_eim_assembly_objects()
{
  return _rb_eim_assembly_objects;
}

void RBEIMConstruction::enrich_RB_space()
{
  LOG_SCOPE("enrich_RB_space()", "RBEIMConstruction");

  // put solution in _ghosted_meshfunction_vector so we can access it from the mesh function
  // this allows us to compute EIM_rhs appropriately
  get_explicit_system().solution->localize(*_ghosted_meshfunction_vector,
                                           get_explicit_system().get_dof_map().get_send_list());

  RBEIMEvaluation & eim_eval = cast_ref<RBEIMEvaluation &>(get_rb_evaluation());

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
          get_explicit_system().solution->add(-eim_eval.RB_solution(i), get_rb_evaluation().get_basis_function(i));
        }
    }

  // need to update since context uses current_local_solution
  get_explicit_system().update();

  // Find the quadrature point at which solution (which now stores
  // the "EIM residual") has maximum absolute value
  // by looping over the mesh
  Point optimal_point;
  Number optimal_value = 0.;
  unsigned int optimal_var = 0;
  dof_id_type optimal_elem_id = DofObject::invalid_id;

  // Initialize largest_abs_value to be negative so that it definitely gets updated.
  Real largest_abs_value = -1.;

  // Compute truth representation via projection
  MeshBase & mesh = this->get_mesh();

  UniquePtr<DGFEMContext> explicit_c(new DGFEMContext( get_explicit_system() ));
  DGFEMContext & explicit_context = cast_ref<DGFEMContext &>(*explicit_c);
  init_context_with_sys(explicit_context, get_explicit_system());

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      explicit_context.pre_fe_reinit(get_explicit_system(), *el);
      explicit_context.elem_fe_reinit();

      for(unsigned int var=0; var<get_explicit_system().n_vars(); var++)
        {
          unsigned int n_qpoints = explicit_context.get_element_qrule().n_points();

          for(unsigned int qp=0; qp<n_qpoints; qp++)
            {
              Number value = explicit_context.interior_value(var, qp);
              Real abs_value = std::abs(value);

              if( abs_value > largest_abs_value )
                {
                  optimal_value = value;
                  largest_abs_value = abs_value;
                  optimal_var = var;
                  optimal_elem_id = (*el)->id();

                  FEBase * elem_fe = libmesh_nullptr;
                  explicit_context.get_element_fe( var, elem_fe );
                  optimal_point = elem_fe->get_xyz()[qp];
                }

            }
        }
    }

  // Find out which processor has the largest of the abs values
  unsigned int proc_ID_index;
  this->comm().maxloc(largest_abs_value, proc_ID_index);

  // Broadcast the optimal point from proc_ID_index
  this->comm().broadcast(optimal_point, proc_ID_index);

  // Also broadcast the corresponding optimal_var, optimal_value, and optimal_elem_id
  this->comm().broadcast(optimal_var, proc_ID_index);
  this->comm().broadcast(optimal_value, proc_ID_index);
  this->comm().broadcast(optimal_elem_id, proc_ID_index);

  // In debug mode, assert that we found an optimal_elem_id
  libmesh_assert_not_equal_to(optimal_elem_id, DofObject::invalid_id);

  // Scale the solution
  get_explicit_system().solution->scale(1./optimal_value);

  // Store optimal point in interpolation_points
  eim_eval.interpolation_points.push_back(optimal_point);
  eim_eval.interpolation_points_var.push_back(optimal_var);
  Elem * elem_ptr = mesh.elem_ptr(optimal_elem_id);
  eim_eval.interpolation_points_elem.push_back( elem_ptr );

  NumericVector<Number> * new_bf = NumericVector<Number>::build(this->comm()).release();
  new_bf->init (get_explicit_system().n_dofs(), get_explicit_system().n_local_dofs(), false, PARALLEL);
  *new_bf = *get_explicit_system().solution;
  get_rb_evaluation().basis_functions.push_back( new_bf );

  if(best_fit_type_flag == PROJECTION_BEST_FIT)
    {
      // In order to speed up dot products, we store the product
      // of the basis function and the inner product matrix

      UniquePtr< NumericVector<Number> > implicit_sys_temp1 = this->solution->zero_clone();
      UniquePtr< NumericVector<Number> > implicit_sys_temp2 = this->solution->zero_clone();
      NumericVector<Number>* matrix_times_new_bf =
        get_explicit_system().solution->zero_clone().release();

      // We must localize new_bf before calling get_explicit_sys_subvector
      UniquePtr<NumericVector<Number> > localized_new_bf =
        NumericVector<Number>::build(this->comm());
      localized_new_bf->init(get_explicit_system().n_dofs(), false, SERIAL);
      new_bf->localize(*localized_new_bf);

      for (unsigned int var=0; var<get_explicit_system().n_vars(); var++)
        {
          get_explicit_sys_subvector(*implicit_sys_temp1,
                                     var,
                                     *localized_new_bf);

          inner_product_matrix->vector_mult(*implicit_sys_temp2, *implicit_sys_temp1);

          set_explicit_sys_subvector(*matrix_times_new_bf,
                                     var,
                                     *implicit_sys_temp2);
        }

      _matrix_times_bfs.push_back(matrix_times_new_bf);
    }
}

void RBEIMConstruction::initialize_parametrized_functions_in_training_set()
{
  if(!serial_training_set)
    libmesh_error_msg("Error: We must have serial_training_set==true in " \
                      << "RBEIMConstruction::initialize_parametrized_functions_in_training_set");

  libMesh::out << "Initializing parametrized functions in training set..." << std::endl;
  // initialize rb_eval's parameters
  get_rb_evaluation().initialize_parameters(*this);

  _parametrized_functions_in_training_set.resize( get_n_training_samples() );
  for(unsigned int i=0; i<get_n_training_samples(); i++)
    {
      set_params_from_training_set(i);
      truth_solve(-1);

      _parametrized_functions_in_training_set[i] = get_explicit_system().solution->clone().release();

      libMesh::out << "Completed solve for training sample " << (i+1) << " of " << get_n_training_samples() << std::endl;
    }

  _parametrized_functions_in_training_set_initialized = true;

  libMesh::out << "Parametrized functions in training set initialized" << std::endl << std::endl;
}

void RBEIMConstruction::plot_parametrized_functions_in_training_set(const std::string& pathname)
{
  libmesh_assert(_parametrized_functions_in_training_set_initialized);

  for (std::size_t i=0; i<_parametrized_functions_in_training_set.size(); i++)
    {
#ifdef LIBMESH_HAVE_EXODUS_API
      *get_explicit_system().solution = *_parametrized_functions_in_training_set[i];

      std::stringstream pathname_i;
      pathname_i << pathname << "_" << i << ".exo";

      std::set<std::string> system_names;
      system_names.insert(get_explicit_system().name());
      ExodusII_IO(get_mesh()).write_equation_systems (pathname_i.str(),
                                                      this->get_equation_systems(),
                                                      &system_names);
      libMesh::out << "Plotted parameterized function " << i << std::endl;
#endif
    }
}



Real RBEIMConstruction::compute_best_fit_error()
{
  LOG_SCOPE("compute_best_fit_error()", "RBEIMConstruction");

  const unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  // load up the parametrized function for the current parameters
  truth_solve(-1);

  switch(best_fit_type_flag)
    {
      // Perform an L2 projection in order to find an approximation to solution (from truth_solve above)
    case(PROJECTION_BEST_FIT):
      {
        // We have pre-stored inner_product_matrix * basis_function[i] for each i
        // so we can just evaluate the dot product here.
        DenseVector<Number> best_fit_rhs(RB_size);
        for(unsigned int i=0; i<RB_size; i++)
          {
            best_fit_rhs(i) = get_explicit_system().solution->dot(*_matrix_times_bfs[i]);
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
      libmesh_error_msg("Should not reach here");
    }

  // load the error into solution
  for (unsigned int i=0; i<get_rb_evaluation().get_n_basis_functions(); i++)
    get_explicit_system().solution->add(-get_rb_evaluation().RB_solution(i),
                                        get_rb_evaluation().get_basis_function(i));

  Real best_fit_error = get_explicit_system().solution->linfty_norm();

  return best_fit_error;
}

Real RBEIMConstruction::truth_solve(int plot_solution)
{
  LOG_SCOPE("truth_solve()", "RBEIMConstruction");

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
      *get_explicit_system().solution =
        *_parametrized_functions_in_training_set[training_parameters_found_index];
      get_explicit_system().update(); // put the solution into current_local_solution as well
    }
  // Otherwise, we have to compute the projection
  else
    {
      if(this->n_vars() != 1)
        {
          libmesh_error_msg("The system that we use to perform EIM L2 solves should have one variable");
        }

      RBEIMEvaluation & eim_eval = cast_ref<RBEIMEvaluation &>(get_rb_evaluation());
      eim_eval.set_parameters( get_parameters() );

      // Compute truth representation via L2 projection
      const MeshBase & mesh = this->get_mesh();

      UniquePtr<DGFEMContext> c(new DGFEMContext( *this ));
      DGFEMContext & context = cast_ref<DGFEMContext &>(*c);
      init_context_with_sys(context, *this);

      // First cache all the element data
      std::vector< std::vector< std::vector<Number> > > parametrized_fn_vals(mesh.n_elem());
      std::vector< std::vector<Real> > JxW_values(mesh.n_elem());
      std::vector< std::vector<std::vector<Real> > > phi_values(mesh.n_elem());

      MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

      for ( ; el != end_el; ++el)
        {
          dof_id_type elem_id = (*el)->id();

          context.pre_fe_reinit(*this, *el);
          context.elem_fe_reinit();

          FEBase * elem_fe = libmesh_nullptr;
          context.get_element_fe( 0, elem_fe );
          unsigned int n_qpoints = context.get_element_qrule().n_points();
          const std::vector<Real> & JxW = elem_fe->get_JxW();
          const std::vector<std::vector<Real> > & phi = elem_fe->get_phi();
          const std::vector<Point> & xyz = elem_fe->get_xyz();

          // Loop over qp before var because parametrized functions often use
          // some caching based on qp.
          parametrized_fn_vals[elem_id].resize(n_qpoints);
          JxW_values[elem_id].resize(n_qpoints);
          phi_values[elem_id].resize(n_qpoints);
          for (unsigned int qp=0; qp<n_qpoints; qp++)
            {
              JxW_values[elem_id][qp] = JxW[qp];

              unsigned int n_var_dofs = cast_int<unsigned int>(context.get_dof_indices().size());
              phi_values[elem_id][qp].resize(n_var_dofs);
              for (unsigned int i=0; i != n_var_dofs; i++)
                {
                  phi_values[elem_id][qp][i] = phi[i][qp];
                }

              parametrized_fn_vals[elem_id][qp].resize(get_explicit_system().n_vars());
              for (unsigned int var=0; var<get_explicit_system().n_vars(); var++)
                {
                  Number eval_result = eim_eval.evaluate_parametrized_function(var, xyz[qp], *(*el));
                  parametrized_fn_vals[elem_id][qp][var] = eval_result;
                }
            }
        }

      // We do a distinct solve for each variable in the ExplicitSystem
      for (unsigned int var=0; var<get_explicit_system().n_vars(); var++)
        {
          rhs->zero();

          MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

          for ( ; el != end_el; ++el)
            {
              dof_id_type elem_id = (*el)->id();

              context.pre_fe_reinit(*this, *el);
              //context.elem_fe_reinit(); <--- skip this because we cached all the FE data

              // Loop over qp before var because parametrized functions often use
              // some caching based on qp.
              for (std::size_t qp=0; qp<JxW_values[elem_id].size(); qp++)
                {
                  unsigned int n_var_dofs = phi_values[elem_id][qp].size();

                  Number eval_result = parametrized_fn_vals[elem_id][qp][var];
                  for (unsigned int i=0; i != n_var_dofs; i++)
                    {
                      context.get_elem_residual()(i) +=
                        JxW_values[elem_id][qp] * eval_result * phi_values[elem_id][qp][i];
                    }
                }

              // Apply constraints, e.g. periodic constraints
              this->get_dof_map().constrain_element_vector(context.get_elem_residual(), context.get_dof_indices() );

              // Add element vector to global vector
              rhs->add_vector(context.get_elem_residual(), context.get_dof_indices() );
            }

          // Solve to find the best fit, then solution stores the truth representation
          // of the function to be approximated
          solve_for_matrix_and_rhs(*inner_product_solver, *inner_product_matrix, *rhs);

          if (assert_convergence)
            check_convergence(*inner_product_solver);

          // Now copy the solution to the explicit system's solution.
          set_explicit_sys_subvector(*get_explicit_system().solution, var, *solution);
        }
      get_explicit_system().update();
    }

  if(plot_solution > 0)
    {
#ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(get_mesh()).write_equation_systems ("truth.exo",
                                                      this->get_equation_systems());
#endif
    }

  return 0.;
}

void RBEIMConstruction::init_context_with_sys(FEMContext & c, System & sys)
{
  // default implementation of init_context
  // for compute_best_fit
  for(unsigned int var=0; var<sys.n_vars(); var++)
    {
      FEBase * elem_fe = libmesh_nullptr;
      c.get_element_fe( var, elem_fe );
      elem_fe->get_JxW();
      elem_fe->get_phi();
      elem_fe->get_xyz();
    }
}

void RBEIMConstruction::update_RB_system_matrices()
{
  LOG_SCOPE("update_RB_system_matrices()", "RBEIMConstruction");

  // First, update the inner product matrix
  {
    unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

    UniquePtr< NumericVector<Number> > explicit_sys_temp =
      get_explicit_system().solution->zero_clone();

    UniquePtr< NumericVector<Number> > temp1 = this->solution->zero_clone();
    UniquePtr< NumericVector<Number> > temp2 = this->solution->zero_clone();

    for(unsigned int i=(RB_size-1); i<RB_size; i++)
      {
        for(unsigned int j=0; j<RB_size; j++)
          {
            // We must localize get_rb_evaluation().get_basis_function(j) before calling
            // get_explicit_sys_subvector
            UniquePtr<NumericVector<Number> > localized_basis_function =
              NumericVector<Number>::build(this->comm());
            localized_basis_function->init(get_explicit_system().n_dofs(), false, SERIAL);
            get_rb_evaluation().get_basis_function(j).localize(*localized_basis_function);

            // Compute reduced inner_product_matrix via a series of matvecs
            for(unsigned int var=0; var<get_explicit_system().n_vars(); var++)
              {
                get_explicit_sys_subvector(*temp1, var, *localized_basis_function);
                inner_product_matrix->vector_mult(*temp2, *temp1);
                set_explicit_sys_subvector(*explicit_sys_temp, var, *temp2);
              }

            Number value = explicit_sys_temp->dot( get_rb_evaluation().get_basis_function(i) );
            get_rb_evaluation().RB_inner_product_matrix(i,j) = value;
            if(i!=j)
              {
                // The inner product matrix is assumed
                // to be hermitian
                get_rb_evaluation().RB_inner_product_matrix(j,i) = libmesh_conj(value);
              }
          }
      }
  }

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  RBEIMEvaluation & eim_eval = cast_ref<RBEIMEvaluation &>(get_rb_evaluation());

  // update the EIM interpolation matrix
  for(unsigned int j=0; j<RB_size; j++)
    {
      // Sample the basis functions at the
      // new interpolation point
      get_rb_evaluation().get_basis_function(j).localize(*_ghosted_meshfunction_vector,
                                                         get_explicit_system().get_dof_map().get_send_list());

      eim_eval.interpolation_matrix(RB_size-1,j) =
        evaluate_mesh_function( eim_eval.interpolation_points_var[RB_size-1],
                                eim_eval.interpolation_points[RB_size-1] );
    }
}

Real RBEIMConstruction::get_RB_error_bound()
{
  Real best_fit_error = compute_best_fit_error();
  return best_fit_error;
}

void RBEIMConstruction::update_system()
{
  libMesh::out << "Updating RB matrices" << std::endl;
  update_RB_system_matrices();
}

bool RBEIMConstruction::greedy_termination_test(Real abs_greedy_error,
                                                Real initial_error,
                                                int)
{
  if(abs_greedy_error < get_abs_training_tolerance())
    {
      libMesh::out << "Absolute error tolerance reached." << std::endl;
      return true;
    }

  Real rel_greedy_error = abs_greedy_error/initial_error;
  if(rel_greedy_error < get_rel_training_tolerance())
    {
      libMesh::out << "Relative error tolerance reached." << std::endl;
      return true;
    }

  if(get_rb_evaluation().get_n_basis_functions() >= this->get_Nmax())
    {
      libMesh::out << "Maximum number of basis functions reached: Nmax = "
                   << get_Nmax() << "." << std::endl;
      return true;
    }

  return false;
}

void RBEIMConstruction::set_explicit_sys_subvector(NumericVector<Number> & dest,
                                                   unsigned int var,
                                                   NumericVector<Number> & source)
{
  LOG_SCOPE("set_explicit_sys_subvector()", "RBEIMConstruction");

  // For convenience we localize the source vector first to make it easier to
  // copy over (no need to do distinct send/receives).
  UniquePtr<NumericVector<Number> > localized_source =
    NumericVector<Number>::build(this->comm());
  localized_source->init(this->n_dofs(), false, SERIAL);
  source.localize(*localized_source);

  for (std::size_t i=0; i<_dof_map_between_systems[var].size(); i++)
    {
      dof_id_type implicit_sys_dof_index = i;
      dof_id_type explicit_sys_dof_index = _dof_map_between_systems[var][i];

      if ((dest.first_local_index() <= explicit_sys_dof_index) &&
          (explicit_sys_dof_index < dest.last_local_index()))
        dest.set(explicit_sys_dof_index,
                 (*localized_source)(implicit_sys_dof_index));
    }

  dest.close();
}

void RBEIMConstruction::get_explicit_sys_subvector(NumericVector<Number> & dest,
                                                   unsigned int var,
                                                   NumericVector<Number> & localized_source)
{
  LOG_SCOPE("get_explicit_sys_subvector()", "RBEIMConstruction");

  for (std::size_t i=0; i<_dof_map_between_systems[var].size(); i++)
    {
      dof_id_type implicit_sys_dof_index = i;
      dof_id_type explicit_sys_dof_index = _dof_map_between_systems[var][i];

      if ((dest.first_local_index() <= implicit_sys_dof_index) &&
          (implicit_sys_dof_index < dest.last_local_index()))
        dest.set(implicit_sys_dof_index,
                 localized_source(explicit_sys_dof_index));
    }

  dest.close();
}

void RBEIMConstruction::init_dof_map_between_systems()
{
  LOG_SCOPE("init_dof_map_between_systems()", "RBEIMConstruction");

  unsigned int n_vars = get_explicit_system().n_vars();
  unsigned int n_sys_dofs = this->n_dofs();

  _dof_map_between_systems.resize(n_vars);
  for(unsigned int var=0; var<n_vars; var++)
    {
      _dof_map_between_systems[var].resize(n_sys_dofs);
    }

  std::vector<dof_id_type> implicit_sys_dof_indices;
  std::vector<dof_id_type> explicit_sys_dof_indices;

  MeshBase::const_element_iterator       el     = get_mesh().active_elements_begin();
  const MeshBase::const_element_iterator end_el = get_mesh().active_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      this->get_dof_map().dof_indices (elem, implicit_sys_dof_indices);

      const unsigned int n_dofs = implicit_sys_dof_indices.size();

      for(unsigned int var=0; var<n_vars; var++)
        {
          get_explicit_system().get_dof_map().dof_indices (elem, explicit_sys_dof_indices, var);

          libmesh_assert(explicit_sys_dof_indices.size() == n_dofs);

          for(unsigned int i=0; i<n_dofs; i++)
            {
              dof_id_type implicit_sys_dof_index = implicit_sys_dof_indices[i];
              dof_id_type explicit_sys_dof_index = explicit_sys_dof_indices[i];

              _dof_map_between_systems[var][implicit_sys_dof_index] =
                explicit_sys_dof_index;
            }
        }
    }
}

} // namespace libMesh
