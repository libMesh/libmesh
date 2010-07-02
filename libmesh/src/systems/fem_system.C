// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "equation_systems.h"
#include "fem_system.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "numeric_vector.h"
#include "parallel.h"
#include "quadrature.h"
#include "sparse_matrix.h"
#include "time_solver.h"
#include "unsteady_solver.h" // For euler_residual

namespace libMesh
{





FEMSystem::FEMSystem (EquationSystems& es,
                      const std::string& name,
                      const unsigned int number)
  : Parent(es, name, number),
    fe_reinit_during_postprocess(true),
    extra_quadrature_order(0),
    numerical_jacobian_h(1.e-6),
    verify_analytic_jacobians(0.0),
    _mesh_sys(NULL),
    _mesh_x_var(libMesh::invalid_uint),
    _mesh_y_var(libMesh::invalid_uint),
    _mesh_z_var(libMesh::invalid_uint)
{
}


FEMSystem::~FEMSystem ()
{
  this->clear();
}



void FEMSystem::clear()
{
  Parent::clear();
}



void FEMSystem::init_data ()
{
  // First initialize LinearImplicitSystem data
  Parent::init_data();
}


void FEMSystem::assembly (bool get_residual, bool get_jacobian)
{
  libmesh_assert(get_residual || get_jacobian);
  std::string log_name;
  if (get_residual && get_jacobian)
    log_name = "assembly()";
  else if (get_residual)
    log_name = "assembly(get_residual)";
  else
    log_name = "assembly(get_jacobian)";

  START_LOG(log_name, "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

//  this->get_vector("_nonlinear_solution").localize
//    (*current_local_nonlinear_solution,
//     dof_map.get_send_list());
  this->update();

  if (print_solution_norms)
    {
//      this->get_vector("_nonlinear_solution").close();
      this->solution->close();
      libMesh::out << "|U| = "
//                    << this->get_vector("_nonlinear_solution").l1_norm()
                    << this->solution->l1_norm()
                    << std::endl;
    }
  if (print_solutions)
    {
      unsigned int old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
//      libMesh::out << "U = [" << this->get_vector("_nonlinear_solution")
      libMesh::out << "U = [" << *(this->solution)
                    << "];" << std::endl;
      libMesh::out.precision(old_precision);
    }

  // Is this definitely necessary? [RHS]
  if (get_jacobian)
    matrix->zero();
  if (get_residual)
    rhs->zero();

  // Stupid C++ lets you set *Real* verify_analytic_jacobians = true!
  if (verify_analytic_jacobians > 0.5)
    {
      libMesh::err << "WARNING!  verify_analytic_jacobians was set "
                    << "to absurdly large value of "
                    << verify_analytic_jacobians << std::endl;
      libMesh::err << "Resetting to 1e-6!" << std::endl;
      verify_analytic_jacobians = 1e-6;
    }

  // In time-dependent problems, the nonlinear function we're trying
  // to solve at each timestep may depend on the particular solver
  // we're using
  libmesh_assert (time_solver.get() != NULL);

  AutoPtr<DiffContext> con = this->build_context();
  FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);

  // Build the residual and jacobian contributions on every active
  // mesh element on this processor
  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      _femcontext.pre_fe_reinit(*this, *el);
      _femcontext.elem_fe_reinit();

      bool jacobian_computed =
	time_solver->element_residual(get_jacobian, _femcontext);

      // Compute a numeric jacobian if we have to
      if (get_jacobian && !jacobian_computed)
        {
          // Make sure we didn't compute a jacobian and lie about it
          libmesh_assert(_femcontext.elem_jacobian.l1_norm() == 0.0);
          // Logging of numerical jacobians is done separately
          this->numerical_elem_jacobian(_femcontext);
        }

      // Compute a numeric jacobian if we're asked to verify the
      // analytic jacobian we got
      if (get_jacobian && jacobian_computed &&
          verify_analytic_jacobians != 0.0)
        {
          DenseMatrix<Number> analytic_jacobian(_femcontext.elem_jacobian);

          _femcontext.elem_jacobian.zero();
          // Logging of numerical jacobians is done separately
          this->numerical_elem_jacobian(_femcontext);

          Real analytic_norm = analytic_jacobian.l1_norm();
          Real numerical_norm = _femcontext.elem_jacobian.l1_norm();

          // If we can continue, we'll probably prefer the analytic jacobian
          analytic_jacobian.swap(_femcontext.elem_jacobian);

          // The matrix "analytic_jacobian" will now hold the error matrix
          analytic_jacobian.add(-1.0, _femcontext.elem_jacobian);
          Real error_norm = analytic_jacobian.l1_norm();

          Real relative_error = error_norm /
                                std::max(analytic_norm, numerical_norm);

          if (relative_error > verify_analytic_jacobians)
            {
              libMesh::err << "Relative error " << relative_error
                            << " detected in analytic jacobian on element "
                            << _femcontext.elem->id() << '!' << std::endl;

              unsigned int old_precision = libMesh::out.precision();
              libMesh::out.precision(16);
	      libMesh::out << "J_analytic " << _femcontext.elem->id() << " = "
                        << _femcontext.elem_jacobian << std::endl;
              analytic_jacobian.add(1.0, _femcontext.elem_jacobian);
	      libMesh::out << "J_numeric " << _femcontext.elem->id() << " = "
                        << analytic_jacobian << std::endl;

              libMesh::out.precision(old_precision);

              libmesh_error();
            }
        }

      for (_femcontext.side = 0;
           _femcontext.side != _femcontext.elem->n_sides();
           ++_femcontext.side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!compute_internal_sides && 
	      _femcontext.elem->neighbor(_femcontext.side) != NULL)
            continue;

          // Any mesh movement has already been done (and restored,
          // if the TimeSolver isn't broken), but
          // reinitializing the side FE objects is still necessary
          _femcontext.side_fe_reinit();

          DenseMatrix<Number> old_jacobian;
          // If we're in DEBUG mode, we should always verify that the
	  // user's side_residual function doesn't alter our existing
          // jacobian and then lie about it
#ifndef DEBUG
	  // Even if we're not in DEBUG mode, when we're verifying
	  // analytic jacobians we'll want to verify each side's
          // jacobian contribution separately
          if (verify_analytic_jacobians != 0.0 && get_jacobian)
#endif // ifndef DEBUG
            {
              old_jacobian = _femcontext.elem_jacobian;
              _femcontext.elem_jacobian.zero();
            }
	  jacobian_computed =
            time_solver->side_residual(get_jacobian, _femcontext);

          // Compute a numeric jacobian if we have to
          if (get_jacobian && !jacobian_computed)
            {
	      // In DEBUG mode, we've already set elem_jacobian == 0,
	      // so we can make sure side_residual didn't compute a
              // jacobian and lie about it
#ifdef DEBUG
              libmesh_assert(_femcontext.elem_jacobian.l1_norm() == 0.0);
#endif
              // Logging of numerical jacobians is done separately
              this->numerical_side_jacobian(_femcontext);

              // If we're in DEBUG mode or if
	      // verify_analytic_jacobians is on, we've moved
              // elem_jacobian's accumulated values into old_jacobian.
              // Now let's add them back.
#ifndef DEBUG
              if (verify_analytic_jacobians != 0.0)
#endif // ifndef DEBUG
                _femcontext.elem_jacobian += old_jacobian;
            }

          // Compute a numeric jacobian if we're asked to verify the
          // analytic jacobian we got
	  if (get_jacobian && jacobian_computed &&
              verify_analytic_jacobians != 0.0)
            {
	      DenseMatrix<Number> analytic_jacobian(_femcontext.elem_jacobian);

              _femcontext.elem_jacobian.zero();
              // Logging of numerical jacobians is done separately
              this->numerical_side_jacobian(_femcontext);

              Real analytic_norm = analytic_jacobian.l1_norm();
              Real numerical_norm = _femcontext.elem_jacobian.l1_norm();

              // If we can continue, we'll probably prefer the analytic jacobian
              analytic_jacobian.swap(_femcontext.elem_jacobian);

              // The matrix "analytic_jacobian" will now hold the error matrix
              analytic_jacobian.add(-1.0, _femcontext.elem_jacobian);
              Real error_norm = analytic_jacobian.l1_norm();

              Real relative_error = error_norm /
                                    std::max(analytic_norm, numerical_norm);

              if (relative_error > verify_analytic_jacobians)
                {
                  libMesh::err << "Relative error " << relative_error
                                << " detected in analytic jacobian on element "
                                << _femcontext.elem->id()
			        << ", side "
                                << _femcontext.side << '!' << std::endl;

                  unsigned int old_precision = libMesh::out.precision();
                  libMesh::out.precision(16);
	          libMesh::out << "J_analytic " << _femcontext.elem->id() << " = "
                                << _femcontext.elem_jacobian << std::endl;
                  analytic_jacobian.add(1.0, _femcontext.elem_jacobian);
	          libMesh::out << "J_numeric " << _femcontext.elem->id() << " = "
                                << analytic_jacobian << std::endl;
                  libMesh::out.precision(old_precision);

                  libmesh_error();
                }
              // Once we've verified a side, we'll want to add back the
              // rest of the accumulated jacobian
              _femcontext.elem_jacobian += old_jacobian;
            }
	  // In DEBUG mode, we've set elem_jacobian == 0, and we
          // may still need to add the old jacobian back
#ifdef DEBUG
	  if (get_jacobian && jacobian_computed &&
              verify_analytic_jacobians == 0.0)
            {
              _femcontext.elem_jacobian += old_jacobian;
            }
#endif // ifdef DEBUG
        }

#ifdef LIBMESH_ENABLE_AMR
      // We turn off the asymmetric constraint application;
      // enforce_constraints_exactly() should be called in the solver
      if (get_residual && get_jacobian)
        this->get_dof_map().constrain_element_matrix_and_vector
          (_femcontext.elem_jacobian, _femcontext.elem_residual,
           _femcontext.dof_indices, false);
      else if (get_residual)
        this->get_dof_map().constrain_element_vector
          (_femcontext.elem_residual, _femcontext.dof_indices, false);
      else if (get_jacobian)
        this->get_dof_map().constrain_element_matrix
          (_femcontext.elem_jacobian, _femcontext.dof_indices, false);
#endif // #ifdef LIBMESH_ENABLE_AMR

      if (get_jacobian && print_element_jacobians)
        {
          unsigned int old_precision = libMesh::out.precision();
          libMesh::out.precision(16);
	  libMesh::out << "J_elem " << _femcontext.elem->id() << " = "
                    << _femcontext.elem_jacobian << std::endl;
          libMesh::out.precision(old_precision);
        }

      if (get_jacobian)
        this->matrix->add_matrix (_femcontext.elem_jacobian,
                                  _femcontext.dof_indices);
      if (get_residual)
        this->rhs->add_vector (_femcontext.elem_residual,
                               _femcontext.dof_indices);
    }


  if (get_residual && print_residual_norms)
    {
      this->rhs->close();
      libMesh::out << "|F| = " << this->rhs->l1_norm() << std::endl;
    }
  if (get_residual && print_residuals)
    {
      this->rhs->close();
      unsigned int old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      libMesh::out << "F = [" << *(this->rhs) << "];" << std::endl;
      libMesh::out.precision(old_precision);
    }
  if (get_jacobian && print_jacobian_norms)
    {
      this->matrix->close();
      libMesh::out << "|J| = " << this->matrix->l1_norm() << std::endl;
    }
  if (get_jacobian && print_jacobians)
    {
      this->matrix->close();
      unsigned int old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      libMesh::out << "J = [" << *(this->matrix) << "];" << std::endl;
      libMesh::out.precision(old_precision);
    }
  STOP_LOG(log_name, "FEMSystem");
}



void FEMSystem::solve()
{
  Parent::solve();

  // On a moving mesh we want the mesh to reflect the new solution
  this->mesh_position_set();
}



void FEMSystem::mesh_position_set()
{
  // If we don't need to move the mesh, we're done
  if (_mesh_sys != this)
    return;

  const MeshBase& mesh = this->get_mesh();

  // This code won't work on a parallelized mesh yet -
  // it won't get ancestor elements right.
  libmesh_assert(mesh.is_serial());

  // In fact it won't work in parallel at all yet - we don't have
  // current_solution() data for non-ghost elements.  That's going to
  // be put on hold until we're debugged in serial and until I figure
  // out how to better abstract the "ask other procs for their local
  // data, respond to others' queries, work on my query's results"
  // pattern that we seem to be using a lot.
  libmesh_assert(libMesh::n_processors() == 1);
  
  AutoPtr<DiffContext> con = this->build_context();
  FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);

  // Move every mesh element we can
  MeshBase::const_element_iterator el =
    mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_elements_end();

  for ( ; el != end_el; ++el)
    {
      // We need the algebraic data
      _femcontext.pre_fe_reinit(*this, *el);
      // And when asserts are on, we also need the FE so
      // we can assert that the mesh data is of the right type.
#ifndef NDEBUG
      _femcontext.elem_fe_reinit();
#endif

      // This code won't handle moving subactive elements
      libmesh_assert(!_femcontext.elem->has_children());

      _femcontext.elem_position_set(1.);
    }
}



void FEMSystem::postprocess ()
{
  START_LOG("postprocess()", "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

  this->update();

  AutoPtr<DiffContext> con = this->build_context();
  FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);

  // Loop over every active mesh element on this processor
  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      _femcontext.pre_fe_reinit(*this, *el);

      // Optionally initialize all the interior FE objects on elem.
      if (fe_reinit_during_postprocess)
        _femcontext.elem_fe_reinit();
      
      this->element_postprocess(_femcontext);

      for (_femcontext.side = 0;
           _femcontext.side != _femcontext.elem->n_sides();
           ++_femcontext.side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!postprocess_sides ||
              (!compute_internal_sides &&
               _femcontext.elem->neighbor(_femcontext.side) != NULL))
            continue;

          // Optionally initialize all the FE objects on this side.
          if (fe_reinit_during_postprocess)
            _femcontext.side_fe_reinit();

          this->side_postprocess(_femcontext);
        }
    }

  STOP_LOG("postprocess()", "FEMSystem");
}



void FEMSystem::assemble_qoi (const QoISet &qoi_indices)
{
  START_LOG("assemble_qoi()", "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

  this->update();

  AutoPtr<DiffContext> con = this->build_context();
  FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);

  // the quantity of interest is assumed to be a sum of element and
  // side terms
  for (unsigned int i=0; i != qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      qoi[i] = 0;

  // Loop over every active mesh element on this processor
  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      _femcontext.pre_fe_reinit(*this, *el);
      _femcontext.elem_fe_reinit();

      this->element_qoi(_femcontext, qoi_indices);

      for (_femcontext.side = 0;
           _femcontext.side != _femcontext.elem->n_sides();
           ++_femcontext.side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!assemble_qoi_sides ||
              (!compute_internal_sides &&
               _femcontext.elem->neighbor(_femcontext.side) != NULL))
            continue;

          _femcontext.side_fe_reinit();

          this->side_qoi(_femcontext, qoi_indices);
        }
    }

  Parallel::sum(_femcontext.elem_qoi);

  this->qoi = _femcontext.elem_qoi;

  STOP_LOG("assemble_qoi()", "FEMSystem");
}



void FEMSystem::assemble_qoi_derivative (const QoISet& qoi_indices)
{
  START_LOG("assemble_qoi_derivative()", "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

  this->update();

  AutoPtr<DiffContext> con = this->build_context();
  FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);

  // The quantity of interest derivative assembly accumulates on
  // initially zero vectors
  for (unsigned int i=0; i != qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      this->add_adjoint_rhs(i).zero();

  // Loop over every active mesh element on this processor
  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      _femcontext.pre_fe_reinit(*this, *el);
      _femcontext.elem_fe_reinit();

      this->element_qoi_derivative(_femcontext, qoi_indices);

      for (_femcontext.side = 0;
           _femcontext.side != _femcontext.elem->n_sides();
           ++_femcontext.side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!assemble_qoi_sides ||
              (!compute_internal_sides &&
               _femcontext.elem->neighbor(_femcontext.side) != NULL))
            continue;

          _femcontext.side_fe_reinit();

          this->side_qoi_derivative(_femcontext, qoi_indices);
        }

      // We need some unmodified indices to use for constraining
      // multiple vector
      // FIXME - there should be a DofMap::constrain_element_vectors
      // to do this more efficiently
      std::vector<unsigned int> original_dofs = _femcontext.dof_indices;

      for (unsigned int i=0; i != qoi.size(); ++i)
        if (qoi_indices.has_index(i))
          {
            _femcontext.dof_indices = original_dofs;
            this->get_dof_map().constrain_element_vector
              (_femcontext.elem_qoi_derivative[i], _femcontext.dof_indices, false);

            this->get_adjoint_rhs(i).add_vector
              (_femcontext.elem_qoi_derivative[i], _femcontext.dof_indices);
          }
    }

  STOP_LOG("assemble_qoi_derivative()", "FEMSystem");
}



void FEMSystem::numerical_jacobian (TimeSolverResPtr res,
                                    FEMContext &context)
{
  // Logging is done by numerical_elem_jacobian
  // or numerical_side_jacobian

  DenseVector<Number> original_residual(context.elem_residual);
  DenseVector<Number> backwards_residual(context.elem_residual);
  DenseMatrix<Number> numerical_jacobian(context.elem_jacobian);
#ifdef DEBUG
  DenseMatrix<Number> old_jacobian(context.elem_jacobian);
#endif

  Real numerical_point_h = 0.;
  if (_mesh_sys == this)
    numerical_point_h = numerical_jacobian_h * context.elem->hmin();

  for (unsigned int j = 0; j != context.dof_indices.size(); ++j)
    {
      // Take the "minus" side of a central differenced first derivative
      Number original_solution = context.elem_solution(j);
      context.elem_solution(j) -= numerical_jacobian_h;

      // Make sure to catch any moving mesh terms
      // FIXME - this could be less ugly
      Real *coord = NULL;
      if (_mesh_sys == this)
        {
          if (_mesh_x_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != context.dof_indices_var[_mesh_x_var].size(); ++k)
              if (context.dof_indices_var[_mesh_x_var][k] ==
                  context.dof_indices[j])
                coord = &(context.elem->point(k)(0));
          if (_mesh_y_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != context.dof_indices_var[_mesh_y_var].size(); ++k)
              if (context.dof_indices_var[_mesh_y_var][k] == 
                  context.dof_indices[j])
                coord = &(context.elem->point(k)(1));
          if (_mesh_z_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != context.dof_indices_var[_mesh_z_var].size(); ++k)
	      if (context.dof_indices_var[_mesh_z_var][k] ==
                  context.dof_indices[j])
                coord = &(context.elem->point(k)(2));
        }
      if (coord)
        {
          // We have enough information to scale the perturbations
          // here appropriately
          context.elem_solution(j) = original_solution - numerical_point_h;
          *coord = libmesh_real(context.elem_solution(j));
        }

      context.elem_residual.zero();
      ((*time_solver).*(res))(false, context);
#ifdef DEBUG
      libmesh_assert(old_jacobian == context.elem_jacobian);
#endif
      backwards_residual = context.elem_residual;

      // Take the "plus" side of a central differenced first derivative
      context.elem_solution(j) = original_solution + numerical_jacobian_h;
      if (coord)
        {
          context.elem_solution(j) = original_solution + numerical_point_h;
          *coord = libmesh_real(context.elem_solution(j));
        }
      context.elem_residual.zero();
      ((*time_solver).*(res))(false, context);
#ifdef DEBUG
      libmesh_assert(old_jacobian == context.elem_jacobian);
#endif

      context.elem_solution(j) = original_solution;
      if (coord)
        {
          *coord = libmesh_real(context.elem_solution(j));
	  for (unsigned int i = 0; i != context.dof_indices.size(); ++i)
            {
              numerical_jacobian(i,j) =
                (context.elem_residual(i) - backwards_residual(i)) /
                2. / numerical_point_h;
            }
        }
      else
        {
          for (unsigned int i = 0; i != context.dof_indices.size(); ++i)
            {
              numerical_jacobian(i,j) =
                (context.elem_residual(i) - backwards_residual(i)) /
                2. / numerical_jacobian_h;
            }
        }
    }

  context.elem_residual = original_residual;
  context.elem_jacobian = numerical_jacobian;
}



void FEMSystem::numerical_elem_jacobian (FEMContext &context)
{
  START_LOG("numerical_elem_jacobian()", "FEMSystem");
  this->numerical_jacobian(&TimeSolver::element_residual, context);
  STOP_LOG("numerical_elem_jacobian()", "FEMSystem");
}



void FEMSystem::numerical_side_jacobian (FEMContext &context)
{
  START_LOG("numerical_side_jacobian()", "FEMSystem");
  this->numerical_jacobian(&TimeSolver::side_residual, context);
  STOP_LOG("numerical_side_jacobian()", "FEMSystem");
}



AutoPtr<DiffContext> FEMSystem::build_context ()
{
  return AutoPtr<DiffContext>(new FEMContext(*this));
}



void FEMSystem::init_context(DiffContext &c)
{
  Parent::init_context(c);

  FEMContext &context = libmesh_cast_ref<FEMContext&>(c);

  // Make sure we're prepared to do mass integration
  for (unsigned int var = 0; var != this->n_vars(); ++var)
    if (_time_evolving[var])
      {
        context.element_fe_var[var]->get_JxW();
        context.element_fe_var[var]->get_phi();
      }
}



void FEMSystem::time_evolving (unsigned int var)
{
  // Call the parent function
  Parent::time_evolving(var);
}



void FEMSystem::mesh_position_get()
{
  // This function makes no sense unless we've already picked out some
  // variable(s) to reflect mesh position coordinates
  if (!_mesh_sys)
    libmesh_error();

  // We currently assume mesh variables are in our own system
  if (_mesh_sys != this)
    libmesh_not_implemented();

  // Loop over every active mesh element on this processor
  const MeshBase& mesh = this->get_mesh();

  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  AutoPtr<DiffContext> con = this->build_context();
  FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);
  this->init_context(_femcontext);

  // Get the solution's mesh variables from every element
  for ( ; el != end_el; ++el)
    {
      _femcontext.pre_fe_reinit(*this, *el);

      _femcontext.elem_position_get();

      if (_mesh_x_var != libMesh::invalid_uint)
        this->solution->insert(*_femcontext.elem_subsolutions[_mesh_x_var],
                               _femcontext.dof_indices_var[_mesh_x_var]);
      if (_mesh_y_var != libMesh::invalid_uint)
        this->solution->insert(*_femcontext.elem_subsolutions[_mesh_y_var],
                               _femcontext.dof_indices_var[_mesh_y_var]);
      if (_mesh_z_var != libMesh::invalid_uint)
        this->solution->insert(*_femcontext.elem_subsolutions[_mesh_z_var],
                               _femcontext.dof_indices_var[_mesh_z_var]);
    }

  this->solution->close();

  // And make sure the current_local_solution is up to date too
  this->update();
}


bool FEMSystem::eulerian_residual (bool request_jacobian,
                                   DiffContext &c)
{
  // Only calculate a mesh movement residual if it's necessary
  if (!_mesh_sys)
    return request_jacobian;

  FEMContext &context = libmesh_cast_ref<FEMContext&>(c);

  // This function only supports fully coupled mesh motion for now
  libmesh_assert(_mesh_sys == this);

  unsigned int n_qpoints = context.element_qrule->n_points();

  const unsigned int n_x_dofs = (_mesh_x_var == libMesh::invalid_uint) ?
                                0 : context.dof_indices_var[_mesh_x_var].size();
  const unsigned int n_y_dofs = (_mesh_y_var == libMesh::invalid_uint) ?
                                0 : context.dof_indices_var[_mesh_y_var].size();
  const unsigned int n_z_dofs = (_mesh_z_var == libMesh::invalid_uint) ?
                                0 : context.dof_indices_var[_mesh_z_var].size();

  const unsigned int mesh_xyz_var = n_x_dofs ? _mesh_x_var :
                                   (n_y_dofs ? _mesh_y_var :
                                   (n_z_dofs ? _mesh_z_var :
                                    libMesh::invalid_uint));

  // If we're our own _mesh_sys, we'd better be in charge of
  // at least one coordinate, and we'd better have the same
  // FE type for all coordinates we are in charge of
  libmesh_assert(mesh_xyz_var != libMesh::invalid_uint);
  libmesh_assert(!n_x_dofs || context.element_fe_var[_mesh_x_var] ==
                              context.element_fe_var[mesh_xyz_var]);
  libmesh_assert(!n_y_dofs || context.element_fe_var[_mesh_y_var] ==
                              context.element_fe_var[mesh_xyz_var]);
  libmesh_assert(!n_z_dofs || context.element_fe_var[_mesh_z_var] ==
                              context.element_fe_var[mesh_xyz_var]);

  const std::vector<std::vector<Real> >     &psi =
    context.element_fe_var[mesh_xyz_var]->get_phi();

  for (unsigned int var = 0; var != this->n_vars(); ++var)
    {
      // Mesh motion only affects time-evolving variables
      if (!_time_evolving[var])
        continue;

      // The mesh coordinate variables themselves are Lagrangian,
      // not Eulerian, and no convective term is desired.
      if (_mesh_sys == this &&
          (var == _mesh_x_var ||
           var == _mesh_y_var ||
           var == _mesh_z_var))
        continue;

      // Some of this code currently relies on the assumption that
      // we can pull mesh coordinate data from our own system
      if (_mesh_sys != this)
        libmesh_not_implemented();

      // This residual should only be called by unsteady solvers:
      // if the mesh is steady, there's no mesh convection term!
      UnsteadySolver *unsteady =
        dynamic_cast<UnsteadySolver *>(this->time_solver.get());
      if (!unsteady)
        return request_jacobian;

      const std::vector<Real> &JxW = 
        context.element_fe_var[var]->get_JxW();

      const std::vector<std::vector<Real> >     &phi =
        context.element_fe_var[var]->get_phi();

      const std::vector<std::vector<RealGradient> > &dphi =
        context.element_fe_var[var]->get_dphi();

      const unsigned int n_u_dofs = context.dof_indices_var[var].size();

      DenseSubVector<Number> &Fu = *context.elem_subresiduals[var];
      DenseSubMatrix<Number> &Kuu = *context.elem_subjacobians[var][var];

      DenseSubMatrix<Number> *Kux = n_x_dofs ?
        context.elem_subjacobians[var][_mesh_x_var] : NULL;
      DenseSubMatrix<Number> *Kuy = n_y_dofs ?
        context.elem_subjacobians[var][_mesh_y_var] : NULL;
      DenseSubMatrix<Number> *Kuz = n_z_dofs ?
        context.elem_subjacobians[var][_mesh_z_var] : NULL;

      std::vector<Real> delta_x(n_x_dofs, 0.);
      std::vector<Real> delta_y(n_y_dofs, 0.);
      std::vector<Real> delta_z(n_z_dofs, 0.);

      for (unsigned int i = 0; i != n_x_dofs; ++i)
        {
          unsigned int j = context.dof_indices_var[_mesh_x_var][i];
          delta_x[i] = libmesh_real(this->current_solution(j)) -
                       libmesh_real(unsteady->old_nonlinear_solution(j));
        }

      for (unsigned int i = 0; i != n_y_dofs; ++i)
        {
          unsigned int j = context.dof_indices_var[_mesh_y_var][i];
          delta_y[i] = libmesh_real(this->current_solution(j)) -
                       libmesh_real(unsteady->old_nonlinear_solution(j));
        }

      for (unsigned int i = 0; i != n_z_dofs; ++i)
        {
          unsigned int j = context.dof_indices_var[_mesh_z_var][i];
          delta_z[i] = libmesh_real(this->current_solution(j)) -
                       libmesh_real(unsteady->old_nonlinear_solution(j));
        }

      for (unsigned int qp = 0; qp != n_qpoints; ++qp)
        {
          Gradient grad_u = context.interior_gradient(var, qp);
          RealGradient convection(0.);

          for (unsigned int i = 0; i != n_x_dofs; ++i)
	    convection(0) += delta_x[i] * psi[i][qp];
          for (unsigned int i = 0; i != n_y_dofs; ++i)
	    convection(1) += delta_y[i] * psi[i][qp];
          for (unsigned int i = 0; i != n_z_dofs; ++i)
	    convection(2) += delta_z[i] * psi[i][qp];

          for (unsigned int i = 0; i != n_u_dofs; ++i)
            {
              Number JxWxPhiI = JxW[qp] * phi[i][qp];
              Fu(i) += (convection * grad_u) * JxWxPhiI;
              if (request_jacobian)
                {
                  Number JxWxPhiI = JxW[qp] * phi[i][qp];
                  for (unsigned int j = 0; j != n_u_dofs; ++j)
                    Kuu(i,j) += JxWxPhiI * (convection * dphi[j][qp]);

                  Number JxWxPhiIoverDT = JxWxPhiI/this->deltat;

                  Number JxWxPhiIxDUDXoverDT = JxWxPhiIoverDT * grad_u(0);
                  for (unsigned int j = 0; j != n_x_dofs; ++j)
                    (*Kux)(i,j) += JxWxPhiIxDUDXoverDT * psi[j][qp];

                  Number JxWxPhiIxDUDYoverDT = JxWxPhiIoverDT * grad_u(1);
                  for (unsigned int j = 0; j != n_y_dofs; ++j)
                    (*Kuy)(i,j) += JxWxPhiIxDUDYoverDT * psi[j][qp];

                  Number JxWxPhiIxDUDZoverDT = JxWxPhiIoverDT * grad_u(2);
                  for (unsigned int j = 0; j != n_z_dofs; ++j)
                    (*Kuz)(i,j) += JxWxPhiIxDUDZoverDT * psi[j][qp];
                }
            }
        }
    }

  return request_jacobian;
}

      

bool FEMSystem::mass_residual (bool request_jacobian,
                               DiffContext &c)
{
  FEMContext &context = libmesh_cast_ref<FEMContext&>(c);

  unsigned int n_qpoints = context.element_qrule->n_points();

  for (unsigned int var = 0; var != this->n_vars(); ++var)
    {
      if (!_time_evolving[var])
        continue;

      const std::vector<Real> &JxW = 
        context.element_fe_var[var]->get_JxW();

      const std::vector<std::vector<Real> > &phi =
        context.element_fe_var[var]->get_phi();

      const unsigned int n_dofs = context.dof_indices_var[var].size();

      DenseSubVector<Number> &Fu = *context.elem_subresiduals[var];
      DenseSubMatrix<Number> &Kuu = *context.elem_subjacobians[var][var];

      for (unsigned int qp = 0; qp != n_qpoints; ++qp)
        {
          Number u = context.interior_value(var, qp);
          Number JxWxU = JxW[qp] * u;
          for (unsigned int i = 0; i != n_dofs; ++i)
            {
              Fu(i) += JxWxU * phi[i][qp];
              if (request_jacobian && context.elem_solution_derivative)
                {
                  libmesh_assert (context.elem_solution_derivative == 1.0);

                  Number JxWxPhiI = JxW[qp] * phi[i][qp];
                  Kuu(i,i) += JxWxPhiI * phi[i][qp];
                  for (unsigned int j = i+1; j < n_dofs; ++j)
                    {
                      Number Kij = JxWxPhiI * phi[j][qp];
                      Kuu(i,j) += Kij;
                      Kuu(j,i) += Kij;
                    }
                }
            }
        }
    }

  return request_jacobian;
}

} // namespace libMesh
