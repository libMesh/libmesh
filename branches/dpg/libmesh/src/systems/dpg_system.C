// $Id: fem_system.C 4784 2011-08-03 15:29:54Z trumanellis $

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
#include "equation_systems.h"
#include "fe_base.h"
#include "dpg_context.h"
#include "dpg_system.h"
#include "libmesh_logging.h"
#include "mesh_base.h"
#include "numeric_vector.h"
#include "parallel.h"
#include "quadrature.h"
#include "sparse_matrix.h"
#include "time_solver.h"
#include "unsteady_solver.h" // For eulerian_residual

namespace libMesh
{





DPGSystem::DPGSystem (EquationSystems& es,
                      const std::string& name,
                      const unsigned int number)
  : Parent(es, name, number)
{
}


DPGSystem::~DPGSystem ()
{
  this->clear();
}



void DPGSystem::clear()
{
  _test_variables.clear();

  _test_variable_numbers.clear();

  Parent::clear();
}



void DPGSystem::init_data ()
{
  // First initialize LinearImplicitSystem data
  Parent::init_data();
}


void DPGSystem::assembly (bool get_residual, bool get_jacobian)
{
  libmesh_assert(get_residual || get_jacobian);
  std::string log_name;
  if (get_residual && get_jacobian)
    log_name = "assembly()";
  else if (get_residual)
    log_name = "assembly(get_residual)";
  else
    log_name = "assembly(get_jacobian)";

  START_LOG(log_name, "DPGSystem");

  const MeshBase& mesh = this->get_mesh();

  this->update();

  if (print_solution_norms)
    {
      this->solution->close();
      libMesh::out << "|U| = "
                    << this->solution->l1_norm()
                    << std::endl;
    }
  if (print_solutions)
    {
      unsigned int old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
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
  this->init_context(_femcontext);

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
                                << static_cast<unsigned int>(_femcontext.side) << '!' << std::endl;

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
  STOP_LOG(log_name, "DPGSystem");
}



void DPGSystem::solve()
{
  Parent::solve();
}



AutoPtr<DiffContext> DPGSystem::build_context ()
{
  AutoPtr<DiffContext> ap(new DPGContext(*this));

  ap->set_deltat_pointer( &deltat );
  
  return ap;
}



void DPGSystem::init_context(DiffContext &c)
{
  Parent::init_context(c);
}



unsigned int DPGSystem::add_test_variable (const std::string& var,
    const FEType& type,
    const std::set<subdomain_id_type> * const active_subdomains)
{  
  // Make sure the test variable isn't there already
  // or if it is, that it's the type we want
  for (unsigned int v=0; v<this->n_test_vars(); v++)
    if (this->test_variable_name(v) == var)
      {
	if (this->test_variable_type(v) == type)
	  return _test_variables[v].number();

	libMesh::err << "ERROR: incompatible test variable "
		      << var
		      << " has already been added for this system!"
		      << std::endl;
	libmesh_error();
      }

  const unsigned int curr_n_test_vars = this->n_test_vars();						

  // Add the test variable to the list
  _test_variables.push_back((active_subdomains == NULL) ?
		       Variable(var, curr_n_test_vars, type) :
		       Variable(var, curr_n_test_vars, type, *active_subdomains));

  libmesh_assert ((curr_n_test_vars+1) == this->n_test_vars());

  _test_variable_numbers[var] = curr_n_test_vars;

  // // Add the test variable to the _dof_map
  // _dof_map->add_test_variable (_test_variables.back());

  // Return the number of the new variable
  return curr_n_test_vars;
}



unsigned int DPGSystem::add_test_variable (const std::string& var,
    const Order order,
    const FEFamily family,
    const std::set<subdomain_id_type> * const active_subdomains)
{
  return this->add_test_variable(var, 
			    FEType(order, family), 
			    active_subdomains);
}



bool DPGSystem::has_test_variable (const std::string& var) const
{
  return _test_variable_numbers.count(var);
}



unsigned short int DPGSystem::test_variable_number (const std::string& var) const
{
  // Make sure the test variable exists
  std::map<std::string, unsigned short int>::const_iterator
    pos = _test_variable_numbers.find(var);
  
  if (pos == _test_variable_numbers.end())
    {
      libMesh::err << "ERROR: test variable "
		    << var
		    << " does not exist in this system!"
		    << std::endl;      
      libmesh_error();
    }
  libmesh_assert (_test_variables[pos->second].name() == var);

  return pos->second;
}

} // namespace libMesh
