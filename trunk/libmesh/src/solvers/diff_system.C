#include "diff_system.h"
#include "dof_map.h"
#include "numeric_vector.h"
#include "time_solver.h"



DifferentiableSystem::DifferentiableSystem
                      (EquationSystems& es,
                       const std::string& name,
                       const unsigned int number) :
  Parent      (es, name, number),
  compute_internal_sides(false),
  postprocess_sides(false),
  time_solver (NULL),
  time(0.),
  deltat(1.),
  print_solution_norms(false),
  print_solutions(false),
  print_residual_norms(false),
  print_residuals(false),
  print_jacobian_norms(false),
  print_jacobians(false),
  current_local_nonlinear_solution(NumericVector<Number>::build())
{
  untested();
}



DifferentiableSystem::~DifferentiableSystem ()
{
}



void DifferentiableSystem::reinit ()
{
  Parent::reinit();

  // Resize the serial nonlinear solution for the current mesh
  current_local_nonlinear_solution->init (this->n_dofs());

  time_solver->reinit();
}



void DifferentiableSystem::init_data ()
{
  // First, allocate a vector for the iterate in our quasi_Newton
  // solver

  // FIXME - there really ought to be a way to just use the System
  // solution vector if the solver doesn't need an extra vector!

  // We don't want to project more solutions than we have to
//  this->add_vector("_nonlinear_solution", false);
  this->add_vector("_nonlinear_solution");
//  this->project_solution_on_reinit() = false;

  // Next, give us flags for every variable that might be time
  // evolving
  _time_evolving.resize(this->n_vars(), false);

  // Do any initialization our solvers need
  assert (time_solver.get() != NULL);
  time_solver->init();

  // Next initialize ImplicitSystem data
  Parent::init_data();

  // Finally initialize solution/residual/jacobian data structures
  unsigned int n_vars = this->n_vars();

  dof_indices_var.resize(n_vars);

  elem_subsolutions.clear();
  elem_subsolutions.reserve(n_vars);
  elem_subresiduals.clear();
  elem_subresiduals.reserve(n_vars);
  elem_subjacobians.clear();
  elem_subjacobians.resize(n_vars);
  for (unsigned int i=0; i != n_vars; ++i)
    {
      elem_subsolutions.push_back(new DenseSubVector<Number>(elem_solution));
      elem_subresiduals.push_back(new DenseSubVector<Number>(elem_residual));
      elem_subjacobians[i].clear();
      elem_subjacobians[i].reserve(n_vars);

      for (unsigned int j=0; j != n_vars; ++j)
        {
          elem_subjacobians[i].push_back
            (new DenseSubMatrix<Number>(elem_jacobian));
        }
    }

  // Resize the serial nonlinear solution for the current mesh
  current_local_nonlinear_solution->init (this->n_dofs());
}



void DifferentiableSystem::assemble ()
{
  this->assembly(true, true);
}



void DifferentiableSystem::solve ()
{
  time_solver->solve();
}


Number DifferentiableSystem::current_nonlinear_solution (const unsigned int global_dof_number) const
{
  // Check the sizes
  assert (global_dof_number < this->get_dof_map().n_dofs());
  assert (global_dof_number < current_local_nonlinear_solution->size());

  return (*current_local_nonlinear_solution)(global_dof_number);
}
