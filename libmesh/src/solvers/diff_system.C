#include "diff_system.h"
#include "time_solver.h"



DifferentiableSystem::DifferentiableSystem
                      (EquationSystems& es,
                       const std::string& name,
                       const unsigned int number) :
  Parent      (es, name, number),
  time_solver (NULL)
{
}



DifferentiableSystem::~DifferentiableSystem ()
{
}



void DifferentiableSystem::init_data ()
{
  // First, allocate a vector for the iterate in our quasi_Newton solver

  // FIXME - there really ought to be a way to just use the System
  // solution vector if the solver doesn't need an extra vector!

  this->add_vector("_nonlinear_solution", false);

  // Next initialize ImplicitSystem data
  Parent::init_data();

}



void DifferentiableSystem::assemble ()
{
  this->assembly(true, true);
}



void DifferentiableSystem::solve ()
{
  assert (time_solver.get() != NULL);
  time_solver->solve();
}
