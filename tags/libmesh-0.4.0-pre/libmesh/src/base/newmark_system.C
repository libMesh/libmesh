// $Id: newmark_system.C,v 1.4 2003-05-04 23:59:00 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



// C++ includes


// Local includes
#include "newmark_system.h"
#include "equation_systems.h"



// ------------------------------------------------------------
// NewmarkSystem implementation
NewmarkSystem::NewmarkSystem (EquationSystems<NewmarkSystem>& es,
			      const std::string&           name,
			      const unsigned int           number) :
  SystemBase             (es.get_mesh(), name, number),
  init_cond_fptr         (NULL),
  assemble_fptr          (NULL),
  solve_fptr             (NULL),
  _equation_systems      (es)
{
  // default values of the newmark parameters
  _equation_systems.set_parameter("Newmark alpha") = .25;
  _equation_systems.set_parameter("Newmark delta") = .5;

  // add additional matrices and vectors that will be used in the
  // newmark algorithm to the data structure
  // functions SystemBase::add_matrix and SystemBase::add_vector
  // are used so we do not have to bother about initialization and
  // dof mapping

  // system matrices
  this->add_matrix ("stiffness");
  this->add_matrix ("damping");
  this->add_matrix ("mass");

  // load vector  
  this->add_vector ("force");

  // the displacement and the time derivatives  
  this->add_vector ("displacement");
  this->add_vector ("velocity");
  this->add_vector ("acceleration");

  // contributions to the rhs
  this->add_vector ("rhs_m");
  this->add_vector ("rhs_c");
  
  // results from the previous time step
  this->add_vector ("old_solution");
  this->add_vector ("old_acceleration");

  // set bool to false
  _finished_assemble =false;
}



NewmarkSystem::~NewmarkSystem ()
{
  //init_system_fptr = assemble_fptr = NULL;
  // clear the parameters and integration constants
  _equation_systems.unset_parameter("Newmark alpha");
  _equation_systems.unset_parameter("Newmark delta");

}



void NewmarkSystem::clear ()
{
  // use parent clear this will also clear the
  // matrices and vectors added in the constructor
  SystemBase::clear();

  // clear the parameters and integration constants
  _equation_systems.unset_parameter("Newmark alpha");
  _equation_systems.unset_parameter("Newmark delta");


  // set bool to false
  _finished_assemble = false;
}




void NewmarkSystem::init ()
{
  assert (_mesh.is_prepared());
  
  // initialize parent data
  SystemBase::init();
}



void NewmarkSystem::reinit ()
{
  assert (_mesh.is_prepared());
  
  // initialize parent data
  SystemBase::reinit();

  error();
}



void NewmarkSystem::assemble ()
{
  assert (assemble_fptr != NULL);
  assert (!_finished_assemble);

  // prepare matrix with the help of the _dof_map, 
  // fill with sparsity pattern
  SystemBase::assemble();

  // Log how long the user's matrix assembly code takes
  START_LOG("assemble()", "NewmarkSystem");
  
  // Call the user-specified matrix assembly function
  this->assemble_fptr (_equation_systems, this->name());

  // Stop logging the user code
  STOP_LOG("assemble()", "NewmarkSystem");

  // compute the effective system matrix
  this->compute_matrix();

  // apply initial conditions
  this->initial_conditions();

  _finished_assemble = true;

}


void NewmarkSystem::initial_conditions ()
{
  // assert (init_cond_fptr != NULL);

  // Log how long the user's matrix assembly code takes
  START_LOG("initial_conditions ()", "NewmarkSystem");
  
  // Call the user-specified function for initial conditions.
  // If this function is not provided by the user
  // we simply assume at rest conditions and zero the vectors.
  if (init_cond_fptr != NULL)
    this->init_cond_fptr (_equation_systems, this->name());
  else
    { 
      this->get_vector("displacement").zero();
      this->get_vector("velocity").zero();
      this->get_vector("acceleration").zero();
    }

  // Stop logging the user code
  STOP_LOG("initial_conditions ()", "NewmarkSystem");


}


void NewmarkSystem::compute_matrix ()
{

  // zero the global matrix
  this->matrix->zero();

  // close the system matrices
  // this->stiffness->close ();
  // this->damping->close   ();
  // this->mass->close      ();

  // add up the matrices
  this->matrix->add (1., this->get_matrix ("stiffness"));
  this->matrix->add (_a_0, this->get_matrix ("mass"));
  this->matrix->add (_a_1, this->get_matrix ("damping"));

}


void NewmarkSystem::update_rhs ()
{
  // zero the rhs-vector
  NumericVector<Number>& rhs = *this->rhs;
  rhs.zero();

  // get writable references to some vectors
  NumericVector<Number>&  rhs_m          = this->get_vector("rhs_m");
  NumericVector<Number>&  rhs_c          = this->get_vector("rhs_c");


  // zero the vectors for matrix-vector product
  rhs_m.zero();
  rhs_c.zero();

  // compute auxiliary vectors rhs_m and rhs_c
  rhs_m.add(_a_0, this->get_vector("displacement"));
  rhs_m.add(_a_2, this->get_vector("velocity"));
  rhs_m.add(_a_3, this->get_vector("acceleration"));

  rhs_c.add(_a_1, this->get_vector("displacement"));
  rhs_c.add(_a_4, this->get_vector("velocity"));
  rhs_c.add(_a_5, this->get_vector("acceleration"));

  // compute rhs
  rhs.add(this->get_vector("force"));
  rhs.add_vector(rhs_m, this->get_matrix("mass"));
  rhs.add_vector(rhs_c, this->get_matrix("damping"));

}


void NewmarkSystem::update_u_v_a ()
{
  // get some references for convenience
  const NumericVector<Number>&  solu     = *this->solution;

  NumericVector<Number>&  disp_vec       = this->get_vector("displacement");
  NumericVector<Number>&  vel_vec        = this->get_vector("velocity");
  NumericVector<Number>&  acc_vec        = this->get_vector("acceleration");
  NumericVector<Number>&  old_acc        = this->get_vector("old_acceleration");
  NumericVector<Number>&  old_solu       = this->get_vector("old_solution");

  // copy data
  old_solu = disp_vec;
  disp_vec = solu;
  old_acc = acc_vec;

  // compute the new acceleration vector
  acc_vec.scale(-_a_3);
  acc_vec.add(_a_0, disp_vec);
  acc_vec.add(-_a_0, old_solu);
  acc_vec.add(-_a_2,vel_vec);

  // compute the new velocity vector
  vel_vec.add(_a_6,old_acc);
  vel_vec.add(_a_7,acc_vec);

}

std::pair<unsigned int, Real>
NewmarkSystem::solve ()
{

  // make sure that the system is assembled
  if(!_finished_assemble)
    {
      std::cerr << "ERROR: Need to assemble the system before calling solve(). " << std::endl;
      error();
    }


  // Assemble the linear system
  // this->assemble ();

  // Call the user specified pre-solve method
  // if it is provided
  // if(this->solve_fptr != NULL)
  //    this->solve_fptr (_equation_systems, this->name());

  // Log how long the linear solve takes.
  START_LOG("solve()", "NewmarkSystem");
  
  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    _equation_systems.parameter("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    static_cast<unsigned int>(_equation_systems.parameter("linear solver maximum iterations"));

  // Solve the linear system
  const std::pair<unsigned int, Real> rval = 
    linear_solver_interface->solve (*matrix, *solution, *rhs, tol, maxits);

  // Stop logging the linear solve
  STOP_LOG("solve()", "NewmarkSystem");

  return rval; 
}



void NewmarkSystem::set_newmark_parameters (const Real delta_T, 
					    const Real alpha,
					    const Real delta)
{
  assert(delta_T != 0.);

  // the newmark parameters
  _equation_systems.set_parameter("Newmark alpha") = alpha;
  _equation_systems.set_parameter("Newmark delta") = delta;

  // time step size.  should be handled at a later stage through EquationSystems?
  _equation_systems.set_parameter("Newmark time step") = delta_T;

  // the constants for time integration
  _a_0 = 1./(alpha*delta_T*delta_T);
  _a_1 = delta/(alpha*delta_T);
  _a_2 = 1./(alpha*delta_T);
  _a_3 = 1./(2.*alpha)-1.;
  _a_4 = delta/alpha -1.;
  _a_5 = delta_T/2.*(delta/alpha-2.);
  _a_6 = delta_T*(1.-delta);
  _a_7 = delta*delta_T;


}




bool NewmarkSystem::compare (const NewmarkSystem& other_system, 
			     const Real threshold,
			     const bool verbose) const
{
  /* provided the parameters alpha and delta are identical
   * (what EquationSystems will check, since it owns the parameters),
   * we have no additional data to compare, so let SystemBase do the job
   */
  const SystemBase& other_system_base = static_cast<const SystemBase&>(other_system);
  return SystemBase::compare (other_system_base, threshold, verbose);
}




void NewmarkSystem::attach_init_cond_function(void fptr(EquationSystems<NewmarkSystem>& es,
							const std::string& name))
{
  assert (fptr != NULL);
  
  init_cond_fptr = fptr;
}



void NewmarkSystem::attach_assemble_function(void fptr(EquationSystems<NewmarkSystem>& es,
						    const std::string& name))
{
  assert (fptr != NULL);
  
  assemble_fptr = fptr;  
}



void NewmarkSystem::attach_solve_function(void fptr(EquationSystems<NewmarkSystem>& es,
						      const std::string& name))
{
  assert (fptr != NULL);
  
  solve_fptr = fptr;
}
