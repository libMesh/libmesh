// $Id: frequency_system.C,v 1.11 2003-04-11 23:57:05 ddreyer Exp $

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
#include <stdio.h>          /* avoid the ostringstream */

// Local includes
#include "frequency_system.h"
#include "equation_systems.h"


/*
 * Require complex arithmetic
 */
#if defined(USE_COMPLEX_NUMBERS)


// ------------------------------------------------------------
// FrequencySystem implementation
FrequencySystem::FrequencySystem (EquationSystems<FrequencySystem>& es,
				  const std::string&  name,
				  const unsigned int  number,
				  const SolverPackage solver_package) :
  SystemBase                (es.get_mesh(), name, number, solver_package),
  _assemble_fptr            (NULL),
  _solve_fptr               (NULL),
  _equation_systems         (es),
  _finished_set_frequencies (false),
  _finished_init            (false),
  _finished_assemble        (false)
{
  // dumb safety: clear the frequency vector
  _frequencies.clear();

  // initialize additional matrices to store the (inherently) real-valued
  // mass and stiffness as real and imaginary part, respectively, and
  // the damping in another separate matrix
//  SystemBase::add_matrix("mass_stiffness");
//  SystemBase::add_matrix("damping");

  // default value for wave speed & fluid density
  _equation_systems.set_parameter("wave speed") = 340.;
  _equation_systems.set_parameter("rho")        = 1.225;
}



FrequencySystem::~FrequencySystem ()
{
  //_assemble_fptr = _solve_fptr = NULL;

  // clear frequencies: 
  // 1. in the parameters section of the 
  //    EquationSystems<FrequencySystem> object
  // 2. in the local vector
  for (unsigned int n=0; n < this->n_frequencies(); n++)
      _equation_systems.unset_parameter(this->form_freq_param_name(n));
  _equation_systems.unset_parameter("current frequency");

  this->_frequencies.clear();

  // the additional matrices and vectors are cleared and zero'ed in SystemBase
}




void FrequencySystem::clear ()
{
  SystemBase::clear();

  //_assemble_fptr = _solve_fptr = NULL;

  _finished_set_frequencies = false;
  _finished_init            = false;
  _finished_assemble        = false;

  // clear frequencies: 
  // 1. in the parameters section of the 
  //    EquationSystems<FrequencySystem> object
  // 2. in the local vector
  for (unsigned int n=0; n < this->n_frequencies(); n++)
      _equation_systems.unset_parameter(this->form_freq_param_name(n));
  _equation_systems.unset_parameter("current frequency");

  this->_frequencies.clear();
}




void FrequencySystem::init ()
{
  // Log how long initializing the system takes
  START_LOG("init()", "FrequencySystem");

  // make sure we have frequencies to solve for
  if (!_finished_set_frequencies)
    {
      std::cerr << "ERROR: Need to set frequencies before calling init(). " << std::endl;
      error();
    }

  // initialize parent data and additional solution vectors
  SystemBase::init();

  _finished_init = true;

  // Stop logging init()
  STOP_LOG("init()", "FrequencySystem");
}



void FrequencySystem::assemble ()
{
  assert (_finished_init);

  if (_finished_assemble)  
    {
      std::cerr << "ERROR: Matrices already assembled." << std::endl;
      error (); 
    }

  // Log how long assemble() takes
  START_LOG("assemble()", "FrequencySystem");

  // prepare matrix with the help of the _dof_map, 
  // fill with sparsity pattern, initialize the
  // additional matrices
  SystemBase::assemble();

  // Optionally call the user-specified matrix assembly function,
  // if the user provided it.  Note that \p FrequencySystem also
  // works without an _assemble_fptr function
  if (_assemble_fptr != NULL)
    this->_assemble_fptr (_equation_systems, this->name());

  //matrix.print ();
  //rhs.print    ();

  _finished_assemble = true;

  // Log how long assemble() takes
  STOP_LOG("assemble()", "FrequencySystem");
}



void FrequencySystem::set_frequencies_by_steps (const Real base_freq,
						const Real freq_step,
						const unsigned int n_freq)
{
  // sanity check
  assert(this->n_frequencies() == 0);

  if (_finished_set_frequencies)
    {
      std::cerr << "ERROR: frequencies already initialized. " << std::endl;
      error();
    }

  _frequencies.resize (n_freq);

  // store number of frequencies as parameter
  _equation_systems.set_parameter("n_frequencies") = n_freq;
  

  for (unsigned int n=0; n<n_freq; n++)
    {
      // local storage of frequencies
      _frequencies[n] = base_freq + n * freq_step;
      // remember frequencies as parameters, so that they
      // are saved, once the EquationSystems object is written
      _equation_systems.set_parameter(this->form_freq_param_name(n)) = _frequencies[n];

      // build storage for solution vector
      SystemBase::add_vector(this->form_solu_vec_name(n));
    }  

  // set the current frequency
  this->set_current_frequency(0);

  _finished_set_frequencies = true;
}



void FrequencySystem::set_frequencies_by_range (const Real min_freq,
						const Real max_freq,
						const unsigned int n_freq)
{
  // sanity checks
  assert(this->n_frequencies() == 0);
  assert(max_freq > min_freq);
  assert(n_freq > 0);

  if (_finished_set_frequencies)
    {
      std::cerr << "ERROR: frequencies already initialized. " << std::endl;
      error();
    }

  _frequencies.resize (n_freq);

  // store number of frequencies as parameter
  _equation_systems.set_parameter("n_frequencies") = n_freq;

  // set frequencies, build solution storage
  for (unsigned int n=0; n<n_freq; n++)
    {
      // local storage of frequencies
      _frequencies[n] = min_freq + n*(max_freq-min_freq)/(n_freq-1);
      // remember frequencies as parameters, so that they
      // are saved, once the EquationSystems object is written
      _equation_systems.set_parameter(this->form_freq_param_name(n)) = _frequencies[n];
      
      // build storage for solution vector
      SystemBase::add_vector(this->form_solu_vec_name(n));
    }  

  // set the current frequency
  this->set_current_frequency(0);

  _finished_set_frequencies = true;
}



std::vector< std::pair<unsigned int, Real> >
FrequencySystem::solve (const unsigned int n_start_in, 
			const unsigned int n_stop_in)
{
  // Assemble the linear system, if not already done
  if (!_finished_assemble)
    this->assemble (); 

  // the user-supplied solve method _has_ to be provided by the user
  assert (_solve_fptr != NULL);

//  Do not call this, otherwise perflog may count the time twice,
//  due to the additional START/STOP_LOGs below
//  // Log how long solve() takes
//  START_LOG("solve()", "FrequencySystem");


  // default values, solve for whole frequency range
  unsigned int n_start = 0;
  unsigned int n_stop  = this->n_frequencies()-1;

  if ((n_start_in!=static_cast<unsigned int>(-1)) &&
      (n_stop_in !=static_cast<unsigned int>(-1)))
    {
      n_start = n_start_in;
      n_stop  = n_stop_in;
    }
  else if (n_stop ==static_cast<unsigned int>(-1))
    {
      std::cerr << "ERROR: Forgot to set n_stop." << std::endl;
      error();
    }

  // existence & range checks
  assert(n_frequencies() > 0);
  assert(n_stop < n_frequencies());


  // Get the user-specified linear solver tolerance,
  //     the user-specified maximum # of linear solver iterations,
  //     the user-specified wave speed
  const Real tol            =
    _equation_systems.parameter("linear solver tolerance");
  const unsigned int maxits =
    static_cast<unsigned int>(_equation_systems.parameter("linear solver maximum iterations"));


  // return values
  std::vector< std::pair<unsigned int, Real> > vec_rval;

  // start solver loop
  for (unsigned int n=n_start; n<= n_stop; n++)
    {
      // set the current frequency
      this->set_current_frequency(n);

      // Call the user-supplied pre-solve method
      START_LOG("user_pre_solve()", "FrequencySystem");
      
      this->_solve_fptr (_equation_systems, this->name());
      
      STOP_LOG("user_pre_solve()", "FrequencySystem");


      START_LOG("linear_equation_solve()", "FrequencySystem");

      // Solve the linear system for this specific frequency
      const std::pair<unsigned int, Real> rval = 
	linear_solver_interface->solve (*matrix, *solution, *rhs, tol, maxits);

      STOP_LOG("linear_equation_solve()", "FrequencySystem");

      vec_rval.push_back(rval);      

      /**
       * store the current solution in the additional vector
       */
      this->get_vector(this->form_solu_vec_name(n)) = *solution;

    }  

  // sanity check
  assert (vec_rval.size() == (n_stop-n_start+1));

//  Do not call this, otherwise perflog may count the time twice,
//  due to the additional START/STOP_LOGs below
//  // Log how long solve() takes
//  STOP_LOG("solve()", "FrequencySystem");

  return vec_rval; 
}




bool FrequencySystem::compare (const FrequencySystem& other_system, 
			       const Real threshold,
			       const bool verbose) const
{
  // we have no additional data to compare,
  // let SystemBase do the job
  const SystemBase& other_system_base = static_cast<const SystemBase&>(other_system);
  return SystemBase::compare (other_system_base, threshold, verbose);
}




void FrequencySystem::attach_assemble_function(void fptr(EquationSystems<FrequencySystem>& es,
							 const std::string& name))
{
  assert (fptr != NULL);
  
  _assemble_fptr = fptr;  
}




void FrequencySystem::attach_solve_function(void fptr(EquationSystems<FrequencySystem>& es,
						      const std::string& name))
{
  assert (fptr != NULL);
  
  _solve_fptr = fptr;
}



void FrequencySystem::set_current_frequency(unsigned int n)
{
  assert(n < _frequencies.size());
  _equation_systems.set_parameter("current frequency") = _frequencies[n];
}



std::string FrequencySystem::form_freq_param_name(const unsigned int n) const
{
  assert (n < 9999);
  char buf[15];
  sprintf(buf, "frequency_%04d", n);
  return (buf);
}




std::string FrequencySystem::form_solu_vec_name(const unsigned int n) const
{
  assert (n < 9999);
  char buf[15];
  sprintf(buf, "solution_%04d", n);
  return (buf);
}


#endif // if defined(USE_COMPLEX_NUMBERS)
