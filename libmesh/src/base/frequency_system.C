// $Id: frequency_system.C,v 1.16 2003-07-07 23:19:28 ddreyer Exp $

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



// Local includes
#include "mesh_config.h"

/*
 * Require complex arithmetic
 */
#if defined(USE_COMPLEX_NUMBERS)


// C++ includes
#include <stdio.h>          /* avoid the ostringstream */

// Local includes
#include "frequency_system.h"
#include "equation_systems.h"
#include "mesh_logging.h"
#include "linear_solver_interface.h"



// ------------------------------------------------------------
// FrequencySystem implementation
FrequencySystem::FrequencySystem (EquationSystems& es,
				  const std::string& name,
				  const unsigned int number) :
  SteadySystem              (es, name, number),
  solve_system              (NULL),
  _finished_set_frequencies (false),
  _finished_init            (false),
  _finished_assemble        (false)
{
  // default value for wave speed & fluid density
  //_equation_systems.set_parameter("wave speed") = 340.;
  //_equation_systems.set_parameter("rho")        = 1.225;
}



FrequencySystem::~FrequencySystem ()
{
  this->clear ();
  
  // the additional matrices and vectors are cleared and zero'ed in SystemBase
}




void FrequencySystem::clear ()
{
  SteadySystem::clear();

  _finished_set_frequencies = false;
  _finished_init            = false;
  _finished_assemble        = false;

  // clear frequencies: 
  // 1. in the parameters section of the 
  //    EquationSystems<FrequencySystem> object
  // 2. in the local vector
  if (_equation_systems.parameter_exists ("n_frequencies"))
    {
      for (unsigned int n=0; n < _equation_systems.parameter("n_frequencies"); n++)
	  _equation_systems.unset_parameter(this->form_freq_param_name(n));
      _equation_systems.unset_parameter("current frequency");
    }
}




void FrequencySystem::init_data ()
{
  // Log how long initializing the system takes
  START_LOG("init()", "FrequencySystem");

  // make sure we have frequencies to solve for
  if (!_finished_set_frequencies)
    {
      /*
       * when this system was read from file, check
       * if this has a "n_frequencies" parameter,
       * and initialize us with these.
       */
      if (_equation_systems.parameter_exists ("n_frequencies"))
        {
	  const unsigned int n_freq = 
	      static_cast<unsigned int>(_equation_systems.parameter("n_frequencies"));

	  assert(this->n_frequencies() == 0);
	  assert(n_freq > 0);

	  _finished_set_frequencies = true;

	  this->set_current_frequency(0);
        }
      else
        {
	  std::cerr << "ERROR: Need to set frequencies before calling init(). " << std::endl;
	  error();
	}
    }

  // initialize parent data and additional solution vectors
  SteadySystem::init_data();

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
  SteadySystem::assemble();

//   // Optionally call the user-specified matrix assembly function,
//   // if the user provided it.  Note that \p FrequencySystem also
//   // works without an _assemble_fptr function
//   if (_assemble_fptr != NULL)
//     this->_assemble_fptr (_equation_systems, this->name());

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
  if (_finished_set_frequencies)
    {
      std::cerr << "ERROR: frequencies already initialized." 
		<< std::endl;
      error();
    }

  // store number of frequencies as parameter
  _equation_systems.set_parameter("n_frequencies") = n_freq;

  for (unsigned int n=0; n<n_freq; n++)
    {
      // remember frequencies as parameters
      _equation_systems.set_parameter(this->form_freq_param_name(n)) = 
	  base_freq + n * freq_step;

      // build storage for solution vector
      SystemBase::add_vector(this->form_solu_vec_name(n));
    }  

  _finished_set_frequencies = true;

  // set the current frequency
  this->set_current_frequency(0);
}



void FrequencySystem::set_frequencies_by_range (const Real min_freq,
						const Real max_freq,
						const unsigned int n_freq)
{
  // sanity checks
  assert(max_freq > min_freq);
  assert(n_freq > 0);

  if (_finished_set_frequencies)
    {
      std::cerr << "ERROR: frequencies already initialized. " << std::endl;
      error();
    }

  // store number of frequencies as parameter
  _equation_systems.set_parameter("n_frequencies") = n_freq;

  // set frequencies, build solution storage
  for (unsigned int n=0; n<n_freq; n++)
    {
      // remember frequencies as parameters
      _equation_systems.set_parameter(this->form_freq_param_name(n)) = 
	  min_freq + n*(max_freq-min_freq)/(n_freq-1);
      
      // build storage for solution vector
      SystemBase::add_vector(this->form_solu_vec_name(n));
    }  

  _finished_set_frequencies = true;

  // set the current frequency
  this->set_current_frequency(0);
}



void FrequencySystem::set_frequencies (const std::vector<Real>& frequencies)
{
  // sanity checks
  assert(!frequencies.empty());

  if (_finished_set_frequencies)
    {
      std::cerr << "ERROR: frequencies already initialized. " << std::endl;
      error();
    }

  // store number of frequencies as parameter
  _equation_systems.set_parameter("n_frequencies") = frequencies.size();

  // set frequencies, build solution storage
  for (unsigned int n=0; n<frequencies.size(); n++)
    {
      // remember frequencies as parameters
      _equation_systems.set_parameter(this->form_freq_param_name(n)) = frequencies[n];
      
      // build storage for solution vector
      SystemBase::add_vector(this->form_solu_vec_name(n));
    }  

  _finished_set_frequencies = true;

  // set the current frequency
  this->set_current_frequency(0);
}




unsigned int FrequencySystem::n_frequencies () const
{
  assert(_finished_set_frequencies);
  return (static_cast<unsigned int>(_equation_systems.parameter("n_frequencies")));
}





std::vector< std::pair<unsigned int, Real> >
FrequencySystem::solve (const unsigned int n_start_in, 
			const unsigned int n_stop_in)
{
  // Assemble the linear system, if not already done
  if (!_finished_assemble)
    this->assemble (); 

  // the user-supplied solve method _has_ to be provided by the user
  assert (solve_system != NULL);

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
      
      this->solve_system (_equation_systems, this->name());
      
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



void FrequencySystem::attach_solve_function(void fptr(EquationSystems& es,
						      const std::string& name))
{
  assert (fptr != NULL);
  
  solve_system = fptr;
}



void FrequencySystem::set_current_frequency(unsigned int n)
{
  assert(n < n_frequencies());
  _equation_systems.set_parameter("current frequency") = 
      _equation_systems.parameter(this->form_freq_param_name(n));
}



std::string FrequencySystem::form_freq_param_name(const unsigned int n) const
{
  assert (n < 9999);
  char buf[15];
  sprintf(buf, "frequency %04d", n);
  return (buf);
}



std::string FrequencySystem::form_solu_vec_name(const unsigned int n) const
{
  assert (n < 9999);
  char buf[15];
  sprintf(buf, "solution %04d", n);
  return (buf);
}


#endif // if defined(USE_COMPLEX_NUMBERS)
