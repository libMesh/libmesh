// $Id: frequency_system.C,v 1.2 2003-02-13 22:56:08 benkirk Exp $

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
#include "frequency_system.h"


/*
 * For the moment, only PETSc provides complex support
 */
#if defined(USE_COMPLEX_NUMBERS) && defined(HAVE_PETSC)



// ------------------------------------------------------------
// FrequencySystem implementation
FrequencySystem::FrequencySystem (EquationSystems& es,
				  const std::string& name,
				  const SolverPackage solver_package) :
  SystemBase       (es.get_mesh(), name, solver_package),
  init_system_fptr (NULL),
  assemble_fptr    (NULL),
  mass             (SparseMatrix::build (solver_package)),
  damping          (SparseMatrix::build (solver_package)),
  stiffness        (SparseMatrix::build (solver_package)),
  equation_systems (es),
  /* set_frequency() initializes the _solutions vectors */
  _is_assembled    (false),
  _have_freq       (false),
  _n_freq          (0)
{
  std::cerr << "ERROR: Not working (yet)." << std::endl;

  error ();
}



FrequencySystem::~FrequencySystem ()
{
  //init_system_fptr = assemble_fptr = NULL;

  mass->clear ();

  damping->clear ();

  stiffness->clear ();

  clear_frequencies ();
}




void  FrequencySystem::clear_frequencies ()
{
  if (_have_freq)
    {
      for (unsigned int n=0; n<_n_freq; n++)
	_solutions[n]->clear ();

      _frequencies.clear ();

      _have_freq = false;
    }
}




void FrequencySystem::clear ()
{
  SystemBase::clear();

  //init_system_fptr = assemble_fptr = NULL;

  mass->clear ();

  damping->clear ();

  stiffness->clear ();

  _is_assembled = false;

  clear_frequencies ();
}




void FrequencySystem::init ()
{
  // initialize parent data
  SystemBase::init();

  // Possibly call a user-supplied initialization
  // method.
  if (init_system_fptr != NULL)
    {
      init_system_fptr (equation_systems, name());
    }
}



void FrequencySystem::assemble ()
{
  assert (assemble_fptr != NULL);

  if (_is_assembled)  
    {
      std::cerr << "ERROR: Matrices already assembled." << std::endl;
      error (); 
    }

  _is_assembled = true;

  // prepare matrix with the help of the _dof_map, 
  // fill with sparsity pattern
  SystemBase::assemble();


  std::cerr << "ERROR: Have to initialize the mass, damping and" << std::endl
	    << " stiffness matrices somehow. PETSc's MatDuplicate() does not " << std::endl
	    << " work, since the matrix from which to duplicate has to be assembled first"
	    << std::endl;
  error ();

  // ...

  // Call the user-specified matrix assembly function
  assemble_fptr (equation_systems, name());

  //matrix.print ();
  //rhs.print    ();
}



void FrequencySystem::set_frequencies (const Real base_freq,
				       const Real freq_step,
				       const unsigned int n_steps)
{
  if (_have_freq)
    {
      std::cerr << "WARNING: frequencies already initialized, " << std::endl
		<< " have to clear up first." << std::endl;

      clear_frequencies ();
    }

  _have_freq = true;

  _n_freq = n_steps;

  _frequencies.resize (n_steps);

  // set frequencies, build solution storage
  for (unsigned int n=0; n<n_steps; n++)
    {
      AutoPtr<NumericVector> ap (NumericVector::build(equation_systems.get_solver_package()));
      _solutions.insert(_solutions.end(), ap.release() );

      _frequencies[n] = base_freq + n * freq_step;
    }  

}



std::vector< std::pair<unsigned int, Real> >
FrequencySystem::solve ()
{
  // Assemble the linear system, if not already done
  if (_is_assembled != true)
    {
      std::cerr << "WARNING: matrices not assembled yet," << std::endl
		<< " have to assemble first." << std::endl;

      assemble (); 
    }

  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    equation_systems.parameter("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    static_cast<unsigned int>(equation_systems.parameter("linear solver maximum iterations"));

  // return values
  std::vector< std::pair<unsigned int, Real> > vec_rval;


  for (unsigned int n=0; n<_n_freq; n++)
    {
      //Real f = _frequencies[n];
      //Real omega = 2*pi???

      std::cerr << "ERROR: missing how to combine the matrices" << std::endl;
      error ();
      // ...

      // Solve the linear system for this specific frequency
      const std::pair<unsigned int, Real> rval = 
	linear_solver_interface->solve (*matrix, *solution, *rhs, tol, maxits);

      vec_rval.insert(vec_rval.end(), rval);
   }  

  // sanity
  assert (vec_rval.size() == _n_freq);

  return vec_rval; 
}



void FrequencySystem::attach_init_function(void fptr(EquationSystems& es,
						   const std::string& name))
{
  assert (fptr != NULL);
  
  init_system_fptr = fptr;
}



void FrequencySystem::attach_assemble_function(void fptr(EquationSystems& es,
						       const std::string& name))
{
  assert (fptr != NULL);
  
  assemble_fptr = fptr;  
}




#endif // if defined(USE_COMPLEX_NUMBERS) && defined(HAVE_PETSC)
