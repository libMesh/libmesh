// $Id: nonlinear.h,v 1.1 2004-01-03 15:37:42 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __nonlinear_h__
#define __nonlinear_h__

// C++ includes

// Local includes
#include "libmesh_common.h"
#include "numeric_vector.h"
#include "linear.h"


/**
 * This is a generic class that defines a nonlinear to be used in a
 * simulation.  A user can define a nonlinear by deriving from this
 * class and implementing certain functions.
 *
 * @author Benjamin S. Kirk, 2003-2004.
 */

// ------------------------------------------------------------
// Nonlinear class definition

template <class T = Linear<> >
class Nonlinear : public T
{
public:
  
  /**
   * Constructor. Requires a reference to a system to be solved.
   */
  Nonlinear (EquationSystems& es);

  /**
   * Constructor.  Requires a referece to the \p EquationSystems object.
   */
  Nonlinear (EquationSystems& es,
	     const std::string& name,
	     const unsigned int number);

  /**
   * Destructor.
   */
  ~Nonlinear ();

  /**
   * Re-implement the solve member to do a fixed number of
   * linear solves
   */
  virtual void solve ();

  /**
   * @returns the maximum number of nonlinear steps to take.
   */
  unsigned int max_nonlinear_steps () const { return _max_nl_steps; }

  /**
   * Sets the maximum number of nonlinear steps to take.
   */
  unsigned int & max_nonlinear_steps () { return _max_nl_steps; }

  /**
   * @returns the nonlinear solver tolerance.
   */
  Real nonlinear_tolerance () const { return _nl_tol; }

  /**
   * Sets the nonlinear solver tolerance.
   */
  Real & nonlinear_tolerance () { return _nl_tol; }

  
private:

  /**
   * The maximum number of nonlinear steps to take.
   */
  unsigned int _max_nl_steps;

  /**
   * The nonlinear solver tolerance.
   */
  Real _nl_tol;
};



// ------------------------------------------------------------
// Nonlinear inline members
template <class T>
Nonlinear<T>::Nonlinear (EquationSystems& es) :
  T             (es),  // Call the base class constructor
  _max_nl_steps (5),   // Default solver attributes
  _nl_tol       (1.e-6)
{
}



template <class T>
Nonlinear<T>::Nonlinear (EquationSystems& es,
			 const std::string& name,
			 const unsigned int number) :
  Nonlinear (es),
  T         (es, name, number)  // Call the base class constructor
{
}



template <class T>
Nonlinear<T>::~Nonlinear ()
{
}



template <class T>
void Nonlinear<T>::solve ()
{
  for (unsigned int l=0; l<this->max_nonlinear_steps(); l++)
    {
      // Get a copy of the solution at the current nonlinear
      // iteration
      AutoPtr<NumericVector<Number> >
	last_nonlinear_soln (this->system().solution->clone());
      
      // Call the base class solver
      T::solve ();

      // Compute the difference between this solution
      // and the last iterate
      last_nonlinear_soln->add (-1., *(this->system().solution));

      // We must close the vector before we ask it for its norm
      last_nonlinear_soln->close();
      
      // Compute the l2 norm of the difference
      const Real norm_delta = last_nonlinear_soln->l2_norm();

      // Print out convergence information
      std::cout << "Nonlinear convergence: ||u - u_old|| = "
		<< norm_delta
		<< std::endl;

      // Terminate the solution iteration if the difference between
      // this iteration and the last is sufficiently small.
      if (norm_delta < this->nonlinear_tolerance())
	{
	  std::cout << " Nonlinear solver converged at step "
		    << l
		    << std::endl;
	  break;
	}
    }
}


#endif // #define __nonlinear_h__
