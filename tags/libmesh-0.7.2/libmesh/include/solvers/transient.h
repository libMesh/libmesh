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



#ifndef __transient_h__
#define __transient_h__

// C++ includes

// Local includes
#include "libmesh_common.h"
#include "linear.h"

namespace libMesh
{


/**
 * This is a generic class that defines a transient to be used in a
 * simulation.  A user can define a transient by deriving from this
 * class and implementing certain functions.
 *
 * @author Benjamin S. Kirk, 2003-2004.
 */

// ------------------------------------------------------------
// Transient class definition

template <class T = Linear<> >
class Transient : public T
{
public:
  
  /**
   * Constructor. Requires a reference to a system to be solved.
   */
  Transient (EquationSystems& es);

  /**
   * Constructor.  Requires a referece to the \p EquationSystems object.
   */
  Transient (EquationSystems& es,
	     const std::string& name,
	     const unsigned int number);

  /**
   * Destructor.
   */
  ~Transient ();

  /**
   * Re-implement the solve member to do solve a transient system.
   */
  virtual void solve ();

  /**
   * @returns the current time.
   */
  Real time () const { return _time; }

  /**
   * @returns the target simulation time.
   */
  Real t_end () const { return _t_end; }

  /**
   * Sets the target simulation time.
   */
  Real & t_end () { return _t_end; }
  
  /**
   * @returns the time step used to advance the solution.
   */
  Real dt () const { return _dt; }

  /**
   * Set the time step.
   */
  Real & dt () { return _dt; }

  /**
   * @returns the current time step (0,1,2...)
   */
  unsigned int time_step () const { return _time_step; }

  /**
   * @returns the maximum number of time steps to take.
   */
  unsigned int max_time_steps () const { return _max_time_steps; }

  /**
   * Sets the maximum number of time steps to take.
   */
  unsigned int & max_time_steps () { return _max_time_steps; }

 
protected:

  /**
   * Set the current time.
   * Only to be used by this and derived classes.
   */
  Real & time () { return _time; }

  /**
   * Set the current time step.
   * Only to be used by this and derived classes.
   */
  unsigned int & time_step () { return _time_step; }

  
private:

  /**
   * The current time.
   */
  Real _time;

  /**
   * The time at which the simulation will be terminated.
   */
  Real _t_end;

  /**
   * The time incriment used to advance the solution.
   */
  Real _dt;

  /**
   * The current time step.
   */
  unsigned int _time_step;
  
  /**
   * The maximum number of time steps to take.
   */
  unsigned int _max_time_steps;
};



// ------------------------------------------------------------
// Transient inline members
template <class T>
inline
Transient<T>::Transient(EquationSystems& es) :
  T               (es),   // Call the base class constructor
  _time           (0.),   // Default solver attributes
  _t_end          (1.e20),
  _dt             (1.e-2),
  _time_step      (0),
  _max_time_steps (5)
{
libmesh_deprecated();
}



template <class T>
inline
Transient<T>::Transient (EquationSystems& es,
			 const std::string& name,
			 const unsigned int number) :
  Transient (es),
  T         (es, name, number)
{
libmesh_deprecated();
}



template <class T>
inline
Transient<T>::~Transient ()
{
}



template <class T>
inline
void Transient<T>::solve ()
{
  do
    {
      // Incriment the time step counter
      this->time_step()++;

      // Incriment the time counter
      this->time() += this->dt();
      
      libMesh::out << "Solving time step "
		   << this->time_step()
		   << std::endl;

      // Call the base class solver
      T::solve ();	
    }
  while ((this->time_step() < this->max_time_steps()) &&
	 (this->time()      < this->t_end()));
}

} // namespace libMesh


#endif // #define __transient_h__
