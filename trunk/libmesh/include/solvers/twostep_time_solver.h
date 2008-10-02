// $Id: adaptive_time_solver.h 3039 2008-09-17 01:26:31Z benkirk $

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



#ifndef __twostep_time_solver_h__
#define __twostep_time_solver_h__

// C++ includes

// Local includes
#include "adaptive_time_solver.h"

// Forward declarations
class System;

// UPDATE THIS DESCRIPTION

/**
 * This class wraps another UnsteadySolver derived class, and compares
 * the results of timestepping with deltat and timestepping with
 * 2*deltat to adjust future timestep lengths.
 *
 * Currently this class only works on fully coupled Systems
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2007
 */

// ------------------------------------------------------------
// Solver class definition
class TwostepTimeSolver : public AdaptiveTimeSolver
{
public:
  /**
   * The parent class
   */
  typedef AdaptiveTimeSolver Parent;
  
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  TwostepTimeSolver (sys_type& s);
  
  /**
   * Destructor.
   */
  ~TwostepTimeSolver ();

  void solve();
};



#endif // #define __twostep_time_solver_
