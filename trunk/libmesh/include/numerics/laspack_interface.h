// $Id: laspack_interface.h,v 1.3 2004-09-22 18:43:01 benkirk Exp $

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



#ifndef __laspack_interface_h__
#define __laspack_interface_h__

#include "libmesh_common.h"

#if defined(HAVE_LASPACK)
//#if defined(HAVE_LASPACK) && !defined(USE_COMPLEX_NUMBERS)


// C++ includes

// Local includes
#include "linear_solver_interface.h"
#include "laspack_vector.h"
#include "laspack_matrix.h"


#include <itersolv.h>
#include <rtc.h>
#include <errhandl.h>



/**
 * This class provides an interface to Laspack
 * iterative solvers that is compatible with the \p libMesh
 * \p LinearSolverInterface<>
 * 
 * @author Benjamin Kirk, 2002-2004
 */
template <typename T>
class LaspackInterface : public LinearSolverInterface<T>
{
 public:
  /**
   *  Constructor. Initializes Laspack data structures
   */
  LaspackInterface ();
    
  /**
   * Destructor.
   */
  ~LaspackInterface ();
  
  /**
   * Release all memory and clear data structures.
   */
  void clear ();

  /**
   * Initialize data structures if not done so already.
   */
  void init ();
  

  /**
   * Call the Laspack solver
   */
    
  std::pair<unsigned int, Real> 
    solve (SparseMatrix<T>  &matrix,
	   NumericVector<T> &solution,
	   NumericVector<T> &rhs,
	   const double tol,
	   const unsigned int m_its);
   
 private:
  
  /**
   * Tells LASPACK to use the user-specified preconditioner stored in
   * \p _preconditioner_type
   */
  void set_laspack_preconditioner_type ();

  /**
   * Preconditioner type
   */
  PrecondProcType _precond_type;
};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
LaspackInterface<T>::LaspackInterface () :
  _precond_type (ILUPrecond)
{
}



template <typename T>
inline
LaspackInterface<T>::~LaspackInterface ()
{
  this->clear ();
}


#endif // #ifdef HAVE_LASPACK
#endif // #ifndef __laspack_interface_h__
