// $Id: laspack_interface.h,v 1.4 2003-02-20 04:59:58 benkirk Exp $

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



#ifndef __laspack_interface_h__
#define __laspack_interface_h__

#include "mesh_common.h"

#if defined(HAVE_LASPACK) && !defined(USE_COMPLEX_NUMBERS)


// C++ includes

// Local includes
#include "linear_solver_interface.h"
#include "laspack_vector.h"
#include "laspack_matrix.h"


namespace Laspack {
#include <itersolv.h>
#include <rtc.h>
#include <errhandl.h>
}




/**
 * This class provides a deal.II interface to the Laspack
 * iterative solver library.
 * Currently Laspack only supports real datatypes, so
 * this class is a full specialization of \p NumericVector<>
 * with \p Tp = \p Real*
 * @author Benjamin Kirk, 2002
 */

class LaspackInterface : public LinearSolverInterface<Real>
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
    solve (SparseMatrix<Real> &matrix,
	   NumericVector<Real> &solution,
	   NumericVector<Real> &rhs,
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
  Laspack::PrecondProcType _precond_type;
};


/*----------------------- functions ----------------------------------*/
inline
LaspackInterface::LaspackInterface () :
  _precond_type (Laspack::ILUPrecond)
{
}



inline
LaspackInterface::~LaspackInterface ()
{
  clear ();
}


#endif // #ifdef HAVE_LASPACK
#endif // #ifndef __laspack_interface_h__
