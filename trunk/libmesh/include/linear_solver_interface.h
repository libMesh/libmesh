// $Id: linear_solver_interface.h,v 1.2 2003-02-13 22:56:07 benkirk Exp $

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



#ifndef __linear_solver_interface_h__
#define __linear_solver_interface_h__


// C++ includes

// Local includes
#include "mesh_common.h"
#include "auto_ptr.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "enum_solver_package.h"
#include "enum_solver_type.h"
#include "enum_preconditioner_type.h"
#include "reference_counted_object.h"


// forward declarations
class LinearSolverInterface;





/**
 * This class provides a deal.II interface to the Solver
 * iterative solver library.
 *
 * @author Benjamin Kirk, 2003
 */

class LinearSolverInterface : public ReferenceCountedObject<LinearSolverInterface>
{
 public:
  
  /**
   *  Constructor. Initializes Solver data structures
   */
  LinearSolverInterface ();
    
  /**
   * Destructor.
   */
  virtual ~LinearSolverInterface ();
  
  /**
   * Builds a \p LinearSolverInterface using the linear solver package specified by
   * \p solver_package
   */
  static AutoPtr<LinearSolverInterface> build(const SolverPackage solver_package);
  
  /**
   * @returns true if the data structures are
   * initialized, false otherwise.
   */
  bool initialized () const { return _is_initialized; }
  
  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () {}

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init () = 0;

  /**
   * Returns the type of solver to use.
   */
  SolverType solver_type () const { return _solver_type; }

  /**
   * Sets the type of solver to use.
   */
  void set_solver_type (const SolverType st)
  { _solver_type = st; }

  /**
   * Returns the type of preconditioner to use.
   */
  PreconditionerType preconditioner_type () const
  { return _preconditioner_type; }

  /**
   * Sets the type of preconditioner to use.
   */
  void set_preconditioner_type (const PreconditionerType  pct)
  { _preconditioner_type = pct; }
  

  /**
   * Call the Solver solver
   */    
  virtual std::pair<unsigned int, Real> 
    solve (SparseMatrix&,
	   NumericVector&,
	   NumericVector&,
	   const double,
	   const unsigned int) = 0;
   
protected:
  
  /**
   * Enum stating which type of iterative solver to use.
   */
  SolverType _solver_type;

  /**
   * Enum statitng with type of preconditioner to use.
   */
  PreconditionerType _preconditioner_type;
  
  /**
   * Flag indicating if the data structures have been initialized.
   */
  bool _is_initialized;
};


/*----------------------- inline functions ----------------------------------*/
inline
LinearSolverInterface::LinearSolverInterface () :
  _solver_type (GMRES),
  _preconditioner_type (ILU_PRECOND),
  _is_initialized (false)
{}



inline
LinearSolverInterface::~LinearSolverInterface ()
{
  clear ();
}



#endif // #ifdef __solver_interface_h__
