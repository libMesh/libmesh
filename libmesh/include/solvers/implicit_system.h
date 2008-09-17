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



#ifndef __implicit_system_h__
#define __implicit_system_h__

// C++ includes

// Local Includes
#include "explicit_system.h"

// Forward declarations
template <typename T> class SparseMatrix;



/**
 * This class provides a specific system class.  It aims
 * at implicit systems, offering nothing more than just
 * the essentials needed to solve a system.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the parent class \p ExplicitSystem.
 */

// ------------------------------------------------------------
// ImplicitSystem class definition

class ImplicitSystem : public ExplicitSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  ImplicitSystem (EquationSystems& es,
		  const std::string& name,
		  const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~ImplicitSystem ();

  /**
   * The type of system.
   */
  typedef ImplicitSystem sys_type;

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }

  /**
   * The type of the parent.
   */
  typedef ExplicitSystem Parent;
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  virtual void clear ();

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit ();
   
  /**
   * Prepares \p matrix and \p _dof_map for matrix assembly.
   * Does not actually assemble anything.  For matrix assembly,
   * use the \p assemble() in derived classes.
   * @e Should be overloaded in derived classes.
   */
  virtual void assemble ();
 
//   /**
//    * Assembles & solves the linear system Ax=b. 
//    */
//   virtual void solve ();
 
  /**
   * @returns \p "Implicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "Implicit"; }

  /**
   * Matrix iterator typedefs.
   */
  typedef std::map<std::string, SparseMatrix<Number>* >::iterator        matrices_iterator;
  typedef std::map<std::string, SparseMatrix<Number>* >::const_iterator  const_matrices_iterator;

  /**
   * Adds the additional matrix \p mat_name to this system.  Only
   * allowed @e prior to \p assemble().  All additional matrices
   * have the same sparsity pattern as the matrix used during
   * solution.  When not \p System but the @e user wants to
   * initialize the mayor matrix, then all the additional matrices,
   * if existent, have to be initialized by the user, too.
   */
  SparseMatrix<Number> & add_matrix (const std::string& mat_name);

  /**
   * @returns \p true if this \p System has a matrix associated with the
   * given name, \p false otherwise.
   */
  bool have_matrix (const std::string& mat_name) const;

  /**
   * @returns a const reference to this system's @e additional matrix
   * named \p mat_name.  @e None of these matrices is involved in the 
   * solution process.  Access is only granted when the matrix is already
   * properly initialized.
   */
  const SparseMatrix<Number> & get_matrix (const std::string& mat_name) const;

  /**
   * @returns a writeable reference to this system's @e additional matrix
   * named \p mat_name.  @e None of these matrices is involved in the 
   * solution process.  Access is only granted when the matrix is already
   * properly initialized.
   */
  SparseMatrix<Number> & get_matrix (const std::string& mat_name);

  /**
   * @returns the number of matrices handled by this system
   */
  unsigned int n_matrices () const;
  
  /**
   * The system matrix.  Implicit systems are characterized by
   * the need to solve the linear system Ax=b.  This is the
   * system matrix A.
   */
  SparseMatrix<Number> * matrix;


  
protected:
  
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();
  
  /**
   * Initializes the matrices associated with this system.
   */
  virtual void init_matrices ();


  
private:

  /**
   * Add the system matrix to the \p _matrices data structure.
   * Useful in initialization.
   */
  void add_system_matrix ();

  /**
   * Some systems need an arbitrary number of matrices.
   */
  std::map<std::string, SparseMatrix<Number>* > _matrices;

  /**
   * \p true when additional matrices may still be added, \p false otherwise.
   */
  bool _can_add_matrices;
};



// ------------------------------------------------------------
// ImplicitSystem inline methods
inline
bool ImplicitSystem::have_matrix (const std::string& mat_name) const
{
  return (_matrices.count(mat_name));
}


inline
unsigned int ImplicitSystem::n_matrices () const
{
 return _matrices.size(); 
}


#endif
