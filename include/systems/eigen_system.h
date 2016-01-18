// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_EIGEN_SYSTEM_H
#define LIBMESH_EIGEN_SYSTEM_H

#include "libmesh/libmesh_config.h"

// Currently, the EigenSystem should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC)

// Local Includes
#include "libmesh/system.h"
#include "libmesh/eigen_solver.h"

// C++ includes

namespace libMesh
{

// Forward Declarations
template <typename T> class SparseMatrix;


/**
 * This class provides a specific system class.  It aims
 * at solving eigenvalue problems.  Currently, this class
 * is able  to handle standard eigenvalue problems
 * \p A*x=lambda*x  and generalited eigenvalue problems
 * \p A*x=lambda*B*x.
 */
class EigenSystem : public System
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  EigenSystem (EquationSystems & es,
               const std::string & name_in,
               const unsigned int number_in);

  /**
   * Destructor.
   */
  virtual ~EigenSystem ();

  /**
   * The type of system.
   */
  typedef EigenSystem sys_type;

  /**
   * The type of the parent
   */
  typedef System Parent;

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear () libmesh_override;

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit () libmesh_override;

  /**
   * Assembles & solves the eigen system.
   */
  virtual void solve () libmesh_override;

  /**
   * Assembles the system matrix.
   */
  virtual void assemble () libmesh_override;

  /**
   * Returns real and imaginary part of the ith eigenvalue and copies
   * the respective eigen vector to the solution vector.
   */
  virtual std::pair<Real, Real> get_eigenpair (unsigned int i);

  /**
   * @returns \p "Eigen".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const libmesh_override { return "Eigen"; }

  /**
   * @returns the number of matrices handled by this system
   */
  virtual unsigned int n_matrices () const libmesh_override;

  /**
   * @returns the number of converged eigenpairs.
   */
  unsigned int get_n_converged () const {return _n_converged_eigenpairs;}

  /**
   * @returns the number of eigen solver iterations.
   */
  unsigned int get_n_iterations () const {return _n_iterations;}

  /**
   * Sets the type of the current eigen problem.
   */
  void set_eigenproblem_type (EigenProblemType ept);

  /**
   * @returns the eigen problem type.
   */
  EigenProblemType get_eigenproblem_type () const {return _eigen_problem_type;}

  /**
   * @returns true if the underlying problem is generalized
   * , false otherwise.
   */
  bool generalized () const { return _is_generalized_eigenproblem; }

  /**
   * The system matrix for standard eigenvalue problems.
   */
  SparseMatrix<Number> * matrix_A;

  /**
   * A second system matrix for generalized eigenvalue problems.
   */
  SparseMatrix<Number> * matrix_B;

  /**
   * The EigenSolver, definig which interface, i.e solver
   * package to use.
   */
  UniquePtr<EigenSolver<Number> > eigen_solver;


protected:


  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data () libmesh_override;

  /**
   * Initializes the matrices associated with the system
   */
  virtual void init_matrices ();

  /**
   * Set the _n_converged_eigenpairs member, useful for
   * subclasses of EigenSystem.
   */
  void set_n_converged (unsigned int nconv)
  { _n_converged_eigenpairs = nconv; }

  /**
   * Set the _n_iterations member, useful for subclasses of
   * EigenSystem.
   */
  void set_n_iterations (unsigned int its)
  { _n_iterations = its;}


private:

  /**
   * The number of converged eigenpairs.
   */
  unsigned int _n_converged_eigenpairs;

  /**
   * The number of iterations of the eigen solver algorithm.
   */
  unsigned int _n_iterations;

  /**
   * A boolean flag to indicate whether we are dealing with
   * a generalized eigenvalue problem.
   */
  bool _is_generalized_eigenproblem;

  /**
   * The type of the eigenvalue problem.
   */
  EigenProblemType _eigen_problem_type;
};



// ------------------------------------------------------------
// EigenSystem inline methods
inline
unsigned int EigenSystem::n_matrices () const
{
  if(_is_generalized_eigenproblem)
    return 2;

  return 1;
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC

#endif // LIBMESH_EIGEN_SYSTEM_H
