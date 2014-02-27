// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// ------------------------------------------------------------
// EigenSystem class definition

class EigenSystem : public System
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  EigenSystem (EquationSystems& es,
               const std::string& name,
               const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~EigenSystem ();

  /**
   * Abstract base class to be used for assembly of sensitivity
   * data for EigenSystem. A user class derived from this class may be used to
   * assemble the sensitivity of system by attaching an object
   * with the method \p attach_eigenproblem_sensitivity_assemble_object.
   */
  class EigenproblemSensitivityAssembly
  {
  public:
    /**
     * Destructor.  Virtual because we will have virtual functions.
     */
    virtual ~EigenproblemSensitivityAssembly () {}
    
    /**
     * Assembly function.  This function will be called
     * to assemble the sensitivity of eigenproblem matrices. 
     * The method provides dA/dp_i and dB/dpi for \par i ^th parameter 
     * in the vector \par parameters.
     *
     * If the routine is not able to provide sensitivity for this parameter,
     * then it should return false, and the system will attempt to use
     * finite differencing.
     */
    virtual bool sensitivity_assemble (const ParameterVector& parameters,
                                       const unsigned int i,
                                       SparseMatrix<Number>* sensitivity_A,
                                       SparseMatrix<Number>* sensitivity_B) = 0;
  };

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
  virtual void clear ();

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit ();

  /**
   * Assembles & solves the eigen system.
   */
  virtual void solve ();

  /**
   * Solves the sensitivity system, for the provided parameters. The return
   * parameters are irrelevant for EigenSystem. Sensitivity of eigenvalues
   * are returned in \p sens.
   *
   * This method is only implemented in some derived classes.
   */
  virtual std::pair<unsigned int, Real>
  sensitivity_solve (const ParameterVector& parameters,
                     std::vector<Number>& sens);

  
  /**
   * Assembles the system matrix.
   */
  virtual void assemble ();

  /*!
   *  Assembles the sensitivity of matrix_A and matrix_B with respect to the 
   *  specified parameter
   */
  virtual void assemble_eigensystem_sensitivity(const ParameterVector& parameters,
                                                const unsigned int p);
  
  /**
   * Returns real and imaginary part of the ith eigenvalue. If the vectors are
   * provided in the function argument through \p vec_re and \p vec_im, this
   * method copies the eigenvector in the given vector(s), else copies the
   * vector to System::solution.
   *
   * Note that with Number = Complex, \p vec_im must be NULL, and for 
   * Number = Real and eigen problem type HEP or GHEP, \p vec_im must be NULL.
   * For Number = Real and eigenproblem type NHEP or GNHEP, the real and imag.
   * parts of the eigenvector are copied to \p vec_re and \p vec_im,
   * respectively. If \p vec_im is not provided, then only the real part will be 
   * copied to either \p vec_re or System::solution depending on the second
   * argument.
   */
  virtual std::pair<Real, Real> get_eigenpair (unsigned int i,
                                               NumericVector<Number>* vec_re = NULL,
                                               NumericVector<Number>* vec_im = NULL);

  /**
   * @returns \p "Eigen".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "Eigen"; }

  /**
   * @returns the number of matrices handled by this system
   */
  virtual unsigned int n_matrices () const;

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
   * Register a user function to use in assembling the system
   * RHS sensitivity. If the routine is unable to provide sensitivity for this
   * parameter, then it should return false.
   */
  void attach_eigenproblem_sensitivity_assemble_function
  (bool fptr(EquationSystems& es,
             const std::string& name,
             const ParameterVector& parameters,
             const unsigned int i,
             SparseMatrix<Number>* sensitivity_A,
             SparseMatrix<Number>* sensitivity_B));
  
  /**
   * Register a user object to use in assembling the system
   * RHS sensitivity.
   */
  void attach_eigenproblem_sensitivity_assemble_object (EigenproblemSensitivityAssembly& assemble);

  /**
   * The system matrix for standard eigenvalue problems.
   */
  SparseMatrix<Number> *matrix_A;

  /**
   * A second system matrix for generalized eigenvalue problems.
   */
  SparseMatrix<Number> *matrix_B;

  /**
   * The EigenSolver, definig which interface, i.e solver
   * package to use.
   */
  AutoPtr<EigenSolver<Number> > eigen_solver;


protected:


  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

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


  /*!
   *   checks if either a user provided function or object is available to calculate
   *   the sensitivity of A & B matrices for this eigenproblem. Returns true
   *   if user provided function/object is able to calculate the sensitivity
   *   for this parameter, otherwise returns false.
   */
  bool user_eigensystem_sensitivity_assemble(const ParameterVector& parameters,
                                             const unsigned int p);

  
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

  /**
   * Function that assembles the sensitivity of eigen_system.
   */
  bool (* _eigenproblem_sensitivity_assemble_system_function) (EquationSystems& es,
                                                               const std::string& name,
                                                               const ParameterVector& parameter,
                                                               const unsigned int i,
                                                               SparseMatrix<Number>* sensitivity_A,
                                                               SparseMatrix<Number>* sensitivity_B);
  
  /**
   * Object that assembles the sensitivity of eigen_system.
   */
  EigenproblemSensitivityAssembly * _eigenproblem_sensitivity_assemble_system_object;

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
