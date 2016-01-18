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



#ifndef LIBMESH_NLOPT_OPTIMIZATION_SOLVER_H
#define LIBMESH_NLOPT_OPTIMIZATION_SOLVER_H

#include "libmesh/libmesh_config.h"

// Petsc include files.
#if defined(LIBMESH_HAVE_NLOPT) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)

// Local includes
#include "libmesh/optimization_solver.h"

// NLopt include
#include "nlopt.h"

// C++ includes

namespace libMesh
{

// NLopt callback to set the objective function.
double __libmesh_nlopt_objective(unsigned n,
                                 const double * x,
                                 double * gradient,
                                 void * data);

// NLopt callback to set equality constraints.
void __libmesh_nlopt_equality_constraints(unsigned m,
                                          double * result,
                                          unsigned n,
                                          const double * x,
                                          double * gradient,
                                          void * data);

// NLopt callback to set inequality constraints.
void __libmesh_nlopt_inequality_constraints(unsigned m,
                                            double * result,
                                            unsigned n,
                                            const double * x,
                                            double * gradient,
                                            void * data);

/**
 * This class provides an interface to the NLopt optimization solvers.
 * http://ab-initio.mit.edu/wiki/index.php/NLopt
 *
 * \author David Knezevic
 * \author John Peterson
 * \date 2015
 */
template <typename T>
class NloptOptimizationSolver : public OptimizationSolver<T>
{
public:

  /**
   * The type of system that we use in conjunction with this solver.
   */
  typedef OptimizationSystem sys_type;

  /**
   * Constructor.
   */
  explicit
  NloptOptimizationSolver (sys_type & system);

  /**
   * Destructor.
   */
  ~NloptOptimizationSolver ();

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () libmesh_override;

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init () libmesh_override;

  /**
   * Returns the raw NLopt object.
   */
  nlopt_opt get_nlopt_object() { this->init(); return _opt; }

  /**
   * Call the NLopt solver.
   */
  virtual void solve () libmesh_override;

  /**
   * Prints a useful message about why the latest optimization solve
   * con(di)verged.
   */
  virtual void print_converged_reason() libmesh_override;

  /**
   * Returns the currently-available (or most recently obtained, if the NLopt object has
   * been destroyed) convergence reason.  Refer to NLopt docs for the meaning of different
   * the value.
   */
  virtual int get_converged_reason() libmesh_override;

  /**
   * Returns a writeable reference to the current iteration count
   * which can be incremented in the objective function.
   */
  unsigned & get_iteration_count() { return _iteration_count; }

protected:

  /**
   * Optimization solver context
   */
  nlopt_opt _opt;

  /**
   * Store the result (i.e. convergence/divergence) for the most recent NLopt solve.
   */
  nlopt_result _result;

  /**
   * Stores the current iteration index (incremented at each call of __libmesh_nlopt_objective).
   */
  unsigned _iteration_count;

  /**
   * NLopt requires us to specify a tolerance for the constraints.
   */
  double _constraints_tolerance;

private:

  // Make NLopt callback functions friends
  friend double __libmesh_nlopt_objective (unsigned n,
                                           const double * x,
                                           double * gradient,
                                           void * data);

  friend void __libmesh_nlopt_equality_constraints(unsigned m,
                                                   double * result,
                                                   unsigned n,
                                                   const double * x,
                                                   double * gradient,
                                                   void * data);

  friend void __libmesh_nlopt_inequality_constraints(unsigned m,
                                                     double * result,
                                                     unsigned n,
                                                     const double * x,
                                                     double * gradient,
                                                     void * data);

  // Map between strings and NLopt algorithms for command line parsing.
  static std::map<std::string, nlopt_algorithm> _nlopt_algorithms;

  // Static function used to initialize the _nlopt_algorithms map, see
  // below for its use.  Naming scheme:
  // G/L == global/local optimization
  // N/D == no/yes gradient required
  // See the full list of algorithms at:
  // http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
  static std::map<std::string, nlopt_algorithm> build_map()
  {
    std::map<std::string, nlopt_algorithm> ret;
    ret["LD_SLSQP"]                   = NLOPT_LD_SLSQP;
    ret["LD_MMA"]                     = NLOPT_LD_MMA;
    ret["LD_CCSAQ"]                   = NLOPT_LD_CCSAQ;
    ret["LD_LBFGS"]                   = NLOPT_LD_LBFGS;
    ret["LD_LBFGS_NOCEDAL"]           = NLOPT_LD_LBFGS_NOCEDAL;
    ret["LD_TNEWTON"]                 = NLOPT_LD_TNEWTON;
    ret["LD_TNEWTON_RESTART"]         = NLOPT_LD_TNEWTON_RESTART;
    ret["LD_TNEWTON_PRECOND"]         = NLOPT_LD_TNEWTON_PRECOND;
    ret["LD_TNEWTON_PRECOND_RESTART"] = NLOPT_LD_TNEWTON_PRECOND_RESTART;
    ret["LD_AUGLAG"]                  = NLOPT_LD_AUGLAG;
    ret["LD_VAR1"]                    = NLOPT_LD_VAR1;
    ret["LD_VAR2"]                    = NLOPT_LD_VAR2;
    ret["LN_COBYLA"]                  = NLOPT_LN_COBYLA;
    ret["LN_BOBYQA"]                  = NLOPT_LN_BOBYQA;
    ret["LN_PRAXIS"]                  = NLOPT_LN_PRAXIS;
    ret["LN_NELDERMEAD"]              = NLOPT_LN_NELDERMEAD;
    ret["LN_SBPLX"]                   = NLOPT_LN_SBPLX;
    ret["GN_ISRES"]                   = NLOPT_GN_ISRES;
    return ret;
  }
};



// Call the class-static function to define the class-static member.
// Since it's a template class, you actually do this in the header,
// not the source file.
template <typename T>
std::map<std::string, nlopt_algorithm>
NloptOptimizationSolver<T>::_nlopt_algorithms = NloptOptimizationSolver<T>::build_map();

} // namespace libMesh


#endif // #if defined(LIBMESH_HAVE_NLOPT) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
#endif // LIBMESH_NLOPT_OPTIMIZATION_SOLVER_H
