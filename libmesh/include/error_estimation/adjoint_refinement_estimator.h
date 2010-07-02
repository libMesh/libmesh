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



#ifndef __adjoint_refinement_estimator_h__
#define __adjoint_refinement_estimator_h__

// C++ includes
#include <vector>

// Local Includes
#include "error_estimator.h"
#include "libmesh.h"

namespace libMesh
{

#ifdef LIBMESH_ENABLE_AMR

/**
 * This class implements a ``brute force'' goal-oriented error
 * estimator which computes an estimate of error in a quantity of
 * interest based on the residual of the current coarse grid primal
 * solution as weighted against an adjoint solution on a uniformly
 * refined (in h and/or p, for an arbitrary number of levels) grid.
 *
 * @author Roy H. Stogner, 2009.
 */
class AdjointRefinementEstimator : public ErrorEstimator
{
public:

  /**
   * Constructor.  Sets the most common default parameter values.
   */
  AdjointRefinementEstimator() : number_h_refinements(1),
                                 number_p_refinements(0),
                                 _qoi_set(QoISet())
  {
    // We're not actually going to use error_norm; our norms are
    // absolute values of QoI error.
    error_norm = INVALID_NORM;
  }
  
  /**
   * Destructor.  
   */
  ~AdjointRefinementEstimator() {}

  /**
   * Access to the QoISet (default: weight all QoIs equally) to use
   * when computing errors
   */
  QoISet &qoi_set() { return _qoi_set; }

  /**
   * This function does uniform refinements and an adjoint
   * solve to get an adjoint solution on each cell,
   * then estimates the error by finding the weighted residual
   * of the coarse solution with the fine adjoint solution.
   *
   * system.solve() and system.assembly() must be called, and so
   * should have no side effects.
   *
   * Only the provided system is solved on the refined mesh;
   * we don't support adjoint solves on loosely coupled
   * collections of Systems.
   *
   * The estimated error is output in the vector
   * \p error_per_cell
   */
  virtual void estimate_error (const System& system,
			       ErrorVector& error_per_cell,
			       const NumericVector<Number>* solution_vector = NULL,
			       bool estimate_parent_error = false);

  /**
   * How many h refinements to perform to get the fine grid
   */
  unsigned char number_h_refinements;
  
  /**
   * How many p refinements to perform to get the fine grid
   */
  unsigned char number_p_refinements;

protected:
  /**
   * The code for estimate_error and both estimate_errors versions is very
   * similar, so we use the same function for all three
   */
  virtual void _estimate_error (const EquationSystems *equation_systems,
                                const System* system,
				ErrorVector* error_per_cell,
			        std::map<std::pair<const System*, unsigned int>, ErrorVector*>* errors_per_cell,
			        const std::map<const System*, const NumericVector<Number>* >* solution_vectors = NULL,
				bool estimate_parent_error = false);

  /**
   * A QoISet to handle cases with multiple QoIs available
   */
  QoISet _qoi_set;
};

#endif // #ifdef LIBMESH_ENABLE_AMR

} // namespace libMesh

#endif

