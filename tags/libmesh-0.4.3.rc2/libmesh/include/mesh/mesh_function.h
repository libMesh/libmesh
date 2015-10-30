// $Id: mesh_function.h,v 1.2 2004-01-03 15:37:42 benkirk Exp $

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



#ifndef __mesh_function_h__
#define __mesh_function_h__

// C++ includes
#include <vector>


// Local Includes
#include "function_base.h"
#include "dense_vector.h"



// Forward Declarations
class EquationSystems;
template <typename T> class NumericVector;
class DofMap;
class PointLocatorBase;



/**
 * This class provides function-like objects for data
 * distributed over a mesh.  
 *
 * @author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// MeshFunction class definition
class MeshFunction : public FunctionBase
{
public:

  /**
   * Constructor for mesh based functions with vectors
   * as return value.  Optionally takes a master function.
   */
  MeshFunction (const EquationSystems& eqn_systems,
		const NumericVector<Number>& vec,
		const DofMap& dof_map,
		const std::vector<unsigned int>& vars,
		const FunctionBase* master=NULL);

  /**
   * Constructor for mesh based functions with a number
   * as return value.  Optionally takes a master function.
   */
  MeshFunction (const EquationSystems& eqn_systems,
		const NumericVector<Number>& vec,
		const DofMap& dof_map,
		const unsigned int var,
		const FunctionBase* master=NULL);

  /**
   * Destructor.
   */
  ~MeshFunction ();



  /**
   * The actual initialization process.
   */
  void init ();

  /**
   * @returns the \f$ 0^{th} \f$ entry of the \p std::vector<Number> at point
   * \p p and for \p time, which defaults to zero.  Creates a
   * \p DenseVector<Number> as input to the user-provided method, 
   * so it may be worth thinking about using \p evaluate().
   */
  Number operator() (const Point& p, 
		     const Real time=0.);

  /**
   * Computes values at coordinate \p p and for time \p time, defaults
   * to zero.  It is up to the user-provided method \p _analytical_fptr
   * whether \p output has to have the correct length or should be
   * resized.
   */
  void operator() (const Point& p,
		   const Real time,
		   DenseVector<Number>& output);

protected:


  /**
   * The equation systems handler, from which
   * the data are gathered.
   */
  const EquationSystems& _eqn_systems;

  /**
   * A reference to the vector that holds the data
   * that is to be interpolated.
   */
  const NumericVector<Number>& _vector;

  /**
   * Need access to the \p DofMap of the other system.
   */
  const DofMap& _dof_map;

  /**
   * The indices of the variables within the other system 
   * for which data are to be gathered.
   */
  const std::vector<unsigned int> _system_vars;

  /**
   * A point locator is needed to locate the
   * points in the mesh.
   */
  PointLocatorBase* _point_locator;

};




// ------------------------------------------------------------
// MeshFunction inline methods



#endif

