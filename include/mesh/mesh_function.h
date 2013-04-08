// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MESH_FUNCTION_H
#define LIBMESH_MESH_FUNCTION_H

// Local Includes
#include "libmesh/function_base.h"
#include "libmesh/dense_vector.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/tree_base.h"
#include "libmesh/parallel_object.h"

// C++ includes
#include <cstddef>
#include <vector>

namespace libMesh
{


// Forward Declarations
template <typename T> class DenseVector;
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
class MeshFunction : public FunctionBase<Number>,
		     public ParallelObject
{
public:

  /**
   * Constructor for mesh based functions with vectors
   * as return value.  Optionally takes a master function.
   * If the MeshFunction is to be evaluated outside of the
   * local partition of the mesh, then both the mesh in
   * \p eqn_systems and the coefficient vector \p vec
   * should be serialized.
   */
  MeshFunction (const EquationSystems& eqn_systems,
		const NumericVector<Number>& vec,
		const DofMap& dof_map,
		const std::vector<unsigned int>& vars,
		const FunctionBase<Number>* master=NULL);

  /**
   * Constructor for mesh based functions with a number
   * as return value.  Optionally takes a master function.
   * If the MeshFunction is to be evaluated outside of the
   * local partition of the mesh, then both the mesh in
   * \p eqn_systems and the coefficient vector \p vec
   * should be serialized.
   */
  MeshFunction (const EquationSystems& eqn_systems,
		const NumericVector<Number>& vec,
		const DofMap& dof_map,
		const unsigned int var,
		const FunctionBase<Number>* master=NULL);

  /**
   * Destructor.
   */
  ~MeshFunction ();



  /**
   * The actual initialization process.
   * specifies the method to use when building a \p PointLocator
   */
  void init () { this->init(Trees::NODES); };


  /**
   * The actual initialization process.  Takes an optional argument which
   * specifies the method to use when building a \p PointLocator
   */
  virtual void init (const Trees::BuildType point_locator_build_type);

  /**
   * Clears the function.
   */
  virtual void clear ();

  /**
   * Returns a new copy of the function.  The new copy uses the
   * original as a master function to enable simultaneous evaluations
   * of the copies in different threads.
   * Note that this implies the copy should not be used after the
   * original is destroyed.
   */
  virtual AutoPtr<FunctionBase<Number> > clone () const;

  /**
   * @returns the value of variable 0 at point
   * \p p and for \p time, which defaults to zero.
   */
  Number operator() (const Point& p,
		     const Real time=0.);

  /**
   * @returns the first derivatives of variable 0 at point
   * \p p and for \p time, which defaults to zero.
   */
  Gradient gradient (const Point& p,
		     const Real time=0.);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * @returns the second derivatives of variable 0 at point
   * \p p and for \p time, which defaults to zero.
   */
  Tensor hessian (const Point& p,
		  const Real time=0.);
#endif

  /**
   * Computes values at coordinate \p p and for time \p time, which
   * defaults to zero.
   */
  void operator() (const Point& p,
		   const Real time,
		   DenseVector<Number>& output);

  /**
   * Computes gradients at coordinate \p p and for time \p time, which
   * defaults to zero.
   */
  void gradient (const Point& p,
		 const Real time,
		 std::vector<Gradient>& output);

  /**
   * Computes gradients at coordinate \p p and for time \p time, which
   * defaults to zero.
   */
  void hessian (const Point& p,
		const Real time,
		std::vector<Tensor>& output);

  /**
   * Returns the current \p PointLocator object, for you might want to
   * use it elsewhere.  The \p MeshFunction object must be initialized
   * before.
   */
  const PointLocatorBase& get_point_locator (void) const;

  /**
   * Enables out-of-mesh mode.  In this mode, if asked for a point
   * that is not contained in any element, the \p MeshFunction will
   * return the given \p value instead of crashing.  This mode is off
   * per default.  If you use a master mesh function and you want to
   * enable this mode, you will have to enable it for the master mesh
   * function as well and for all mesh functions that have the same
   * master mesh function.  You may, however, specify different
   * values.
   */
  void enable_out_of_mesh_mode(const DenseVector<Number>& value);

  /**
   * Enables out-of-mesh mode.  In this mode, if asked for a point
   * that is not contained in any element, the \p MeshFunction will
   * return the given \p value instead of crashing.  This mode is off
   * per default.  If you use a master mesh function and you want to
   * enable this mode, you will have to enable it for the master mesh
   * function as well and for all mesh functions that have the same
   * master mesh function.  You may, however, specify different
   * values.
   */
  void enable_out_of_mesh_mode(const Number& value);

  /**
   * Disables out-of-mesh mode.  This is also the default.
   */
  void disable_out_of_mesh_mode(void);

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

  /**
   * \p true if out-of-mesh mode is enabled.  See \p
   * enable_out_of_mesh_mode() for more details.  Default is \p false.
   */
  bool _out_of_mesh_mode;

  /**
   * Value to return outside the mesh if out-of-mesh mode is enabled.
   * See \p enable_out_of_mesh_mode() for more details.
   */
  DenseVector<Number> _out_of_mesh_value;
};




// ------------------------------------------------------------
// MeshFunction inline methods


} // namespace libMesh


#endif // LIBMESH_MESH_FUNCTION_H
