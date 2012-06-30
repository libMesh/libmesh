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



#ifndef __fem_context_h__
#define __fem_context_h__

// Local Includes
#include "diff_context.h"
#include "vector_value.h"
#include "fe_type.h"
#include "fe_base.h"

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#include "tensor_value.h"
#endif

// C++ includes
#include <map>

namespace libMesh
{

  // Forward Declarations
  class Elem;
  template <typename T> class FEGenericBase;
  typedef FEGenericBase<Real> FEBase;
  class QBase;
  class Point;
  template <typename T> class NumericVector;

/**
 * This class provides all data required for a physics package
 * (e.g. an FEMSystem subclass) to perform local element residual
 * and jacobian integrations.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2009
 */

// ------------------------------------------------------------
// FEMContext class definition

class FEMContext : public DiffContext
{
public:

  /**
   * Constructor.  Allocates some but fills no data structures.
   */
  explicit
  FEMContext (const System &sys);

  /**
   * Destructor.
   */
  virtual ~FEMContext ();

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   * This API currently present for backward compatibility.
   */
  Number interior_value(unsigned int var, unsigned int qp) const;

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element side.
   * This API currently present for backward compatibility.
   */
  Number side_value(unsigned int var, unsigned int qp) const;

  /**
   * Returns the value of the solution variable \p var at the physical
   * point \p p on the current element.
   * This API currently present for backward compatibility.
   */
  Number point_value(unsigned int var, const Point &p) const;

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   * This API currently present for backward compatibility.
   */
  Gradient interior_gradient(unsigned int var, unsigned int qp) const;

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element side.
   * This API currently present for backward compatibility.
   */
  Gradient side_gradient(unsigned int var, unsigned int qp) const;

  /**
   * Returns the gradient of the solution variable \p var at the physical
   * point \p p on the current element.
   * This API currently present for backward compatibility.
   */
  Gradient point_gradient(unsigned int var, const Point &p) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Returns the hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   * This API currently present for backward compatibility.
   */
  Tensor interior_hessian(unsigned int var, unsigned int qp) const;

  /**
   * Returns the hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element side.
   * This API currently present for backward compatibility.
   */
  Tensor side_hessian(unsigned int var, unsigned int qp) const;

  /**
   * Returns the hessian of the solution variable \p var at the physical
   * point \p p on the current element.
   * This API currently present for backward compatibility.
   */
  Tensor point_hessian(unsigned int var, const Point &p) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Returns the value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   * This API currently present for backward compatibility.
   */
  Number fixed_interior_value(unsigned int var, unsigned int qp) const;

  /**
   * Returns the value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side.
   * This API currently present for backward compatibility.
   */
  Number fixed_side_value(unsigned int var, unsigned int qp) const;

  /**
   * Returns the value of the fixed_solution variable \p var at the physical
   * point \p p on the current element.
   * This API currently present for backward compatibility.
   */
  Number fixed_point_value(unsigned int var, const Point &p) const;

  /**
   * Returns the gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   * This API currently present for backward compatibility.
   */
  Gradient fixed_interior_gradient(unsigned int var, unsigned int qp) const;

  /**
   * Returns the gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side.
   * This API currently present for backward compatibility.
   */
  Gradient fixed_side_gradient(unsigned int var, unsigned int qp) const;

  /**
   * Returns the gradient of the fixed_solution variable \p var at the physical
   * point \p p on the current element.
   * This API currently present for backward compatibility.
   */
  Gradient fixed_point_gradient(unsigned int var, const Point &p) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Returns the hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   * This API currently present for backward compatibility.
   */
  Tensor fixed_interior_hessian(unsigned int var, unsigned int qp) const;

  /**
   * Returns the hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side.
   * This API currently present for backward compatibility.
   */
  Tensor fixed_side_hessian(unsigned int var, unsigned int qp) const;

  /**
   * Returns the hessian of the fixed_solution variable \p var at the physical
   * point \p p on the current element.
   * This API currently present for backward compatibility.
   */
  Tensor fixed_point_hessian (unsigned int var, const Point &p) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Accessor for interior finite element object for variable var.
   */
  template<typename OutputShape>
  void get_element_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const;

  /**
   * Accessor for edge/face (2D/3D) finite element object for variable var.
   */
  template<typename OutputShape>
  void get_side_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const;

  /**
   * Accessor for edge (3D only!) finite element object for variable var.
   */
  template<typename OutputShape>
  void get_edge_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const;

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputShape> 
  void interior_value(unsigned int var, unsigned int qp, OutputShape& u) const;

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputShape> 
  void side_value(unsigned int var, unsigned int qp, OutputShape& u) const;

  /**
   * Returns the value of the solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputShape>
  void point_value(unsigned int var, const Point &p, OutputShape& u) const;

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputShape>
  void interior_gradient(unsigned int var, unsigned int qp, 
			 typename FEGenericBase<OutputShape>::OutputGradient& du) const;

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputShape> 
  void side_gradient(unsigned int var, unsigned int qp, 
		     typename FEGenericBase<OutputShape>::OutputGradient& du) const;

   /**
   * Returns the gradient of the solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputShape>
  void point_gradient(unsigned int var, const Point &p, 
		      typename FEGenericBase<OutputShape>::OutputGradient& grad_u) const;
  
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Returns the hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputShape>
  void interior_hessian(unsigned int var, unsigned int qp,
			typename FEGenericBase<OutputShape>::OutputTensor& d2u) const;

  /**
   * Returns the hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputShape>
  void side_hessian(unsigned int var, unsigned int qp, 
		    typename FEGenericBase<OutputShape>::OutputTensor& d2u) const;

  /**
   * Returns the hessian of the solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputShape>
  void point_hessian(unsigned int var, const Point &p, 
		      typename FEGenericBase<OutputShape>::OutputTensor& hess_u) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Returns the value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputShape>
  void fixed_interior_value(unsigned int var, unsigned int qp, OutputShape& u) const;

  /**
   * Returns the value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputShape>
  void fixed_side_value(unsigned int var, unsigned int qp, OutputShape& u) const;

  /**
   * Returns the value of the fixed_solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputShape>
  void fixed_point_value(unsigned int var, const Point &p, OutputShape& u) const;

  /**
   * Returns the gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputShape>
  void fixed_interior_gradient(unsigned int var, unsigned int qp,
			       typename FEGenericBase<OutputShape>::OutputGradient& grad_u) const;

  /**
   * Returns the gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputShape>
  void fixed_side_gradient(unsigned int var, unsigned int qp,
			   typename FEGenericBase<OutputShape>::OutputGradient& grad_u) const;

  /**
   * Returns the gradient of the fixed_solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputShape>
  void fixed_point_gradient(unsigned int var, const Point &p,
			    typename FEGenericBase<OutputShape>::OutputGradient& grad_u) const;
  
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Returns the hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputShape>
  void fixed_interior_hessian(unsigned int var, unsigned int qp,
			      typename FEGenericBase<OutputShape>::OutputTensor& hess_u) const;
  
  /**
   * Returns the hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputShape>
  void fixed_side_hessian(unsigned int var, unsigned int qp,
			  typename FEGenericBase<OutputShape>::OutputTensor& hess_u) const;
  
  /**
   * Returns the hessian of the fixed_solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputShape>
  void fixed_point_hessian(unsigned int var, const Point &p,
			   typename FEGenericBase<OutputShape>::OutputTensor& hess_u) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Reinitialize all my FEM context data on a given
   * element for the given system
   */
//  virtual void reinit(const FEMSystem&, Elem*);

// should be protected:
  /**
   * Reinitialize Elem and FE objects if necessary for integration at a new
   * point in time: specifically, handle moving elements in moving mesh
   * schemes.
   */
  virtual void elem_reinit(Real theta);

  /**
   * Reinitialize Elem and side FE objects if necessary for integration at a
   * new point in time: specifically, handle moving elements in moving mesh
   * schemes.
   */
  virtual void elem_side_reinit(Real theta);

  /**
   * Reinitialize Elem and edge FE objects if necessary for
   * integration at a new point in time: specifically, handle moving
   * elements in moving mesh schemes.
   */
  virtual void elem_edge_reinit(Real theta);

  /**
   * Reinitializes local data vectors/matrices on the current geometric element
   */
  virtual void pre_fe_reinit(const System&, const Elem *e);

  /**
   * Reinitializes interior FE objects on the current geometric element
   */
  void elem_fe_reinit();

  /**
   * Reinitializes side FE objects on the current geometric element
   */
  void side_fe_reinit();

  /**
   * Reinitializes edge FE objects on the current geometric element
   */
  void edge_fe_reinit();

// should be protected?:
  /**
   * Finite element objects for each variable's interior, sides and edges.
   */
  std::map<FEType, FEBase *> element_fe;
  std::map<FEType, FEBase *> side_fe;
  std::map<FEType, FEBase *> edge_fe;

  /**
   * Pointers to the same finite element objects, but indexed
   * by variable number
   */
  std::vector<FEBase *> element_fe_var;
  std::vector<FEBase *> side_fe_var;
  std::vector<FEBase *> edge_fe_var;

  /**
   * Quadrature rule for element interior.
   * The FEM context will try to find a quadrature rule that
   * correctly integrates all variables
   */
  QBase *element_qrule;

  /**
   * Quadrature rules for element sides
   * The FEM context will try to find a quadrature rule that
   * correctly integrates all variables
   */
  QBase *side_qrule;

  /**
   * Quadrature rules for element edges.  If the FEM context is told
   * to prepare for edge integration on 3D elements, it will try to
   * find a quadrature rule that correctly integrates all variables
   */
  QBase *edge_qrule;

  /**
   * Uses the coordinate data specified by mesh_*_position configuration
   * to set the geometry of \p elem to the value it would take after a fraction
   * \p theta of a timestep.
   */
  void elem_position_set(Real theta);

  /**
   * Uses the geometry of \p elem to set the coordinate data specified
   * by mesh_*_position configuration.
   */
  void elem_position_get();

  /**
   * System from which to acquire moving mesh information
   */
  const System *_mesh_sys;

  /**
   * Variables from which to acquire moving mesh information
   */
  unsigned int _mesh_x_var, _mesh_y_var, _mesh_z_var;

  /**
   * Current element for element_* to examine
   */
  const Elem *elem;

  /**
   * Current side for side_* to examine
   */
  unsigned char side;

  /**
   * Current edge for edge_* to examine
   */
  unsigned char edge;

  /**
   * Cached dimension of elements in this mesh
   */
  unsigned char dim;

protected:

  /**
   * Helper function to reduce some code duplication in the *_point_* methods.
   */
  template<typename OutputShape>
  AutoPtr<FEGenericBase<OutputShape> > build_new_fe( const FEGenericBase<OutputShape>* fe, const Point &p ) const;

  /**
   * Finite element objects for each variable's interior, sides and edges.
   */
  std::map<FEType, FEAbstract*> _element_fe;
  std::map<FEType, FEAbstract*> _side_fe;
  std::map<FEType, FEAbstract*> _edge_fe;
  

  /**
   * Pointers to the same finite element objects, but indexed
   * by variable number
   */
  std::vector<FEAbstract*> _element_fe_var;
  std::vector<FEAbstract*> _side_fe_var;
  std::vector<FEAbstract*> _edge_fe_var;

private:
  /**
   * Uses the coordinate data specified by mesh_*_position configuration
   * to set the geometry of \p elem to the value it would take after a fraction
   * \p theta of a timestep.
   *
   * This does the work of elem_position_set, but isn't safe to call
   * without _mesh_sys/etc. defined first.
   */
  void _do_elem_position_set(Real theta);

  /**
   * Update the time in the context object for the given value of
   * theta, based on the values of "time" and "deltat" stored in the
   * system which created this context.
   */
  void _update_time_from_system(Real theta);
};



// ------------------------------------------------------------
// FEMContext inline methods

inline
void FEMContext::elem_position_set(Real theta)
{
  if (_mesh_sys)
    this->_do_elem_position_set(theta);
}

template<typename OutputShape>
inline
void FEMContext::get_element_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const
{
  libmesh_assert( var < _element_fe_var.size() );
  fe = libmesh_cast_ptr<FEGenericBase<OutputShape>*>( _element_fe_var[var] );
}

template<typename OutputShape>
inline
void FEMContext::get_side_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const
{
  libmesh_assert( var < _side_fe_var.size() );
  fe = libmesh_cast_ptr<FEGenericBase<OutputShape>*>( _side_fe_var[var] );
}

template<typename OutputShape>
inline
void FEMContext::get_edge_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const
{ 
  libmesh_assert( var < _edge_fe_var.size() );
  fe = libmesh_cast_ptr<FEGenericBase<OutputShape>*>( _edge_fe_var[var] );
}


} // namespace libMesh

#endif //__fem_context_h__
