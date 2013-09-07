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



#ifndef LIBMESH_FEM_CONTEXT_H
#define LIBMESH_FEM_CONTEXT_H

// Local Includes
#include "libmesh/diff_context.h"
#include "libmesh/id_types.h"
#include "libmesh/fe_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/vector_value.h"

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#include "libmesh/tensor_value.h"
#endif

// C++ includes
#include <map>

namespace libMesh
{

  // Forward Declarations
  class BoundaryInfo;
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
   * Reports if the boundary id is found on the current side
   */
  bool has_side_boundary_id(boundary_id_type id) const;

  /**
   * Lists the boundary ids found on the current side
   */
  std::vector<boundary_id_type> side_boundary_ids() const;

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
   * Accessor for interior finite element object for scalar-valued variable var.
   */
  FEBase* get_element_fe( unsigned int var ) const;

  /**
   * Accessor for edge/face (2D/3D) finite element object for variable var.
   */
  template<typename OutputShape>
  void get_side_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const;

  /**
   * Accessor for side finite element object for scalar-valued variable var.
   */
  FEBase* get_side_fe( unsigned int var ) const;

  /**
   * Accessor for edge (3D only!) finite element object for variable var.
   */
  template<typename OutputShape>
  void get_edge_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const;

  /**
   * Accessor for edge finite element object for scalar-valued variable var.
   */
  FEBase* get_edge_fe( unsigned int var ) const;

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputType>
  void interior_value(unsigned int var, unsigned int qp,
                      OutputType& u) const;

  /**
   * Fills a vector of values of the _system_vector at the all the quadrature
   * points in the current element interior.
   */
  template<typename OutputType>
  void interior_values(unsigned int var, const NumericVector<Number> & _system_vector,
		       std::vector<OutputType>& interior_values_vector) const;

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputType>
  void side_value(unsigned int var, unsigned int qp,
                  OutputType& u) const;

  /**
   * Fills a vector of values of the _system_vector at the all the quadrature
   * points on the current element side.
   */
  template<typename OutputType>
  void side_values(unsigned int var, const NumericVector<Number> & _system_vector,
		   std::vector<OutputType>& side_values_vector) const;

  /**
   * Returns the value of the solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputType>
  void point_value(unsigned int var, const Point &p,
                   OutputType& u) const;

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputType>
  void interior_gradient(unsigned int var, unsigned int qp,
			 OutputType& du) const;

  /**
   * Fills a vector with the gradient of the solution variable \p var at all the quadrature
   * points in the current element interior. This is the preferred API.
   */
  template<typename OutputType>
  void interior_gradients(unsigned int var, const NumericVector<Number> & _system_vector,
			  std::vector<OutputType>& interior_gradients_vector) const;

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputType>
  void side_gradient(unsigned int var, unsigned int qp,
		     OutputType & du) const;

  /**
   * Fills a vector with the gradient of the solution variable \p var at all the quadrature
   * points on the current element side. This is the preferred API.
   */
  template<typename OutputType>
  void side_gradients(unsigned int var, const NumericVector<Number> & _system_vector,
		      std::vector<OutputType>& side_gradients_vector) const;

   /**
   * Returns the gradient of the solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputType>
  void point_gradient(unsigned int var, const Point &p,
		      OutputType& grad_u) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Returns the hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputType>
  void interior_hessian(unsigned int var, unsigned int qp,
			OutputType& d2u) const;

  /**
   * Fills a vector of hessians of the _system_vector at the all the
   * quadrature points in the current element interior. This is the
   * preferred API.
   */
  template<typename OutputType>
  void interior_hessians(unsigned int var, const NumericVector<Number> & _system_vector,
		         std::vector<OutputType>& d2u_vals) const;

  /**
   * Returns the hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputType>
  void side_hessian(unsigned int var, unsigned int qp,
		    OutputType& d2u) const;

  /**
   * Fills a vector of hessians of the _system_vector at the all the
   * quadrature points on the current element side.  This is the
   * preferred API.
   */
  template<typename OutputType>
  void side_hessians(unsigned int var, const NumericVector<Number> & _system_vector,
		     std::vector<OutputType>& d2u_vals) const;

  /**
   * Returns the hessian of the solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputType>
  void point_hessian(unsigned int var, const Point &p,
		     OutputType& hess_u) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Returns the value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputType>
  void fixed_interior_value(unsigned int var, unsigned int qp,
                            OutputType& u) const;

  /**
   * Returns the value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputType>
  void fixed_side_value(unsigned int var, unsigned int qp,
                        OutputType& u) const;

  /**
   * Returns the value of the fixed_solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputType>
  void fixed_point_value(unsigned int var, const Point &p,
                         OutputType& u) const;

  /**
   * Returns the gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputType>
  void fixed_interior_gradient(unsigned int var, unsigned int qp,
			       OutputType& grad_u) const;

  /**
   * Returns the gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputType>
  void fixed_side_gradient(unsigned int var, unsigned int qp,
			   OutputType& grad_u) const;

  /**
   * Returns the gradient of the fixed_solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputType>
  void fixed_point_gradient(unsigned int var, const Point &p,
			    OutputType& grad_u) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Returns the hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior. This is the preferred API.
   */
  template<typename OutputType>
  void fixed_interior_hessian(unsigned int var, unsigned int qp,
			      OutputType& hess_u) const;

  /**
   * Returns the hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side. This is the preferred API.
   */
  template<typename OutputType>
  void fixed_side_hessian(unsigned int var, unsigned int qp,
			  OutputType& hess_u) const;

  /**
   * Returns the hessian of the fixed_solution variable \p var at the physical
   * point \p p on the current element. This is the preferred API.
   */
  template<typename OutputType>
  void fixed_point_hessian(unsigned int var, const Point &p,
			   OutputType& hess_u) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Returns the curl of the solution variable \p var at the physical
   * point \p p on the current element.
   */
  template<typename OutputType>
  void interior_curl(unsigned int var, unsigned int qp,
		     OutputType& curl_u) const;

  /**
   * Returns the curl of the solution variable \p var at the physical
   * point \p p on the current element.
   */
  template<typename OutputType>
  void point_curl(unsigned int var, const Point &p,
		  OutputType& curl_u) const;

  /**
   * Returns the divergence of the solution variable \p var at the physical
   * point \p p on the current element.
   */
  template<typename OutputType>
  void interior_div(unsigned int var, unsigned int qp,
		    OutputType& div_u) const;

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
  virtual void side_fe_reinit();

  /**
   * Reinitializes edge FE objects on the current geometric element
   */
  void edge_fe_reinit();

  /**
   * Accessor for element interior quadrature rule.
   */
  const QBase& get_element_qrule() const
  { return *(this->element_qrule); }

  /**
   * Accessor for element side quadrature rule.
   */
  const QBase& get_side_qrule() const
  { return *(this->side_qrule); }

  /**
   * Accessor for element edge quadrature rule.
   */
  const QBase& get_edge_qrule() const
  { return *(this->edge_qrule); }

  /**
   * Tells the FEMContext that system \p sys contains the
   * isoparametric Lagrangian variables which correspond to the
   * coordinates of mesh nodes, in problems where the mesh itself is
   * expected to move in time.
   *
   * This should be set automatically if the FEMPhysics requires it.
   */
  virtual void set_mesh_system(System* sys)
  { this->_mesh_sys = sys; }

  /**
   * Accessor for moving mesh System
   */
  const System* get_mesh_system() const
  { return this->_mesh_sys; }

  /**
   * Accessor for moving mesh System
   */
  System* get_mesh_system()
  { return this->_mesh_sys; }

  /**
   * Accessor for x-variable of moving mesh System
   */
  unsigned int get_mesh_x_var() const
  { return _mesh_x_var; }

  /**
   * Accessor for x-variable of moving mesh System
   *
   * This should be set automatically if the FEMPhysics requires it.
   */
  void set_mesh_x_var(unsigned int x_var)
  { _mesh_x_var = x_var; }

  /**
   * Accessor for y-variable of moving mesh System
   */
  unsigned int get_mesh_y_var() const
  { return _mesh_y_var; }

  /**
   * Accessor for y-variable of moving mesh System
   *
   * This should be set automatically if the FEMPhysics requires it.
   */
  void set_mesh_y_var(unsigned int y_var)
  { _mesh_y_var = y_var; }

  /**
   * Accessor for z-variable of moving mesh System
   */
  unsigned int get_mesh_z_var() const
  { return _mesh_z_var; }

  /**
   * Accessor for z-variable of moving mesh System
   *
   * This should be set automatically if the FEMPhysics requires it.
   */
  void set_mesh_z_var(unsigned int z_var)
  { _mesh_z_var = z_var; }

  /**
   * Accessor for current Elem object
   */
  const Elem& get_elem() const
  { return *elem; }

  /**
   * Accessor for current Elem object
   */
  Elem& get_elem()
  { return *(const_cast<Elem*>(elem)); }

  /**
   * Accessor for current side of Elem object
   */
  unsigned char get_side() const
  { return side; }

  /**
   * Accessor for current edge of Elem object
   */
  unsigned char get_edge() const
  { return edge; }

  /**
   * Accessor for cached element dimension
   */
  unsigned char get_dim() const
  { return dim; }

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
  System *_mesh_sys;

  /**
   * Variables from which to acquire moving mesh information
   */
  unsigned int _mesh_x_var, _mesh_y_var, _mesh_z_var;

  /**
   * Current side for side_* to examine
   */
  unsigned char side;

  /**
   * Current edge for edge_* to examine
   */
  unsigned char edge;

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

  /**
   * Saved pointer to BoundaryInfo on the mesh for this System.  Used
   * to answer boundary id requests.
   */
  BoundaryInfo* _boundary_info;

  /**
   * Current element for element_* to examine
   */
  const Elem *elem;

  /**
   * Cached dimension of elements in this mesh
   */
  unsigned char dim;

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
  libmesh_assert_less ( var, _element_fe_var.size() );
  fe = libmesh_cast_ptr<FEGenericBase<OutputShape>*>( _element_fe_var[var] );
}

inline
FEBase* FEMContext::get_element_fe( unsigned int var ) const
{
  libmesh_assert_less ( var, _element_fe_var.size() );
  return libmesh_cast_ptr<FEBase*>( _element_fe_var[var] );
}

template<typename OutputShape>
inline
void FEMContext::get_side_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const
{
  libmesh_assert_less ( var, _side_fe_var.size() );
  fe = libmesh_cast_ptr<FEGenericBase<OutputShape>*>( _side_fe_var[var] );
}

inline
FEBase* FEMContext::get_side_fe( unsigned int var ) const
{
  libmesh_assert_less ( var, _side_fe_var.size() );
  return libmesh_cast_ptr<FEBase*>( _side_fe_var[var] );
}

template<typename OutputShape>
inline
void FEMContext::get_edge_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const
{
  libmesh_assert_less ( var, _edge_fe_var.size() );
  fe = libmesh_cast_ptr<FEGenericBase<OutputShape>*>( _edge_fe_var[var] );
}

inline
FEBase* FEMContext::get_edge_fe( unsigned int var ) const
{
  libmesh_assert_less ( var, _edge_fe_var.size() );
  return libmesh_cast_ptr<FEBase*>( _edge_fe_var[var] );
}


} // namespace libMesh

#endif // LIBMESH_FEM_CONTEXT_H
