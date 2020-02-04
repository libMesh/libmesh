// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <set>

namespace libMesh
{

// Forward Declarations
template <typename> class BoundaryInfoTempl;
typedef BoundaryInfoTempl<Real> BoundaryInfo;
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
class QBase;
template <typename> class PointTempl;
typedef PointTempl<Real> Point;
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
 * \author Roy H. Stogner
 * \date 2009
 */
class FEMContext : public DiffContext
{
public:

  /**
   * Constructor.  Allocates some but fills no data structures.
   */
  explicit
  FEMContext (const System & sys);

  /**
   * Constructor.  Specify the extra quadrature order instead
   * of getting it from \p sys.
   */
  explicit
  FEMContext (const System & sys, int extra_quadrature_order);

  /**
   * Destructor.
   */
  virtual ~FEMContext ();

  /**
   * Use quadrature rules designed to over-integrate a mass matrix,
   * plus \p extra_quadrature_order.
   */
  void use_default_quadrature_rules(int extra_quadrature_order=0);

  /**
   * Use quadrature rules designed to exactly integrate unweighted
   * undistorted basis functions, plus \p extra_quadrature_order.
   */
  void use_unweighted_quadrature_rules(int extra_quadrature_order=0);

  /**
   * Reports if the boundary id is found on the current side
   */
  bool has_side_boundary_id(boundary_id_type id) const;

  /**
   * Lists the boundary ids found on the current side
   *
   * \deprecated Instead, use the version that takes a reference to a
   * std::set.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  std::vector<boundary_id_type> side_boundary_ids() const;
#endif

  /**
   * As above, but fills in the std::set provided by the user.
   */
  void side_boundary_ids(std::vector<boundary_id_type> & vec_to_fill) const;

  /**
   * \returns The value of the solution variable \p var at the
   * quadrature point \p qp on the current element interior.
   *
   * \note This API is currently present for backward compatibility.
   */
  Number interior_value(unsigned int var, unsigned int qp) const;

  /**
   * \returns The value of the solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This API currently is present for backward compatibility.
   */
  Number side_value(unsigned int var, unsigned int qp) const;

  /**
   * \returns The value of the solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This API is currently present for backward compatibility.
   */
  Number point_value(unsigned int var, const Point & p) const;

  /**
   * \returns The gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This API is currently present for backward compatibility.
   */
  Gradient interior_gradient(unsigned int var, unsigned int qp) const;

  /**
   * \returns The gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This API is currently present for backward compatibility.
   */
  Gradient side_gradient(unsigned int var, unsigned int qp) const;

  /**
   * \returns The gradient of the solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This API is currently present for backward compatibility.
   */
  Gradient point_gradient(unsigned int var, const Point & p) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * \returns The hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This API is currently present for backward compatibility.
   */
  Tensor interior_hessian(unsigned int var, unsigned int qp) const;

  /**
   * \returns The hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This API is currently present for backward compatibility.
   */
  Tensor side_hessian(unsigned int var, unsigned int qp) const;

  /**
   * \returns The hessian of the solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This API currently present for backward compatibility.
   */
  Tensor point_hessian(unsigned int var, const Point & p) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This API is currently present for backward compatibility.
   */
  Number fixed_interior_value(unsigned int var, unsigned int qp) const;

  /**
   * \returns The value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This API is currently present for backward compatibility.
   */
  Number fixed_side_value(unsigned int var, unsigned int qp) const;

  /**
   * \returns The value of the fixed_solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This API is currently present for backward compatibility.
   */
  Number fixed_point_value(unsigned int var, const Point & p) const;

  /**
   * \returns The gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This API is currently present for backward compatibility.
   */
  Gradient fixed_interior_gradient(unsigned int var, unsigned int qp) const;

  /**
   * \returns The gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This API is currently present for backward compatibility.
   */
  Gradient fixed_side_gradient(unsigned int var, unsigned int qp) const;

  /**
   * \returns The gradient of the fixed_solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This API is currently present for backward compatibility.
   */
  Gradient fixed_point_gradient(unsigned int var, const Point & p) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * \returns The hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This API is currently present for backward compatibility.
   */
  Tensor fixed_interior_hessian(unsigned int var, unsigned int qp) const;

  /**
   * \returns The hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This API is currently present for backward compatibility.
   */
  Tensor fixed_side_hessian(unsigned int var, unsigned int qp) const;

  /**
   * \returns The hessian of the fixed_solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This API is currently present for backward compatibility.
   */
  Tensor fixed_point_hessian (unsigned int var, const Point & p) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Accessor for interior finite element object for variable var for
   * the largest dimension in the mesh. We default to the largest mesh dim
   * because this method may be called before the Elem * is set in the FEMContext,
   * e.g. in FEMSystem::init_context (or a subclass).
   * If you have lower dimensional elements in the mesh and need to query for
   * those FE objects, use the alternative get_element_fe method.
   */
  template<typename OutputShape>
  void get_element_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const
  { this->get_element_fe<OutputShape>(var,fe,this->get_dim()); }

  /**
   * Accessor for interior finite element object for scalar-valued variable var
   * for the largest dimension in the mesh. We default to the largest mesh dim
   * because this method may be called before the Elem * is set in the FEMContext,
   * e.g. in FEMSystem::init_context (or a subclass).
   * If you have lower dimensional elements in the mesh and need to query for
   * those FE objects, use the alternative get_element_fe method.
   */
  FEBase * get_element_fe( unsigned int var ) const
  { return this->get_element_fe(var,this->get_dim()); }

  /**
   * Accessor for interior finite element object for variable var for
   * dimension dim.
   */
  template<typename OutputShape>
  void get_element_fe( unsigned int var, FEGenericBase<OutputShape> *& fe,
                       unsigned short dim ) const;

  /**
   * Accessor for interior finite element object for scalar-valued variable var for
   * dimension dim.
   */
  FEBase * get_element_fe( unsigned int var, unsigned short dim ) const;

  /**
   * Accessor for edge/face (2D/3D) finite element object for variable var
   * for the largest dimension in the mesh. We default to the largest mesh dim
   * because this method may be called before the Elem * is set in the FEMContext,
   * e.g. in FEMSystem::init_context (or a subclass).
   * If you have lower dimensional elements in the mesh and need to query for
   * those FE objects, use the alternative get_side_fe method.
   */
  template<typename OutputShape>
  void get_side_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const
  { this->get_side_fe<OutputShape>(var,fe,this->get_dim()); }

  /**
   * Accessor for side finite element object for scalar-valued variable var
   * for the largest dimension in the mesh. We default to the largest mesh dim
   * because this method may be called before the Elem * is set in the FEMContext,
   * e.g. in FEMSystem::init_context (or a subclass).
   * If you have lower dimensional elements in the mesh and need to query for
   * those FE objects, use the alternative get_side_fe method.
   */
  FEBase * get_side_fe( unsigned int var ) const
  { return this->get_side_fe(var,this->get_dim()); }

  /**
   * Accessor for edge/face (2D/3D) finite element object for variable var
   * for dimension dim.
   */
  template<typename OutputShape>
  void get_side_fe( unsigned int var, FEGenericBase<OutputShape> *& fe,
                    unsigned short dim ) const;

  /**
   * Accessor for side finite element object for scalar-valued variable var
   * for dimension dim.
   */
  FEBase * get_side_fe( unsigned int var, unsigned short dim ) const;

  /**
   * Accessor for edge (3D only!) finite element object for variable var.
   */
  template<typename OutputShape>
  void get_edge_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const;

  /**
   * Accessor for edge (3D only!) finite element object for scalar-valued variable var.
   */
  FEBase * get_edge_fe( unsigned int var ) const;

  /**
   * \returns The value of the solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void interior_value(unsigned int var,
                      unsigned int qp,
                      OutputType & u) const;

  /**
   * Fills a vector of values of the _system_vector at the all the quadrature
   * points in the current element interior.
   */
  template<typename OutputType>
  void interior_values(unsigned int var,
                       const NumericVector<Number> & _system_vector,
                       std::vector<OutputType> & interior_values_vector) const;

  /**
   * \returns The value of the solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void side_value(unsigned int var,
                  unsigned int qp,
                  OutputType & u) const;

  /**
   * Fills a vector of values of the _system_vector at the all the quadrature
   * points on the current element side.
   */
  template<typename OutputType>
  void side_values(unsigned int var,
                   const NumericVector<Number> & _system_vector,
                   std::vector<OutputType> & side_values_vector) const;

  /**
   * \returns The value of the solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This is the preferred API.
   *
   * Allows evaluation of points within a relative tolerance outside
   * the element.
   */
  template<typename OutputType>
  void point_value(unsigned int var,
                   const Point & p,
                   OutputType & u,
                   const Real tolerance = TOLERANCE) const;

  /**
   * \returns The gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void interior_gradient(unsigned int var, unsigned int qp,
                         OutputType & du) const;

  /**
   * Fills a vector with the gradient of the solution variable \p var at all the quadrature
   * points in the current element interior.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void interior_gradients(unsigned int var,
                          const NumericVector<Number> & _system_vector,
                          std::vector<OutputType> & interior_gradients_vector) const;

  /**
   * \returns The gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void side_gradient(unsigned int var,
                     unsigned int qp,
                     OutputType & du) const;

  /**
   * Fills a vector with the gradient of the solution variable \p var at all the quadrature
   * points on the current element side.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void side_gradients(unsigned int var,
                      const NumericVector<Number> & _system_vector,
                      std::vector<OutputType> & side_gradients_vector) const;

  /**
   * \returns The gradient of the solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This is the preferred API.
   *
   * Allows evaluation of points within a relative tolerance outside
   * the element.
   */
  template<typename OutputType>
  void point_gradient(unsigned int var,
                      const Point & p,
                      OutputType & grad_u,
                      const Real tolerance = TOLERANCE) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * \returns The hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void interior_hessian(unsigned int var,
                        unsigned int qp,
                        OutputType & d2u) const;

  /**
   * Fills a vector of hessians of the _system_vector at the all the
   * quadrature points in the current element interior. This is the
   * preferred API.
   */
  template<typename OutputType>
  void interior_hessians(unsigned int var,
                         const NumericVector<Number> & _system_vector,
                         std::vector<OutputType> & d2u_vals) const;

  /**
   * \returns The hessian of the solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void side_hessian(unsigned int var,
                    unsigned int qp,
                    OutputType & d2u) const;

  /**
   * Fills a vector of hessians of the _system_vector at the all the
   * quadrature points on the current element side.  This is the
   * preferred API.
   */
  template<typename OutputType>
  void side_hessians(unsigned int var,
                     const NumericVector<Number> & _system_vector,
                     std::vector<OutputType> & d2u_vals) const;

  /**
   * \returns The hessian of the solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This is the preferred API.
   *
   * Allows evaluation of points within a relative tolerance outside
   * the element.
   */
  template<typename OutputType>
  void point_hessian(unsigned int var,
                     const Point & p,
                     OutputType & hess_u,
                     const Real tolerance = TOLERANCE) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The time derivative (rate) of the solution variable
   * \p var at the quadrature point \p qp on the current element
   * interior.
   */
  template<typename OutputType>
  void interior_rate(unsigned int var,
                     unsigned int qp,
                     OutputType & u) const;


  /**
   * \returns The time derivative (rate) of the solution gradient
   * of variable \p var at the quadrature point \p qp on the current
   * element interior.
   */
  template<typename OutputType>
  void interior_rate_gradient(unsigned int var,
                              unsigned int qp,
                              OutputType & u) const;


  /**
   * \returns The time derivative (rate) of the solution variable
   * \p var at the quadrature point \p qp on the current element side.
   */
  template<typename OutputType>
  void side_rate(unsigned int var,
                 unsigned int qp,
                 OutputType & u) const;

  /**
   * \returns The time derivative (rate) of the solution variable
   * \p var at the physical point \p p on the current element.
   */
  template<typename OutputType>
  void point_rate(unsigned int var,
                  const Point & p,
                  OutputType & u) const;

  /**
   * \returns The second time derivative (acceleration) of the solution variable
   * \p var at the quadrature point \p qp on the current element
   * interior.
   */
  template<typename OutputType>
  void interior_accel(unsigned int var,
                      unsigned int qp,
                      OutputType & u) const;

  /**
   * \returns The second time derivative (acceleration) of the solution variable
   * \p var at the quadrature point \p qp on the current element side.
   */
  template<typename OutputType>
  void side_accel(unsigned int var,
                  unsigned int qp,
                  OutputType & u) const;

  /**
   * \returns The second time derivative (acceleration) of the solution variable
   * \p var at the physical point \p p on the current element.
   */
  template<typename OutputType>
  void point_accel(unsigned int var,
                   const Point & p,
                   OutputType & u) const;

  /**
   * \returns The value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void fixed_interior_value(unsigned int var,
                            unsigned int qp,
                            OutputType & u) const;

  /**
   * \returns The value of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void fixed_side_value(unsigned int var,
                        unsigned int qp,
                        OutputType & u) const;

  /**
   * \returns The value of the fixed_solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This is the preferred API.
   *
   * Allows evaluation of points within a relative tolerance outside
   * the element.
   */
  template<typename OutputType>
  void fixed_point_value(unsigned int var,
                         const Point & p,
                         OutputType & u,
                         const Real tolerance = TOLERANCE) const;

  /**
   * \returns The gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void fixed_interior_gradient(unsigned int var,
                               unsigned int qp,
                               OutputType & grad_u) const;

  /**
   * \returns The gradient of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void fixed_side_gradient(unsigned int var,
                           unsigned int qp,
                           OutputType & grad_u) const;

  /**
   * \returns The gradient of the fixed_solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This is the preferred API.
   *
   * Allows evaluation of points within a relative tolerance outside
   * the element.
   */
  template<typename OutputType>
  void fixed_point_gradient(unsigned int var,
                            const Point & p,
                            OutputType & grad_u,
                            const Real tolerance = TOLERANCE) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * \returns The hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element interior.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void fixed_interior_hessian(unsigned int var,
                              unsigned int qp,
                              OutputType & hess_u) const;

  /**
   * \returns The hessian of the fixed_solution variable \p var at the quadrature
   * point \p qp on the current element side.
   *
   * \note This is the preferred API.
   */
  template<typename OutputType>
  void fixed_side_hessian(unsigned int var,
                          unsigned int qp,
                          OutputType & hess_u) const;

  /**
   * \returns The hessian of the fixed_solution variable \p var at the physical
   * point \p p on the current element.
   *
   * \note This is the preferred API.
   *
   * Allows evaluation of points within a relative tolerance outside
   * the element.
   */
  template<typename OutputType>
  void fixed_point_hessian(unsigned int var,
                           const Point & p,
                           OutputType & hess_u,
                           const Real tolerance = TOLERANCE) const;

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * \returns The curl of the solution variable \p var at the physical
   * point \p p on the current element.
   */
  template<typename OutputType>
  void interior_curl(unsigned int var,
                     unsigned int qp,
                     OutputType & curl_u) const;

  /**
   * \returns The curl of the solution variable \p var at the physical
   * point \p p on the current element.
   *
   * Allows evaluation of points within a relative tolerance outside
   * the element.
   */
  template<typename OutputType>
  void point_curl(unsigned int var,
                  const Point & p,
                  OutputType & curl_u,
                  const Real tolerance = TOLERANCE) const;

  /**
   * \returns The divergence of the solution variable \p var at the physical
   * point \p p on the current element.
   */
  template<typename OutputType>
  void interior_div(unsigned int var,
                    unsigned int qp,
                    OutputType & div_u) const;

  // should be protected:
  /**
   * Resets the current time in the context. Additionally, reinitialize Elem
   * and FE objects if there's a moving mesh present in the system such that
   * the mesh is deformed to its position at \f$ t_{\theta} \f$.
   */
  virtual void elem_reinit(Real theta) override;

  /**
   * Resets the current time in the context. Additionally, reinitialize Elem
   * and FE objects if there's a moving mesh present in the system such that
   * the mesh is deformed to its position at \f$ t_{\theta} \f$.
   */
  virtual void elem_side_reinit(Real theta) override;

  /**
   * Resets the current time in the context. Additionally, reinitialize Elem
   * and FE objects if there's a moving mesh present in the system such that
   * the mesh is deformed to its position at \f$ t_{\theta} \f$.
   */
  virtual void elem_edge_reinit(Real theta) override;

  /**
   * Gives derived classes the opportunity to reinitialize data needed
   * for nonlocal calculations at a new point within a timestep
   */
  virtual void nonlocal_reinit(Real theta) override;

  /**
   * Reinitializes local data vectors/matrices on the current geometric element
   */
  virtual void pre_fe_reinit(const System &,
                             const Elem * e);

  /**
   * Reinitializes interior FE objects on the current geometric element
   */
  virtual void elem_fe_reinit(const std::vector<Point> * const pts = nullptr);

  /**
   * Reinitializes side FE objects on the current geometric element
   */
  virtual void side_fe_reinit();

  /**
   * Reinitializes edge FE objects on the current geometric element
   */
  virtual void edge_fe_reinit();

  /**
   * Accessor for element interior quadrature rule for the dimension of the
   * current _elem.
   */
  const QBase & get_element_qrule() const
  { return this->get_element_qrule(this->get_elem_dim()); }

  /**
   * Accessor for element side quadrature rule for the dimension of the
   * current _elem.
   */
  const QBase & get_side_qrule() const
  { return this->get_side_qrule(this->get_elem_dim()); }

  /**
   * Accessor for element interior quadrature rule.
   */
  const QBase & get_element_qrule( unsigned short dim ) const
  { libmesh_assert(_element_qrule[dim]);
    return *(this->_element_qrule[dim]); }

  /**
   * Accessor for element side quadrature rule.
   */
  const QBase & get_side_qrule( unsigned short dim ) const
  {
    libmesh_assert(_side_qrule[dim]);
    return *(this->_side_qrule[dim]);
  }

  /**
   * Accessor for element edge quadrature rule.
   */
  const QBase & get_edge_qrule() const
  { return *(this->_edge_qrule); }

  /**
   * Tells the FEMContext that system \p sys contains the
   * isoparametric Lagrangian variables which correspond to the
   * coordinates of mesh nodes, in problems where the mesh itself is
   * expected to move in time.
   *
   * This should be set automatically if the FEMPhysics requires it.
   */
  virtual void set_mesh_system(System * sys)
  { this->_mesh_sys = sys; }

  /**
   * Accessor for moving mesh System
   */
  const System * get_mesh_system() const
  { return this->_mesh_sys; }

  /**
   * Accessor for moving mesh System
   */
  System * get_mesh_system()
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
   * Test for current Elem object
   */
  bool has_elem() const
  { return (this->_elem != nullptr); }

  /**
   * Accessor for current Elem object
   */
  const Elem & get_elem() const
  { libmesh_assert(this->_elem);
    return *(this->_elem); }

  /**
   * Accessor for current Elem object
   */
  Elem & get_elem()
  { libmesh_assert(this->_elem);
    return *(const_cast<Elem *>(this->_elem)); }

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
   * Accessor for cached mesh dimension. This is the largest dimension
   * of the elements in the mesh. For the dimension of this->_elem, use
   * get_elem_dim();
   */
  unsigned char get_dim() const
  { return this->_dim; }

  /**
   * \returns The dimension of this->_elem. For mixed dimension meshes, this
   * may be different from get_dim().
   */
  unsigned char get_elem_dim() const
  { return _elem_dim; }

  /**
   * \returns Set of dimensions of elements present in the mesh at
   * context initialization.
   */
  const std::set<unsigned char> & elem_dimensions() const
  { return _elem_dims; }

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
   * Enum describing what data to use when initializing algebraic
   * structures on each element.
   */
  enum AlgebraicType { NONE = 0,  // Do not reinitialize dof_indices
                       DOFS_ONLY, // Reinitialize dof_indices, not
                                  // algebraic structures
                       CURRENT,   // Use dof_indices, current solution
                       OLD,       // Use old_dof_indices, custom solution
                       OLD_DOFS_ONLY}; // Reinitialize old_dof_indices, not
                                       // algebraic structures

  /**
   * Setting which determines whether to initialize algebraic
   * structures (elem_*) on each element and set their values from
   * current_local_solution.  Algebraic initialization may be disabled
   * for efficiency in cases where FEMContext is only used as a
   * convenient container of FE objects.
   */
  void set_algebraic_type(const AlgebraicType atype)
  { _atype = atype; }

  /*
   * Get the current AlgebraicType setting
   */
  AlgebraicType algebraic_type() const { return _atype; }

  /**
   * Set a NumericVector to be used in place of current_local_solution
   * for calculating elem_solution.  Set to nullptr to restore the
   * current_local_solution behavior.  Advanced DifferentiableSystem
   * specific capabilities will only be enabled in the
   * current_local_solution case.
   */
  void set_custom_solution(const NumericVector<Number> * custom_sol)
  { _custom_solution = custom_sol; }

  /**
   * Calls set_jacobian_tolerance() on all the FE objects controlled
   * by this class. (Actually, it calls this on the underlying)
   */
  void set_jacobian_tolerance(Real tol);

  /**
   * System from which to acquire moving mesh information
   */
  System * _mesh_sys;

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

  /**
   * Helper function for creating quadrature rules
   */
  FEType find_hardest_fe_type();

  /**
   * Helper function for attaching quadrature rules
   */
  void attach_quadrature_rules();

  /**
   * Helper function to reduce some code duplication in the *_point_* methods.
   */
  template<typename OutputShape>
  FEGenericBase<OutputShape> * build_new_fe(const FEGenericBase<OutputShape> * fe,
                                            const Point & p,
                                            const Real tolerance = TOLERANCE) const;

protected:

  /**
   * Keep track of what type of algebra reinitialization is to be done
   */
  AlgebraicType _atype;

  /**
   * Data with which to do algebra reinitialization
   */
  const NumericVector<Number> * _custom_solution;

  mutable std::unique_ptr<FEGenericBase<Real>>         _real_fe;
  mutable std::unique_ptr<FEGenericBase<RealGradient>> _real_grad_fe;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  mutable bool _real_fe_is_inf;
  mutable bool _real_grad_fe_is_inf;
#endif

  template<typename OutputShape>
  FEGenericBase<OutputShape> * cached_fe( const unsigned int elem_dim,
                                          const FEType fe_type ) const;

  /**
   * Helper function to promote accessor usage
   */
  void set_elem( const Elem * e );

  // gcc-3.4, oracle 12.3 require this typedef to be public
  // in order to use it in a return type
public:

  /**
   * Helper typedef to simplify refactoring
   */
  typedef const DenseSubVector<Number> & (DiffContext::*diff_subsolution_getter)(unsigned int) const;

protected:
  /**
   * Helper nested class for C++03-compatible "template typedef"
   */
  template <typename OutputType>
  struct FENeeded
  {
    // Rank decrementer helper types
    typedef typename TensorTools::DecrementRank<OutputType>::type Rank1Decrement;
    typedef typename TensorTools::DecrementRank<Rank1Decrement>::type Rank2Decrement;

    // Typedefs for "Value getter" function pointer
    typedef typename TensorTools::MakeReal<OutputType>::type value_shape;
    typedef FEGenericBase<value_shape> value_base;
    typedef void (FEMContext::*value_getter) (unsigned int, value_base *&, unsigned short) const;

    // Typedefs for "Grad getter" function pointer
    typedef typename TensorTools::MakeReal<Rank1Decrement>::type grad_shape;
    typedef FEGenericBase<grad_shape> grad_base;
    typedef void (FEMContext::*grad_getter) (unsigned int, grad_base *&, unsigned short) const;

    // Typedefs for "Hessian getter" function pointer
    typedef typename TensorTools::MakeReal<Rank2Decrement>::type hess_shape;
    typedef FEGenericBase<hess_shape> hess_base;
    typedef void (FEMContext::*hess_getter) (unsigned int, hess_base *&, unsigned short) const;
  };



  /**
   * Helper function to reduce some code duplication in the
   * *interior_value methods.
   */
  template<typename OutputType,
           typename FENeeded<OutputType>::value_getter fe_getter,
           diff_subsolution_getter subsolution_getter>
  void some_value(unsigned int var, unsigned int qp, OutputType & u) const;

  /**
   * Helper function to reduce some code duplication in the
   * *interior_gradient methods.
   */
  template<typename OutputType,
           typename FENeeded<OutputType>::grad_getter fe_getter,
           diff_subsolution_getter subsolution_getter>
  void some_gradient(unsigned int var, unsigned int qp, OutputType & u) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  /**
   * Helper function to reduce some code duplication in the
   * *interior_hessian methods.
   */
  template<typename OutputType,
           typename FENeeded<OutputType>::hess_getter fe_getter,
           diff_subsolution_getter subsolution_getter>
  void some_hessian(unsigned int var, unsigned int qp, OutputType & u) const;
#endif

  /**
   * Finite element objects for each variable's interior, sides and edges.
   * We store FE objects for each element dimension present in the mesh,
   * except for edge_fe which only applies to 3D elements.
   */
  std::vector<std::map<FEType, std::unique_ptr<FEAbstract<>>>> _element_fe;
  std::vector<std::map<FEType, std::unique_ptr<FEAbstract<>>>> _side_fe;
  std::map<FEType, std::unique_ptr<FEAbstract<>>> _edge_fe;


  /**
   * Pointers to the same finite element objects, but indexed
   * by variable number. We store FE objects for each element dimension
   * present in the mesh, except for edge_fe_var which only applies
   * for 3D elements.
   */
  std::vector<std::vector<FEAbstract<> *>> _element_fe_var;
  std::vector<std::vector<FEAbstract<> *>> _side_fe_var;
  std::vector<FEAbstract<> *> _edge_fe_var;

  /**
   * Saved reference to BoundaryInfo on the mesh for this System.
   * Used to answer boundary id requests.
   */
  const BoundaryInfo & _boundary_info;

  /**
   * Current element for element_* to examine
   */
  const Elem * _elem;

  /**
   * Cached dimension of largest dimension element in this mesh
   */
  unsigned char _dim;

  /**
   * Cached dimension of this->_elem.
   */
  unsigned char _elem_dim;

  /**
   * Cached dimensions of elements in the mesh, plus dimension 0 if
   * SCALAR variables are in use.
   */
  std::set<unsigned char> _elem_dims;

  /**
   * Quadrature rule for element interior.
   * The FEM context will try to find a quadrature rule that
   * correctly integrates all variables. We prepare quadrature
   * rules for each element dimension in the mesh.
   */
  std::vector<std::unique_ptr<QBase>> _element_qrule;

  /**
   * Quadrature rules for element sides
   * The FEM context will try to find a quadrature rule that
   * correctly integrates all variables. We prepare quadrature
   * rules for each element dimension in the mesh.
   */
  std::vector<std::unique_ptr<QBase>> _side_qrule;

  /**
   * Quadrature rules for element edges.  If the FEM context is told
   * to prepare for edge integration on 3D elements, it will try to
   * find a quadrature rule that correctly integrates all variables.
   * Because edge rules only apply to 3D elements, we don't need to
   * worry about multiple dimensions
   */
  std::unique_ptr<QBase> _edge_qrule;

  /**
   * The extra quadrature order for this context.
   */
  int _extra_quadrature_order;

private:
  /**
   * Helper function used in constructors to set up internal data.
   */
  void init_internal_data(const System & sys);

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
void FEMContext::get_element_fe( unsigned int var, FEGenericBase<OutputShape> *& fe,
                                 unsigned short dim ) const
{
  libmesh_assert( !_element_fe_var[dim].empty() );
  libmesh_assert_less ( var, (_element_fe_var[dim].size() ) );
  fe = cast_ptr<FEGenericBase<OutputShape> *>( (_element_fe_var[dim][var] ) );
}

inline
FEBase * FEMContext::get_element_fe( unsigned int var, unsigned short dim ) const
{
  libmesh_assert( !_element_fe_var[dim].empty() );
  libmesh_assert_less ( var, (_element_fe_var[dim].size() ) );
  return cast_ptr<FEBase *>( (_element_fe_var[dim][var] ) );
}

template<typename OutputShape>
inline
void FEMContext::get_side_fe( unsigned int var, FEGenericBase<OutputShape> *& fe,
                              unsigned short dim ) const
{
  libmesh_assert( !_side_fe_var[dim].empty() );
  libmesh_assert_less ( var, (_side_fe_var[dim].size() ) );
  fe = cast_ptr<FEGenericBase<OutputShape> *>( (_side_fe_var[dim][var] ) );
}

inline
FEBase * FEMContext::get_side_fe( unsigned int var, unsigned short dim ) const
{
  libmesh_assert( !_side_fe_var[dim].empty() );
  libmesh_assert_less ( var, (_side_fe_var[dim].size() ) );
  return cast_ptr<FEBase *>( (_side_fe_var[dim][var] ) );
}

template<typename OutputShape>
inline
void FEMContext::get_edge_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const
{
  libmesh_assert_less ( var, _edge_fe_var.size() );
  fe = cast_ptr<FEGenericBase<OutputShape> *>( _edge_fe_var[var] );
}

inline
FEBase * FEMContext::get_edge_fe( unsigned int var ) const
{
  libmesh_assert_less ( var, _edge_fe_var.size() );
  return cast_ptr<FEBase *>( _edge_fe_var[var] );
}


} // namespace libMesh

#endif // LIBMESH_FEM_CONTEXT_H
