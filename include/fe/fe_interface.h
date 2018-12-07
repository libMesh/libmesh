// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FE_INTERFACE_H
#define LIBMESH_FE_INTERFACE_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"

// C++ includes
#include <map>
#include <vector>

namespace libMesh
{


// forward declarations
class BoundaryInfo;
class DofConstraints;
class DofMap;
class Elem;
class FEType;
class FEComputeData;
class Point;
class MeshBase;
enum FEFamily : int;
enum Order : int;
enum FEFieldType : int;
enum ElemType : int;

#ifdef LIBMESH_ENABLE_PERIODIC
class PeriodicBoundaries;
class PointLocatorBase;
#endif

/**
 * This class provides an encapsulated access to all static
 * public member functions of finite element classes.
 * Using this class, one need not worry about the correct
 * finite element class.
 *
 * \author Daniel Dreyer
 * \date 2002-2007
 * \brief Interface class which provides access to FE functions.
 */
class FEInterface
{
private:

  /**
   * Empty constructor. Do not create an object of this type.
   */
  FEInterface();

public:

  /**
   * Destructor.
   */
  virtual ~FEInterface() {}

  /**
   * \returns The number of shape functions associated with this
   * finite element of type \p fe_t.
   * Automatically decides which finite element class to use.
   *
   * On a p-refined element, \p fe_t.order should be the total order of the element.
   */
  static unsigned int n_shape_functions(const unsigned int dim,
                                        const FEType & fe_t,
                                        const ElemType t);

  /**
   * \returns The number of shape functions associated with this
   * finite element.
   * Automatically decides which finite element class to use.
   *
   * On a p-refined element, \p fe_t.order should be the total order of the element.
   */
  static unsigned int n_dofs(const unsigned int dim,
                             const FEType & fe_t,
                             const ElemType t);

  /**
   * \returns The number of dofs at node n for a finite element
   * of type \p fe_t.
   * Automatically decides which finite element class to use.
   *
   * On a p-refined element, \p fe_t.order should be the total order of the element.
   */
  static unsigned int n_dofs_at_node(const unsigned int dim,
                                     const FEType & fe_t,
                                     const ElemType t,
                                     const unsigned int n);

  typedef unsigned int (*n_dofs_at_node_ptr) (const ElemType,
                                              const Order,
                                              const unsigned int);

  /**
   * \returns A function which evaluates n_dofs_at_node for the
   * requested FE type and dimension.
   */
  static n_dofs_at_node_ptr
  n_dofs_at_node_function(const unsigned int dim,
                          const FEType & fe_t);

  /**
   * \returns The number of dofs interior to the element,
   * not associated with any interior nodes.
   * Automatically decides which finite element class to use.
   *
   * On a p-refined element, \p fe_t.order should be the total order of the element.
   */
  static unsigned int n_dofs_per_elem(const unsigned int dim,
                                      const FEType & fe_t,
                                      const ElemType t);

  /**
   * Fills the vector di with the local degree of freedom indices
   * associated with side \p s of element \p elem
   * Automatically decides which finite element class to use.
   *
   * On a p-refined element, \p fe_t.order should be the base order of the element.
   */
  static void dofs_on_side(const Elem * const elem,
                           const unsigned int dim,
                           const FEType & fe_t,
                           unsigned int s,
                           std::vector<unsigned int> & di);

  /**
   * Fills the vector di with the local degree of freedom indices
   * associated with edge \p e of element \p elem
   * Automatically decides which finite element class to use.
   *
   * On a p-refined element, \p fe_t.order should be the base order of the element.
   */
  static void dofs_on_edge(const Elem * const elem,
                           const unsigned int dim,
                           const FEType & fe_t,
                           unsigned int e,
                           std::vector<unsigned int> & di);

  /**
   * Build the nodal soln from the element soln.
   * This is the solution that will be plotted.
   * Automatically passes the request to the appropriate
   * finite element class member.  To indicate that
   * results from this specific implementation of
   * \p nodal_soln should not be used, the vector
   * \p nodal_soln is returned empty.
   *
   * On a p-refined element, \p fe_t.order should be the base order of the element.
   */
  static void nodal_soln(const unsigned int dim,
                         const FEType & fe_t,
                         const Elem * elem,
                         const std::vector<Number> & elem_soln,
                         std::vector<Number> & nodal_soln);

  /**
   * \returns The point in physical space corresponding to the
   * reference point \p p which is passed in.
   */
  static Point map(unsigned int dim,
                   const FEType & fe_t,
                   const Elem * elem,
                   const Point & p);

  /**
   * \returns The location (on the reference element) of the
   * point \p p located in physical space.  This function requires
   * inverting the (probably nonlinear) transformation map, so
   * it is not trivial. The optional parameter \p tolerance defines
   * how close is "good enough."  The map inversion iteration
   * computes the sequence \f$ \{ p_n \} \f$, and the iteration is
   * terminated when \f$ \|p - p_n\| < \mbox{\texttt{tolerance}} \f$
   */
  static Point inverse_map (const unsigned int dim,
                            const FEType & fe_t,
                            const Elem * elem,
                            const Point & p,
                            const Real tolerance = TOLERANCE,
                            const bool secure = true);

  /**
   * \returns The location (on the reference element) of the points \p
   * physical_points located in physical space.  This function
   * requires inverting the (probably nonlinear) transformation map,
   * so it is not trivial. The location of each point on the reference
   * element is returned in the vector \p reference_points. The
   * optional parameter \p tolerance defines how close is "good
   * enough."  The map inversion iteration computes the sequence \f$
   * \{ p_n \} \f$, and the iteration is terminated when \f$ \|p -
   * p_n\| < \mbox{\texttt{tolerance}} \f$
   */
  static void  inverse_map (const unsigned int dim,
                            const FEType & fe_t,
                            const Elem * elem,
                            const std::vector<Point> & physical_points,
                            std::vector<Point> &       reference_points,
                            const Real tolerance = TOLERANCE,
                            const bool secure = true);

  /**
   * \returns \p true if the point \p p is located on the reference element
   * for element type t, \p false otherwise.
   *
   * Since we are doing floating point comparisons, the parameter
   * \p eps can be specified to indicate a tolerance.  For example,
   * \f$ \xi \le 1 \f$ becomes \f$ \xi \le 1 + \epsilon \f$.
   */
  static bool on_reference_element(const Point & p,
                                   const ElemType t,
                                   const Real eps=TOLERANCE);
  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p. This method allows you to specify the dimension,
   * element type, and order directly. Automatically passes the
   * request to the appropriate finite element class member.
   *
   * \note On a p-refined element, \p fe_t.order should be the total
   * order of the element.
   */
  static Real shape(const unsigned int dim,
                    const FEType & fe_t,
                    const ElemType t,
                    const unsigned int i,
                    const Point & p);

  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p. This method allows you to specify the dimension,
   * element type, and order directly. Automatically passes the
   * request to the appropriate finite element class member.
   *
   * \note On a p-refined element, \p fe_t.order should be the base
   * order of the element.
   */
  static Real shape(const unsigned int dim,
                    const FEType & fe_t,
                    const Elem * elem,
                    const unsigned int i,
                    const Point & p);

  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p. This method allows you to specify the dimension,
   * element type, and order directly. Automatically passes the
   * request to the appropriate *scalar* finite element class member.
   *
   * \note On a p-refined element, \p fe_t.order should be the total
   * order of the element.
   */
  template<typename OutputType>
  static void shape(const unsigned int dim,
                    const FEType & fe_t,
                    const ElemType t,
                    const unsigned int i,
                    const Point & p,
                    OutputType & phi);

  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p. This method allows you to specify the dimension,
   * element type, and order directly. Automatically passes the
   * request to the appropriate *scalar* finite element class member.
   *
   * \note On a p-refined element, \p fe_t.order should be the total
   * order of the element.
   */
  template<typename OutputType>
  static void shape(const unsigned int dim,
                    const FEType & fe_t,
                    const Elem * elem,
                    const unsigned int i,
                    const Point & p,
                    OutputType & phi);

  /**
   * \returns The \f$ j^{th} \f$ coordinate of the gradient of
   * the \f$ i^{th} \f$ shape function at point \p p.
   * This method allows you to specify the dimension,
   * element type, and order directly. Automatically passes the
   * request to the appropriate *scalar* finite element class member.
   *
   * \note On a p-refined element, \p fe_t.order should be the total
   * order of the element.
   */
  static Real shape_deriv(const unsigned int dim,
                          const FEType & fe_t,
                          const ElemType t,
                          const unsigned int i,
                          const unsigned int j,
                          const Point & p);

  /**
   * \returns The \f$ j^{th} \f$ coordinate of the gradient of
   * the \f$ i^{th} \f$ shape function at point \p p.
   * This method allows you to specify the dimension,
   * element type, and order directly. Automatically passes the
   * request to the appropriate *scalar* finite element class member.
   *
   * \note On a p-refined element, \p fe_t.order should be the total
   * order of the element.
   */
  static Real shape_deriv (const unsigned int dim,
                           const FEType & fe_t,
                           const Elem *elem,
                           const unsigned int i,
                           const unsigned int j,
                           const Point & p);


  /**
   * Lets the appropriate child of \p FEBase compute the requested
   * data for the input specified in \p data, and sets the values in
   * \p data.  See this as a generalization of \p shape().  With
   * infinite elements disabled, computes values for all shape
   * functions of \p elem evaluated at \p p.
   *
   * \note On a p-refined element, \p fe_t.order should be the base
   * order of the element.
   */
  static void compute_data(const unsigned int dim,
                           const FEType & fe_t,
                           const Elem * elem,
                           FEComputeData & data);

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Computes the constraint matrix contributions (for
   * non-conforming adapted meshes) corresponding to
   * variable number \p var_number.
   */
  static void compute_constraints (DofConstraints & constraints,
                                   DofMap & dof_map,
                                   const unsigned int variable_number,
                                   const Elem * elem);
#endif // #ifdef LIBMESH_ENABLE_AMR

#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * Computes the constraint matrix contributions (for
   * periodic boundary conditions) corresponding to
   * variable number \p var_number.
   */
  static void compute_periodic_constraints (DofConstraints & constraints,
                                            DofMap & dof_map,
                                            const PeriodicBoundaries & boundaries,
                                            const MeshBase & mesh,
                                            const PointLocatorBase * point_locator,
                                            const unsigned int variable_number,
                                            const Elem * elem);
#endif // #ifdef LIBMESH_ENABLE_PERIODIC

  /**
   * \returns The maximum polynomial degree that the given finite
   * element family can support on the given geometric element.
   */
  static unsigned int max_order (const FEType & fe_t,
                                 const ElemType & el_t);

  /**
   * \returns \p true if separate degrees of freedom must be allocated for
   * vertex DoFs and edge/face DoFs at a hanging node.
   */
  static bool extra_hanging_dofs (const FEType & fe_t);

  /**
   * \returns The number of components of a vector-valued element.
   * Scalar-valued elements return 1.
   */
  static FEFieldType field_type (const FEType & fe_type);

  /**
   * \returns The number of components of a vector-valued element.
   * Scalar-valued elements return 1.
   */
  static FEFieldType field_type (const FEFamily & fe_family);

  /**
   * \returns The number of components of a vector-valued element.
   * Scalar-valued elements return 1.
   */
  static unsigned int n_vec_dim (const MeshBase & mesh,
                                 const FEType & fe_type);

private:


  /**
   * \returns \p true if \p et is an element to be processed by class
   * \p InfFE, false otherwise or when !LIBMESH_ENABLE_INFINITE_ELEMENTS.
   */
  static bool is_InfFE_elem(const ElemType et);


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /*
   * All these private members do the same as their public
   * counterparts, except for infinite elements. This disentangles
   * the calls to \p FE and \p InfFE.
   */

  static unsigned int ifem_n_shape_functions(const unsigned int dim,
                                             const FEType & fe_t,
                                             const ElemType t);

  static unsigned int ifem_n_dofs(const unsigned int dim,
                                  const FEType & fe_t,
                                  const ElemType t);

  static unsigned int ifem_n_dofs_at_node(const unsigned int dim,
                                          const FEType & fe_t,
                                          const ElemType t,
                                          const unsigned int n);

  static unsigned int ifem_n_dofs_per_elem(const unsigned int dim,
                                           const FEType & fe_t,
                                           const ElemType t);

  static void ifem_nodal_soln(const unsigned int dim,
                              const FEType & fe_t,
                              const Elem * elem,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln);

  static Point ifem_map (const unsigned int dim,
                         const FEType & fe_t,
                         const Elem * elem,
                         const Point & p);

  static Point ifem_inverse_map (const unsigned int dim,
                                 const FEType & fe_t,
                                 const Elem * elem,
                                 const Point & p,
                                 const Real tolerance = TOLERANCE,
                                 const bool secure = true);

  static void ifem_inverse_map (const unsigned int dim,
                                const FEType & fe_t,
                                const Elem * elem,
                                const std::vector<Point> & physical_points,
                                std::vector<Point> &       reference_points,
                                const Real tolerance = TOLERANCE,
                                const bool secure = true);


  static bool ifem_on_reference_element(const Point & p,
                                        const ElemType t,
                                        const Real eps);

  static Real ifem_shape(const unsigned int dim,
                         const FEType & fe_t,
                         const ElemType t,
                         const unsigned int i,
                         const Point & p);

  static Real ifem_shape(const unsigned int dim,
                         const FEType & fe_t,
                         const Elem * elem,
                         const unsigned int i,
                         const Point & p);


  static Real ifem_shape_deriv(const unsigned int dim,
                               const FEType & fe_t,
                               const ElemType t,
                               const unsigned int i,
                               const unsigned int j,
                               const Point & p);

  static Real ifem_shape_deriv(const unsigned int dim,
                               const FEType & fe_t,
                               const Elem * elem,
                               const unsigned int i,
                               const unsigned int j,
                               const Point & p);

  static void ifem_compute_data(const unsigned int dim,
                                const FEType & fe_t,
                                const Elem * elem,
                                FEComputeData & data);

#endif


};

} // namespace libMesh

#endif // LIBMESH_FE_INTERFACE_H
