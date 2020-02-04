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



#ifndef LIBMESH_INF_FE_H
#define LIBMESH_INF_FE_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/fe_base.h"
#include "libmesh/inf_fe_map.h"

// C++ includes
#include <cstddef>

namespace libMesh
{


// forward declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
class FEComputeData;


/**
 * Infinite elements are in some sense directional, compared
 * to conventional finite elements.  All methods related
 * to the radial part, which extends perpendicular from the base,
 * are collected in this class.  This class offers static methods for
 * use in \p InfFE and \p InfFEMap
 *
 * \author Daniel Dreyer
 * \date 2003
 */
class InfFERadial
{
private:

  /**
   * Never use an object of this type.
   */
  InfFERadial() {}

public:

  /**
   * \returns The decay in the radial direction of
   * the \p Dim dimensional infinite element.
   */
  static Real decay (const unsigned int dim, const Real v);

  /**
   * \returns The first (local) derivative of the
   * decay in radial direction of the infinite element.
   *
   * This is only valid for 3D?? - RHS
   */
  static Real decay_deriv (const Real) { return -.5; }

  /**
   * \returns The radial weight D, used as an additional weight
   * for the test function, evaluated at local radial coordinate \p v.
   */
  static Real D (const Real v) { return (1.-v)*(1.-v)/4.; }

  /**
   * \returns The first (local) radial derivative of the radial weight D.
   */
  static Real D_deriv (const Real v) { return (v-1.)/2.; }

  /**
   * \returns The Order of the mapping functions
   * in the radial direction. Currently, this is always \p FIRST.
   */
  static Order mapping_order() { return FIRST; }

  /**
   * \returns The number of shape functions in radial direction
   * associated with this infinite element.
   * Either way, if the modes are stored as nodal dofs (\p n_dofs_at_node)
   * or as element dofs (\p n_dofs_per_elem), in each case we have the
   * same number of modes in radial direction.
   *
   * \note For the case of 1D infinite elements, in the base the
   * dof-per-node scheme is used.
   *
   * From the formulation of the infinite elements, we have
   * 1 mode, when \p o_radial=CONST.
   * Therefore, we have a total of \p o_radial+1 modes in radial direction.
   */
  static unsigned int n_dofs (const Order o_radial)
  { return static_cast<unsigned int>(o_radial)+1; }

  /**
   * \returns The number of dofs in radial direction on "onion slice"
   * \p n (either 0 or 1) for an infinite element of type \p inf_elem_type and
   * radial order \p o_radial.
   *
   * Currently, the first radial mode is associated with the nodes in
   * the base.  All higher radial modes are associated with
   * the physically existing nodes further out.
   */
  static unsigned int n_dofs_at_node (const Order o_radial,
                                      const unsigned int n_onion);

  /**
   * \returns The number of modes in radial direction interior to the element,
   * not associated with any interior nodes.
   *
   * \note These modes are a discontinuous approximation, therefore
   * we have no special formulation for coupling in the base, like in the
   * case of associating (possibly) multiple dofs per (outer) node.
   */
  static unsigned int n_dofs_per_elem (const Order o_radial)
  { return static_cast<unsigned int>(o_radial)+1; }

};



/**
 * This nested class contains most of the static methods related
 * to the base part of an infinite element.  Only static members
 * are provided, for use within \p InfFE and \p InfFEMap.
 *
 * \author Daniel Dreyer
 * \date 2003
 */
class InfFEBase
{
private:

  /**
   * Never use an object of this type.
   */
  InfFEBase() {}

public:

  /**
   * Build the base element of an infinite element.  Be careful,
   * this method allocates memory!  So be sure to delete the
   * new element afterward.
   */
  static Elem * build_elem (const Elem * inf_elem);

  /**
   * \returns The base element associated to
   * \p type.  This is, for example, \p TRI3 for
   * \p INFPRISM6.
   */
  static ElemType get_elem_type (const ElemType type);

  /**
   * \returns The number of shape functions used in the
   * mapping in the base element of type \p base_elem_type
   * mapped with order \p base_mapping_order
   */
  static unsigned int n_base_mapping_sf (const Elem & base_elem,
                                         const Order base_mapping_order);

};



/**
 * A specific instantiation of the \p FEBase class. This
 * class is templated, and specific template instantiations
 * will result in different Infinite Element families, similar
 * to the \p FE class.  \p InfFE builds a \p FE<Dim-1,T_base>,
 * and most of the requests related to the base are handed over
 * to this object.  All methods related to the radial part
 * are collected in the class \p InfFERadial.  Similarly,
 * most of the static methods concerning base approximation
 * are contained in \p InfFEBase.
 *
 * Having different shape approximation families in radial direction
 * introduces the requirement for an additional \p Order in this
 * class. Therefore, the \p FEType internals change when infinite
 * elements are enabled.
 * When the specific infinite element type is not known at compile
 * time, use the \p FEBase::build() member to create abstract
 * (but still optimized) infinite elements at run time.
 *
 * The node numbering scheme is the one from the current
 * infinite element.  Each node in the base holds exactly
 * the same number of dofs as an adjacent conventional \p FE
 * would contain.  The nodes further out hold the additional
 * dof necessary for radial approximation.  The order of the outer nodes'
 * components is such that the radial shapes have highest
 * priority, followed by the base shapes.
 *
 * \author Daniel Dreyer
 * \date 2003
 * \brief Base class for all the infinite geometric element types.
 */
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
class InfFE : public FEBase
{

  /*
   * Protect the nested class
   */
protected:


public:

  // InfFE continued

  /**
   * Constructor and empty destructor.
   * Initializes some data structures.  Builds a \p FE<Dim-1,T_base>
   * object to handle  approximation in the base, so that
   * there is no need to template \p InfFE<Dim,T_radial,T_map> also with
   * respect to the base approximation \p T_base.
   *
   * The same remarks concerning compile-time optimization for
   * \p FE also hold for \p InfFE.  Use the
   * \p FEBase::build_InfFE(const unsigned int, const FEType &)
   * method to build specific instantiations of \p InfFE at
   * run time.
   */
  explicit
  InfFE(const FEType & fet);
  ~InfFE() {}

  // The static public members for access from FEInterface etc

  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method lets you specify the relevant
   * data directly, and is therefore allowed to be static.
   *
   * \note This class member is not as efficient as its counterpart in
   * \p FE<Dim,T>, and is not employed in the \p reinit() cycle.
   *
   * \note This method does not return physically correct shapes,
   * instead use \p compute_data().  The \p shape() methods should
   * only be used for mapping.
   */
  static Real shape(const FEType & fet,
                    const ElemType t,
                    const unsigned int i,
                    const Point & p);

  /**
   * \returns The value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method lets you specify the relevant
   * data directly, and is therefore allowed to be static.
   *
   * \note This class member is not as efficient as its counterpart in
   * \p FE<Dim,T>, and is not employed in the \p reinit() cycle.
   *
   * \note This method does not return physically correct shapes,
   * instead use \p compute_data().  The \p shape() methods should
   * only be used for mapping.
   */
  static Real shape(const FEType & fet,
                    const Elem * elem,
                    const unsigned int i,
                    const Point & p);

  /**
   * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at point \p p. This method lets you specify the relevant
   * data directly, and is therefore allowed to be static.
   *
   * \note This class member is not as efficient as its counterpart in
   * \p FE<Dim,T>, and is not employed in the \p reinit() cycle.
   *
   * \note This method does not return physically correct shape gradients,
   * instead use \p compute_data().  The \p shape_deriv() methods should
   * only be used for mapping.
   */
  static Real shape_deriv (const FEType & fet,
                           const Elem * inf_elem,
                           const unsigned int i,
                           const unsigned int j,
                           const Point & p);

  /**
   * \returns The \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at point \p p. This method lets you specify the relevant
   * data directly, and is therefore allowed to be static.
   *
   * \note This class member is not as efficient as its counterpart in
   * \p FE<Dim,T>, and is not employed in the \p reinit() cycle.
   *
   * \note This method does not return physically correct shape gradients,
   * instead use \p compute_data().  The \p shape_deriv() methods should
   * only be used for mapping.
   */
  static Real shape_deriv (const FEType & fet,
                           const ElemType inf_elem_type,
                           const unsigned int i,
                           const unsigned int j,
                           const Point & p);

  /**
   * Generalized version of \p shape(), takes an \p Elem *.  The \p data
   * contains both input and output parameters.  For frequency domain
   * simulations, the complex-valued shape is returned.  In time domain
   * both the computed shape, and the phase is returned.
   *
   * \note The phase (proportional to the distance of the \p Point \p
   * data.p from the envelope) is actually a measure how far into the
   * future the results are.
   */
  static void compute_data(const FEType & fe_t,
                           const Elem * inf_elem,
                           FEComputeData & data);

  /**
   * \returns The number of shape functions associated with
   * a finite element of type \p t and approximation order \p o.
   */
  static unsigned int n_shape_functions (const FEType & fet,
                                         const ElemType t)
  { return n_dofs(fet, t); }

  /**
   * \returns The number of shape functions associated with this
   * infinite element.  Currently, we have \p o_radial+1 modes in
   * radial direction, and \code FE<Dim-1,T>::n_dofs(...) \endcode
   * in the base.
   */
  static unsigned int n_dofs(const FEType & fet,
                             const ElemType inf_elem_type);

  /**
   * \returns The number of dofs at infinite element node \p n
   * (not dof!) for an element of type \p t and order \p o.
   */
  static unsigned int n_dofs_at_node(const FEType & fet,
                                     const ElemType inf_elem_type,
                                     const unsigned int n);

  /**
   * \returns The number of dofs interior to the element,
   * not associated with any interior nodes.
   */
  static unsigned int n_dofs_per_elem(const FEType & fet,
                                      const ElemType inf_elem_type);

  /**
   * \returns The continuity of the element.
   */
  virtual FEContinuity get_continuity() const override
  { return C_ZERO; }  // FIXME - is this true??

  /**
   * \returns \p true if the element's higher order shape functions are
   * hierarchic
   */
  virtual bool is_hierarchic() const override
  { return false; }  // FIXME - Inf FEs don't handle p elevation yet

  /**
   * Usually, this method would build the nodal soln from the
   * element soln.  But infinite elements require additional
   * simulation-specific data to compute physically correct
   * results.  Use \p compute_data() to compute results.  For
   * compatibility an empty vector is returned.
   */
  static void nodal_soln(const FEType & fet,
                         const Elem * elem,
                         const std::vector<Number> & elem_soln,
                         std::vector<Number> & nodal_soln);


  static Point map (const Elem * inf_elem,
                    const Point & reference_point) {
    // libmesh_deprecated(); // soon
    return InfFEMap::map(Dim, inf_elem, reference_point);
  }


  static Point inverse_map (const Elem * elem,
                            const Point & p,
                            const Real tolerance = TOLERANCE,
                            const bool secure = true) {
    // libmesh_deprecated(); // soon
    return InfFEMap::inverse_map(Dim, elem, p, tolerance, secure);
  }


  static void inverse_map (const Elem * elem,
                           const std::vector<Point> & physical_points,
                           std::vector<Point> &       reference_points,
                           const Real tolerance = TOLERANCE,
                           const bool secure = true) {
    // libmesh_deprecated(); // soon
    return InfFEMap::inverse_map(Dim, elem, physical_points,
                                 reference_points, tolerance, secure);
  }


  // The workhorses of InfFE. These are often used during matrix assembly.

  /**
   * This is at the core of this class. Use this for each
   * new element in the mesh.  Reinitializes all the physical
   * element-dependent data based on the current element
   * \p elem.
   * \note: pts need to be in reference space coordinates, not physical ones.
   */
  virtual void reinit (const Elem * elem,
                       const std::vector<Point> * const pts = nullptr,
                       const std::vector<Real> * const weights = nullptr) override;

  /**
   * Reinitializes all the physical
   * element-dependent data based on the \p side of an infinite
   * element.
   */
  virtual void reinit (const Elem * inf_elem,
                       const unsigned int s,
                       const Real tolerance = TOLERANCE,
                       const std::vector<Point> * const pts = nullptr,
                       const std::vector<Real> * const weights = nullptr) override;

  /**
   * Not implemented yet.  Reinitializes all the physical
   * element-dependent data based on the \p edge of an infinite
   * element.
   */
  virtual void edge_reinit (const Elem * elem,
                            const unsigned int edge,
                            const Real tolerance = TOLERANCE,
                            const std::vector<Point> * const pts = nullptr,
                            const std::vector<Real> * const weights = nullptr) override;

  /**
   * Computes the reference space quadrature points on the side of
   * an element based on the side quadrature points.
   */
  virtual void side_map (const Elem * /* elem */,
                         const Elem * /* side */,
                         const unsigned int /* s */,
                         const std::vector<Point> & /* reference_side_points */,
                         std::vector<Point> & /* reference_points */) override
  {
    libmesh_not_implemented();
  }

  /**
   * The use of quadrature rules with the \p InfFE class is somewhat
   * different from the approach of the \p FE class.  While the
   * \p FE class requires an appropriately initialized quadrature
   * rule object, and simply uses it, the \p InfFE class requires only
   * the quadrature rule object of the current \p FE class.
   * From this \p QBase *, it determines the necessary data,
   * and builds two appropriate quadrature classes, one for radial,
   * and another for base integration, using the convenient
   * \p QBase::build() method.
   */
  virtual void attach_quadrature_rule (QBase * q) override;

  /**
   * \returns The number of shape functions associated with
   * this infinite element.
   */
  virtual unsigned int n_shape_functions () const override
  { return _n_total_approx_sf; }

  /**
   * \returns The total number of quadrature points.  Call this
   * to get an upper bound for the \p for loop in your simulation
   * for matrix assembly of the current element.
   */
  virtual unsigned int n_quadrature_points () const override
  { libmesh_assert(radial_qrule); return _n_total_qp; }


protected:

  // static members used by the workhorses

  /**
   * \returns The value of the \f$ i^{th} \f$ polynomial evaluated
   * at \p v.  This method provides the approximation
   * in radial direction for the overall shape functions,
   * which is defined in \p InfFE::shape().
   * This method is allowed to be static, since it is independent
   * of dimension and base_family.  It is templated, though,
   * w.r.t. to radial \p FEFamily.
   *
   * \returns The value of the \f$ i^{th} \f$ mapping shape function
   * in radial direction evaluated at \p v when T_radial ==
   * INFINITE_MAP.  Currently, only one specific mapping shape is
   * used.  Namely the one by Marques JMMC, Owen DRJ: Infinite
   * elements in quasi-static materially nonlinear problems, Computers
   * and Structures, 1984.
   */
  static Real eval(Real v,
                   Order o_radial,
                   unsigned int i);

  /**
   * \returns The value of the first derivative of the
   * \f$ i^{th} \f$ polynomial at coordinate \p v.
   * See \p eval for details.
   */
  static Real eval_deriv(Real v,
                         Order o_radial,
                         unsigned int i);



  // Non-static members used by the workhorses

  /**
   * Updates the protected member \p base_elem to the appropriate base element
   * for the given \p inf_elem.
   */
  void update_base_elem (const Elem * inf_elem);

  /**
   * Do not use this derived member in \p InfFE<Dim,T_radial,T_map>.
   */
  virtual void init_base_shape_functions(const std::vector<Point> &,
                                         const Elem *) override
  { libmesh_not_implemented(); }

  /**
   * Some of the member data only depend on the radial part of the
   * infinite element.  The parts that only change when the radial
   * order changes, are initialized here.
   */
  void init_radial_shape_functions(const Elem * inf_elem,
                                   const std::vector<Point> * radial_pts = nullptr);

  /**
   * Initialize all the data fields like \p weight, \p mode,
   * \p phi, \p dphidxi, \p dphideta, \p dphidzeta, etc.
   * for the current element.  This method prepares the data
   * related to the base part, and some of the combined fields.
   */
  void init_shape_functions(const std::vector<Point> & radial_qp,
                            const std::vector<Point> & base_qp,
                            const Elem * inf_elem);

  /**
   * Initialize all the data fields like \p weight,
   * \p phi, etc for the side \p s.
   */
  void init_face_shape_functions (const std::vector<Point> &,
                                  const Elem * inf_side);

  /**
   * Combines the shape functions, which were formed in
   * \p init_shape_functions(Elem *), with geometric data.
   * Has to be called every time the geometric configuration
   * changes.  Afterward, the fields are ready to be used
   * to compute global derivatives, the jacobian etc, see
   * \p FEAbstract::compute_map().
   */
  void combine_base_radial(const Elem * inf_elem);

  /**
   * After having updated the jacobian and the transformation
   * from local to global coordinates in FEAbstract::compute_map(),
   * the first derivatives of the shape functions are
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx/y/z, \p dphasedx/y/z, \p dweight. This method
   * should barely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected.
   */
  virtual void compute_shape_functions(const Elem *, const std::vector<Point> &) override;



  // Miscellaneous static members

  /**
   * Computes the indices in the base \p base_node and in radial
   * direction \p radial_node (either 0 or 1) associated to the
   * node \p outer_node_index of an infinite element of type
   * \p inf_elem_type.
   */
  static void compute_node_indices (const ElemType inf_elem_type,
                                    const unsigned int outer_node_index,
                                    unsigned int & base_node,
                                    unsigned int & radial_node);

  /**
   * Does the same as \p compute_node_indices(), but stores
   * the maps for the current element type.  Provided the
   * infinite element type changes seldom, this is probably
   * faster than using \p compute_node_indices () alone.
   * This is possible since the number of nodes is not likely
   * to change.
   */
  static void compute_node_indices_fast (const ElemType inf_elem_type,
                                         const unsigned int outer_node_index,
                                         unsigned int & base_node,
                                         unsigned int & radial_node);

  /**
   * Computes the indices of shape functions in the base \p base_shape and
   * in radial direction \p radial_shape (0 in the base, \f$ \ge 1 \f$ further
   * out) associated to the shape with global index \p i of an infinite element
   * of type \p inf_elem_type.
   */
  static void compute_shape_indices (const FEType & fet,
                                     const ElemType inf_elem_type,
                                     const unsigned int i,
                                     unsigned int & base_shape,
                                     unsigned int & radial_shape);

  /**
   * the radial distance of the base nodes from the origin
   */
  std::vector<Real> dist;

  /**
   * the additional radial weight \f$ 1/{r^2} \f$ in local coordinates,
   * over all quadrature points. The weight does not vary in base
   * direction.  However, for uniform access to the data fields from the
   * outside, this data field is expanded to all quadrature points.
   */
  std::vector<Real> dweightdv;

  /**
   * the radial decay \f$ 1/r \f$ in local coordinates.  Needed when
   * setting up the overall shape functions.
   *
   * \note It is this decay which ensures that the Sommerfeld
   * radiation condition is satisfied in advance.
   */
  std::vector<Real> som;
  /**
   * the first local derivative of the radial decay \f$ 1/r \f$ in local
   * coordinates.  Needed when setting up the overall shape functions.
   */
  std::vector<Real> dsomdv;

  /**
   * the radial approximation shapes in local coordinates
   * Needed when setting up the overall shape functions.
   */
  std::vector<std::vector<Real>> mode;

  /**
   * the first local derivative of the radial approximation shapes.
   * Needed when setting up the overall shape functions.
   */
  std::vector<std::vector<Real>> dmodedv;

  /**
   * the radial mapping shapes in local coordinates
   */
  std::vector<std::vector<Real>> radial_map;


  /**
   * the first local derivative of the radial mapping shapes
   */
  std::vector<std::vector<Real>> dradialdv_map;

  /**
   * the first local derivative (for 3D, the first in the base)
   * of the phase term in local coordinates.
   * Needed in the overall weak form of infinite element formulations.
   */
  std::vector<Real> dphasedxi;

  /**
   * the second local derivative (for 3D, the second in the base)
   * of the phase term in local coordinates.
   * Needed in the overall weak form of infinite element formulations.
   */
  std::vector<Real> dphasedeta;

  /**
   * the third local derivative (for 3D, the derivative in radial
   * direction) of the phase term in local coordinates.
   * Needed in the overall weak form of infinite element formulations.
   */
  std::vector<Real> dphasedzeta;




  // numbering scheme maps

  /**
   * The internal structure of the \p InfFE
   * -- tensor product of base element times radial
   * nodes -- has to be determined from the node numbering
   * of the current infinite element.  This vector
   * maps the infinite \p Elem node number to the
   * radial node (either 0 or 1).
   */
  std::vector<unsigned int> _radial_node_index;

  /**
   * The internal structure of the \p InfFE
   * -- tensor product of base element times radial
   * nodes -- has to be determined from the node numbering
   * of the current element.  This vector
   * maps the infinite \p Elem node number to the
   * associated node in the base element.
   */
  std::vector<unsigned int> _base_node_index;

  /**
   * The internal structure of the \p InfFE
   * -- tensor product of base element shapes times radial
   * shapes -- has to be determined from the dof numbering
   * scheme of the current infinite element.  This vector
   * maps the infinite \p Elem dof index to the radial
   * \p InfFE shape index (\p 0..radial_order+1 ).
   */
  std::vector<unsigned int> _radial_shape_index;

  /**
   * The internal structure of the \p InfFE
   * -- tensor product of base element shapes times radial
   * shapes -- has to be determined from the dof numbering
   * scheme of the current infinite element.  This vector
   * maps the infinite \p Elem dof index to the associated
   * dof in the base \p FE.
   */
  std::vector<unsigned int> _base_shape_index;




  // some more protected members

  /**
   * The number of total approximation shape functions for
   * the current configuration
   */
  unsigned int _n_total_approx_sf;

  /**
   * The total number of quadrature points
   * for the current configuration
   */
  unsigned int _n_total_qp;

  /**
   * this vector contains the combined integration weights, so
   * that \p FEAbstract::compute_map() can still be used
   */
  std::vector<Real> _total_qrule_weights;

  /**
   * The quadrature rule for the base element associated
   * with the current infinite element
   */
  std::unique_ptr<QBase> base_qrule;

  /**
   * The quadrature rule for the base element associated
   * with the current infinite element
   */
  std::unique_ptr<QBase> radial_qrule;

  /**
   * The base element associated with the
   * current infinite element
   */
  std::unique_ptr<Elem> base_elem;

  /**
   * Have a \p FE<Dim-1,T_base> handy for base approximation.
   * Since this one is created using the \p FEBase::build() method,
   * the \p InfFE class is not required to be templated w.r.t.
   * to the base approximation shape.
   */
  std::unique_ptr<FEBase> base_fe;

  /**
   * This \p FEType stores the characteristics for which the data
   * structures \p phi, \p phi_map etc are currently initialized.
   * This avoids re-initializing the radial part.
   *
   * \note Currently only \p order may change, both the FE families
   * and \p base_order must remain constant.
   */
  FEType current_fe_type;


private:

  /**
   * \returns \p false, currently not required.
   */
  virtual bool shapes_need_reinit() const override;

  /**
   * When \p compute_node_indices_fast() is used, this static
   * variable remembers the element type for which the
   * static variables in  \p compute_node_indices_fast()
   * are currently set.  Using a class member for the
   * element type helps initializing it to a default value.
   */
  static ElemType _compute_node_indices_fast_current_elem_type;


#ifdef DEBUG

  /**
   * static members that are used to issue warning messages only once.
   */
  static bool _warned_for_nodal_soln;
  static bool _warned_for_shape;
  static bool _warned_for_dshape;

#endif

  /**
   * Make all \p InfFE<Dim,T_radial,T_map> classes
   * friends of each other, so that the protected
   * \p eval() may be accessed.
   */
  template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
  friend class InfFE;

  friend class InfFEMap;
};



// InfFERadial class inline members
inline
Real InfFERadial::decay(const unsigned int dim, const Real v)
{
  switch (dim)
    //TODO:[DD] What decay do i have in 2D and 1D?
    {
    case 3:
      return (1.-v)/2.;

    case 2:
      return 0.;

    case 1:
      return 0.;

    default:
      libmesh_error_msg("Invalid dim = " << dim);
    }
}

} // namespace libMesh


#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


#endif // LIBMESH_INF_FE_H
