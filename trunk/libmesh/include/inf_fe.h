// $Id: inf_fe.h,v 1.5 2003-01-24 17:24:38 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __inf_fe_h__
#define __inf_fe_h__

// C++ includes

// Local includes
#include "fe_base.h"
#include "elem.h"
#include "mesh_config.h"


// forward declarations




#ifdef ENABLE_INFINITE_ELEMENTS


/**
 * A specific instatiation of the \p FEBase class. This
 * class is templated, and specific template instantiations
 * will result in different Infinite Element families, similar
 * to the \p FE class.  \p InfFE builds a \p FE<Dim-1,T_base>,
 * and most of the requests related to the base are handed over
 * to this object.  All methods related to the radial part 
 * are collected in the nested class \p Radial.  Similarly,
 * most of the static methods concerning base approximation
 * are contained in \p Base.
 * 
 * Having different shape approximation families in radial direction 
 * introduces the requirement for an additional \p Order in this
 * class. Therefore, the \p FEType internals change when infinite
 * elements are used. 
 * When the specific infinite element type is not known at compile
 * time, use the \p FEBase::build() member to create abstract 
 * (but still optimized) infinite elements at run time.
 *
 * This is the numbering scheme, also applicable to the quadrature points:
   \verbatim

                      to infinity
                           ^                                   
                           |                                   

    o 5         ( )             (o)           (o)           (o)
    |            |               |                          / 2      
    o 4          |               o            o            o        
    |            |               |                        /         
    o 3          |               o 9         o ...       o          
    |            |               |                      /           
    o 2          |               o 6        o 8        o 7          
    |            |               |                    /             
    o 1          |               o 3       o 5       o 4            
    |            |               |                  /               
    o 0         ( )            ({#})----({#})----({#})              
                                  0        2        1                
  radial       radial 
  modes        mapping

                                {X}------{X}------{X}              
                                  0        2        1  


                                        base
                                      element


  o        radial approximation dof
 ( )       radial mapping node
  #        radial and base approximation dof
  X        base approximation dof
 { }       base mapping node
   \endverbatim
 *
 * \author Daniel Dreyer
 * \date 2003
 * \version $Revision: 1.5 $
 */

//-------------------------------------------------------------
// InfFE class definition
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
class InfFE : public FEBase
{

/*
 * Protect the nested class
 */
protected:

  /**
   * Infinite elements are in some sense directional, compared
   * to conventional finite elements.  All methods related 
   * to the radial part, which extends perpendicular from the base,
   * are collected in this nested class.  This class offers
   * static methods, which are only available to \p InfFE members.
   *
   * \author Daniel Dreyer
   * \date 2003
   * \version $Revision: 1.5 $
   */
  //-------------------------------------------------------------
  // InfFE::Radial class definition
  class Radial
  {
  public:

    /**
     * Constructor.  Don't use an object of this type.
     */
    Radial();

    /**
     * @returns the decay in radial direction of
     * the \p Dim dimensional infinite element.
     */
    static real decay(const real v);

    /**
     * @returns the first (local) derivative of the
     * decay in radial direction of the infinite element.
     */
    static real decay_deriv(const real) { return .5; };

    /**
     * @returns the radial weight D, used as an additional weight
     * for the test function, evaluated at local radial coordinate \p v.
     */
    static real D(const real v) { return (1.-v)*(1.-v)/4.; };

    /**
     * @returns the first (local) radial derivative of the radial weight D.
     */
    static real D_deriv(const real v) { return (v-1.)/2.; };

    /**
     * @returns the index (0 for the base, 1,2,... for the outer shells) 
     * in @e radial direction of the infinite element.  \p i is the index in
     * the whole infinite element, and \p n_base_dofs is the number of dofs
     * in the base of the infinite element.
     * Note that for  \p i<Elem::n_nodes(), i.e. a node, not a dof,
     * either 0 (base) or 1 (outer node) is returned.
     */
    static unsigned int index(const unsigned int n_base_dofs,
			      const unsigned int i)
	{ return i / n_base_dofs; };

    /**
     * @returns the index (0 for the base, 1,2,... for the outer shells) 
     * in @e radial direction for the infinite element with the @e base \p FEType
     * \p base_fe_type and the @e base element \p base_elem_type.  Used by the 
     * public static members of \p InfFE.
     */
    static unsigned int index(const FEType& base_fe_type,
			      const ElemType base_elem_type,
			      const unsigned int i);

    /**
     * @returns the Order of the mapping functions
     * in radial direction. Currently, this is @e always \p FIRST.
     */
    static Order mapping_order() { return FIRST; };

    /**
     * @returns the number of shape functions in radial direction
     * associated with this infinite element.
     * Either way, if the modes are stored as nodal dofs (\p n_dofs_at_node) 
     * or as element dofs (\p n_dofs_per_elem), in each case we have the
     * same number of modes in radial direction. Note that for the case of 1D
     * infinite elements, in the base the dof-per-node scheme is used, 
     * as defined in \p FE<0,T>.
     * 
     * From the formulation of the infinite elements, we have
     * 1 mode, when \p o_radial=CONST.
     * Therefore, we have a total of \p o_radial+1 modes in radial direction.
     */
    static unsigned int n_dofs(const ElemType /* t */,
			       const Order o_radial)
	{ return static_cast<unsigned int>(o_radial)+1; };

    /**
     * @returns the number of dofs in radial direction on "onion slice" 
     * \p n (either 0 or 1) for an infinite element of type \p inf_elem_type and 
     * radial order \p o_radial.
     *
     * Currently, the first radial mode is associated with the nodes in
     * the base.  All higher radial modes are associated with
     * the physically existing nodes further out.
     */
    static unsigned int n_dofs_at_node(const ElemType /* inf_elem_type */,
				       const Order o_radial,
				       const unsigned int n_onion);

    /**
     * @returns the number of modes in radial direction interior to the element,
     * not associated with any interior nodes.
     * Note that these modes are a discontinuous approximation, therefore
     * we have no special formulation for coupling in the base, like in the 
     * case of associating (possibly) multiple dofs per (outer) node.
     */
    static unsigned int n_dofs_per_elem(const ElemType /* inf_elem */,
					const Order    o_radial)
	{ return static_cast<unsigned int>(o_radial)+1; };
				       
    /**
     * @returns the location in radial direction (on the reference axis) 
     * of the point \p p located in physical space.  This function requires
     * inverting the (possibly nonlinear) transformation map, so
     * it is not trivial.
     */
    static Point inverse_map (const Elem* inf_elem,
			      const real dist_origin);

  };



  /**
   * This nested class contains most of the static methods related 
   * to the base part of an infinite element.  Only static members
   * are provided, and these should only be accessible from within \p InfFE.
   *
   * \author Daniel Dreyer
   * \date 2003
   * \version $Revision: 1.5 $
   */
  //-------------------------------------------------------------
  // InfFE::Base class definition
  class Base
  {
  public:

    /**
     * Constructor.  Don't use an object of this type.
     */
    Base();

    /**
     * Build the base element of an infinite element.  Be careful,
     * this method allocates memory!  So be sure to delete the
     * new element afterwards.
     */
    static Elem* build_elem (const Elem* inf_elem)
	{ AutoPtr<Elem> ape = inf_elem->build_side(0); return ape.release(); };

    /**
     * Build the \p FEType for the base element of an infinite element,
     * where \p inf_fe_type is the \p FEType of the whole infinite element.
     */
    static FEType build_fe_type (const FEType& inf_fe_type);
 
    /**
     * @returns the index in the @e base element. \p i
     * is the index in the whole infinite element,
     * \p n_base_dofs is the number of dofs in the @e base
     * of the infinite element.
     */
    static unsigned int index(const unsigned int n_base_dofs,
			      const unsigned int i)
	{ return i % n_base_dofs; };

    /**
     * @returns the index in the @e base element. \p i is the index in the 
     * whole infinite element, \p base_fe_type is the  @e base \p FEType
     * and \p base_elem_type is the the @e base element.  Used by the 
     * public static members of \p InfFE.
     */
    static unsigned int index(const FEType& base_fe_type,
			      const ElemType base_elem_type,
			      const unsigned int i);

    /**
     * @returns the base element associated to
     * \p type.  This is, for example, \p TRI3 for
     * \p INFPRISM6.
     */
    static ElemType get_elem_type(const ElemType type);


  };







public:

  //-------------------------------------------------------------
  // InfFE continued

  /**
   * Constructor.
   * Initializes some data structures, and builds \p base_fe with 
   * \p fet.base_order as \p order, and \p fet.base_family as 
   * \p family.  Through this we get the appropriate \p base_fe* 
   * for handling approximation in the base, avoiding to template
   * InfFE also with respect to the base approximation \p FEFamily.
   * This is only possible through the \p FEBase::build() method
   * and some protected additional virtual members in \p FEBase.
   */
  InfFE(const FEType& fet);

  /**
   * Desctructor.  Clean up.
   */
  ~InfFE();





  //-------------------------------------------------------------
  // The static public members for access from FEInterface etc
  /**
   * @returns the value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method lets you specify the relevant
   * data directly, and is therefore allowed to be static.
   * Note that this class member is by far not as efficient as 
   * its counterpart in \p FE<Dim,T>, and is @e not employed
   * in the \p reinit() cycle.
   */
  static real shape(const FEType& fet,
		    const ElemType t,
		    const unsigned int i,
		    const Point& p);

  /**
   * @returns the value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method lets you specify the relevant
   * data directly, and is therefore allowed to be static.
   * Note that this class member is not as efficient as its 
   * counterpart in \p FE<Dim,T>, and is @e not employed 
   * in the \p reinit() cycle.
   */
  static real shape(const FEType& fet,
		    const Elem* elem,
		    const unsigned int i,
		    const Point& p);
  
  /**
   * @returns the number of shape functions associated with
   * a finite element of type \p t and approximation order \p o.
   */
  static unsigned int n_shape_functions (const FEType& fet,
					 const ElemType t);

  /**
   * @returns the number of shape functions associated with this
   * infinite element.  Currently, we have \p o_radial+1 modes in 
   * radial direction, and \p FE<Dim-1,T>::n_dofs(...) in the base.
   */
  static unsigned int n_dofs(const FEType& fet,
			     const ElemType inf_elem_type);

  /**
   * @returns the number of dofs at infinite element @e node \p n 
   * (not dof!) for an element of type \p t and order \p o.  
   */
  static unsigned int n_dofs_at_node(const FEType& fet,
				     const ElemType inf_elem_type,
				     const unsigned int n);

  /**
   * @returns the number of dofs interior to the element,
   * not associated with any interior nodes.
   */
  static unsigned int n_dofs_per_elem(const FEType& fet,
				      const ElemType inf_elem_type);

  /**
   * Build the nodal soln from the element soln.
   * This is the solution that will be plotted.
   */
  static void nodal_soln(const FEType& fet,
			 const Elem* elem, 
			 const std::vector<number>& elem_soln,
			 std::vector<number>& nodal_soln);

  /**
   * @returns the location (on the reference element) of the
   * point \p p located in physical space.  This function requires
   * inverting the (possibly nonlinear) transformation map, so
   * it is not trivial.
   */
  static Point inverse_map (const Elem* elem,
			    const Point& p);

  /**
   * @returns true if \p et is an element to be processed by
   * class \p InfFE.  Otherwise, it returns false, and
   * this element should be processed using \p FE.
   * This method is particularly helpful during the actual
   * matrix assembly process.
   */
  static bool is_InfFE_elem(const ElemType et);




  //-------------------------------------------------------------
  // The work-horses of InfFE. These are often used during matrix assembly
  /**
   * This is at the core of this class. Use this for each
   * new element in the mesh.  Reinitializes all the physical 
   * element-dependent data based on the current element 
   * \p elem.
   */
  void reinit (const Elem* elem);
    
  /**
   * Not implemented yet.  Reinitializes all the physical 
   * element-dependent data based on the \p side of an infinite 
   * element.
   */
  void reinit (QBase* qside,
	       const Elem* elem,
	       const unsigned int side);

  /**
   * The use of quadrature rules with the \p InfFE class is somewhat
   * different from the approach of the \p FE class.  While the
   * \p FE class requires an appropriately initialized quadrature
   * rule object, and simply uses it, the \p InfFE class requires only
   * the quadrature rule object of the current \p FE class.
   * From this \p QBase*, it determines the necessary data,
   * and @e builds two appropriate quadrature classes, one for radial,
   * and another for base integration.  Currently, for both
   * base and radial integration only Gaussian quadrature is 
   * available, but may be changed easily.
   */
  void attach_quadrature_rule (QBase* q);
  
  /**
   * @returns the total number of quadrature points.  Call this
   * to get an upper bound for the \p for loop in your simulation
   * for matrix assembly of the current element.
   */
  unsigned int n_quadrature_points () const
      { assert (radial_qrule != NULL); return _n_total_qp; };

  /**
   * @returns the global first derivative of the phase term in
   * infinite elements, evaluated at the quadrature points.
   */
  const std::vector<Point>& get_dphase() const
      { return dphase; };

  /**
   * @returns the multiplicative weight at each quadrature point.
   * This weight is used for certain infinite element weak 
   * formulations, so that @e weighted Sobolev spaces are
   * used for the trial function space.  This renders the
   * variational form easily computable.
   */
  const std::vector<real>& get_Sobolev_weight() const
      { return weight; };

  /**
   * @returns the first global derivative of the multiplicative 
   * weight at each quadrature point. See \p get_Sobolev_weight()
   * for details.
   */
  const std::vector<Point>& get_Sobolev_dweight() const
      { return dweight; };

  /**
   * @returns the number of shape functions associated with
   * this infinite element.
   */
  unsigned int n_shape_functions () const { return _n_total_approx_sf; };




protected:

  //-------------------------------------------------------------
  // static members used by the "work-horses"
  /**
   * @returns the value of the \f$ i^{th} \f$ polynomial evaluated
   * at \p v.  This method provides the approximation
   * in radial direction for the overall shape functions,
   * which is defined in \p InfFE::shape().
   * This method is allowed to be static, since it is independent
   * of dimension and base_family.  It is templated, though,
   * w.r.t. to radial \p FEFamily.
   *
   * Specialized for \p T_radial=INFINITE_MAP, this function returns 
   * the value of the \f$ i^{th} \f$ @e mapping shape function
   * in radial direction evaluated at \p v.  Currently, only one specific 
   * mapping shape is used.  Namely the one by Marques JMMC, Owen DRJ:
   * Infinite elements in quasi-static materially nonlinear problems,
   * @e Computers @e and @e Structures, 1984.
   */
  static real eval(const real v,
		   const Order o_radial,
		   const unsigned int i);
  
  /**
   * @returns the value of the first derivative of the
   * \f$ i^{th} \f$ polynomial at coordinate \p v.
   * See \p eval for details.
   */
  static real eval_deriv(const real v,
			 const Order o_radial,
			 const unsigned int i);



//protected:

  //-------------------------------------------------------------
  // Non-static members used by the "work-horses"
  /**
   * Updates the protected member \p base_elem to the appropriate base element
   * for the given \p inf_elem.
   */
  void update_base_elem (const Elem* inf_elem);

  /**
   * Do not use this derived member in \p InfFE<Dim,T_radial,T_map>.
   */
  void init_base_shape_functions(const QBase*, const Elem*)
      { error(); };



  //-------------------------------------------------------------
  // The bigger internal methods used during assembly  
  /** 
   * Initialize all the data fields like \p weight, \p mode, 
   * \p dist, \p phi, \p dphidxi, \p dphideta, \p dphidzeta, etc.
   * for the current element.  Of these data fields, only
   * the ones that are independent of base approximation
   * are evaluated.  For constant radial \p order in the mesh, 
   * this method only has to be called once.
   */
  void init_shape_functions(const Elem* inf_elem);

  /** 
   * Not implemented yet.  Initialize all the data fields like \p weight, 
   * \p phi, etc for the side \p s.
   */  
  void init_shape_functions(const QBase* q,
			    const Elem* e,
			    const unsigned int s);

  /** 
   * Combines the base approximation, mapping etc. with
   * the radial counterparts.  Has to be used every time
   * the base approximation changes (so rather often).
   */
  void combine_base_radial();

  /** 
   * After having updated the jacobian and the transformation
   * from local to global coordinates in FEBase::compute_map(),
   * the first derivatives of the shape functions are 
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx/y/z, \p dphasedx/y/z, \p dweight. This method
   * should barely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected.
   * Do not give any values to this method, from this
   * it is even more visible how we use the \p InfFE 
   * and not the \p FEBase version.
   */
  void compute_shape_functions();



  //TODO:[DD] overload compute_map?
  
/*   NOT YET IMPLEMENTED */
/*   /\** */
/*    * @returns the location (in physical space) of the point */
/*    * \p p located on the reference element. */
/*    *\/ */
/*   static Point map (const Elem* elem, */
/* 		    const Point& reference_point); */
  
/*   /\** */
/*    * @returns d(xyz)/dxi (in physical space) of the point */
/*    * \p p located on the reference element. */
/*    *\/ */
/*   static Point map_xi (const Elem* elem, */
/* 		       const Point& reference_point); */
  
/*   /\** */
/*    * @returns d(xyz)/deta (in physical space) of the point */
/*    * \p p located on the reference element. */
/*    *\/ */
/*   static Point map_eta (const Elem* elem, */
/* 			const Point& reference_point); */

/*   /\** */
/*    * @returns d(xyz)/dzeta (in physical space) of the point */
/*    * \p p located on the reference element. */
/*    *\/ */
/*   static Point map_zeta (const Elem* elem, */
/* 			 const Point& reference_point); */





  //--------------------------------------------------------------
  // protected members, which may be accessed from the outside through some inline functions
  /**
   * the first derivatives of the phase term in global coordinates,
   * over @e all quadrature points.
   */
  std::vector<Point> dphase;

  /**
   * the global derivative of the additional radial weight \f$ 1/{r^2} \f$,
   * over @e all quadrature points.
   */
  std::vector<Point> dweight;

  /**
   * the additional radial weight \f$ 1/{r^2} \f$ in local coordinates,
   * over @e all quadrature points.
   */
  std::vector<real>  weight;



  //--------------------------------------------------------------
  // protected members, which are not to be accessed from outside
  /**
   * the radial distance of the base nodes from the origin
   */
  std::vector<real>  dist;

  /**
   * the additional radial weight \f$ 1/{r^2} \f$ in local coordinates,
   * over all quadrature points. The weight does not vary in base
   * direction.  However, for uniform access to the data fields from the 
   * outside, this data field is expanded to @e all quadrature points.
   */
  std::vector<real>  dweightdv;

  /**
   * the radial decay \f$ 1/r \f$ in local coordinates.
   * Needed when setting up the overall shape functions.
   * Note that it is this decay which assures to satisfy
   * the Sommerfeld radiation condition in advance.
   */
  std::vector<real>  som;
  /**
   * the first local derivative of the radial decay \f$ 1/r \f$ in local 
   * coordinates.  Needed when setting up the overall shape functions.
   */
  std::vector<real>  dsomdv;

  /**
   * the radial approximation shapes in local coordinates
   * Needed when setting up the overall shape functions.
   */
  std::vector<std::vector<real> >   mode;

  /**
   * the first local derivative of the radial approximation shapes.
   * Needed when setting up the overall shape functions.
   */
  std::vector<std::vector<real> >   dmodedv;
  
  /**
   * the radial mapping shapes in local coordinates
   */
  std::vector<std::vector<real> >   radial_map;


  /**
   * the first local derivative of the radial mapping shapes
   */
  std::vector<std::vector<real> >   dradialdv_map;

  /**
   * the first local derivative (for 3D, the first in the base) 
   * of the phase term in local coordinates.
   * Needed in the overall weak form of infinite element formulations.
   */
  std::vector<real>  dphasedxi;

  /**
   * the second local derivative (for 3D, the second in the base) 
   * of the phase term in local coordinates.
   * Needed in the overall weak form of infinite element formulations.
   */
  std::vector<real>  dphasedeta;

  /**
   * the third local derivative (for 3D, the derivative in radial
   * direction) of the phase term in local coordinates.
   * Needed in the overall weak form of infinite element formulations.
   */
  std::vector<real>  dphasedzeta;




  //--------------------------------------------------------------
  // some protected members

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
   * The quadrature rule for the base element associated 
   * with the current infinite element
   */
  QBase* base_qrule;

  /**
   * The quadrature rule for the base element associated 
   * with the current infinite element
   */
  QBase* radial_qrule;

  /**
   * The base element associated with the
   * current infinite element
   */
  Elem* base_elem;
//  AutoPtr<Elem> base_elem;

  /**
   * Have a \p FE<Dim-1,T_base> handy for base approximation.
   * Since this one is created using the \p FEBase::build() method,
   * the \p InfFE class is not required to be templated w.r.t.
   * to the base approximation shape.
   */
  FEBase* base_fe;

  /**
   * This \p FEType stores the characteristics for which
   * the data structures \p phi, \p phi_map etc are currently
   * initialized.  This avoids re-initializing the radial
   * part.  But note that currently @e only \p order may change,
   * neither the FE families nor \p base_order!
   */
  FEType current_fe_type;


private:
 
  /**
   * Make the \p FEBase class a friend so that its
   * \p FEBase::build() member will work.
   */
  friend class FEBase;

  /**
   * Make all \p \InfFE<Dim,T_radial,T_map> classes
   * friends of each other, so that the protected
   * \p eval() may be accessed.
   */
  template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
  friend class InfFE;

};




// ------------------------------------------------------------
// InfFE class inline members
template <unsigned int Dim, FEFamily T_radial, FEFamily T_base>
inline
bool InfFE<Dim,T_radial,T_base>::is_InfFE_elem(const ElemType et)
{

  switch (et)
  {
    case INFEDGE2:
    case INFQUAD4:
    case INFQUAD6:
    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
    case INFPRISM6:
    case INFPRISM12:
      {
        return true;
      };

    default:
      { 
	return false;
      };

  };

};



template <unsigned int Dim, FEFamily T_radial, FEFamily T_base>
inline
void InfFE<Dim,T_radial,T_base>::update_base_elem (const Elem* inf_elem)
{
  if (base_elem != NULL)
    delete base_elem;
  base_elem = Base::build_elem(inf_elem);
};



// ------------------------------------------------------------
// InfFE::Radial class inline members
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
inline
real InfFE<Dim,T_radial,T_map>::Radial::decay(const real v)
{
  switch (Dim)
  //TODO:[DD] What decay do i have in 2D and 1D?
  {
    case 3:
      return (1.-v)/2.;
    case 2:
      return 0.;
    case 1:
      return 0.;
    default:
      error();
      return 0.;
  };
};



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
inline
unsigned int InfFE<Dim,T_radial,T_map>::Radial::n_dofs_at_node(const ElemType,
							       const Order o_radial,
							       const unsigned int n_onion)
{
  assert (n_onion < 2);

  if (n_onion == 0)
    /*
     * in the base, no matter what, we have 1 node associated 
     * with radial direction
     */
    return 1;
  else
    /*
     * this works, since for Order o_radial=CONST=0, we still
     * have the (1-v)/2 mode, associated to the base
     */
    return static_cast<unsigned int>(o_radial);
};





// ------------------------------------------------------------
// InfFE::Base class inline members
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
inline
FEType InfFE<Dim,T_radial,T_map>::Base::build_fe_type (const FEType& inf_fe_type)
{
  FEType tmp(INVALID_ORDER, 
	     INVALID_FE, 
	     inf_fe_type.base_order, 
	     inf_fe_type.base_family, 
	     INVALID_INF_MAP );
  return tmp; 
};





#endif //ifdef ENABLE_INFINITE_ELEMENTS


#endif
