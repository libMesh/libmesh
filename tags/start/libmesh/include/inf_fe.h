// $Id: inf_fe.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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


// forward declarations



// Local includes
#include "fe.h"
//#include "polynomial.h"




#ifdef ENABLE_INFINITE_ELEMENTS



/**
 * A specific instatiation of the \p FEBase class.  
 * This creates infinite elements used in exterior acoustics. 
 * Currently, the Astley-Leis, or Mapped Wave Envelope Elements
 * are implemented. Different polynomials for approximation
 * in radial direction are provided by the polynomial class
 * attached through attach_polynomial().
 * 
 * @author Daniel Dreyer, 2002
 */

//-------------------------------------------------------------
// InfFE class definition

class InfFE : public FEBase
{
public:

  /**
   * Constructor.  Optionally initializes required data
   * structures.
   */
  InfFE(const MeshBase& m,
	const unsigned int d,
	const FEType& fet,
	const ElemType t=INVALID_ELEM,
	QBase* q=NULL);//,
//     PolynomialBase* p=NULL);
  

  /**
   * Destructor.
   */
  virtual ~InfFE() {return;};

  /**
   * @returns the number of shape functions associated with this
   * finite element of type t.
   */
  unsigned int n_shape_functions() const;

  /**
   * @returns the number of shape functions associated with this
   * finite element of type t, order o.
   */
  static unsigned int n_shape_functions(const ElemType t,
					const Order o);

  /**
   * @returns the number of shape functions associated with this
   * finite element.
   */
  static unsigned int n_dofs(const ElemType t,
			     const Order o);

  /**
   * @returns the number of dofs at node n for a finite element
   * of type t and order o.
   */
  static unsigned int n_dofs_at_node(const ElemType t,
				     const Order o,
				     const unsigned int n);

  /**
   * @returns the number of dofs interior to the element,
   * not associated with any interior nodes.
   */
  static unsigned int n_dofs_per_elem(const ElemType t,
				      const Order o);
				     
  
  /**
   * Build the nodal soln from the element soln.
   * This is the solution that will be plotted.
   */
  static void nodal_soln(const MeshBase& mesh,
			 const Elem* elem, const Order o,
			 const std::vector<number>& elem_soln,
			 std::vector<number>& nodal_soln);

  /**
   * @returns true if the point p is located on the reference element
   * for element type t, false otherwise.  Since we are doing floating
   * point comparisons here the parameter eps can be specified to
   * indicate a tolerance.  For example, (x<=1) becomes (x<=1+eps). 
   */
  static bool on_reference_element(const Point& p,
				   const ElemType t,
				   const real eps=1.e-6);

  /**
   * @returns the value of the \f$ i^{th} \f$ shape
   * function at point \p p.
   */
  static real shape(const unsigned int d,
		    const ElemType t,
		    const Order o,
		    const unsigned int i,
		    const Point& p);

  /**
   * @returns the value of the \f$ i^{th} \f$ shape
   * function at point \p p.
   */
  static real shape(const unsigned int d,
		    const Elem* elem,
		    const Order o,
		    const unsigned int i,
		    const Point& p);

  /**
   * @returns the first derivative in \f$ j^{th} \f$
   * direction of the \f$ i^{th} \f$ shape function
   * at point \p p.
   */
  static real shape_deriv(const unsigned int d,
			  const ElemType t,
			  const Order o,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p);

  /**
   * @returns the first derivative in \f$ j^{th} \f$
   * direction of the \f$ i^{th} \f$ shape function
   * at point \p p.
   */
  static real shape_deriv(const unsigned int d,
			  const Elem* elem,
			  const Order o,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p);


  /**
   * This is at the core of this class. Use this for each
   * new element in the mesh.  Reinitializes all the physical 
   * element-dependent data based on the current element 
   * \p elem.
   */
  void reinit (const Elem* elem);
  
  
  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of \p face.
   */
  void reinit (QBase* qside,
	       const Elem* elem,
	       const unsigned int side);
  

//  /**
//   *Use this method to attach a polynomial base for radial approximation.
//  */
//  void attach_polynomial(PolynomialBase* p)
//  { poly = p; return; };

  /**
   * Prints relevant information about the current configuration.
   */ 
  void print_info() const;




protected:

//  PolynomialBase* poly;


  /** 
   * Update the various member data fields \p phi,
   * \p dphidxi, \p dphideta, \p dphidzeta, etc.
   * for the current element.
   */
  void init_shape_functions(const QBase* q,
			    const Elem* e);

  /** 
   * Same as before, but for a side.
   */  
  void init_shape_functions(const QBase* q,
			    const Elem* e,
			    const unsigned int s);


  /**
   * Compute the jacobian and some other additional
   * data fields. Currently, only conventionally mapped
   * elements are used. However, using spherical/ellipsoidal
   * mappings is quite doable. Just modify this method.
   */
  void compute_map(const QBase* q,
		   const Elem* e);

  
  /** 
   * Same as before, but for a side.
   */  
  void compute_map(const QBase* q,
		   const Elem* e,
		   const unsigned int s);


private:



};




// ------------------------------------------------------------
// InfFE class inline members
inline
InfFE::InfFE(const MeshBase& m,
	     const unsigned int d,
	     const FEType& fet,
	     const ElemType t,
	     QBase* q) : //,
//	     PolynomialBase* p) :
    FEBase(m,d,fet)//,
//    poly(p)
{
};




inline
void InfFE::print_info() const
{
//  std::cout << "Type of radial approximation." << std::endl;
//  std::cout << poly->get_name() << std::endl;
  FEBase::print_info();
}



// endif for ENABLE_INFINITE_ELEMENTS
#endif


// endif for __inf_fe_h__
#endif

