// $Id: inf_fe_static.C,v 1.13 2003-03-29 12:45:33 ddreyer Exp $

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



// Local includes
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS
#include "inf_fe.h"
#include "fe_interface.h"
#include "elem.h"




// ------------------------------------------------------------
// InfFE::Radial class members
template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
InfFE<Dim,T_radial,T_base>::Radial::Radial () 
{ 
  std::cerr << "Do not define an object of this type." 
	    << std::endl;  
  error(); 
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::Radial::index(const FEType& fe_type,
						      const ElemType base_elem_type,
						      const unsigned int i)
{
  if (Dim > 1)
    return Radial::index( FEInterface::n_dofs(Dim-1, fe_type, base_elem_type), i);
  else
    return Radial::index(1, i);
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Point InfFE<Dim,T_radial,T_map>::Radial::inverse_map(const Elem*,
						     const Real,
						     const Real)

{
  std::cerr << "ERROR: Radial::inverse_map() not yet implemented." << std::endl;
  error();

  Point p;
  return p;
}








// ------------------------------------------------------------
// InfFE::Base class members
template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
InfFE<Dim,T_radial,T_base>::Base::Base () 
{ 
  std::cerr << "Do not define an object of this type." 
	    << std::endl;  
  error(); 
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
Elem* InfFE<Dim,T_radial,T_base>::Base::build_elem (const Elem* inf_elem)
{ 
  AutoPtr<Elem> ape(inf_elem->build_side(0)); 
  return ape.release(); 
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::Base::index(const FEType& fe_type,
						    const ElemType base_elem_type,
						    const unsigned int i)
{
  if (Dim > 1)
    return Base::index(FEInterface::n_dofs(Dim-1, fe_type, base_elem_type), i);
  else
    return Base::index(1, i);
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
ElemType InfFE<Dim,T_radial,T_base>::Base::get_elem_type(const ElemType type)
{
  switch (type)
    {
      // 3D infinite elements:
      // with Dim=3 -> infinite elements on their own
      case INFHEX8:
	  return QUAD4;

      case INFHEX16:
	  return QUAD8;

      case INFHEX18:
	  return QUAD9;
		 
      case INFPRISM6:
	  return TRI3;

      case INFPRISM12:
	  return TRI6;

      // 2D infinite elements:
      // with Dim=3 -> used as boundary condition,
      // with Dim=2 -> infinite elements on their own
      case INFQUAD4:
	  return EDGE2;

      case INFQUAD6:
	  return EDGE3;

      // 1D infinite elements:
      // with Dim=2 -> used as boundary condition,
      // with Dim=1 -> infinite elements on their own,
      //               but no base element!
      case INFEDGE2:
	  return INVALID_ELEM;

      default:
	{
	  std::cerr << "ERROR: Unsupported element type!: " << type
		    << std::endl;
	  error();
	}
    }


  error();
  return INVALID_ELEM;
}







// ------------------------------------------------------------
// InfFE public static class members; mostly for interaction 
// from FEInterface



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_shape_functions (const FEType& fet,
							   const ElemType t)
{
  return InfFE<Dim,T_radial,T_map>::n_dofs(fet, t);
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs(const FEType& fet,
					       const ElemType inf_elem_type)
{
  const ElemType base_et(Base::get_elem_type(inf_elem_type));
    
  if (Dim > 1)
    return FEInterface::n_dofs(Dim-1, fet, base_et) * Radial::n_dofs(fet.radial_order);
  else
    return Radial::n_dofs(fet.radial_order);
}
		




template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs_at_node(const FEType& fet,
						       const ElemType inf_elem_type,
						       const unsigned int n)
{
  const ElemType     base_et  ( Base::get_elem_type(inf_elem_type) );

  const unsigned int n_base   ( Base::index  (fet, base_et, n) );
  const unsigned int n_radial ( Radial::index(fet, base_et, n) );

  if (Dim > 1)
    return FEInterface::n_dofs_at_node(Dim-1, fet, base_et, n_base) 
        * Radial::n_dofs_at_node(fet.radial_order, n_radial);
  else
    return Radial::n_dofs_at_node(fet.radial_order, n_radial);
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs_per_elem(const FEType& fet,
							const ElemType inf_elem_type)
{
  const ElemType     base_et  ( Base::get_elem_type(inf_elem_type) );

  if (Dim > 1)
    return FEInterface::n_dofs_per_elem(Dim-1, fet, base_et) 
        * Radial::n_dofs_per_elem(fet.radial_order);
  else
    return Radial::n_dofs_per_elem(fet.radial_order);
}







template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::nodal_soln(const FEType&,
					   const Elem*,
					   const std::vector<Number>&,
					   std::vector<Number>&)
{
  std::cerr << "ERROR: The concept of a nodal solution is not "
	    << "applicable to infinite elements!" << std::endl;
  error();  
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Point InfFE<Dim,T_radial,T_map>::inverse_map (const Elem*,
					      const Point& physical_point,
					      const Real)
{
//To do: fix inverse_map() for all three dimensions, specialize for each Dim, so that this gets more effective (for field point processing...)

/*
determine origin by simply back-computing...
_or_, put origin in the constructor...
For sure: 
1. #ifdef DEBUG (see fe_map.C)
   check with distance, whether this one is in my element (needed?),
   #endif
2. project this point down to the base, and make inverse_map() for
   this 2D-element
3. get the accurate distance afterwards (or through some hand-calculated
   inversion of the 1/r-map?
*/

  return physical_point;
}



//To do: Probably have to fix on_reference_element() also for InfFE
/*
bool FEBase::on_reference_element(const Point& p, const ElemType t, const Real eps)
*/






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape(const FEType& fet,
				      const ElemType type,
				      const unsigned int i,
				      const Point& p)
{
  assert (Dim != 0);

  const ElemType     base_et  ( Base::get_elem_type(type) );
  const unsigned int i_base   ( Base::index  (fet, base_et, i) );

  const Order        o_radial ( fet.radial_order );
  const unsigned int i_radial ( Radial::index(fet, base_et, i) );
  const Real         v        ( p(Dim-1) );  // holds for all Dim, except for 0, but we don't inst Dim=0 ;-)


  //TODO:[SP/DD]  exp(ikr) is still missing here!
  if (Dim > 1)
    return FEInterface::shape (Dim-1, fet, base_et, i_base, p)
        * InfFE<Dim,T_radial,T_map>::eval (v, o_radial, i_radial)
        * InfFE<Dim,T_radial,T_map>::Radial::decay(v);
  else
    return InfFE<Dim,T_radial,T_map>::eval (v, o_radial, i_radial)
        * InfFE<Dim,T_radial,T_map>::Radial::decay(v);
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape(const FEType& fet,
				      const Elem* elem,
				      const unsigned int i,
				      const Point& p)
{
  assert (elem != NULL);
  assert (Dim != 0);

  const ElemType     base_et (Base::get_elem_type(elem->type()));
  const unsigned int i_base  (Base::index  (fet, base_et, i));

  const Order        o_radial(fet.radial_order);
  const unsigned int i_radial(Radial::index(fet, base_et, i));
  const Real         v       (p(Dim-1));  // holds for all Dim, except for 0, but we don't inst Dim=0 ;-)

  AutoPtr<Elem>      base_el(elem->build_side(0));

  if (Dim > 1)
    return FEInterface::shape (Dim-1, fet, base_el.get(), i_base, p)
        * InfFE<Dim,T_radial,T_map>::eval (v, o_radial, i_radial)
        * InfFE<Dim,T_radial,T_map>::Radial::decay(v);
  else
    return InfFE<Dim,T_radial,T_map>::eval (v, o_radial, i_radial)
        * InfFE<Dim,T_radial,T_map>::Radial::decay(v);
}





//--------------------------------------------------------------
// Explicit instantiations
#include "inf_fe_instantiate_1D.h"
#include "inf_fe_instantiate_2D.h"
#include "inf_fe_instantiate_3D.h"



#endif //ifdef ENABLE_INFINITE_ELEMENTS

