// $Id: fe_interface.C,v 1.2 2003-01-20 16:31:33 jwpeterson Exp $

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
#include "fe_interface.h"
#include "elem.h"
#include "fe.h"




//------------------------------------------------------------
//FEInterface class members
FEInterface::FEInterface()
{
  std::cerr << "Currently, this class is intended only" << std::endl
            << "to provide static member functions." << std::endl
            << "Do not define an object of this type." << std::endl;
  error();
};




#ifdef ENABLE_INFINITE_ELEMENTS

// this one decides whether we have an infinite element or not
inline
bool FEInterface::is_InfFE_elem(const ElemType t)
{

  switch (t)
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
      }

    default:
      { 
	return false;
      }

    }

};

#endif




unsigned int FEInterface::n_shape_functions(const unsigned int dim,
					    const FEType& fe_t,
					    const ElemType t)
{
  const Order o = fe_t.order;
  
  // #ifdef ENABLE_INFINITE_ELEMENTS
  //   if (FEInterface::is_InfFE_elem(t))
  //     {
  //       return InfFE::n_shape_functions(t, o);
  //     };
  // #endif

  // by default, use this
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::n_shape_functions(t, o);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::n_shape_functions(t, o);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::n_shape_functions(t, o);

	  default:
	    error();
	  };
      };

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::n_shape_functions(t, o);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::n_shape_functions(t, o);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::n_shape_functions(t, o);

	  default:
	    error();
	  };
      };

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::n_shape_functions(t, o);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::n_shape_functions(t, o);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::n_shape_functions(t, o);

	  default:
	    error();
	  };
      };


    default:
      error();
    };

  
  error();
  return 0;
};





unsigned int FEInterface::n_dofs(const unsigned int dim,
				 const FEType& fe_t,
				 const ElemType t)
{
  const Order o = fe_t.order;
  
  // #ifdef ENABLE_INFINITE_ELEMENTS
  //   if (FEInterface::is_InfFE_elem(t))
  //     {
  //       return InfFE::n_dofs(t, o);
  //     };
  // #endif

  // by default, use this
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::n_dofs(t, o);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::n_dofs(t, o);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::n_dofs(t, o);

	  default:
	    error();
	  };
      };

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::n_dofs(t, o);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::n_dofs(t, o);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::n_dofs(t, o);

	  default:
	    error();
	  };
      };

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::n_dofs(t, o);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::n_dofs(t, o);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::n_dofs(t, o);

	  default:
	    error();
	  };
      };


    default:
      error();
    };

  
  error();
  return 0;
};

		


unsigned int FEInterface::n_dofs_at_node(const unsigned int dim,
					 const FEType& fe_t,
					 const ElemType t,
					 const unsigned int n)
{
  const Order o = fe_t.order;
  
  // #ifdef ENABLE_INFINITE_ELEMENTS
  //   if (FEInterface::is_InfFE_elem(t))
  //     {
  //       return InfFE::n_dofs_at_node(t, o, n);
  //     };
  // #endif

  // by default, use this
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::n_dofs_at_node(t, o, n);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::n_dofs_at_node(t, o, n);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::n_dofs_at_node(t, o, n);

	  default:
	    error();
	  };
      };

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::n_dofs_at_node(t, o, n);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::n_dofs_at_node(t, o, n);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::n_dofs_at_node(t, o, n);

	  default:
	    error();
	  };
      };

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::n_dofs_at_node(t, o, n);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::n_dofs_at_node(t, o, n);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::n_dofs_at_node(t, o, n);

	  default:
	    error();
	  };
      };


    default:
      error();
    };

  
  error();
  return 0;
};





unsigned int FEInterface::n_dofs_per_elem(const unsigned int dim,
					  const FEType& fe_t,
					  const ElemType t)
{
  const Order o = fe_t.order;
  
  // #ifdef ENABLE_INFINITE_ELEMENTS
  //   if (FEInterface::is_InfFE_elem(t))
  //     {
  //       return InfFE::n_dofs_per_elem(t, o);
  //     };
  // #endif

  // by default, use this
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::n_dofs_per_elem(t, o);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::n_dofs_per_elem(t, o);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::n_dofs_per_elem(t, o);

	  default:
	    error();
	  };
      };

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::n_dofs_per_elem(t, o);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::n_dofs_per_elem(t, o);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::n_dofs_per_elem(t, o);

	  default:
	    error();
	  };
      };

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::n_dofs_per_elem(t, o);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::n_dofs_per_elem(t, o);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::n_dofs_per_elem(t, o);

	  default:
	    error();
	  };
      };


    default:
      error();
    };

  
  error();
  return 0;
};






void FEInterface::nodal_soln(const unsigned int dim,
			     const FEType& fe_t,
			     const Elem* elem,
			     const std::vector<number>& elem_soln,
			     std::vector<number>&       nodal_soln)
{
  const Order order = fe_t.order;
  
  // #ifdef ENABLE_INFINITE_ELEMENTS
  //   if (FEInterface::is_InfFE_elem(elem->type()))
  //     {
  //       InfFE::nodal_soln(mesh, elem, order, elem_soln, nodal_soln);
  //     }
  //   else
  // #endif

  // by default, use this
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    FE<1,HIERARCHIC>::nodal_soln(elem, order,
					 elem_soln, nodal_soln);
	    return;

	  case LAGRANGE:
	    FE<1,LAGRANGE>::nodal_soln(elem, order,
				       elem_soln, nodal_soln);
	    return;
   
	  case MONOMIAL:
	    FE<1,MONOMIAL>::nodal_soln(elem, order,
				       elem_soln, nodal_soln);
	    return;

	  default:
	    error();
	  };
      };

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    FE<2,HIERARCHIC>::nodal_soln(elem, order,
					 elem_soln, nodal_soln);
	    return;
	    
	  case LAGRANGE:
	    FE<2,LAGRANGE>::nodal_soln(elem, order,
				       elem_soln, nodal_soln);
	    return;
   
	  case MONOMIAL:
	    FE<2,MONOMIAL>::nodal_soln(elem, order,
				       elem_soln, nodal_soln);
	    return;

	  default:
	    error();
	  };
      };

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    FE<3,HIERARCHIC>::nodal_soln(elem, order,
					 elem_soln, nodal_soln);
	    return;
	    
	  case LAGRANGE:
	    FE<3,LAGRANGE>::nodal_soln(elem, order,
				       elem_soln, nodal_soln);
	    return;
	    
	  case MONOMIAL:
	    FE<3,MONOMIAL>::nodal_soln(elem, order,
				       elem_soln, nodal_soln);
	    return;

	  default:
	    error();
	  };
      };


    default:
      error();
    };

  
  error();
  return;
};




Point FEInterface::inverse_map (const unsigned int dim,
				const FEType& fe_t,
				const Elem* elem,
				const Point& p)
{
  // #ifdef ENABLE_INFINITE_ELEMENTS
  //   if (FEInterface::is_InfFE_elem(elem->type()))
  //     {
  //       return 0.;//InfFE::inverse_map(elem, p);
  //     };
  // #endif

  // by default, use this
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::inverse_map(elem, p);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::inverse_map(elem, p);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::inverse_map(elem, p);

	  default:
	    error();
	  };
      };

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::inverse_map(elem, p);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::inverse_map(elem, p);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::inverse_map(elem, p);

	  default:
	    error();
	  };
      };

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::inverse_map(elem, p);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::inverse_map(elem, p);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::inverse_map(elem, p);

	  default:
	    error();
	  };
      };


    default:
      error();
    };

  
  error();
  Point pt;
  return pt;
};



bool FEInterface::on_reference_element(const Point& p,
				       const ElemType t,
				       const real eps)
{
  return FEBase::on_reference_element(p,t,eps);
};




real FEInterface::shape(const unsigned int dim,
			const FEType& fe_t,
			const ElemType t,
			const unsigned int i,
			const Point& p)
{
  const Order o = fe_t.order;
  
  // #ifdef ENABLE_INFINITE_ELEMENTS
  //   if (FEInterface::is_InfFE_elem(t))
  //     {
  //       return 0.;//InfFE::shape(d, t, o, i, p);
  //     };
  // #endif

  // by default, use this
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::shape(t,o,i,p);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::shape(t,o,i,p);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::shape(t,o,i,p);

	  default:
	    error();
	  };
      };

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::shape(t,o,i,p);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::shape(t,o,i,p);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::shape(t,o,i,p);

	  default:
	    error();
	  };
      };

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::shape(t,o,i,p);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::shape(t,o,i,p);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::shape(t,o,i,p);

	  default:
	    error();
	  };
      };


    default:
      error();
    };

  
  error();
  return 0.;
};




real FEInterface::shape(const unsigned int dim,
			const FEType& fe_t,
			const Elem* elem,
			const unsigned int i,
			const Point& p)
{
  const Order o = fe_t.order;

  // #ifdef ENABLE_INFINITE_ELEMENTS
  //   if (FEInterface::is_InfFE_elem(elem->type()))
  //     {
  //       return 0.;//InfFE::shape(d, elem, o, i, p);
  //     };
  // #endif

  // by default, use this
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::shape(elem,o,i,p);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::shape(elem,o,i,p);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::shape(elem,o,i,p);

	  default:
	    error();
	  };
      };

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::shape(elem,o,i,p);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::shape(elem,o,i,p);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::shape(elem,o,i,p);

	  default:
	    error();
	  };
      };

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::shape(elem,o,i,p);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::shape(elem,o,i,p);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::shape(elem,o,i,p);

	  default:
	    error();
	  };
      };


    default:
      error();
    };

  
  error();
  return 0.;
};




  
// real FEInterface::shape_deriv(const unsigned int d,
// 			      const ElemType t,
// 			      const Order o,
// 			      const unsigned int i,
// 			      const unsigned int j,
// 			      const Point& p)
// {

// #ifdef ENABLE_INFINITE_ELEMENTS
//   if (FEInterface::is_InfFE_elem(t))
//     {
//       return InfFE::shape_deriv(d, t, o, i, j, p);
//     };
// #endif

//   // by default, use this
//   return FEBase::shape_deriv(d, t, o, i, j, p);
// };





// real FEInterface::shape_deriv(const unsigned int d,
// 			      const Elem* elem,
// 			      const Order o,
// 			      const unsigned int i,
// 			      const unsigned int j,
// 			      const Point& p)
// {

// #ifdef ENABLE_INFINITE_ELEMENTS
//   if (FEInterface::is_InfFE_elem(elem->type()))
//     {
//       return InfFE::shape_deriv(d, elem, o, i, j, p);
//     };
// #endif

//   // by default, use this
//   return FEBase::shape_deriv(d, elem, o, i, j, p);
// };

