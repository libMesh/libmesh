// $Id: fe_interface.C,v 1.7 2003-02-05 20:51:43 ddreyer Exp $

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
  std::cerr << "ERROR: Do not define an object of this type." 
	    << std::endl;
  error();
};




unsigned int FEInterface::n_shape_functions(const unsigned int dim,
					    const FEType& fe_t,
					    const ElemType t)
{

#ifdef ENABLE_INFINITE_ELEMENTS
  /*
   * Since the FEType, stored in DofMap/SystemData, has to
   * be the _same_ for InfFE and FE, we have to catch calls
   * to infinite elements through the element type.
   */

  if ( is_InfFE_elem(t) )
    return ifem_n_shape_functions(dim, fe_t, t);

#endif

  const Order o = fe_t.order;
  
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
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(t) )
    return ifem_n_dofs(dim, fe_t, t);

#endif

  const Order o = fe_t.order;

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
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(t) )
    return ifem_n_dofs_at_node(dim, fe_t, t, n);

#endif

  const Order o = fe_t.order;
  
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
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(t) )
    return ifem_n_dofs_per_elem(dim, fe_t, t);

#endif

  const Order o = fe_t.order;

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
			     const std::vector<Complex>& elem_soln,
			     std::vector<Complex>&       nodal_soln)
{
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
  {
    ifem_nodal_soln(dim, fe_t, elem, elem_soln, nodal_soln);
    return;
  };

#endif

  const Order order = fe_t.order;
  
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
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    return ifem_inverse_map(dim, fe_t, elem, p);

#endif

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
				       const Real eps)
{
  return FEBase::on_reference_element(p,t,eps);
};




Real FEInterface::shape(const unsigned int dim,
			const FEType& fe_t,
			const ElemType t,
			const unsigned int i,
			const Point& p)
{
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(t) )
    return ifem_shape(dim, fe_t, t, i, p);

#endif

  const Order o = fe_t.order;
  
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




Real FEInterface::shape(const unsigned int dim,
			const FEType& fe_t,
			const Elem* elem,
			const unsigned int i,
			const Point& p)
{
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    return ifem_shape(dim, fe_t, elem, i, p);

#endif

  const Order o = fe_t.order;

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


