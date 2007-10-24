// $Id: fe_interface.C,v 1.47 2006-08-30 18:06:29 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "fe_compute_data.h"
#include "dof_map.h"

//------------------------------------------------------------
//FEInterface class members
FEInterface::FEInterface()
{
  std::cerr << "ERROR: Do not define an object of this type." 
	    << std::endl;
  error();
}




unsigned int FEInterface::n_shape_functions(const unsigned int dim,
					    const FEType& fe_t,
					    const ElemType t)
{

#ifdef ENABLE_INFINITE_ELEMENTS
  /*
   * Since the FEType, stored in DofMap/(some System child), has to
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
	  case CLOUGH:
	    return FE<1,CLOUGH>::n_shape_functions(t, o);
	    
	  case HERMITE:
	    return FE<1,HERMITE>::n_shape_functions(t, o);
	    
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::n_shape_functions(t, o);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::n_shape_functions(t, o);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::n_shape_functions(t, o);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<1,BERNSTEIN>::n_shape_functions(t, o);

	  case SZABAB:
	    return FE<1,SZABAB>::n_shape_functions(t, o);

#endif
	    
	  case XYZ:
	    return FEXYZ<1>::n_shape_functions(t, o);


	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    return FE<2,CLOUGH>::n_shape_functions(t, o);
	    
	  case HERMITE:
	    return FE<2,HERMITE>::n_shape_functions(t, o);
	    
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::n_shape_functions(t, o);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::n_shape_functions(t, o);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::n_shape_functions(t, o);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<2,BERNSTEIN>::n_shape_functions(t, o);

	  case SZABAB:
	    return FE<2,SZABAB>::n_shape_functions(t, o);

#endif
	    
	  case XYZ:
	    return FEXYZ<2>::n_shape_functions(t, o);

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    return FE<3,HERMITE>::n_shape_functions(t, o);
	    
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::n_shape_functions(t, o);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::n_shape_functions(t, o);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::n_shape_functions(t, o);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<3,BERNSTEIN>::n_shape_functions(t, o);

	  case SZABAB:
	    return FE<3,SZABAB>::n_shape_functions(t, o);

#endif
	    
	  case XYZ:
	    return FEXYZ<3>::n_shape_functions(t, o);

	  default:
	    error();
	  }
      }


    default:
      error();
    }

  
  error();
  return 0;
}





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
	  case CLOUGH:
	    return FE<1,CLOUGH>::n_dofs(t, o);
	    
	  case HERMITE:
	    return FE<1,HERMITE>::n_dofs(t, o);
	    
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::n_dofs(t, o);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::n_dofs(t, o);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::n_dofs(t, o);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<1,BERNSTEIN>::n_dofs(t, o);

	  case SZABAB:
	    return FE<1,SZABAB>::n_dofs(t, o);

#endif
	    
	  case XYZ:
	    return FEXYZ<1>::n_dofs(t, o);

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    return FE<2,CLOUGH>::n_dofs(t, o);
	    
	  case HERMITE:
	    return FE<2,HERMITE>::n_dofs(t, o);
	    
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::n_dofs(t, o);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::n_dofs(t, o);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::n_dofs(t, o);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<2,BERNSTEIN>::n_dofs(t, o);

	  case SZABAB:
	    return FE<2,SZABAB>::n_dofs(t, o);

#endif
	    
	  case XYZ:
	    return FEXYZ<2>::n_dofs(t, o);

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    return FE<3,HERMITE>::n_dofs(t, o);
	    
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::n_dofs(t, o);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::n_dofs(t, o);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::n_dofs(t, o);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<3,BERNSTEIN>::n_dofs(t, o);

	  case SZABAB:
	    return FE<3,SZABAB>::n_dofs(t, o);

#endif
	    
	  case XYZ:
	    return FEXYZ<3>::n_dofs(t, o);

	  default:
	    error();
	  }
      }


    default:
      error();
    }

  
  error();
  return 0;
}

		


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
	  case CLOUGH:
	    return FE<1,CLOUGH>::n_dofs_at_node(t, o, n);
	    
	  case HERMITE:
	    return FE<1,HERMITE>::n_dofs_at_node(t, o, n);
	    
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::n_dofs_at_node(t, o, n);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::n_dofs_at_node(t, o, n);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::n_dofs_at_node(t, o, n);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<1,BERNSTEIN>::n_dofs_at_node(t, o, n);

	  case SZABAB:
	    return FE<1,SZABAB>::n_dofs_at_node(t, o, n);

#endif
	    
	  case XYZ:
	    return FEXYZ<1>::n_dofs_at_node(t, o, n);

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    return FE<2,CLOUGH>::n_dofs_at_node(t, o, n);
	    
	  case HERMITE:
	    return FE<2,HERMITE>::n_dofs_at_node(t, o, n);
	    
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::n_dofs_at_node(t, o, n);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::n_dofs_at_node(t, o, n);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::n_dofs_at_node(t, o, n);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<2,BERNSTEIN>::n_dofs_at_node(t, o, n);

	  case SZABAB:
	    return FE<2,SZABAB>::n_dofs_at_node(t, o, n);

#endif
	    
	  case XYZ:
	    return FEXYZ<2>::n_dofs_at_node(t, o, n);

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    return FE<3,HERMITE>::n_dofs_at_node(t, o, n);
	    
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::n_dofs_at_node(t, o, n);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::n_dofs_at_node(t, o, n);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::n_dofs_at_node(t, o, n);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<3,BERNSTEIN>::n_dofs_at_node(t, o, n);

	  case SZABAB:
	    return FE<3,SZABAB>::n_dofs_at_node(t, o, n);

#endif
	    
	  case XYZ:
	    return FEXYZ<3>::n_dofs_at_node(t, o, n);

	  default:
	    error();
	  }
      }


    default:
      error();
    }

  
  error();
  return 0;
}





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
	  case CLOUGH:
	    return FE<1,CLOUGH>::n_dofs_per_elem(t, o);
	    
	  case HERMITE:
	    return FE<1,HERMITE>::n_dofs_per_elem(t, o);
	    
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::n_dofs_per_elem(t, o);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::n_dofs_per_elem(t, o);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::n_dofs_per_elem(t, o);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<1,BERNSTEIN>::n_dofs_per_elem(t, o);

	  case SZABAB:
	    return FE<1,SZABAB>::n_dofs_per_elem(t, o);

#endif
	    
	  case XYZ:
	    return FEXYZ<1>::n_dofs_per_elem(t, o);

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    return FE<2,CLOUGH>::n_dofs_per_elem(t, o);
	    
	  case HERMITE:
	    return FE<2,HERMITE>::n_dofs_per_elem(t, o);
	    
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::n_dofs_per_elem(t, o);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::n_dofs_per_elem(t, o);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::n_dofs_per_elem(t, o);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<2,BERNSTEIN>::n_dofs_per_elem(t, o);

	  case SZABAB:
	    return FE<2,SZABAB>::n_dofs_per_elem(t, o);

#endif
	    
	  case XYZ:
	    return FEXYZ<2>::n_dofs_per_elem(t, o);

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    return FE<3,HERMITE>::n_dofs_per_elem(t, o);
	    
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::n_dofs_per_elem(t, o);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::n_dofs_per_elem(t, o);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::n_dofs_per_elem(t, o);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<3,BERNSTEIN>::n_dofs_per_elem(t, o);

	  case SZABAB:
	    return FE<3,SZABAB>::n_dofs_per_elem(t, o);

#endif
	    
	  case XYZ:
	    return FEXYZ<3>::n_dofs_per_elem(t, o);

	  default:
	    error();
	  }
      }


    default:
      error();
    }

  
  error();
  return 0;
}




void FEInterface::dofs_on_side(const Elem* const elem,
			       const unsigned int dim,
			       const FEType& fe_t,
			       unsigned int s,
			       std::vector<unsigned int>& di)
{
  const Order o = fe_t.order;

  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<1,CLOUGH>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case HERMITE:
	    FE<1,HERMITE>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case HIERARCHIC:
	    FE<1,HIERARCHIC>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case LAGRANGE:
	    FE<1,LAGRANGE>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case MONOMIAL:
	    FE<1,MONOMIAL>::dofs_on_side(elem, o, s, di);
            return;

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<1,BERNSTEIN>::dofs_on_side(elem, o, s, di);
            return;

	  case SZABAB:
	    FE<1,SZABAB>::dofs_on_side(elem, o, s, di);
            return;

#endif
	    
	  case XYZ:
	    FEXYZ<1>::dofs_on_side(elem, o, s, di);
            return;

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<2,CLOUGH>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case HERMITE:
	    FE<2,HERMITE>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case HIERARCHIC:
	    FE<2,HIERARCHIC>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case LAGRANGE:
	    FE<2,LAGRANGE>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case MONOMIAL:
	    FE<2,MONOMIAL>::dofs_on_side(elem, o, s, di);
            return;

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<2,BERNSTEIN>::dofs_on_side(elem, o, s, di);
            return;

	  case SZABAB:
	    FE<2,SZABAB>::dofs_on_side(elem, o, s, di);
            return;

#endif
	    
	  case XYZ:
	    FEXYZ<2>::dofs_on_side(elem, o, s, di);
            return;

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    FE<3,HERMITE>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case HIERARCHIC:
	    FE<3,HIERARCHIC>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case LAGRANGE:
	    FE<3,LAGRANGE>::dofs_on_side(elem, o, s, di);
            return;
	    
	  case MONOMIAL:
	    FE<3,MONOMIAL>::dofs_on_side(elem, o, s, di);
            return;

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<3,BERNSTEIN>::dofs_on_side(elem, o, s, di);
            return;

	  case SZABAB:
	    FE<3,SZABAB>::dofs_on_side(elem, o, s, di);
            return;

#endif
	    
	  case XYZ:
	    FEXYZ<3>::dofs_on_side(elem, o, s, di);
            return;

	  default:
	    error();
	  }
      }


    default:
      error();
    }

  error();
}



void FEInterface::dofs_on_edge(const Elem* const elem,
			       const unsigned int dim,
			       const FEType& fe_t,
			       unsigned int e,
			       std::vector<unsigned int>& di)
{
  const Order o = fe_t.order;

  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<1,CLOUGH>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case HERMITE:
	    FE<1,HERMITE>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case HIERARCHIC:
	    FE<1,HIERARCHIC>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case LAGRANGE:
	    FE<1,LAGRANGE>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case MONOMIAL:
	    FE<1,MONOMIAL>::dofs_on_edge(elem, o, e, di);
            return;

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<1,BERNSTEIN>::dofs_on_edge(elem, o, e, di);
            return;

	  case SZABAB:
	    FE<1,SZABAB>::dofs_on_edge(elem, o, e, di);
            return;

#endif
	    
	  case XYZ:
	    FEXYZ<1>::dofs_on_edge(elem, o, e, di);
            return;

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<2,CLOUGH>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case HERMITE:
	    FE<2,HERMITE>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case HIERARCHIC:
	    FE<2,HIERARCHIC>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case LAGRANGE:
	    FE<2,LAGRANGE>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case MONOMIAL:
	    FE<2,MONOMIAL>::dofs_on_edge(elem, o, e, di);
            return;

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<2,BERNSTEIN>::dofs_on_edge(elem, o, e, di);
            return;

	  case SZABAB:
	    FE<2,SZABAB>::dofs_on_edge(elem, o, e, di);
            return;

#endif
	    
	  case XYZ:
	    FEXYZ<2>::dofs_on_edge(elem, o, e, di);
            return;

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    FE<3,HERMITE>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case HIERARCHIC:
	    FE<3,HIERARCHIC>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case LAGRANGE:
	    FE<3,LAGRANGE>::dofs_on_edge(elem, o, e, di);
            return;
	    
	  case MONOMIAL:
	    FE<3,MONOMIAL>::dofs_on_edge(elem, o, e, di);
            return;

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<3,BERNSTEIN>::dofs_on_edge(elem, o, e, di);
            return;

	  case SZABAB:
	    FE<3,SZABAB>::dofs_on_edge(elem, o, e, di);
            return;

#endif
	    
	  case XYZ:
	    FEXYZ<3>::dofs_on_edge(elem, o, e, di);
            return;

	  default:
	    error();
	  }
      }


    default:
      error();
    }

  error();
}




void FEInterface::nodal_soln(const unsigned int dim,
			     const FEType& fe_t,
			     const Elem* elem,
			     const std::vector<Number>& elem_soln,
			     std::vector<Number>&       nodal_soln)
{
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
  {
    ifem_nodal_soln(dim, fe_t, elem, elem_soln, nodal_soln);
    return;
  }

#endif

  const Order order = fe_t.order;
  
  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<1,CLOUGH>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;
	    
	  case HERMITE:
	    FE<1,HERMITE>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;
	    
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

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<1,BERNSTEIN>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;

	  case SZABAB:
	    FE<1,SZABAB>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;

#endif

	  case XYZ:
	    FEXYZ<1>::nodal_soln(elem, order,
				 elem_soln, nodal_soln);
	    return;

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<2,CLOUGH>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;
	    
	  case HERMITE:
	    FE<2,HERMITE>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;
	    
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

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<2,BERNSTEIN>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;

	  case SZABAB:
	    FE<2,SZABAB>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;

#endif

	  case XYZ:
	    FEXYZ<2>::nodal_soln(elem, order,
				 elem_soln, nodal_soln);
	    return;

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    FE<3,HERMITE>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;
	    
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

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<3,BERNSTEIN>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;

	  case SZABAB:
	    FE<3,SZABAB>::nodal_soln(elem, order,
				     elem_soln, nodal_soln);
	    return;

#endif

	  case XYZ:
	    FEXYZ<3>::nodal_soln(elem, order,
				 elem_soln, nodal_soln);
	    return;

	  default:
	    error();
	  }
      }


    default:
      error();
    }
}




Point FEInterface::inverse_map (const unsigned int dim,
				const FEType& fe_t,
				const Elem* elem,
				const Point& p,
				const Real tolerance,
				const bool secure)
{
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    return ifem_inverse_map(dim, fe_t, elem, p,tolerance, secure);

#endif

  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    return FE<1,CLOUGH>::inverse_map(elem, p, tolerance, secure);
	    
	  case HERMITE:
	    return FE<1,HERMITE>::inverse_map(elem, p, tolerance, secure);
	    
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::inverse_map(elem, p, tolerance, secure);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::inverse_map(elem, p, tolerance, secure);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::inverse_map(elem, p, tolerance, secure);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<1,BERNSTEIN>::inverse_map(elem, p, tolerance, secure);

	  case SZABAB:
	    return FE<1,SZABAB>::inverse_map(elem, p, tolerance, secure);

#endif
	  case XYZ:
	    return FEXYZ<1>::inverse_map(elem, p, tolerance, secure);


	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    return FE<2,CLOUGH>::inverse_map(elem, p, tolerance, secure);
	    
	  case HERMITE:
	    return FE<2,HERMITE>::inverse_map(elem, p, tolerance, secure);
	    
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::inverse_map(elem, p, tolerance, secure);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::inverse_map(elem, p, tolerance, secure);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::inverse_map(elem, p, tolerance, secure);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<2,BERNSTEIN>::inverse_map(elem, p, tolerance, secure);

	  case SZABAB:
	    return FE<2,SZABAB>::inverse_map(elem, p, tolerance, secure);

#endif
	  case XYZ:
	    return FEXYZ<2>::inverse_map(elem, p, tolerance, secure);

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    return FE<3,HERMITE>::inverse_map(elem, p, tolerance, secure);
	    
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::inverse_map(elem, p, tolerance, secure);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::inverse_map(elem, p, tolerance, secure);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::inverse_map(elem, p, tolerance, secure);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<3,BERNSTEIN>::inverse_map(elem, p, tolerance, secure);

	  case SZABAB:
	    return FE<3,SZABAB>::inverse_map(elem, p, tolerance, secure);

#endif
	  case XYZ:
	    return FEXYZ<3>::inverse_map(elem, p, tolerance, secure);

	  default:
	    error();
	  }
      }


    default:
      error();
    }

  
  error();
  Point pt;
  return pt;
}




void FEInterface::inverse_map (const unsigned int dim,
			       const FEType& fe_t,
			       const Elem* elem,
			       const std::vector<Point>& physical_points,
			       std::vector<Point>&       reference_points,
			       const Real tolerance,
			       const bool secure)
{
  const unsigned int n_pts = physical_points.size();

  // Resize the vector
  reference_points.resize(n_pts);
  
  if (n_pts == 0)
    {
      std::cerr << "WARNING: empty vector physical_points!"
		<< std::endl;
      here();
      return;
    }
  



  
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    {
      ifem_inverse_map(dim, fe_t, elem, physical_points, reference_points, tolerance, secure);
      return;

//       std::cerr << "ERROR: Not implemented!"
// 		<< std::endl;
//       error();
    }
  
#endif

  switch (dim)
    {
      // 1D
    case 1:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<1,CLOUGH>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case HERMITE:
	    FE<1,HERMITE>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case HIERARCHIC:
	    FE<1,HIERARCHIC>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case LAGRANGE:
	    FE<1,LAGRANGE>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case MONOMIAL:
	    FE<1,MONOMIAL>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<1,BERNSTEIN>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;

	  case SZABAB:
	    FE<1,SZABAB>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;

#endif
	    
	  case XYZ:
	    FEXYZ<1>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<2,CLOUGH>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case HERMITE:
	    FE<2,HERMITE>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case HIERARCHIC:
	    FE<2,HIERARCHIC>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case LAGRANGE:
	    FE<2,LAGRANGE>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case MONOMIAL:
	    FE<2,MONOMIAL>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<2,BERNSTEIN>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;

	  case SZABAB:
	    FE<2,SZABAB>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;

#endif
	    
	  case XYZ:
	    FEXYZ<2>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    FE<3,HERMITE>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case HIERARCHIC:
	    FE<3,HIERARCHIC>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case LAGRANGE:
	    FE<3,LAGRANGE>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  case MONOMIAL:
	    FE<3,MONOMIAL>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    FE<3,BERNSTEIN>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;

	  case SZABAB:
	    FE<3,SZABAB>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;

#endif
	    
	  case XYZ:
	    FEXYZ<3>::inverse_map(elem, physical_points, reference_points, tolerance, secure);
	    return;
	    
	  default:
	    error();
	  }
      }


    default:
      error();
    }

  
  error();
  return;
}



bool FEInterface::on_reference_element(const Point& p,
				       const ElemType t,
				       const Real eps)
{
  return FEBase::on_reference_element(p,t,eps);
}




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
	  case CLOUGH:
	    return FE<1,CLOUGH>::shape(t,o,i,p);
	    
	  case HERMITE:
	    return FE<1,HERMITE>::shape(t,o,i,p);
	    
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::shape(t,o,i,p);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::shape(t,o,i,p);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::shape(t,o,i,p);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<1,BERNSTEIN>::shape(t,o,i,p);

	  case SZABAB:
	    return FE<1,SZABAB>::shape(t,o,i,p);

#endif
	    
	  case XYZ:
	    return FEXYZ<1>::shape(t,o,i,p);

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    return FE<2,CLOUGH>::shape(t,o,i,p);
	    
	  case HERMITE:
	    return FE<2,HERMITE>::shape(t,o,i,p);
	    
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::shape(t,o,i,p);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::shape(t,o,i,p);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::shape(t,o,i,p);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<2,BERNSTEIN>::shape(t,o,i,p);

	  case SZABAB:
	    return FE<2,SZABAB>::shape(t,o,i,p);

#endif
	    
	  case XYZ:
	    return FEXYZ<2>::shape(t,o,i,p);

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    return FE<3,HERMITE>::shape(t,o,i,p);
	    
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::shape(t,o,i,p);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::shape(t,o,i,p);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::shape(t,o,i,p);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<3,BERNSTEIN>::shape(t,o,i,p);

	  case SZABAB:
	    return FE<3,SZABAB>::shape(t,o,i,p);

#endif
	    
	  case XYZ:
	    return FEXYZ<3>::shape(t,o,i,p);

	  default:
	    error();
	  }
      }


    default:
      error();
    }

  
  error();
  return 0.;
}




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
	  case CLOUGH:
	    return FE<1,CLOUGH>::shape(elem,o,i,p);
	    
	  case HERMITE:
	    return FE<1,HERMITE>::shape(elem,o,i,p);
	    
	  case HIERARCHIC:
	    return FE<1,HIERARCHIC>::shape(elem,o,i,p);
	    
	  case LAGRANGE:
	    return FE<1,LAGRANGE>::shape(elem,o,i,p);
	    
	  case MONOMIAL:
	    return FE<1,MONOMIAL>::shape(elem,o,i,p);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<1,BERNSTEIN>::shape(elem,o,i,p);

	  case SZABAB:
	    return FE<1,SZABAB>::shape(elem,o,i,p);

#endif
	    
	  case XYZ:
	    return FEXYZ<1>::shape(elem,o,i,p);

	  default:
	    error();
	  }
      }

      
      // 2D
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    return FE<2,CLOUGH>::shape(elem,o,i,p);
	    
	  case HERMITE:
	    return FE<2,HERMITE>::shape(elem,o,i,p);
	    
	  case HIERARCHIC:
	    return FE<2,HIERARCHIC>::shape(elem,o,i,p);
	    
	  case LAGRANGE:
	    return FE<2,LAGRANGE>::shape(elem,o,i,p);
	    
	  case MONOMIAL:
	    return FE<2,MONOMIAL>::shape(elem,o,i,p);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<2,BERNSTEIN>::shape(elem,o,i,p);

	  case SZABAB:
	    return FE<2,SZABAB>::shape(elem,o,i,p);

#endif
	    
	  case XYZ:
	    return FEXYZ<2>::shape(elem,o,i,p);

	  default:
	    error();
	  }
      }

      
      // 3D
    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    return FE<3,HERMITE>::shape(elem,o,i,p);
	    
	  case HIERARCHIC:
	    return FE<3,HIERARCHIC>::shape(elem,o,i,p);
	    
	  case LAGRANGE:
	    return FE<3,LAGRANGE>::shape(elem,o,i,p);
	    
	  case MONOMIAL:
	    return FE<3,MONOMIAL>::shape(elem,o,i,p);

#ifdef ENABLE_HIGHER_ORDER_SHAPES

	  case BERNSTEIN:
	    return FE<3,BERNSTEIN>::shape(elem,o,i,p);

	  case SZABAB:
	    return FE<3,SZABAB>::shape(elem,o,i,p);

#endif
	    
	  case XYZ:
	    return FEXYZ<3>::shape(elem,o,i,p);

	  default:
	    error();
	  }
      }


    default:
      error();
    }

  
  error();
  return 0.;
}




void FEInterface::compute_data(const unsigned int dim,
			       const FEType& fe_t,
			       const Elem* elem,
			       FEComputeData& data)
{
#ifdef ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    {
      data.init();
      ifem_compute_data(dim, fe_t, elem, data);
      return;
    }

#endif

  FEType p_refined = fe_t;
  p_refined.order = static_cast<Order>(p_refined.order + elem->p_level());

  const unsigned int n_dof = n_dofs (dim, p_refined, elem->type());
  const Point&       p     = data.p;
  data.shape.resize(n_dof);

  // set default values for all the output fields
  data.init();

  for (unsigned int n=0; n<n_dof; n++)
      data.shape[n] = shape(dim, p_refined, elem, n, p);

   return;
}




void FEInterface::compute_constraints (DofConstraints &constraints,
				       DofMap &dof_map,
				       const unsigned int variable_number,
				       const Elem* elem)
{
  assert (elem != NULL);

  const FEType& fe_t = dof_map.variable_type(variable_number);
  
  switch (elem->dim())
    {
    case 1:
      {
	// No constraints in 1D.
	return;
      }

      
    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<2,CLOUGH>::compute_constraints (constraints,
					       dof_map,
					       variable_number,
					       elem); return;

	  case HERMITE:
	    FE<2,HERMITE>::compute_constraints (constraints,
					        dof_map,
					        variable_number,
					        elem); return;

	  case LAGRANGE:
	    FE<2,LAGRANGE>::compute_constraints (constraints,
						 dof_map,
						 variable_number,
						 elem); return;

	  case HIERARCHIC:
	    FE<2,HIERARCHIC>::compute_constraints (constraints,
						   dof_map,
						   variable_number,
						   elem); return;

	  default:
	    return;
	  }
      }


    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    FE<3,HERMITE>::compute_constraints (constraints,
					        dof_map,
					        variable_number,
					        elem); return;

	  case LAGRANGE:
	    FE<3,LAGRANGE>::compute_constraints (constraints,
					         dof_map,
						 variable_number,
						 elem); return;      

	  case HIERARCHIC:
	    FE<3,HIERARCHIC>::compute_constraints (constraints,
						   dof_map,
						   variable_number,
						   elem); return;
	  default:
	    return;
	  }
      }

      
    default:
      error();
    }
}
  


unsigned int FEInterface::max_order(const FEType& fe_t,
			            const ElemType& el_t)
{
  // Yeah, I know, infinity is much larger than 11, but our
  // solvers don't seem to like high degree polynomials, and our
  // quadrature rules and number_lookups tables
  // need to go up higher.
  const unsigned int unlimited = 11;

  // If we used 0 as a default, then elements missing from this 
  // table (e.g. infinite elements) would be considered broken.
  const unsigned int unknown = unlimited;

  switch (fe_t.family)
    {
      case LAGRANGE:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	      return 3;
	    case TRI3:
	      return 1;
	    case TRI6:
	      return 2;
	    case QUAD4:
	      return 1;
	    case QUAD8:
	    case QUAD9:
	      return 2;
	    case TET4:
	      return 1;
	    case TET10:
	      return 2;
	    case HEX8:
	      return 1;
	    case HEX20:
	    case HEX27:
	      return 2;
	    case PRISM6:
	    case PRISM15:
	      return 1;
	    case PRISM18:
	      return 2;
	    case PYRAMID5:
	      return 1;
	    default:
	      return unknown;
	  }
	break;
      case MONOMIAL:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	    case TRI3:
	    case TRI6:
	    case QUAD4:
	    case QUAD8:
	    case QUAD9:
	    case TET4:
	    case TET10:
	    case HEX8:
	    case HEX20:
	    case HEX27:
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 4;
	    default:
	      return unknown;
	  }
	break;
#ifdef ENABLE_HIGHER_ORDER_SHAPES
      case BERNSTEIN:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	      return unlimited;
	    case TRI3:
	      return 0;
	    case TRI6:
	      return 6;
	    case QUAD4:
	      return 0;
	    case QUAD8:
	    case QUAD9:
	      return unlimited;
	    case TET4:
	      return 1;
	    case TET10:
	      return 2;
	    case HEX8:
	      return 0;
	    case HEX20:
	      return 2;
	    case HEX27:
	      return 4;
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
      case SZABAB:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	      return 7;
	    case TRI3:
	      return 0;
	    case TRI6:
	      return 7;
	    case QUAD4:
	      return 0;
	    case QUAD8:
	    case QUAD9:
	      return 7;
	    case TET4:
	    case TET10:
	    case HEX8:
	    case HEX20:
	    case HEX27:
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
#endif
      case XYZ:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	    case TRI3:
	    case TRI6:
	    case QUAD4:
	    case QUAD8:
	    case QUAD9:
	    case TET4:
	    case TET10:
	    case HEX8:
	    case HEX20:
	    case HEX27:
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 4;
	    default:
	      return unknown;
	  }
	break;
      case CLOUGH:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	      return 3;
	    case EDGE4:
	    case TRI3:
	      return 0;
	    case TRI6:
	      return 3;
	    case QUAD4:
	    case QUAD8:
	    case QUAD9:
	    case TET4:
	    case TET10:
	    case HEX8:
	    case HEX20:
	    case HEX27:
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
      case HERMITE:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	      return unlimited;
	    case EDGE4:
	    case TRI3:
	    case TRI6:
	      return 0;
	    case QUAD4:
	      return 3;
	    case QUAD8:
	    case QUAD9:
	      return unlimited;
	    case TET4:
	    case TET10:
	      return 0;
	    case HEX8:
	      return 3;
	    case HEX20:
	    case HEX27:
	      return unlimited;
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
      case HIERARCHIC:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	      return unlimited;
	    case TRI3:
	      return 1;
	    case TRI6:
	      return unlimited;
	    case QUAD4:
	      return 1;
	    case QUAD8:
	    case QUAD9:
	      return unlimited;
	    case TET4:
	    case TET10:
	      return 0;
	    case HEX8:
	    case HEX20:
	      return 1;
	    case HEX27:
	      return unlimited;
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
      default:
	return 0;
	break;
    }
}

  

bool FEInterface::extra_hanging_dofs(const FEType& fe_t)
{
  switch (fe_t.family)
    {
      case LAGRANGE:
      case MONOMIAL:
#ifdef ENABLE_HIGHER_ORDER_SHAPES
      case BERNSTEIN:
      case SZABAB:
#endif
      case XYZ:
	return false;
      case CLOUGH:
      case HERMITE:
      case HIERARCHIC:
      default:
	return true;
    }
}
