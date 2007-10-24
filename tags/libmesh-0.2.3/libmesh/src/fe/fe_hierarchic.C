// $Id: fe_hierarchic.C,v 1.6 2003-02-03 03:51:49 ddreyer Exp $

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
#include "fe.h"
#include "elem.h"




// ------------------------------------------------------------
// Hierarchic-specific implementations
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::nodal_soln(const Elem* elem,
			   const Order order,
			   const std::vector<Complex>& elem_soln,
			   std::vector<Complex>&       nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  
  const ElemType type = elem->type();

  nodal_soln.resize(n_nodes);


  
  switch (order)
    {
      // Constant shape functions
    case CONST:
      {
	assert (elem_soln.size() == 1);
	
	const Complex val = elem_soln[0];
	
	for (unsigned int n=0; n<n_nodes; n++)
	  nodal_soln[n] = val;
	
	return;
      };


      // For other bases do interpolation at the nodes
      // explicitly.
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
    case FIFTH:
      {

	const unsigned int n_sf =
	  FE<Dim,T>::n_shape_functions(type, order);
	
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    const Point mapped_point = FE<Dim,T>::inverse_map(elem,
							      elem->point(n));

	    assert (elem_soln.size() == n_sf);

	    // Zero before summation
	    nodal_soln[n] = 0;

	    // u_i = Sum (alpha_i phi_i)
	    for (unsigned int i=0; i<n_sf; i++)
	      nodal_soln[n] += elem_soln[i]*FE<Dim,T>::shape(elem,
							     order,
							     i,
							     mapped_point);	    
	  };

	return;
      };
      
    default:
      {
	error();
	return;
      };
    };

  
  // How did we get here?
  error();  
  return;
};



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {
      // Hierarchic 1st-order polynomials.
    case FIRST:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 2;

	  case TRI6:
	    return 3;

	  case QUAD8:
	  case QUAD9:
	    return 4;

	  case HEX27:
	    return 8;

	  default:
	    error();
	  };
      };


      // Hierarchic 2nd-order polynomials.
    case SECOND:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 3;

	  case TRI6:
	    return 6;

	  case QUAD8:
	  case QUAD9:
	    return 9;
	    
	  case HEX27:
	    return 27;

	  default:
	    error();
	  };
      };


      // Hierarchic 3rd-order polynomials.
    case THIRD:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 4;

	  case TRI6:
	    return 10;
	    
	  case QUAD8:
	  case QUAD9:
	    return 16;

	  case HEX27:
	    return 64;

	  default:
	    error();
	  };
      };


      // Hierarchic 4th-order polynomials.
    case FOURTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 5;

	  case TRI6:
	    return 15;
	    
	  case QUAD8:
	  case QUAD9:
	    return 25;

	  case HEX27:
	    return 125;

	  default:
	    error();
	  };
      };


      // Hierarchic 5th-order polynomials.
    case FIFTH:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 6;

	  case TRI6:
	    return 21;
	      
	  case QUAD8:
	  case QUAD9:
	    return 36;

	  case HEX27:
	    return 216;

	  default:
	    error();
	  };
      };

      
    default:
      {
	error();
      };
    };
  
  error();  
  return 0;
};



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_at_node(const ElemType t,
				       const Order o,
				       const unsigned int n)
{
  switch (o)
    {
      // The first-order hierarchic shape functions
    case FIRST:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 0;
		  
		default:
		  error();		  
		};
	    };

	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 0;

		default:
		  error();
		};
	    };
	    

	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		  return 1;

		case 4:
		case 5:
		case 6:
		case 7:
		  return 0;
		  
		case 8:
		  return 0;
		  
		default:
		  error();
		};
	    };

	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		  return 1;

		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
		case 14:
		case 15:
		case 16:
		case 17:
		case 18:
		case 19:
		  return 0;
		  
		case 20:
		case 21:
		case 22:
		case 23:
		case 24:
		case 25:
		  return 0;
		  
		case 26:
		  return 0;
		};
	    };
	    
	  default:
	    error();
	    
	  };
      };



      // The second-order hierarchic shape functions
    case SECOND:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 1;
		  
		default:
		  error();		  
		};
	    };

	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 1;

		default:
		  error();
		};
	    };


	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		  return 1;

		case 4:
		case 5:
		case 6:
		case 7:
		  return 1;
		  
		case 8:
		  return 1;
		  
		default:
		  error();
		};
	    };

	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		  return 1;

		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
		case 14:
		case 15:
		case 16:
		case 17:
		case 18:
		case 19:
		  return 1;
		  
		case 20:
		case 21:
		case 22:
		case 23:
		case 24:
		case 25:
		  return 1;
		  
		case 26:
		  return 1;
		};
	    };
	    
	  default:
	    error();
	    
	  };
      };



      // The third-order hierarchic shape functions
    case THIRD:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 2;
		  
		default:
		  error();		  
		};
	    };

	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 2;

		default:
		  error();
		};
	    };


	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		  return 1;

		case 4:
		case 5:
		case 6:
		case 7:
		  return 2;
		  
		case 8:
		  return 4;
		  
		default:
		  error();
		};
	    };

	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		  return 1;

		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
		case 14:
		case 15:
		case 16:
		case 17:
		case 18:
		case 19:
		  return 2;
		  
		case 20:
		case 21:
		case 22:
		case 23:
		case 24:
		case 25:
		  return 4;
		  
		case 26:
		  return 8;
		};
	    };
	    
	  default:
	    error();
	    
	  };
      };



      // The fourth-order hierarchic shape functions
    case FOURTH:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 3;
		  
		default:
		  error();		  
		};
	    };

	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 3;

		default:
		  error();
		};
	    };


	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		  return 1;

		case 4:
		case 5:
		case 6:
		case 7:
		  return 3;
		  
		case 8:
		  return 9;
		  
		default:
		  error();
		};
	    };

	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		  return 1;

		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
		case 14:
		case 15:
		case 16:
		case 17:
		case 18:
		case 19:
		  return 3;
		  
		case 20:
		case 21:
		case 22:
		case 23:
		case 24:
		case 25:
		  return 9;
		  
		case 26:
		  return 27;
		};
	    };
	    
	  default:
	    error();
	    
	  };
      };



      // The fifth-order hierarchic shape functions
    case FIFTH:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;

		case 2:
		  return 4;
		  
		default:
		  error();		  
		};
	    };

	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;

		case 3:
		case 4:
		case 5:
		  return 4;

		default:
		  error();
		};
	    };


	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
	  case QUAD9:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		  return 1;

		case 4:
		case 5:
		case 6:
		case 7:
		  return 4;
		  
		case 8:
		  return 16;
		  
		default:
		  error();
		};
	    };

	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		  return 1;

		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
		case 14:
		case 15:
		case 16:
		case 17:
		case 18:
		case 19:
		  return 4;
		  
		case 20:
		case 21:
		case 22:
		case 23:
		case 24:
		case 25:
		  return 16;
		  
		case 26:
		  return 64;
		};
	    };
	    
	  default:
	    error();
	    
	  };
      };


      
    default:
      {
	error();
      };
    };
  
  error();
  
  return 0;
};



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_per_elem(const ElemType t,
					const Order o)
{
  switch (o)
    {
      // The first-order hierarchic shape functions
    case FIRST:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a two-noded edge
	  case EDGE2:
	    return 0;
	    
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    return 0;

	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD8:
	    return 0;
	    
	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	    
	  case QUAD9:
	    return 0;
	    
	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    return 0;
	    
	    
	  default:
	    error();
	    
	  };
      };



      // The second-order hierarchic shape functions
    case SECOND:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a two-noded edge
	  case EDGE2:
	    return 1;
	    
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    return 0;

	    // The 2D tensor-product hierarchics defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 1;

	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;
	    
	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    return 0;


	    
	  default:
	    error();
	    
	  };
      };



      // The third-order hierarchic shape functions
    case THIRD:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a two-noded edge
	  case EDGE2:
	    return 2;
	    
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    return 1;

	    // The 2D tensor-product hierarchics defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 4;

	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;
	    
	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    return 0;


	    
	  default:
	    error();
	    
	  };
      };



      // The fourth-order hierarchic shape functions
    case FOURTH:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a two-noded edge
	  case EDGE2:
	    return 3;
	    
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    return 3;

	    // The 2D tensor-product hierarchics defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 9;

	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;
	    
	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    return 0;


	    
	  default:
	    error();
	    
	  };
      };



      // The fifth-order hierarchic shape functions
    case FIFTH:
      {
	switch (t)
	  {
	    // The 1D heirarchic defined on a two-noded edge
	  case EDGE2:
	    return 4;
	    
	    // The 1D heirarchic defined on a three-noded edge
	  case EDGE3:
	    return 0;
	    
	    // The 2D hierarchic defined on a 6-noded triangle
	  case TRI6:
	    return 6;

	    // The 2D tensor-product hierarchics defined on a
	    // eight-noded quadrilateral.
	  case QUAD8:
	    return 16;

	    // The 2D tensor-product hierarchics defined on a
	    // nine-noded quadrilateral.
	  case QUAD9:
	    return 0;
	    
	    // The 3D tensor-product hierarchics defined on a
	    // twenty-seven noded hexahedral.
	  case HEX27:
	    return 0;


	    
	  default:
	    error();
	    
	  };
      };



      
      // Otherwise no DOFS per element
    default:
      return 0;
    };
};




//--------------------------------------------------------------
// Explicit instantiations
template class FE<1,HIERARCHIC>;
template class FE<2,HIERARCHIC>;
template class FE<3,HIERARCHIC>;
