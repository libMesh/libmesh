// $Id: fe_lagrange.C,v 1.3 2003-01-20 17:06:22 jwpeterson Exp $

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
// Lagrange-specific implementations
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::nodal_soln(const Elem* elem,
			   const Order order,
			   const std::vector<number>& elem_soln,
			   std::vector<number>&       nodal_soln)
{

  const unsigned int n_nodes = elem->n_nodes();
  const ElemType type        = elem->type();

  nodal_soln.resize(n_nodes);


  
  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
	switch (type)
	  {
	  case EDGE3:
	    {
	      assert (elem_soln.size()  == 2);
	      assert (nodal_soln.size() == 3);

	      nodal_soln[0] = elem_soln[0];
	      nodal_soln[1] = elem_soln[1];
	      nodal_soln[2] = .5*(elem_soln[0] + elem_soln[1]);

	      return;
	    };

	    
	  case TRI6:
	    {
	      assert (elem_soln.size()  == 3);
	      assert (nodal_soln.size() == 6);
	      
	      nodal_soln[0] = elem_soln[0];
	      nodal_soln[1] = elem_soln[1];
	      nodal_soln[2] = elem_soln[2];
	      nodal_soln[3] = .5*(elem_soln[0] + elem_soln[1]);	      
	      nodal_soln[4] = .5*(elem_soln[1] + elem_soln[2]);
	      nodal_soln[5] = .5*(elem_soln[2] + elem_soln[0]);

	      return;
	    };
	    

	  case QUAD8:
	  case QUAD9:
	    {
	      assert (elem_soln.size()  == 4);
	      
	      if (type == QUAD8)
		assert (nodal_soln.size() == 8);
	      else
		assert (nodal_soln.size() == 9);

	      
	      nodal_soln[0] = elem_soln[0];
	      nodal_soln[1] = elem_soln[1];
	      nodal_soln[2] = elem_soln[2];
	      nodal_soln[3] = elem_soln[3];
	      nodal_soln[4] = .5*(elem_soln[0] + elem_soln[1]);
	      nodal_soln[5] = .5*(elem_soln[1] + elem_soln[2]);
	      nodal_soln[6] = .5*(elem_soln[2] + elem_soln[3]);
	      nodal_soln[7] = .5*(elem_soln[3] + elem_soln[0]);
					       
	      if (type == QUAD9)
		nodal_soln[8] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3]);

	      return;
	    };

	    
	  case TET10:
	    {
	      assert (elem_soln.size()  == 4);
	      assert (nodal_soln.size() == 10);
	      
	      nodal_soln[0] = elem_soln[0];
	      nodal_soln[1] = elem_soln[1];
	      nodal_soln[2] = elem_soln[2];
	      nodal_soln[3] = elem_soln[3];
	      nodal_soln[4] = .5*(elem_soln[0] + elem_soln[1]);
	      nodal_soln[5] = .5*(elem_soln[1] + elem_soln[2]);
	      nodal_soln[6] = .5*(elem_soln[2] + elem_soln[0]);
	      nodal_soln[7] = .5*(elem_soln[3] + elem_soln[0]);
	      nodal_soln[8] = .5*(elem_soln[3] + elem_soln[1]);
	      nodal_soln[9] = .5*(elem_soln[3] + elem_soln[2]);
	      
	      return;
	    };

	    
	  case HEX20:
	  case HEX27:
	    {
	      assert (elem_soln.size()  == 8);
	      
	      if (type == HEX20)
		assert (nodal_soln.size() == 20);
	      else
		assert (nodal_soln.size() == 27);
	      
	      nodal_soln[0]  = elem_soln[0];
	      nodal_soln[1]  = elem_soln[1];
	      nodal_soln[2]  = elem_soln[2];
	      nodal_soln[3]  = elem_soln[3];
	      nodal_soln[4]  = elem_soln[4];
	      nodal_soln[5]  = elem_soln[5];
	      nodal_soln[6]  = elem_soln[6];
	      nodal_soln[7]  = elem_soln[7];
	      nodal_soln[8]  = .5*(elem_soln[0] + elem_soln[1]);
	      nodal_soln[9]  = .5*(elem_soln[1] + elem_soln[2]);
	      nodal_soln[10] = .5*(elem_soln[2] + elem_soln[3]);
	      nodal_soln[11] = .5*(elem_soln[3] + elem_soln[0]);
	      nodal_soln[12] = .5*(elem_soln[0] + elem_soln[4]);
	      nodal_soln[13] = .5*(elem_soln[1] + elem_soln[5]);
	      nodal_soln[14] = .5*(elem_soln[2] + elem_soln[6]);
	      nodal_soln[15] = .5*(elem_soln[3] + elem_soln[7]);
	      nodal_soln[16] = .5*(elem_soln[4] + elem_soln[5]);
	      nodal_soln[17] = .5*(elem_soln[5] + elem_soln[6]);
	      nodal_soln[18] = .5*(elem_soln[6] + elem_soln[7]);
	      nodal_soln[19] = .5*(elem_soln[4] + elem_soln[7]);

	      if (type == HEX27)
		{
		  nodal_soln[20] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3]);
		  nodal_soln[21] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[4] + elem_soln[5]);
		  nodal_soln[22] = .25*(elem_soln[1] + elem_soln[2] + elem_soln[5] + elem_soln[6]);
		  nodal_soln[23] = .25*(elem_soln[2] + elem_soln[3] + elem_soln[6] + elem_soln[7]);
		  nodal_soln[24] = .25*(elem_soln[3] + elem_soln[0] + elem_soln[7] + elem_soln[4]);
		  nodal_soln[25] = .25*(elem_soln[4] + elem_soln[5] + elem_soln[6] + elem_soln[7]);

		  nodal_soln[26] = .125*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3] +
					 elem_soln[4] + elem_soln[5] + elem_soln[6] + elem_soln[7]);
		};
	      
	      return;
	    };
	    
	    
	  case PRISM18:
	    {
	      error();
	      
	      return;
	    };


	    
	  default:
	    {
	      // By default the element solution _is_ nodal,
	      // so just copy it.
	      nodal_soln = elem_soln;
	      
	      return;
	    };
	  };
      };


      
    default:
      {
	// By default the element solution _is_ nodal,
	// so just copy it.
	nodal_soln = elem_soln;
	
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

      // linear Lagrange shape functions
    case FIRST:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    return 2;

	  case TRI3:
	  case TRI6:
	    return 3;

	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    return 4;

	  case TET4:
	  case TET10:
	    return 4;

	  case HEX8:
	  case HEX20:
	  case HEX27:
	    return 8;

	  case PRISM6:
	  case PRISM18:
	    return 6;

	  case PYRAMID5:
	    return 5;

	  default:
	    error();	    
	  };
      };

      
      // quadratic Lagrange shape functions
    case SECOND:
      {
	switch (t)
	  {
	  case EDGE3:
	    return 3;

	  case TRI6:
	    return 6;

	  case QUAD8:
	    return 8;
	    
	  case QUAD9:
	    return 9;

	  case TET10:
	    return 10;

	  case HEX20:
	    return 20;
	    
	  case HEX27:
	    return 27;

	  case PRISM18:
	    return 18;

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



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_at_node(const ElemType t,
				       const Order o,
				       const unsigned int n)
{
  switch (o)
    {
      
      // linear Lagrange shape functions
    case FIRST:
      {
	switch (t)
	  {
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		  return 1;
		  
		default:
		  return 0;
		};
	    };

	  case TRI3:
	  case TRI6:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		  return 1;
		  
		default:
		  return 0;
		};
	    };

	  case QUAD4:
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
		  
		default:
		  return 0;
		};
	    };


	  case TET4:
	  case TET10:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		  return 1;
		  
		default:
		  return 0;
		};
	    };
		    
	  case HEX8:
	  case HEX20:
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
		  
		default:
		  return 0;
		};
	    };

	  case PRISM6:
	  case PRISM18:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		  return 1;
		  
		default:
		  return 0;
		};
	    };

	  case PYRAMID5:
	    {
	      switch (n)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
		  return 1;
		  
		default:
		  return 0;
		};
	    };

	  default:
	    error();	    
	  };
      };

      // quadratic Lagrange shape functions
    case SECOND:
      {
	switch (t)
	  {
	    // quadratic lagrange has one dof at each node
	  case EDGE3:
	  case TRI6:
	  case QUAD8:
	  case QUAD9:
	  case TET10:
	  case HEX20:
	  case HEX27:
	  case PRISM18:
	    return 1;

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



template <unsigned int Dim, FEFamily T>
unsigned int FE<Dim,T>::n_dofs_per_elem(const ElemType,
					const Order)
{
  // Lagrange elements have no dofs per element
  // (just at the nodes)
  
  return 0;
};




//--------------------------------------------------------------
// Explicit instantiations
template class FE<1,LAGRANGE>;
template class FE<2,LAGRANGE>;
template class FE<3,LAGRANGE>;
