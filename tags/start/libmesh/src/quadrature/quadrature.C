// $Id: quadrature.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include "quadrature.h"
#include "point.h"



void QGauss::init(const ElemType t)
{
  // check to see if we have already
  // done the work for this quadrature rule
  if (t == _type)
    return;
  else
    _type = t;
    
  
  
  switch(_dim)
    {
    case 1:
      init_1D(_type);

      return;
      
    case 2:
      init_2D(_type);

      return;

    case 3:
      init_3D(_type);

      return;

    default:
      error();
    };

  error();
  
  return;
};




void QGauss::init(const ElemType t,
		  const unsigned int s)
{
  static unsigned int side = static_cast<unsigned int>(-1);

  
  // check to see if we have already
  // done the work for this quadrature rule
  if ((t == _type) &&
      (s == side)) 
    {
      //here();
      return;
    }
  else
    {
      _type = t;
      side = s;
    }


  switch (_dim)
    {
    case 2:
      init_2D(_type,side);

      return;

    case 3:
      init_3D(_type,side);

      return;


    default:
      error();

    };

  error();

  return;
};



void QTrap::init(const ElemType t)
{
  // check to see if we have already
  // done the work for this quadrature rule
  if (t == _type)
    return;
  else
    _type = t;
    
  
  
  switch(_dim)
    {
    case 1:
      init_1D(_type);

      return;
      
    case 2:
      init_2D(_type);

      return;

    case 3:
      init_3D(_type);

      return;

    default:
      error();
    };

  error();
  
  return;
};




void QTrap::init(const ElemType t,
		 const unsigned int s)
{
  static unsigned int side = static_cast<unsigned int>(-1);

  
  // check to see if we have already
  // done the work for this quadrature rule
  if ((t == _type) &&
      (s == side)) 
    {
      //here();
      return;
    }
  else
    {
      _type = t;
      side = s;
    }


  switch (_dim)
    {
    case 2:
      init_2D(_type,side);

      return;

    case 3:
      init_3D(_type,side);

      return;


    default:
      error();

    };

  error();

  return;
};
