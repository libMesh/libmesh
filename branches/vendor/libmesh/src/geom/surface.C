// $Id: surface.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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



// C++ includes
#include <iostream>


// Local includes
#include "surface.h"



// ------------------------------------------------------------
// Surface class member functions
Surface::Surface ()
{
};



Surface::Surface (const Surface& other_surface)
{
};



Surface::~Surface ()
{
};



bool Surface::above_surface (const Point& p) const 
{
  std::cerr << "ERROR: This method must be implemented in the derived class!"
	    << std::endl;
		       
  error();

  return false;
};



bool Surface::below_surface (const Point& p) const 
{
  std::cerr << "ERROR: This method must be implemented in the derived class!"
	    << std::endl;
		       
  error();

  return false;
};



bool Surface::on_surface (const Point& p) const 
{
  std::cerr << "ERROR: This method must be implemented in the derived class!"
	    << std::endl;
		       
  error();

  return false;
};



Point Surface::closest_point (const Point& p) const
{
  std::cerr << "ERROR: This method must be implemented in the derived class!"
	    << std::endl;
		       
  error();

  Point cp;

  return cp;
};



Point Surface::unit_normal (const Point& p) const
{
  std::cerr << "ERROR: This method must be implemented in the derived class!"
	    << std::endl;
		       
  error();

  Point n;

  return n;
};
