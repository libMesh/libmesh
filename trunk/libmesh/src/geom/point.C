// $Id: point.C,v 1.5 2003-01-24 17:24:43 jwpeterson Exp $

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
#include "point.h"
#include <iomanip>




// ------------------------------------------------------------
// Point class member funcions
Point Point::cross(const Point& p) const
{
  assert(DIM == 3);

  // |     i          j          k    |
  // |(*this)(0) (*this)(1) (*this)(2)|
  // |   p(0)       p(1)       p(2)   |

  
  // Point val;
  // val.coords[0] =  (coords[1]*p.coords[2] - coords[2]*p.coords[1]);
  // val.coords[1] = -(coords[0]*p.coords[2] - coords[2]*p.coords[0]);
  // val.coords[2] =  (coords[0]*p.coords[1] - coords[1]*p.coords[0]);
  // return val;

  return Point(  coords[1]*p.coords[2] - coords[2]*p.coords[1],
	        -coords[0]*p.coords[2] + coords[2]*p.coords[0],
	         coords[0]*p.coords[1] - coords[1]*p.coords[0]);
};



Point Point::unit() const
{

  const real length = size();
  
  assert (length != static_cast<real>(0.));
  
#if DIM == 1
  return Point(coords[0]/length);
#endif
  
#if DIM == 2 
  return Point(coords[0]/length,
	       coords[1]/length);
#endif
  
#if DIM == 3
  return Point(coords[0]/length,
	       coords[1]/length, 
	       coords[2]/length);
#endif
  
};



void Point::clear()
{
  for (unsigned int i=0; i<DIM; i++)
    coords[i] = 0;
};



void Point::print() const
{
#if DIM == 1
  
  std::cout << "x=" << (*this)(0) << std::endl;
  
#endif
#if DIM == 2
  
  std::cout << "(x,y)=("
	    << std::setw(8) << (*this)(0) << ", "
	    << std::setw(8) << (*this)(1) << ")"
	    << std::endl;

#endif
#if DIM == 3
  
  std::cout <<  "(x,y,z)=("
	    << std::setw(8) << (*this)(0) << ", "
	    << std::setw(8) << (*this)(1) << ", "
	    << std::setw(8) << (*this)(2) << ")"
	    << std::endl;
#endif
};




void Point::write_unformatted (std::ostream &out) const
{
  assert (out);

  out << std::setiosflags(std::ios::showpoint)
      << (*this)(0) << " "
      << (*this)(1) << " "
      << (*this)(2) << std::endl;
};




bool Point::operator == (const Point& rhs) const
{
  Point res = (*this) - rhs;

  if (res.size_sq() < 1.e-12)
    return true;
  
  return false;
};



bool Point::operator < (const Point& rhs) const
{
  const unsigned int big     = static_cast<unsigned int>(1e5);
  const unsigned int smaller = static_cast<unsigned int>(1e2);
  
  int lhs_x = static_cast<int>(ceil(((*this)(0))*big));
  int lhs_y = static_cast<int>(ceil(((*this)(1))*big));
  int lhs_z = static_cast<int>(ceil(((*this)(2))*big));
  
  int rhs_x = static_cast<int>(ceil(rhs(0)*big));
  int rhs_y = static_cast<int>(ceil(rhs(1)*big));
  int rhs_z = static_cast<int>(ceil(rhs(2)*big));
  
  lhs_x -= lhs_x % smaller;
  lhs_y -= lhs_y % smaller;
  lhs_z -= lhs_z % smaller;

  rhs_x -= rhs_x % smaller;
  rhs_y -= rhs_y % smaller;
  rhs_z -= rhs_z % smaller;

  if (lhs_z < rhs_z)
    return true;
  else
    {
      if (lhs_y < rhs_y)
	return true;
      else
	{
	  if (lhs_x < rhs_x)
	    return true;
	  else
	    return false;
	}
    }

  return false;
};



unsigned int Point::key() const
{
  unsigned int tempx,tempy,tempz;

  int i,j=2,cnt=0;
  unsigned int index[3];
  const real deg = 1.e12;
  
  tempx = static_cast<unsigned int>(((*this)(0)*deg));
  tempy = static_cast<unsigned int>(((*this)(1)*deg));
  tempz = static_cast<unsigned int>(((*this)(2)*deg));

  index[0]=0;
  index[1]=0;
  index[2]=0;

  for(i=sizeof(unsigned int)*8-1;i>=0;i--)
    {
      index[j] += (tempx >> i) & 01;
      index[j]  = index[j] << 01;

      if (( cnt % (sizeof(unsigned int)*8) == 0) && (cnt !=0 ) )
	{
	  cnt = 0;
	  j--;
	}
      else
	cnt++;
    
      index[j] += (tempy >> i) & 01;
      index[j]  = index[j] << 01;

      if (( cnt % (sizeof(unsigned int)*8) == 0) && (cnt !=0 ) )
	{
	  cnt = 0;
	  j--;
	}
      else
	cnt++;
    
      index[j] += (tempz >> i) & 01;
      index[j]  = index[j] << 01;

      if (( cnt % (sizeof(unsigned int)*8) == 0) && (cnt !=0 ) )
	{
	  cnt = 0;
	  j--;
	}
      else
	cnt++;
    };

  return index[2];
};

