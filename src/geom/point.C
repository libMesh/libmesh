// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
// #include "libmesh/point.h"




// ------------------------------------------------------------
// Point class member functions
// unsigned int Point::key() const
// {
//   unsigned int tempx,tempy,tempz;

//   int i,j=2,cnt=0;
//   unsigned int index[3];
//   const Real deg = 1.e12;

//   tempx = static_cast<unsigned int>(((*this)(0)*deg));
//   tempy = static_cast<unsigned int>(((*this)(1)*deg));
//   tempz = static_cast<unsigned int>(((*t!his)(2)*deg));

//   index[0]=0;
//   index[1]=0;
//   index[2]=0;

//   for (i=sizeof(unsigned int)*8-1;i>=0;i--)
//     {
//       index[j] += (tempx >> i) & 01;
//       index[j]  = index[j] << 01;

//       if (( cnt % (sizeof(unsigned int)*8) == 0) && (cnt !=0 ) )
// {
//   cnt = 0;
//   j--;
// }
//       else
// cnt++;

//       index[j] += (tempy >> i) & 01;
//       index[j]  = index[j] << 01;

//       if (( cnt % (sizeof(unsigned int)*8) == 0) && (cnt !=0 ) )
// {
//   cnt = 0;
//   j--;
// }
//       else
// cnt++;

//       index[j] += (tempz >> i) & 01;
//       index[j]  = index[j] << 01;

//       if (( cnt % (sizeof(unsigned int)*8) == 0) && (cnt !=0 ) )
// {
//   cnt = 0;
//   j--;
// }
//       else
// cnt++;
//     }

//   return index[2];
// }
