// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// libMesh includes
#include "libmesh/libmesh_C_isnan.h"

// C includes
#include <math.h>
#include <stdio.h>

#ifdef isnan
int libmesh_C_isnan_float(float a) {
  printf("libmesh_C_isnan_float is deprecated and should not be used!\n");
  return isnan(a);
}
int libmesh_C_isnan_double(double a) {
  printf("libmesh_C_isnan_double is deprecated and should not be used!\n");
  return isnan(a);
}
int libmesh_C_isnan_longdouble(long double a) {
  printf("libmesh_C_isnan_longdouble is deprecated and should not be used!\n");
  return isnan(a);
}

#else
int libmesh_C_isnan_float(float a) {
  printf("libmesh_C_isnan_float is deprecated and should not be used!\n");
  return (a != a);
}
int libmesh_C_isnan_double(double a) {
  printf("libmesh_C_isnan_double is deprecated and should not be used!\n");
  return (a != a);
}
int libmesh_C_isnan_longdouble(long double a) {
  printf("libmesh_C_isnan_longdouble is deprecated and should not be used!\n");
  return (a != a);
}

#endif
