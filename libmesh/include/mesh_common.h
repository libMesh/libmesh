// $Id: mesh_common.h,v 1.3 2003-01-20 17:06:12 jwpeterson Exp $

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



#ifndef __mesh_common_h__
#define __mesh_common_h__


// C++ includes everyone should know about
#include <iostream>
#include <assert.h>
#include <stdlib.h>


#include "mesh_config.h"

/**
 * An anonymous namespace to hold \p typedefs and the like
 * for the library
 */
namespace
{
  
#ifdef REAL
#  undef REAL
#endif

#ifdef NUMBER
#  undef NUMBER
#endif


#ifndef SINGLE_PRECISION
  typedef double real;
  typedef double REAL;
#  ifdef USE_COMPLEX_NUMBERS
  // note that PETSc also uses <complex>, check petscmath.h
#    include <complex>
     typedef std::complex<double> number;
     typedef std::complex<double> NUMBER;
#  else
     typedef double number;
     typedef double NUMBER;
#  endif
#else
  typedef float real;
  typedef float REAL;
#  ifdef USE_COMPLEX_NUMBERS
     // this is _not_ supported by PETSc!
     CHOKE_THIS!
#    include <complex>
     typedef std::complex<float> number;
     typedef std::complex<float> NUMBER;
#  else
     typedef float number;
     typedef float NUMBER;
#  endif
#endif


#undef here
#undef error
#define here()     { std::cout << __FILE__ << ", line " << __LINE__ << std::endl; }
#define error()    { here(); abort(); }
#define untested() { std::cout << "*** Using untested code: " << __FILE__ << ", line " << __LINE__ << std::endl; }


// 3D spatial dimension unless otherwise specified
#ifndef DIM
#  define DIM 3
#endif

}; // end anonymous namespace

#endif // #define __mesh_common_h__
