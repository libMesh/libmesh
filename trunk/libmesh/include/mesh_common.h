// $Id: mesh_common.h,v 1.7 2003-02-04 00:49:09 benkirk Exp $

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

// The library configuration options
#include "mesh_config.h"


  
#ifdef Real
#  undef Real
#endif

#ifdef REAL
#  undef REAL
#endif

#ifdef Complex
#  undef Complex
#endif

#ifdef COMPLEX
#  undef COMPLEX
#endif


#ifndef SINGLE_PRECISION
  typedef double Real;
  typedef double REAL;
#  ifdef USE_COMPLEX_NUMBERS
#    include <complex>
     typedef std::complex<double> Complex;
     typedef std::complex<double> COMPLEX;
#  else
     typedef double Complex;
     typedef double COMPLEX;
#  endif
#else
  typedef float Real;
  typedef float REAL;
#  ifdef USE_COMPLEX_NUMBERS
     // this is _not_ supported by PETSc!
     CHOKE_THIS!
#  else
     typedef float Complex;
     typedef float COMPLEX;
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

#endif // #define __mesh_common_h__
