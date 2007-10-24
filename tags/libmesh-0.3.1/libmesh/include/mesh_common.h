// $Id: mesh_common.h,v 1.8 2003-02-20 04:59:58 benkirk Exp $

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
#include <complex>
#include <assert.h>
#include <stdlib.h>

// The library configuration options
#include "mesh_config.h"



// Undefine any existing macros
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


// Define the type to use for real numbers
#ifndef SINGLE_PRECISION
  typedef double Real;
  typedef double REAL;
#else
  typedef float Real;
  typedef float REAL;
#endif

// Define the type to use for complex numbers
// Always use std::complex<double>, as required by Petsc
typedef std::complex<double> Complex;
typedef std::complex<double> COMPLEX;


// Define the value type for unknowns in simulations.
// This is either Real or Complex, depending on how
// the library was configures
#if   defined (USE_REAL_NUMBERS)
  typedef Real Number;
#elif defined (USE_COMPLEX_NUMBERS)
  typedef Complex Number;
#else
  DIE A HORRIBLE DEATH HERE...
#endif



#undef here
#undef error
#define here()     { std::cout << __FILE__ << ", line " << __LINE__ << std::endl; }
#define error()    { here(); abort(); }
#define untested() { std::cout << "*** Using untested code: " << __FILE__ << ", line " << __LINE__ << " ***" << std::endl; }


// 3D spatial dimension unless otherwise specified
#ifndef DIM
#  define DIM 3
#endif

#endif // #define __mesh_common_h__
