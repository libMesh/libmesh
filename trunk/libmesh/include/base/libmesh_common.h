// $Id: libmesh_common.h,v 1.3 2004-08-06 16:48:39 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __libmesh_common_h__
#define __libmesh_common_h__


// The library configuration options
#include "libmesh_config.h"

// C++ includes everyone should know about
#include <complex>

// Include the MPI definition
#ifdef HAVE_MPI
# include <mpi.h>
#endif

// _basic_ library functionality
#include "libmesh_base.h"





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

#ifdef MPI_REAL
#  undef MPI_REAL
#endif

// Check to see if TOLERANCE has been defined by another
// package, if so we might want to change the name...
#ifdef TOLERANCE
   DIE A HORRIBLE DEATH HERE...
#  undef TOLERANCE
#endif


   
// Define the type to use for real numbers
#ifndef SINGLE_PRECISION
  typedef double Real;
  typedef double REAL;
# define MPI_REAL MPI_DOUBLE
#else
  typedef float Real;
  typedef float REAL;
# define MPI_REAL MPI_FLOAT
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

  


// These are useful macros that behave lilke functions in the code.
// If you want to make sure you are accessing a section of code just
// stick a here(); in it, for example
#undef here
#define here()     { std::cout << "[" << libMesh::processor_id() << "] " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; }

#undef error
#ifdef HAVE_MPI
#  define error()    { std::cerr << "[" << libMesh::processor_id() << "] " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; abort(); }
#else
#  define error()    { std::cerr << "[" << libMesh::processor_id() << "] " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; abort(); }
#endif

#undef untested
#define untested() { std::cout << "*** Using untested code: " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << " ***" << std::endl; }




// 3D spatial dimension unless otherwise specified
#ifndef DIM
#  define DIM 3
#endif



// Define a tolerance.  This is what should be considered "good enough"
// when doing floating point comparisons.  For example, v == 0 is
// changed to fabs(v) < TOLERANCE.
#define TOLERANCE 1.e-6





#endif // #define __libmesh_common_h__
