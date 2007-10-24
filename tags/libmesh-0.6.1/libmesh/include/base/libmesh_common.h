// $Id: libmesh_common.h,v 1.22 2007-10-21 20:48:40 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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

// C/C++ includes everyone should know about
#include <iostream> // needed for std::cout, std::cerr
#include <complex>
#include <cassert>
#ifdef HAVE_STDLIB_H
# include <cstdlib>
#endif

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

//#ifdef REAL
//#  undef REAL
//#endif

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

#define DOUBLE_PRECISION

// Define a corresponding tolerance.  This is what should be
// considered "good enough" when doing floating point comparisons.
// For example, v == 0 is changed to std::abs(v) < TOLERANCE.

#ifndef SINGLE_PRECISION
  #ifdef TRIPLE_PRECISION
    typedef long double Real;
  # define TOLERANCE 1.e-8
  # define MPI_REAL MPI_LONG_DOUBLE
  #else
    typedef double Real;
  # define TOLERANCE 1.e-6
  # define MPI_REAL MPI_DOUBLE
  #endif
#else
  typedef float Real;
  # define TOLERANCE 1.e-3
  # define MPI_REAL MPI_FLOAT
#endif

// For some reason the real std::max, std::min
// don't handle mixed compatible types
namespace std {
  inline long double max(long double a, double b)
  { return (a>b?a:b); }
  inline long double min(long double a, double b)
  { return (a<b?a:b); }

  inline long double max(double a, long double b)
  { return (a>b?a:b); }
  inline long double min(double a, long double b)
  { return (a<b?a:b); }

  inline double max(double a, float b)
  { return (a>b?a:b); }
  inline double min(double a, float b)
  { return (a<b?a:b); }

  inline double max(float a, double b)
  { return (a>b?a:b); }
  inline double min(float a, double b)
  { return (a<b?a:b); }

  inline long double max(long double a, float b)
  { return (a>b?a:b); }
  inline long double min(long double a, float b)
  { return (a<b?a:b); }

  inline long double max(float a, long double b)
  { return (a>b?a:b); }
  inline long double min(float a, long double b)
  { return (a<b?a:b); }
}

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


// Define the value type for error estimates.
// Since AMR/C decisions don't have to be precise,
// we default to float for memory efficiency.
typedef float ErrorVectorReal;
#define MPI_ERRORVECTORREAL MPI_FLOAT


#ifdef HAVE_MPI
namespace libMesh
{
  /**
   * MPI Communicator to be used in the library.
   */
  extern MPI_Comm COMM_WORLD;
}
#endif

// These are useful macros that behave like functions in the code.
// If you want to make sure you are accessing a section of code just
// stick a here(); in it, for example
#undef here
#define here()     { std::cout << "[" << libMesh::processor_id() << "] " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; }

// the stop() macro will stop the code until a SIGCONT signal is recieved.  This is useful, for example, when
// determining the memory used by a given operation.  A stop() could be instered before and after a questionable
// operation and the delta memory can be obtained from a ps or top.  This macro only works for serial cases.
#undef stop
#ifdef HAVE_CSIGNAL
#  include <csignal>
#  define stop()     { if (libMesh::n_processors() == 1) { here(); std::cout << "Stopping process " << getpid() << "..." << std::endl; std::raise(SIGSTOP); std::cout << "Continuing process " << getpid() << "..." << std::endl; } }
#else
#  define stop()     { if (libMesh::n_processors() == 1) { here(); std::cerr << "WARNING:  stop() does not work without the <csignal> header file!" << std::endl; } }
#endif


// The error() macro prints a message and aborts the code
#undef error
#ifdef HAVE_MPI
#  define error()    { std::cerr << "[" << libMesh::processor_id() << "] " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; if (libMesh::n_processors() > 1) MPI_Abort(libMesh::COMM_WORLD,1); std::abort(); }
#else
#  define error()    { std::cerr << "[" << libMesh::processor_id() << "] " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; std::abort(); }
#endif

// The untested macro warns that you are using untested code
#undef untested
#define untested() { std::cout << "*** Warning, This code is untested, experimental, or likely to see future API changes: " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << " ***" << std::endl; }


// The deprecated macro warns that you are using deprecated code
#undef deprecated
#define deprecated() { std::cout << "*** Warning, This Code is Deprecated! " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << " ***" << std::endl; }




// 3D spatial dimension unless otherwise specified
#ifndef DIM
#  define DIM 3
#endif






#endif // #define __libmesh_common_h__
