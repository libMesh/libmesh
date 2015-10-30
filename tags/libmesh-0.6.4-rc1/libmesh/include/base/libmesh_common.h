// $Id$

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



#ifndef __libmesh_common_h__
#define __libmesh_common_h__


// The library configuration options
#include "libmesh_config.h"

// Include the MPI definition
#ifdef LIBMESH_HAVE_MPI
# include <mpi.h>
#endif

// C/C++ includes everyone should know about
#include <iostream> // needed for std::cout, std::cerr
#include <fstream> // needed for argument to print_trace()
#include <unistd.h>  // needed for getpid()
#include <complex>
// #include <cassert>  // Use libmesh_assert() now
#ifdef LIBMESH_HAVE_STDLIB_H
# include <cstdlib>
#endif
#include <typeinfo> // std::bad_cast

// _basic_ library functionality
#include "libmesh_base.h"
#include "libmesh_exceptions.h"

#ifdef LIBMESH_ENABLE_TRACEFILES
#  include "print_trace.h"
#endif





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


// Helper functions for complex/real numbers
// to clean up #ifdef LIBMESH_USE_COMPLEX_NUMBERS elsewhere
template<typename T> inline T libmesh_real(T a) { return a; }
template<typename T> inline T libmesh_norm(T a) { return a*a; }

template<typename T>
inline T libmesh_real(std::complex<T> a) { return std::real(a); }

template<typename T>
inline T libmesh_norm(std::complex<T> a) { return std::norm(a); }


// Define the value type for unknowns in simulations.
// This is either Real or Complex, depending on how
// the library was configures
#if   defined (LIBMESH_USE_REAL_NUMBERS)
  typedef Real Number;
#elif defined (LIBMESH_USE_COMPLEX_NUMBERS)
  typedef Complex Number;
#else
  DIE A HORRIBLE DEATH HERE...
#endif


// Define the value type for error estimates.
// Since AMR/C decisions don't have to be precise,
// we default to float for memory efficiency.
typedef float ErrorVectorReal;
#define MPI_ERRORVECTORREAL MPI_FLOAT


#ifdef LIBMESH_HAVE_MPI
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
// stick a libmesh_here(); in it, for example

#define libmesh_here()     do { std::cerr << "[" << libMesh::processor_id() << "] " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; } while (0)

// the libmesh_stop() macro will stop the code until a SIGCONT signal
// is recieved.  This is useful, for example, when determining the
// memory used by a given operation.  A libmesh_stop() could be
// instered before and after a questionable operation and the delta
// memory can be obtained from a ps or top.  This macro only works for
// serial cases.
#ifdef LIBMESH_HAVE_CSIGNAL
#  include <csignal>
#  define libmesh_stop()     do { if (libMesh::n_processors() == 1) { libmesh_here(); std::cout << "Stopping process " << getpid() << "..." << std::endl; std::raise(SIGSTOP); std::cout << "Continuing process " << getpid() << "..." << std::endl; } } while(0)
#else
#  define libmesh_stop()     do { if (libMesh::n_processors() == 1) { libmesh_here(); std::cerr << "WARNING:  libmesh_stop() does not work without the <csignal> header file!" << std::endl; } } while(0)
#endif

// The libmesh_assert() macro acts like C's assert(), but throws a 
// libmesh_error() (including stack trace, etc) instead of just exiting
#ifdef NDEBUG
#define libmesh_assert(asserted) 
#else
#define libmesh_assert(asserted)  do { if (!(asserted)) { std::cerr << "Assertion `" #asserted "' failed." << std::endl; libmesh_error(); } } while(0)
#endif

// The libmesh_write_traceout() macro writes stack trace files, if
// that feature has been configured.
#ifdef LIBMESH_ENABLE_TRACEFILES
  #define libmesh_write_traceout()   do { std::stringstream outname; outname << "traceout_" << libMesh::processor_id() << '_' << getpid() << ".txt"; std::ofstream traceout(outname.str().c_str()); print_trace(traceout); } while(0)
#else
  #define libmesh_write_traceout()   do { } while (0)
#endif

// The libmesh_error() macro prints a message and throws a LogicError
// exception
//
// The libmesh_not_implemented() macro prints a message and throws a
// NotImplemented exception
//
// The libmesh_file_error(const std::string& filename) macro prints a message
// and throws a FileError exception
//
// The libmesh_convergence_failure() macro prints a message
// and throws a ConvergenceFailure exception
#define libmesh_error()    do { libmesh_write_traceout(); libmesh_here(); LIBMESH_THROW(libMesh::LogicError()); } while(0)
#define libmesh_not_implemented()    do { libmesh_write_traceout(); libmesh_here(); LIBMESH_THROW(libMesh::NotImplemented()); } while(0)
#define libmesh_file_error(filename)    do { libmesh_write_traceout(); libmesh_here(); LIBMESH_THROW(libMesh::FileError(filename)); } while(0)
#define libmesh_convergence_failure()    do { libmesh_here(); LIBMESH_THROW(libMesh::ConvergenceFailure()); } while(0)

// The libmesh_cast functions do a dynamic cast and assert the result,
// if we're in debug or development modes, but just do a faster static
// cast if we're in optimized mode.  Use these casts when you're
// certain that a cast will succeed in correct code but you want to be
// able to double-check.
template <typename Tnew, typename Told>
inline Tnew libmesh_cast_ref(Told& oldvar)
{
#ifndef NDEBUG
  try
    {
      Tnew newvar = dynamic_cast<Tnew>(oldvar);
      return newvar;
    }
  catch (std::bad_cast)
    {
      libmesh_assert (false);
    }
#else
  return(static_cast<Tnew>(oldvar));
#endif
}

// We use two different function names to avoid an odd overloading
// ambiguity bug with icc 10.1.008
template <typename Tnew, typename Told>
inline Tnew libmesh_cast_ptr (Told* oldvar)
{
#ifndef NDEBUG
  Tnew newvar = dynamic_cast<Tnew>(oldvar);
  libmesh_assert (newvar);
  return newvar;
#else
  return(static_cast<Tnew>(oldvar));
#endif
}


// The untested macro warns that you are using untested code
#undef untested
#define untested() do { std::cout << "*** Warning, This code is untested, experimental, or likely to see future API changes: " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << " ***" << std::endl; } while (0)


// The deprecated macro warns that you are using deprecated code
#undef deprecated
#define deprecated() do { std::cout << "*** Warning, This Code is Deprecated! " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << " ***" << std::endl; } while(0)




// 3D spatial dimension unless otherwise specified
#ifndef LIBMESH_DIM
#  define LIBMESH_DIM 3
#endif






#endif // #define __libmesh_common_h__
