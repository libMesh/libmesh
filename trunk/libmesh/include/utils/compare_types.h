
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


#ifndef __compare_types_h__
#define __compare_types_h__

// System includes
#include <complex>

namespace libMesh
{

// Copy of boost's enable_if_c

namespace boostcopy {
  template <bool B, class T = void>
    struct enable_if_c {
      typedef T type;
    };

  template <class T>
    struct enable_if_c<false, T> {};
}



// Complete list of scalar classes, needed for disambiguation
template <typename T>
struct ScalarTraits {
      static const bool value = false;
};

template<>
struct ScalarTraits<signed char> { static const bool value = true; };

template<>
struct ScalarTraits<char> { static const bool value = true; };

template<>
struct ScalarTraits<short> { static const bool value = true; };

template<>
struct ScalarTraits<int> { static const bool value = true; };

template<>
struct ScalarTraits<long> { static const bool value = true; };

template<>
struct ScalarTraits<unsigned char> { static const bool value = true; };

template<>
struct ScalarTraits<unsigned short> { static const bool value = true; };

template<>
struct ScalarTraits<unsigned int> { static const bool value = true; };

template<>
struct ScalarTraits<unsigned long> { static const bool value = true; };

template<>
struct ScalarTraits<float> { static const bool value = true; };

template<>
struct ScalarTraits<double> { static const bool value = true; };

template<>
struct ScalarTraits<long double> { static const bool value = true; };

template<>
struct ScalarTraits<std::complex<float> > { static const bool value = true; };

template<>
struct ScalarTraits<std::complex<double> > { static const bool value = true; };

template<>
struct ScalarTraits<std::complex<long double> > { static const bool value = true; };



// Operators using different but compatible types need a return value
// based on whichever type the other can be upconverted into.  For
// instance, the proper return type for
// TypeVector<float>::operator*(double) is TypeVector<double>.  In
// general, an operation using types S and T should return a value
// based on CompareTypes<S,T>::supertype

template<typename S, typename T>
struct CompareTypes {
  typedef void supertype;
};

// There's got to be some magic template way to do these better - but the best
// thing on the net requires a bunch of Alexandrescu's code and doesn't work
// with older compilers

#define SUPERTYPE(mysub,mysuper) \
	template<> \
	struct CompareTypes<mysub, mysuper> { \
	  typedef mysuper supertype; \
	}; \
	template<> \
	struct CompareTypes<mysuper, mysub> { \
	  typedef mysuper supertype; \
	}; \
	template<> \
	struct CompareTypes<std::complex<mysub>, mysuper> { \
	  typedef std::complex<mysuper> supertype; \
	}; \
	template<> \
	struct CompareTypes<mysuper, std::complex<mysub> > { \
	  typedef std::complex<mysuper> supertype; \
	}; \
	template<> \
	struct CompareTypes<mysub, std::complex<mysuper> > { \
	  typedef std::complex<mysuper> supertype; \
	}; \
	template<> \
	struct CompareTypes<std::complex<mysuper>, mysub> { \
	  typedef std::complex<mysuper> supertype; \
	}; \
	template<> \
	struct CompareTypes<std::complex<mysub>, std::complex<mysuper> > { \
	  typedef std::complex<mysuper> supertype; \
	}; \
	template<> \
	struct CompareTypes<std::complex<mysuper>, std::complex<mysub> > { \
	  typedef std::complex<mysuper> supertype; \
	}

template<typename T>
struct CompareTypes<T, T> {
  typedef T supertype;
};

template<typename T>
struct CompareTypes<T, std::complex<T> > {
  typedef std::complex<T> supertype;
};

template<typename T>
struct CompareTypes<std::complex<T>, T> {
  typedef std::complex<T> supertype;
};

SUPERTYPE(unsigned char, short);
SUPERTYPE(unsigned char, int);
SUPERTYPE(unsigned char, float);
SUPERTYPE(unsigned char, double);
SUPERTYPE(unsigned char, long double);
SUPERTYPE(unsigned short, int);
SUPERTYPE(unsigned short, float);
SUPERTYPE(unsigned short, double);
SUPERTYPE(unsigned short, long double);
SUPERTYPE(unsigned int, float);
SUPERTYPE(unsigned int, double);
SUPERTYPE(unsigned int, long double);
SUPERTYPE(char, short);
SUPERTYPE(char, int);
SUPERTYPE(char, float);
SUPERTYPE(char, double);
SUPERTYPE(char, long double);
SUPERTYPE(short, int);
SUPERTYPE(short, float);
SUPERTYPE(short, double);
SUPERTYPE(short, long double);
SUPERTYPE(int, float);
SUPERTYPE(int, double);
SUPERTYPE(int, long double);
SUPERTYPE(float, double);
SUPERTYPE(float, long double);
SUPERTYPE(double, long double);

// gcc can't tell which of the following is the most specialized?  Weak.
/*
template<typename S, typename T>
struct CompareTypes<std::complex<S>, std::complex<T> > {
  typedef std::complex<typename CompareTypes<S, T>::supertype> supertype;
};

template<typename S, typename T>
struct CompareTypes<std::complex<S>, T> {
  typedef std::complex<typename CompareTypes<S, T>::supertype> supertype;
};

template<typename S, typename T>
struct CompareTypes<S, std::complex<T> > {
  typedef std::complex<typename CompareTypes<S, T>::supertype> supertype;
};
*/

} // namespace libMesh

#endif // __compare_types_h__
