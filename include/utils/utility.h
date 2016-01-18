// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_UTILITY_H
#define LIBMESH_UTILITY_H

// Local includes
#include "libmesh/libmesh_common.h" // for Real

// System includes
#include <string>
#include <vector>
#include <algorithm> // for std::lower_bound

namespace libMesh
{


// ------------------------------------------------------------
// The Utility namespace is for functions
// which are useful but don't necessarily belong anywhere else.

namespace Utility
{
//-------------------------------------------------------------------
/**
 * The \p system_info function returns information about the system
 * you are running on.
 */
std::string system_info();



//-------------------------------------------------------------------
/**
 * \p Utility::iota is a duplication of the SGI STL extension
 * \p std::iota.  It simply assigns sequentially increasing values
 * to a range. That is, it assigns \p value to \p *first, \p value + 1
 * to \p *(first + 1) and so on. In general, each iterator \p i in the
 * range [first, last) is assigned \p value + (i - \p first).
 */
template <typename ForwardIter, typename T>
void iota (ForwardIter first, ForwardIter last, T value)
{
  while (first != last)
    {
      *first = value++;
      ++first;
    }
}


/**
 * Utility::is_sorted mimics the behavior of the SGI STL extension
 * std::is_sorted.  Checks to see if the range [first,last) is
 * sorted in non-decreasing order, ie. for each "i" in
 * [first,last) *i <= *(i+1).
 */
template< class InputIterator >
bool is_sorted(InputIterator first, InputIterator last)
{
  if ( first == last )
    return true;

  // "prev" always points to the entry just to the left of "first"
  //  [-    -    -    -    -    -]
  //   ^    ^
  // prev first
  //
  //  [-    -    -    -    -    -]
  //        ^    ^
  //      prev first
  //
  //  [-    -    -    -    -    -]
  //             ^    ^
  //           prev first
  InputIterator prev( first );
  for ( ++first; first != last; ++prev, ++first )
    if ( *first < *prev  )    // Note: this is the same as *prev > *first,
      return false;        // but we only require op< to be defined.

  // If we haven't returned yet, it's sorted!
  return true;


  // A one-liner version using adjacent_find.  This doesn't work for
  // C-style arrays, since their pointers do not have a value_type.
  //
  // Works by checking to see if adjacent entries satisfy *i >
  // *(i+1) and returns the first one which does.  If "last" is
  // returned, no such pair was found, and therefore the range must
  // be in non-decreasing order.
  //
  // return (last ==
  // std::adjacent_find(first, last,
  // std::greater< typename InputIterator::value_type >()));

  // A second one-linear attempt.  This one checks for a **strictly
  // increasing** (no duplicate entries) range.  Also doesn't work
  // with C-style arrays.
  //
  // return (last ==
  // std::adjacent_find(first, last,
  // std::not2(std::less<typename InputIterator::value_type>())));
}


/**
 * The STL provides binary_search() which returns true/false depending
 * on whether the searched-for value is found.  Utility::binary_find() uses a
 * binary search on a sorted range to return an iterator to the searched-for
 * element, or "last" if the element is not found.
 */
template<class ForwardIterator, class T>
ForwardIterator binary_find(ForwardIterator first, ForwardIterator last, const T & value)
{
  ForwardIterator it = std::lower_bound(first, last, value);
  return (it == last || value < *it) ? last : it;
}

/**
 * As above, but takes a custom comparison object.
 */
template<class ForwardIterator, class T, class Compare>
ForwardIterator binary_find(ForwardIterator first, ForwardIterator last, const T & value, Compare comp)
{
  ForwardIterator it = std::lower_bound(first, last, value, comp);
  return (it == last || comp(value,*it)) ? last : it;
}


//-------------------------------------------------------------------
/**
 * An efficient template instantiation for raising
 * to an arbitrary integer power.
 */
template <int N, typename T>
struct do_pow {
  static inline T apply (const T & x)
  {
    libmesh_assert(N>1);

    if (N%2) // odd exponent
      return x * do_pow<N-1,T>::apply(x);

    const T xNover2 = do_pow<N/2,T>::apply(x);

    return xNover2*xNover2;
  }
};

// An efficient compiler would distill N=6 down to 3
// multiplications, but an inefficient one (or a complicated
// T::operator*) might do worse, so we'll specialize here.
template <typename T>
struct do_pow<6,T> {
  static inline T apply (const T & x)
  {
    const T x2 = x*x,
      x4 = x2*x2;

    return x4*x2;
  }
};

template <typename T>
struct do_pow<1,T> {
  static inline T apply (const T & x) { return x; }
};

template <typename T>
struct do_pow<0,T> {
  static inline T apply (const T &) { return 1; }
};


template <int N, typename T>
inline
T pow(const T & x)
{
  return do_pow<N,T>::apply(x);
}

//-------------------------------------------------------------------
/**
 * A simple implementation of the factorial.
 */
inline
unsigned int factorial(unsigned int n)
{

  unsigned int factorial_n = 1;

  if (n==0)
    return factorial_n;

  for (unsigned int i=1; i<n; i++)
    factorial_n *= i+1;

  return factorial_n;
}


//-------------------------------------------------------------------
/**
 * A convenient method to truly empty a vector using the "swap trick"
 */
template <typename T>
void deallocate (std::vector<T> & vec)
{
  std::vector<T>().swap(vec);
}


//-------------------------------------------------------------------
// Utility functions useful when dealing with complex numbers.

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

/**
 * @returns for \p r_o_c = 0 the filename for output of the real part
 * of complex data, and for  \p r_o_c = 1 the filename for the imaginary
 * part.
 */
std::string complex_filename (const std::string & basename,
                              unsigned int r_o_c=0);

/**
 * Prepare complex data for writing.
 */
void prepare_complex_data (const std::vector<Complex> & source,
                           std::vector<Real> & real_part,
                           std::vector<Real> & imag_part);

#endif // #ifdef LIBMESH_USE_COMPLEX_NUMBERS



//-------------------------------------------------------------------
/**
 * This Functor simply takes an object and reverses its byte
 * representation.  This is useful for changing endian-ness
 * for file IO.  This class has been tested on x86 architectures
 * with 4-byte words.
 *
 *
 */
class ReverseBytes
{
public:

  /**
   * Constructor.  Takes a bool, determines if we will actually
   * do byte reversing.
   */
  explicit
  ReverseBytes (const bool dr);

  /**
   * Functor.  Takes the data to reverse and performs the
   * byte-ordering reversal.
   */
  template <typename T>
  T operator () (T & data) const;

private:

  /**
   * Returns the value of the reverse flag.
   */
  bool reverse () const { return _do_reverse; }

  /**
   * flag
   */
  const bool _do_reverse;
};



//---------------------------------------------------------
// ReverseBytes inline members
inline
ReverseBytes::ReverseBytes (const bool rb) :
  _do_reverse (rb)
{}


template <typename T>
inline
T ReverseBytes::operator() (T & data) const
{
  // Possibly reverse the byte ordering
  if (this->reverse())
    {
      unsigned char * b = (unsigned char *) &data;

      int i=0;
      int j=(sizeof(T) - 1);

      while (i < j)
        {
          std::swap (b[i], b[j]);
          i++; j--;
        }
    }

  return data;
}


}

} // namespace libMesh

#endif // LIBMESH_UTILITY_H
