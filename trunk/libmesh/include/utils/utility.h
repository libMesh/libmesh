// $Id: utility.h,v 1.3 2004-02-09 17:12:28 jwpeterson Exp $

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


#ifndef __utility_h__
#define __utility_h__

// Local includes

// System includes
#include <string>


// ------------------------------------------------------------
// The Utility namespace is for functions
// which are useful but don't necessarily belong anywhere else.

namespace Utility
{
  /**
   * The \p system_info function returns information about the system
   * you are running on.
   */
  std::string system_info();

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
   * An efficient, temporary-free swapping algorithm
   * on machines where this is possible.
   */
  template <typename T>
  void swap (T& a, T& b)
  {
#ifdef __HP_aCC
    
    const T buf = a;
    a = b;
    b = buf;
    
#else

    a^=b^=a^=b;
      
#endif    
  }

  /**
   * An efficient template instantiation for raising
   * to an arbitrary integer power.
   */
  template <int N>
  inline
  Real pow(const Real x) { return x * pow<N-1>(x); }

  /**
   * You have to also provide a full specialization for
   * raising to the zero power which returns 1.  Otherwise,
   * the function above will simply expand to the maximum
   * template depth on your machine and cause a compilation
   * error.
   */
  template <>
  inline
  Real pow<0>(const Real) { return 1.; }
}

#endif
