// $Id: utility.h,v 1.6 2004-03-20 15:16:56 benkirk Exp $

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

// System includes
#include <string>
#include <vector>

// Local includes
#include "libmesh_common.h" // for Real


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
   * An efficient template instantiation for raising
   * to an arbitrary integer power.
   */
  template <int N>
  inline
  Real pow(const Real x) { assert(N>0); return x * pow<N-1>(x); }

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


  //-------------------------------------------------------------------
  // Utility functions useful when dealing with complex numbers.
  
#ifdef USE_COMPLEX_NUMBERS

  /**
   * @returns for \p r_o_c = 0 the filename for output of the real part
   * of complex data, and for  \p r_o_c = 1 the filename for the imaginary 
   * part.
   */
  const char* complex_filename (std::string name,
				unsigned int r_o_c=0);

  /**
   * Prepare complex data for writing.
   */
  void prepare_complex_data (const std::vector<Complex>& source,
			     std::vector<Real>& real_part,
			     std::vector<Real>& imag_part);

#endif // #ifdef USE_COMPLEX_NUMBERS

}

#endif // #define __utility_h__
