// $Id: o_f_stream.h,v 1.1 2003-03-22 21:04:30 ddreyer Exp $

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



#ifndef __o_f_stream_h__
#define __o_f_stream_h__

// C++ includes
#include <fstream>
#include <string>
#include <vector>

// Local Includes
#include "mesh_common.h"


/*
 * Need sprintf for scientific format, while other 
 * compilers need iomanip
 */
#ifdef BROKEN_IOSTREAM
# include <stdio.h>
#else
# include <iomanip>
#endif


// Forward Declarations



/**
 * This class provides a compatibility class for broken
 * features in the \p std::ofstream of the older \p GCC
 * versions.  For other compilers, this class is simply 
 * a \p std::ofstream.
 */

// ------------------------------------------------------------
// OFStream class definition
class OFStream : public std::ofstream
{
public:

  /**
   * Default constructor.
   */
  OFStream () {};

  /**
   * Default destructor.
   */
  ~OFStream () {};

  /**
   * convenient typedef
   */
  typedef std::string::size_type sizetype;

  /**
   * @returns a \p OFStream, where \p r
   * was directed in scientific style with
   * size \p w.
   */
  OFStream& scientific (const sizetype w,
			const Real r);

  /**
   * @returns a \p OFStream, where all entries of
   * \p v were directed in scientific style with
   * size \p w, each float separated from the other
   * by \p sep.  Works also for complex numbers,
   * if enabled.  Note that the separator should 
   * contain at most 10 characters, and defaults
   * to just one whitespace.
   */
  template <typename T>
  OFStream& scientific (const sizetype w,
			const std::vector<T>& v,
			const char* sep = " ");

};



// ------------------------------------------------------------
// OFStream inline methods
inline
OFStream& OFStream::scientific (const sizetype w,
				const Real r)
{
#ifndef BROKEN_IOSTREAM
  *this << std::setw(w)
	<< std::scientific
	<< r;
#else
  assert (w < 30);
  char buf[30];  
  char format[8];
  // form the format for r
  sprintf (format, "%%%de", w);
  // form string as desired
  sprintf (buf, format, r);
  *this << buf;
#endif
  return *this;
}


#endif // ifndef __o_f_stream_h__

