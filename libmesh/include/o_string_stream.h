// $Id: o_string_stream.h,v 1.1 2003-03-22 21:04:30 ddreyer Exp $

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



#ifndef __o_string_stream_h__
#define __o_string_stream_h__

// C++ includes
#include <sstream>

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
 * features in the \p ostringstream of the older \p GCC
 * versions.  For other compilers, this class is simply 
 * a \p std::ostringstream.
 */

// ------------------------------------------------------------
// OStringStream class definition
class OStringStream : public std::ostringstream
{
public:

  /**
   * Default constructor.
   */
  OStringStream () {};

  /**
   * Default destructor.
   */
  ~OStringStream () {};

  /**
   * convenient typedef
   */
  typedef std::string::size_type sizetype;

  /**
   * @returns a \p OStringStream, where \p s
   * was directed in ragged left style with
   * size \p w.
   */
  OStringStream& left (const sizetype w,
		       const std::string& s);

  /**
   * @returns a \p OStringStream, where \p r
   * was directed in ragged left style with
   * size \p w and precision \p prec.
   */
  OStringStream& left (const sizetype w,
		       const sizetype prec,
		       const Real r);

  /**
   * @returns a \p OStringStream, where \p r
   * was directed in ragged left style with
   * size \p w.
   */
  OStringStream& left (const sizetype w,
		       const int n);

  /**
   * @returns a \p OStringStream, where \p s
   * was directed in ragged right style with
   * size \p w.
   */
  OStringStream& right (const sizetype w,
			const std::string& s);


  /**
   * @returns a \p OStringStream, where \p r
   * was directed in scientific style with
   * size \p w.
   */
  OStringStream& scientific (const sizetype w,
			     const Real r);

protected:

#ifdef BROKEN_IOSTREAM
  /**
   * Appends \p n whitespaces to the string \p s.
   */
  void print_ws (const sizetype n);
#endif

};



// ------------------------------------------------------------
// OStringStream inline methods
inline
OStringStream& OStringStream::left (const sizetype w,
				    const std::string& s)
{
#ifndef BROKEN_IOSTREAM
  this->width(w);
  *this << std::left << s;
#else
  *this << s;
  // pad with whitespaces afterwards
  print_ws ((w-s.size()));
#endif
  return *this;
}



inline
OStringStream& OStringStream::left (const sizetype w,
				    const sizetype prec,
				    const Real r)
{
#ifndef BROKEN_IOSTREAM
  this->width(w);
  this->precision(prec);
  *this << r;
#else
  assert (w < 30);
  char buf[30];  
  char format[8];
  // form the format for r
//  sprintf (format, "%%%d.%df", int((w-prec)/2), prec);
  sprintf (format, "%%.%df", prec);
  // form string as desired
  sprintf (buf, format, r);
  *this << buf;
  // pad with whitespaces afterwards
  print_ws (w-std::string(buf).size());
//  print_ws (w-int((w-prec)/2));
#endif
  return *this;
}




inline
OStringStream& OStringStream::left (const sizetype w,
				    const int n)
{
#ifndef BROKEN_IOSTREAM
  this->width(w);
  *this << n;
#else
  assert (w < 30);
  char buf[30];  
  // form string as desired
  sprintf (buf, "%d", n);
  *this << buf;
  // pad with whitespaces afterwards
  print_ws (w-std::string(buf).size());
#endif
  return *this;

}



inline
OStringStream& OStringStream::right (const sizetype w,
				     const std::string& s)
{
#ifndef BROKEN_IOSTREAM
  this->width(w);
  *this << std::right << s;
#else
  // first pad with whitespaces
  print_ws ((w-s.size()));
  *this << s;
#endif
  return *this;
}



inline
OStringStream& OStringStream::scientific (const sizetype w,
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



#ifdef BROKEN_IOSTREAM
inline
void OStringStream::print_ws (const sizetype n)
{
  for (sizetype i = 0; i < n; i++)
    *this << ' ';
}
#endif



#endif // ifndef __o_string_stream_h__

