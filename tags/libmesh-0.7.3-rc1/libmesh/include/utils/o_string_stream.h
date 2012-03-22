// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Mesh configuration options
#include "libmesh_config.h"



#ifdef LIBMESH_HAVE_SSTREAM

// C++ includes
#include <sstream>

// Local Includes
#include "libmesh_common.h"


// Forward Declarations



/*
 * Some compilers, at least HP \p aCC do not even
 * accept empty classes derived from \p std::ostringstream.
 * Therefore, resort to preprocessor definitions.
 */

#ifndef LIBMESH_BROKEN_IOSTREAM
 /*
  * ---------------------------------------------------------------------------------
  * Everything for a clean iostream
  */

# include <iomanip>

 /*
  * Outputs \p std::string \p d with width
  * \p v left-formatted to stream \p o.
  */
# define OSSStringleft(o,v,d)       (o).width(v);  (o) << std::left << (d)

 /*
  * Outputs \p std::string \p d with width
  * \p v right-formatted to stream \p o.
  */
# define OSSStringright(o,v,d)      (o).width(v);  (o) << std::right << (d)

 /*
  * Outputs \p Real \p d with width \p v and
  * precision \p p to stream \p o (padded with
  * whitespaces).
  */
# define OSSRealleft(o,v,p,d)       (o).width(v);  (o).precision(p); (o).fill(' '); (o) << (d)

 /*
  * Outputs \p Real \p d with width \p v and
  * precision \p p to stream \p o (padded with
  * zeros).
  */
# define OSSRealzeroleft(o,v,p,d)   (o).width(v); (o).precision(p); (o).fill('0'); (o) << std::showpoint << std::left << (d)

 /*
  * Outputs \p Real \p d with width \p v and
  * precision \p p to stream \p o (padded with
  * whitespaces).
  */
# define OSSRealright(o,v,p,d)       (o).width(v);  (o).precision(p); (o).fill(' '); (o) << std::right << (d)

 /*
  * Outputs \p Real \p d with width \p v and
  * precision \p p to stream \p o (padded with
  * zeros).
  */
# define OSSRealzeroright(o,v,p,d)   (o).width(v);  (o).precision(p); (o).fill('0'); (o) << std::right << (d)

 /*
  * Outputs \p Real \p d with width \p v
  * in scientific format to stream \p o.
  */
# define OSSRealscientific(o,v,d)   (o) << std::setw(v) << std::scientific << (d)

 /*
  * Outputs \p int \p d with width
  * \p v to stream \p o.
  */
# define OSSInt(o,v,d)              (o).width(v);  (o) << (d)

 /*
  * class alias
  */
//# define OStringStream           std::ostringstream
typedef std::ostringstream OStringStream;




#else
 /*
  * ---------------------------------------------------------------------------------
  * Everything for broken iostream
  */

# include <cstdio>

 /*
  * Outputs \p std::string \p d with width
  * \p v left-formatted to stream \p o.
  */
# define OSSStringleft(o,v,d)       (o).left( (v), (d) )

 /*
  * Outputs \p std::string \p d with width
  * \p v right-formatted to stream \p o.
  */
# define OSSStringright(o,v,d)      (o).right( (v), (d) )

 /*
  * Outputs \p Real \p d with width \p v and
  * precision \p p to stream \p o (padded with
  * whitespaces).
  */
# define OSSRealleft(o,v,p,d)       (o).left( (v), (p), (d) , ' ')

 /*
  * Outputs \p Real \p d with width \p v and
  * precision \p p to stream \p o (padded with
  * zeros).
  */
# define OSSRealzeroleft(o,v,p,d)       (o).left( (v), (p), (d) , '0')

 /*
  * Outputs \p Real \p d with width \p v and
  * precision \p p to stream \p o (padded with
  * whitespaces).
  */
# define OSSRealright(o,v,p,d)      (o).right( (v), (p), (d) ,' ')

 /*
  * Outputs \p Real \p d with width \p v and
  * precision \p p to stream \p o (padded with
  * zeros).
  */
# define OSSRealzeroright(o,v,p,d)      (o).right( (v), (p), (d) ,'0')

 /*
  * Outputs \p Real \p d with width \p v
  * in scientific format to stream \p o.
  */
# define OSSRealscientific(o,v,d)   (o).scientific( (v), (d) )

 /*
  * Outputs \p int \p d with width
  * \p v to stream \p o.
  */
# define OSSInt(o,v,d)              (o).left( (v), (d) )

//  /*
//   * class alias
//   */
// # define OSSOStringStream           OStringStream


namespace libMesh
{


 /**
  * This class provides a compatibility class for broken
  * features in the \p ostringstream of the older \p GCC
  * versions.  Other compilers do not see this class.
  */

 // ------------------------------------------------------------
 // OStringStream class definition
 class OStringStream : public std::ostringstream
 {
 public:

  /**
   * Default constructor.
   */
  OStringStream () : std::ostringstream() {};

  /**
   * Default destructor.
   */
  ~OStringStream () {};

  /**
   * convenient typedef
   */
  typedef std::string::size_type sizetype;

  /**
   * Outputs in a \p OStringStream, where \p s
   * was directed in ragged left style with
   * size \p w.
   */
  void left (const sizetype w,
	     const std::string& s);

  /**
   * Outputs in a \p OStringStream, where \p r
   * was directed in ragged left style with
   * size \p w and precision \p prec.
   */
  void left (const sizetype w,
	     const sizetype prec,
	     const Real r,
	     const char c = ' ');

  /**
   * Outputs in a \p OStringStream, where \p r
   * was directed in ragged left style with
   * size \p w, padded with chars \p c, defaults
   * to whitespace.
   */
  void left (const sizetype w,
	     const int n);

  /**
   * Outputs in a \p OStringStream, where \p s
   * was directed in ragged right style with
   * size \p w.
   */
  void right (const sizetype w,
	      const std::string& s);


  /**
   * Outputs in a \p OStringStream, where \p r
   * was directed in ragged right style with
   * size \p w and precision \p prec.  The leading
   * characters are \p c, defaults to whitespaces.
   */
  void right (const sizetype w,
	      const sizetype prec,
	      const Real r,
	      const char c = ' ');

  /**
   * Outputs in a \p OStringStream, where \p r
   * was directed in scientific style with
   * size \p w.
   */
  void scientific (const sizetype w,
		   const Real r);

 protected:
  /**
   * Appends \p n chars (defaults to whitespaces)
   * to the string \p s.
   */
  void print_ws (const sizetype n,
		 const char c = ' ');

 };



 // ------------------------------------------------------------
 // OStringStream inline methods
 inline
 void OStringStream::left (const sizetype w,
			   const std::string& s)
 {
   *this << s;
   // pad with whitespaces afterwards
   print_ws ((w-s.size()));
 }


 inline
 void OStringStream::left (const sizetype w,
			   const sizetype prec,
			   const Real r,
			   const char c)
 {
   libmesh_assert (w < 30);
   char buf[30];
   char format[8];
   // form the format for r
   // ALTERNATIVE: sprintf (format, "%%%d.%df", int((w-prec)/2), prec);
   sprintf (format, "%%.%df", prec);
   // form string as desired
   sprintf (buf, format, r);
   *this << buf;
   // pad with whitespaces afterwards
   print_ws (w-std::string(buf).size(), c);
   // ALTERNATIVE: print_ws (w-int((w-prec)/2));
 }


 inline
 void OStringStream::left (const sizetype w,
			   const int n)
 {
   libmesh_assert (w < 30);
   char buf[30];
   // form string as desired
   sprintf (buf, "%d", n);
   *this << buf;
   // pad with whitespaces afterwards
   print_ws (w-std::string(buf).size());
 }


 inline
 void OStringStream::right (const sizetype w,
			    const std::string& s)
 {
   // first pad with whitespaces
   print_ws ((w-s.size()));
   *this << s;
 }


 inline
 void OStringStream::right (const sizetype w,
			    const sizetype prec,
			    const Real r,
			    const char c)
 {
   libmesh_assert (w < 30);
   char buf[30];
   char format[8];
   // form the format for r
   sprintf (format, "%%.%df", prec);
   // form string as desired
   sprintf (buf, format, r);
   // first pad with the user-defined char
   print_ws (w-std::string(buf).size(), c);
   // then print the float
   *this << buf;
 }


 inline
 void OStringStream::scientific (const sizetype w,
				 const Real r)
 {
   libmesh_assert (w < 30);
   char buf[30];
   char format[8];
   // form the format for r
   sprintf (format, "%%%de", w);
   // form string as desired
   sprintf (buf, format, r);
   *this << buf;
 }


 inline
 void OStringStream::print_ws (const sizetype n,
			       const char c)
 {
   for (sizetype i = 0; i < n; i++)
     *this << c;
 }


} // namespace libMesh

#endif // ifndef ... else ... LIBMESH_BROKEN_IOSTREAM




#endif // ifdef LIBMESH_HAVE_SSTREAM

#endif // ifndef __o_string_stream_h__

