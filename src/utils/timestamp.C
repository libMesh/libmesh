// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_config.h"

// C++ includes
#include <ctime>
#include <sstream>
#ifdef LIBMESH_HAVE_LOCALE
#include <locale>
#endif

// Local includes
#include "libmesh/timestamp.h"

namespace libMesh
{

namespace Utility
{
// If the locale header is available, we will use the "C++" way of
// creating a timestamp, otherwise, we'll fall back on the C way.
std::string get_timestamp()
{
#ifdef LIBMESH_HAVE_LOCALE
  // Create time_put "facet"
  std::locale loc;
  const std::time_put<char> & tp = std::use_facet <std::time_put<char>> (loc);

  // Call C-style time getting functions
  time_t now = time(NULL);
  tm * tm_struct = localtime(&now);

  // Date will eventually be stored in this ostringstream's string
  std::ostringstream date_stream;

  // See below for documentation on the use of the
  // std::time_put::put() function
  tp.put(date_stream,        /*s*/
         date_stream,        /*str*/
         date_stream.fill(), /*fill*/
         tm_struct,          /*tm*/
         'c');               /*format*/

  // Another way to use it is to totally customize the format...
  //    char pattern[]="%d %B %Y %I:%M:%S %p";
  //    tp.put(date_stream,                /*s*/
  //   date_stream,                /*str*/
  //   date_stream.fill(),         /*fill*/
  //   tm_struct,                  /*tm*/
  //   pattern,                    /*format begin*/
  //   pattern+sizeof(pattern)-1); /*format end  */

  return date_stream.str();
#else
  // C-style code originally found here:
  // http://people.sc.fsu.edu/~burkardt/cpp_src/timestamp/timestamp.C
  // Author: John Burkardt, 24 September 2003
  const unsigned int time_size = 40;
  char time_buffer[time_size];

  time_t now = time (NULL);
  tm * tm_struct = localtime (&now);

  // No more than time_size characters will be placed into the array.  If the
  // total number of resulting characters, including the terminating
  // NUL character, is not more than time_size, strftime() returns the
  // number of characters in the array, not counting the terminating
  // NUL.  Otherwise, zero is returned and the buffer contents are
  // indeterminate.
  size_t len = strftime ( time_buffer, time_size, "%c", tm_struct );

  if (len != 0)
    return std::string(time_buffer);
  else
    {
      libMesh::out << "Error formatting time buffer, returning empty string!" << std::endl;
      return std::string("");
    }

#endif // LIBMESH_HAVE_LOCALE
}

// std::time_put::put() documentation
// s=Iterator pointing to the first character of the output sequence.
//
// str=Object of a class derived from ios_base (generally an output
// stream object). It is used to obtain formatting information.
//
// fill=Fill character. The fill character is used in output insertion
// operations to fill spaces when the format requires some character padding.
//
// tm=Pointer to an object of type struct tm
//
// format=The final argument to time_put::put is an individual format character.
// The function will format some of the information pointed by the tm struct
// into a sequence of characters as specified by this character, just as if
// it was preceded by a percentage sign in a format string passed to strftime.
// (see 'man strftime' for more...)
// Example: 'c' is national representation of time and date
// 'c' = Thu Feb  4 12:24:11 2010
//  tp.put(date_stream /*s*/,
// date_stream /*str*/,
// date_stream.fill() /*fill*/,
// tm_struct, /*tm*/
// 'c'/*format*/);

// We can also pass to char * to the beginning and end of the desired format:
// const charT * pattern, const charT * pat_end.  This allows us to have the full
// flexibility of strftime.
// Example: "%d %B %Y %I:%M:%S %p"
// 04 February 2010 01:44:10 PM
}

} // namespace libMesh
