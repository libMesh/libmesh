// $Id: o_f_stream.C,v 1.1 2003-03-22 21:04:31 ddreyer Exp $

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



// C++ includes


// Local Includes
#include "o_f_stream.h"



// the template for reals
template <typename T>
OFStream& OFStream::scientific (const sizetype w,
				const std::vector<T>& v,
				const char* sep)
{
#ifndef BROKEN_IOSTREAM
  typename std::vector<T>::const_iterator pos = v.begin();
  for (; pos != v.end(); ++pos)
      *this << std::setw(w)
	    << std::scientific
	    << *pos 
	    << sep;
#else
  assert (w < 30);
  char buf[30];  
  char format[8];
  // form the format for r
  sprintf (format, "%%%de", w);
  typename std::vector<T>::const_iterator pos = v.begin();
  for (; pos != v.end(); ++pos)
    {
      sprintf (buf, format, *pos);
      *this << buf << sep;
    }
#endif
  return *this;
}





// full specialization in case of complex numbers
#if defined(USE_COMPLEX_NUMBERS) 

template <>
OFStream& OFStream::scientific (const sizetype w,
				const std::vector<Complex>& v,
				const char* sep)
{
#ifndef BROKEN_IOSTREAM
  std::vector<Complex>::const_iterator pos = v.begin();
  for (; pos != v.end(); ++pos)
      *this << std::setw(w)
	    << std::scientific
	    << (*pos).real()
	    << sep
	    << std::setw(w)
	    << std::scientific
	    << (*pos).imag()
	    << sep;
#else
  assert (w < 30);
  char buf[60];  
  char format[20];
  // form the format for r
  sprintf (format, "%%%de%s%%%de%s", w, sep, w, sep);
  std::vector<Complex>::const_iterator pos = v.begin();
  unsigned int cnt=0;
  for (; pos != v.end(); ++pos)
    {
      sprintf (buf, format, (*pos).real(), (*pos).imag());
      *this << buf;
    }
#endif
  return *this;
}

#endif // if defined(USE_COMPLEX_NUMBERS) 





//--------------------------------------------------------------
// Explicit instantiations for reals
template OFStream& OFStream::scientific (const sizetype w,
					 const std::vector<double>& v,
					 const char* sep);

template OFStream& OFStream::scientific (const sizetype w,
					 const std::vector<float>& v,
					 const char* sep);
