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



// C++ includes


// Local Includes
#include "o_f_stream.h"

/*
 * This class is only alive when iostream is broken
 */
#ifdef LIBMESH_BROKEN_IOSTREAM


namespace libMesh
{


// the template for reals
template <typename T>
OFStream& OFStream::scientific (const sizetype w,
				const T r)
{
  libmesh_assert (w < 30);
  char buf[30];  
  char format[8];
  // form the format for r
  sprintf (format, "%%%de", w);
  // form string as desired
  sprintf (buf, format, r);
  *this << buf;
  return *this;
}



// full specialization in case of complex numbers
#if defined(LIBMESH_USE_COMPLEX_NUMBERS) 

template <>
OFStream& OFStream::scientific (const sizetype w,
				const Complex r)
{
  libmesh_assert (w < 30);
  char buf[60];  
  char format[16];
  // form the format for r
  sprintf (format, "%%%de %%%de", w, w);
  // form string as desired
  sprintf (buf, format, r.real(), r.imag());
  *this << buf;
  return *this;
}

#endif // if defined(LIBMESH_USE_COMPLEX_NUMBERS) 



//--------------------------------------------------------------
// Explicit instantiations for reals
template OFStream& OFStream::scientific (const sizetype w,
					 const double r);

template OFStream& OFStream::scientific (const sizetype w,
					 const float r);

} // namespace libMesh



#endif // ifdef LIBMESH_BROKEN_IOSTREAM
