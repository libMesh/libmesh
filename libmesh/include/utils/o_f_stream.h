// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include "libmesh_common.h"



/*
 * Some compilers, at least HP \p aCC do not even
 * accept empty classes derived from \p std::ostringstream.
 * Therefore, resort to preprocessor definitions. 
 */

#ifndef BROKEN_IOSTREAM
 /*
  * ---------------------------------------------------------------------------------
  * Everything for a clean iostream
  */

# include <iomanip>

 /*
  * Outputs \p Real \p d  with width
  * \p v in scientific format to stream \p o.
  */
# define OFSRealscientific(o,v,d)       (o) << std::setprecision(v) << std::scientific << (d)

 /*
  * Outputs \p Number \p d, (note that \p Number
  * may be either real or complex) with width
  * \p v in scientific format to stream \p o.
  */
# if defined(USE_COMPLEX_NUMBERS) 
#  define OFSNumberscientific(o,v,d)    (o) << std::setprecision(v) << std::scientific << (d).real() << " " \
                                            << std::setprecision(v) << std::scientific << (d).imag()
# else
#  define OFSNumberscientific(o,v,d)    (o) << std::setprecision(v) << std::scientific << (d)
# endif

 /*
  * class alias
  */
# define OFStream                    std::ofstream


#else
 /*
  * ---------------------------------------------------------------------------------
  * Everything for broken iostream
  */

# include <cstdio>

 /*
  * Outputs \p Real \p d with width 
  * \p v in scientific format to stream \p o.
  */
# define OFSRealscientific(o,v,d)       (o).scientific( (v), (d) )

 /*
  * Outputs \p Number \p d, (note that \p Number
  * may be either real or complex) with width
  * \p v in scientific format to stream \p o.
  */
# define OFSNumberscientific(o,v,d)     (o).scientific( (v), (d) )

//  /*
//   * class alias
//   */
// # define OFSOFStream                    OFStream



 /**
  * This class provides a compatibility class for broken
  * features in the \p std::ofstream of the older \p GCC
  * versions.  Other compilers do not see this class.
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
  template <typename T>
  OFStream& scientific (const sizetype w,
			const T r);

 };


 // ------------------------------------------------------------
 // OFStream inline methods
 

#endif // ifndef ... else ... BROKEN_IOSTREAM



#endif // ifndef __o_f_stream_h__

