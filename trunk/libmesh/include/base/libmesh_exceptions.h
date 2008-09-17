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



#ifndef __libmesh_exceptions_h__
#define __libmesh_exceptions_h__

#include "libmesh_config.h"

#ifdef ENABLE_EXCEPTIONS
#include <stdexcept>


namespace libMesh {

  /**
   * A class to represent the internal "this should never happen"
   * errors to be thrown by "libmesh_error();"
   */
  class LogicError : public std::logic_error
    {
    public:
      LogicError() : std::logic_error( "Error in libMesh internal logic" ) {}
    };


  /**
   *
   */
  class NotImplemented : public std::logic_error
    {
    public:
      NotImplemented() : std::logic_error( "Error: not implemented!" ) {}
    };

}

#define LIBMESH_THROW(e) do { throw e; } while (0)

#else

#define LIBMESH_THROW(e) do { std::abort(); } while (0)

#endif // ENABLE_EXCEPTIONS

#endif // #define __libmesh_exceptions_h__
