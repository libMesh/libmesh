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

#ifdef LIBMESH_ENABLE_EXCEPTIONS
#include <stdexcept>


namespace libMesh {

  /**
   * A class to represent the internal "this should never happen"
   * errors; to be thrown by "libmesh_error();"
   */
  class LogicError : public std::logic_error
    {
    public:
      LogicError() : std::logic_error( "Error in libMesh internal logic" ) {}
    };


  /**
   * A class to stub for features that should be in libMesh, but
   * haven't been written yet; to be thrown by
   * "libmesh_not_implemented();"
   */
  class NotImplemented : public std::logic_error
    {
    public:
      NotImplemented() : std::logic_error( "Error: not implemented!" ) {}
    };


  /**
   * A class representing a failed attempt by the library to open a
   * file (or construct an fstream, etc); to be thrown by
   * "libmesh_file_error(filename);" For ease of debugging, "filename"
   * should include any (absolute or relative or implicit) pathname
   * that was part of the failed open.
   */
  class FileError : public std::runtime_error
    {
    public:
      FileError(const std::string& filename) : std::runtime_error( "Error accessing file: " + filename ) {}
    };

}

#define LIBMESH_THROW(e) do { throw e; } while (0)

#else

#define LIBMESH_THROW(e) do { std::abort(); } while (0)

#endif // LIBMESH_ENABLE_EXCEPTIONS

#endif // #define __libmesh_exceptions_h__
