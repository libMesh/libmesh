// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_version.h"

#include <iostream>
#include <iomanip>

void libMesh::libmesh_version_stdout()
{
  std::cout << "--------------------------------------------------------" << std::endl;
  std::cout << "libMesh Library: Version = " << LIBMESH_LIB_VERSION;
  std::cout << " (" << get_libmesh_version() << ")" << std::endl << std::endl;

  std::cout << LIBMESH_LIB_RELEASE << std::endl << std::endl;

  std::cout << "Build Date   = " << LIBMESH_BUILD_DATE     << std::endl;
  std::cout << "Build Host   = " << LIBMESH_BUILD_HOST     << std::endl;
  std::cout << "Build User   = " << LIBMESH_BUILD_USER     << std::endl;
  std::cout << "Build Arch   = " << LIBMESH_BUILD_ARCH     << std::endl;
  std::cout << "Build Rev    = " << LIBMESH_BUILD_VERSION  << std::endl << std::endl;

  // CXXFLAGS is ambiguous wth multiple methods - could add all three but why not libmesh-config?
  //std::cout << "C++ Config   = " << LIBMESH_CXX << " " << LIBMESH_CXXFLAGS << std::endl;
  std::cout << "--------------------------------------------------------" << std::endl;

  return;
}



int libMesh::get_libmesh_version()
{
  /* Note: return format follows the versioning convention xx.yy.zz where

     xx = major version number
     yy = minor version number
     zz = micro version number

     For example:
     v.   0.23  -> 002300 = 2300
     v   0.23.1 -> 002301 = 2301
     v. 10.23.2 -> 102302         */

  int major_version = 0;
  int minor_version = 0;
  int micro_version = 0;

#ifdef LIBMESH_MAJOR_VERSION
  major_version = LIBMESH_MAJOR_VERSION;
#endif

#ifdef LIBMESH_MINOR_VERSION
  minor_version = LIBMESH_MINOR_VERSION;
#endif

#ifdef LIBMESH_MICRO_VERSION
  micro_version = LIBMESH_MICRO_VERSION;
#endif

  return major_version*10000 + minor_version*100 + micro_version;
}



std::string libMesh::get_io_compatibility_version ()
{
  std::string retval(LIBMESH_IO_COMPATIBILITY_VERSION);
  return retval;
}
