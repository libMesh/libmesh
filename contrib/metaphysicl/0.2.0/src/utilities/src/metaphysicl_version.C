//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// MetaPhysicL - A metaprogramming library for physics calculations
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "metaphysicl/metaphysicl_version.h"

namespace MetaPhysicL
{

  void metaphysicl_version_stdout()
  {
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "metaphysicl Package: Version = " << METAPHYSICL_LIB_VERSION;
    std::cout << " (" << get_metaphysicl_version() << ")" << std::endl << std::endl;
  
    std::cout << METAPHYSICL_LIB_RELEASE << std::endl << std::endl;
  
    std::cout << "Build Date   = " << METAPHYSICL_BUILD_DATE     << std::endl;
    std::cout << "Build Host   = " << METAPHYSICL_BUILD_HOST     << std::endl;
    std::cout << "Build User   = " << METAPHYSICL_BUILD_USER     << std::endl;
    std::cout << "Build Arch   = " << METAPHYSICL_BUILD_ARCH     << std::endl;
    std::cout << "Build Rev    = " << METAPHYSICL_BUILD_VERSION  << std::endl << std::endl;
  
    std::cout << "C++ Config   = " << METAPHYSICL_CXX << " " << METAPHYSICL_CXXFLAGS << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
  
    return;
  }

  int get_metaphysicl_version()
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

#ifdef METAPHYSICL_MAJOR_VERSION
    major_version = METAPHYSICL_MAJOR_VERSION;
#endif

#ifdef METAPHYSICL_MINOR_VERSION
    minor_version = METAPHYSICL_MINOR_VERSION;
#endif

#ifdef METAPHYSICL_MICRO_VERSION
    micro_version = METAPHYSICL_MICRO_VERSION;
#endif
      
    return major_version*10000 + minor_version*100 + micro_version;
  }

} // end namespace MetaPhysicL
