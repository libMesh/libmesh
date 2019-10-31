// The PassMess Message-Passing Parallelism Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "passmess/passmess_version.h"

namespace PassMess
{

  void passmess_version_stdout()
  {
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "PassMess Package: Version = " << PASSMESS_LIB_VERSION;
    std::cout << " (" << get_passmess_version() << ")" << std::endl << std::endl;
  
    std::cout << PASSMESS_LIB_RELEASE << std::endl << std::endl;
  
    std::cout << "Build Date   = " << PASSMESS_BUILD_DATE     << std::endl;
    std::cout << "Build Host   = " << PASSMESS_BUILD_HOST     << std::endl;
    std::cout << "Build User   = " << PASSMESS_BUILD_USER     << std::endl;
    std::cout << "Build Arch   = " << PASSMESS_BUILD_ARCH     << std::endl;
    std::cout << "Build Rev    = " << PASSMESS_BUILD_VERSION  << std::endl << std::endl;
  
    std::cout << "C++ Config   = " << PASSMESS_CXX << " " << PASSMESS_CXXFLAGS << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
  
    return;
  }

  int get_passmess_version()
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

#ifdef PASSMESS_MAJOR_VERSION
    major_version = PASSMESS_MAJOR_VERSION;
#endif

#ifdef PASSMESS_MINOR_VERSION
    minor_version = PASSMESS_MINOR_VERSION;
#endif

#ifdef PASSMESS_MICRO_VERSION
    micro_version = PASSMESS_MICRO_VERSION;
#endif
      
    return major_version*10000 + minor_version*100 + micro_version;
  }

} // end namespace PassMess
