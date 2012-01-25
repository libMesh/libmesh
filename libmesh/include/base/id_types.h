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



#ifndef __id_types_h__
#define __id_types_h__

#include <limits>
#include "libmesh_config.h"
#include <stdint.h>

namespace libMesh
{

// A useful way to debug:
#if 0
class TestClass {
//int _c;
unsigned int _c;
public:
TestClass() : _c(0) {}
TestClass(unsigned int c) : _c(c) {}
TestClass& operator=(unsigned int c) { _c = c; return *this; }
bool operator<(const TestClass &l) const { return _c < l._c; }
operator int() const { return _c; }
};
typedef TestClass subdomain_id_type;
#endif

#if LIBMESH_SUBDOMAIN_ID_BYTES == 1
typedef uint8_t subdomain_id_type;
#elif LIBMESH_SUBDOMAIN_ID_BYTES == 4
/**
 * Note: subdomain_id_types are positive integers - however limitation in the exodusII
 * API force us to use a signed integer here to represent subdomains.  This gives us 2^31
 * possible unique blocks
 */
typedef int32_t subdomain_id_type;
#else // LIBMESH_SUBDOMAIN_ID_BYTES = 2 (default)
typedef uint16_t subdomain_id_type;
#endif

} // namespace libMesh

#endif // end #ifndef __id_types_h__
