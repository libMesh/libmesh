
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



#ifndef __location_maps_h__
#define __location_maps_h__

#include "libmesh_config.h"

// C++ Includes   -----------------------------------
#if   defined(LIBMESH_HAVE_UNORDERED_MAP)
# include <unordered_map>
#elif defined(LIBMESH_HAVE_TR1_UNORDERED_MAP)
# include <tr1/unordered_map>
#elif defined(LIBMESH_HAVE_HASH_MAP)
# include <hash_map>
#elif defined(LIBMESH_HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#else
# include <map>
#endif
#include <vector>

// Local Includes -----------------------------------
#include "libmesh_common.h"
#include "point.h"

namespace libMesh
{

// Forward Declarations -----------------------------
class Elem;
class MeshBase;
class Node;


  /**
   * Data structures that enable location-based lookups
   * The key is a hash of the Point location.
   * For efficiency we will use a hashed multimap if it is
   * available, otherwise a regular multimap.
   */
template <typename T>
class LocationMap
{
#if   defined(LIBMESH_HAVE_UNORDERED_MAP)
  typedef std::unordered_multimap<unsigned int, T*> map_type;
#elif defined(LIBMESH_HAVE_TR1_UNORDERED_MAP)
  typedef std::tr1::unordered_multimap<unsigned int, T*> map_type;
#elif defined(LIBMESH_HAVE_HASH_MAP)
  typedef std::hash_multimap<unsigned int, T*> map_type;
#elif defined(LIBMESH_HAVE_EXT_HASH_MAP)
# if   (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
  typedef std::hash_multimap<unsigned int, T*> map_type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
  typedef __gnu_cxx::hash_multimap<unsigned int, T*> map_type;
# else
// XLC and who knows what other compilers get here.
// Try the most standard thing we can:
  typedef std::multimap<unsigned int, T*> map_type;
# endif
#else
  typedef std::multimap<unsigned int, T*> map_type;
#endif
public: 
  void init(MeshBase&);

  void clear() { _map.clear(); }

  void insert(T&);

  bool empty() const { return _map.empty(); }

  T* find(const Point&,
	  const Real tol = TOLERANCE);

  Point point_of(const T&) const;

protected:
  unsigned int key(const Point&);

  void fill(MeshBase&);

private: 
  map_type          _map;
  std::vector<Real> _lower_bound;
  std::vector<Real> _upper_bound;
};

} // namespace libMesh


#endif // #define __location_maps_h__
