// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_INF_FE_MACRO_H
#define LIBMESH_INF_FE_MACRO_H


// Local includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS



/**
 * This macro helps in instantiating specific versions
 * of the \p InfFE class.  Better do not use this macro
 * directly, but instantiate through #include`ing the
 * file(s) \p inf_fe_instantiate_1D.h, \p inf_fe_instantiate_2D.h,
 * and \p inf_fe_instantiate_3D.h for 1D, 2D and 3D, respectively.
 */
#define INSTANTIATE_INF_FE(_dim,_map_type) template  class InfFE< _dim, INFINITE_MAP, _map_type >; \
  template  class InfFE< _dim, JACOBI_20_00, _map_type >;               \
  template  class InfFE< _dim, JACOBI_30_00, _map_type >;               \
  template  class InfFE< _dim, LEGENDRE,     _map_type >;               \
  template  class InfFE< _dim, LAGRANGE,     _map_type >

#define INSTANTIATE_INF_FE_MBRF(_dim,_map_type,_return,_function)       \
  template _return InfFE< _dim,INFINITE_MAP,_map_type>::_function;      \
  template _return InfFE< _dim,JACOBI_20_00,_map_type>::_function;      \
  template _return InfFE< _dim,JACOBI_30_00,_map_type>::_function;      \
  template _return InfFE< _dim,LEGENDRE,_map_type>::_function;          \
  template _return InfFE< _dim,LAGRANGE,_map_type>::_function

#else



#define INSTANTIATE_INF_FE(_dim,_map_type)      \



#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_INF_FE_MACRO_H
