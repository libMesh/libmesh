// $Id: inf_fe_macro.h,v 1.1 2003-01-24 17:24:38 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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


#ifndef __inf_fe_macro_h__
#define __inf_fe_macro_h__


// Local includes
#include "mesh_config.h"

#ifdef ENABLE_INFINITE_ELEMENTS



/**
 * Include this file when you either 
 * (a) want to instantiate the template classes 
 *     \p InfFE<Dim,T_radial,T_map>, or
 * (b) declare these classes as friends.  
 *
 * The classes are instantiated in \p inf_fe_instantiate_1D.h, 
 * \p inf_fe_instantiate_2D.h, and \p inf_fe_instantiate_3D.h
 * for 1D, 2D and 3D, respectively.
 */
#define SAY_THAT_INF_FE_IS(_say_what,_dim,_map_type) _say_what  class InfFE< _dim, INFINITE_MAP, _map_type >; \
                                                     _say_what  class InfFE< _dim, JACOBI_20_00, _map_type >; \
                                                     _say_what  class InfFE< _dim, JACOBI_30_00, _map_type >; \
                                                     _say_what  class InfFE< _dim, LEGENDRE,     _map_type >; \
                                                     _say_what  class InfFE< _dim, INF_LAGRANGE, _map_type >; 



#else



#define SAY_THAT_INF_FE_IS(_say_what,_dim,_map_type) \



#endif //ifdef ENABLE_INFINITE_ELEMENTS


#endif
