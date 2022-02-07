// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FE_INTERFACE_MACROS_H
#define LIBMESH_FE_INTERFACE_MACROS_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#include "libmesh/enum_to_string.h"

#define inf_fe_switch(func_and_args)                                \
  do {                                                              \
    switch (dim)                                                    \
      {                                                             \
      case 1:                                                       \
        inf_fe_family_mapping_switch (1, func_and_args, return, ;); \
      case 2:                                                       \
        inf_fe_family_mapping_switch (2, func_and_args, return, ;); \
      case 3:                                                       \
        inf_fe_family_mapping_switch (3, func_and_args, return, ;); \
      default:                                                      \
        libmesh_error_msg("Invalid dim = " << dim);                 \
      }                                                             \
  } while (0)

#define inf_fe_family_mapping_switch(dim, func_and_args, prefix, suffix)    \
   do{                                                                      \
     switch(fe_t.inf_map)                                                   \
       {                                                                    \
       case CARTESIAN:                                                      \
         {                                                                  \
          switch (fe_t.radial_family)                                       \
            {                                                               \
            case INFINITE_MAP:                                              \
              prefix InfFE<dim,INFINITE_MAP,CARTESIAN>::func_and_args suffix\
            case JACOBI_20_00:                                              \
              prefix InfFE<dim,JACOBI_20_00,CARTESIAN>::func_and_args suffix\
            case JACOBI_30_00:                                              \
              prefix InfFE<dim,JACOBI_30_00,CARTESIAN>::func_and_args suffix\
            case LEGENDRE:                                                  \
              prefix InfFE<dim,LEGENDRE,CARTESIAN>::func_and_args suffix    \
            case LAGRANGE:                                                  \
              prefix InfFE<dim,LAGRANGE,CARTESIAN>::func_and_args suffix    \
            default:                                                        \
              libmesh_error_msg("Invalid radial family = " << Utility::enum_to_string(fe_t.radial_family)); \
            }                                                               \
            suffix                                                          \
         }                                                                  \
       case SPHERICAL:                                                      \
       case ELLIPSOIDAL:                                                    \
         libmesh_not_implemented();                                         \
       default:                                                             \
         libmesh_error_msg("Invalid radial mapping " << Utility::enum_to_string(fe_t.inf_map)); \
      }                                                                     \
  } while (0)

#endif //LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // define LIBMESH_FE_INTERFACE_MACROS_H
