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



#ifndef LIBMESH_ENUM_ORDER_H
#define LIBMESH_ENUM_ORDER_H

// ------------------------------------------------------------
// enum Order definition
namespace libMesh {

/**
 * \enum libMesh::Order defines an \p enum for polynomial orders.
 * Fixing each label to a specific int, since \p InfFE and p refinement
 * may cast between them
 */
enum Order {CONSTANT     =  0,
            FIRST        =  1,
            SECOND       =  2,
            THIRD        =  3,
            FOURTH       =  4,
            FIFTH        =  5,
            SIXTH        =  6,
            SEVENTH      =  7,
            EIGHTH       =  8,
            NINTH        =  9,
            TENTH        = 10,

            ELEVENTH     = 11,
            TWELFTH      = 12,
            THIRTEENTH   = 13,
            FOURTEENTH   = 14,
            FIFTEENTH    = 15,
            SIXTEENTH    = 16,
            SEVENTEENTH  = 17,
            EIGHTTEENTH  = 18,
            NINTEENTH    = 19, // deprecated
            NINETEENTH   = 19,
            TWENTIETH    = 20,

            TWENTYFIRST   = 21,
            TWENTYSECOND  = 22,
            TWENTYTHIRD   = 23,
            TWENTYFOURTH  = 24,
            TWENTYFIFTH   = 25,
            TWENTYSIXTH   = 26,
            TWENTYSEVENTH = 27,
            TWENTYEIGHTH  = 28,
            TWENTYNINTH   = 29,
            THIRTIETH     = 30,

            THIRTYFIRST   = 31,
            THIRTYSECOND  = 32,
            THIRTYTHIRD   = 33,
            THIRTYFOURTH  = 34,
            THIRTYFIFTH   = 35,
            THIRTYSIXTH   = 36,
            THIRTYSEVENTH = 37,
            THIRTYEIGHTH  = 38,
            THIRTYNINTH   = 39,
            FORTIETH     = 40,

            FORTYFIRST   = 41,
            FORTYSECOND  = 42,
            FORTYTHIRD   = 43,

            INVALID_ORDER};

}

#endif // LIBMESH_ENUM_ORDER_H
