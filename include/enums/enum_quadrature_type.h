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



#ifndef LIBMESH_ENUM_QUADRATURE_TYPE_H
#define LIBMESH_ENUM_QUADRATURE_TYPE_H

// ------------------------------------------------------------
// enum QuadratureType definition
namespace libMesh {

/**
 * Defines an \p enum for currently available quadrature rules.
 */
enum QuadratureType {QGAUSS            = 0,

                     QJACOBI_1_0       = 1,
                     QJACOBI_2_0       = 2,

                     QSIMPSON          = 3,
                     QTRAP             = 4,
                     QGRID             = 5,
                     QGRUNDMANN_MOLLER = 6,
                     QMONOMIAL         = 7,
                     QCONICAL          = 8,
                     QGAUSS_LOBATTO    = 9,

                     QCLOUGH           = 21,

                     QCOMPOSITE        = 31,

                     INVALID_Q_RULE    = 127};
}

#endif // LIBMESH_ENUM_QUADRATURE_TYPE_H
