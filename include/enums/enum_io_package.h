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



#ifndef LIBMESH_ENUM_IO_PACKAGE_H
#define LIBMESH_ENUM_IO_PACKAGE_H


namespace libMesh {

/**
 * libMesh interfaces with several different software packages
 * for the purposes of creating, reading, and writing mesh files.
 * These enumerations give an easy way of selecting one or the
 * other.
 *
 * The fixed type, i.e. ": int", enumeration syntax used here allows
 * this enum to be forward declared as
 * enum IOPackage : int;
 * reducing header file dependencies.
 */
enum IOPackage : int
  {
    TECPLOT,
    GMV,
    GMSH,
    VTK,
    DIVA,
    TETGEN,
    UCD,
    LIBMESH,
    // Invalid
    INVALID_IO_PACKAGE
  };
}

#endif
