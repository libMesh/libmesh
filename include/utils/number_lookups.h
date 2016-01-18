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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA


#ifndef LIBMESH_NUMBER_LOOKUPS_H
#define LIBMESH_NUMBER_LOOKUPS_H

namespace libMesh
{

// Lookup tables for hierarchic numbering of basis functions
extern const unsigned char triangular_number_row[];
extern const unsigned char triangular_number_column[];
extern const unsigned char square_number_row[];
extern const unsigned char square_number_column[];
extern const unsigned char cube_number_row[];
extern const unsigned char cube_number_column[];
extern const unsigned char cube_number_page[];

} // namespace libMesh

#endif // LIBMESH_NUMBER_LOOKUPS_H
