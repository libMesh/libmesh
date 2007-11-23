// $Id: number_lookups.h,v 1.2 2007-10-21 20:48:45 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson

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


#ifndef __number_lookups_h__
#define __number_lookups_h__

// Lookup tables for hierarchic numbering of basis functions
extern const unsigned char triangular_number_row[];
extern const unsigned char triangular_number_column[];
extern const unsigned char square_number_row[];
extern const unsigned char square_number_column[];
extern const unsigned char cube_number_row[];
extern const unsigned char cube_number_column[];
extern const unsigned char cube_number_page[];

#endif // ifndef __number_lookups_h__
