// $Id: enum_io_package.h,v 1.1 2004-07-14 19:26:06 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __enum_io_package_h__
#define __enum_io_package_h__

// C++ includes

// Local includes
#include "libmesh_config.h"




/*
 * The \p libMeshEnums namespace is the namespace all \p enum definitions
 * should be put into.
 */

// ------------------------------------------------------------
// enum IOPackage definition
namespace libMeshEnums {
  
  /**
   * libMesh interfaces with several different software packages
   * for the purposes of creating, reading, and writing mesh files.
   * These enumerations give an easy way of selecting one or the
   * other.
   */
  enum IOPackage
    { 
      TECPLOT,
      GMV,
      GMSH,
      VTK,
      DIVA,
      TETGEN,
      INVALID_IO_PACKAGE
    };
}

using namespace libMeshEnums;



#endif




