// $Id: mesh_generation.h,v 1.1 2004-11-15 22:09:12 benkirk Exp $

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



#ifndef __mesh_generation_h__
#define __mesh_generation_h__



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "libmesh_common.h"
#include "enum_elem_type.h"

// forward declarations
class Mesh;



/**
 * Tools for \p Mesh generation.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 * \version $Revision: 1.1 $
 */


// ------------------------------------------------------------
// MeshTools::Generation namespace
namespace MeshTools
{
  namespace Generation
  {
    /**
     * Builds a \f$ nx \times ny \times nz \f$ (elements) cube.
     * Defaults to a unit cube (or line in 1D, square in 2D),
     * but the dimensions can be specified through the optional
     * arguments.
     */  
    void build_cube (Mesh& mesh,
		     const unsigned int nx=0,
		     const unsigned int ny=0,
		     const unsigned int nz=0,
		     const Real xmin=0., const Real xmax=1.,
		     const Real ymin=0., const Real ymax=1.,
		     const Real zmin=0., const Real zmax=1.,
		     const ElemType type=INVALID_ELEM);

    /**
     * A specialized \p build_cube() for 2D meshes.
     */
    void build_square (Mesh& mesh,
		       const unsigned int nx,
		       const unsigned int ny,
		       const Real xmin=0., const Real xmax=1.,
		       const Real ymin=0., const Real ymax=1.,
		       const ElemType type=INVALID_ELEM);

    /**
     * Meshes a spherical or mapped-spherical domain.
     */
    void build_sphere (Mesh& mesh,
		       const Real rad=1,
		       const unsigned int nr=2,
		       const ElemType type=INVALID_ELEM);
    
  } // end namespace Meshtools::Generation
} // end namespace MeshTools


#endif // #define __mesh_generation_h__
