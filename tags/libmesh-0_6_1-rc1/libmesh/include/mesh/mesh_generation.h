// $Id: mesh_generation.h,v 1.14 2007-09-25 19:54:42 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include <vector>

// Local Includes -----------------------------------
// #include "libmesh_common.h" // needed for Real
#include "enum_elem_type.h" // needed for ElemType enum
#include "libmesh.h"        // needed for libMesh::invalid_uint
#include "point.h"

// forward declarations
class UnstructuredMesh;



// ------------------------------------------------------------
// MeshTools::Generation namespace
namespace MeshTools
{
  /**
   * Tools for \p Mesh generation.
   *
   * \author Benjamin S. Kirk
   * \date 2004
   * \version $Revision: 1.14 $
   */
  namespace Generation
  {
    /**
     * Builds a \f$ nx \times ny \times nz \f$ (elements) cube.
     * Defaults to a unit cube (or line in 1D, square in 2D),
     * but the dimensions can be specified through the optional
     * arguments.
     */  
    void build_cube (UnstructuredMesh& mesh,
		     const unsigned int nx=0,
		     const unsigned int ny=0,
		     const unsigned int nz=0,
		     const Real xmin=0., const Real xmax=1.,
		     const Real ymin=0., const Real ymax=1.,
		     const Real zmin=0., const Real zmax=1.,
		     const ElemType type=INVALID_ELEM,
		     const bool gauss_lobatto_grid=false);

    /**
     * A specialized \p build_cube() for 1D meshes
     */
    void build_line (UnstructuredMesh& mesh,
                     const unsigned int nx,
                     const Real xmin=0., const Real xmax=1.,
                     const ElemType type=INVALID_ELEM,
                     const bool gauss_lobatto_grid=false);

    /**
     * A specialized \p build_cube() for 2D meshes.
     */
    void build_square (UnstructuredMesh& mesh,
		       const unsigned int nx,
		       const unsigned int ny,
		       const Real xmin=0., const Real xmax=1.,
		       const Real ymin=0., const Real ymax=1.,
		       const ElemType type=INVALID_ELEM,
		       const bool gauss_lobatto_grid=false);

    /**
     * Meshes a spherical or mapped-spherical domain.
     */
    void build_sphere (UnstructuredMesh& mesh,
		       const Real rad=1,
		       const unsigned int nr=2,
		       const ElemType type=INVALID_ELEM);

    
    namespace Private
    {
      /**
       * A useful inline function which replaces the #defines
       * used previously.  Not private since this is a namespace,
       * but would be if this were a class.  The first one returns
       * the proper node number for 2D elements while the second
       * one returns the node number for 3D elements.
       */
      inline
      unsigned int idx(const ElemType type,
		       const unsigned int nx,
		       const unsigned int i,
		       const unsigned int j)
      {
	switch(type)
	  {
	  case INVALID_ELEM:
	  case QUAD4:
	  case TRI3:
	    {
	      return i + j*(nx+1);
	      break;
	    }

	  case QUAD8:
	  case QUAD9:
	  case TRI6:
	    {
	      return i + j*(2*nx+1);
	      break;
	    }
	  
	  default:
	    {
	      std::cerr << "ERROR: Unrecognized 2D element type." << std::endl;
	      error();
	    }
	  }

	return libMesh::invalid_uint;
      }


    
      // Same as the function above, but for 3D elements
      inline
      unsigned int idx(const ElemType type,
		       const unsigned int nx,
		       const unsigned int ny,
		       const unsigned int i,
		       const unsigned int j,
		       const unsigned int k)
      {
	switch(type)
	  {
	  case INVALID_ELEM:
	  case HEX8:
	  case PRISM6:
	    {
	      return i + (nx+1)*(j + k*(ny+1));
	      break;
	    }

	  case HEX20:
	  case HEX27:
	  case TET4:  // TET4's are created from an initial HEX27 discretization
	  case TET10: // TET10's are created from an initial HEX27 discretization
	  case PRISM15:
	  case PRISM18:
	    {
	      return i + (2*nx+1)*(j + k*(2*ny+1));
	      break;
	    }
	  
	  default:
	    {
	      std::cerr << "ERROR: Unrecognized element type." << std::endl;
	      error();
	    }
	  }
      
	return libMesh::invalid_uint;
      }
    } // end namespace Meshtools::Generation::Private
  } // end namespace Meshtools::Generation
} // end namespace MeshTools


#endif // #define __mesh_generation_h__
