// $Id: mesh_modification.h,v 1.4 2004-12-17 20:55:07 benkirk Exp $

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



#ifndef __mesh_modification_h__
#define __mesh_modification_h__



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "libmesh_common.h"

// forward declarations
class MeshBase;



// ------------------------------------------------------------
// MeshTools::Modification namespace
namespace MeshTools
{
  /**
   * Tools for \p Mesh modification.
   *
   * \author Benjamin S. Kirk
   * \date 2004
   * \version $Revision: 1.4 $
   */  
  namespace Modification
  {
    /**
     * Randomly perturb the nodal locations.  This function will
     * move each node \p factor fraction of its minimum neighboring
     * node separation distance.  Nodes on the boundary are not moved
     * by default, however they may be by setting the flag
     * \p perturb_boundary true.
     */
    void distort (MeshBase& mesh,
		  const Real factor, const bool perturb_boundary=false);
    
    /**
     * Translates the mesh.  The grid points are translated in the
     * \p x direction by \p xt, in the \p y direction by \p yt,
     * etc...
     */
    void translate (MeshBase& mesh,
		    const Real xt=0., const Real yt=0., const Real zt=0.); 

//     /**
//      * Rotates the mesh in the xy plane. The rotation is
//      * counter-clock-wise (mathematical definition).
//      * The angle is in degrees (360 make a full circle)
//      */
//     void rotate2D (MeshBase& mesh,
//                    const Real alpha=0.); 

    /**
     * Rotates the mesh in 3D space. 
     * Here the standard Euler angles are adopted
     * (http://mathworld.wolfram.com/EulerAngles.html)
     * The angles are in degrees (360 make a full circle)
     */
    void rotate (MeshBase& mesh,
		 const Real phi, const Real theta=0., const Real psi=0.); 

    /**
     * Scales the mesh.  The grid points are scaled in the
     * \p x direction by \p xs, in the \p y direction by \p ys,
     * etc...  If only \p xs is specified then the scaling is
     * assumed uniform in all directions.
     */
    void scale (MeshBase& mesh,
		const Real xs, const Real ys=0., const Real zs=0.);

    /**
     * Converts the 2D quadrilateral elements of a Mesh into
     * triangular elements.
     * Note: Only works for 2D elements!  3D elements are ignored.
     * Note: Probably won't do the right thing for meshes which
     * have been refined previously.
     */
    void all_tri (MeshBase& mesh);
    
  } // end namespace Meshtools::Modification
} // end namespace MeshTools


#endif // #define __mesh_modification_h__
