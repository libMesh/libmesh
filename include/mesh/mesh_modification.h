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



#ifndef LIBMESH_MESH_MODIFICATION_H
#define LIBMESH_MESH_MODIFICATION_H



// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/id_types.h" // for boundary_id_type, subdomain_id_type

// C++ Includes   -----------------------------------

namespace libMesh
{

// forward declarations
template <typename Output> class FunctionBase;
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
void distort (MeshBase & mesh,
              const Real factor, const bool perturb_boundary=false);

/**
 * Deterministically perturb the nodal locations.  This function will
 * move each node from it's current x/y/z coordinates to a new x/y/z
 * coordinate given by the first LIBMESH_DIM components of the
 * specified function \p mapfunc
 *
 * Nodes on the boundary are also moved.
 *
 * Currently, non-vertex nodes are moved in the same way as vertex
 * nodes, according to (newx,newy,newz) = mapfunc(x,y,z).  This
 * behavior is often suboptimal for higher order geometries and may be
 * subject to change in future libMesh versions.
 */
void redistribute (MeshBase & mesh,
                   const FunctionBase<Real> & mapfunc);


/**
 * Translates the mesh.  The grid points are translated in the
 * \p x direction by \p xt, in the \p y direction by \p yt,
 * etc...
 */
void translate (MeshBase & mesh,
                const Real xt=0., const Real yt=0., const Real zt=0.);

//     /**
//      * Rotates the mesh in the xy plane. The rotation is
//      * counter-clock-wise (mathematical definition).
//      * The angle is in degrees (360 make a full circle)
//      */
//     void rotate2D (MeshBase & mesh,
//                    const Real alpha=0.);

/**
 * Rotates the mesh in 3D space.
 * Here the standard Euler angles are adopted
 * (http://mathworld.wolfram.com/EulerAngles.html)
 * The angles are in degrees (360 make a full circle)
 */
void rotate (MeshBase & mesh,
             const Real phi, const Real theta=0., const Real psi=0.);

/**
 * Scales the mesh.  The grid points are scaled in the
 * \p x direction by \p xs, in the \p y direction by \p ys,
 * etc...  If only \p xs is specified then the scaling is
 * assumed uniform in all directions.
 */
void scale (MeshBase & mesh,
            const Real xs, const Real ys=0., const Real zs=0.);

/**
 * Converts the 2D quadrilateral elements of a Mesh into
 * triangular elements.
 * Note: Only works for 2D elements!  3D elements are ignored.
 * Note: Probably won't do the right thing for meshes which
 * have been refined previously.
 */
void all_tri (MeshBase & mesh);

/**
 * Smooth the mesh with a simple Laplace smoothing algorithm.  The mesh is
 * smoothed \p n_iterations times.  If the parameter \p power is 0, each
 * node is moved to the average postition of the neighboring connected
 * nodes. If \p power > 0, the node positions are weighted by their
 * distance.  The positions of higher order nodes, and nodes living in
 * refined elements, are calculated from the vertex positions of their
 * parent nodes.  Only works in 2D.
 *
 * \author Martin Luthi (luthi@gi.alaska.edu)
 * \date 2005
 */
void smooth(MeshBase &, unsigned int, Real);

#ifdef LIBMESH_ENABLE_AMR
/**
 * Removes all the refinement tree structure of Mesh, leaving
 * only the highest-level (most-refined) elements.  This is useful
 * when you want to write out a uniformly-refined grid to be treated later
 * as an initial mesh.  Note that many functions in LibMesh assume a
 * conforming (with no hanging nodes) grid exists at some level, so
 * you probably only want to do this on meshes which have been uniformly
 * refined.
 */
void flatten(MeshBase & mesh);
#endif // #ifdef LIBMESH_ENABLE_AMR

/**
 * Finds any boundary ids that are currently old_id,
 * changes them to new_id
 */
void change_boundary_id (MeshBase & mesh,
                         const boundary_id_type old_id,
                         const boundary_id_type new_id);

/**
 * Finds any subdomain ids that are currently old_id,
 * changes them to new_id
 */
void change_subdomain_id (MeshBase & mesh,
                          const subdomain_id_type old_id,
                          const subdomain_id_type new_id);

} // end namespace Meshtools::Modification
} // end namespace MeshTools

} // namespace libMesh


#endif // LIBMESH_MESH_MODIFICATION_H
