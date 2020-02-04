// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MESH_GENERATION_H
#define LIBMESH_MESH_GENERATION_H

// Local Includes
#include "libmesh/libmesh.h"
#include "libmesh/enum_elem_type.h" // INVALID_ELEM
#include "libmesh/vector_value.h"
#ifdef LIBMESH_HAVE_TRIANGLE
#include "libmesh/mesh_triangle_interface.h"
#endif

// C++ Includes
#include <cstddef>
#include <vector>

namespace libMesh
{

// forward declarations
template <typename> class MeshBaseTempl;
template <typename> class UnstructuredMeshTempl;
template <typename> class ElemTempl;

// ------------------------------------------------------------
// MeshTools::Generation namespace
namespace MeshTools
{
/**
 * Tools for \p Mesh generation.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
namespace Generation
{

// forward declaration
template <typename> class QueryElemSubdomainIDBaseTempl;

/**
 * Builds a \f$ nx \times ny \times nz \f$ (elements) cube.
 * Defaults to a unit cube (or line in 1D, square in 2D),
 * but the dimensions can be specified through the optional
 * arguments.
 *
 * Boundary ids are set to be equal to the side indexing on a
 * master hex
 */
template <typename RealType>
void build_cube (UnstructuredMeshTempl<RealType> & mesh,
                 const unsigned int nx=0,
                 const unsigned int ny=0,
                 const unsigned int nz=0,
                 const Real xmin=0., const Real xmax=1.,
                 const Real ymin=0., const Real ymax=1.,
                 const Real zmin=0., const Real zmax=1.,
                 const ElemType type=INVALID_ELEM,
                 const bool gauss_lobatto_grid=false);

/**
 * A specialized \p build_cube() for 0D meshes.  The resulting
 * mesh is a single NodeElem suitable for ODE tests
 */
template <typename RealType>
void build_point (UnstructuredMeshTempl<RealType> & mesh,
                  const ElemType type=INVALID_ELEM,
                  const bool gauss_lobatto_grid=false);

/**
 * A specialized \p build_cube() for 1D meshes
 *
 * Boundary ids are set to be equal to the side indexing on a
 * master edge
 */
template <typename RealType>
void build_line (UnstructuredMeshTempl<RealType> & mesh,
                 const unsigned int nx,
                 const Real xmin=0., const Real xmax=1.,
                 const ElemType type=INVALID_ELEM,
                 const bool gauss_lobatto_grid=false);

/**
 * A specialized \p build_cube() for 2D meshes.
 *
 * Boundary ids are set to be equal to the side indexing on a
 * master quad
 */
template <typename RealType>
void build_square (UnstructuredMeshTempl<RealType> & mesh,
                   const unsigned int nx,
                   const unsigned int ny,
                   const Real xmin=0., const Real xmax=1.,
                   const Real ymin=0., const Real ymax=1.,
                   const ElemType type=INVALID_ELEM,
                   const bool gauss_lobatto_grid=false);

/**
 * Meshes a spherical or mapped-spherical domain.
 */
template <typename RealType>
void build_sphere (UnstructuredMeshTempl<RealType> & mesh,
                   const Real rad=1,
                   const unsigned int nr=2,
                   const ElemType type=INVALID_ELEM,
                   const unsigned int n_smooth=2,
                   const bool flat=true);

/**
 * Meshes the tensor product of a 1D and a 1D-or-2D domain.
 */
template <typename RealType>
void build_extrusion (UnstructuredMeshTempl<RealType> & mesh,
                      const MeshBaseTempl<RealType> & cross_section,
                      const unsigned int nz,
                      RealVectorValue extrusion_vector,
                      QueryElemSubdomainIDBaseTempl<RealType> * elem_subdomain = nullptr);

#if defined(LIBMESH_HAVE_TRIANGLE) && LIBMESH_DIM > 1
/**
 * Meshes a rectangular (2D) region (with or without holes) with a
 * Delaunay triangulation.  This function internally calls the
 * triangle library written by J.R. Shewchuk.
 */
template <typename RealType>
void build_delaunay_square(UnstructuredMeshTempl<RealType> & mesh,
                           const unsigned int nx, // num. of elements in x-dir
                           const unsigned int ny, // num. of elements in y-dir
                           const Real xmin, const Real xmax,
                           const Real ymin, const Real ymax,
                           const ElemType type,
                           const std::vector<TriangleInterface::Hole*> * holes=nullptr);
#endif // LIBMESH_HAVE_TRIANGLE && LIBMESH_DIM > 1

/**
 * Class for receiving the callback during extrusion generation and providing user-defined
 * subdomains based on the old (existing) element id and the current layer.
 */
template <typename RealType>
class QueryElemSubdomainIDBaseTempl
{
public:
  virtual ~QueryElemSubdomainIDBaseTempl() {}

  virtual subdomain_id_type get_subdomain_for_layer(const ElemTempl<RealType> * old_elem, unsigned int layer) = 0;
};

typedef QueryElemSubdomainIDBaseTempl<Real> QueryElemSubdomainIDBase;

} // end namespace Meshtools::Generation
} // end namespace MeshTools


} // namespace libMesh

#endif // LIBMESH_MESH_GENERATION_H
