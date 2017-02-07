// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::sqrt


// Local includes
#include "libmesh/mesh_generation.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_edge4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_quad8.h"
#include "libmesh/face_quad9.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex20.h"
#include "libmesh/cell_hex27.h"
#include "libmesh/cell_prism6.h"
#include "libmesh/cell_prism15.h"
#include "libmesh/cell_prism18.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_pyramid5.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/boundary_info.h"
#include "libmesh/sphere.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_smoother_laplace.h"
#include "libmesh/node_elem.h"
#include "libmesh/vector_value.h"
#include "libmesh/function_base.h"

namespace libMesh
{

namespace MeshTools {
namespace Generation {
namespace Private {
/**
 * A useful inline function which replaces the macros
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
      }

    case QUAD8:
    case QUAD9:
    case TRI6:
      {
        return i + j*(2*nx+1);
      }

    default:
      libmesh_error_msg("ERROR: Unrecognized 2D element type.");
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
      }

    case HEX20:
    case HEX27:
    case TET4:  // TET4's are created from an initial HEX27 discretization
    case TET10: // TET10's are created from an initial HEX27 discretization
    case PYRAMID5: // PYRAMID5's are created from an initial HEX27 discretization
    case PYRAMID13:
    case PYRAMID14:
    case PRISM15:
    case PRISM18:
      {
        return i + (2*nx+1)*(j + k*(2*ny+1));
      }

    default:
      libmesh_error_msg("ERROR: Unrecognized element type.");
    }

  return libMesh::invalid_uint;
}


/**
 * This object is passed to MeshTools::Modification::redistribute() to
 * redistribute the points on a uniform grid into the Gauss-Lobatto
 * points on the actual grid.
 */
class GaussLobattoRedistributionFunction : public FunctionBase<Real>
{
public:
  /**
   * Constructor class base class ctor with NULL master.
   */
  GaussLobattoRedistributionFunction(unsigned int nx,
                                     Real xmin,
                                     Real xmax,
                                     unsigned int ny=0,
                                     Real ymin=0,
                                     Real ymax=0,
                                     unsigned int nz=0,
                                     Real zmin=0,
                                     Real zmax=0) :
    FunctionBase<Real>(libmesh_nullptr)
  {
    _nelem.resize(3);
    _nelem[0] = nx;
    _nelem[1] = ny;
    _nelem[2] = nz;

    _mins.resize(3);
    _mins[0] = xmin;
    _mins[1] = ymin;
    _mins[2] = zmin;

    _widths.resize(3);
    _widths[0] = xmax - xmin;
    _widths[1] = ymax - ymin;
    _widths[2] = zmax - zmin;

    // Precompute the cosine values.
    _cosines.resize(3);
    for (unsigned dir=0; dir<3; ++dir)
      if (_nelem[dir] != 0)
        {
          _cosines[dir].resize(_nelem[dir]+1);
          for (std::size_t i=0; i<_cosines[dir].size(); ++i)
            _cosines[dir][i] = std::cos(libMesh::pi * i / _nelem[dir]);
        }
  }

  /**
   * We must provide a way to clone ourselves to satisfy the pure
   * virtual interface.  We use the autogenerated copy constructor.
   */
  virtual UniquePtr<FunctionBase<Real> > clone () const libmesh_override
  {
    return UniquePtr<FunctionBase<Real> > (new GaussLobattoRedistributionFunction(*this));
  }

  /**
   * This is the actual function that
   * MeshTools::Modification::redistribute() calls.  Moves the points
   * of the grid to the Gauss-Lobatto points.
   */
  virtual void operator() (const Point & p,
                           const Real /*time*/,
                           DenseVector<Real> & output) libmesh_override
  {
    output.resize(3);

    for (unsigned dir=0; dir<3; ++dir)
      if (_nelem[dir] != 0)
        {
          // Figure out the index of the current point.
          Real float_index = (p(dir) - _mins[dir]) * _nelem[dir] / _widths[dir];

          // std::modf separates the fractional and integer parts of the index.
          Real integer_part = 0;
          Real fractional_part = std::modf(float_index, &integer_part);

          // Vertex node?
          if (std::abs(fractional_part) < TOLERANCE || std::abs(fractional_part - 1.0) < TOLERANCE)
            {
              int index = round(float_index);

              // Move node to the Gauss-Lobatto position.
              output(dir) = _mins[dir] + _widths[dir] * 0.5 * (1.0 - _cosines[dir][index]);
            }

          // Mid-edge (quadratic) node?
          else if (std::abs(fractional_part - 0.5) < TOLERANCE)
            {
              // Move node to the Gauss-Lobatto position, which is the average of
              // the node to the left and the node to the right.
              output(dir) = _mins[dir] + _widths[dir] * 0.5 *
                (1.0 - 0.5*(_cosines[dir][integer_part] + _cosines[dir][integer_part+1]));
            }

          // 1D only: Left interior (cubic) node?
          else if (std::abs(fractional_part - 1./3.) < TOLERANCE)
            {
              // Move node to the Gauss-Lobatto position, which is
              // 2/3*left_vertex + 1/3*right_vertex.
              output(dir) = _mins[dir] + _widths[dir] * 0.5 *
                (1.0 - 2./3.*_cosines[dir][integer_part] - 1./3.*_cosines[dir][integer_part+1]);
            }

          // 1D only: Right interior (cubic) node?
          else if (std::abs(fractional_part - 2./3.) < TOLERANCE)
            {
              // Move node to the Gauss-Lobatto position, which is
              // 1/3*left_vertex + 2/3*right_vertex.
              output(dir) = _mins[dir] + _widths[dir] * 0.5 *
                (1.0 - 1./3.*_cosines[dir][integer_part] - 2./3.*_cosines[dir][integer_part+1]);
            }

          else
            libmesh_error_msg("Cannot redistribute node: " << p);
        }
  }

  /**
   * We must also override operator() which returns a Real, but this function
   * should never be called, so it's left unimplemented.
   */
  virtual Real operator() (const Point & /*p*/,
                           const Real /*time*/) libmesh_override
  {
    libmesh_not_implemented();
  }

protected:
  // Stored data
  std::vector<Real> _mins;
  std::vector<unsigned int> _nelem;
  std::vector<Real> _widths;

  // Precomputed values
  std::vector<std::vector<Real> > _cosines;
};


} // namespace Private
} // namespace Generation
} // namespace MeshTools

// ------------------------------------------------------------
// MeshTools::Generation function for mesh generation
void MeshTools::Generation::build_cube(UnstructuredMesh & mesh,
                                       const unsigned int nx,
                                       const unsigned int ny,
                                       const unsigned int nz,
                                       const Real xmin, const Real xmax,
                                       const Real ymin, const Real ymax,
                                       const Real zmin, const Real zmax,
                                       const ElemType type,
                                       const bool gauss_lobatto_grid)
{
  START_LOG("build_cube()", "MeshTools::Generation");

  // Declare that we are using the indexing utility routine
  // in the "Private" part of our current namespace.  If this doesn't
  // work in GCC 2.95.3 we can either remove it or stop supporting
  // 2.95.3 altogether.
  // Changing this to import the whole namespace... just importing idx
  // causes an internal compiler error for Intel Compiler 11.0 on Linux
  // in debug mode.
  using namespace MeshTools::Generation::Private;

  // Clear the mesh and start from scratch
  mesh.clear();

  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  if (nz != 0)
    {
      mesh.set_mesh_dimension(3);
      mesh.set_spatial_dimension(3);
    }
  else if (ny != 0)
    {
      mesh.set_mesh_dimension(2);
      mesh.set_spatial_dimension(2);
    }
  else if (nx != 0)
    {
      mesh.set_mesh_dimension(1);
      mesh.set_spatial_dimension(1);
    }
  else
    {
      // Will we get here?
      mesh.set_mesh_dimension(0);
      mesh.set_spatial_dimension(0);
    }

  switch (mesh.mesh_dimension())
    {
      //---------------------------------------------------------------------
      // Build a 0D point
    case 0:
      {
        libmesh_assert_equal_to (nx, 0);
        libmesh_assert_equal_to (ny, 0);
        libmesh_assert_equal_to (nz, 0);

        libmesh_assert (type == INVALID_ELEM || type == NODEELEM);

        // Build one nodal element for the mesh
        mesh.add_point (Point(0, 0, 0), 0);
        Elem * elem = mesh.add_elem (new NodeElem);
        elem->set_node(0) = mesh.node_ptr(0);

        break;
      }



      //---------------------------------------------------------------------
      // Build a 1D line
    case 1:
      {
        libmesh_assert_not_equal_to (nx, 0);
        libmesh_assert_equal_to (ny, 0);
        libmesh_assert_equal_to (nz, 0);
        libmesh_assert_less (xmin, xmax);

        // Reserve elements
        switch (type)
          {
          case INVALID_ELEM:
          case EDGE2:
          case EDGE3:
          case EDGE4:
            {
              mesh.reserve_elem (nx);
              break;
            }

          default:
            libmesh_error_msg("ERROR: Unrecognized 1D element type.");
          }

        // Reserve nodes
        switch (type)
          {
          case INVALID_ELEM:
          case EDGE2:
            {
              mesh.reserve_nodes(nx+1);
              break;
            }

          case EDGE3:
            {
              mesh.reserve_nodes(2*nx+1);
              break;
            }

          case EDGE4:
            {
              mesh.reserve_nodes(3*nx+1);
              break;
            }

          default:
            libmesh_error_msg("ERROR: Unrecognized 1D element type.");
          }


        // Build the nodes, depends on whether we're using linears,
        // quadratics or cubics and whether using uniform grid or Gauss-Lobatto
        unsigned int node_id = 0;
        switch(type)
          {
          case INVALID_ELEM:
          case EDGE2:
            {
              for (unsigned int i=0; i<=nx; i++)
                mesh.add_point (Point(static_cast<Real>(i)/nx, 0, 0), node_id++);

              break;
            }

          case EDGE3:
            {
              for (unsigned int i=0; i<=2*nx; i++)
                mesh.add_point (Point(static_cast<Real>(i)/(2*nx), 0, 0), node_id++);
              break;
            }

          case EDGE4:
            {
              for (unsigned int i=0; i<=3*nx; i++)
                mesh.add_point (Point(static_cast<Real>(i)/(3*nx), 0, 0), node_id++);

              break;
            }

          default:
            libmesh_error_msg("ERROR: Unrecognized 1D element type.");

          }

        // Build the elements of the mesh
        switch(type)
          {
          case INVALID_ELEM:
          case EDGE2:
            {
              for (unsigned int i=0; i<nx; i++)
                {
                  Elem * elem = mesh.add_elem (new Edge2);
                  elem->set_node(0) = mesh.node_ptr(i);
                  elem->set_node(1) = mesh.node_ptr(i+1);

                  if (i == 0)
                    boundary_info.add_side(elem, 0, 0);

                  if (i == (nx-1))
                    boundary_info.add_side(elem, 1, 1);
                }
              break;
            }

          case EDGE3:
            {
              for (unsigned int i=0; i<nx; i++)
                {
                  Elem * elem = mesh.add_elem (new Edge3);
                  elem->set_node(0) = mesh.node_ptr(2*i);
                  elem->set_node(2) = mesh.node_ptr(2*i+1);
                  elem->set_node(1) = mesh.node_ptr(2*i+2);

                  if (i == 0)
                    boundary_info.add_side(elem, 0, 0);

                  if (i == (nx-1))
                    boundary_info.add_side(elem, 1, 1);
                }
              break;
            }

          case EDGE4:
            {
              for (unsigned int i=0; i<nx; i++)
                {
                  Elem * elem = mesh.add_elem (new Edge4);
                  elem->set_node(0) = mesh.node_ptr(3*i);
                  elem->set_node(2) = mesh.node_ptr(3*i+1);
                  elem->set_node(3) = mesh.node_ptr(3*i+2);
                  elem->set_node(1) = mesh.node_ptr(3*i+3);

                  if (i == 0)
                    boundary_info.add_side(elem, 0, 0);

                  if (i == (nx-1))
                    boundary_info.add_side(elem, 1, 1);
                }
              break;
            }

          default:
            libmesh_error_msg("ERROR: Unrecognized 1D element type.");
          }

        // Move the nodes to their final locations.
        if (gauss_lobatto_grid)
          {
            GaussLobattoRedistributionFunction func(nx, xmin, xmax);
            MeshTools::Modification::redistribute(mesh, func);
          }
        else // !gauss_lobatto_grid
          {
            for (unsigned int p=0; p<mesh.n_nodes(); p++)
              mesh.node_ref(p)(0) = (mesh.node_ref(p)(0))*(xmax-xmin) + xmin;
          }

        // Add sideset names to boundary info
        boundary_info.sideset_name(0) = "left";
        boundary_info.sideset_name(1) = "right";

        // Add nodeset names to boundary info
        boundary_info.nodeset_name(0) = "left";
        boundary_info.nodeset_name(1) = "right";

        break;
      }










      //---------------------------------------------------------------------
      // Build a 2D quadrilateral
    case 2:
      {
        libmesh_assert_not_equal_to (nx, 0);
        libmesh_assert_not_equal_to (ny, 0);
        libmesh_assert_equal_to (nz, 0);
        libmesh_assert_less (xmin, xmax);
        libmesh_assert_less (ymin, ymax);

        // Reserve elements.  The TRI3 and TRI6 meshes
        // have twice as many elements...
        switch (type)
          {
          case INVALID_ELEM:
          case QUAD4:
          case QUAD8:
          case QUAD9:
            {
              mesh.reserve_elem (nx*ny);
              break;
            }

          case TRI3:
          case TRI6:
            {
              mesh.reserve_elem (2*nx*ny);
              break;
            }

          default:
            libmesh_error_msg("ERROR: Unrecognized 2D element type.");
          }



        // Reserve nodes.  The quadratic element types
        // need to reserve more nodes than the linear types.
        switch (type)
          {
          case INVALID_ELEM:
          case QUAD4:
          case TRI3:
            {
              mesh.reserve_nodes( (nx+1)*(ny+1) );
              break;
            }

          case QUAD8:
          case QUAD9:
          case TRI6:
            {
              mesh.reserve_nodes( (2*nx+1)*(2*ny+1) );
              break;
            }


          default:
            libmesh_error_msg("ERROR: Unrecognized 2D element type.");
          }



        // Build the nodes. Depends on whether you are using a linear
        // or quadratic element, and whether you are using a uniform
        // grid or the Gauss-Lobatto grid points.
        unsigned int node_id = 0;
        switch (type)
          {
          case INVALID_ELEM:
          case QUAD4:
          case TRI3:
            {
              for (unsigned int j=0; j<=ny; j++)
                for (unsigned int i=0; i<=nx; i++)
                  mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(nx),
                                        static_cast<Real>(j)/static_cast<Real>(ny),
                                        0.), node_id++);

              break;
            }

          case QUAD8:
          case QUAD9:
          case TRI6:
            {
              for (unsigned int j=0; j<=(2*ny); j++)
                for (unsigned int i=0; i<=(2*nx); i++)
                  mesh.add_point (Point(static_cast<Real>(i)/static_cast<Real>(2*nx),
                                        static_cast<Real>(j)/static_cast<Real>(2*ny),
                                        0), node_id++);

              break;
            }


          default:
            libmesh_error_msg("ERROR: Unrecognized 2D element type.");
          }






        // Build the elements.  Each one is a bit different.
        switch (type)
          {

          case INVALID_ELEM:
          case QUAD4:
            {
              for (unsigned int j=0; j<ny; j++)
                for (unsigned int i=0; i<nx; i++)
                  {
                    Elem * elem = mesh.add_elem(new Quad4);

                    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
                    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j)  );
                    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+1,j+1));
                    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i,j+1)  );

                    if (j == 0)
                      boundary_info.add_side(elem, 0, 0);

                    if (j == (ny-1))
                      boundary_info.add_side(elem, 2, 2);

                    if (i == 0)
                      boundary_info.add_side(elem, 3, 3);

                    if (i == (nx-1))
                      boundary_info.add_side(elem, 1, 1);
                  }
              break;
            }


          case TRI3:
            {
              for (unsigned int j=0; j<ny; j++)
                for (unsigned int i=0; i<nx; i++)
                  {
                    Elem * elem = libmesh_nullptr;

                    // Add first Tri3
                    elem = mesh.add_elem(new Tri3);

                    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
                    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j)  );
                    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+1,j+1));

                    if (j == 0)
                      boundary_info.add_side(elem, 0, 0);

                    if (i == (nx-1))
                      boundary_info.add_side(elem, 1, 1);

                    // Add second Tri3
                    elem = mesh.add_elem(new Tri3);

                    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
                    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+1,j+1));
                    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i,j+1)  );

                    if (j == (ny-1))
                      boundary_info.add_side(elem, 1, 2);

                    if (i == 0)
                      boundary_info.add_side(elem, 2, 3);
                  }
              break;
            }



          case QUAD8:
          case QUAD9:
            {
              for (unsigned int j=0; j<(2*ny); j += 2)
                for (unsigned int i=0; i<(2*nx); i += 2)
                  {
                    Elem * elem = (type == QUAD8) ?
                      mesh.add_elem(new Quad8) :
                      mesh.add_elem(new Quad9);


                    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
                    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+2,j)  );
                    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+2,j+2));
                    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i,j+2)  );
                    elem->set_node(4) = mesh.node_ptr(idx(type,nx,i+1,j)  );
                    elem->set_node(5) = mesh.node_ptr(idx(type,nx,i+2,j+1));
                    elem->set_node(6) = mesh.node_ptr(idx(type,nx,i+1,j+2));
                    elem->set_node(7) = mesh.node_ptr(idx(type,nx,i,j+1)  );
                    if (type == QUAD9)
                      elem->set_node(8) = mesh.node_ptr(idx(type,nx,i+1,j+1));


                    if (j == 0)
                      boundary_info.add_side(elem, 0, 0);

                    if (j == 2*(ny-1))
                      boundary_info.add_side(elem, 2, 2);

                    if (i == 0)
                      boundary_info.add_side(elem, 3, 3);

                    if (i == 2*(nx-1))
                      boundary_info.add_side(elem, 1, 1);
                  }
              break;
            }


          case TRI6:
            {
              for (unsigned int j=0; j<(2*ny); j += 2)
                for (unsigned int i=0; i<(2*nx); i += 2)
                  {
                    Elem * elem = libmesh_nullptr;

                    // Add first Tri6
                    elem = mesh.add_elem(new Tri6);

                    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
                    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+2,j)  );
                    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i+2,j+2));
                    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i+1,j)  );
                    elem->set_node(4) = mesh.node_ptr(idx(type,nx,i+2,j+1));
                    elem->set_node(5) = mesh.node_ptr(idx(type,nx,i+1,j+1));

                    if (j == 0)
                      boundary_info.add_side(elem, 0, 0);

                    if (i == 2*(nx-1))
                      boundary_info.add_side(elem, 1, 1);

                    // Add second Tri6
                    elem = mesh.add_elem(new Tri6);

                    elem->set_node(0) = mesh.node_ptr(idx(type,nx,i,j)    );
                    elem->set_node(1) = mesh.node_ptr(idx(type,nx,i+2,j+2));
                    elem->set_node(2) = mesh.node_ptr(idx(type,nx,i,j+2)  );
                    elem->set_node(3) = mesh.node_ptr(idx(type,nx,i+1,j+1));
                    elem->set_node(4) = mesh.node_ptr(idx(type,nx,i+1,j+2));
                    elem->set_node(5) = mesh.node_ptr(idx(type,nx,i,j+1)  );

                    if (j == 2*(ny-1))
                      boundary_info.add_side(elem, 1, 2);

                    if (i == 0)
                      boundary_info.add_side(elem, 2, 3);

                  }
              break;
            };


          default:
            libmesh_error_msg("ERROR: Unrecognized 2D element type.");
          }




        // Scale the nodal positions
        if (gauss_lobatto_grid)
          {
            GaussLobattoRedistributionFunction func(nx, xmin, xmax,
                                                    ny, ymin, ymax);
            MeshTools::Modification::redistribute(mesh, func);
          }
        else // !gauss_lobatto_grid
          {
            for (unsigned int p=0; p<mesh.n_nodes(); p++)
              {
                mesh.node_ref(p)(0) = (mesh.node_ref(p)(0))*(xmax-xmin) + xmin;
                mesh.node_ref(p)(1) = (mesh.node_ref(p)(1))*(ymax-ymin) + ymin;
              }
          }

        // Add sideset names to boundary info
        boundary_info.sideset_name(0) = "bottom";
        boundary_info.sideset_name(1) = "right";
        boundary_info.sideset_name(2) = "top";
        boundary_info.sideset_name(3) = "left";

        // Add nodeset names to boundary info
        boundary_info.nodeset_name(0) = "bottom";
        boundary_info.nodeset_name(1) = "right";
        boundary_info.nodeset_name(2) = "top";
        boundary_info.nodeset_name(3) = "left";

        break;
      }











      //---------------------------------------------------------------------
      // Build a 3D mesh using hexes, tets, prisms, or pyramids.
    case 3:
      {
        libmesh_assert_not_equal_to (nx, 0);
        libmesh_assert_not_equal_to (ny, 0);
        libmesh_assert_not_equal_to (nz, 0);
        libmesh_assert_less (xmin, xmax);
        libmesh_assert_less (ymin, ymax);
        libmesh_assert_less (zmin, zmax);


        // Reserve elements.  Meshes with prismatic elements require
        // twice as many elements.
        switch (type)
          {
          case INVALID_ELEM:
          case HEX8:
          case HEX20:
          case HEX27:
          case TET4:  // TET4's are created from an initial HEX27 discretization
          case TET10: // TET10's are created from an initial HEX27 discretization
          case PYRAMID5: // PYRAMIDs are created from an initial HEX27 discretization
          case PYRAMID13:
          case PYRAMID14:
            {
              mesh.reserve_elem(nx*ny*nz);
              break;
            }

          case PRISM6:
          case PRISM15:
          case PRISM18:
            {
              mesh.reserve_elem(2*nx*ny*nz);
              break;
            }

          default:
            libmesh_error_msg("ERROR: Unrecognized 3D element type.");
          }





        // Reserve nodes.  Quadratic elements need twice as many nodes as linear elements.
        switch (type)
          {
          case INVALID_ELEM:
          case HEX8:
          case PRISM6:
            {
              mesh.reserve_nodes( (nx+1)*(ny+1)*(nz+1) );
              break;
            }

          case HEX20:
          case HEX27:
          case TET4: // TET4's are created from an initial HEX27 discretization
          case TET10: // TET10's are created from an initial HEX27 discretization
          case PYRAMID5: // PYRAMIDs are created from an initial HEX27 discretization
          case PYRAMID13:
          case PYRAMID14:
          case PRISM15:
          case PRISM18:
            {
              // FYI: The resulting TET4 mesh will have exactly
              // 5*(nx*ny*nz) + 2*(nx*ny + nx*nz + ny*nz) + (nx+ny+nz) + 1
              // nodes once the additional mid-edge nodes for the HEX27 discretization
              // have been deleted.
              mesh.reserve_nodes( (2*nx+1)*(2*ny+1)*(2*nz+1) );
              break;
            }

          default:
            libmesh_error_msg("ERROR: Unrecognized 3D element type.");
          }




        // Build the nodes.
        unsigned int node_id = 0;
        switch (type)
          {
          case INVALID_ELEM:
          case HEX8:
          case PRISM6:
            {
              for (unsigned int k=0; k<=nz; k++)
                for (unsigned int j=0; j<=ny; j++)
                  for (unsigned int i=0; i<=nx; i++)
                    mesh.add_point(Point(static_cast<Real>(i)/static_cast<Real>(nx),
                                         static_cast<Real>(j)/static_cast<Real>(ny),
                                         static_cast<Real>(k)/static_cast<Real>(nz)), node_id++);

              break;
            }

          case HEX20:
          case HEX27:
          case TET4: // TET4's are created from an initial HEX27 discretization
          case TET10: // TET10's are created from an initial HEX27 discretization
          case PYRAMID5: // PYRAMIDs are created from an initial HEX27 discretization
          case PYRAMID13:
          case PYRAMID14:
          case PRISM15:
          case PRISM18:
            {
              for (unsigned int k=0; k<=(2*nz); k++)
                for (unsigned int j=0; j<=(2*ny); j++)
                  for (unsigned int i=0; i<=(2*nx); i++)
                    mesh.add_point(Point(static_cast<Real>(i)/static_cast<Real>(2*nx),
                                         static_cast<Real>(j)/static_cast<Real>(2*ny),
                                         static_cast<Real>(k)/static_cast<Real>(2*nz)), node_id++);

              break;
            }


          default:
            libmesh_error_msg("ERROR: Unrecognized 3D element type.");
          }




        // Build the elements.
        switch (type)
          {
          case INVALID_ELEM:
          case HEX8:
            {
              for (unsigned int k=0; k<nz; k++)
                for (unsigned int j=0; j<ny; j++)
                  for (unsigned int i=0; i<nx; i++)
                    {
                      Elem * elem = mesh.add_elem(new Hex8);

                      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i,j,k)      );
                      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
                      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
                      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
                      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i,j,k+1)    );
                      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
                      elem->set_node(6) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
                      elem->set_node(7) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );

                      if (k == 0)
                        boundary_info.add_side(elem, 0, 0);

                      if (k == (nz-1))
                        boundary_info.add_side(elem, 5, 5);

                      if (j == 0)
                        boundary_info.add_side(elem, 1, 1);

                      if (j == (ny-1))
                        boundary_info.add_side(elem, 3, 3);

                      if (i == 0)
                        boundary_info.add_side(elem, 4, 4);

                      if (i == (nx-1))
                        boundary_info.add_side(elem, 2, 2);
                    }
              break;
            }




          case PRISM6:
            {
              for (unsigned int k=0; k<nz; k++)
                for (unsigned int j=0; j<ny; j++)
                  for (unsigned int i=0; i<nx; i++)
                    {
                      // First Prism
                      Elem * elem = libmesh_nullptr;
                      elem = mesh.add_elem(new Prism6);

                      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i,j,k)      );
                      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
                      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
                      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i,j,k+1)    );
                      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
                      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );

                      // Add sides for first prism to boundary info object
                      if (i==0)
                        boundary_info.add_side(elem, 3, 4);

                      if (j==0)
                        boundary_info.add_side(elem, 1, 1);

                      if (k==0)
                        boundary_info.add_side(elem, 0, 0);

                      if (k == (nz-1))
                        boundary_info.add_side(elem, 4, 5);

                      // Second Prism
                      elem = mesh.add_elem(new Prism6);

                      elem->set_node(0) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k)    );
                      elem->set_node(1) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
                      elem->set_node(2) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k)    );
                      elem->set_node(3) = mesh.node_ptr(idx(type,nx,ny,i+1,j,k+1)  );
                      elem->set_node(4) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
                      elem->set_node(5) = mesh.node_ptr(idx(type,nx,ny,i,j+1,k+1)  );

                      // Add sides for second prism to boundary info object
                      if (i == (nx-1))
                        boundary_info.add_side(elem, 1, 2);

                      if (j == (ny-1))
                        boundary_info.add_side(elem, 2, 3);

                      if (k==0)
                        boundary_info.add_side(elem, 0, 0);

                      if (k == (nz-1))
                        boundary_info.add_side(elem, 4, 5);
                    }
              break;
            }






          case HEX20:
          case HEX27:
          case TET4: // TET4's are created from an initial HEX27 discretization
          case TET10: // TET10's are created from an initial HEX27 discretization
          case PYRAMID5: // PYRAMIDs are created from an initial HEX27 discretization
          case PYRAMID13:
          case PYRAMID14:
            {
              for (unsigned int k=0; k<(2*nz); k += 2)
                for (unsigned int j=0; j<(2*ny); j += 2)
                  for (unsigned int i=0; i<(2*nx); i += 2)
                    {
                      Elem * elem = (type == HEX20) ?
                        mesh.add_elem(new Hex20) :
                        mesh.add_elem(new Hex27);

                      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k)  );
                      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k)  );
                      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k)  );
                      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k)  );
                      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+2));
                      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+2));
                      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+2));
                      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+2));
                      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k)  );
                      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k)  );
                      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k)  );
                      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k)  );
                      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+1));
                      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+1));
                      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+1));
                      elem->set_node(15) = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+1));
                      elem->set_node(16) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+2));
                      elem->set_node(17) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+2));
                      elem->set_node(18) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+2));
                      elem->set_node(19) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+2));
                      if ((type == HEX27) || (type == TET4) || (type == TET10) ||
                          (type == PYRAMID5) || (type == PYRAMID13) || (type == PYRAMID14))
                        {
                          elem->set_node(20) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
                          elem->set_node(21) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+1));
                          elem->set_node(22) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+1));
                          elem->set_node(23) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+1));
                          elem->set_node(24) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+1));
                          elem->set_node(25) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
                          elem->set_node(26) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
                        }


                      if (k == 0)
                        boundary_info.add_side(elem, 0, 0);

                      if (k == 2*(nz-1))
                        boundary_info.add_side(elem, 5, 5);

                      if (j == 0)
                        boundary_info.add_side(elem, 1, 1);

                      if (j == 2*(ny-1))
                        boundary_info.add_side(elem, 3, 3);

                      if (i == 0)
                        boundary_info.add_side(elem, 4, 4);

                      if (i == 2*(nx-1))
                        boundary_info.add_side(elem, 2, 2);
                    }
              break;
            }




          case PRISM15:
          case PRISM18:
            {
              for (unsigned int k=0; k<(2*nz); k += 2)
                for (unsigned int j=0; j<(2*ny); j += 2)
                  for (unsigned int i=0; i<(2*nx); i += 2)
                    {
                      // First Prism
                      Elem * elem = libmesh_nullptr;
                      elem = ((type == PRISM15) ?
                              mesh.add_elem(new Prism15) :
                              mesh.add_elem(new Prism18));

                      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k)  );
                      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k)  );
                      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k)  );
                      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+2));
                      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+2));
                      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+2));
                      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k)  );
                      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
                      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k)  );
                      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i,  j,  k+1));
                      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+2,j,  k+1));
                      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,  j+2,k+1));
                      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+2));
                      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
                      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+2));
                      if (type == PRISM18)
                        {
                          elem->set_node(15) = mesh.node_ptr(idx(type,nx,ny,i+1,j,  k+1));
                          elem->set_node(16) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
                          elem->set_node(17) = mesh.node_ptr(idx(type,nx,ny,i,  j+1,k+1));
                        }

                      // Add sides for first prism to boundary info object
                      if (i==0)
                        boundary_info.add_side(elem, 3, 4);

                      if (j==0)
                        boundary_info.add_side(elem, 1, 1);

                      if (k==0)
                        boundary_info.add_side(elem, 0, 0);

                      if (k == 2*(nz-1))
                        boundary_info.add_side(elem, 4, 5);


                      // Second Prism
                      elem = ((type == PRISM15) ?
                              mesh.add_elem(new Prism15) :
                              mesh.add_elem(new Prism18));

                      elem->set_node(0)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k)     );
                      elem->set_node(1)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k)   );
                      elem->set_node(2)  = mesh.node_ptr(idx(type,nx,ny,i,j+2,k)     );
                      elem->set_node(3)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k+2)   );
                      elem->set_node(4)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+2) );
                      elem->set_node(5)  = mesh.node_ptr(idx(type,nx,ny,i,j+2,k+2)   );
                      elem->set_node(6)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k)  );
                      elem->set_node(7)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k)  );
                      elem->set_node(8)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k)  );
                      elem->set_node(9)  = mesh.node_ptr(idx(type,nx,ny,i+2,j,k+1)  );
                      elem->set_node(10) = mesh.node_ptr(idx(type,nx,ny,i+2,j+2,k+1));
                      elem->set_node(11) = mesh.node_ptr(idx(type,nx,ny,i,j+2,k+1)  );
                      elem->set_node(12) = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+2));
                      elem->set_node(13) = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+2));
                      elem->set_node(14) = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+2));
                      if (type == PRISM18)
                        {
                          elem->set_node(15)  = mesh.node_ptr(idx(type,nx,ny,i+2,j+1,k+1));
                          elem->set_node(16)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+2,k+1));
                          elem->set_node(17)  = mesh.node_ptr(idx(type,nx,ny,i+1,j+1,k+1));
                        }

                      // Add sides for second prism to boundary info object
                      if (i == 2*(nx-1))
                        boundary_info.add_side(elem, 1, 2);

                      if (j == 2*(ny-1))
                        boundary_info.add_side(elem, 2, 3);

                      if (k==0)
                        boundary_info.add_side(elem, 0, 0);

                      if (k == 2*(nz-1))
                        boundary_info.add_side(elem, 4, 5);

                    }
              break;
            }





          default:
            libmesh_error_msg("ERROR: Unrecognized 3D element type.");
          }




        //.......................................
        // Scale the nodal positions
        if (gauss_lobatto_grid)
          {
            GaussLobattoRedistributionFunction func(nx, xmin, xmax,
                                                    ny, ymin, ymax,
                                                    nz, zmin, zmax);
            MeshTools::Modification::redistribute(mesh, func);
          }
        else // !gauss_lobatto_grid
          {
            for (unsigned int p=0; p<mesh.n_nodes(); p++)
              {
                mesh.node_ref(p)(0) = (mesh.node_ref(p)(0))*(xmax-xmin) + xmin;
                mesh.node_ref(p)(1) = (mesh.node_ref(p)(1))*(ymax-ymin) + ymin;
                mesh.node_ref(p)(2) = (mesh.node_ref(p)(2))*(zmax-zmin) + zmin;
              }
          }



        // Additional work for tets and pyramids: we take the existing
        // HEX27 discretization and split each element into 24
        // sub-tets or 6 sub-pyramids.
        //
        // 24 isn't the minimum-possible number of tets, but it
        // obviates any concerns about the edge orientations between
        // the various elements.
        if ((type == TET4) ||
            (type == TET10) ||
            (type == PYRAMID5) ||
            (type == PYRAMID13) ||
            (type == PYRAMID14))
          {
            // Temporary storage for new elements. (24 tets per hex, 6 pyramids)
            std::vector<Elem *> new_elements;

            if ((type == TET4) || (type == TET10))
              new_elements.reserve(24*mesh.n_elem());
            else
              new_elements.reserve(6*mesh.n_elem());

            // Create tetrahedra or pyramids
            {
              MeshBase::element_iterator       el     = mesh.elements_begin();
              const MeshBase::element_iterator end_el = mesh.elements_end();

              for ( ; el != end_el;  ++el)
                {
                  // Get a pointer to the HEX27 element.
                  Elem * base_hex = *el;

                  // Get a pointer to the node located at the HEX27 centroid
                  Node * apex_node = base_hex->node_ptr(26);

                  // Container to catch ids handed back from BoundaryInfo
                  std::vector<boundary_id_type> ids;

                  for (unsigned short s=0; s<base_hex->n_sides(); ++s)
                    {
                      // Get the boundary ID(s) for this side
                      boundary_info.boundary_ids(*el, s, ids);

                      // We're creating this Mesh, so there should be 0 or 1 boundary IDs.
                      libmesh_assert(ids.size() <= 1);

                      // A convenient name for the side's ID.
                      boundary_id_type b_id = ids.empty() ? BoundaryInfo::invalid_id : ids[0];

                      // Need to build the full-ordered side!
                      UniquePtr<Elem> side = base_hex->build_side_ptr(s);

                      if ((type == TET4) || (type == TET10))
                        {
                          // Build 4 sub-tets per side
                          for (unsigned int sub_tet=0; sub_tet<4; ++sub_tet)
                            {
                              new_elements.push_back( new Tet4 );
                              Elem * sub_elem = new_elements.back();
                              sub_elem->set_node(0) = side->node_ptr(sub_tet);
                              sub_elem->set_node(1) = side->node_ptr(8);                           // centroid of the face
                              sub_elem->set_node(2) = side->node_ptr(sub_tet==3 ? 0 : sub_tet+1 ); // wrap-around
                              sub_elem->set_node(3) = apex_node;                                   // apex node always used!

                              // If the original hex was a boundary hex, add the new sub_tet's side
                              // 0 with the same b_id.  Note: the tets are all aligned so that their
                              // side 0 is on the boundary.
                              if (b_id != BoundaryInfo::invalid_id)
                                boundary_info.add_side(sub_elem, 0, b_id);
                            }
                        } // end if ((type == TET4) || (type == TET10))

                      else // type==PYRAMID5 || type==PYRAMID13 || type==PYRAMID14
                        {
                          // Build 1 sub-pyramid per side.
                          new_elements.push_back(new Pyramid5);
                          Elem * sub_elem = new_elements.back();

                          // Set the base.  Note that since the apex is *inside* the base_hex,
                          // and the pyramid uses a counter-clockwise base numbering, we need to
                          // reverse the [1] and [3] node indices.
                          sub_elem->set_node(0) = side->node_ptr(0);
                          sub_elem->set_node(1) = side->node_ptr(3);
                          sub_elem->set_node(2) = side->node_ptr(2);
                          sub_elem->set_node(3) = side->node_ptr(1);

                          // Set the apex
                          sub_elem->set_node(4) = apex_node;

                          // If the original hex was a boundary hex, add the new sub_pyr's side
                          // 4 (the square base) with the same b_id.
                          if (b_id != BoundaryInfo::invalid_id)
                            boundary_info.add_side(sub_elem, 4, b_id);
                        } // end else type==PYRAMID5 || type==PYRAMID13 || type==PYRAMID14
                    }
                }
            }


            // Delete the original HEX27 elements from the mesh, and the boundary info structure.
            {
              MeshBase::element_iterator       el     = mesh.elements_begin();
              const MeshBase::element_iterator end_el = mesh.elements_end();

              for ( ; el != end_el;  ++el)
                {
                  boundary_info.remove(*el); // Safe even if *el has no boundary info.
                  mesh.delete_elem(*el);
                }
            }

            // Add the new elements
            for (std::size_t i=0; i<new_elements.size(); ++i)
              mesh.add_elem(new_elements[i]);

          } // end if (type == TET4,TET10,PYRAMID5,PYRAMID13,PYRAMID14


        // Use all_second_order to convert the TET4's to TET10's or PYRAMID5's to PYRAMID14's
        if ((type == TET10) || (type == PYRAMID14))
          mesh.all_second_order();

        else if (type == PYRAMID13)
          mesh.all_second_order(/*full_ordered=*/false);


        // Add sideset names to boundary info (Z axis out of the screen)
        boundary_info.sideset_name(0) = "back";
        boundary_info.sideset_name(1) = "bottom";
        boundary_info.sideset_name(2) = "right";
        boundary_info.sideset_name(3) = "top";
        boundary_info.sideset_name(4) = "left";
        boundary_info.sideset_name(5) = "front";

        // Add nodeset names to boundary info
        boundary_info.nodeset_name(0) = "back";
        boundary_info.nodeset_name(1) = "bottom";
        boundary_info.nodeset_name(2) = "right";
        boundary_info.nodeset_name(3) = "top";
        boundary_info.nodeset_name(4) = "left";
        boundary_info.nodeset_name(5) = "front";

        break;
      } // end case dim==3

    default:
      libmesh_error_msg("Unknown dimension " << mesh.mesh_dimension());
    }

  STOP_LOG("build_cube()", "MeshTools::Generation");



  // Done building the mesh.  Now prepare it for use.
  mesh.prepare_for_use (/*skip_renumber =*/ false);
}



void MeshTools::Generation::build_point (UnstructuredMesh & mesh,
                                         const ElemType type,
                                         const bool gauss_lobatto_grid)
{
  // This method only makes sense in 0D!
  // But we now just turn a non-0D mesh into a 0D mesh
  //libmesh_assert_equal_to (mesh.mesh_dimension(), 1);

  build_cube(mesh,
             0, 0, 0,
             0., 0.,
             0., 0.,
             0., 0.,
             type,
             gauss_lobatto_grid);
}


void MeshTools::Generation::build_line (UnstructuredMesh & mesh,
                                        const unsigned int nx,
                                        const Real xmin, const Real xmax,
                                        const ElemType type,
                                        const bool gauss_lobatto_grid)
{
  // This method only makes sense in 1D!
  // But we now just turn a non-1D mesh into a 1D mesh
  //libmesh_assert_equal_to (mesh.mesh_dimension(), 1);

  build_cube(mesh,
             nx, 0, 0,
             xmin, xmax,
             0., 0.,
             0., 0.,
             type,
             gauss_lobatto_grid);
}



void MeshTools::Generation::build_square (UnstructuredMesh & mesh,
                                          const unsigned int nx,
                                          const unsigned int ny,
                                          const Real xmin, const Real xmax,
                                          const Real ymin, const Real ymax,
                                          const ElemType type,
                                          const bool gauss_lobatto_grid)
{
  // This method only makes sense in 2D!
  // But we now just turn a non-2D mesh into a 2D mesh
  //libmesh_assert_equal_to (mesh.mesh_dimension(), 2);

  // Call the build_cube() member to actually do the work for us.
  build_cube (mesh,
              nx, ny, 0,
              xmin, xmax,
              ymin, ymax,
              0., 0.,
              type,
              gauss_lobatto_grid);
}









#ifndef LIBMESH_ENABLE_AMR
void MeshTools::Generation::build_sphere (UnstructuredMesh &,
                                          const Real,
                                          const unsigned int,
                                          const ElemType,
                                          const unsigned int,
                                          const bool)
{
  libmesh_error_msg("Building a circle/sphere only works with AMR.");
}

#else

void MeshTools::Generation::build_sphere (UnstructuredMesh & mesh,
                                          const Real rad,
                                          const unsigned int nr,
                                          const ElemType type,
                                          const unsigned int n_smooth,
                                          const bool flat)
{
  libmesh_assert_greater (rad, 0.);
  //libmesh_assert_greater (nr, 0); // must refine at least once otherwise will end up with a square/cube

  START_LOG("build_sphere()", "MeshTools::Generation");

  // Clear the mesh and start from scratch, but save the original
  // mesh_dimension, since the original intent of this function was to
  // allow the geometric entity (line, circle, ball, sphere)
  // constructed to be determined by the mesh's dimension.
  unsigned int orig_mesh_dimension = mesh.mesh_dimension();
  mesh.clear();
  mesh.set_mesh_dimension(orig_mesh_dimension);

  // If mesh.mesh_dimension()==1, it *could* be because the user
  // constructed a Mesh without specifying a dimension (since this is
  // allowed now) and hence it got the default dimension of 1.  In
  // this case, we will try to infer the dimension they *really*
  // wanted from the requested ElemType, and if they don't match, go
  // with the ElemType.
  if (mesh.mesh_dimension() == 1)
    {
      if (type==HEX8 || type==HEX27)
        mesh.set_mesh_dimension(3);

      else if (type==TRI3 || type==QUAD4)
        mesh.set_mesh_dimension(2);

      else if (type==EDGE2 || type==EDGE3 || type==EDGE4 || type==INVALID_ELEM)
        mesh.set_mesh_dimension(1);

      else
        libmesh_error_msg("build_sphere(): Please specify a mesh dimension or a valid ElemType (EDGE{2,3,4}, TRI3, QUAD4, HEX{8,27})");
    }

  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  // Sphere is centered at origin by default
  const Point cent;

  const Sphere sphere (cent, rad);

  switch (mesh.mesh_dimension())
    {
      //-----------------------------------------------------------------
      // Build a line in one dimension
    case 1:
      {
        build_line (mesh, 3, -rad, rad, type);
      }




      //-----------------------------------------------------------------
      // Build a circle or hollow sphere in two dimensions
    case 2:
      {
        // For DistributedMesh, if we don't specify node IDs the Mesh
        // will try to pick an appropriate (unique) one for us.  But
        // since we are adding these nodes on all processors, we want
        // to be sure they have consistent IDs across all processors.
        unsigned node_id = 0;

        if (flat)
          {
            const Real sqrt_2     = std::sqrt(2.);
            const Real rad_2      = .25*rad;
            const Real rad_sqrt_2 = rad/sqrt_2;

            // (Temporary) convenient storage for node pointers
            std::vector<Node *> nodes(8);

            // Point 0
            nodes[0] = mesh.add_point (Point(-rad_2,-rad_2, 0.), node_id++);

            // Point 1
            nodes[1] = mesh.add_point (Point( rad_2,-rad_2, 0.), node_id++);

            // Point 2
            nodes[2] = mesh.add_point (Point( rad_2, rad_2, 0.), node_id++);

            // Point 3
            nodes[3] = mesh.add_point (Point(-rad_2, rad_2, 0.), node_id++);

            // Point 4
            nodes[4] = mesh.add_point (Point(-rad_sqrt_2,-rad_sqrt_2, 0.), node_id++);

            // Point 5
            nodes[5] = mesh.add_point (Point( rad_sqrt_2,-rad_sqrt_2, 0.), node_id++);

            // Point 6
            nodes[6] = mesh.add_point (Point( rad_sqrt_2, rad_sqrt_2, 0.), node_id++);

            // Point 7
            nodes[7] = mesh.add_point (Point(-rad_sqrt_2, rad_sqrt_2, 0.), node_id++);

            // Build the elements & set node pointers

            // Element 0
            {
              Elem * elem0 = mesh.add_elem (new Quad4);
              elem0->set_node(0) = nodes[0];
              elem0->set_node(1) = nodes[1];
              elem0->set_node(2) = nodes[2];
              elem0->set_node(3) = nodes[3];
            }

            // Element 1
            {
              Elem * elem1 = mesh.add_elem (new Quad4);
              elem1->set_node(0) = nodes[4];
              elem1->set_node(1) = nodes[0];
              elem1->set_node(2) = nodes[3];
              elem1->set_node(3) = nodes[7];
            }

            // Element 2
            {
              Elem * elem2 = mesh.add_elem (new Quad4);
              elem2->set_node(0) = nodes[4];
              elem2->set_node(1) = nodes[5];
              elem2->set_node(2) = nodes[1];
              elem2->set_node(3) = nodes[0];
            }

            // Element 3
            {
              Elem * elem3 = mesh.add_elem (new Quad4);
              elem3->set_node(0) = nodes[1];
              elem3->set_node(1) = nodes[5];
              elem3->set_node(2) = nodes[6];
              elem3->set_node(3) = nodes[2];
            }

            // Element 4
            {
              Elem * elem4 = mesh.add_elem (new Quad4);
              elem4->set_node(0) = nodes[3];
              elem4->set_node(1) = nodes[2];
              elem4->set_node(2) = nodes[6];
              elem4->set_node(3) = nodes[7];
            }

          }
        else
          {
            // Create the 12 vertices of a regular unit icosahedron
            Real t = 0.5 * (1 + std::sqrt(5.0));
            Real s = rad / std::sqrt(1 + t*t);
            t *= s;

            mesh.add_point (Point(-s,  t,  0), node_id++);
            mesh.add_point (Point( s,  t,  0), node_id++);
            mesh.add_point (Point(-s, -t,  0), node_id++);
            mesh.add_point (Point( s, -t,  0), node_id++);

            mesh.add_point (Point( 0, -s,  t), node_id++);
            mesh.add_point (Point( 0,  s,  t), node_id++);
            mesh.add_point (Point( 0, -s, -t), node_id++);
            mesh.add_point (Point( 0,  s, -t), node_id++);

            mesh.add_point (Point( t,  0, -s), node_id++);
            mesh.add_point (Point( t,  0,  s), node_id++);
            mesh.add_point (Point(-t,  0, -s), node_id++);
            mesh.add_point (Point(-t,  0,  s), node_id++);

            // Create the 20 triangles of the icosahedron
            static const unsigned int idx1 [6] = {11, 5, 1, 7, 10, 11};
            static const unsigned int idx2 [6] = {9, 4, 2, 6, 8, 9};
            static const unsigned int idx3 [6] = {1, 5, 11, 10, 7, 1};

            for (unsigned int i = 0; i < 5; ++i)
              {
                // 5 elems around point 0
                Elem * new_elem = mesh.add_elem (new Tri3);
                new_elem->set_node(0) = mesh.node_ptr(0);
                new_elem->set_node(1) = mesh.node_ptr(idx1[i]);
                new_elem->set_node(2) = mesh.node_ptr(idx1[i+1]);

                // 5 adjacent elems
                new_elem = mesh.add_elem (new Tri3);
                new_elem->set_node(0) = mesh.node_ptr(idx3[i]);
                new_elem->set_node(1) = mesh.node_ptr(idx3[i+1]);
                new_elem->set_node(2) = mesh.node_ptr(idx2[i]);

                // 5 elems around point 3
                new_elem = mesh.add_elem (new Tri3);
                new_elem->set_node(0) = mesh.node_ptr(3);
                new_elem->set_node(1) = mesh.node_ptr(idx2[i]);
                new_elem->set_node(2) = mesh.node_ptr(idx2[i+1]);

                // 5 adjacent elems
                new_elem = mesh.add_elem (new Tri3);
                new_elem->set_node(0) = mesh.node_ptr(idx2[i+1]);
                new_elem->set_node(1) = mesh.node_ptr(idx2[i]);
                new_elem->set_node(2) = mesh.node_ptr(idx3[i+1]);
              }
          }

        break;
      } // end case 2





      //-----------------------------------------------------------------
      // Build a sphere in three dimensions
    case 3:
      {
        // (Currently) supported types
        if (!((type == HEX8) || (type == HEX27)))
          {
            // FIXME: We'd need an all_tet() routine (which could also be used by
            // build_square()) to do Tets, or Prisms for that matter.
            libmesh_error_msg("Error: Only HEX8/27 currently supported.");
          }


        // 3D analog of 2D initial grid:
        const Real
          r_small = 0.25*rad,                      //  0.25 *radius
          r_med   = (0.125*std::sqrt(2.)+0.5)*rad; // .67677*radius

        // (Temporary) convenient storage for node pointers
        std::vector<Node *> nodes(16);

        // For DistributedMesh, if we don't specify node IDs the Mesh
        // will try to pick an appropriate (unique) one for us.  But
        // since we are adding these nodes on all processors, we want
        // to be sure they have consistent IDs across all processors.
        unsigned node_id = 0;

        // Points 0-7 are the initial HEX8
        nodes[0] = mesh.add_point (Point(-r_small,-r_small, -r_small), node_id++);
        nodes[1] = mesh.add_point (Point( r_small,-r_small, -r_small), node_id++);
        nodes[2] = mesh.add_point (Point( r_small, r_small, -r_small), node_id++);
        nodes[3] = mesh.add_point (Point(-r_small, r_small, -r_small), node_id++);
        nodes[4] = mesh.add_point (Point(-r_small,-r_small,  r_small), node_id++);
        nodes[5] = mesh.add_point (Point( r_small,-r_small,  r_small), node_id++);
        nodes[6] = mesh.add_point (Point( r_small, r_small,  r_small), node_id++);
        nodes[7] = mesh.add_point (Point(-r_small, r_small,  r_small), node_id++);

        //  Points 8-15 are for the outer hexes, we number them in the same way
        nodes[8]  = mesh.add_point (Point(-r_med,-r_med, -r_med), node_id++);
        nodes[9]  = mesh.add_point (Point( r_med,-r_med, -r_med), node_id++);
        nodes[10] = mesh.add_point (Point( r_med, r_med, -r_med), node_id++);
        nodes[11] = mesh.add_point (Point(-r_med, r_med, -r_med), node_id++);
        nodes[12] = mesh.add_point (Point(-r_med,-r_med,  r_med), node_id++);
        nodes[13] = mesh.add_point (Point( r_med,-r_med,  r_med), node_id++);
        nodes[14] = mesh.add_point (Point( r_med, r_med,  r_med), node_id++);
        nodes[15] = mesh.add_point (Point(-r_med, r_med,  r_med), node_id++);

        // Now create the elements and add them to the mesh
        // Element 0 - center element
        {
          Elem * elem0 = mesh.add_elem (new Hex8);
          elem0->set_node(0) = nodes[0];
          elem0->set_node(1) = nodes[1];
          elem0->set_node(2) = nodes[2];
          elem0->set_node(3) = nodes[3];
          elem0->set_node(4) = nodes[4];
          elem0->set_node(5) = nodes[5];
          elem0->set_node(6) = nodes[6];
          elem0->set_node(7) = nodes[7];
        }

        // Element 1 - "bottom"
        {
          Elem * elem1 = mesh.add_elem (new Hex8);
          elem1->set_node(0) = nodes[8];
          elem1->set_node(1) = nodes[9];
          elem1->set_node(2) = nodes[10];
          elem1->set_node(3) = nodes[11];
          elem1->set_node(4) = nodes[0];
          elem1->set_node(5) = nodes[1];
          elem1->set_node(6) = nodes[2];
          elem1->set_node(7) = nodes[3];
        }

        // Element 2 - "front"
        {
          Elem * elem2 = mesh.add_elem (new Hex8);
          elem2->set_node(0) = nodes[8];
          elem2->set_node(1) = nodes[9];
          elem2->set_node(2) = nodes[1];
          elem2->set_node(3) = nodes[0];
          elem2->set_node(4) = nodes[12];
          elem2->set_node(5) = nodes[13];
          elem2->set_node(6) = nodes[5];
          elem2->set_node(7) = nodes[4];
        }

        // Element 3 - "right"
        {
          Elem * elem3 = mesh.add_elem (new Hex8);
          elem3->set_node(0) = nodes[1];
          elem3->set_node(1) = nodes[9];
          elem3->set_node(2) = nodes[10];
          elem3->set_node(3) = nodes[2];
          elem3->set_node(4) = nodes[5];
          elem3->set_node(5) = nodes[13];
          elem3->set_node(6) = nodes[14];
          elem3->set_node(7) = nodes[6];
        }

        // Element 4 - "back"
        {
          Elem * elem4 = mesh.add_elem (new Hex8);
          elem4->set_node(0) = nodes[3];
          elem4->set_node(1) = nodes[2];
          elem4->set_node(2) = nodes[10];
          elem4->set_node(3) = nodes[11];
          elem4->set_node(4) = nodes[7];
          elem4->set_node(5) = nodes[6];
          elem4->set_node(6) = nodes[14];
          elem4->set_node(7) = nodes[15];
        }

        // Element 5 - "left"
        {
          Elem * elem5 = mesh.add_elem (new Hex8);
          elem5->set_node(0) = nodes[8];
          elem5->set_node(1) = nodes[0];
          elem5->set_node(2) = nodes[3];
          elem5->set_node(3) = nodes[11];
          elem5->set_node(4) = nodes[12];
          elem5->set_node(5) = nodes[4];
          elem5->set_node(6) = nodes[7];
          elem5->set_node(7) = nodes[15];
        }

        // Element 6 - "top"
        {
          Elem * elem6 = mesh.add_elem (new Hex8);
          elem6->set_node(0) = nodes[4];
          elem6->set_node(1) = nodes[5];
          elem6->set_node(2) = nodes[6];
          elem6->set_node(3) = nodes[7];
          elem6->set_node(4) = nodes[12];
          elem6->set_node(5) = nodes[13];
          elem6->set_node(6) = nodes[14];
          elem6->set_node(7) = nodes[15];
        }

        break;
      } // end case 3

    default:
      libmesh_error_msg("Unknown dimension " << mesh.mesh_dimension());



    } // end switch (dim)

  // Now we have the beginnings of a sphere.
  // Add some more elements by doing uniform refinements and
  // popping nodes to the boundary.
  MeshRefinement mesh_refinement (mesh);

  // Loop over the elements, refine, pop nodes to boundary.
  for (unsigned int r=0; r<nr; r++)
    {
      mesh_refinement.uniformly_refine(1);

      MeshBase::element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::element_iterator end = mesh.active_elements_end();

      for (; it != end; ++it)
        {
          Elem * elem = *it;

          for (unsigned int s=0; s<elem->n_sides(); s++)
            if (elem->neighbor_ptr(s) == libmesh_nullptr || (mesh.mesh_dimension() == 2 && !flat))
              {
                UniquePtr<Elem> side(elem->build_side_ptr(s));

                // Pop each point to the sphere boundary
                for (unsigned int n=0; n<side->n_nodes(); n++)
                  side->point(n) =
                    sphere.closest_point(side->point(n));
              }
        }
    }

  // The mesh now contains a refinement hierarchy due to the refinements
  // used to generate the grid.  In order to call other support functions
  // like all_tri() and all_second_order, you need a "flat" mesh file (with no
  // refinement trees) so
  MeshTools::Modification::flatten(mesh);

  // In 2D, convert all the quads to triangles if requested
  if (mesh.mesh_dimension()==2)
    {
      if ((type == TRI6) || (type == TRI3))
        {
          MeshTools::Modification::all_tri(mesh);
        }
    }


  // Convert to second-order elements if the user requested it.
  if (Elem::build(type)->default_order() != FIRST)
    {
      // type is already a second-order, determine if it is the
      // "full-ordered" second-order element, or the "serendipity"
      // second order element.  Note also that all_second_order
      // can't be called once the mesh has been refined.
      bool full_ordered = !((type==QUAD8) || (type==HEX20));
      mesh.all_second_order(full_ordered);

      // And pop to the boundary again...
      MeshBase::element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::element_iterator end = mesh.active_elements_end();

      for (; it != end; ++it)
        {
          Elem * elem = *it;

          for (unsigned int s=0; s<elem->n_sides(); s++)
            if (elem->neighbor_ptr(s) == libmesh_nullptr)
              {
                UniquePtr<Elem> side(elem->build_side_ptr(s));

                // Pop each point to the sphere boundary
                for (unsigned int n=0; n<side->n_nodes(); n++)
                  side->point(n) =
                    sphere.closest_point(side->point(n));
              }
        }
    }


  // The meshes could probably use some smoothing.
  LaplaceMeshSmoother smoother(mesh);
  smoother.smooth(n_smooth);

  // We'll give the whole sphere surface a boundary id of 0
  {
    MeshBase::element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::element_iterator end = mesh.active_elements_end();

    for (; it != end; ++it)
      {
        Elem * elem = *it;
        for (unsigned short s=0; s != elem->n_sides(); ++s)
          if (!elem->neighbor_ptr(s))
            boundary_info.add_side(elem, s, 0);
      }
  }

  STOP_LOG("build_sphere()", "MeshTools::Generation");


  // Done building the mesh.  Now prepare it for use.
  mesh.prepare_for_use(/*skip_renumber =*/ false);
}

#endif // #ifndef LIBMESH_ENABLE_AMR


// Meshes the tensor product of a 1D and a 1D-or-2D domain.
void MeshTools::Generation::build_extrusion (UnstructuredMesh & mesh,
                                             const MeshBase & cross_section,
                                             const unsigned int nz,
                                             RealVectorValue extrusion_vector,
                                             QueryElemSubdomainIDBase * elem_subdomain)
{
  if (!cross_section.n_elem())
    return;

  START_LOG("build_extrusion()", "MeshTools::Generation");

  dof_id_type orig_elem = cross_section.n_elem();
  dof_id_type orig_nodes = cross_section.n_nodes();

  unsigned int order = 1;

  BoundaryInfo & boundary_info = mesh.get_boundary_info();
  const BoundaryInfo & cross_section_boundary_info = cross_section.get_boundary_info();

  // If cross_section is distributed, so is its extrusion
  if (!cross_section.is_serial())
    mesh.delete_remote_elements();

  // We know a priori how many elements we'll need
  mesh.reserve_elem(nz*orig_elem);

  // For straightforward meshes we need one or two additional layers per
  // element.
  if ((*cross_section.elements_begin())->default_order() == SECOND)
    order = 2;

  mesh.reserve_nodes((order*nz+1)*orig_nodes);

  // Container to catch the boundary IDs handed back by the BoundaryInfo object
  std::vector<boundary_id_type> ids_to_copy;

  MeshBase::const_node_iterator       nd  = cross_section.nodes_begin();
  const MeshBase::const_node_iterator nend = cross_section.nodes_end();
  for (; nd!=nend; ++nd)
    {
      const Node * node = *nd;

      for (unsigned int k=0; k != order*nz+1; ++k)
        {
          Node * new_node =
            mesh.add_point(*node +
                           (extrusion_vector * k / nz / order),
                           node->id() + (k * orig_nodes),
                           node->processor_id());

          cross_section_boundary_info.boundary_ids(node, ids_to_copy);
          boundary_info.add_node(new_node, ids_to_copy);
        }
    }

  const std::set<boundary_id_type> & side_ids =
    cross_section_boundary_info.get_side_boundary_ids();
  const boundary_id_type next_side_id = side_ids.empty() ?
    0 : cast_int<boundary_id_type>(*side_ids.rbegin() + 1);

  MeshBase::const_element_iterator       el  = cross_section.elements_begin();
  const MeshBase::const_element_iterator end = cross_section.elements_end();
  for (; el!=end; ++el)
    {
      const Elem * elem = *el;
      const ElemType etype = elem->type();

      // build_extrusion currently only works on coarse meshes
      libmesh_assert (!elem->parent());

      for (unsigned int k=0; k != nz; ++k)
        {
          Elem * new_elem;
          switch (etype)
            {
            case EDGE2:
              {
                new_elem = new Quad4;
                new_elem->set_node(0) = mesh.node_ptr(elem->node_ptr(0)->id() + (k * orig_nodes));
                new_elem->set_node(1) = mesh.node_ptr(elem->node_ptr(1)->id() + (k * orig_nodes));
                new_elem->set_node(2) = mesh.node_ptr(elem->node_ptr(1)->id() + ((k+1) * orig_nodes));
                new_elem->set_node(3) = mesh.node_ptr(elem->node_ptr(0)->id() + ((k+1) * orig_nodes));
                break;
              }
            case EDGE3:
              {
                new_elem = new Quad9;
                new_elem->set_node(0) = mesh.node_ptr(elem->node_ptr(0)->id() + (2*k * orig_nodes));
                new_elem->set_node(1) = mesh.node_ptr(elem->node_ptr(1)->id() + (2*k * orig_nodes));
                new_elem->set_node(2) = mesh.node_ptr(elem->node_ptr(1)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(3) = mesh.node_ptr(elem->node_ptr(0)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(4) = mesh.node_ptr(elem->node_ptr(2)->id() + (2*k * orig_nodes));
                new_elem->set_node(5) = mesh.node_ptr(elem->node_ptr(1)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(6) = mesh.node_ptr(elem->node_ptr(2)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(7) = mesh.node_ptr(elem->node_ptr(0)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(8) = mesh.node_ptr(elem->node_ptr(2)->id() + ((2*k+1) * orig_nodes));
                break;
              }
            case TRI3:
              {
                new_elem = new Prism6;
                new_elem->set_node(0) = mesh.node_ptr(elem->node_ptr(0)->id() + (k * orig_nodes));
                new_elem->set_node(1) = mesh.node_ptr(elem->node_ptr(1)->id() + (k * orig_nodes));
                new_elem->set_node(2) = mesh.node_ptr(elem->node_ptr(2)->id() + (k * orig_nodes));
                new_elem->set_node(3) = mesh.node_ptr(elem->node_ptr(0)->id() + ((k+1) * orig_nodes));
                new_elem->set_node(4) = mesh.node_ptr(elem->node_ptr(1)->id() + ((k+1) * orig_nodes));
                new_elem->set_node(5) = mesh.node_ptr(elem->node_ptr(2)->id() + ((k+1) * orig_nodes));
                break;
              }
            case TRI6:
              {
                new_elem = new Prism18;
                new_elem->set_node(0) = mesh.node_ptr(elem->node_ptr(0)->id() + (2*k * orig_nodes));
                new_elem->set_node(1) = mesh.node_ptr(elem->node_ptr(1)->id() + (2*k * orig_nodes));
                new_elem->set_node(2) = mesh.node_ptr(elem->node_ptr(2)->id() + (2*k * orig_nodes));
                new_elem->set_node(3) = mesh.node_ptr(elem->node_ptr(0)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(4) = mesh.node_ptr(elem->node_ptr(1)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(5) = mesh.node_ptr(elem->node_ptr(2)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(6) = mesh.node_ptr(elem->node_ptr(3)->id() + (2*k * orig_nodes));
                new_elem->set_node(7) = mesh.node_ptr(elem->node_ptr(4)->id() + (2*k * orig_nodes));
                new_elem->set_node(8) = mesh.node_ptr(elem->node_ptr(5)->id() + (2*k * orig_nodes));
                new_elem->set_node(9) = mesh.node_ptr(elem->node_ptr(0)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(10) = mesh.node_ptr(elem->node_ptr(1)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(11) = mesh.node_ptr(elem->node_ptr(2)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(12) = mesh.node_ptr(elem->node_ptr(3)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(13) = mesh.node_ptr(elem->node_ptr(4)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(14) = mesh.node_ptr(elem->node_ptr(5)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(15) = mesh.node_ptr(elem->node_ptr(3)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(16) = mesh.node_ptr(elem->node_ptr(4)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(17) = mesh.node_ptr(elem->node_ptr(5)->id() + ((2*k+1) * orig_nodes));
                break;
              }
            case QUAD4:
              {
                new_elem = new Hex8;
                new_elem->set_node(0) = mesh.node_ptr(elem->node_ptr(0)->id() + (k * orig_nodes));
                new_elem->set_node(1) = mesh.node_ptr(elem->node_ptr(1)->id() + (k * orig_nodes));
                new_elem->set_node(2) = mesh.node_ptr(elem->node_ptr(2)->id() + (k * orig_nodes));
                new_elem->set_node(3) = mesh.node_ptr(elem->node_ptr(3)->id() + (k * orig_nodes));
                new_elem->set_node(4) = mesh.node_ptr(elem->node_ptr(0)->id() + ((k+1) * orig_nodes));
                new_elem->set_node(5) = mesh.node_ptr(elem->node_ptr(1)->id() + ((k+1) * orig_nodes));
                new_elem->set_node(6) = mesh.node_ptr(elem->node_ptr(2)->id() + ((k+1) * orig_nodes));
                new_elem->set_node(7) = mesh.node_ptr(elem->node_ptr(3)->id() + ((k+1) * orig_nodes));
                break;
              }
            case QUAD9:
              {
                new_elem = new Hex27;
                new_elem->set_node(0) = mesh.node_ptr(elem->node_ptr(0)->id() + (2*k * orig_nodes));
                new_elem->set_node(1) = mesh.node_ptr(elem->node_ptr(1)->id() + (2*k * orig_nodes));
                new_elem->set_node(2) = mesh.node_ptr(elem->node_ptr(2)->id() + (2*k * orig_nodes));
                new_elem->set_node(3) = mesh.node_ptr(elem->node_ptr(3)->id() + (2*k * orig_nodes));
                new_elem->set_node(4) = mesh.node_ptr(elem->node_ptr(0)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(5) = mesh.node_ptr(elem->node_ptr(1)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(6) = mesh.node_ptr(elem->node_ptr(2)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(7) = mesh.node_ptr(elem->node_ptr(3)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(8) = mesh.node_ptr(elem->node_ptr(4)->id() + (2*k * orig_nodes));
                new_elem->set_node(9) = mesh.node_ptr(elem->node_ptr(5)->id() + (2*k * orig_nodes));
                new_elem->set_node(10) = mesh.node_ptr(elem->node_ptr(6)->id() + (2*k * orig_nodes));
                new_elem->set_node(11) = mesh.node_ptr(elem->node_ptr(7)->id() + (2*k * orig_nodes));
                new_elem->set_node(12) = mesh.node_ptr(elem->node_ptr(0)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(13) = mesh.node_ptr(elem->node_ptr(1)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(14) = mesh.node_ptr(elem->node_ptr(2)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(15) = mesh.node_ptr(elem->node_ptr(3)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(16) = mesh.node_ptr(elem->node_ptr(4)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(17) = mesh.node_ptr(elem->node_ptr(5)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(18) = mesh.node_ptr(elem->node_ptr(6)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(19) = mesh.node_ptr(elem->node_ptr(7)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(20) = mesh.node_ptr(elem->node_ptr(8)->id() + (2*k * orig_nodes));
                new_elem->set_node(21) = mesh.node_ptr(elem->node_ptr(4)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(22) = mesh.node_ptr(elem->node_ptr(5)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(23) = mesh.node_ptr(elem->node_ptr(6)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(24) = mesh.node_ptr(elem->node_ptr(7)->id() + ((2*k+1) * orig_nodes));
                new_elem->set_node(25) = mesh.node_ptr(elem->node_ptr(8)->id() + ((2*k+2) * orig_nodes));
                new_elem->set_node(26) = mesh.node_ptr(elem->node_ptr(8)->id() + ((2*k+1) * orig_nodes));
                break;
              }
            default:
              {
                libmesh_not_implemented();
                break;
              }
            }

          new_elem->set_id(elem->id() + (k * orig_elem));
          new_elem->processor_id() = elem->processor_id();

          if (!elem_subdomain)
            // maintain the subdomain_id
            new_elem->subdomain_id() = elem->subdomain_id();
          else
            // Allow the user to choose new subdomain_ids
            new_elem->subdomain_id() = elem_subdomain->get_subdomain_for_layer(elem, k);

          new_elem = mesh.add_elem(new_elem);

          // Copy any old boundary ids on all sides
          for (unsigned short s = 0; s != elem->n_sides(); ++s)
            {
              cross_section_boundary_info.boundary_ids(elem, s, ids_to_copy);

              if (new_elem->dim() == 3)
                {
                  // For 2D->3D extrusion, we give the boundary IDs
                  // for side s on the old element to side s+1 on the
                  // new element.  This is just a happy coincidence as
                  // far as I can tell...
                  boundary_info.add_side(new_elem,
                                         cast_int<unsigned short>(s+1),
                                         ids_to_copy);
                }
              else
                {
                  // For 1D->2D extrusion, the boundary IDs map as:
                  // Old elem -> New elem
                  // 0        -> 3
                  // 1        -> 1
                  libmesh_assert_less(s, 2);
                  const unsigned short sidemap[2] = {3, 1};
                  boundary_info.add_side(new_elem, sidemap[s], ids_to_copy);
                }
            }

          // Give new boundary ids to bottom and top
          if (k == 0)
            boundary_info.add_side(new_elem, 0, next_side_id);
          if (k == nz-1)
            {
              // For 2D->3D extrusion, the "top" ID is 1+the original
              // element's number of sides.  For 1D->2D extrusion, the
              // "top" ID is side 2.
              const unsigned short top_id = new_elem->dim() == 3 ?
                cast_int<unsigned short>(elem->n_sides()+1) : 2;
              boundary_info.add_side
                (new_elem, top_id,
                 cast_int<boundary_id_type>(next_side_id+1));
            }
        }
    }

  STOP_LOG("build_extrusion()", "MeshTools::Generation");

  // Done building the mesh.  Now prepare it for use.
  mesh.prepare_for_use(/*skip_renumber =*/ false);
}




#ifdef LIBMESH_HAVE_TRIANGLE

// Triangulates a 2D rectangular region with or without holes
void MeshTools::Generation::build_delaunay_square(UnstructuredMesh & mesh,
                                                  const unsigned int nx, // num. of elements in x-dir
                                                  const unsigned int ny, // num. of elements in y-dir
                                                  const Real xmin, const Real xmax,
                                                  const Real ymin, const Real ymax,
                                                  const ElemType type,
                                                  const std::vector<TriangleInterface::Hole*> * holes)
{
  // Check for reasonable size
  libmesh_assert_greater_equal (nx, 1); // need at least 1 element in x-direction
  libmesh_assert_greater_equal (ny, 1); // need at least 1 element in y-direction
  libmesh_assert_less (xmin, xmax);
  libmesh_assert_less (ymin, ymax);

  // Clear out any data which may have been in the Mesh
  mesh.clear();

  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  // Make sure the new Mesh will be 2D
  mesh.set_mesh_dimension(2);

  // The x and y spacing between boundary points
  const Real delta_x = (xmax-xmin) / static_cast<Real>(nx);
  const Real delta_y = (ymax-ymin) / static_cast<Real>(ny);

  // Bottom
  for (unsigned int p=0; p<=nx; ++p)
    mesh.add_point(Point(xmin + p*delta_x, ymin));

  // Right side
  for (unsigned int p=1; p<ny; ++p)
    mesh.add_point(Point(xmax, ymin + p*delta_y));

  // Top
  for (unsigned int p=0; p<=nx; ++p)
    mesh.add_point(Point(xmax - p*delta_x, ymax));

  // Left side
  for (unsigned int p=1; p<ny; ++p)
    mesh.add_point(Point(xmin,  ymax - p*delta_y));

  // Be sure we added as many points as we thought we did
  libmesh_assert_equal_to (mesh.n_nodes(), 2*(nx+ny));

  // Construct the Triangle Interface object
  TriangleInterface t(mesh);

  // Set custom variables for the triangulation
  t.desired_area()       = 0.5 * (xmax-xmin)*(ymax-ymin) / static_cast<Real>(nx*ny);
  t.triangulation_type() = TriangleInterface::PSLG;
  t.elem_type()          = type;

  if (holes != libmesh_nullptr)
    t.attach_hole_list(holes);

  // Triangulate!
  t.triangulate();

  // The mesh is now generated, but we still need to mark the boundaries
  // to be consistent with the other build_square routines.  Note that all
  // hole boundary elements get the same ID, 4.
  MeshBase::element_iterator       el     = mesh.elements_begin();
  const MeshBase::element_iterator end_el = mesh.elements_end();
  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      for (unsigned int s=0; s<elem->n_sides(); s++)
        if (elem->neighbor_ptr(s) == libmesh_nullptr)
          {
            UniquePtr<const Elem> side (elem->build_side_ptr(s));

            // Check the location of the side's midpoint.  Since
            // the square has straight sides, the midpoint is not
            // on the corner and thus it is uniquely on one of the
            // sides.
            Point side_midpoint= 0.5f*( side->point(0) + side->point(1) );

            // The boundary ids are set following the same convention as Quad4 sides
            // bottom = 0
            // right  = 1
            // top = 2
            // left = 3
            // hole = 4
            boundary_id_type bc_id=4;

            // bottom
            if      (std::fabs(side_midpoint(1) - ymin) < TOLERANCE)
              bc_id=0;

            // right
            else if (std::fabs(side_midpoint(0) - xmax) < TOLERANCE)
              bc_id=1;

            // top
            else if (std::fabs(side_midpoint(1) - ymax) < TOLERANCE)
              bc_id=2;

            // left
            else if (std::fabs(side_midpoint(0) - xmin) < TOLERANCE)
              bc_id=3;

            // If the point is not on any of the external boundaries, it
            // is on one of the holes....

            // Finally, add this element's information to the boundary info object.
            boundary_info.add_side(elem->id(), s, bc_id);
          }
    }

} // end build_delaunay_square

#endif // LIBMESH_HAVE_TRIANGLE



} // namespace libMesh
