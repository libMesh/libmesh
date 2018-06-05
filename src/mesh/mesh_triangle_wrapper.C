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


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_TRIANGLE

// Local includes
#include "libmesh/mesh_triangle_wrapper.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/point.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/enum_elem_type.h"

namespace libMesh
{

void TriangleWrapper::init(TriangleWrapper::triangulateio & t)
{
  t.pointlist                    = static_cast<REAL*>(libmesh_nullptr);
  t.pointattributelist           = static_cast<REAL*>(libmesh_nullptr);
  t.pointmarkerlist              = static_cast<int *>(libmesh_nullptr);
  t.numberofpoints               = 0 ;
  t.numberofpointattributes      = 0 ;

  t.trianglelist                 = static_cast<int *>(libmesh_nullptr);
  t.triangleattributelist        = static_cast<REAL*>(libmesh_nullptr);
  t.trianglearealist             = static_cast<REAL*>(libmesh_nullptr);
  t.neighborlist                 = static_cast<int *>(libmesh_nullptr);
  t.numberoftriangles            = 0;
  t.numberofcorners              = 0;
  t.numberoftriangleattributes   = 0;

  t.segmentlist                  = static_cast<int *>(libmesh_nullptr);
  t.segmentmarkerlist            = static_cast<int *>(libmesh_nullptr);
  t.numberofsegments             = 0;

  t.holelist                     = static_cast<REAL*>(libmesh_nullptr);
  t.numberofholes                = 0;

  t.regionlist                   = static_cast<REAL*>(libmesh_nullptr);
  t.numberofregions              = 0;

  t.edgelist                     = static_cast<int *>(libmesh_nullptr);
  t.edgemarkerlist               = static_cast<int *>(libmesh_nullptr);
  t.normlist                     = static_cast<REAL*>(libmesh_nullptr);
  t.numberofedges                = 0;
}






void TriangleWrapper::destroy(TriangleWrapper::triangulateio & t, TriangleWrapper::IO_Type io_type)
{
  std::free (t.pointlist            );
  std::free (t.pointattributelist   );
  std::free (t.pointmarkerlist      );
  std::free (t.trianglelist         );
  std::free (t.triangleattributelist);
  std::free (t.trianglearealist     );
  std::free (t.neighborlist         );
  std::free (t.segmentlist          );
  std::free (t.segmentmarkerlist    );

  // Only attempt to free these when t was used as an input struct!
  if (io_type==INPUT)
    {
      std::free (t.holelist  );
      std::free (t.regionlist);
    }

  std::free (t.edgelist      );
  std::free (t.edgemarkerlist);
  std::free (t.normlist      );

  // Reset
  // TriangleWrapper::init(t);
}






void TriangleWrapper::copy_tri_to_mesh(const triangulateio & triangle_data_input,
                                       UnstructuredMesh & mesh_output,
                                       const ElemType type)
{
  // Transfer the information into the LibMesh mesh.
  mesh_output.clear();

  // Make sure the new Mesh will be 2D
  mesh_output.set_mesh_dimension(2);

  // Node information
  for (int i=0, c=0; c<triangle_data_input.numberofpoints; i+=2, ++c)
    {
      // Specify ID when adding point, otherwise, if this is DistributedMesh,
      // it might add points with a non-sequential numbering...
      mesh_output.add_point( Point(triangle_data_input.pointlist[i],
                                   triangle_data_input.pointlist[i+1]),
                             /*id=*/c);
    }

  // Element information
  for (int i=0; i<triangle_data_input.numberoftriangles; ++i)
    {
      switch (type)
        {
        case TRI3:
          {
            Elem * elem = mesh_output.add_elem (new Tri3);

            for (unsigned int n=0; n<3; ++n)
              elem->set_node(n) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*3 + n]);

            break;
          }

        case TRI6:
          {
            Elem * elem = mesh_output.add_elem (new Tri6);

            // Triangle number TRI6 nodes in a different way to libMesh
            elem->set_node(0) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 0]);
            elem->set_node(1) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 1]);
            elem->set_node(2) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 2]);
            elem->set_node(3) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 5]);
            elem->set_node(4) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 3]);
            elem->set_node(5) = mesh_output.node_ptr(triangle_data_input.trianglelist[i*6 + 4]);

            break;
          }

        default:
          libmesh_error_msg("ERROR: Unrecognized triangular element type.");
        }
    }

  // Note: If the input mesh was a parallel one, calling
  // prepare_for_use() now will re-parallelize it by a call to
  // delete_remote_elements()... We do not actually want to
  // reparallelize it here though: the triangulate() function may
  // still do some Mesh smoothing.  The main thing needed (for
  // smoothing) is the neighbor information, so let's just find
  // neighbors...
  //mesh_output.prepare_for_use(/*skip_renumber =*/false);
  mesh_output.find_neighbors();
}


}

#endif // LIBMESH_HAVE_TRIANGLE
