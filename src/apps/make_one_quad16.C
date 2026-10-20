// make_one_quad16.C
//
// Writes a single-element QUAD16 mesh on [-1,1]^2 to one_quad16.xda,
// mirroring reference_elements/2D/one_quad9.xda (nodes in local order,
// identity connectivity, sides 0-3 tagged with boundary ids 0-3).
//
// Build against your libMesh, e.g. (from build dir, dbg method):
//   $(libmesh-config --cxx) make_one_quad16.C -o make_one_quad16 \
//       $(libmesh-config --cppflags --cxxflags --include --libs)
// or drop it in src/apps/ style and build with the usual app rules.
//
// Run:  ./make_one_quad16     -> one_quad16.xda

#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"
#include "libmesh/point.h"
#include "libmesh/boundary_info.h"
#include "libmesh/enum_elem_type.h"

#include <array>
#include <cmath>
#include <iostream>

using namespace libMesh;

int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  ReplicatedMesh mesh(init.comm(), /*dim=*/2);

  // Cubic Lagrange nodes sit at -1, -1/3, +1/3, +1 along each axis.
  const Real a = 1.0 / 3.0;

  // Local node order (must match Quad16's numbering / build_cube):
  //   0-3   corners (counter-clockwise from bottom-left)
  //   4,5   bottom edge, left -> right
  //   6,7   right edge,  bottom -> top
  //   8,9   top edge,    right -> left
  //   10,11 left edge,   top -> bottom
  //   12-15 interior: (-a,-a), (a,-a), (a,a), (-a,a)
  const std::array<Point, 16> pts =
    {{
      Point(-1, -1), Point( 1, -1), Point( 1,  1), Point(-1,  1),
      Point(-a, -1), Point( a, -1),
      Point( 1, -a), Point( 1,  a),
      Point( a,  1), Point(-a,  1),
      Point(-1,  a), Point(-1, -a),
      Point(-a, -a), Point( a, -a), Point( a,  a), Point(-a,  a)
    }};

  Elem * elem = mesh.add_elem(Elem::build_with_id(QUAD16, 0));

  if (elem->n_nodes() != 16)
    libmesh_error_msg("QUAD16 reports n_nodes() = " << elem->n_nodes()
                      << ", expected 16");

  for (unsigned int i = 0; i < 16; ++i)
    {
      Node * n = mesh.add_point(pts[i], /*id=*/i);
      elem->set_node(i, n);
    }

  // Sanity check against the class's own reference coordinates.
  // If master_point() is not implemented yet this may throw; comment out then.
  bool mismatch = false;
  for (unsigned int i = 0; i < 16; ++i)
    {
      const Point mp = elem->master_point(i);
      if ((mp - pts[i]).norm() > 1e-12)
        {
          libMesh::err << "Node " << i << ": master_point = " << mp
                       << " but file will use " << pts[i] << std::endl;
          mismatch = true;
        }
    }
  if (mismatch)
    libMesh::err << "WARNING: node order differs from Quad16::master_point()."
                 << " Fix one or the other before using this file.\n";

  // Same boundary convention as one_quad9.xda: side s gets id s.
  BoundaryInfo & bi = mesh.get_boundary_info();
  for (unsigned short s = 0; s < 4; ++s)
    bi.add_side(elem, s, s);

  mesh.prepare_for_use();
  mesh.write("one_quad16.xda");

  libMesh::out << "Wrote one_quad16.xda ("
               << mesh.n_nodes() << " nodes, "
               << mesh.n_elem() << " element)\n";
  return 0;
}
