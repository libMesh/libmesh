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


#ifndef LIBMESH_MESH_TETGEN_INTERFACE_H
#define LIBMESH_MESH_TETGEN_INTERFACE_H

#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_TETGEN


// Local includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/point.h" // used for specifying holes

// C++ includes
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace libMesh
{
// Forward Declarations
class UnstructuredMesh;
class TetGenWrapper;
class Elem;

/**
 * Class \p TetGenMeshInterface provides an interface for
 * tetrahedrization of meshes using the TetGen library.  For
 * information about TetGen cf.
 * <a href="http://tetgen.org/">TetGen home page</a>.
 *
 * \author Steffen Petersen
 * \date 2004
 * \author John W. Peterson
 * \date 2011
 */
class TetGenMeshInterface
{
public:

  /**
   * Constructor. Takes a reference to the mesh.
   */
  explicit
  TetGenMeshInterface (UnstructuredMesh & mesh);

  /**
   * Empty destructor.
   */
  ~TetGenMeshInterface() {}

  /**
   * Method to set swithes to tetgen, allowing for different behaviours
   */
  void set_switches(const std::string &);

  /**
   * Method invokes TetGen library to compute a Delaunay tetrahedrization
   * from the nodes point set.
   */
  void triangulate_pointset ();

  /**
   * Method invokes TetGen library to compute a Delaunay tetrahedrization
   * from the nodes point set. Stores only 2D hull surface elements.
   */
  void pointset_convexhull ();

  /**
   * Method invokes TetGen library to compute a Delaunay tetrahedrization
   * from the nodes point set. Boundary constraints are taken from
   * elements array.
   */
  void triangulate_conformingDelaunayMesh (double quality_constraint=0.,
                                           double volume_constraint=0.);

  /**
   * Method invokes TetGen library to compute a Delaunay tetrahedrization
   * from the nodes point set. Boundary constraints are taken from
   * elements array. Include carve-out functionality.
   */
  void triangulate_conformingDelaunayMesh_carvehole (const std::vector<Point> & holes,
                                                     double quality_constraint=0.,
                                                     double volume_constraint=0.);



protected:
  /**
   * This function copies nodes from the _mesh into TetGen's
   * pointlist.  Takes some pains to ensure that non-sequential
   * node numberings (which can happen with e.g. DistributedMesh)
   * are handled.
   */
  void fill_pointlist(TetGenWrapper & wrapper);

  /**
   * Assigns the node IDs contained in the 'node_labels'
   * array to 'elem'.
   */
  void assign_nodes_to_elem(unsigned * node_labels, Elem * elem);

  /**
   * This function checks the integrity of the current set of
   * elements in the Mesh to see if they comprise a convex hull,
   * that is:
   * .) If they are all TRI3 elements
   * .) They all have non-NULL neighbors
   *
   * Returns:
   * .) 0 if the mesh forms a valid convex hull
   * .) 1 if a non-TRI3 element is found
   * .) 2 if an element with a NULL-neighbor is found
   */
  unsigned check_hull_integrity();

  /**
   * This function prints an informative message and
   * crashes based on the output of the check_hull_integrity()
   * function.  It is a separate function so that you
   * can check hull integrity without crashing if you desire.
   */
  void process_hull_integrity_result(unsigned result);

  /**
   * Delete original convex hull elements from the Mesh
   * after performing a Delaunay tetrahedralization.
   */
  void delete_2D_hull_elements();

  /**
   * Local reference to the mesh we are working with.
   */
  UnstructuredMesh & _mesh;

  /**
   * We should not assume libmesh nodes are numbered sequentially...
   * This is not the default behavior of DistributedMesh, for example,
   * unless you specify node IDs explicitly.  So this array allows us
   * to keep a mapping between the sequential numbering in
   * tetgen_data.pointlist.
   */
  std::vector<unsigned> _sequential_to_libmesh_node_map;

  /**
   * Tetgen only operates on serial meshes.
   */
  MeshSerializer _serializer;

  /**
   * Parameter controling the behaviour of tetgen.
   * By default quiet.
   */
  std::string _switches;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_TETGEN

#endif // LIBMESH_MESH_TETGEN_INTERFACE_H
