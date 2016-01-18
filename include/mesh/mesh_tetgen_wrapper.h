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

#ifndef LIBMESH_MESH_TETGEN_WRAPPER_H
#define LIBMESH_MESH_TETGEN_WRAPPER_H

#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_TETGEN

// Local includes

// TetGen include file
#include "tetgen.h"  // Defines REAL and other Tetgen types

// C++ includes
#include <string>

namespace libMesh
{
/**
 * The \p TetGenWrapper provides an interface for basic
 * access to TetGen data structures and methods.
 *
 * \author Steffen Petersen
 * \date 2004
 * \author John W. Peterson
 * \date 2011
 */
class TetGenWrapper
{
public:

  /**
   * Constructor.
   */
  TetGenWrapper ();

  /**
   * Destructor.  Empty.
   */
  ~TetGenWrapper ();

  /**
   * Method to set TetGen commandline switches
   * -p Tetrahedralizes a piecewise linear complex (.poly or .smesh file).
   * -q Quality mesh generation. A minimum radius-edge ratio may be specified (default 2.0).
   * -a Applies a maximum tetrahedron volume constraint.
   * -A Assigns attributes to identify tetrahedra in certain regions.
   * -r Reconstructs and Refines a previously generated mesh.
   * -Y Suppresses boundary facets/segments splitting.
   * -i Inserts a list of additional points into mesh.
   * -M Does not merge coplanar facets.
   * -T Set a tolerance for coplanar test (default 1e-8).
   * -d Detect intersections of PLC facets.
   * -z Numbers all output items starting from zero.
   * -o2 Generates second-order subparametric elements.
   * -f Outputs faces (including non-boundary faces) to .face file.
   * -e Outputs subsegments to .edge file.
   * -n Outputs tetrahedra neighbors to .neigh file.
   * -g Outputs mesh to .mesh file for viewing by Medit.
   * -G Outputs mesh to .msh file for viewing by Gid.
   * -O Outputs mesh to .off file for viewing by Geomview.
   * -J No jettison of unused vertices from output .node file.
   * -B Suppresses output of boundary information.
   * -N Suppresses output of .node file.
   * -E Suppresses output of .ele file.
   * -F Suppresses output of .face file.
   * -I Suppresses mesh iteration numbers.
   * -C Checks the consistency of the final mesh.
   * -Q Quiet: No terminal output except errors.
   * -V Verbose: Detailed information, more terminal output.
   * -v Prints the version information.
   * -h Help: A brief instruction for using TetGen.
   */
  void set_switches(const std::string & s);

  /**
   * Method starts triangulization.
   */
  void run_tetgen();

  /**
   * Method returns number of tetrahedra in TetGen output.
   */
  int  get_numberoftetrahedra();

  /**
   * Method returns number of triangle surface elts. in TetGen output.
   */
  int  get_numberoftrifaces();

  /**
   * Method sets number of nodes in TetGen input.
   */
  void set_numberofpoints(int i);

  /**
   * Method returns number of nodes in TetGen output.
   */
  int get_numberofpoints();

  /**
   * Method sets number of facets in TetGen input.
   */
  void set_numberoffacets(int i);

  /**
   * Method sets number of holes in TetGen input.
   */
  void set_numberofholes(int i);

  /**
   * Method sets number of regions in TetGen input.
   */
  void set_numberofregions(int i);

  /**
   * Method allocates memory, sets number of nodes in TetGen input.
   */
  void allocate_pointlist(int numofpoints);

  /**
   * Method allocates memory, sets number of facets, holes in TetGen input.
   */
  void allocate_facetlist(int numoffacets, int numofholes);

  /**
   * Method allocates memory, sets number of regions in TetGen input.
   */
  void allocate_regionlist(int numofregions);

  /**
   * Method sets coordinates of point i in TetGen input.
   */
  void set_node(unsigned i, REAL x, REAL y, REAL z);

  /**
   * Method returns coordinates of point i in TetGen output.
   */
  void get_output_node(unsigned i, REAL & x, REAL & y, REAL & z);

  /**
   * Method returns index of jth node from element i in TetGen output.
   */
  int  get_element_node(unsigned i, unsigned j);

  /**
   * Method returns index of jth node from surface triangle i in TetGen output.
   */
  int  get_triface_node(unsigned i, unsigned j);

  /**
   * Method returns attribute of element i in TetGen output.
   */
  REAL get_element_attribute(unsigned i);

  /**
   * Method sets coordinates of hole i in TetGen input.
   */
  void set_hole(unsigned i, REAL x, REAL y, REAL z);

  /**
   * Method sets number of polygons for facet i in TetGen input.
   */
  void set_facet_numberofpolygons(unsigned i, int num);

  /**
   * Method sets number of holes for facet i in TetGen input.
   */
  void set_facet_numberofholes(unsigned i, int num);

  /**
   * Method allocates memory, sets number of polygons for facet i
   * in TetGen input.
   */
  void allocate_facet_polygonlist(unsigned i, int numofpolygons);

  /**
   * Method sets number of vertices for polygon j, facet i in TetGen input.
   */
  void set_polygon_numberofvertices(unsigned i, unsigned j, int num);

  /**
   * Method allocates memory, sets number of vertices for polygon j,
   * facet i in TetGen input.
   */
  void allocate_polygon_vertexlist(unsigned i, unsigned j, int numofvertices);

  /**
   * Method sets index of ith facet, jth polygon, kth vertex in
   * TetGen input.
   */
  void set_vertex(unsigned i, unsigned j, unsigned k, int nodeindex);

  /**
   * Method sets coordinates, attribute and volume constraint for
   * region i in TetGen input.  Note that coordinates and attributes
   * will only be considered if the corresponding switches are
   * enabled.  See TetGen documentation for more details.
   */
  void set_region(unsigned i, REAL x, REAL y, REAL z,
                  REAL attribute, REAL vol_constraint);

  /**
   * TetGen input structure.
   */
  tetgenio   tetgen_data;

  /**
   * TetGen output structure.
   */
  tetgenio *  tetgen_output;

  /**
   * TetGen mesh structure (from the TetGen library).
   */
  tetgenmesh      tetgen_mesh;

  /**
   * TetGen control class (from the TetGen library).
   */
  tetgenbehavior  tetgen_be;
};



} // namespace libMesh


#endif // LIBMESH_HAVE_TETGEN
#endif // LIBMESH_MESH_TETGEN_WRAPPER_H
