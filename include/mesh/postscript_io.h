// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_POSTSCRIPT_IO_H
#define LIBMESH_POSTSCRIPT_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_output.h"
//#include "libmesh/dense_matrix.h"
#include "libmesh/point.h"

// C++ includes
#include <fstream>
#include <sstream>
#include <vector>

namespace libMesh
{

// Forward declarations
class MeshBase;
class Elem;

/**
 * This class implements writing 2D meshes in Postscript.  It borrows
 * several ideas from, and is a more simple-minded version of, the
 * DataOutBase::write_eps() function from Deal II.  Only output is
 * supported here, and only the Mesh (none of the data) is written.
 * The main use I imagined for this class is creating nice Mesh
 * images for publications, since I didn't find/don't know of a free
 * visualization program which would do this.
 *
 * @author John W. Peterson, 2008
 */
class PostscriptIO : public MeshOutput<MeshBase>
{
 public:
  /**
   * Constructor.
   */
  explicit
  PostscriptIO (const MeshBase& mesh);

  /**
   * Destructor.
   */
  virtual ~PostscriptIO ();

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string& );

  /**
   * Controls greyscale shading of cells.  By default this value
   * is 0.0 (which actually corresponds to black) and this indicates
   * "no shading" i.e. only mesh lines will be drawn.  Any other
   * value in (0,1] will cause the cells to be grey-shaded to some
   * degree, with higher values being lighter.  A value of 0.75
   * gives decent results.
   */
  Real shade_value;

  /**
   * Control the thickness of the lines used.  0.5 is a reasonable default
   * for printed images, but you may need to decrease this value (or
   * choose it adaptively) when there are very slim cells present in
   * the mesh.
   */
  Real line_width;

  /**
   * Draws an element with Bezier curves
   */
  void plot_quadratic_elem(const Elem* elem);

  /**
   * Draws an element with straight lines
   */
  void plot_linear_elem(const Elem* elem);

private:
  /**
   * Given a quadratic edge Elem which lies in the x-y plane,
   * computes the Bezier coefficients.  These may be passed to
   * the Postscript routine "curveto".
   */
  void _compute_edge_bezier_coeffs(const Elem* elem);

  /**
   * Coefficients of the transformation from physical-space
   * edge coordinates to Bezier basis coefficients.  Transforms
   * x and y separately.
   */
  //DenseMatrix<float> _M;
  static const float _bezier_transform[3][3];

  /**
   * Vector containing 3 points corresponding to Bezier coefficients,
   * as computed by _compute_edge_bezier_coeffs.
   */
  std::vector<Point> _bezier_coeffs;

  /**
   * Amount to add to every (x,y) point to place it in Postscript coordinates.
   */
  Point _offset;

  /**
   * Amount by which to stretch each point to place it in Postscript coordinates.
   */
  Real _scale;

  /**
   * A point object used for temporary calculations
   */
  Point _current_point;

  /**
   * Drawing style-independent data for a single cell.  This can be
   * used as a temporary buffer for storing data which may be sent to
   * the output stream multiple times.
   */
  std::ostringstream _cell_string;

  /**
   * Output file stream which will be opened when the file name is known
   */
  std::ofstream _out;
};

} // namespace libMesh

#endif // LIBMESH_POSTSCRIPT_IO_H
