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



#ifndef LIBMESH_GNUPLOT_IO_H
#define LIBMESH_GNUPLOT_IO_H

// Local includes
#include "libmesh/mesh_output.h"

// C++ includes
#include <cstddef>


namespace libMesh
{

// forward declaration
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;

/**
 * This class implements writing meshes using GNUplot, designed for
 * use only with 1D meshes.
 *
 * \author David Knezevic
 * \date 2005
 */
class GnuPlotIO : public MeshOutput<MeshBase>
{
public:

  /**
   * Define enumerations to set plotting properties on construction
   */
  enum PlottingProperties { GRID_ON    = 1,
                            PNG_OUTPUT = 2};

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * To set the properties, we input a bitwise OR of the
   * GnuPlotIO::PlottingProperties enumerations,
   * e.g. GnuPlotIO::GRID_ON | GnuPlotIO::PNG_OUTPUT
   */
  explicit
  GnuPlotIO (const MeshBase &,
             const std::string & = std::string("FE 1D Solution"),
             int properties=0);

  /**
   * Write the mesh to the specified file.
   */
  virtual void write(const std::string &) override;

  /**
   * Bring in base class functionality for name resolution and to
   * avoid warnings about hidden overloaded virtual functions.
   */
  using MeshOutput<MeshBase>::write_nodal_data;

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string &,
                                 const std::vector<Number> &,
                                 const std::vector<std::string> &) override;

  /**
   * Set title of plot
   */
  void set_title(const std::string & title) { _title = title; }

  /**
   * Turn grid on or off.
   */
  void use_grid(bool grid) { _grid = grid; }


  /**
   * Write output to a .png file using gnuplot
   */
  void set_png_output(bool png_output) { _png_output = png_output; }

  /**
   * GNUplot automatically adjusts the x and y-axes of 2D plots
   * to "zoom in" on the data.  You can set this string to force
   * GNUplot to maintain a fixed set of axes.
   * Example: axes_limits = "[0:1] [0:1]" would force x and y
   * to be plotted on the range 0<=x<=1 and 0<=y<=1 regardless
   * of where the data lie.
   */
  std::string axes_limits;

private:
  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write an ASCII file.
   */
  void write_solution (const std::string &,
                       const std::vector<Number> * = nullptr,
                       const std::vector<std::string> * = nullptr);

  std::string _title;

  bool _grid;
  bool _png_output;
};

} // namespace libMesh


#endif // LIBMESH_GNUPLOT_IO_H
