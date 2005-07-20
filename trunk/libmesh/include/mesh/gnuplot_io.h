// $Id: gnuplot_io.h,v 1.4 2005-07-20 20:20:05 knezed01 Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __gnuplot_io_h__
#define __gnuplot_io_h__

#include "mesh_output.h"


// forward declaration
class MeshBase;

/**
 * This class implements writing meshes using GNUplot, designed for use only
 * with 1D meshes.
 *
 * @author David Knezevic, 2005
 */

// ------------------------------------------------------------
// GnuPlotIO class definition
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
  GnuPlotIO (const MeshBase&, const std::string& = "FE 1D Solution", 
             int properties);

  /**
   * Write the mesh to the specified file.
   */
  virtual void write(const std::string&);

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string&,
				 const std::vector<Number>&,
				 const std::vector<std::string>&);

  /**
   * Set title of plot
   */
  void set_title(const std::string& title) { _title = title; }

  /**
   * Turn grid on or off.
   */
  void use_grid(bool grid) { _grid = grid; }


  /**
   * Write output to a .png file useing gnuplot
   */
  void set_png_output(bool png_output) { _png_output = png_output; }

 private:
  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write an ASCII file.
   */
  void write_solution (const std::string&,
		       const std::vector<Number>* = NULL,
		       const std::vector<std::string>* = NULL);

  std::string _title;

  bool _grid;
  bool _png_output;
};

    
#endif
