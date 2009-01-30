// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include <ctime>

// Local includes
#include "postscript_io.h"
#include "mesh_tools.h"
#include "elem.h"

PostscriptIO::PostscriptIO (const MeshBase& mesh)
  : MeshOutput<MeshBase> (mesh),
    shade_value(0.0)
{
  // This code is still undergoing some development.
  untested();
}



PostscriptIO::~PostscriptIO ()
{
}



void PostscriptIO::write (const std::string& fname)
{
  if (libMesh::processor_id() == 0)
    {
      // Get a constant reference to the mesh.
      const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

      // Only works in 2D
      libmesh_assert(mesh.mesh_dimension() == 2);
      
      // Create output file stream
      std::ofstream out(fname.c_str());

      // Make sure it opened correctly
      if (!out.good())
        libmesh_file_error(fname.c_str());

      // The mesh bounding box gives us info about what the
      // Postscript bounding box should be.
      MeshTools::BoundingBox bbox = MeshTools::bounding_box(mesh);

      // Add a little extra padding to the "true" bounding box so
      // that we can still see the boundary
      const Real percent_padding = 0.01;
      const Real dx=bbox.second(0)-bbox.first(0); libmesh_assert(dx > 0.0);
      const Real dy=bbox.second(1)-bbox.first(1); libmesh_assert(dy > 0.0);
      
      const Real x_min = bbox.first(0)  - percent_padding*dx;
      const Real y_min = bbox.first(1)  - percent_padding*dy;
      const Real x_max = bbox.second(0) + percent_padding*dx;
      const Real y_max = bbox.second(1) + percent_padding*dy;

      // Width of the output as given in postscript units.
      // This usually is given by the strange unit 1/72 inch. 
      // A width of 300 represents a size of roughly 10 cm.
      const Real width = 300;
      const Real scale = width / (x_max-x_min);
      const Point offset(x_min, y_min);
      
      // Header writing stuff stolen from Deal.II
      std::time_t  time1= std::time (0);
      std::tm     *time = std::localtime(&time1);
      out << "%!PS-Adobe-2.0 EPSF-1.2" << '\n'
	//<< "%!PS-Adobe-1.0" << '\n' // Lars' PS version
	  << "%%Filename: " << fname << '\n'
	  << "%%Title: LibMesh Output" << '\n'
	  << "%%Creator: LibMesh: A C++ finite element library" << '\n'
	  << "%%Creation Date: "
	  << time->tm_year+1900 << "/"
	  << time->tm_mon+1 << "/"
	  << time->tm_mday << " - "
	  << time->tm_hour << ":"
	  << std::setw(2) << time->tm_min << ":"
	  << std::setw(2) << time->tm_sec << '\n'
	  << "%%BoundingBox: "
	// lower left corner
	  << "0 0 "
	// upper right corner
	  << static_cast<unsigned int>( rint((x_max-x_min) * scale ))
	  << ' '
	  << static_cast<unsigned int>( rint((y_max-y_min) * scale ))
	  << '\n';

      // define some abbreviations to keep
      // the output small:
      // m=move turtle to
      // l=define a line
      // s=set rgb color
      // sg=set gray value
      // lx=close the line and plot the line
      // lf=close the line and fill the interior
      out << "/m {moveto} bind def"      << '\n'
	  << "/l {lineto} bind def"      << '\n'
	  << "/s {setrgbcolor} bind def" << '\n'
	  << "/sg {setgray} bind def"    << '\n'
	  << "/lx {lineto closepath stroke} bind def" << '\n'
	  << "/lf {lineto closepath fill} bind def"   << '\n';

      out << "%%EndProlog" << '\n';
      //	  << '\n';
      
      // set line width.  0.5 is a reasonable default for printed images
      const Real line_width = 0.5;
      out << line_width << " setlinewidth" << '\n';

      // Set line cap and join options
      out << "1 setlinecap" << '\n';
      out << "1 setlinejoin" << '\n';
  
      // allow only five digits for output (instead of the default
      // six); this should suffice even for fine grids, but reduces
      // the file size significantly
      out << std::setprecision (5);

      // By default, our "stream <<" function for points generates formatted output
      // for x, y, and z.  Here we just want to print x and y so we'll
      // use this temporary Point object.
      Point current_point;

      // Shading-independent data for drawing a single cell.  This can be
      // used to draw each cell twice, once for a filled cell and once for
      // black edges.
      std::ostringstream cell_string;
	
      // Loop over the active elements, draw lines for the edges.  We
      // draw even quadratic elements with straight sides, i.e. a straight
      // line sits between each pair of vertices.  Also we draw every edge
      // for an element regardless of the fact that it may overlap with
      // another.  This would probably be a useful optimization...
      MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_elements_end(); 
      for ( ; el != end_el; ++el)
	{
	  const Elem* elem = *el;

	  // Clear the string contents.  Yes, this really is how you do that...
	  cell_string.str("");
	  
	  // The general strategy is:
	  // 1.) Use m  := {moveto} to go to vertex 0.
	  // 2.) Use l  := {lineto} commands to draw lines to vertex 1, 2, ... N-1.
	  // 3a.) Use lx := {lineto closepath stroke} command at  vertex N to draw the last line.
	  // 3b.)     lf := {lineto closepath fill} command to shade the cell just drawn 
	  // All of our 2D elements' vertices are numbered in counterclockwise order,
	  // so we can just draw them in the same order.

	  // 1.)
	  current_point = (elem->point(0) - offset) * scale;
	  cell_string << current_point(0) << " " << current_point(1) << " "; // write x y 
	  cell_string << "m ";

	  // 2.)
	  const unsigned int nv=elem->n_vertices();
	  for (unsigned int v=1; v<nv-1; ++v)
	    {
	      current_point = (elem->point(v) - offset) * scale;
	      cell_string << current_point(0) << " " << current_point(1) << " "; // write x y 
	      cell_string << "l ";
	    }
	  
	  // 3.)
	  current_point = (elem->point(nv-1) - offset) * scale;
	  cell_string << current_point(0) << " " << current_point(1) << " "; // write x y 

	  // We draw the shaded (interior) parts first, if applicable.
	  if (shade_value > 0.0)
	    out << shade_value << " sg " << cell_string.str() << "lf\n";
	  
	  // Draw the black lines (I guess we will always do this)
	  out << "0 sg " << cell_string.str() << "lx\n";
	}
      
      // Issue the showpage command, and we're done.
      out << "showpage" << std::endl;
      
    } // end if (libMesh::processor_id() == 0)
}
