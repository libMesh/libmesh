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
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

// Local includes
#include "libmesh/postscript_io.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/elem.h"

namespace libMesh
{


// Transformation map between monomial (physical space) and Bezier bases.
const float PostscriptIO::_bezier_transform[3][3] =
  {
    {-1.f/6.f, 1.f/6.f, 1.},
    {-1.f/6.f, 0.5,     1.f/6.f},
    {0.,       1.,      0.}
  };


PostscriptIO::PostscriptIO (const MeshBase & mesh_in) :
  MeshOutput<MeshBase> (mesh_in),
  shade_value(0.0),
  line_width(0.5),
  //_M(3,3),
  _offset(0., 0.),
  _scale(1.0),
  _current_point(0., 0.)
{
  // This code is still undergoing some development.
  libmesh_experimental();

  // Entries of transformation matrix from physical to Bezier coords.
  // _M(0,0) = -1./6.;    _M(0,1) = 1./6.;    _M(0,2) = 1.;
  // _M(1,0) = -1./6.;    _M(1,1) = 0.5  ;    _M(1,2) = 1./6.;
  // _M(2,0) = 0.    ;    _M(2,1) = 1.   ;    _M(2,2) = 0.;

  // Make sure there is enough room to store Bezier coefficients.
  _bezier_coeffs.resize(3);
}



PostscriptIO::~PostscriptIO ()
{
}



void PostscriptIO::write (const std::string & fname)
{
  // We may need to gather a DistributedMesh to output it, making that
  // const qualifier in our constructor a dirty lie
  MeshSerializer serialize(const_cast<MeshBase &>(this->mesh()), !_is_parallel_format);

  if (this->mesh().processor_id() == 0)
    {
      // Get a constant reference to the mesh.
      const MeshBase & the_mesh = MeshOutput<MeshBase>::mesh();

      // Only works in 2D
      libmesh_assert_equal_to (the_mesh.mesh_dimension(), 2);

      // Create output file stream.
      // _out is now a private member of the class.
      _out.open(fname.c_str());

      // Make sure it opened correctly
      if (!_out.good())
        libmesh_file_error(fname.c_str());

      // The mesh bounding box gives us info about what the
      // Postscript bounding box should be.
      MeshTools::BoundingBox bbox = MeshTools::bounding_box(the_mesh);

      // Add a little extra padding to the "true" bounding box so
      // that we can still see the boundary
      const Real percent_padding = 0.01;
      const Real dx=bbox.second(0)-bbox.first(0); libmesh_assert_greater (dx, 0.0);
      const Real dy=bbox.second(1)-bbox.first(1); libmesh_assert_greater (dy, 0.0);

      const Real x_min = bbox.first(0)  - percent_padding*dx;
      const Real y_min = bbox.first(1)  - percent_padding*dy;
      const Real x_max = bbox.second(0) + percent_padding*dx;
      const Real y_max = bbox.second(1) + percent_padding*dy;

      // Width of the output as given in postscript units.
      // This usually is given by the strange unit 1/72 inch.
      // A width of 300 represents a size of roughly 10 cm.
      const Real width = 300;
      _scale = width / (x_max-x_min);
      _offset(0) = x_min;
      _offset(1) = y_min;

      // Header writing stuff stolen from Deal.II
      std::time_t  time1= std::time (0);
      std::tm     * time = std::localtime(&time1);
      _out << "%!PS-Adobe-2.0 EPSF-1.2" << '\n'
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
           << static_cast<unsigned int>( rint((x_max-x_min) * _scale ))
           << ' '
           << static_cast<unsigned int>( rint((y_max-y_min) * _scale ))
           << '\n';

      // define some abbreviations to keep
      // the output small:
      // m=move turtle to
      // l=define a line
      // s=set rgb color
      // sg=set gray value
      // lx=close the line and plot the line
      // lf=close the line and fill the interior
      _out << "/m {moveto} bind def"      << '\n'
           << "/l {lineto} bind def"      << '\n'
           << "/s {setrgbcolor} bind def" << '\n'
           << "/sg {setgray} bind def"    << '\n'
           << "/cs {curveto stroke} bind def" << '\n'
           << "/lx {lineto closepath stroke} bind def" << '\n'
           << "/lf {lineto closepath fill} bind def"   << '\n';

      _out << "%%EndProlog" << '\n';
      //  << '\n';

      // Set line width in the postscript file.
      _out << line_width << " setlinewidth" << '\n';

      // Set line cap and join options
      _out << "1 setlinecap" << '\n';
      _out << "1 setlinejoin" << '\n';

      // allow only five digits for output (instead of the default
      // six); this should suffice even for fine grids, but reduces
      // the file size significantly
      _out << std::setprecision (5);

      // Loop over the active elements, draw lines for the edges.  We
      // draw even quadratic elements with straight sides, i.e. a straight
      // line sits between each pair of vertices.  Also we draw every edge
      // for an element regardless of the fact that it may overlap with
      // another.  This would probably be a useful optimization...
      MeshBase::const_element_iterator       el     = the_mesh.active_elements_begin();
      const MeshBase::const_element_iterator end_el = the_mesh.active_elements_end();
      for ( ; el != end_el; ++el)
        {
          this->plot_linear_elem(*el);
          //this->plot_quadratic_elem(*el); // Experimental
        }

      // Issue the showpage command, and we're done.
      _out << "showpage" << std::endl;

    } // end if (this->mesh().processor_id() == 0)
}






void PostscriptIO::plot_linear_elem(const Elem * elem)
{
  // Clear the string contents.  Yes, this really is how you do that...
  _cell_string.str("");

  // The general strategy is:
  // 1.) Use m  := {moveto} to go to vertex 0.
  // 2.) Use l  := {lineto} commands to draw lines to vertex 1, 2, ... N-1.
  // 3a.) Use lx := {lineto closepath stroke} command at  vertex N to draw the last line.
  // 3b.)     lf := {lineto closepath fill} command to shade the cell just drawn
  // All of our 2D elements' vertices are numbered in counterclockwise order,
  // so we can just draw them in the same order.

  // 1.)
  _current_point = (elem->point(0) - _offset) * _scale;
  _cell_string << _current_point(0) << " " << _current_point(1) << " "; // write x y
  _cell_string << "m ";

  // 2.)
  const unsigned int nv=elem->n_vertices();
  for (unsigned int v=1; v<nv-1; ++v)
    {
      _current_point = (elem->point(v) - _offset) * _scale;
      _cell_string << _current_point(0) << " " << _current_point(1) << " "; // write x y
      _cell_string << "l ";
    }

  // 3.)
  _current_point = (elem->point(nv-1) - _offset) * _scale;
  _cell_string << _current_point(0) << " " << _current_point(1) << " "; // write x y

  // We draw the shaded (interior) parts first, if applicable.
  if (shade_value > 0.0)
    _out << shade_value << " sg " << _cell_string.str() << "lf\n";

  // Draw the black lines (I guess we will always do this)
  _out << "0 sg " << _cell_string.str() << "lx\n";
}




void PostscriptIO::plot_quadratic_elem(const Elem * elem)
{
  for (unsigned int ns=0; ns<elem->n_sides(); ++ns)
    {
      // Build the quadratic side
      UniquePtr<const Elem> side = elem->build_side_ptr(ns);

      // Be sure it's quadratic (Edge2).  Eventually we could
      // handle cubic elements as well...
      libmesh_assert_equal_to ( side->type(), EDGE3 );

      _out << "0 sg ";

      // Move to the first point on this side.
      _current_point = (side->point(0) - _offset) * _scale;
      _out << _current_point(0) << " " << _current_point(1) << " "; // write x y
      _out << "m ";

      // Compute _bezier_coeffs for this edge.  This fills up
      // the _bezier_coeffs vector.
      this->_compute_edge_bezier_coeffs(side.get());

      // Print curveto path to file
      for (std::size_t i=0; i<_bezier_coeffs.size(); ++i)
        _out << _bezier_coeffs[i](0) << " " << _bezier_coeffs[i](1) << " ";
      _out << " cs\n";
    }
}




void PostscriptIO::_compute_edge_bezier_coeffs(const Elem * elem)
{
  // I only know how to do this for an Edge3!
  libmesh_assert_equal_to (elem->type(), EDGE3);

  // Get x-coordinates into an array, transform them,
  // and repeat for y.
  float phys_coords[3] = {0., 0., 0.};
  float bez_coords[3]  = {0., 0., 0.};

  for (unsigned int i=0; i<2; ++i)
    {
      // Initialize vectors.  Physical coordinates are initialized
      // by their postscript-scaled values.
      for (unsigned int j=0; j<3; ++j)
        {
          phys_coords[j] = static_cast<float>
            ((elem->point(j)(i) - _offset(i)) * _scale);
          bez_coords[j] = 0.; // zero out result vector
        }

      // Multiply matrix times vector
      for (unsigned int j=0; j<3; ++j)
        for (unsigned int k=0; k<3; ++k)
          bez_coords[j] += _bezier_transform[j][k]*phys_coords[k];

      // Store result in _bezier_coeffs
      for (unsigned int j=0; j<3; ++j)
        _bezier_coeffs[j](i) = phys_coords[j];
    }
}

} // namespace libMesh
