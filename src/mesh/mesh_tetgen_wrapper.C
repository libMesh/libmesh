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

#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_TETGEN

// C++ includes
#include <iostream>

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_tetgen_wrapper.h"

namespace libMesh
{

TetGenWrapper::TetGenWrapper() :
  tetgen_output(new tetgenio)
{
  this->tetgen_data.mesh_dim                = 3;
  this->tetgen_data.numberofpointattributes = 0;
  this->tetgen_data.firstnumber             = 0;
}



TetGenWrapper::~TetGenWrapper()
{
}



void TetGenWrapper::set_node(unsigned i, REAL x, REAL y, REAL z)
{
  unsigned index = i*3;
  tetgen_data.pointlist[index++] = x;
  tetgen_data.pointlist[index++] = y;
  tetgen_data.pointlist[index++] = z;
}



void TetGenWrapper::set_hole(unsigned i, REAL x, REAL y, REAL z)
{
  unsigned index = i*3;
  tetgen_data.holelist[index++] = x;
  tetgen_data.holelist[index++] = y;
  tetgen_data.holelist[index++] = z;
}



void TetGenWrapper::set_numberofpoints(int i)
{
  // This is an int in tetgen, so use an int here even though it should be unsigned
  tetgen_data.numberofpoints = i;
}



void TetGenWrapper::get_output_node(unsigned i, REAL & x, REAL & y, REAL & z)
{
  // Bounds checking...
  if (i >= static_cast<unsigned>(tetgen_output->numberofpoints))
    libmesh_error_msg("Error, requested point "        \
                      << i                             \
                      << ", but there are only "       \
                      << tetgen_output->numberofpoints \
                      << " points available.");

  x = tetgen_output->pointlist[3*i];
  y = tetgen_output->pointlist[3*i+1];
  z = tetgen_output->pointlist[3*i+2];
}



int TetGenWrapper::get_numberoftetrahedra()
{
  return tetgen_output->numberoftetrahedra;
}



int TetGenWrapper::get_numberoftrifaces()
{
  return tetgen_output->numberoftrifaces;
}



int TetGenWrapper::get_numberofpoints()
{
  return tetgen_output->numberofpoints;
}



int TetGenWrapper::get_element_node(unsigned i, unsigned j)
{
  return tetgen_output->tetrahedronlist[i*4+j];
}



int TetGenWrapper::get_triface_node(unsigned i, unsigned j)
{
  return tetgen_output->trifacelist[i*3+j];
}



REAL TetGenWrapper::get_element_attribute(unsigned i)
{
  libmesh_assert(tetgen_output->numberoftetrahedronattributes>0);
  return tetgen_output->tetrahedronattributelist[tetgen_output->numberoftetrahedronattributes*i];
}



void TetGenWrapper::allocate_pointlist(int numofpoints)
{
  // This is stored as an int in tetgen, so we store it that way as well.
  this->set_numberofpoints(numofpoints);

  // Don't try to allocate an array of size zero, this is not portable...
  if (this->tetgen_data.numberofpoints > 0)
    {
      // Is there previously-allocated memory here?
      if (this->tetgen_data.pointlist != nullptr)
        libmesh_error_msg("Cannot allocate on top of previously allocated memory!");

      // We allocate memory here, the tetgenio destructor will delete it.
      this->tetgen_data.pointlist = new REAL[this->tetgen_data.numberofpoints * 3];
    }
}



void TetGenWrapper::set_switches(const std::string & s)
{
  // A temporary buffer for passing to the C API, it requires
  // a char *, not a const char *...
  char buffer[256];

  // Make sure char buffer has enough room
  if (s.size() >= sizeof(buffer)-1)
    libmesh_error_msg("Fixed size buffer of length "                  \
                      << sizeof(buffer)                               \
                      << " not large enough to hold TetGen switches.");

  // Copy the string, don't forget to terminate!
  buffer[ s.copy( buffer , sizeof( buffer ) - 1 ) ] = '\0' ;

  if (!tetgen_be.parse_commandline(buffer))
    libMesh::out << "TetGen replies: Wrong switches!" << std::endl;
}



void TetGenWrapper::run_tetgen()
{
  // Call tetrahedralize from the TetGen library.
  tetrahedralize(&tetgen_be, &tetgen_data, tetgen_output.get());
}



void TetGenWrapper::set_numberoffacets(int i)
{
  // This is stored as an int in TetGen
  this->tetgen_data.numberoffacets = i;
}



void TetGenWrapper::set_numberofholes(int i)
{
  // This is stored as an int in TetGen
  this->tetgen_data.numberofholes = i;
}



void TetGenWrapper::set_numberofregions(int i)
{
  // This is stored as an int in TetGen
  this->tetgen_data.numberofregions = i;
}



void TetGenWrapper::allocate_facetlist(int numoffacets, int numofholes)
{
  // These are both stored as ints in TetGen
  this->set_numberoffacets(numoffacets);
  this->set_numberofholes(numofholes);

  // Don't try to allocate an array of size zero, this is not portable...
  if (this->tetgen_data.numberoffacets > 0)
    {
      // Is there previously-allocated memory here?
      if (this->tetgen_data.facetlist != nullptr)
        libmesh_error_msg("Cannot allocate on top of previously allocated memory!");

      // We allocate memory here, the tetgenio destructor cleans it up.
      this->tetgen_data.facetlist = new tetgenio::facet[this->tetgen_data.numberoffacets];

      for (int i=0; i<numoffacets; i++)
        this->tetgen_data.init(&(this->tetgen_data.facetlist[i]));
    }


  // Don't try to allocate an array of size zero, this is not portable...
  if (this->tetgen_data.numberofholes > 0)
    {
      // Is there previously-allocated memory here?
      if (this->tetgen_data.holelist != nullptr)
        libmesh_error_msg("Cannot allocate on top of previously allocated memory!");

      this->tetgen_data.holelist = new REAL[this->tetgen_data.numberofholes * 3];
    }
}



void TetGenWrapper::allocate_regionlist(int numofregions)
{
  this->set_numberofregions(numofregions);

  // Don't try to allocate an array of size zero, this is not portable...
  if (this->tetgen_data.numberofregions > 0)
    {
      // Is there previously-allocated memory here?
      if (this->tetgen_data.regionlist != nullptr)
        libmesh_error_msg("Cannot allocate on top of previously allocated memory!");

      // We allocate memory here, the tetgenio destructor cleans it up.
      this->tetgen_data.regionlist = new REAL[this->tetgen_data.numberofregions * 5];
    }
}



void TetGenWrapper::set_facet_numberofpolygons(unsigned i, int num)
{
  // numberofpolygons is stored as an int in TetGen
  this->tetgen_data.facetlist[i].numberofpolygons = num;
}



void TetGenWrapper::set_facet_numberofholes(unsigned i, int num)
{
  // numberofholes is stored as an int in TetGen
  this->tetgen_data.facetlist[i].numberofholes = num;
}




void TetGenWrapper::allocate_facet_polygonlist(unsigned i, int numofpolygons)
{
  this->set_facet_numberofpolygons(i, numofpolygons);
  this->set_facet_numberofholes(i, 0);

  // Don't try to create an array of size zero, this isn't portable
  if (numofpolygons > 0)
    {
      // Is there previously-allocated memory here?
      if (this->tetgen_data.facetlist[i].polygonlist != nullptr)
        libmesh_error_msg("Cannot allocate on top of previously allocated memory!");

      // We allocate memory here, the tetgenio destructor cleans it up.
      this->tetgen_data.facetlist[i].polygonlist = new tetgenio::polygon[numofpolygons];

      for (int j=0; j<this->tetgen_data.facetlist[i].numberofpolygons; j++)
        this->tetgen_data.init(&(this->tetgen_data.facetlist[i].polygonlist[j]));
    }
}



void TetGenWrapper::set_polygon_numberofvertices(unsigned i, unsigned j, int num)
{
  // numberofvertices is stored as an int in TetGen
  this->tetgen_data.facetlist[i].polygonlist[j].numberofvertices = num;
}



void TetGenWrapper::allocate_polygon_vertexlist(unsigned i, unsigned j, int numofvertices)
{
  this->set_polygon_numberofvertices(i, j, numofvertices);

  // Don't try to create an array of size zero, this isn't portable
  if (numofvertices > 0)
    {
      // Is there previously-allocated memory here?
      if (this->tetgen_data.facetlist[i].polygonlist[j].vertexlist != nullptr)
        libmesh_error_msg("Cannot allocate on top of previously allocated memory!");

      // We allocate memory here, the tetgenio destructor cleans it up.
      this->tetgen_data.facetlist[i].polygonlist[j].vertexlist = new int[numofvertices];
    }
}




void TetGenWrapper::set_vertex(unsigned i, unsigned j, unsigned k, int nodeindex)
{
  // vertexlist entries are stored as ints in TetGen
  this->tetgen_data.facetlist[i].polygonlist[j].vertexlist[k] = nodeindex;
}



void TetGenWrapper::set_region(unsigned i, REAL x, REAL y, REAL z,
                               REAL attribute, REAL vol_constraint)
{
  unsigned index = i*5;
  tetgen_data.regionlist[index++] = x;
  tetgen_data.regionlist[index++] = y;
  tetgen_data.regionlist[index++] = z;
  tetgen_data.regionlist[index++] = attribute;
  tetgen_data.regionlist[index++] = vol_constraint;
}

} // namespace libMesh


#endif // LIBMESH_HAVE_TETGEN
