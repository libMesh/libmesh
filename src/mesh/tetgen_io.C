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
#include <fstream>

// Local includes
#include "libmesh/tetgen_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_tet10.h"

namespace libMesh
{

// ------------------------------------------------------------
// TetgenIO class members
void TetGenIO::read (const std::string & name)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  std::string name_node, name_ele, dummy;

  // tetgen only works in 3D
  MeshInput<MeshBase>::mesh().set_mesh_dimension(3);

#if LIBMESH_DIM < 3
  libmesh_error_msg("Cannot open dimension 3 mesh file when configured without 3D support.");
#endif

  // Check name for *.node or *.ele extension.
  // Set std::istream for node_stream and ele_stream.
  //
  if (name.rfind(".node") < name.size())
    {
      name_node            = name;
      dummy                = name;
      std::size_t position = dummy.rfind(".node");
      name_ele             = dummy.replace(position, 5, ".ele");
    }
  else if (name.rfind(".ele") < name.size())
    {
      name_ele = name;
      dummy    = name;
      std::size_t position = dummy.rfind(".ele");
      name_node    = dummy.replace(position, 4, ".node");
    }
  else
    libmesh_error_msg("ERROR: Unrecognized file name: " << name);



  // Set the streams from which to read in
  std::ifstream node_stream (name_node.c_str());
  std::ifstream ele_stream  (name_ele.c_str());

  if ( !node_stream.good() || !ele_stream.good() )
    libmesh_error_msg("Error while opening either "     \
                      << name_node                      \
                      << " or "                         \
                      << name_ele);

  libMesh::out<< "TetGenIO found the tetgen files to read " <<std::endl;

  // Skip the comment lines at the beginning
  this->skip_comment_lines (node_stream, '#');
  this->skip_comment_lines (ele_stream, '#');

  // Read the nodes and elements from the streams
  this->read_nodes_and_elem (node_stream, ele_stream);
  libMesh::out<< "TetGenIO read in nodes and elements " <<std::endl;
}



void TetGenIO::read_nodes_and_elem (std::istream & node_stream,
                                    std::istream & ele_stream)
{
  _num_nodes    = 0;
  _num_elements = 0;

  // Read all the datasets.
  this->node_in    (node_stream);
  this->element_in (ele_stream);

  // some more clean-up
  _assign_nodes.clear();
}



//----------------------------------------------------------------------
// Function to read in the node table.
void TetGenIO::node_in (std::istream & node_stream)
{
  // Check input buffer
  libmesh_assert (node_stream.good());

  // Get a reference to the mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  unsigned int dimension=0, nAttributes=0, BoundaryMarkers=0;

  node_stream >> _num_nodes       // Read the number of nodes from the stream
              >> dimension        // Read the dimension from the stream
              >> nAttributes      // Read the number of attributes from stream
              >> BoundaryMarkers; // Read if or not boundary markers are included in *.node (0 or 1)

  // Read the nodal coordinates from the node_stream (*.node file).
  unsigned int node_lab=0;
  Point xyz;
  Real dummy;

  // If present, make room for node attributes to be stored.
  this->node_attributes.resize(nAttributes);
  for (unsigned i=0; i<nAttributes; ++i)
    this->node_attributes[i].resize(_num_nodes);


  for (unsigned int i=0; i<_num_nodes; i++)
    {
      // Check input buffer
      libmesh_assert (node_stream.good());

      node_stream >> node_lab  // node number
                  >> xyz(0)    // x-coordinate value
                  >> xyz(1)    // y-coordinate value
                  >> xyz(2);   // z-coordinate value

      // Read and store attributes from the stream.
      for (unsigned int j=0; j<nAttributes; j++)
        node_stream >> node_attributes[j][i];

      // Read (and discard) boundary marker if BoundaryMarker=1.
      // TODO: should we store this somehow?
      if (BoundaryMarkers == 1)
        node_stream >> dummy;

      // Store the new position of the node under its label.
      //_assign_nodes.insert (std::make_pair(node_lab,i));
      _assign_nodes[node_lab] = i;

      // Add this point to the Mesh.
      mesh.add_point(xyz, i);
    }
}



//----------------------------------------------------------------------
// Function to read in the element table.
void TetGenIO::element_in (std::istream & ele_stream)
{
  // Check input buffer
  libmesh_assert (ele_stream.good());

  // Get a reference to the mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Read the elements from the ele_stream (*.ele file).
  unsigned int element_lab=0, n_nodes=0, nAttri=0;

  ele_stream >> _num_elements // Read the number of tetrahedrons from the stream.
             >> n_nodes       // Read the number of nodes per tetrahedron from the stream (defaults to 4).
             >> nAttri;       // Read the number of attributes from stream.

  // Vector that assigns element nodes to their correct position.
  // TetGen is normaly 0-based
  // (right now this is strictly not necessary since it is the identity map,
  //  but in the future TetGen could change their numbering scheme.)
  static const unsigned int assign_elm_nodes[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  // If present, make room for element attributes to be stored.
  this->element_attributes.resize(nAttri);
  for (unsigned i=0; i<nAttri; ++i)
    this->element_attributes[i].resize(_num_elements);

  for (dof_id_type i=0; i<_num_elements; i++)
    {
      libmesh_assert (ele_stream.good());

      // TetGen only supports Tet4 and Tet10 elements.
      Elem * elem;

      if (n_nodes==4)
        elem = new Tet4;

      else if (n_nodes==10)
        elem = new Tet10;

      else
        libmesh_error_msg("Elements with " << n_nodes << " nodes are not supported in the LibMesh tetgen module.");

      elem->set_id(i);

      mesh.add_elem (elem);

      libmesh_assert(elem);
      libmesh_assert_equal_to (elem->n_nodes(), n_nodes);

      // Read the element label
      ele_stream >> element_lab;

      // Read node labels
      for (dof_id_type j=0; j<n_nodes; j++)
        {
          dof_id_type node_label;
          ele_stream >> node_label;

          // Assign node to element
          elem->set_node(assign_elm_nodes[j]) =
            mesh.node_ptr(_assign_nodes[node_label]);
        }

      // Read and store attributes from the stream.
      for (unsigned int j=0; j<nAttri; j++)
        ele_stream >> this->element_attributes[j][i];
    }
}


/**
 * This method implements writing a mesh to a specified ".poly" file.
 * ".poly" files defines so called Piecewise Linear Complex (PLC).
 */
void TetGenIO::write (const std::string & fname)
{
  // libmesh_assert three dimensions (should be extended later)
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().mesh_dimension(), 3);

  if (!(fname.rfind(".poly") < fname.size()))
    libmesh_error_msg("ERROR: Unrecognized file name: " << fname);

  // Open the output file stream
  std::ofstream out_stream (fname.c_str());

  // Make sure it opened correctly
  if (!out_stream.good())
    libmesh_file_error(fname.c_str());

  // Get a reference to the mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // Begin interfacing with the .poly file
  {
    // header:
    out_stream << "# poly file output generated by libmesh\n"
               << mesh.n_nodes() << " 3 0 0\n";

    // write the nodes:
    for (dof_id_type v=0; v<mesh.n_nodes(); v++)
      out_stream << v << " "
                 << mesh.point(v)(0) << " "
                 << mesh.point(v)(1) << " "
                 << mesh.point(v)(2) << "\n";
  }

  {
    // write the connectivity:
    out_stream << "# Facets:\n"
               << mesh.n_elem() << " 0\n";

    //     const_active_elem_iterator       it (mesh.elements_begin());
    //     const const_active_elem_iterator end(mesh.elements_end());

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    for ( ; it != end; ++it)
      out_stream << "1\n3 " // no. of facet polygons
        //  << (*it)->n_nodes() << " "
                 << (*it)->node_id(0)   << " "
                 << (*it)->node_id(1)   << " "
                 << (*it)->node_id(2)   << "\n";
  }

  // end of the file
  out_stream << "0\n"; // no holes output!
  out_stream << "\n\n# end of file\n";
}

} // namespace libMesh
