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


// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/unv_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/face_quad8.h"
#include "libmesh/face_quad9.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex20.h"
#include "libmesh/cell_tet10.h"
#include "libmesh/cell_prism6.h"
#include LIBMESH_INCLUDE_UNORDERED_MAP

// C++ includes
#include <iomanip>
#include <algorithm> // for std::sort
#include <fstream>
#include <ctype.h> // isspace
#include <sstream> // std::istringstream

#ifdef LIBMESH_HAVE_GZSTREAM
# include "gzstream.h" // For reading/writing compressed streams
#endif



namespace libMesh
{

//-----------------------------------------------------------------------------
// UNVIO class static members
const std::string UNVIO::_nodes_dataset_label    = "2411";
const std::string UNVIO::_elements_dataset_label = "2412";
const std::string UNVIO::_groups_dataset_label   = "2467";



// ------------------------------------------------------------
// UNVIO class members

UNVIO::UNVIO (MeshBase & mesh) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _verbose (false)
{
}



UNVIO::UNVIO (const MeshBase & mesh) :
  MeshOutput<MeshBase> (mesh),
  _verbose (false)
{
}



UNVIO::~UNVIO ()
{
}



bool & UNVIO::verbose ()
{
  return _verbose;
}



void UNVIO::read (const std::string & file_name)
{
  if (file_name.rfind(".gz") < file_name.size())
    {
#ifdef LIBMESH_HAVE_GZSTREAM

      igzstream in_stream (file_name.c_str());
      this->read_implementation (in_stream);

#else

      libmesh_error_msg("ERROR:  You must have the zlib.h header files and libraries to read and write compressed streams.");

#endif
      return;
    }

  else
    {
      std::ifstream in_stream (file_name.c_str());
      this->read_implementation (in_stream);
      return;
    }
}


void UNVIO::read_implementation (std::istream & in_stream)
{
  // Keep track of what kinds of elements this file contains
  elems_of_dimension.clear();
  elems_of_dimension.resize(4, false);

  {
    if ( !in_stream.good() )
      libmesh_error_msg("ERROR: Input file not good.");

    // Flags to be set when certain sections are encountered
    bool
      found_node  = false,
      found_elem  = false,
      found_group = false;

    // strings for reading the file line by line
    std::string
      old_line,
      current_line;

    while (true)
      {
        // Save off the old_line.  This will provide extra reliability
        // for detecting the beginnings of the different sections.
        old_line = current_line;

        // Try to read something.  This may set EOF!
        std::getline(in_stream, current_line);

        // If the stream is still "valid", parse the line
        if (in_stream)
          {
            // UNV files always have some amount of leading
            // whitespace, let's not rely on exactly how much...  This
            // command deletes it.
            current_line.erase(std::remove_if(current_line.begin(), current_line.end(), isspace),
                               current_line.end());

            // Parse the nodes section
            if (current_line == _nodes_dataset_label &&
                old_line == "-1")
              {
                found_node = true;
                this->nodes_in(in_stream);
              }

            // Parse the elements section
            else if (current_line == _elements_dataset_label &&
                     old_line == "-1")
              {
                // The current implementation requires the nodes to
                // have been read before reaching the elements
                // section.
                if (!found_node)
                  libmesh_error_msg("ERROR: The Nodes section must come before the Elements section of the UNV file!");

                found_elem = true;
                this->elements_in(in_stream);
              }

            // Parse the groups section
            else if (current_line == _groups_dataset_label &&
                     old_line == "-1")
              {
                // The current implementation requires the nodes and
                // elements to have already been read before reaching
                // the groups section.
                if (!found_node || !found_elem)
                  libmesh_error_msg("ERROR: The Nodes and Elements sections must come before the Groups section of the UNV file!");

                found_group = true;
                this->groups_in(in_stream);
              }

            // We can stop reading once we've found the nodes, elements,
            // and group sections.
            if (found_node && found_elem && found_group)
              break;

            // If we made it here, we're not done yet, so keep reading
            continue;
          }

        // if (!in_stream) check to see if EOF was set.  If so, break out of while loop.
        if (in_stream.eof())
          break;

        // If !in_stream and !in_stream.eof(), stream is in a bad state!
        libmesh_error_msg("Stream is bad! Perhaps the file does not exist?");
      } // end while (true)

    // By now we better have found the datasets for nodes and elements,
    // otherwise the unv file is bad!
    if (!found_node)
      libmesh_error_msg("ERROR: Could not find nodes!");

    if (!found_elem)
      libmesh_error_msg("ERROR: Could not find elements!");
  }



  {
    // Set the mesh dimension to the largest element dimension encountered
    MeshInput<MeshBase>::mesh().set_mesh_dimension(max_elem_dimension_seen());

#if LIBMESH_DIM < 3
    if (MeshInput<MeshBase>::mesh().mesh_dimension() > LIBMESH_DIM)
      libmesh_error_msg("Cannot open dimension "                        \
                        << MeshInput<MeshBase>::mesh().mesh_dimension() \
                        << " mesh file when configured without "        \
                        << MeshInput<MeshBase>::mesh().mesh_dimension() \
                        << "D support." );
#endif

    // Delete any lower-dimensional elements that might have been
    // added to the mesh stricly for setting BCs.
    {
      // Grab reference to the Mesh
      MeshBase & mesh = MeshInput<MeshBase>::mesh();

      unsigned char max_dim = this->max_elem_dimension_seen();

      MeshBase::const_element_iterator       el     = mesh.elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.elements_end();

      for (; el != end_el; ++el)
        {
          Elem * elem = *el;

          if (elem->dim() < max_dim)
            mesh.delete_elem(elem);
        }
    }

    if (this->verbose())
      libMesh::out << "  Finished." << std::endl << std::endl;
  }
}





void UNVIO::write (const std::string & file_name)
{
  if (file_name.rfind(".gz") < file_name.size())
    {
#ifdef LIBMESH_HAVE_GZSTREAM

      ogzstream out_stream(file_name.c_str());
      this->write_implementation (out_stream);

#else

      libmesh_error_msg("ERROR:  You must have the zlib.h header files and libraries to read and write compressed streams.");

#endif

      return;
    }

  else
    {
      std::ofstream out_stream (file_name.c_str());
      this->write_implementation (out_stream);
      return;
    }
}




void UNVIO::write_implementation (std::ostream & out_file)
{
  if ( !out_file.good() )
    libmesh_error_msg("ERROR: Output file not good.");

  // write the nodes, then the elements
  this->nodes_out    (out_file);
  this->elements_out (out_file);
}




void UNVIO::nodes_in (std::istream & in_file)
{
  LOG_SCOPE("nodes_in()","UNVIO");

  if (this->verbose())
    libMesh::out << "  Reading nodes" << std::endl;

  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // node label, we use an int here so we can read in a -1
  int node_label;

  // always 3 coordinates in the UNV file, no matter
  // what LIBMESH_DIM is.
  Point xyz;

  // We'll just read the floating point values as strings even when
  // there are no "D" characters in the file.  This will make parsing
  // the file a bit slower but the logic to do so will all be
  // simplified and in one place...
  std::string line;
  std::istringstream coords_stream;

  // Continue reading nodes until there are none left
  unsigned ctr = 0;
  while (true)
    {
      // Read the node label
      in_file >> node_label;

      // Break out of the while loop when we hit -1
      if (node_label == -1)
        break;

      // Read and discard the the rest of the node data on this line
      // which we do not currently use:
      // .) exp_coord_sys_num
      // .) disp_coord_sys_num
      // .) color
      std::getline(in_file, line);

      // read the floating-point data
      std::getline(in_file, line);

      // Replace any "D" characters used for exponents with "E"
      size_t last_pos = 0;
      while (true)
        {
          last_pos = line.find("D", last_pos);

          if (last_pos != std::string::npos)
            line.replace(last_pos, 1, "E");
          else
            break;
        }

      // Construct a stream object from the current line and then
      // stream in the xyz values.
      coords_stream.str(line);
      coords_stream.clear(); // clear iostate bits!  Important!
      coords_stream >> xyz(0) >> xyz(1) >> xyz(2);

      // Add node to the Mesh, bump the counter.
      Node * added_node = mesh.add_point(xyz, ctr++);

      // Maintain the mapping between UNV node ids and libmesh Node
      // pointers.
      _unv_node_id_to_libmesh_node_ptr[node_label] = added_node;
    }
}



unsigned char UNVIO::max_elem_dimension_seen ()
{
  unsigned char max_dim = 0;

  unsigned char elem_dimensions_size = cast_int<unsigned char>
    (elems_of_dimension.size());
  // The elems_of_dimension array is 1-based in the UNV reader
  for (unsigned char i=1; i<elem_dimensions_size; ++i)
    if (elems_of_dimension[i])
      max_dim = i;

  return max_dim;
}



bool UNVIO::need_D_to_e (std::string & number)
{
  // find "D" in string, start looking at 6th element since dont expect a "D" earlier.
  std::string::size_type position = number.find("D", 6);

  if (position != std::string::npos)
    {
      // replace "D" in string
      number.replace(position, 1, "e");
      return true;
    }
  else
    return false;
}



void UNVIO::groups_in (std::istream & in_file)
{
  // Grab reference to the Mesh, so we can add boundary info data to it
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Record the max and min element dimension seen while reading the file.
  unsigned char max_dim = this->max_elem_dimension_seen();

  // Container which stores lower-dimensional elements (based on
  // Elem::key()) for later assignment of boundary conditions.  We
  // could use e.g. an unordered_set with Elems as keys for this, but
  // this turns out to be much slower on the search side, since we
  // have to build an entire side in order to search, rather than just
  // calling elem->key(side) to compute a key.
  typedef LIBMESH_BEST_UNORDERED_MULTIMAP<dof_id_type, Elem *> map_type;
  map_type provide_bcs;

  // Read groups until there aren't any more to read...
  while (true)
    {
      // If we read a -1, it means there is nothing else to read in this section.
      int group_number;
      in_file >> group_number;

      if (group_number == -1)
        break;

      // The first record consists of 8 entries:
      // Field 1 -- group number (that we just read)
      // Field 2 -- active constraint set no. for group
      // Field 3 -- active restraint set no. for group
      // Field 4 -- active load set no. for group
      // Field 5 -- active dof set no. for group
      // Field 6 -- active temperature set no. for group
      // Field 7 -- active contact set no. for group
      // Field 8 -- number of entities in group
      // Only the first and last of these are relevant to us...
      unsigned num_entities;
      std::string group_name;
      {
        unsigned dummy;
        in_file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
                >> num_entities;

        // The second record has 1 field, the group name
        in_file >> group_name;
      }

      // The dimension of the elements in the group will determine
      // whether this is a sideset group or a subdomain group.
      bool
        is_subdomain_group = false,
        is_sideset_group = false;

      // Each subsequent line has 8 entries, there are two entity type
      // codes and two tags per line that we need to read.  In all my
      // examples, the entity type code was always "8", which stands for
      // "finite elements" but it's possible that we could eventually
      // support other entity type codes...
      // Field 1 -- entity type code
      // Field 2 -- entity tag
      // Field 3 -- entity node leaf id.
      // Field 4 -- entity component/ ham id.
      // Field 5 -- entity type code
      // Field 6 -- entity tag
      // Field 7 -- entity node leaf id.
      // Field 8 -- entity component/ ham id.
      {
        unsigned entity_type_code, entity_tag, dummy;
        for (unsigned entity=0; entity<num_entities; ++entity)
          {
            in_file >> entity_type_code >> entity_tag >> dummy >> dummy;

            if (entity_type_code != 8)
              libMesh::err << "Warning, unrecognized entity type code = "
                           << entity_type_code
                           << " in group "
                           << group_name
                           << std::endl;

            // Try to find this ID in the map from UNV element ids to libmesh ids.
            std::map<unsigned, unsigned>::iterator it =
              _unv_elem_id_to_libmesh_elem_id.find(entity_tag);

            if (it != _unv_elem_id_to_libmesh_elem_id.end())
              {
                unsigned libmesh_elem_id = it->second;

                // Attempt to get a pointer to the elem listed in the group
                Elem * group_elem = mesh.elem_ptr(libmesh_elem_id);

                // dim < max_dim means the Elem defines a boundary condition
                if (group_elem->dim() < max_dim)
                  {
                    is_sideset_group = true;

                    // We can only handle elements that are *exactly*
                    // one dimension lower than the max element
                    // dimension.  Not sure if "edge" BCs in 3D
                    // actually make sense/are required...
                    if (group_elem->dim()+1 != max_dim)
                      libmesh_error_msg("ERROR: Expected boundary element of dimension " \
                                        << max_dim-1 << " but got " << group_elem->dim());

                    // Set the current group number as the lower-dimensional element's subdomain ID.
                    // We will use this later to set a boundary ID.
                    group_elem->subdomain_id() =
                      cast_int<subdomain_id_type>(group_number);

                    // Store the lower-dimensional element in the provide_bcs container.
                    provide_bcs.insert(std::make_pair(group_elem->key(), group_elem));
                  }

                // dim == max_dim means this group defines a subdomain ID
                else if (group_elem->dim() == max_dim)
                  {
                    is_subdomain_group = true;
                    group_elem->subdomain_id() =
                      cast_int<subdomain_id_type>(group_number);
                  }

                else
                  libmesh_error_msg("ERROR: Found an elem with dim=" \
                                    << group_elem->dim() << " > " << "max_dim=" << +max_dim);
              }
            else
              libMesh::err << "WARNING: UNV Element " << entity_tag << " not found while parsing groups." << std::endl;
          } // end for (entity)
      } // end scope

      // Associate this group_number with the group_name in the BoundaryInfo object.
      if (is_sideset_group)
        mesh.get_boundary_info().sideset_name
          (cast_int<boundary_id_type>(group_number)) = group_name;

      if (is_subdomain_group)
        mesh.subdomain_name
          (cast_int<subdomain_id_type>(group_number)) = group_name;

    } // end while (true)

  // Loop over elements and try to assign boundary information
  {
    MeshBase::element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::element_iterator end = mesh.active_elements_end();
    for ( ; it != end; ++it)
      {
        Elem * elem = *it;

        if (elem->dim() == max_dim)
          {
            // This is a max-dimension element that may require BCs.
            // For each of its sides, including internal sides, we'll
            // see if a lower-dimensional element provides boundary
            // information for it.  Note that we have not yet called
            // find_neighbors(), so we can't use elem->neighbor(sn) in
            // this algorithm...
            for (unsigned short sn=0; sn<elem->n_sides(); sn++)
              {
                // Look for this key in the provide_bcs map
                std::pair<map_type::const_iterator,
                          map_type::const_iterator>
                  range = provide_bcs.equal_range (elem->key(sn));

                // Add boundary information for each side in the range.
                for (map_type::const_iterator iter = range.first;
                     iter != range.second; ++iter)
                  {
                    // Build a side to confirm the hash mapped to the correct side.
                    UniquePtr<Elem> side (elem->build_side_ptr(sn));

                    // Get a pointer to the lower-dimensional element
                    Elem * lower_dim_elem = iter->second;

                    // This was a hash, so it might not be perfect.  Let's verify...
                    if (*lower_dim_elem == *side)
                      mesh.get_boundary_info().add_side(elem, sn, lower_dim_elem->subdomain_id());
                  }
              }
          }
      }
  }
}



void UNVIO::elements_in (std::istream & in_file)
{
  LOG_SCOPE("elements_in()","UNVIO");

  if (this->verbose())
    libMesh::out << "  Reading elements" << std::endl;

  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // node label, we use an int here so we can read in a -1
  int element_label;

  unsigned
    n_nodes,           // number of nodes on element
    fe_descriptor_id,  // FE descriptor id
    phys_prop_tab_num, // physical property table number (not supported yet)
    mat_prop_tab_num,  // material property table number (not supported yet)
    color;             // color (not supported yet)

  // vector that temporarily holds the node labels defining element
  std::vector<unsigned int> node_labels (21);

  // vector that assigns element nodes to their correct position
  // for example:
  // 44:plane stress      | QUAD4
  // linear quadrilateral |
  // position in UNV-file | position in libmesh
  // assign_elem_node[1]   = 0
  // assign_elem_node[2]   = 3
  // assign_elem_node[3]   = 2
  // assign_elem_node[4]   = 1
  //
  // UNV is 1-based, we leave the 0th element of the vectors unused in order
  // to prevent confusion, this way we can store elements with up to 20 nodes
  unsigned int assign_elem_nodes[21];

  // Read elements until there are none left
  unsigned ctr = 0;
  while (true)
    {
      // read element label, break out when we read -1
      in_file >> element_label;

      if (element_label == -1)
        break;

      in_file >> fe_descriptor_id        // read FE descriptor id
              >> phys_prop_tab_num       // (not supported yet)
              >> mat_prop_tab_num        // (not supported yet)
              >> color                   // (not supported yet)
              >> n_nodes;                // read number of nodes on element

      // For "beam" type elements, the next three numbers are:
      // .) beam orientation node number
      // .) beam fore-end cross section number
      // .) beam aft-end cross section number
      // which we discard in libmesh.  The "beam" type elements:
      // 11  Rod
      // 21  Linear beam
      // 22  Tapered beam
      // 23  Curved beam
      // 24  Parabolic beam
      // all have fe_descriptor_id < 25.
      // http://www.sdrl.uc.edu/universal-file-formats-for-modal-analysis-testing-1/file-format-storehouse/unv_2412.htm
      if (fe_descriptor_id < 25)
        {
          unsigned dummy;
          in_file >> dummy >> dummy >> dummy;
        }

      // read node labels (1-based)
      for (unsigned int j=1; j<=n_nodes; j++)
        in_file >> node_labels[j];

      // element pointer, to be allocated
      Elem * elem = libmesh_nullptr;

      switch (fe_descriptor_id)
        {
        case 11: // Rod
          {
            elem = new Edge2;

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=1;
            break;
          }

        case 41: // Plane Stress Linear Triangle
        case 91: // Thin Shell   Linear Triangle
          {
            elem = new Tri3;  // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=2;
            assign_elem_nodes[3]=1;
            break;
          }

        case 42: // Plane Stress Quadratic Triangle
        case 92: // Thin Shell   Quadratic Triangle
          {
            elem = new Tri6;  // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=5;
            assign_elem_nodes[3]=2;
            assign_elem_nodes[4]=4;
            assign_elem_nodes[5]=1;
            assign_elem_nodes[6]=3;
            break;
          }

        case 43: // Plane Stress Cubic Triangle
          libmesh_error_msg("ERROR: UNV-element type 43: Plane Stress Cubic Triangle not supported.");

        case 44: // Plane Stress Linear Quadrilateral
        case 94: // Thin Shell   Linear Quadrilateral
          {
            elem = new Quad4; // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=3;
            assign_elem_nodes[3]=2;
            assign_elem_nodes[4]=1;
            break;
          }

        case 45: // Plane Stress Quadratic Quadrilateral
        case 95: // Thin Shell   Quadratic Quadrilateral
          {
            elem = new Quad8; // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=7;
            assign_elem_nodes[3]=3;
            assign_elem_nodes[4]=6;
            assign_elem_nodes[5]=2;
            assign_elem_nodes[6]=5;
            assign_elem_nodes[7]=1;
            assign_elem_nodes[8]=4;
            break;
          }

        case 300: // Thin Shell   Quadratic Quadrilateral (nine nodes)
          {
            elem = new Quad9; // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=7;
            assign_elem_nodes[3]=3;
            assign_elem_nodes[4]=6;
            assign_elem_nodes[5]=2;
            assign_elem_nodes[6]=5;
            assign_elem_nodes[7]=1;
            assign_elem_nodes[8]=4;
            assign_elem_nodes[9]=8;
            break;
          }

        case 46: // Plane Stress Cubic Quadrilateral
          libmesh_error_msg("ERROR: UNV-element type 46: Plane Stress Cubic Quadrilateral not supported.");

        case 111: // Solid Linear Tetrahedron
          {
            elem = new Tet4;  // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=1;
            assign_elem_nodes[3]=2;
            assign_elem_nodes[4]=3;
            break;
          }

        case 112: // Solid Linear Prism
          {
            elem = new Prism6;  // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=1;
            assign_elem_nodes[3]=2;
            assign_elem_nodes[4]=3;
            assign_elem_nodes[5]=4;
            assign_elem_nodes[6]=5;
            break;
          }

        case 115: // Solid Linear Brick
          {
            elem = new Hex8;  // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=4;
            assign_elem_nodes[3]=5;
            assign_elem_nodes[4]=1;
            assign_elem_nodes[5]=3;
            assign_elem_nodes[6]=7;
            assign_elem_nodes[7]=6;
            assign_elem_nodes[8]=2;
            break;
          }

        case 116: // Solid Quadratic Brick
          {
            elem = new Hex20; // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=12;
            assign_elem_nodes[3]=4;
            assign_elem_nodes[4]=16;
            assign_elem_nodes[5]=5;
            assign_elem_nodes[6]=13;
            assign_elem_nodes[7]=1;
            assign_elem_nodes[8]=8;

            assign_elem_nodes[9]=11;
            assign_elem_nodes[10]=19;
            assign_elem_nodes[11]=17;
            assign_elem_nodes[12]=9;

            assign_elem_nodes[13]=3;
            assign_elem_nodes[14]=15;
            assign_elem_nodes[15]=7;
            assign_elem_nodes[16]=18;
            assign_elem_nodes[17]=6;
            assign_elem_nodes[18]=14;
            assign_elem_nodes[19]=2;
            assign_elem_nodes[20]=10;
            break;
          }

        case 117: // Solid Cubic Brick
          libmesh_error_msg("Error: UNV-element type 117: Solid Cubic Brick not supported.");

        case 118: // Solid Quadratic Tetrahedron
          {
            elem = new Tet10; // create new element

            assign_elem_nodes[1]=0;
            assign_elem_nodes[2]=4;
            assign_elem_nodes[3]=1;
            assign_elem_nodes[4]=5;
            assign_elem_nodes[5]=2;
            assign_elem_nodes[6]=6;
            assign_elem_nodes[7]=7;
            assign_elem_nodes[8]=8;
            assign_elem_nodes[9]=9;
            assign_elem_nodes[10]=3;
            break;
          }

        default: // Unrecognized element type
          libmesh_error_msg("ERROR: UNV-element type " << fe_descriptor_id << " not supported.");
        }

      // nodes are being stored in element
      for (dof_id_type j=1; j<=n_nodes; j++)
        {
          // Map the UNV node ID to the libmesh node ID
          std::map<dof_id_type, Node *>::iterator it =
            _unv_node_id_to_libmesh_node_ptr.find(node_labels[j]);

          if (it != _unv_node_id_to_libmesh_node_ptr.end())
            elem->set_node(assign_elem_nodes[j]) = it->second;
          else
            libmesh_error_msg("ERROR: UNV node " << node_labels[j] << " not found!");
        }

      elems_of_dimension[elem->dim()] = true;

      // Set the element's ID
      elem->set_id(ctr);

      // Maintain a map from the libmesh (0-based) numbering to the
      // UNV numbering.
      _unv_elem_id_to_libmesh_elem_id[element_label] = ctr;

      // Add the element to the Mesh
      mesh.add_elem(elem);

      // Increment the counter for the next iteration
      ctr++;
    } // end while(true)
}






void UNVIO::nodes_out (std::ostream & out_file)
{
  if (this->verbose())
    libMesh::out << "  Writing " << MeshOutput<MeshBase>::mesh().n_nodes() << " nodes" << std::endl;

  // Write beginning of dataset
  out_file << "    -1\n"
           << "  "
           << _nodes_dataset_label
           << '\n';

  unsigned int
    exp_coord_sys_dummy  = 0, // export coordinate sys. (not supported yet)
    disp_coord_sys_dummy = 0, // displacement coordinate sys. (not supp. yet)
    color_dummy          = 0; // color(not supported yet)

  // A reference to the parent class's mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // Use scientific notation with captial E and 16 digits for printing out the coordinates
  out_file << std::scientific << std::setprecision(16) << std::uppercase;

  MeshBase::const_node_iterator       nd  = mesh.nodes_begin();
  const MeshBase::const_node_iterator end = mesh.nodes_end();

  for (; nd != end; ++nd)
    {
      const Node * current_node = *nd;

      dof_id_type node_id = current_node->id();

      out_file << std::setw(10) << node_id
               << std::setw(10) << exp_coord_sys_dummy
               << std::setw(10) << disp_coord_sys_dummy
               << std::setw(10) << color_dummy
               << '\n';

      // The coordinates - always write out three coords regardless of LIBMESH_DIM
      Real x = (*current_node)(0);

#if LIBMESH_DIM > 1
      Real y = (*current_node)(1);
#else
      Real y = 0.;
#endif

#if LIBMESH_DIM > 2
      Real z = (*current_node)(2);
#else
      Real z = 0.;
#endif

      out_file << std::setw(25) << x
               << std::setw(25) << y
               << std::setw(25) << z
               << '\n';
    }

  // Write end of dataset
  out_file << "    -1\n";
}






void UNVIO::elements_out(std::ostream & out_file)
{
  if (this->verbose())
    libMesh::out << "  Writing elements" << std::endl;

  // Write beginning of dataset
  out_file << "    -1\n"
           << "  "
           << _elements_dataset_label
           << "\n";

  unsigned
    fe_descriptor_id = 0,    // FE descriptor id
    phys_prop_tab_dummy = 2, // physical property (not supported yet)
    mat_prop_tab_dummy = 1,  // material property (not supported yet)
    color_dummy = 7;         // color (not supported yet)

  // vector that assigns element nodes to their correct position
  // currently only elements with up to 20 nodes
  //
  // Example:
  // QUAD4               | 44:plane stress
  //                     | linear quad
  // position in libMesh | UNV numbering
  // (note: 0-based)     | (note: 1-based)
  //
  // assign_elem_node[0]  = 0
  // assign_elem_node[1]  = 3
  // assign_elem_node[2]  = 2
  // assign_elem_node[3]  = 1
  unsigned int assign_elem_nodes[20];

  unsigned int n_elem_written=0;

  // A reference to the parent class's mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  MeshBase::const_element_iterator it  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  for (; it != end; ++it)
    {
      const Elem * elem = *it;

      elem->n_nodes();

      switch (elem->type())
        {

        case TRI3:
          {
            fe_descriptor_id = 41; // Plane Stress Linear Triangle
            assign_elem_nodes[0] = 0;
            assign_elem_nodes[1] = 2;
            assign_elem_nodes[2] = 1;
            break;
          }

        case TRI6:
          {
            fe_descriptor_id = 42; // Plane Stress Quadratic Triangle
            assign_elem_nodes[0] = 0;
            assign_elem_nodes[1] = 5;
            assign_elem_nodes[2] = 2;
            assign_elem_nodes[3] = 4;
            assign_elem_nodes[4] = 1;
            assign_elem_nodes[5] = 3;
            break;
          }

        case QUAD4:
          {
            fe_descriptor_id = 44; // Plane Stress Linear Quadrilateral
            assign_elem_nodes[0] = 0;
            assign_elem_nodes[1] = 3;
            assign_elem_nodes[2] = 2;
            assign_elem_nodes[3] = 1;
            break;
          }

        case QUAD8:
          {
            fe_descriptor_id = 45; // Plane Stress Quadratic Quadrilateral
            assign_elem_nodes[0] = 0;
            assign_elem_nodes[1] = 7;
            assign_elem_nodes[2] = 3;
            assign_elem_nodes[3] = 6;
            assign_elem_nodes[4] = 2;
            assign_elem_nodes[5] = 5;
            assign_elem_nodes[6] = 1;
            assign_elem_nodes[7] = 4;
            break;
          }

        case QUAD9:
          {
            fe_descriptor_id = 300; // Plane Stress Quadratic Quadrilateral
            assign_elem_nodes[0] = 0;
            assign_elem_nodes[1] = 7;
            assign_elem_nodes[2] = 3;
            assign_elem_nodes[3] = 6;
            assign_elem_nodes[4] = 2;
            assign_elem_nodes[5] = 5;
            assign_elem_nodes[6] = 1;
            assign_elem_nodes[7] = 4;
            assign_elem_nodes[8] = 8;
            break;
          }

        case TET4:
          {
            fe_descriptor_id = 111; // Solid Linear Tetrahedron
            assign_elem_nodes[0] = 0;
            assign_elem_nodes[1] = 1;
            assign_elem_nodes[2] = 2;
            assign_elem_nodes[3] = 3;
            break;
          }

        case PRISM6:
          {
            fe_descriptor_id = 112; // Solid Linear Prism
            assign_elem_nodes[0] = 0;
            assign_elem_nodes[1] = 1;
            assign_elem_nodes[2] = 2;
            assign_elem_nodes[3] = 3;
            assign_elem_nodes[4] = 4;
            assign_elem_nodes[5] = 5;
            break;
          }

        case HEX8:
          {
            fe_descriptor_id = 115; // Solid Linear Brick
            assign_elem_nodes[0] = 0;
            assign_elem_nodes[1] = 4;
            assign_elem_nodes[2] = 5;
            assign_elem_nodes[3] = 1;
            assign_elem_nodes[4] = 3;
            assign_elem_nodes[5] = 7;
            assign_elem_nodes[6] = 6;
            assign_elem_nodes[7] = 2;
            break;
          }

        case HEX20:
          {
            fe_descriptor_id = 116; // Solid Quadratic Brick
            assign_elem_nodes[ 0] = 0;
            assign_elem_nodes[ 1] = 12;
            assign_elem_nodes[ 2] = 4;
            assign_elem_nodes[ 3] = 16;
            assign_elem_nodes[ 4] = 5;
            assign_elem_nodes[ 5] = 13;
            assign_elem_nodes[ 6] = 1;
            assign_elem_nodes[ 7] = 8;

            assign_elem_nodes[ 8] = 11;
            assign_elem_nodes[ 9] = 19;
            assign_elem_nodes[10] = 17;
            assign_elem_nodes[11] = 9;

            assign_elem_nodes[12] = 3;
            assign_elem_nodes[13] = 15;
            assign_elem_nodes[14] = 7;
            assign_elem_nodes[15] = 18;
            assign_elem_nodes[16] = 6;
            assign_elem_nodes[17] = 14;
            assign_elem_nodes[18] = 2;
            assign_elem_nodes[19] = 10;


            break;
          }

        case TET10:
          {
            fe_descriptor_id = 118; // Solid Quadratic Tetrahedron
            assign_elem_nodes[0] = 0;
            assign_elem_nodes[1] = 4;
            assign_elem_nodes[2] = 1;
            assign_elem_nodes[3] = 5;
            assign_elem_nodes[4] = 2;
            assign_elem_nodes[5] = 6;
            assign_elem_nodes[6] = 7;
            assign_elem_nodes[7] = 8;
            assign_elem_nodes[8] = 9;
            assign_elem_nodes[9] = 3;
            break;
          }

        default:
          libmesh_error_msg("ERROR: Element type = " << elem->type() << " not supported in " << "UNVIO!");
        }

      dof_id_type elem_id = elem->id();

      out_file << std::setw(10) << elem_id             // element ID
               << std::setw(10) << fe_descriptor_id    // type of element
               << std::setw(10) << phys_prop_tab_dummy // not supported
               << std::setw(10) << mat_prop_tab_dummy  // not supported
               << std::setw(10) << color_dummy         // not supported
               << std::setw(10) << elem->n_nodes()     // No. of nodes per element
               << '\n';

      for (unsigned int j=0; j<elem->n_nodes(); j++)
        {
          // assign_elem_nodes[j]-th node: i.e., j loops over the
          // libMesh numbering, and assign_elem_nodes[j] over the
          // UNV numbering.
          const Node * node_in_unv_order = elem->node_ptr(assign_elem_nodes[j]);

          // new record after 8 id entries
          if (j==8 || j==16)
            out_file << '\n';

          // Write foreign label for this node
          dof_id_type node_id = node_in_unv_order->id();

          out_file << std::setw(10) << node_id;
        }

      out_file << '\n';

      n_elem_written++;
    }

  if (this->verbose())
    libMesh::out << "  Finished writing " << n_elem_written << " elements" << std::endl;

  // Write end of dataset
  out_file << "    -1\n";
}



void UNVIO::read_dataset(std::string file_name)
{
  std::ifstream in_stream(file_name.c_str());

  if (!in_stream.good())
    libmesh_error_msg("Error opening UNV data file.");

  std::string olds, news, dummy;

  while (true)
    {
      in_stream >> olds >> news;

      // A "-1" followed by a number means the beginning of a dataset.
      while (((olds != "-1") || (news == "-1")) && !in_stream.eof())
        {
          olds = news;
          in_stream >> news;
        }

      if (in_stream.eof())
        break;

      // We only support reading datasets in the "2414" format.
      if (news == "2414")
        {
          // Ignore the rest of the current line and the next two
          // lines that contain analysis dataset label and name.
          for (unsigned int i=0; i<3; i++)
            std::getline(in_stream, dummy);

          // Read the dataset location, where
          // 1: Data at nodes
          // 2: Data on elements
          // Currently only data on nodes is supported.
          unsigned int dataset_location;
          in_stream >> dataset_location;

          // Currently only nodal datasets are supported.
          if (dataset_location != 1)
            libmesh_error_msg("ERROR: Currently only Data at nodes is supported.");

          // Ignore the rest of this line and the next five records.
          for (unsigned int i=0; i<6; i++)
            std::getline(in_stream, dummy);

          // These data are all of no interest to us...
          unsigned int
            model_type,
            analysis_type,
            data_characteristic,
            result_type;

          // The type of data (complex, real, float, double etc, see
          // below).
          unsigned int data_type;

          // The number of floating-point values per entity.
          unsigned int num_vals_per_node;

          in_stream >> model_type           // not used here
                    >> analysis_type        // not used here
                    >> data_characteristic  // not used here
                    >> result_type          // not used here
                    >> data_type
                    >> num_vals_per_node;

          // Ignore the rest of record 9, and records 10-13.
          for (unsigned int i=0; i<5; i++)
            std::getline(in_stream, dummy);

          // Now get the foreign (aka UNV node) node number and
          // the respective nodal data.
          int f_n_id;
          std::vector<Number> values;

          while (true)
            {
              in_stream >> f_n_id;

              // If -1 then we have reached the end of the dataset.
              if (f_n_id == -1)
                break;

              // Resize the values vector (usually data in three
              // principle directions, i.e. num_vals_per_node = 3).
              values.resize(num_vals_per_node);

              // Read the values for the respective node.
              for (unsigned int data_cnt=0; data_cnt<num_vals_per_node; data_cnt++)
                {
                  // Check what data type we are reading.
                  // 2,4: Real
                  // 5,6: Complex
                  // other data types are not supported yet.
                  // As again, these floats may also be written
                  // using a 'D' instead of an 'e'.
                  if (data_type == 2 || data_type == 4)
                    {
                      std::string buf;
                      in_stream >> buf;
                      this->need_D_to_e(buf);
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
                      values[data_cnt] = Complex(std::atof(buf.c_str()), 0.);
#else
                      values[data_cnt] = std::atof(buf.c_str());
#endif
                    }

                  else if (data_type == 5 || data_type == 6)
                    {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
                      Real re_val, im_val;

                      std::string buf;
                      in_stream >> buf;

                      if (this->need_D_to_e(buf))
                        {
                          re_val = std::atof(buf.c_str());
                          in_stream >> buf;
                          this->need_D_to_e(buf);
                          im_val = std::atof(buf.c_str());
                        }
                      else
                        {
                          re_val = std::atof(buf.c_str());
                          in_stream >> im_val;
                        }

                      values[data_cnt] = Complex(re_val,im_val);
#else

                      libmesh_error_msg("ERROR: Complex data only supported when libMesh is configured with --enable-complex!");
#endif
                    }

                  else
                    libmesh_error_msg("ERROR: Data type not supported.");

                } // end loop data_cnt

              // Get a pointer to the Node associated with the UNV node id.
              std::map<dof_id_type, Node *>::const_iterator it =
                _unv_node_id_to_libmesh_node_ptr.find(f_n_id);

              if (it == _unv_node_id_to_libmesh_node_ptr.end())
                libmesh_error_msg("UNV node id " << f_n_id << " was not found.");

              // Store the nodal values in our map which uses the
              // libMesh Node* as the key.  We use operator[] here
              // because we want to create an empty vector if the
              // entry does not already exist.
              _node_data[it->second] = values;
            } // end while (true)
        } // end if (news == "2414")
    } // end while (true)
}



const std::vector<Number> *
UNVIO::get_data (Node * node) const
{
  std::map<Node *, std::vector<Number> >::const_iterator
    it = _node_data.find(node);

  if (it == _node_data.end())
    return libmesh_nullptr;
  else
    return &(it->second);
}


} // namespace libMesh
