// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <iomanip>
#include <cstdio>   // for std::sprintf
#include <algorithm> // for std::sort
#include <fstream>
#include <ctype.h> // isspace

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/unv_io.h"
#include "libmesh/mesh_data.h"
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

#ifdef LIBMESH_HAVE_GZSTREAM
# include "gzstream.h" // For reading/writing compressed streams
#endif



namespace libMesh
{



//-----------------------------------------------------------------------------
// UNVIO class static members
const std::string UNVIO::_label_dataset_nodes    = "2411";
const std::string UNVIO::_label_dataset_elements = "2412";
const std::string UNVIO::_label_dataset_groups   = "2467";



// ------------------------------------------------------------
// UNVIO class members

UNVIO::UNVIO (MeshBase& mesh, MeshData& mesh_data) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _verbose (false),
  _mesh_data (mesh_data)
{
}



UNVIO::UNVIO (const MeshBase& mesh, MeshData& mesh_data) :
  MeshOutput<MeshBase> (mesh),
  _verbose (false),
  _mesh_data (mesh_data)
{
}



UNVIO::~UNVIO ()
{
  this->clear ();
}



bool & UNVIO::verbose ()
{
  return _verbose;
}



bool UNVIO::beginning_of_dataset (std::istream& in_file,
                                  const std::string& ds_name) const
{
  libmesh_assert (in_file.good());
  libmesh_assert (!ds_name.empty());

  std::string olds, news;

  while (true)
    {
      in_file >> olds >> news;

      /*
       * a "-1" followed by a number means the beginning of a dataset
       * stop combing at the end of the file
       */
      while( ((olds != "-1") || (news == "-1") ) && !in_file.eof() )
        {
          olds = news;
          in_file >> news;
        }

      if (in_file.eof())
        return false;

      if (news == ds_name)
        return true;
    }

  // should never end up here
  libmesh_error_msg("We'll never get here!");
  return false;
}



inline
Real UNVIO::D_to_e (std::string& number) const
{
  /* find "D" in string, start looking at
   * 6th element, to improve speed.
   * We dont expect a "D" earlier
   */
  const std::string::size_type position = number.find("D",6);

  libmesh_assert (position != std::string::npos);
  number.replace(position, 1, "e");

  return std::atof (number.c_str());
}



void UNVIO::clear ()
{
  /*
   * Initialize these to dummy values
   */
  this->_n_nodes     = 0;
  this->_n_elements  = 0;
  this->_need_D_to_e = true;

  this->_assign_nodes.clear();
  this->_ds_position.clear();
}





void UNVIO::read (const std::string& file_name)
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


void UNVIO::read_implementation (std::istream& in_stream)
{
  // clear everything, so that
  // we can start from scratch
  this->clear ();

  // Keep track of what kinds of elements this file contains
  elems_of_dimension.clear();
  elems_of_dimension.resize(4, false);

  // Note that we read this file
  // @e twice.  First time to
  // detect the number of nodes
  // and elements (and possible
  // conversion tasks like D_to_e)
  // and the order of datasets
  // (nodes first, then elements,
  // or the other way around),
  // and second to do the actual
  // read.
  std::vector<std::string> order_of_datasets;
  order_of_datasets.reserve(2);

  {
    // the first time we read the file,
    // merely to obtain overall info
    if ( !in_stream.good() )
      libmesh_error_msg("ERROR: Input file not good.");

    // Count nodes and elements, then let
    // other methods read the element and
    // node data.  Also remember which
    // dataset comes first: nodes or elements
    if (this->verbose())
      libMesh::out << "  Counting nodes and elements" << std::endl;

    bool found_node  = false;
    bool found_elem  = false;
    bool found_group = false;

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
            current_line.erase(std::remove_if(current_line.begin(), current_line.end(), isspace), current_line.end());

            // Parse the nodes section
            if (current_line == _label_dataset_nodes &&
                old_line == "-1")
              {
                found_node = true;
                order_of_datasets.push_back (_label_dataset_nodes);
                this->count_nodes (in_stream);
              }

            // Parse the elements section
            else if (current_line == _label_dataset_elements &&
                     old_line == "-1")
              {
                found_elem = true;
                order_of_datasets.push_back (_label_dataset_elements);
                this->count_elements (in_stream);
              }

            // Parse the groups section
            else if (current_line == _label_dataset_groups &&
                     old_line == "-1")
              {
                found_group = true;
                order_of_datasets.push_back (_label_dataset_groups);
                // TODO: We can't currently read the groups until we've
                // read the nodes and elements, unless we want to store
                // the boundary info independently of the BoundaryInfo
                // object.  When we evnetually switch over to not reading
                // the input file twice, we can call groups_in() here.
                // this->groups_in(in_stream);
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
    if (!found_elem)
      libmesh_error_msg("ERROR: Could not find elements!");

    if (!found_node)
      libmesh_error_msg("ERROR: Could not find nodes!");

    // We only support reading the groups *after* the nodes and elements sections
    std::vector<std::string>::iterator
      nodes_pos  = std::find(order_of_datasets.begin(), order_of_datasets.end(), _label_dataset_nodes),
      elems_pos  = std::find(order_of_datasets.begin(), order_of_datasets.end(), _label_dataset_elements),
      groups_pos = std::find(order_of_datasets.begin(), order_of_datasets.end(), _label_dataset_groups);

    if (groups_pos != order_of_datasets.end() &&
        ((nodes_pos > groups_pos) || (elems_pos > groups_pos)))
      libmesh_error_msg("ERROR: Group section must come *after* both nodes and elements sections.");

    // Don't close, just seek to the beginning
    in_stream.seekg(0, std::ios::beg);

    if (!in_stream.good())
      libmesh_error_msg("ERROR: Cannot re-read input file.");
  }





  // We finished scanning the file,
  // and our member data
  // \p this->_n_nodes,
  // \p this->_n_elements,
  // \p this->_need_D_to_e
  // should be properly initialized.
  {
    // Read the datasets in the order that
    // we already know
    libmesh_assert_greater_equal (order_of_datasets.size(), 2);

    for (unsigned int ds=0; ds < order_of_datasets.size(); ds++)
      {
        if (order_of_datasets[ds] == _label_dataset_nodes)
          this->node_in (in_stream);

        else if (order_of_datasets[ds] == _label_dataset_elements)
          this->element_in (in_stream);

        else if (order_of_datasets[ds] == _label_dataset_groups)
          // Reading the groups currently *requires* the nodes and
          // elements to have already been read in.
          this->groups_in(in_stream);

        else
          libmesh_error_msg("Unrecognized dataset type " << order_of_datasets[ds]);
      }

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

    // tell the MeshData object that we are finished
    // reading data
    this->_mesh_data.close_foreign_id_maps ();

    // Delete any lower-dimensional elements that might have been
    // added to the mesh stricly for setting BCs.
    {
      // Grab reference to the Mesh, so we can add boundary info data to it
      MeshBase& mesh = MeshInput<MeshBase>::mesh();

      unsigned max_dim = this->max_elem_dimension_seen();

      MeshBase::const_element_iterator       el     = mesh.elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.elements_end();

      for (; el != end_el; ++el)
        {
          Elem* elem = *el;

          if (elem->dim() < max_dim)
            mesh.delete_elem(elem);
        }
    }


    if (this->verbose())
      libMesh::out << "  Finished." << std::endl << std::endl;
  }

  // save memory
  this->_assign_nodes.clear();
  this->_ds_position.clear();
}





void UNVIO::write (const std::string& file_name)
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




void UNVIO::write_implementation (std::ostream& out_file)
{
  if ( !out_file.good() )
    libmesh_error_msg("ERROR: Output file not good.");


  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  // already know these data, so initialize
  // them.  Does not hurt.
  this->_n_nodes      = mesh.n_nodes();
  this->_n_elements   = mesh.n_elem();
  this->_need_D_to_e  = false;



  // we need the MeshData, otherwise we do not
  // know the foreign node id
  if (!this->_mesh_data.active())
    if (!this->_mesh_data.compatibility_mode())
      {
        libMesh::err << std::endl
                     << "*************************************************************************" << std::endl
                     << "* WARNING: MeshData neither active nor in compatibility mode.           *" << std::endl
                     << "*          Enable compatibility mode for MeshData.  Use this Universal  *" << std::endl
                     << "*          file with caution: libMesh node and element ids are used.    *" << std::endl
                     << "*************************************************************************" << std::endl
                     << std::endl;
        this->_mesh_data.enable_compatibility_mode();
      }



  // write the nodes,  then the elements
  this->node_out    (out_file);
  this->element_out (out_file);
}





void UNVIO::count_nodes (std::istream& in_file)
{
  START_LOG("count_nodes()","UNVIO");

  // if this->_n_nodes is not 0 the dataset
  // has already been scanned
  if (this->_n_nodes != 0)
    libmesh_error_msg("Error: Trying to scan nodes twice!");


  // Read from file, count nodes,
  // check if floats have to be converted
  std::string data;

  in_file >> data; // read the first node label


  if (data == "-1")
    libmesh_error_msg("ERROR: Bad, already reached end of dataset before even starting to read nodes!");


  // ignore the misc data for this node
  in_file.ignore(256,'\n');



  // Now we are there to verify whether we need
  // to convert from D to e or not
  in_file >> data;

  // When this "data" contains a "D", then
  // we have to convert each and every float...
  // But also assume when _this_ specific
  // line does not contain a "D", then the
  // other lines won't, too.
  {
    std::string::size_type position = data.find("D",6);

    if (position!=std::string::npos) // npos means no position
      {
        this->_need_D_to_e = true;

        if (this->verbose())
          libMesh::out << "  Convert from \"D\" to \"e\"" << std::endl;
      }
    else
      this->_need_D_to_e = false;
  }

  // read the remaining two coordinates
  in_file >> data;
  in_file >> data;


  // this was our first node
  this->_n_nodes++;



  // proceed _counting_ the remaining
  // nodes.
  while (in_file.good())
    {
      // read the node label
      in_file >> data;

      if (data == "-1")
        // end of dataset is reached
        break;

      // ignore the remaining data (coord_sys_label, color etc)
      in_file.ignore (256, '\n');
      // ignore the coordinates
      in_file.ignore (256, '\n');

      this->_n_nodes++;
    }


  if (in_file.eof())
    libmesh_error_msg("ERROR: File ended before end of node dataset!");

  if (this->verbose())
    libMesh::out << "  Nodes   : " << this->_n_nodes << std::endl;

  STOP_LOG("count_nodes()","UNVIO");
}






void UNVIO::count_elements (std::istream& in_file)
{
  START_LOG("count_elements()","UNVIO");

  if (this->_n_elements != 0)
    libmesh_error_msg("Error: Trying to scan elements twice!");


  // Simply read the element
  // dataset for the @e only
  // purpose to count nodes!

  std::string data;
  unsigned int fe_id;

  // A dummy string for ignoring lines
  std::string dummy;

  while (!in_file.eof())
    {
      // read element label
      in_file >> data;

      // end of dataset?
      if (data == "-1")
        break;

      // read fe_id
      in_file >> fe_id;

      // Skip related data,
      // and node number list
      std::getline(in_file, dummy);
      std::getline(in_file, dummy);

      // For some elements the node numbers
      // are given more than one record

      // "beam" elements (fe_id < 25) have an extra line of IDs we need to ignore
      // TET10 or QUAD9 - have so many nodes we need to skip an extra line
      if (fe_id < 25 || fe_id == 118 || fe_id == 300)
        std::getline(in_file, dummy);

      // HEX20 - has so many nodes we have to skip two more lines
      if (fe_id == 116)
        {
          std::getline(in_file, dummy);
          std::getline(in_file, dummy);
        }

      this->_n_elements++;
    }


  if (in_file.eof())
    libmesh_error_msg("ERROR: File ended before end of element dataset!");

  if (this->verbose())
    libMesh::out << "  Elements: " << this->_n_elements << std::endl;

  STOP_LOG("count_elements()","UNVIO");
}



void UNVIO::node_in (std::istream& in_file)
{
  START_LOG("node_in()","UNVIO");

  if (this->verbose())
    libMesh::out << "  Reading nodes" << std::endl;

  // adjust the \p istream to our position
  bool ok = this->beginning_of_dataset(in_file, _label_dataset_nodes);

  if (!ok)
    libmesh_error_msg("ERROR: Could not find node dataset!");

  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  // node label, we use an int here so we can read in a -1
  int node_label;

  unsigned int
    exp_coord_sys_num,  // export coordinate system number       (not used yet)
    disp_coord_sys_num, // displacement coordinate system number (not used yet)
    color;              // color                                 (not used yet)

  // always 3 coordinates in the UNV file, no matter
  // which dimensionality libMesh is in
  Point xyz;

  // Continue reading nodes until there are none left
  unsigned ctr = 0;
  while (true)
    {
      // Read the node label
      in_file >> node_label;

      // Break out of the while loop when we hit -1
      if (node_label == -1)
        break;

      // Read the rest of the node data
      in_file >> exp_coord_sys_num       // (not used yet)
              >> disp_coord_sys_num      // (not used yet)
              >> color;                  // (not used yet)

      // read the floating-point data
      for (unsigned int d=0; d<3; d++)
        {
          if (this->_need_D_to_e)
            {
              // Note that \p count_nodes() already determined
              // whether this file uses "D" of "e".
              std::string num_buf;
              in_file >> num_buf;
              xyz(d) = this->D_to_e (num_buf);
            }
          else
            in_file >> xyz(d);
        }

      // set up the id map
      this->_assign_nodes.push_back (node_label);

      // add node to the Mesh &
      // tell the MeshData object the foreign node id
      // (note that mesh.add_point() returns a pointer to the new node)
      this->_mesh_data.add_foreign_node_id (mesh.add_point(xyz, ctr++), node_label);
    }

  // now we need to sort the _assign_nodes vector so we can
  // search it efficiently like a map
  std::sort (this->_assign_nodes.begin(), this->_assign_nodes.end());

  STOP_LOG("node_in()","UNVIO");
}



unsigned UNVIO::max_elem_dimension_seen ()
{
  unsigned max_dim = 0;

  // The elems_of_dimension array is 1-based in the UNV reader
  for (unsigned i=1; i<elems_of_dimension.size(); ++i)
    if (elems_of_dimension[i])
      max_dim = i;

  // Debugging:
  // libMesh::out << "max_dim=" << max_dim << std::endl;

  return max_dim;
}



void UNVIO::groups_in (std::istream& in_file)
{
  // Grab reference to the Mesh, so we can add boundary info data to it
  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  // Read the "-1" strings that divide each section, and the dataset
  // group identifier, so that we are ready to read the actual groups.
  bool ok = this->beginning_of_dataset(in_file, _label_dataset_groups);

  if (!ok)
    libmesh_error_msg("ERROR: Could not find groups dataset!");

  // Record the max and min element dimension seen while reading the file.
  unsigned max_dim = this->max_elem_dimension_seen();

  // map from (node ids) to elem of lower dimensional elements that can provide boundary conditions
  typedef std::map<std::vector<dof_id_type>, Elem*> provide_bcs_t;
  provide_bcs_t provide_bcs;

  // Read groups until there aren't any more to read...
  while (true)
    {
      // If we read a -1, it means there is nothing else to read in this section.
      int group_number;
      in_file >> group_number;

      if (group_number == -1)
        {
          // Do we need to finish reading the rest of the line to the carriage return?
          break;
        }

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

        // Debugging:
        // libMesh::out << "group_number=" << group_number << std::endl;
        // libMesh::out << "num_entities=" << num_entities << std::endl;

        // The second record has 1 field, the group name
        in_file >> group_name;

        // libMesh::out << "group_name=" << group_name << std::endl;
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
                Elem* group_elem = mesh.elem(libmesh_elem_id);
                if (!group_elem)
                  libmesh_error_msg("Group referred to non-existent element with ID " << libmesh_elem_id);

                // dim < max_dim means the Elem defines a boundary condition
                if (group_elem->dim() < max_dim)
                  {
                    is_sideset_group = true;

                    // We can only handle elements that are *exactly*
                    // one dimension lower than the max element
                    // dimension.  Not sure if "edge" BCs in 3D
                    // actually make sense/are required...
                    if (group_elem->dim() != max_dim-1)
                      libmesh_error_msg("ERROR: Expected boundary element of dimension " << max_dim-1 << " but got " << group_elem->dim());

                    // libMesh::out << "UNV Element "
                    //              << entity_tag
                    //              << "(libmesh element "
                    //              << libmesh_elem_id
                    //              << ") dim=="
                    //              << group_elem->dim()
                    //              << " is in group "
                    //              << group_name
                    //              << std::endl;

                    // To be pushed into the provide_bcs data container
                    std::vector<dof_id_type> group_elem_node_ids(group_elem->n_nodes());

                    // Save node IDs in a local vector which will be used as a key for the map.
                    for (unsigned n=0; n<group_elem->n_nodes(); n++)
                      group_elem_node_ids[n] = group_elem->node(n);

                    // Set the current group number as the lower-dimensional element's subdomain ID.
                    // We will use this later to set a boundary ID.
                    group_elem->subdomain_id() = group_number;

                    // Sort before putting into the map
                    std::sort(group_elem_node_ids.begin(), group_elem_node_ids.end());

                    // Unchecked insert:
                    // provide_bcs[group_elem_node_ids] = group_elem;

                    // Ensure that this key doesn't already exist when
                    // inserting it.  We would need to use a multimap if
                    // the same element appears in multiple group numbers!
                    // This actually seems like it could be relatively
                    // common, but I don't have a test for it at the
                    // moment...
                    std::pair<provide_bcs_t::iterator, bool> result =
                      provide_bcs.insert(std::make_pair(group_elem_node_ids, group_elem));

                    if (!result.second)
                      libmesh_error_msg("Boundary element " << group_elem->id() << " was not inserted, it was already in the map!");
                  }

                // dim == max_dim means this group defines a subdomain ID
                else if (group_elem->dim() == max_dim)
                  {
                    is_subdomain_group = true;
                    group_elem->subdomain_id() = group_number;
                  }

                else
                  libmesh_error_msg("ERROR: Found an elem with dim=" << group_elem->dim() << " > " << "max_dim=" << max_dim);
              }
            else
              libMesh::err << "WARNING: UNV Element " << entity_tag << " not found while parsing groups." << std::endl;
          } // end for (entity)
      } // end scope

      // Associate this group_number with the group_name in the BoundaryInfo object.
      if (is_sideset_group)
        mesh.boundary_info->sideset_name(group_number) = group_name;

      if (is_subdomain_group)
        mesh.subdomain_name(group_number) = group_name;

    } // end while (true)

  // Debugging: What did we put in the provide_bcs data structure?
  // {
  //   provide_bcs_t::iterator
  //     provide_it = provide_bcs.begin(),
  //     provide_end = provide_bcs.end();
  //
  //   for ( ; provide_it != provide_end; ++provide_it)
  //     {
  //       const std::vector<dof_id_type> & node_list = provide_it->first;
  //       Elem* elem = provide_it->second;
  //
  //       libMesh::out << "Elem " << elem->id() << " provides BCs for the face with nodes: ";
  //       for (unsigned i=0; i<node_list.size(); ++i)
  //         libMesh::out << node_list[i] << " ";
  //       libMesh::out << std::endl;
  //     }
  // }


  // Loop over elements and try to assign boundary information
  {
    MeshBase::element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::element_iterator end = mesh.active_elements_end();
    for ( ; it != end; ++it)
      {
        Elem* elem = *it;

        if (elem->dim() == max_dim)
          {
            // This is a max-dimension element that
            // may require BCs.  For each of its
            // sides, including internal sides, we'll
            // see if a lower-dimensional element
            // provides boundary information for it.
            // Note that we have not yet called
            // find_neighbors(), so we can't use
            // elem->neighbor(sn) in this algorithm...

            for (unsigned int sn=0; sn<elem->n_sides(); sn++)
              {
                AutoPtr<Elem> side (elem->build_side(sn));

                // Build up a node_ids vector, which is the key
                std::vector<dof_id_type> node_ids(side->n_nodes());
                for (unsigned n=0; n<side->n_nodes(); n++)
                  node_ids[n] = side->node(n);

                // Sort the vector before using it as a key
                std::sort(node_ids.begin(), node_ids.end());

                // Look for this key in the provide_bcs map
                provide_bcs_t::iterator iter = provide_bcs.find(node_ids);

                if (iter != provide_bcs.end())
                  {
                    Elem* lower_dim_elem = iter->second;

                    // Debugging
                    // libMesh::out << "Elem "
                    //              << lower_dim_elem->id()
                    //              << " provides BCs for side "
                    //              << sn
                    //              << " of Elem "
                    //              << elem->id()
                    //              << std::endl;

                    // Add boundary information based on the lower-dimensional element's subdomain id.
                    mesh.boundary_info->add_side(elem, sn, lower_dim_elem->subdomain_id());
                  }
              }
          }
      }
  }

}



void UNVIO::element_in (std::istream& in_file)
{
  START_LOG("element_in()","UNVIO");

  if (this->verbose())
    libMesh::out << "  Reading elements" << std::endl;

  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  // adjust the \p istream to our
  // position
  const bool ok = this->beginning_of_dataset(in_file, _label_dataset_elements);

  if (!ok)
    libmesh_error_msg("ERROR: Could not find element dataset!");

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

  // Get the beginning and end of the _assign_nodes vector
  // to eliminate repeated function calls
  const std::vector<dof_id_type>::const_iterator
    it_begin = this->_assign_nodes.begin(),
    it_end = this->_assign_nodes.end();

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
      Elem* elem = NULL;

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
          // Find the position of node_labels[j] in the _assign_nodes vector.
          const std::pair<std::vector<dof_id_type>::const_iterator, std::vector<dof_id_type>::const_iterator>
            range = std::equal_range (it_begin, it_end, node_labels[j]);

          // it better be there, so libmesh_assert that it was found.
          libmesh_assert (range.first != range.second);
          libmesh_assert_equal_to (*(range.first), node_labels[j]);

          // Now, the distance between this UNV id and the beginning of
          // the _assign_nodes vector will give us a unique id in the
          // range [0,n_nodes) that we can use for defining a contiguous
          // connectivity.
          const dof_id_type assigned_node = libmesh_cast_int<dof_id_type>(std::distance (it_begin, range.first));

          // Make sure we didn't get an out-of-bounds id
          libmesh_assert_less (assigned_node, this->_n_nodes);

          elem->set_node(assign_elem_nodes[j]) = mesh.node_ptr(assigned_node);
        }

      elems_of_dimension[elem->dim()] = true;

      // Set the element's ID
      elem->set_id(ctr);

      // Maintain a map from the libmesh (0-based) numbering to the
      // UNV numbering.  This probably duplicates what the MeshData
      // object does, but hopefully the MeshData object will be going
      // away at some point...
      //_libmesh_elem_id_to_unv_elem_id[i] = element_label;
      _unv_elem_id_to_libmesh_elem_id[element_label] = ctr;

      // Tell the MeshData object the foreign elem id
      this->_mesh_data.add_foreign_elem_id (mesh.add_elem(elem), element_label);

      // Increment the counter for the next iteration
      ctr++;
    } // end while(true)

  STOP_LOG("element_in()","UNVIO");
}






void UNVIO::node_out (std::ostream& out_file)
{

  libmesh_assert (this->_mesh_data.active() ||
                  this->_mesh_data.compatibility_mode());


  if (this->verbose())
    libMesh::out << "  Writing " << this->_n_nodes << " nodes" << std::endl;

  // Write beginning of dataset
  out_file << "    -1\n"
           << "  "
           << _label_dataset_nodes
           << '\n';


  unsigned int exp_coord_sys_dummy  = 0; // export coordinate sys. (not supported yet)
  unsigned int disp_coord_sys_dummy = 0; // displacement coordinate sys. (not supp. yet)
  unsigned int color_dummy          = 0; // color(not supported yet)

  // A reference to the parent class's mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  MeshBase::const_node_iterator       nd  = mesh.nodes_begin();
  const MeshBase::const_node_iterator end = mesh.nodes_end();

  for (; nd != end; ++nd)
    {
      const Node* current_node = *nd;

      char buf[78];
      std::sprintf(buf, "%10d%10u%10u%10u\n",
                   this->_mesh_data.node_to_foreign_id(current_node),
                   exp_coord_sys_dummy,
                   disp_coord_sys_dummy,
                   color_dummy);
      out_file << buf;

      // the coordinates
      if (mesh.spatial_dimension() == 3)
        std::sprintf(buf, "%25.16E%25.16E%25.16E\n",
                     static_cast<double>((*current_node)(0)),
                     static_cast<double>((*current_node)(1)),
                     static_cast<double>((*current_node)(2)));
      else if (mesh.spatial_dimension() == 2)
        std::sprintf(buf, "%25.16E%25.16E\n",
                     static_cast<double>((*current_node)(0)),
                     static_cast<double>((*current_node)(1)));
      else
        std::sprintf(buf, "%25.16E\n",
                     static_cast<double>((*current_node)(0)));

      out_file << buf;
    }


  // Write end of dataset
  out_file << "    -1\n";
}






void UNVIO::element_out(std::ostream& out_file)
{
  libmesh_assert (this->_mesh_data.active() ||
                  this->_mesh_data.compatibility_mode());

  if (this->verbose())
    libMesh::out << "  Writing elements" << std::endl;

  // Write beginning of dataset
  out_file << "    -1\n"
           << "  "
           << _label_dataset_elements
           << "\n";

  unsigned long int fe_descriptor_id = 0;    // FE descriptor id
  unsigned long int phys_prop_tab_dummy = 2; // physical property (not supported yet)
  unsigned long int mat_prop_tab_dummy = 1;  // material property (not supported yet)
  unsigned long int color_dummy = 7;         // color (not supported yet)


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
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  MeshBase::const_element_iterator it  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  for (; it != end; ++it)
    {
      const Elem* elem = *it;

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


      out_file << std::setw(10) << this->_mesh_data.elem_to_foreign_id(elem)  // element ID
               << std::setw(10) << fe_descriptor_id                           // type of element
               << std::setw(10) << phys_prop_tab_dummy                        // not supported
               << std::setw(10) << mat_prop_tab_dummy                         // not supported
               << std::setw(10) << color_dummy                                // not supported
               << std::setw(10) << elem->n_nodes()                            // No. of nodes per element
               << '\n';

      for (unsigned int j=0; j<elem->n_nodes(); j++)
        {
          // assign_elem_nodes[j]-th node: i.e., j loops over the
          // libMesh numbering, and assign_elem_nodes[j] over the
          // UNV numbering.
          const Node* node_in_unv_order = elem->get_node(assign_elem_nodes[j]);

          // new record after 8 id entries
          if (j==8 || j==16)
            out_file << '\n';

          // write foreign label for this node
          out_file << std::setw(10) << this->_mesh_data.node_to_foreign_id(node_in_unv_order);


        }

      out_file << '\n';

      n_elem_written++;
    }

  if (this->verbose())
    libMesh::out << "  Finished writing " << n_elem_written << " elements" << std::endl;

  // Write end of dataset
  out_file << "    -1\n";
}

} // namespace libMesh
