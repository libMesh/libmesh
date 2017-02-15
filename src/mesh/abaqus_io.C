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
#include <string>
#include <cstdlib> // std::strtol
#include <sstream>
#include <ctype.h> // isspace

// Local includes
#include "libmesh/abaqus_io.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h"
#include LIBMESH_INCLUDE_UNORDERED_MAP

// Anonymous namespace to hold mapping Data for Abaqus/libMesh element types
namespace
{
using namespace libMesh;

/**
 * Data structure used for mapping Abaqus IDs to libMesh IDs, and
 * eventually (possibly) vice-versa.
 */
struct ElementDefinition
{
  // Maps (zero-based!) Abaqus local node numbers to libmesh local node numbers
  std::vector<unsigned> abaqus_zero_based_node_id_to_libmesh_node_id;

  // Maps (zero-based!) Abaqus side numbers to libmesh side numbers
  std::vector<unsigned short> abaqus_zero_based_side_id_to_libmesh_side_id;
};

/**
 * Locally-available map containing all element data.
 */
std::map<ElemType, ElementDefinition> eletypes;

/**
 * Helper function to fill up eletypes map
 */
void add_eletype_entry(ElemType libmesh_elem_type,
                       const unsigned * node_map,
                       unsigned node_map_size,
                       const unsigned short * side_map,
                       unsigned side_map_size)
{
  // If map entry does not exist, this will create it
  ElementDefinition & map_entry = eletypes[libmesh_elem_type];


  // Use the "swap trick" from Scott Meyer's "Effective STL" to swap
  // an unnamed temporary vector into the map_entry's vector.  Note:
  // the vector(iter, iter) constructor is used.
  std::vector<unsigned>
    (node_map, node_map+node_map_size).swap
    (map_entry.abaqus_zero_based_node_id_to_libmesh_node_id);

  std::vector<unsigned short>
    (side_map, side_map+side_map_size).swap
    (map_entry.abaqus_zero_based_side_id_to_libmesh_side_id);
}


/**
 * Helper function to initialize the eletypes map.
 */
void init_eletypes ()
{
  // This should happen only once.  The first time this method is
  // called the eletypes data struture will be empty, and we will
  // fill it.  Any subsequent calls will find an initialized
  // eletypes map and will do nothing.
  if (eletypes.empty())
    {
      {
        // EDGE2
        const unsigned int   node_map[] = {0,1}; // identity
        const unsigned short side_map[] = {0,1}; // identity
        add_eletype_entry(EDGE2, node_map, 2, side_map, 2);
      }

      {
        // TRI3
        const unsigned int   node_map[] = {0,1,2}; // identity
        const unsigned short side_map[] = {0,1,2}; // identity
        add_eletype_entry(TRI3, node_map, 3, side_map, 3);
      }

      {
        // QUAD4
        const unsigned int   node_map[] = {0,1,2,3}; // identity
        const unsigned short side_map[] = {0,1,2,3}; // identity
        add_eletype_entry(QUAD4, node_map, 4, side_map, 4);
      }

      {
        // TET4
        const unsigned int   node_map[] = {0,1,2,3}; // identity
        const unsigned short side_map[] = {0,1,2,3}; // identity
        add_eletype_entry(TET4, node_map, 4, side_map, 4);
      }

      {
        // TET10
        const unsigned int   node_map[] = {0,1,2,3,4,5,6,7,8,9}; // identity
        const unsigned short side_map[] = {0,1,2,3};             // identity
        add_eletype_entry(TET10, node_map, 10, side_map, 4);
      }

      {
        // HEX8
        const unsigned int   node_map[] = {0,1,2,3,4,5,6,7}; // identity
        const unsigned short side_map[] = {0,5,1,2,3,4};     // inverse = 0,2,3,4,5,1
        add_eletype_entry(HEX8, node_map, 8, side_map, 6);
      }

      {
        // HEX20
        const unsigned int   node_map[] = // map is its own inverse
          {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
        const unsigned short side_map[] = // inverse = 0,2,3,4,5,1
          {0,5,1,2,3,4};
        add_eletype_entry(HEX20, node_map, 20, side_map, 6);
      }

      {
        // HEX27
        const unsigned int   node_map[] = // inverse = ...,21,23,24,25,26,22,20
          {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,26,20,25,21,22,23,24};
        const unsigned short side_map[] = // inverse = 0,2,3,4,5,1
          {0,5,1,2,3,4};
        add_eletype_entry(HEX27, node_map, 27, side_map, 6);
      }

      {
        // PRISM6
        const unsigned int   node_map[] = {0,1,2,3,4,5}; // identity
        const unsigned short side_map[] = {0,4,1,2,3};   // inverse = 0,2,3,4,1
        add_eletype_entry(PRISM6, node_map, 6, side_map, 5);
      }

      {
        // PRISM15
        const unsigned int   node_map[] = // map is its own inverse
          {0,1,2,3,4,5,6,7,8,12,13,14,9,10,11};
        const unsigned short side_map[] = // inverse = 0,2,3,4,1
          {0,4,1,2,3};
        add_eletype_entry(PRISM15, node_map, 15, side_map, 5);
      }

      {
        // PRISM18
        const unsigned int   node_map[] = // map is its own inverse
          {0,1,2,3,4,5,6,7,8,12,13,14,9,10,11,15,16,17};
        const unsigned short side_map[] = // inverse = 0,2,3,4,1
          {0,4,1,2,3};
        add_eletype_entry(PRISM18, node_map, 18, side_map, 5);
      }
    } // if (eletypes.empty())
}
} // anonymous namespace





namespace libMesh
{

AbaqusIO::AbaqusIO (MeshBase & mesh_in) :
  MeshInput<MeshBase> (mesh_in),
  build_sidesets_from_nodesets(false),
  _already_seen_part(false)
{
}




AbaqusIO::~AbaqusIO ()
{
}




void AbaqusIO::read (const std::string & fname)
{
  // Get a reference to the mesh we are reading
  MeshBase & the_mesh = MeshInput<MeshBase>::mesh();

  // Clear any existing mesh data
  the_mesh.clear();

  // Open stream for reading
  _in.open(fname.c_str());
  libmesh_assert(_in.good());

  // Initialize the elems_of_dimension array.  We will use this in a
  // "1-based" manner so that elems_of_dimension[d]==true means
  // elements of dimension d have been seen.
  elems_of_dimension.resize(4, false);

  // Read file line-by-line... this is based on a set of different
  // test input files.  I have not looked at the full input file
  // specs for Abaqus.
  std::string s;
  while (true)
    {
      // Try to read something.  This may set EOF!
      std::getline(_in, s);

      if (_in)
        {
          // Process s...
          //
          // There are many sections in Abaqus files, we read some
          // but others are just ignored...  Some sections may occur
          // more than once.  For example for a hybrid grid, you
          // will have multiple *Element sections...

          // Some Abaqus files use all upper-case for section names,
          // so we will just convert s to uppercase
          std::string upper(s);
          std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

          // 0.) Look for the "*Part" section
          if (upper.find("*PART") == static_cast<std::string::size_type>(0))
            {
              if (_already_seen_part)
                libmesh_error_msg("We currently don't support reading Abaqus files with multiple PART sections");

              _already_seen_part = true;
            }

          // 1.) Look for the "*Nodes" section
          if (upper.find("*NODE") == static_cast<std::string::size_type>(0))
            {
              // Some sections that begin with *NODE are actually
              // "*NODE OUTPUT" sections which we want to skip.  I
              // have only seen this with a single space, but it would
              // probably be more robust to remove whitespace before
              // making this check.
              if (upper.find("*NODE OUTPUT") == static_cast<std::string::size_type>(0))
                continue;

              // Some *Node sections also specify an Nset name on the same line.
              // Look for one here.
              std::string nset_name = this->parse_label(s, "nset");

              // Process any lines of comments that may be present
              this->process_and_discard_comments();

              // Read a block of nodes
              this->read_nodes(nset_name);
            }



          // 2.) Look for the "*Element" section
          else if (upper.find("*ELEMENT,") == static_cast<std::string::size_type>(0))
            {
              // Some sections that begin with *ELEMENT are actually
              // "*ELEMENT OUTPUT" sections which we want to skip.  I
              // have only seen this with a single space, but it would
              // probably be more robust to remove whitespace before
              // making this check.
              if (upper.find("*ELEMENT OUTPUT") == static_cast<std::string::size_type>(0))
                continue;

              // Some *Element sections also specify an Elset name on the same line.
              // Look for one here.
              std::string elset_name = this->parse_label(s, "elset");

              // Process any lines of comments that may be present
              this->process_and_discard_comments();

              // Read a block of elements
              this->read_elements(upper, elset_name);
            }



          // 3.) Look for a Nodeset section
          else if (upper.find("*NSET") == static_cast<std::string::size_type>(0))
            {
              std::string nset_name = this->parse_label(s, "nset");

              // I haven't seen an unnamed nset yet, but let's detect it
              // just in case...
              if (nset_name == "")
                libmesh_error_msg("Unnamed nset encountered!");

              // Is this a "generated" nset, i.e. one which has three
              // entries corresponding to (first, last, stride)?
              bool is_generated = this->detect_generated_set(upper);

              // Process any lines of comments that may be present
              this->process_and_discard_comments();

              // Read the IDs, storing them in _nodeset_ids
              if (is_generated)
                this->generate_ids(nset_name, _nodeset_ids);
              else
                this->read_ids(nset_name, _nodeset_ids);
            } // *Nodeset



          // 4.) Look for an Elset section
          else if (upper.find("*ELSET") == static_cast<std::string::size_type>(0))
            {
              std::string elset_name = this->parse_label(s, "elset");

              // I haven't seen an unnamed elset yet, but let's detect it
              // just in case...
              if (elset_name == "")
                libmesh_error_msg("Unnamed elset encountered!");

              // Is this a "generated" elset, i.e. one which has three
              // entries corresponding to (first, last, stride)?
              bool is_generated = this->detect_generated_set(upper);

              // Process any lines of comments that may be present
              this->process_and_discard_comments();

              // Read the IDs, storing them in _elemset_ids
              if (is_generated)
                this->generate_ids(elset_name, _elemset_ids);
              else
                this->read_ids(elset_name, _elemset_ids);
            } // *Elset



          // 5.) Look for a Surface section.  Need to be a little
          // careful, since there are also "surface interaction"
          // sections we don't want to read here.
          else if (upper.find("*SURFACE,") == static_cast<std::string::size_type>(0))
            {
              // Get the name from the Name=Foo label.  This will be the map key.
              std::string sideset_name = this->parse_label(s, "name");

              // Process any lines of comments that may be present
              this->process_and_discard_comments();

              // Read the sideset IDs
              this->read_sideset(sideset_name, _sideset_ids);
            }

          continue;
        } // if (_in)

      // If !file, check to see if EOF was set.  If so, break out
      // of while loop.
      if (_in.eof())
        break;

      // If !in and !in.eof(), stream is in a bad state!
      libmesh_error_msg("Stream is bad! Perhaps the file: " << fname << " does not exist?");
    } // while

  // Set the Mesh dimension based on the highest dimension element seen.
  the_mesh.set_mesh_dimension(this->max_elem_dimension_seen());

  // Set element IDs based on the element sets.
  this->assign_subdomain_ids();

  // Assign nodeset values to the BoundaryInfo object
  this->assign_boundary_node_ids();

  // Assign sideset values in the BoundaryInfo object
  this->assign_sideset_ids();

  // Abaqus files only contain nodesets by default.  To be useful in
  // applying most types of BCs in libmesh, we will definitely need
  // sidesets.  So we can call the new BoundaryInfo function which
  // generates sidesets from nodesets.
  if (build_sidesets_from_nodesets)
    the_mesh.get_boundary_info().build_side_list_from_node_list();

  // Delete lower-dimensional elements from the Mesh.  We assume these
  // were only used for setting BCs, and aren't part of the actual
  // Mesh.
  {
    unsigned char max_dim = this->max_elem_dimension_seen();

    MeshBase::element_iterator       el     = the_mesh.elements_begin();
    const MeshBase::element_iterator end_el = the_mesh.elements_end();

    for (; el != end_el; ++el)
      {
        Elem * elem = *el;

        if (elem->dim() < max_dim)
          the_mesh.delete_elem(elem);
      }
  }
}







void AbaqusIO::read_nodes(std::string nset_name)
{
  // Get a reference to the mesh we are reading
  MeshBase & the_mesh = MeshInput<MeshBase>::mesh();

  // In the input files I have, Abaqus neither tells what
  // the mesh dimension is nor how many nodes it has...
  //
  // The node line format is:
  // id, x, y, z
  // and you do have to parse out the commas.
  // The z-coordinate will only be present for 3D meshes

  // Temporary variables for parsing lines of text
  char c;
  std::string line;

  // Defines the sequential node numbering used by libmesh.  Since
  // there can be multiple *NODE sections in an Abaqus file, we always
  // start our numbering with the number of nodes currently in the
  // Mesh.
  dof_id_type libmesh_node_id = the_mesh.n_nodes();

  // We need to duplicate some of the read_ids code if this *NODE
  // section also defines an NSET.  We'll set up the id_storage
  // pointer and push back IDs into this vector in the loop below...
  std::vector<dof_id_type> * id_storage = libmesh_nullptr;
  if (nset_name != "")
    id_storage = &(_nodeset_ids[nset_name]);

  // We will read nodes until the next line begins with *, since that will be the
  // next section.
  // TODO: Is Abaqus guaranteed to start the line with '*' or can there be leading white space?
  while (_in.peek() != '*' && _in.peek() != EOF)
    {
      // Read an entire line which corresponds to a single point's id
      // and (x,y,z) values.
      std::getline(_in, line);

      // Remove all whitespace characters from the line.  This way we
      // can do the remaining parsing without worrying about tabs,
      // different numbers of spaces, etc.
      line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());

      // Make a stream out of the modified line so we can stream values
      // from it in the usual way.
      std::stringstream ss(line);

      // Values to be read in from file
      dof_id_type abaqus_node_id=0;
      Real x=0, y=0, z=0;

      // Note: we assume *at least* 2D points here, should we worry about
      // trying to read 1D Abaqus meshes?
      ss >> abaqus_node_id >> c >> x >> c >> y;

      // Peek at the next character.  If it is a comma, then there is another
      // value to read!
      if (ss.peek() == ',')
        ss >> c >> z;

      // If this *NODE section defines an NSET, also store the abaqus ID in id_storage
      if (id_storage)
        id_storage->push_back(abaqus_node_id);

      // Set up the abaqus -> libmesh node mapping.  This is usually just the
      // "off-by-one" map, but it doesn't have to be.
      _abaqus_to_libmesh_node_mapping[abaqus_node_id] = libmesh_node_id;

      // Add the point to the mesh using libmesh's numbering,
      // and post-increment the libmesh node counter.
      the_mesh.add_point(Point(x,y,z), libmesh_node_id++);
    } // while
}





void AbaqusIO::read_elements(std::string upper, std::string elset_name)
{
  // Get a reference to the mesh we are reading
  MeshBase & the_mesh = MeshInput<MeshBase>::mesh();

  // initialize the eletypes map (eletypes is a file-global variable)
  init_eletypes();

  ElemType elem_type = INVALID_ELEM;
  unsigned n_nodes_per_elem = 0;

  // Within s, we should have "type=XXXX"
  if (upper.find("T3D2") != std::string::npos)
    {
      elem_type = EDGE2;
      n_nodes_per_elem = 2;
      elems_of_dimension[1] = true;
    }
  else if (upper.find("CPE4") != std::string::npos ||
           upper.find("CPS4") != std::string::npos)
    {
      elem_type = QUAD4;
      n_nodes_per_elem = 4;
      elems_of_dimension[2] = true;
    }
  else if (upper.find("CPS3") != std::string::npos ||
           upper.find("S3") != std::string::npos)
    {
      elem_type = TRI3;
      n_nodes_per_elem = 3;
      elems_of_dimension[2] = true;
    }
  else if (upper.find("C3D8") != std::string::npos)
    {
      elem_type = HEX8;
      n_nodes_per_elem = 8;
      elems_of_dimension[3] = true;
    }
  else if (upper.find("C3D4") != std::string::npos)
    {
      elem_type = TET4;
      n_nodes_per_elem = 4;
      elems_of_dimension[3] = true;
    }
  else if (upper.find("C3D20") != std::string::npos)
    {
      elem_type = HEX20;
      n_nodes_per_elem = 20;
      elems_of_dimension[3] = true;
    }
  else if (upper.find("C3D6") != std::string::npos)
    {
      elem_type = PRISM6;
      n_nodes_per_elem = 6;
      elems_of_dimension[3] = true;
    }
  else if (upper.find("C3D15") != std::string::npos)
    {
      elem_type = PRISM15;
      n_nodes_per_elem = 15;
      elems_of_dimension[3] = true;
    }
  else if (upper.find("C3D10") != std::string::npos)
    {
      elem_type = TET10;
      n_nodes_per_elem = 10;
      elems_of_dimension[3] = true;
    }
  else
    libmesh_error_msg("Unrecognized element type: " << upper);

  // Insert the elem type we detected into the set of all elem types for this mesh
  _elem_types.insert(elem_type);

  // Grab a reference to the element definition for this element type
  const ElementDefinition & eledef = eletypes[elem_type];

  // If the element definition was not found, the call above would have
  // created one with an uninitialized struct.  Check for that here...
  if (eledef.abaqus_zero_based_node_id_to_libmesh_node_id.size() == 0)
    libmesh_error_msg("No Abaqus->LibMesh mapping information for ElemType " \
                      << Utility::enum_to_string(elem_type)             \
                      << "!");

  // We will read elements until the next line begins with *, since that will be the
  // next section.
  while (_in.peek() != '*' && _in.peek() != EOF)
    {
      // Read the element ID, it is the first number on each line.  It is
      // followed by a comma, so read that also.  We will need this ID later
      // when we try to assign subdomain IDs
      dof_id_type abaqus_elem_id = 0;
      char c;
      _in >> abaqus_elem_id >> c;

      // Add an element of the appropriate type to the Mesh.
      Elem * elem = the_mesh.add_elem(Elem::build(elem_type).release());

      // Associate the ID returned from libmesh with the abaqus element ID
      //_libmesh_to_abaqus_elem_mapping[elem->id()] = abaqus_elem_id;
      _abaqus_to_libmesh_elem_mapping[abaqus_elem_id] = elem->id();

      // The count of the total number of IDs read for the current element.
      unsigned id_count=0;

      // Continue reading line-by-line until we have read enough nodes for this element
      while (id_count < n_nodes_per_elem)
        {
          // Read entire line (up to carriage return) of comma-separated values
          std::string csv_line;
          std::getline(_in, csv_line);

          // Create a stream object out of the current line
          std::stringstream line_stream(csv_line);

          // Process the comma-separated values
          std::string cell;
          while (std::getline(line_stream, cell, ','))
            {
              // FIXME: factor out this strtol stuff into a utility function.
              char * endptr;
              dof_id_type abaqus_global_node_id = cast_int<dof_id_type>
                (std::strtol(cell.c_str(), &endptr, /*base=*/10));

              if (abaqus_global_node_id!=0 || cell.c_str() != endptr)
                {
                  // Use the global node number mapping to determine the corresponding libmesh global node id
                  dof_id_type libmesh_global_node_id = _abaqus_to_libmesh_node_mapping[abaqus_global_node_id];

                  // Grab the node pointer from the mesh for this ID
                  Node * node = the_mesh.node_ptr(libmesh_global_node_id);

                  // If node_ptr() returns NULL, it may mean we have not yet read the
                  // *Nodes section, though I assumed that always came before the *Elements section...
                  if (node == libmesh_nullptr)
                    libmesh_error_msg("Error!  Mesh returned NULL Node pointer.  Either no node exists with ID " \
                                      << libmesh_global_node_id         \
                                      << " or perhaps this input file has *Elements defined before *Nodes?");

                  // Note: id_count is the zero-based abaqus (elem local) node index.  We therefore map
                  // it to a libmesh elem local node index using the element definition map
                  unsigned libmesh_elem_local_node_id =
                    eledef.abaqus_zero_based_node_id_to_libmesh_node_id[id_count];

                  // Set this node pointer within the element.
                  elem->set_node(libmesh_elem_local_node_id) = node;

                  // Increment the count of IDs read for this element
                  id_count++;
                } // end if strtol success
            } // end while getline(',')
        } // end while (id_count)

      // Ensure that we read *exactly* as many nodes as we were expecting to, no more.
      if (id_count != n_nodes_per_elem)
        libmesh_error_msg("Error: Needed to read " \
                          << n_nodes_per_elem      \
                          << " nodes, but read "   \
                          << id_count              \
                          << " instead!");

      // If we are recording Elset IDs, add this element to the correct set for later processing.
      // Make sure to add it with the Abaqus ID, not the libmesh one!
      if (elset_name != "")
        _elemset_ids[elset_name].push_back(abaqus_elem_id);
    } // end while (peek)
}




std::string AbaqusIO::parse_label(std::string line, std::string label_name) const
{
  // Handle files which have weird line endings from e.g. windows.
  // You can check what kind of line endings you have with 'cat -vet'.
  // For example, some files may have two kinds of line endings like:
  //
  // 4997,^I496,^I532,^I487,^I948^M$
  //
  // and we don't want to deal with this when extracting a label, so
  // just remove all the space characters, which should include all
  // kinds of remaining newlines.  (I don't think Abaqus allows
  // whitespace in label names.)
  line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());

  // Do all string comparisons in upper-case
  std::string
    upper_line(line),
    upper_label_name(label_name);
  std::transform(upper_line.begin(), upper_line.end(), upper_line.begin(), ::toupper);
  std::transform(upper_label_name.begin(), upper_label_name.end(), upper_label_name.begin(), ::toupper);

  // Get index of start of "label="
  size_t label_index = upper_line.find(upper_label_name + "=");

  if (label_index != std::string::npos)
    {
      // Location of the first comma following "label="
      size_t comma_index = upper_line.find(",", label_index);

      // Construct iterators from which to build the sub-string.
      // Note: The +1 while initializing beg is to skip past the "=" which follows the label name
      std::string::iterator
        beg = line.begin() + label_name.size() + 1 + label_index,
        end = (comma_index == std::string::npos) ? line.end() : line.begin() + comma_index;

      return std::string(beg, end);
    }

  // The label index was not found, return the empty string
  return std::string("");
}




bool AbaqusIO::detect_generated_set(std::string upper) const
{
  // Avoid issues with weird line endings, spaces before commas, etc.
  upper.erase(std::remove_if(upper.begin(), upper.end(), isspace), upper.end());

  // Check each comma-separated value in "upper" to see if it is the generate flag.
  std::string cell;
  std::stringstream line_stream(upper);
  while (std::getline(line_stream, cell, ','))
    if (cell == "GENERATE")
      return true;

  return false;
}



void AbaqusIO::read_ids(std::string set_name, container_t & container)
{
  // Grab a reference to a vector that will hold all the IDs
  std::vector<dof_id_type> & id_storage = container[set_name];

  // Read until the start of another section is detected, or EOF is encountered
  while (_in.peek() != '*' && _in.peek() != EOF)
    {
      // Read entire comma-separated line into a string
      std::string csv_line;
      std::getline(_in, csv_line);

      // On that line, use std::getline again to parse each
      // comma-separated entry.
      std::string cell;
      std::stringstream line_stream(csv_line);
      while (std::getline(line_stream, cell, ','))
        {
          // If no conversion can be performed by strtol, 0 is returned.
          //
          // If endptr is not NULL, strtol() stores the address of the
          // first invalid character in *endptr.  If there were no
          // digits at all, however, strtol() stores the original
          // value of str in *endptr.
          char * endptr;

          // FIXME - this needs to be updated for 64-bit inputs
          dof_id_type id = cast_int<dof_id_type>
            (std::strtol(cell.c_str(), &endptr, /*base=*/10));

          // Note that lists of comma-separated values in abaqus also
          // *end* with a comma, so the last call to getline on a given
          // line will get an empty string, which we must detect.
          if (id != 0 || cell.c_str() != endptr)
            {
              // 'cell' is now a string with an integer id in it
              id_storage.push_back( id );
            }
        }
    }
}




void AbaqusIO::generate_ids(std::string set_name, container_t & container)
{
  // Grab a reference to a vector that will hold all the IDs
  std::vector<dof_id_type> & id_storage = container[set_name];

  // Read until the start of another section is detected, or EOF is
  // encountered.  "generate" sections seem to only have one line,
  // although I suppose it's possible they could have more.
  while (_in.peek() != '*' && _in.peek() != EOF)
    {
      // Read entire comma-separated line into a string
      std::string csv_line;
      std::getline(_in, csv_line);

      // Remove all whitespaces from csv_line.
      csv_line.erase(std::remove_if(csv_line.begin(), csv_line.end(), isspace), csv_line.end());

      // Create a new stringstream object from the string, and stream
      // in the comma-separated values.
      char c;
      dof_id_type start, end, stride;
      std::stringstream line_stream(csv_line);
      line_stream >> start >> c >> end >> c >> stride;

      // Generate entries in the id_storage.  Note: each element can
      // only belong to a single Elset (since this corresponds to the
      // subdomain_id) so if an element appears in multiple Elsets,
      // the "last" one (alphabetically, based on set name) in the
      // _elemset_ids map will "win".
      for (dof_id_type current = start; current <= end; current += stride)
        id_storage.push_back(current);
    }
}




void AbaqusIO::read_sideset(std::string sideset_name, sideset_container_t & container)
{
  // Grab a reference to a vector that will hold all the IDs
  std::vector<std::pair<dof_id_type, unsigned> > & id_storage = container[sideset_name];

  // Variables for storing values read in from file
  dof_id_type elem_id=0;
  unsigned side_id=0;
  char c;
  std::string dummy;

  // Read until the start of another section is detected, or EOF is encountered
  while (_in.peek() != '*' && _in.peek() != EOF)
    {
      // The strings are of the form: "391, S2"

      // Read the element ID and the leading comma
      _in >> elem_id >> c;

      // Read another character (the 'S') and finally the side ID
      _in >> c >> side_id;

      // Store this pair of data in the vector
      id_storage.push_back( std::make_pair(elem_id, side_id) );

      // Extract remaining characters on line including newline
      std::getline(_in, dummy);
    } // while
}




void AbaqusIO::assign_subdomain_ids()
{
  // Get a reference to the mesh we are reading
  MeshBase & the_mesh = MeshInput<MeshBase>::mesh();

  // The number of elemsets we've found while reading
  std::size_t n_elemsets = _elemset_ids.size();

  // Fill in a temporary map with (ElemType, index) pairs based on the _elem_types set.  This
  // will allow us to easily look up this index in the loop below.
  std::map<ElemType, unsigned> elem_types_map;
  {
    unsigned ctr=0;
    for (std::set<ElemType>::iterator it=_elem_types.begin(); it!=_elem_types.end(); ++it)
      elem_types_map[*it] = ctr++;
  }

  // Loop over each Elemset and assign subdomain IDs to Mesh elements
  {
    // The maximum element dimension seen while reading the Mesh
    unsigned char max_dim = this->max_elem_dimension_seen();

    // The elemset_id counter assigns a logical numbering to the _elemset_ids keys
    container_t::iterator it = _elemset_ids.begin();
    for (unsigned elemset_id=0; it != _elemset_ids.end(); ++it, ++elemset_id)
      {
        // Grab a reference to the vector of IDs
        std::vector<dof_id_type> & id_vector = it->second;

        // Loop over this vector
        for (std::size_t i=0; i<id_vector.size(); ++i)
          {
            // Map the id_vector[i]'th element ID (Abaqus numbering) to LibMesh numbering
            dof_id_type libmesh_elem_id = _abaqus_to_libmesh_elem_mapping[ id_vector[i] ];

            // Get reference to that element
            Elem & elem = the_mesh.elem_ref(libmesh_elem_id);

            // We won't assign subdomain ids to lower-dimensional
            // elements, as they are assumed to represent boundary
            // conditions.  Since lower-dimensional elements can
            // appear in multiple sidesets, it doesn't make sense to
            // assign them a single subdomain id... only the last one
            // assigned would actually "stick".
            if (elem.dim() < max_dim)
              break;

            // Compute the proper subdomain ID, based on the formula in the
            // documentation for this function.
            subdomain_id_type computed_id = cast_int<subdomain_id_type>
              (elemset_id + (elem_types_map[elem.type()] * n_elemsets));

            // Assign this ID to the element in question
            elem.subdomain_id() = computed_id;

            // We will also assign a unique name to the computed_id,
            // which is created by appending the geometric element
            // name to the elset name provided by the user in the
            // Abaqus file.
            std::string computed_name = it->first + "_" + Utility::enum_to_string(elem.type());
            the_mesh.subdomain_name(computed_id) = computed_name;
          }
      }
  }
}




void AbaqusIO::assign_boundary_node_ids()
{
  // Get a reference to the mesh we are reading
  MeshBase & the_mesh = MeshInput<MeshBase>::mesh();

  // Iterate over the container of nodesets
  container_t::iterator it = _nodeset_ids.begin();
  for (unsigned short current_id=0; it != _nodeset_ids.end(); ++it, ++current_id)
    {
      // Associate current_id with the name we determined earlier
      the_mesh.get_boundary_info().nodeset_name(current_id) = it->first;

      // Get a reference to the current vector of nodeset ID values
      std::vector<dof_id_type> & nodeset_ids = it->second;

      for (std::size_t i=0; i<nodeset_ids.size(); ++i)
        {
          // Map the Abaqus global node ID to the libmesh node ID
          dof_id_type libmesh_global_node_id = _abaqus_to_libmesh_node_mapping[nodeset_ids[i]];

          // Get node pointer from the mesh
          Node * node = the_mesh.node_ptr(libmesh_global_node_id);

          if (node == libmesh_nullptr)
            libmesh_error_msg("Error! Mesh returned NULL node pointer!");

          // Add this node with the current_id (which is determined by the
          // alphabetical ordering of the map) to the BoundaryInfo object
          the_mesh.get_boundary_info().add_node(node, current_id);
        }
    }
}




void AbaqusIO::assign_sideset_ids()
{
  // Get a reference to the mesh we are reading
  MeshBase & the_mesh = MeshInput<MeshBase>::mesh();

  // initialize the eletypes map (eletypes is a file-global variable)
  init_eletypes();

  // Iterate over the container of sidesets
  {
    sideset_container_t::iterator it = _sideset_ids.begin();
    for (unsigned short current_id=0; it != _sideset_ids.end(); ++it, ++current_id)
      {
        // Associate current_id with the name we determined earlier
        the_mesh.get_boundary_info().sideset_name(current_id) = it->first;

        // Get a reference to the current vector of nodeset ID values
        std::vector<std::pair<dof_id_type,unsigned> > & sideset_ids = it->second;

        for (std::size_t i=0; i<sideset_ids.size(); ++i)
          {
            // sideset_ids is a vector of pairs (elem id, side id).  Pull them out
            // now to make the code below more readable.
            dof_id_type  abaqus_elem_id = sideset_ids[i].first;
            unsigned abaqus_side_number = sideset_ids[i].second;

            // Map the Abaqus element ID to LibMesh numbering
            dof_id_type libmesh_elem_id = _abaqus_to_libmesh_elem_mapping[ abaqus_elem_id ];

            // Get a reference to that element
            Elem & elem = the_mesh.elem_ref(libmesh_elem_id);

            // Grab a reference to the element definition for this element type
            const ElementDefinition & eledef = eletypes[elem.type()];

            // If the element definition was not found, the call above would have
            // created one with an uninitialized struct.  Check for that here...
            if (eledef.abaqus_zero_based_side_id_to_libmesh_side_id.size() == 0)
              libmesh_error_msg("No Abaqus->LibMesh mapping information for ElemType " \
                                << Utility::enum_to_string(elem.type())  \
                                << "!");

            // Add this node with the current_id (which is determined by the
            // alphabetical ordering of the map).  Side numbers in Abaqus are 1-based,
            // so we subtract 1 here before passing the abaqus side number to the
            // mapping array
            the_mesh.get_boundary_info().add_side
              (&elem,
               eledef.abaqus_zero_based_side_id_to_libmesh_side_id[abaqus_side_number-1],
               current_id);
          }
      }
  }


  // Some elsets (if they contain lower-dimensional elements) also
  // define sidesets.  So loop over them and build a searchable data
  // structure we can use to assign sidesets.
  {
    unsigned char max_dim = this->max_elem_dimension_seen();

    // multimap from lower-dimensional-element-hash-key to
    // pair(lower-dimensional-element, boundary_id).  The
    // lower-dimensional element is used to verify the results of the
    // hash table search.  The boundary_id will be used to set a
    // boundary ID on a higher-dimensional element.  We use a multimap
    // because the lower-dimensional elements can belong to more than
    // 1 sideset, and multiple lower-dimensional elements can hash to
    // the same value, but this is very rare.
    typedef LIBMESH_BEST_UNORDERED_MULTIMAP<dof_id_type,
                                            std::pair<Elem *, boundary_id_type> > provide_bcs_t;
    provide_bcs_t provide_bcs;

    // The elemset_id counter assigns a logical numbering to the
    // _elemset_ids keys.  We are going to use these ids as boundary
    // ids, so elemset_id is of type boundary_id_type.
    container_t::iterator it = _elemset_ids.begin();
    for (boundary_id_type elemset_id=0; it != _elemset_ids.end(); ++it, ++elemset_id)
      {
        // Grab a reference to the vector of IDs
        std::vector<dof_id_type> & id_vector = it->second;

        // Loop over this vector
        for (std::size_t i=0; i<id_vector.size(); ++i)
          {
            // Map the id_vector[i]'th element ID (Abaqus numbering) to LibMesh numbering
            dof_id_type libmesh_elem_id = _abaqus_to_libmesh_elem_mapping[ id_vector[i] ];

            // Get a reference to that element
            Elem & elem = the_mesh.elem_ref(libmesh_elem_id);

            // If the element dimension is equal to the maximum
            // dimension seen, we can break out of this for loop --
            // this elset does not contain sideset information.
            if (elem.dim() == max_dim)
              break;

            // We can only handle elements that are *exactly*
            // one dimension lower than the max element
            // dimension.  Not sure if "edge" BCs in 3D
            // actually make sense/are required...
            if (elem.dim()+1 != max_dim)
              libmesh_error_msg("ERROR: Expected boundary element of dimension " << max_dim-1 << " but got " << elem.dim());

            // Insert the current (key, pair(elem,id)) into the multimap.
            provide_bcs.insert(std::make_pair(elem.key(),
                                              std::make_pair(&elem,
                                                             elemset_id)));

            // Associate the name of this sideset with the ID we've
            // chosen.  It's not necessary to do this for every
            // element in the set, but it's convenient to do it here
            // since we have all the necessary information...
            the_mesh.get_boundary_info().sideset_name(elemset_id) = it->first;
          }
      }

    // Loop over elements and try to assign boundary information
    {
      MeshBase::element_iterator       e_it  = the_mesh.active_elements_begin();
      const MeshBase::element_iterator end = the_mesh.active_elements_end();
      for ( ; e_it != end; ++e_it)
        {
          Elem * elem = *e_it;

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
                  std::pair<provide_bcs_t::const_iterator,
                            provide_bcs_t::const_iterator>
                    range = provide_bcs.equal_range (elem->key(sn));

                  // Add boundary information for each side in the range.
                  for (provide_bcs_t::const_iterator s_it = range.first;
                       s_it != range.second; ++s_it)
                    {
                      // We'll need to compare the lower dimensional element against the current side.
                      UniquePtr<Elem> side (elem->build_side_ptr(sn));

                      // Get the value mapped by the iterator.
                      std::pair<Elem *, boundary_id_type> p = s_it->second;

                      // Extract the relevant data from the iterator.
                      Elem * lower_dim_elem = p.first;
                      boundary_id_type bid = p.second;

                      // This was a hash, so it might not be perfect.  Let's verify...
                      if (*lower_dim_elem == *side)
                        the_mesh.get_boundary_info().add_side(elem, sn, bid);
                    }
                }
            }
        }
    }
  }
}



void AbaqusIO::process_and_discard_comments()
{
  std::string dummy;
  while (true)
    {
      // We assume we are at the beginning of a line that may be
      // comments or may be data.  We need to only discard the line if
      // it begins with **, but we must avoid calling std::getline()
      // since there's no way to put that back.
      if (_in.peek() == '*')
        {
          // The first character was a star, so actually read it from the stream.
          _in.get();

          // Peek at the next character...
          if (_in.peek() == '*')
            {
              // OK, second character was star also, by definition this
              // line must be a comment!  Read the rest of the line and discard!
              std::getline(_in, dummy);
            }
          else
            {
              // The second character was _not_ a star, so put back the first star
              // we pulled out so that the line can be parsed correctly by somebody
              // else!
              _in.unget();

              // Finally, break out of the while loop, we are done parsing comments
              break;
            }
        }
      else
        {
          // First character was not *, so this line must be data! Break out of the
          // while loop!
          break;
        }
    }
}



unsigned char AbaqusIO::max_elem_dimension_seen ()
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


} // namespace
