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
#include <string>

// Local includes
#include "abaqus_io.h"
#include "point.h"
#include "elem.h"
#include "string_to_enum.h"
#include "boundary_info.h"
#include "utility.h"

// Anonymous namespace to hold mapping Data for Abaqus/libMesh element types
namespace
{
  /**
   * Data structure used for mapping Abaqus IDs to libMesh IDs, and
   * eventually (possibly) vice-versa.
   */
  struct ElementDefinition
  {
    // Maps (zero-based!) Abaqus node numbers to libmesh node numbers
    std::vector<unsigned> abaqus_zero_based_node_id_to_libmesh_node_id;

    // Maps (zero-based!) Abaqus side numbers to libmesh side numbers
    std::vector<unsigned> abaqus_zero_based_side_id_to_libmesh_side_id;
  };

  /**
   * Locally-available map containing all element data.
   */
  std::map<ElemType, ElementDefinition> eletypes;

  /**
   * Helper function to fill up eletypes map
   */
  void add_eletype_entry(ElemType libmesh_elem_type,
			 const unsigned* node_map,
			 unsigned node_map_size,
			 const unsigned* side_map,
			 unsigned side_map_size)
  {
    // If map entry does not exist, this will create it
    ElementDefinition& map_entry = eletypes[libmesh_elem_type];


    // Use the "swap trick" from Scott Meyer's "Effective STL" to swap
    // an unnamed temporary vector into the map_entry's vector.  Note:
    // the vector(iter, iter) constructor is used.
    std::vector<unsigned>(node_map,
			  node_map+node_map_size).swap(map_entry.abaqus_zero_based_node_id_to_libmesh_node_id);

    std::vector<unsigned>(side_map,
			  side_map+side_map_size).swap(map_entry.abaqus_zero_based_side_id_to_libmesh_side_id);
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
	  // TRI3
	  const unsigned int node_map[] = {0,1,2}; // identity
	  const unsigned int side_map[] = {0,1,2}; // identity
	  add_eletype_entry(TRI3, node_map, 3, side_map, 3);
	}

	{
	  // QUAD4
	  const unsigned int node_map[] = {0,1,2,3}; // identity
	  const unsigned int side_map[] = {0,1,2,3}; // identity
	  add_eletype_entry(QUAD4, node_map, 4, side_map, 4);
	}

	{
	  // TET4
	  const unsigned int node_map[] = {0,1,2,3}; // identity
	  const unsigned int side_map[] = {0,1,2,3}; // identity
	  add_eletype_entry(TET4, node_map, 4, side_map, 4);
	}

	{
	  // TET10
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,9}; // identity
	  const unsigned int side_map[] = {0,1,2,3};             // identity
	  add_eletype_entry(TET10, node_map, 10, side_map, 4);
	}

	{
	  // HEX8
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7}; // identity
	  const unsigned int side_map[] = {0,5,1,2,3,4};     // inverse = 0,2,3,4,5,1
	  add_eletype_entry(HEX8, node_map, 8, side_map, 6);
	}

	{
	  // HEX20
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15}; // map is its own inverse
	  const unsigned int side_map[] = {0,5,1,2,3,4};                                       // inverse = 0,2,3,4,5,1
	  add_eletype_entry(HEX20, node_map, 20, side_map, 6);
	}

	{
	  // HEX27
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,26,20,25,21,22,23,24}; // inverse = ...,21,23,24,25,26,22,20
	  const unsigned int side_map[] = {0,5,1,2,3,4};                                                            // inverse = 0,2,3,4,5,1
	  add_eletype_entry(HEX27, node_map, 27, side_map, 6);
	}

	{
	  // PRISM6
	  const unsigned int node_map[] = {0,1,2,3,4,5}; // identity
	  const unsigned int side_map[] = {0,4,1,2,3};   // inverse = 0,2,3,4,1
	  add_eletype_entry(PRISM6, node_map, 6, side_map, 5);
	}

	{
	  // PRISM15
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,12,13,14,9,10,11}; // map is its own inverse
	  const unsigned int side_map[] = {0,4,1,2,3};                          // inverse = 0,2,3,4,1
	  add_eletype_entry(PRISM15, node_map, 15, side_map, 5);
	}

	{
	  // PRISM18
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,12,13,14,9,10,11,15,16,17}; // map is its own inverse
	  const unsigned int side_map[] = {0,4,1,2,3};                                   // inverse = 0,2,3,4,1
	  add_eletype_entry(PRISM18, node_map, 18, side_map, 5);
	}



      } // if (eletypes.empty())
  } // init_eletypes()
} // anonymous namespace





namespace libMesh
{

  Abaqus_IO::Abaqus_IO (MeshBase& mesh) :
    MeshInput<MeshBase> (mesh),
    build_sidesets_from_nodesets(false),
    _already_seen_part(false),
    _max_csv_per_line(16)
  {
  }




  Abaqus_IO::~Abaqus_IO ()
  {
  }




  void Abaqus_IO::read (const std::string& fname)
  {
    // Get a reference to the mesh we are reading
    MeshBase& mesh = MeshInput<MeshBase>::mesh();

    // Clear any existing mesh data
    mesh.clear();

    // Open stream for reading
    _in.open(fname.c_str());
    libmesh_assert(_in.good());

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
	    if (upper.find("*PART") == 0)
	      {
		// libMesh::out << "Found parts section!" << std::endl;

		if (_already_seen_part)
		  {
		    libMesh::err << "We currently don't support reading Abaqus files with multiple PART sections" << std::endl;
		    libmesh_error();
		  }

		_already_seen_part = true;
	      }

	    // 1.) Look for the "*Nodes" section
	    if (upper.find("*NODE") == 0)
	      {
		// Process any lines of comments that may be present
		this->process_and_discard_comments();

		// Read a block of nodes
		this->read_nodes();
	      }



	    // 2.) Look for the "*Element" section
	    else if (upper.find("*ELEMENT,") == 0)
	      {
		// Process any lines of comments that may be present
		this->process_and_discard_comments();

		// Read a block of elements
		this->read_elements(upper);
	      }



	    // 3.) Look for a Nodeset section
	    else if (upper.find("*NSET") == 0)
	      {
		std::string nset_name = this->parse_label(s, "nset");

		// I haven't seen an unnamed elset yet, but let's detect it
		// just in case...
		if (nset_name == "")
		  {
		    libMesh::err << "Unnamed nset encountered!" << std::endl;
		    libmesh_error();
		  }

		// Process any lines of comments that may be present
		this->process_and_discard_comments();

		// Read the IDs, storing them in _nodeset_ids
		this->read_ids(nset_name, _nodeset_ids);
	      } // *Nodeset



	    // 4.) Look for an Elset section
	    else if (upper.find("*ELSET") == 0)
	      {
		std::string elset_name = this->parse_label(s, "elset");

		// I haven't seen an unnamed elset yet, but let's detect it
		// just in case...
		if (elset_name == "")
		  {
		    libMesh::err << "Unnamed elset encountered!" << std::endl;
		    libmesh_error();
		  }

		// Process any lines of comments that may be present
		this->process_and_discard_comments();

		// Read the IDs, storing them in _elemset_ids
		this->read_ids(elset_name, _elemset_ids);
	      } // *Elset



	    // 5.) Look for a Surface section.  Need to be a little
	    // careful, since there are also "surface interaction"
	    // sections we don't want to read here.
	    else if (upper.find("*SURFACE,") == 0)
	      {
		// libMesh::out << "Found SURFACE section: " << s << std::endl;

		// Get the name from the Name=Foo label.  This will be the map key.
		std::string sideset_name = this->parse_label(s, "name");

		// Print name of section we just found
		// libMesh::out << "Found surface section named: " << sideset_name << std::endl;

		// Process any lines of comments that may be present
		this->process_and_discard_comments();

		// Read the sideset IDs
		this->read_sideset(sideset_name, _sideset_ids);

		// Debugging: print status of most recently read sideset
		// libMesh::out << "Read " << _sideset_ids[sideset_name].size() << " sides in " << sideset_name << std::endl;
	      }


	    continue;
	  } // if (_in)

	// If !file, check to see if EOF was set.  If so, break out
	// of while loop.
	if (_in.eof())
	  break;

	// If !in and !in.eof(), stream is in a bad state!
	libMesh::err << "Stream is bad!" << std::endl;
	libmesh_error();
      } // while


    //
    // We've read everything we can from the file at this point.  Now
    // do some more processing.
    //
    libMesh::out << "Mesh contains "
		 << mesh.n_elem()
		 << " elements, and "
		 << mesh.n_nodes()
		 << " nodes." << std::endl;

    // TODO: Remove these or write a function to do it?
//    {
//      container_t::iterator it=_nodeset_ids.begin();
//      for (; it != _nodeset_ids.end(); ++it)
//	{
//	  libMesh::out << "Node set '" << (*it).first << "' contains " << (*it).second.size() << " ID(s)." << std::endl;
//	}
//    }
//
//    {
//      container_t::iterator it=_elemset_ids.begin();
//      for (; it != _elemset_ids.end(); ++it)
//	{
//	  libMesh::out << "Elem set '" << (*it).first << "' contains " << (*it).second.size() << " ID(s)." << std::endl;
//	}
//    }


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
      mesh.boundary_info->build_side_list_from_node_list();

  } // read()







  void Abaqus_IO::read_nodes()
  {
    // Get a reference to the mesh we are reading
    MeshBase& mesh = MeshInput<MeshBase>::mesh();

    // In the input file I have, Abaqus neither tells what
    // the mesh dimension is nor how many nodes it has...

    // The node line format is:
    // id, x, y, z
    // and you do have to parse out the commas.
    // The z-coordinate will only be present for 3D meshes

    // Temporary variables for parsing lines of text
    unsigned node_id=0;
    Real x=0, y=0, z=0;
    char c;
    std::string dummy;

    // We will read nodes until the next line begins with *, since that will be the
    // next section.
    // TODO: Is Abaqus guaranteed to start the line with '*' or can there be leading white space?
    while (_in.peek() != '*' && _in.peek() != EOF)
      {
	node_id=0;
	x = y = z = 0.;

	// Note: we assume *at least* 2D points here, should we worry about
	// trying to read 1D Abaqus meshes?
	_in >> node_id >> c >> x >> c >> y;

	// Peek at the next character.  If it is a comma, then there is another
	// value to read!
	if (_in.peek() == ',')
	  _in >> c >> z;

	// Print what we just read in.
	// libMesh::out << "node_id=" << node_id
	// 	     << ", x=" << x
	// 	     << ", y=" << y
	// 	     << ", z=" << z
	// 	     << std::endl;

	// Read (and discard) the rest of the line, including the newline.
	// This is required so that our 'peek()' at the beginning of this
	// loop doesn't read the newline character, for example.
	std::getline(_in, dummy);

	// Add the point to the mesh using the Abaqus numbering.
	// Note that the Abaqus node IDs are not in general in order
	// in the input file!
	mesh.add_point(Point(x,y,z), node_id);
      } // while
  } // read_nodes()





  void Abaqus_IO::read_elements(std::string upper)
  {
    // Some *Element sections also specify an Elset name on the same line.
    // Look for one here.
    std::string elset_name = this->parse_label(upper, "elset");

    // Get a reference to the mesh we are reading
    MeshBase& mesh = MeshInput<MeshBase>::mesh();

    // initialize the eletypes map (eletypes is a file-global variable)
    init_eletypes();

    ElemType elem_type = INVALID_ELEM;
    unsigned n_nodes_per_elem = 0;

    // Within s, we should have "type=XXXX"
    if (upper.find("CPE4") != std::string::npos ||
	upper.find("CPS4") != std::string::npos)
      {
	// libMesh::out << "Element type is QUAD4!" << std::endl;
	elem_type = QUAD4;
	n_nodes_per_elem = 4;
	_elem_types.insert(QUAD4);
	mesh.set_mesh_dimension(2);
      }
    else if (upper.find("CPS3") != std::string::npos)
      {
	// libMesh::out << "Element type is TRI3!" << std::endl;
	elem_type = TRI3;
	n_nodes_per_elem = 3;
	_elem_types.insert(TRI3);
	mesh.set_mesh_dimension(2);
      }
    else if (upper.find("C3D8") != std::string::npos)
      {
	// libMesh::out << "Element type is HEX8!" << std::endl;
	elem_type = HEX8;
	n_nodes_per_elem = 8;
	_elem_types.insert(HEX8);
	mesh.set_mesh_dimension(3);
      }
    else if (upper.find("C3D4") != std::string::npos)
      {
	// libMesh::out << "Element type is TET4!" << std::endl;
	elem_type = TET4;
	n_nodes_per_elem = 4;
	_elem_types.insert(TET4);
	mesh.set_mesh_dimension(3);
      }
    else if (upper.find("C3D20") != std::string::npos)
      {
	// libMesh::out << "Element type is HEX20!" << std::endl;
	elem_type = HEX20;
	n_nodes_per_elem = 20;
	_elem_types.insert(HEX20);
	mesh.set_mesh_dimension(3);
      }
    else
      {
	libMesh::err << "Unrecognized element type: " << upper << std::endl;
	libmesh_error();
      }

    // For reading in line endings
    std::string dummy;

    // Grab a reference to the element definition for this element type
    const ElementDefinition& eledef = eletypes[elem_type];

    // If the element definition was not found, the call above would have
    // created one with an uninitialized struct.  Check for that here...
    if (eledef.abaqus_zero_based_node_id_to_libmesh_node_id.size() == 0)
      {
	libMesh::err << "No Abaqus->LibMesh mapping information for ElemType " << Utility::enum_to_string(elem_type) << "!" << std::endl;
	libmesh_error();
      }

    // We will read elements until the next line begins with *, since that will be the
    // next section.
    while (_in.peek() != '*' && _in.peek() != EOF)
      {
	// Read the element ID, it is the first number on each line.  It is
	// followed by a comma, so read that also.  We will need this ID later
	// when we try to assign subdomain IDs
	unsigned abaqus_elem_id = 0;
	char c;
	_in >> abaqus_elem_id >> c;
	// libMesh::out << "Reading data for element " << abaqus_elem_id << std::endl;

	// Add an element of the appropriate type to the Mesh.
	Elem* elem = mesh.add_elem(Elem::build(elem_type).release());

	// Associate the ID returned from libmesh with the abaqus element ID
	//_libmesh_to_abaqus_elem_mapping[elem->id()] = abaqus_elem_id;
	_abaqus_to_libmesh_elem_mapping[abaqus_elem_id] = elem->id();

	// The count of the total number of IDs read for the current element.
	unsigned id_count=0;

	// The count of IDs read from the current line.  This must always be
	// less than or equatl to _max_csv_per_line!  This is initialized to
	// 1 since we have already read the elem_id on this line...
	unsigned current_line_id_count=1;

	// Read all the node IDs for this element, which may be spread across > 1 line.
	// 1.) 9999
	// 2.) 9999, 9999
	// 3.) 9999, 9999, 9999,
	while (true)
	  {
	    // Read in the value from file
	    unsigned node_id = 0;
	    _in >> node_id;

	    // Grab the node pointer from the mesh for this ID
	    Node* node = mesh.node_ptr(node_id);

	    // If node_ptr() returns NULL, it may mean we have not yet read the
	    // *Nodes section, though I assumed that always came before the *Elements section...
	    if (node == NULL)
	      {
		libMesh::err << "Error!  Mesh returned NULL Node pointer.\n";
		libMesh::err << "Either no node exists with ID " << node_id
			     << " or perhaps this input file has *Elements defined before *Nodes?" << std::endl;
		libmesh_error();
	      }

	    // Note: id_count is the zero-based abaqus node index.  We therefore map
	    // it to a libmesh local node index using the element definition map
	    unsigned libmesh_node_id =
	      eledef.abaqus_zero_based_node_id_to_libmesh_node_id[id_count];

	    // Set this node pointer within the element.
	    elem->set_node(libmesh_node_id) = node;

	    // Increment the count of IDs read for this element
	    id_count++;

	    // Increment the number of values read from the current line
	    current_line_id_count++;

	    // Check to see if the next character is a comma.  If not, this could be
	    // the last number on the line.
	    bool found_comma = (_in.peek() == ',');

	    // Read comma from stream if found.
	    if (found_comma)
	      _in >> c;

	    // If we didn't read a comma, then this is the last node
	    // ID to be read, and we need to discard the rest of this
	    // line and break out of the while!
	    if (!found_comma)
	      {
		std::getline(_in, dummy);
		break;
	      }

	    // On the other hand, if we've already read _max_csv_per_line entries from
	    // this line, it's time to reset the current_line_id_count to zero, extract
	    // and discard the rest of this line from the file, and then continue with
	    // the while loop.
	    if (current_line_id_count == _max_csv_per_line)
	      {
		current_line_id_count = 0;
		std::getline(_in, dummy);
	      }
	  } // while

	// Ensure that we read as many nodes as we were expecting to.
	if (id_count != n_nodes_per_elem)
	  {
	    libMesh::err << "Error: Needed to read "
			 << n_nodes_per_elem
			 << " nodes, but read "
			 << id_count
			 << " instead!" << std::endl;
	    libmesh_error();
	  }

	// If we are recording Elset IDs, add this element to the correct set for later processing.
	// Make sure to add it with the Abaqus ID, not the libmesh one!
	if (elset_name != "")
	  {
	    // Debugging:
	    // libMesh::out << "Adding Elem " << abaqus_elem_id << " to Elmset " << elset_name << std::endl;
	    _elemset_ids[elset_name].push_back(abaqus_elem_id);
	  }

      } // while
  } // read_elements()




  std::string Abaqus_IO::parse_label(std::string line, std::string label_name)
  {
    // Do all string comparisons in upper-case
    std::string upper_line(line), upper_label_name(label_name);
    std::transform(upper_line.begin(), upper_line.end(), upper_line.begin(), ::toupper);
    std::transform(upper_label_name.begin(), upper_label_name.end(), upper_label_name.begin(), ::toupper);

    // Get index of start of "label="
    size_t label_index = upper_line.find(upper_label_name + "=");

    if (label_index != std::string::npos)
      {
	// Location of the first comma following "label="
	size_t comma_index = upper_line.find(",", label_index);

	// Construct iterators from which to build the sub-string.
	// Note the +1 is to skip past the "=" which follows the label name
	std::string::iterator
	  beg = line.begin() + label_name.size() + 1 + label_index,
	  end = (comma_index == std::string::npos) ? line.end() : line.begin()+comma_index;

	return std::string(beg, end);
      }

    // The label index was not found, return the empty string
    return std::string("");
  } // parse_label()




  void Abaqus_IO::read_ids(std::string set_name, container_t& container)
  {
    // Grab a reference to a vector that will hold all the IDs
    std::vector<unsigned>& id_storage = container[set_name];

    // Read until the start of another section is detected, or EOF is encountered
    while (_in.peek() != '*' && _in.peek() != EOF)
      {
	unsigned id=0;
	char c;
	unsigned id_count=0;
	std::string dummy;

	// Read a single line of IDs.  There are 3 possibilities to handle:
	// 1.) 9999
	// 2.) 9999, 9999
	// 3.) 9999, 9999, 9999,
	while (true)
	  {
	    // Read in the actual value, store in array
	    _in >> id;
	    id_storage.push_back(id);

	    // Increment the count of IDs read
	    id_count++;

	    // Check to see if the next character is a comma.  If not, this could be
	    // the last number on the line.
	    bool found_comma = (_in.peek() == ',');

	    // Read comma from stream if found.
	    if (found_comma)
	      _in >> c;

	    // If we didn't read a comma or we have
	    // already read max_id_count values, read any
	    // remaining whitespace and the newline
	    // character, and then break out of this while
	    // loop
	    if (!found_comma || id_count == _max_csv_per_line)
	      {
		std::getline(_in, dummy);
		break;
	      }
	  } // while
      } // while

    // Status message
    // libMesh::out << "Read " << id_storage.size() << " ID(s) for the set " << set_name << std::endl;
  } // read_ids()




  void Abaqus_IO::read_sideset(std::string sideset_name, sideset_container_t& container)
  {
    // Grab a reference to a vector that will hold all the IDs
    std::vector<std::pair<unsigned, unsigned> >& id_storage = container[sideset_name];

    // Variables for storing values read in from file
    unsigned elem_id=0, side_id=0;
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

	// Debugging: print status
	// std::cout << "Read elem_id=" << elem_id << ", side_id=" << side_id << std::endl;

	// Store this pair of data in the vector
	id_storage.push_back( std::make_pair(elem_id, side_id) );

	// Extract remaining characters on line including newline
	std::getline(_in, dummy);
      } // while
  }




  void Abaqus_IO::assign_subdomain_ids()
  {
    // Get a reference to the mesh we are reading
    MeshBase& mesh = MeshInput<MeshBase>::mesh();

    // The number of elemsets we've found while reading
    unsigned n_elemsets = _elemset_ids.size();

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
      // The elemset_id counter assigns a logical numbering to the _elemset_ids keys
      container_t::iterator it=_elemset_ids.begin();
      for (unsigned elemset_id=0; it != _elemset_ids.end(); ++it, ++elemset_id)
	{
	  // Grab a reference to the vector of IDs
	  std::vector<unsigned>& id_vector = (*it).second;

	  // Loop over this vector
	  for (unsigned i=0; i<id_vector.size(); ++i)
	    {
	      // Map the id_vector[i]'th element ID (Abaqus numbering) to LibMesh numbering
	      unsigned libmesh_elem_id = _abaqus_to_libmesh_elem_mapping[ id_vector[i] ];

	      // Get pointer to that element
	      Elem* elem = mesh.elem(libmesh_elem_id);

	      if (elem == NULL)
		{
		  libMesh::err << "Mesh returned NULL pointer for Elem " << libmesh_elem_id << std::endl;
		  libmesh_error();
		}

	      // Compute the proper subdomain ID, based on the formula in the
	      // documentation for this function.
	      subdomain_id_type computed_id = elemset_id + (elem_types_map[elem->type()] * n_elemsets);

	      // Assign this ID to the element in question
	      elem->subdomain_id() = computed_id;
	    }
	}
    }
  } // assign_subdomain_ids()




  void Abaqus_IO::assign_boundary_node_ids()
  {
    // Get a reference to the mesh we are reading
    MeshBase& mesh = MeshInput<MeshBase>::mesh();

    // Iterate over the container of nodesets
    container_t::iterator it=_nodeset_ids.begin();
    for (unsigned current_id=0; it != _nodeset_ids.end(); ++it, ++current_id)
      {
	libMesh::out << "Assigning node boundary ID " << current_id << " to nodeset '"
		     << (*it).first
		     << "'." << std::endl;

	// Get a reference to the current vector of nodeset ID values
	std::vector<unsigned>& nodeset_ids = (*it).second;

	for (unsigned i=0; i<nodeset_ids.size(); ++i)
	  {
	    // Get node pointer from the mesh
	    Node* node = mesh.node_ptr(nodeset_ids[i]);

	    if (node == NULL)
	      {
		libMesh::err << "Error! Mesh returned NULL node pointer!" << std::endl;
		libmesh_error();
	      }

	    // Add this node with the current_id (which is determined by the
	    // alphabetical ordering of the map)
	    mesh.boundary_info->add_node(node, current_id);
	  }
      }

  } // assign_boundary_node_ids()




  void Abaqus_IO::assign_sideset_ids()
  {
    // Get a reference to the mesh we are reading
    MeshBase& mesh = MeshInput<MeshBase>::mesh();

    // initialize the eletypes map (eletypes is a file-global variable)
    init_eletypes();

    // Iterate over the container of sidesets
    sideset_container_t::iterator it=_sideset_ids.begin();
    for (unsigned current_id=0; it != _sideset_ids.end(); ++it, ++current_id)
      {
	libMesh::out << "Assigning sideset ID " << current_id << " to sideset '"
		     << (*it).first
		     << "'." << std::endl;

	// Get a reference to the current vector of nodeset ID values
	std::vector<std::pair<unsigned,unsigned> >& sideset_ids = (*it).second;

	for (unsigned i=0; i<sideset_ids.size(); ++i)
	  {
	    // sideset_ids is a vector of pairs (elem id, side id).  Pull them out
	    // now to make the code below more readable.
	    unsigned abaqus_elem_id = sideset_ids[i].first;
	    unsigned abaqus_side_number = sideset_ids[i].second;

	    // Map the Abaqus element ID to LibMesh numbering
	    unsigned libmesh_elem_id = _abaqus_to_libmesh_elem_mapping[ abaqus_elem_id ];

	    // Get pointer to that element
	    Elem* elem = mesh.elem(libmesh_elem_id);

	    // Check that the pointer returned from the Mesh is non-NULL
	    if (elem == NULL)
	      {
		libMesh::err << "Mesh returned NULL pointer for Elem " << libmesh_elem_id << std::endl;
		libmesh_error();
	      }

	    // Grab a reference to the element definition for this element type
	    const ElementDefinition& eledef = eletypes[elem->type()];

	    // If the element definition was not found, the call above would have
	    // created one with an uninitialized struct.  Check for that here...
	    if (eledef.abaqus_zero_based_side_id_to_libmesh_side_id.size() == 0)
	      {
		libMesh::err << "No Abaqus->LibMesh mapping information for ElemType " << Utility::enum_to_string(elem->type()) << "!" << std::endl;
		libmesh_error();
	      }

	    // Add this node with the current_id (which is determined by the
	    // alphabetical ordering of the map).  Side numbers in Abaqus are 1-based,
	    // so we subtract 1 here before passing the abaqus side number to the
	    // mapping array
	    mesh.boundary_info->add_side(elem,
					 eledef.abaqus_zero_based_side_id_to_libmesh_side_id[abaqus_side_number-1],
					 current_id);
	  }
      }
  } // assign_sideset_ids()



  void Abaqus_IO::process_and_discard_comments()
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

		// Debugging:
		// libMesh::out << "Read comment line: " << dummy << std::endl;
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

  } // process_and_discard_comments()


} // namespace
