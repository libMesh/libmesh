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


#ifndef LIBMESH_ABAQUS_IO_H
#define LIBMESH_ABAQUS_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"

// C++ includes
#include <fstream>
#include <set>

namespace libMesh
{

/**
 * The AbaqusIO class is a preliminary implementation for reading
 * Abaqus mesh files in ASCII format.
 *
 * \author John W. Peterson
 * \date 2011
 */
class AbaqusIO : public MeshInput<MeshBase>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   */
  explicit
  AbaqusIO (MeshBase & mesh);

  /**
   * Destructor.
   */
  virtual ~AbaqusIO ();

  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string & name) libmesh_override;

  /**
   * Default false.  Abaqus files have only nodesets in them by
   * default.  Set this flag to true if you want libmesh to automatically
   * generate sidesets from Abaqus' nodesets.
   */
  bool build_sidesets_from_nodesets;

private:
  /**
   * The type of data structure used to store Node and Elemset IDs.
   */
  typedef std::map<std::string, std::vector<dof_id_type> > container_t;

  /**
   * Type of the data structure for storing the (elem ID, side) pairs
   * defining sidesets.  These come from the *Surface sections of the
   * input file.
   */
  typedef std::map<std::string, std::vector<std::pair<dof_id_type, unsigned> > > sideset_container_t;

  /**
   * This function parses a block of nodes in the Abaqus file once
   * such a block has been found.  If the *NODE section specifies an
   * NSET name, also pass that to this function.
   */
  void read_nodes(std::string nset_name);

  /**
   * This function parses a block of elements in the Abaqus file.
   * You must pass it an upper-cased version of the string declaring
   * this section, which is typically something like:
   * *ELEMENT, TYPE=CPS3
   * so that it can determine the type of elements to read.
   */
  void read_elements(std::string upper, std::string elset_name);

  /**
   * This function parses a label of the form foo=bar from a
   * comma-delimited line of the form
   * ..., foo=bar, ...
   * The input to the function in this case would be foo, the
   * output would be bar
   */
  std::string parse_label(std::string line, std::string label_name) const;

  /**
   * Returns true if the input string is a generated elset or nset,
   * false otherwise. The input string is assumed to already be in all
   * upper case. Generated nsets are assumed to have the following
   * format:
   * *Nset, nset=Set-1, generate
   */
  bool detect_generated_set(std::string upper) const;

  /**
   * This function reads all the IDs for the current node or element
   * set of the given name, storing them in the passed map using the
   * name as key.
   */
  void read_ids(std::string set_name, container_t & container);

  /**
   * This function handles "generated" nset and elset sections.  These
   * are denoted by having the comma-separated "GENERATE" keyword in
   * their definition, e.g.
   * *Nset, nset=Set-1, generate
   */
  void generate_ids(std::string set_name, container_t & container);

  /**
   * This function is called after all the elements have been
   * read and assigns element subdomain IDs.
   *
   * The IDs are simply chosen in the order in which the elset
   * labels are stored in the map (roughly alphabetically).  To make
   * this easy on people who are planning to use Exodus output,
   * we'll assign different geometric elements to different (but
   * related) subdomains, i.e. assuming there are E elemsets:
   *
   * Elemset 0, Geometric Type 0: ID 0
   * Elemset 0, Geometric Type 1: ID 0+E
   * ...
   * Elemset 0, Geometric Type N: ID 0+N*E
   * --------------------------------------
   * Elemset 1, Geometric Type 0: ID 1
   * Elemset 1, Geometric Type 1: ID 1+E
   * ...
   * Elemset 1, Geometric Type N: ID 1+N*E
   * etc.
   */
  void assign_subdomain_ids();

  /**
   * This function reads a sideset from the input file.  This is defined
   * by a "*Surface" section in the file, and then a list of element ID
   * and side IDs for the set.
   */
  void read_sideset(std::string sideset_name, sideset_container_t & container);

  /**
   * This function assigns boundary IDs to node sets based on the
   * alphabetical order in which the sets are labelled in the Abaqus
   * file.  We choose the alphabetical ordering simply because
   * Abaqus does not provide a numerical one within the file.
   */
  void assign_boundary_node_ids();

  /**
   * Called at the end of the read() function, assigns any sideset IDs
   * found when reading the file to the BoundaryInfo object.
   */
  void assign_sideset_ids();

  /**
   * Any of the various sections can start with some number of lines
   * of comments, which start with "**".  This function discards
   * any lines of comments that it finds from the stream, leaving
   * trailing data intact.
   */
  void process_and_discard_comments();

  /**
   * Returns the maximum geometric element dimension encountered while
   * reading the Mesh.  Only valid after the elements have been read
   * in and the elems_of_dimension array has been populated.
   */
  unsigned char max_elem_dimension_seen();

  /**
   * Abaqus writes nodesets and elemsets with labels.  As we read
   * them in, we'll use these maps to provide a natural ordering for
   * them.
   */
  container_t _nodeset_ids;
  container_t _elemset_ids;
  sideset_container_t _sideset_ids;

  /**
   * Stream object used to interact with the file
   */
  std::ifstream _in;

  /**
   * A set of the different geometric element types detected when reading the
   * mesh.
   */
  std::set<ElemType> _elem_types;

  /**
   * Map from libmesh element number -> abaqus element number,
   * and the converse.
   */
  // std::map<dof_id_type, dof_id_type> _libmesh_to_abaqus_elem_mapping;
  std::map<dof_id_type, dof_id_type> _abaqus_to_libmesh_elem_mapping;

  /**
   * Map from abaqus node number -> sequential, 0-based libmesh node numbering.
   * Note that in every Abaqus file I've ever seen the node numbers were 1-based,
   * sequential, and all in order, so that this map is probably overkill.
   * Nevertheless, it is the most general solution in case we come across a
   * weird Abaqus file some day.
   */
  std::map<dof_id_type, dof_id_type> _abaqus_to_libmesh_node_mapping;

  /**
   * This flag gets set to true after the first "*PART" section
   * we see.  If it is still true when we see a second PART
   * section, we will print an error message... we don't currently
   * handle input files with multiple parts.
   */
  bool _already_seen_part;
};

} // namespace

#endif // LIBMESH_ABAQUS_IO_H
