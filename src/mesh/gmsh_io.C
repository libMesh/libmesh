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
#include <set>
#include <cstring> // std::memcpy
#include <numeric>

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/mesh_base.h"
#include LIBMESH_INCLUDE_UNORDERED_MAP

namespace libMesh
{

// Initialize the static data member
GmshIO::ElementMaps GmshIO::_element_maps = GmshIO::build_element_maps();



// Definition of the static function which constructs the ElementMaps object.
GmshIO::ElementMaps GmshIO::build_element_maps()
{
  // Object to be filled up
  ElementMaps em;

  // POINT (import only)
  em.in.insert(std::make_pair(15, ElementDefinition(NODEELEM, 15, 0, 1)));

  // Add elements with trivial node mappings
  em.add_def(ElementDefinition(EDGE2, 1, 1, 2));
  em.add_def(ElementDefinition(EDGE3, 8, 1, 3));
  em.add_def(ElementDefinition(TRI3, 2, 2, 3));
  em.add_def(ElementDefinition(TRI6, 9, 2, 6));
  em.add_def(ElementDefinition(QUAD4, 3, 2, 4));
  em.add_def(ElementDefinition(QUAD8, 16, 2, 8));
  em.add_def(ElementDefinition(QUAD9, 10, 2, 9));
  em.add_def(ElementDefinition(HEX8, 5, 3, 8));
  em.add_def(ElementDefinition(TET4, 4, 3, 4));
  em.add_def(ElementDefinition(PRISM6, 6, 3, 6));
  em.add_def(ElementDefinition(PYRAMID5, 7, 3, 5));

  // Add elements with non-trivial node mappings

  // HEX20
  {
    ElementDefinition eledef(HEX20, 17, 3, 20);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,15,16,19,17,18};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  // HEX27
  {
    ElementDefinition eledef(HEX27, 12, 3, 27);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,
                                  15,16,19,17,18,20,21,24,22,23,25,26};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  // TET10
  {
    ElementDefinition eledef(TET10, 11, 3, 10);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,7,9,8};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  // PRISM15
  {
    ElementDefinition eledef(PRISM15, 18, 3, 15);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,12,14,13};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  // PRISM18
  {
    ElementDefinition eledef(PRISM18, 13, 3, 18);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,12,14,13,15,17,16};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  return em;
}



GmshIO::GmshIO (const MeshBase & mesh) :
  MeshOutput<MeshBase>(mesh),
  _binary(false),
  _write_lower_dimensional_elements(true)
{
}



GmshIO::GmshIO (MeshBase & mesh) :
  MeshInput<MeshBase>  (mesh),
  MeshOutput<MeshBase> (mesh),
  _binary (false),
  _write_lower_dimensional_elements(true)
{
}



bool & GmshIO::binary ()
{
  return _binary;
}



bool & GmshIO::write_lower_dimensional_elements ()
{
  return _write_lower_dimensional_elements;
}



void GmshIO::read (const std::string & name)
{
  std::ifstream in (name.c_str());
  this->read_mesh (in);
}



void GmshIO::read_mesh(std::istream & in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  libmesh_assert(in.good());

  // clear any data in the mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();

  // some variables
  int format=0, size=0;
  Real version = 1.0;

  // Keep track of lower-dimensional blocks which are not BCs, but
  // actually blocks of lower-dimensional elements.
  std::set<subdomain_id_type> lower_dimensional_blocks;

  // Mapping from physical id -> (physical dim, physical name) pairs.
  // These can refer to either "sidesets" or "subdomains"; we need to
  // wait until the Mesh has been read to know which is which.  Note
  // that we are using 'int' as the key here rather than
  // subdomain_id_type or boundary_id_type, since at this point, it
  // could be either.
  typedef std::pair<unsigned, std::string> GmshPhysical;
  std::map<int, GmshPhysical> gmsh_physicals;

  // map to hold the node numbers for translation
  // note the the nodes can be non-consecutive
  std::map<unsigned int, unsigned int> nodetrans;

  // For reading the file line by line
  std::string s;

  while (true)
    {
      // Try to read something.  This may set EOF!
      std::getline(in, s);

      if (in)
        {
          // Process s...

          if (s.find("$MeshFormat") == static_cast<std::string::size_type>(0))
            {
              in >> version >> format >> size;
              if ((version != 2.0) && (version != 2.1) && (version != 2.2))
                {
                  // Some notes on gmsh mesh versions:
                  //
                  // Mesh version 2.0 goes back as far as I know.  It's not explicitly
                  // mentioned here: http://www.geuz.org/gmsh/doc/VERSIONS.txt
                  //
                  // As of gmsh-2.4.0:
                  // bumped mesh version format to 2.1 (small change in the $PhysicalNames
                  // section, where the group dimension is now required);
                  // [Since we don't even parse the PhysicalNames section at the time
                  //  of this writing, I don't think this change affects us.]
                  //
                  // Mesh version 2.2 tested by Manav Bhatia; no other
                  // libMesh code changes were required for support
                  libmesh_error_msg("Error: Unknown msh file version " << version);
                }

              if (format)
                libmesh_error_msg("Error: Unknown data format for mesh in Gmsh reader.");
            }

          // Read and process the "PhysicalNames" section.
          else if (s.find("$PhysicalNames") == static_cast<std::string::size_type>(0))
            {
              // The lines in the PhysicalNames section should look like the following:
              // 2 1 "frac" lower_dimensional_block
              // 2 3 "top"
              // 2 4 "bottom"
              // 3 2 "volume"

              // Read in the number of physical groups to expect in the file.
              unsigned int num_physical_groups = 0;
              in >> num_physical_groups;

              // Read rest of line including newline character.
              std::getline(in, s);

              for (unsigned int i=0; i<num_physical_groups; ++i)
                {
                  // Read an entire line of the PhysicalNames section.
                  std::getline(in, s);

                  // Use an istringstream to extract the physical
                  // dimension, physical id, and physical name from
                  // this line.
                  std::istringstream s_stream(s);
                  unsigned phys_dim;
                  int phys_id;
                  std::string phys_name;
                  s_stream >> phys_dim >> phys_id >> phys_name;

                  // Not sure if this is true for all Gmsh files, but
                  // my test file has quotes around the phys_name
                  // string.  So let's erase any quotes now...
                  phys_name.erase(std::remove(phys_name.begin(), phys_name.end(), '"'), phys_name.end());

                  // Record this ID for later assignment of subdomain/sideset names.
                  gmsh_physicals[phys_id] = std::make_pair(phys_dim, phys_name);

                  // If 's' also contains the libmesh-specific string
                  // "lower_dimensional_block", add this block ID to
                  // the list of blocks which are not boundary
                  // conditions.
                  if (s.find("lower_dimensional_block") != std::string::npos)
                    {
                      lower_dimensional_blocks.insert(cast_int<subdomain_id_type>(phys_id));

                      // The user has explicitly told us that this
                      // block is a subdomain, so set that association
                      // in the Mesh.
                      mesh.subdomain_name(cast_int<subdomain_id_type>(phys_id)) = phys_name;
                    }
                }
            }

          // read the node block
          else if (s.find("$NOD") == static_cast<std::string::size_type>(0) ||
                   s.find("$NOE") == static_cast<std::string::size_type>(0) ||
                   s.find("$Nodes") == static_cast<std::string::size_type>(0))
            {
              unsigned int num_nodes = 0;
              in >> num_nodes;
              mesh.reserve_nodes (num_nodes);

              // read in the nodal coordinates and form points.
              Real x, y, z;
              unsigned int id;

              // add the nodal coordinates to the mesh
              for (unsigned int i=0; i<num_nodes; ++i)
                {
                  in >> id >> x >> y >> z;
                  mesh.add_point (Point(x, y, z), i);
                  nodetrans[id] = i;
                }

              // read the $ENDNOD delimiter
              std::getline(in, s);
            }


          // Read the element block
          else if (s.find("$ELM") == static_cast<std::string::size_type>(0) ||
                   s.find("$Elements") == static_cast<std::string::size_type>(0))
            {
              // For reading the number of elements and the node ids from the stream
              unsigned int
                num_elem = 0,
                node_id = 0;

              // read how many elements are there, and reserve space in the mesh
              in >> num_elem;
              mesh.reserve_elem (num_elem);

              // As of version 2.2, the format for each element line is:
              // elm-number elm-type number-of-tags < tag > ... node-number-list
              // From the Gmsh docs:
              // * the first tag is the number of the
              //   physical entity to which the element belongs
              // * the second is the number of the elementary geometrical
              //   entity to which the element belongs
              // * the third is the number of mesh partitions to which the element
              //   belongs
              // * The rest of the tags are the partition ids (negative
              //   partition ids indicate ghost cells). A zero tag is
              //   equivalent to no tag. Gmsh and most codes using the
              //   MSH 2 format require at least the first two tags
              //   (physical and elementary tags).

              // Keep track of all the element dimensions seen
              std::vector<unsigned> elem_dimensions_seen(3);

              // read the elements
              for (unsigned int iel=0; iel<num_elem; ++iel)
                {
                  unsigned int
                    id, type,
                    physical=1, elementary=1,
                    nnodes=0, ntags;

                  // Note: tag has to be an int because it could be negative,
                  // see above.
                  int tag;

                  if (version <= 1.0)
                    in >> id >> type >> physical >> elementary >> nnodes;

                  else
                    {
                      in >> id >> type >> ntags;

                      if (ntags > 2)
                        libmesh_do_once(libMesh::err << "Warning, ntags=" << ntags << ", but we currently only support reading 2 flags." << std::endl;);

                      for (unsigned int j = 0; j < ntags; j++)
                        {
                          in >> tag;
                          if (j == 0)
                            physical = tag;
                          else if (j == 1)
                            elementary = tag;
                        }
                    }

                  // Consult the import element table to determine which element to build
                  std::map<unsigned int, GmshIO::ElementDefinition>::iterator eletypes_it = _element_maps.in.find(type);

                  // Make sure we actually found something
                  if (eletypes_it == _element_maps.in.end())
                    libmesh_error_msg("Element type " << type << " not found!");

                  // Get a reference to the ElementDefinition
                  const GmshIO::ElementDefinition & eletype = eletypes_it->second;

                  // If we read nnodes, make sure it matches the number in eletype.nnodes
                  if (nnodes != 0 && nnodes != eletype.nnodes)
                    libmesh_error_msg("nnodes = " << nnodes << " and eletype.nnodes = " << eletype.nnodes << " do not match.");

                  // Assign the value from the eletype object.
                  nnodes = eletype.nnodes;

                  // Don't add 0-dimensional "point" elements to the
                  // Mesh.  They should *always* be treated as boundary
                  // "nodeset" data.
                  if (eletype.dim > 0)
                    {
                      // Record this element dimension as being "seen".
                      // We will treat all elements with dimension <
                      // max(dimension) as specifying boundary conditions,
                      // but we won't know what max_elem_dimension_seen is
                      // until we read the entire file.
                      elem_dimensions_seen[eletype.dim-1] = 1;

                      // Add the element to the mesh
                      {
                        Elem * elem = Elem::build(eletype.type).release();
                        elem->set_id(iel);
                        elem = mesh.add_elem(elem);

                        // Make sure that the libmesh element we added has nnodes nodes.
                        if (elem->n_nodes() != nnodes)
                          libmesh_error_msg("Number of nodes for element " \
                                            << id \
                                            << " of type " << eletype.type \
                                            << " (Gmsh type " << type \
                                            << ") does not match Libmesh definition. " \
                                            << "I expected " << elem->n_nodes() \
                                            << " nodes, but got " << nnodes);

                        // Add node pointers to the elements.
                        // If there is a node translation table, use it.
                        if (eletype.nodes.size() > 0)
                          for (unsigned int i=0; i<nnodes; i++)
                            {
                              in >> node_id;
                              elem->set_node(eletype.nodes[i]) = mesh.node_ptr(nodetrans[node_id]);
                            }
                        else
                          {
                            for (unsigned int i=0; i<nnodes; i++)
                              {
                                in >> node_id;
                                elem->set_node(i) = mesh.node_ptr(nodetrans[node_id]);
                              }
                          }

                        // Finally, set the subdomain ID to physical.  If this is a lower-dimension element, this ID will
                        // eventually go into the Mesh's BoundaryInfo object.
                        elem->subdomain_id() = static_cast<subdomain_id_type>(physical);
                      }
                    }

                  // Handle 0-dimensional elements (points) by adding
                  // them to the BoundaryInfo object with
                  // boundary_id == physical.
                  else
                    {
                      // This seems like it should always be the same
                      // number as the 'id' we already read in on this
                      // line.  At least it was in the example gmsh
                      // file I had...
                      in >> node_id;
                      mesh.get_boundary_info().add_node
                        (nodetrans[node_id],
                         static_cast<boundary_id_type>(physical));
                    }
                } // element loop

              // read the $ENDELM delimiter
              std::getline(in, s);

              // Record the max and min element dimension seen while reading the file.
              unsigned char
                max_elem_dimension_seen=1,
                min_elem_dimension_seen=3;

              for (std::size_t i=0; i<elem_dimensions_seen.size(); ++i)
                if (elem_dimensions_seen[i])
                  {
                    // Debugging
                    // libMesh::out << "Seen elements of dimension " << i+1 << std::endl;
                    max_elem_dimension_seen =
                      std::max(max_elem_dimension_seen, cast_int<unsigned char>(i+1));
                    min_elem_dimension_seen =
                      std::min(min_elem_dimension_seen, cast_int<unsigned char>(i+1));
                  }

              // Debugging:
              // libMesh::out << "max_elem_dimension_seen=" << max_elem_dimension_seen << std::endl;
              // libMesh::out << "min_elem_dimension_seen=" << min_elem_dimension_seen << std::endl;

              // If the difference between the max and min element dimension seen is larger than
              // 1, (e.g. the file has 1D and 3D elements only) we don't handle this case.
              if (max_elem_dimension_seen - min_elem_dimension_seen > 1)
                libmesh_error_msg("Cannot handle meshes with dimension mismatch greater than 1.");

              // How many different element dimensions did we see while reading from file?
              unsigned n_dims_seen = std::accumulate(elem_dimensions_seen.begin(),
                                                     elem_dimensions_seen.end(),
                                                     static_cast<unsigned>(0),
                                                     std::plus<unsigned>());

              // Have not yet tested a case where 1, 2, and 3D elements are all in the same Mesh,
              // though it should theoretically be possible to handle.
              if (n_dims_seen == 3)
                libmesh_error_msg("Reading meshes with 1, 2, and 3D elements not currently supported.");

              // Set mesh_dimension based on the largest element dimension seen.
              mesh.set_mesh_dimension(max_elem_dimension_seen);

              // Now that we know the maximum element dimension seen,
              // we know whether the physical names are subdomain
              // names or sideset names.
              {
                std::map<int, GmshPhysical>::iterator it = gmsh_physicals.begin();
                for (; it != gmsh_physicals.end(); ++it)
                  {
                    // Extract data
                    int phys_id = it->first;
                    unsigned phys_dim = it->second.first;
                    std::string phys_name = it->second.second;

                    // If the physical's dimension matches the largest
                    // dimension we've seen, it's a subdomain name.
                    if (phys_dim == max_elem_dimension_seen)
                      mesh.subdomain_name(cast_int<subdomain_id_type>(phys_id)) = phys_name;

                    // Otherwise, if it's not a lower-dimensional
                    // block, it's a sideset name.
                    else if (phys_dim < max_elem_dimension_seen &&
                             !lower_dimensional_blocks.count(cast_int<boundary_id_type>(phys_id)))
                      mesh.boundary_info->sideset_name(cast_int<boundary_id_type>(phys_id)) = phys_name;
                  }
              }

              if (n_dims_seen > 1)
                {
                  // Store lower-dimensional elements in a map sorted
                  // by Elem::key().  Bob Jenkins' hash functions are
                  // very good, but it's not possible for them to be
                  // perfect... so we use a multimap.
                  typedef LIBMESH_BEST_UNORDERED_MULTIMAP<dof_id_type, Elem *> provide_container_t;
                  provide_container_t provide_bcs;

                  // 1st loop over active elements - get info about lower-dimensional elements.
                  {
                    MeshBase::element_iterator       it  = mesh.active_elements_begin();
                    const MeshBase::element_iterator end = mesh.active_elements_end();
                    for ( ; it != end; ++it)
                      {
                        Elem * elem = *it;

                        if (elem->dim() < max_elem_dimension_seen &&
                            !lower_dimensional_blocks.count(elem->subdomain_id()))
                          {
                            // To be consistent with the previous
                            // GmshIO behavior, add all the
                            // lower-dimensional elements' nodes to
                            // the Mesh's BoundaryInfo object with the
                            // lower-dimensional element's subdomain
                            // ID.
                            for (unsigned n=0; n<elem->n_nodes(); n++)
                              mesh.get_boundary_info().add_node(elem->node_id(n),
                                                                elem->subdomain_id());

                            // Store this elem in a quickly-searchable
                            // container to use it to assign boundary
                            // conditions later.
                            provide_bcs.insert(std::make_pair(elem->key(), elem));
                          }
                      }
                  }

                  // 2nd loop over active elements - use lower dimensional element data to set BCs for higher dimensional elements
                  {
                    MeshBase::element_iterator       it  = mesh.active_elements_begin();
                    const MeshBase::element_iterator end = mesh.active_elements_end();

                    for ( ; it != end; ++it)
                      {
                        Elem * elem = *it;

                        if (elem->dim() == max_elem_dimension_seen)
                          {
                            // This is a max-dimension element that
                            // may require BCs.  For each of its
                            // sides, including internal sides, we'll
                            // see if a lower-dimensional element
                            // provides boundary information for it.
                            // Note that we have not yet called
                            // find_neighbors(), so we can't use
                            // elem->neighbor(sn) in this algorithm...
                            for (unsigned short sn=0; sn<elem->n_sides(); sn++)
                              {
                                // Look for the current side in the provide_bcs multimap.
                                std::pair<provide_container_t::iterator,
                                          provide_container_t::iterator>
                                  rng = provide_bcs.equal_range(elem->key(sn));

                                for (provide_container_t::iterator iter = rng.first;
                                     iter != rng.second; ++iter)
                                  {
                                    // Construct the side for hash verification.
                                    UniquePtr<Elem> side (elem->build_side_ptr(sn));

                                    // Construct the lower-dimensional element to compare to the side.
                                    Elem * lower_dim_elem = iter->second;

                                    // This was a hash, so it might not be perfect.  Let's verify...
                                    if (*lower_dim_elem == *side)
                                      {
                                        // Add the lower-dimensional
                                        // element's subdomain_id as a
                                        // boundary_id for the
                                        // higher-dimensional element.
                                        boundary_id_type bid = cast_int<boundary_id_type>(lower_dim_elem->subdomain_id());
                                        mesh.get_boundary_info().add_side(elem, sn, bid);

                                        // We only allow one match, so break out of for loop.
                                        break;
                                      }
                                  }
                              }
                          }
                      }
                  } // end 2nd loop over active elements

                  // 3rd loop over active elements - Remove the lower-dimensional elements
                  {
                    MeshBase::element_iterator       it  = mesh.active_elements_begin();
                    const MeshBase::element_iterator end = mesh.active_elements_end();
                    for ( ; it != end; ++it)
                      {
                        Elem * elem = *it;

                        if (elem->dim() < max_elem_dimension_seen &&
                            !lower_dimensional_blocks.count(elem->subdomain_id()))
                          mesh.delete_elem(elem);
                      }
                  } // end 3rd loop over active elements
                } // end if (n_dims_seen > 1)
            } // if $ELM

          continue;
        } // if (in)


      // If !in, check to see if EOF was set.  If so, break out
      // of while loop.
      if (in.eof())
        break;

      // If !in and !in.eof(), stream is in a bad state!
      libmesh_error_msg("Stream is bad! Perhaps the file does not exist?");

    } // while true
}



void GmshIO::write (const std::string & name)
{
  if (MeshOutput<MeshBase>::mesh().processor_id() == 0)
    {
      // Open the output file stream
      std::ofstream out_stream (name.c_str());

      // Make sure it opened correctly
      if (!out_stream.good())
        libmesh_file_error(name.c_str());

      this->write_mesh (out_stream);
    }
}



void GmshIO::write_nodal_data (const std::string & fname,
                               const std::vector<Number> & soln,
                               const std::vector<std::string> & names)
{
  LOG_SCOPE("write_nodal_data()", "GmshIO");

  if (MeshOutput<MeshBase>::mesh().processor_id() == 0)
    this->write_post  (fname, &soln, &names);
}



void GmshIO::write_mesh (std::ostream & out_stream)
{
  // Be sure that the stream is valid.
  libmesh_assert (out_stream.good());

  // Get a const reference to the mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // If requested, write out lower-dimensional elements for
  // element-side-based boundary conditions.
  unsigned int n_boundary_faces = 0;
  if (this->write_lower_dimensional_elements())
    n_boundary_faces = mesh.boundary_info->n_boundary_conds();

  // Note: we are using version 2.0 of the gmsh output format.

  // Write the file header.
  out_stream << "$MeshFormat\n";
  out_stream << "2.0 0 " << sizeof(Real) << '\n';
  out_stream << "$EndMeshFormat\n";

  // write the nodes in (n x y z) format
  out_stream << "$Nodes\n";
  out_stream << mesh.n_nodes() << '\n';

  for (unsigned int v=0; v<mesh.n_nodes(); v++)
    out_stream << mesh.node_ref(v).id()+1 << " "
               << mesh.node_ref(v)(0) << " "
               << mesh.node_ref(v)(1) << " "
               << mesh.node_ref(v)(2) << '\n';
  out_stream << "$EndNodes\n";

  {
    // write the connectivity
    out_stream << "$Elements\n";
    out_stream << mesh.n_active_elem() + n_boundary_faces << '\n';

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    // loop over the elements
    for ( ; it != end; ++it)
      {
        const Elem * elem = *it;

        // Make sure we have a valid entry for
        // the current element type.
        libmesh_assert (_element_maps.out.count(elem->type()));

        // consult the export element table
        std::map<ElemType, ElementDefinition>::iterator def_it =
          _element_maps.out.find(elem->type());

        // Assert that we found it
        if (def_it == _element_maps.out.end())
          libmesh_error_msg("Element type " << elem->type() << " not found in _element_maps.");

        // Get a reference to the ElementDefinition object
        const ElementDefinition & eletype = def_it->second;

        // The element mapper better not require any more nodes
        // than are present in the current element!
        libmesh_assert_less_equal (eletype.nodes.size(), elem->n_nodes());

        // elements ids are 1 based in Gmsh
        out_stream << elem->id()+1 << " ";
        // element type
        out_stream << eletype.gmsh_type;

        // write the number of tags (3) and their values:
        // 1 (physical entity)
        // 2 (geometric entity)
        // 3 (partition entity)
        out_stream << " 3 "
                   << static_cast<unsigned int>(elem->subdomain_id())
                   << " 0 "
                   << elem->processor_id()+1
                   << " ";

        // if there is a node translation table, use it
        if (eletype.nodes.size() > 0)
          for (unsigned int i=0; i < elem->n_nodes(); i++)
            out_stream << elem->node_id(eletype.nodes[i])+1 << " "; // gmsh is 1-based
        // otherwise keep the same node order
        else
          for (unsigned int i=0; i < elem->n_nodes(); i++)
            out_stream << elem->node_id(i)+1 << " ";                  // gmsh is 1-based
        out_stream << "\n";
      } // element loop
  }

  {
    // A counter for writing surface elements to the Gmsh file
    // sequentially.  We start numbering them with a number strictly
    // larger than the largest element ID in the mesh.  Note: the
    // MeshBase docs say "greater than or equal to" the maximum
    // element id in the mesh, so technically we might need a +1 here,
    // but all of the implementations return an ID strictly greater
    // than the largest element ID in the Mesh.
    unsigned int e_id = mesh.max_elem_id();

    // loop over the elements, writing out boundary faces
    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    if (n_boundary_faces)
      {
        // Construct the list of boundary sides
        std::vector<dof_id_type> element_id_list;
        std::vector<unsigned short int> side_list;
        std::vector<boundary_id_type> bc_id_list;

        mesh.boundary_info->build_side_list(element_id_list, side_list, bc_id_list);

        // Loop over these lists, writing data to the file.
        for (std::size_t idx=0; idx<element_id_list.size(); ++idx)
          {
            const Elem & elem = mesh.elem_ref(element_id_list[idx]);

            UniquePtr<const Elem> side = elem.build_side_ptr(side_list[idx]);

            // Map from libmesh elem type to gmsh elem type.
            std::map<ElemType, ElementDefinition>::iterator def_it =
              _element_maps.out.find(side->type());

            // If we didn't find it, that's an error
            if (def_it == _element_maps.out.end())
              libmesh_error_msg("Element type " << side->type() << " not found in _element_maps.");

            // consult the export element table
            const GmshIO::ElementDefinition & eletype = def_it->second;

            // The element mapper better not require any more nodes
            // than are present in the current element!
            libmesh_assert_less_equal (eletype.nodes.size(), side->n_nodes());

            // elements ids are 1-based in Gmsh
            out_stream << e_id+1 << " ";

            // element type
            out_stream << eletype.gmsh_type;

            // write the number of tags:
            // 1 (physical entity)
            // 2 (geometric entity)
            // 3 (partition entity)
            out_stream << " 3 "
                       << bc_id_list[idx]
                       << " 0 "
                       << elem.processor_id()+1
                       << " ";

            // if there is a node translation table, use it
            if (eletype.nodes.size() > 0)
              for (unsigned int i=0; i < side->n_nodes(); i++)
                out_stream << side->node_id(eletype.nodes[i])+1 << " "; // gmsh is 1-based

            // otherwise keep the same node order
            else
              for (unsigned int i=0; i < side->n_nodes(); i++)
                out_stream << side->node_id(i)+1 << " ";                // gmsh is 1-based

            // Go to the next line
            out_stream << "\n";

            // increment this index too...
            ++e_id;
          }
      }

    out_stream << "$EndElements\n";
  }
}



void GmshIO::write_post (const std::string & fname,
                         const std::vector<Number> * v,
                         const std::vector<std::string> * solution_names)
{

  // Should only do this on processor 0!
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  // Create an output stream
  std::ofstream out_stream(fname.c_str());

  // Make sure it opened correctly
  if (!out_stream.good())
    libmesh_file_error(fname.c_str());

  // create a character buffer
  char buf[80];

  // Get a constant reference to the mesh.
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  //  write the data
  if ((solution_names != libmesh_nullptr) && (v != libmesh_nullptr))
    {
      const unsigned int n_vars =
        cast_int<unsigned int>(solution_names->size());

      if (!(v->size() == mesh.n_nodes()*n_vars))
        libMesh::err << "ERROR: v->size()=" << v->size()
                     << ", mesh.n_nodes()=" << mesh.n_nodes()
                     << ", n_vars=" << n_vars
                     << ", mesh.n_nodes()*n_vars=" << mesh.n_nodes()*n_vars
                     << "\n";

      libmesh_assert_equal_to (v->size(), mesh.n_nodes()*n_vars);

      // write the header
      out_stream << "$PostFormat\n";
      if (this->binary())
        out_stream << "1.2 1 " << sizeof(double) << "\n";
      else
        out_stream << "1.2 0 " << sizeof(double) << "\n";
      out_stream << "$EndPostFormat\n";

      // Loop over the elements to see how much of each type there are
      unsigned int n_points=0, n_lines=0, n_triangles=0, n_quadrangles=0,
        n_tetrahedra=0, n_hexahedra=0, n_prisms=0, n_pyramids=0;
      unsigned int n_scalar=0, n_vector=0, n_tensor=0;
      unsigned int nb_text2d=0, nb_text2d_chars=0, nb_text3d=0, nb_text3d_chars=0;

      {
        MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
        const MeshBase::const_element_iterator end = mesh.active_elements_end();


        for ( ; it != end; ++it)
          {
            const ElemType elemtype = (*it)->type();

            switch (elemtype)
              {
              case EDGE2:
              case EDGE3:
              case EDGE4:
                {
                  n_lines += 1;
                  break;
                }
              case TRI3:
              case TRI6:
                {
                  n_triangles += 1;
                  break;
                }
              case QUAD4:
              case QUAD8:
              case QUAD9:
                {
                  n_quadrangles += 1;
                  break;
                }
              case TET4:
              case TET10:
                {
                  n_tetrahedra += 1;
                  break;
                }
              case HEX8:
              case HEX20:
              case HEX27:
                {
                  n_hexahedra += 1;
                  break;
                }
              case PRISM6:
              case PRISM15:
              case PRISM18:
                {
                  n_prisms += 1;
                  break;
                }
              case PYRAMID5:
                {
                  n_pyramids += 1;
                  break;
                }
              default:
                libmesh_error_msg("ERROR: Nonexistent element type " << (*it)->type());
              }
          }
      }

      // create a view for each variable
      for (unsigned int ivar=0; ivar < n_vars; ivar++)
        {
          std::string varname = (*solution_names)[ivar];

          // at the moment, we just write out scalar quantities
          // later this should be made configurable through
          // options to the writer class
          n_scalar = 1;

          // write the variable as a view, and the number of time steps
          out_stream << "$View\n" << varname << " " << 1 << "\n";

          // write how many of each geometry type are written
          out_stream << n_points * n_scalar << " "
                     << n_points * n_vector << " "
                     << n_points * n_tensor << " "
                     << n_lines * n_scalar << " "
                     << n_lines * n_vector << " "
                     << n_lines * n_tensor << " "
                     << n_triangles * n_scalar << " "
                     << n_triangles * n_vector << " "
                     << n_triangles * n_tensor << " "
                     << n_quadrangles * n_scalar << " "
                     << n_quadrangles * n_vector << " "
                     << n_quadrangles * n_tensor << " "
                     << n_tetrahedra * n_scalar << " "
                     << n_tetrahedra * n_vector << " "
                     << n_tetrahedra * n_tensor << " "
                     << n_hexahedra * n_scalar << " "
                     << n_hexahedra * n_vector << " "
                     << n_hexahedra * n_tensor << " "
                     << n_prisms * n_scalar << " "
                     << n_prisms * n_vector << " "
                     << n_prisms * n_tensor << " "
                     << n_pyramids * n_scalar << " "
                     << n_pyramids * n_vector << " "
                     << n_pyramids * n_tensor << " "
                     << nb_text2d << " "
                     << nb_text2d_chars << " "
                     << nb_text3d << " "
                     << nb_text3d_chars << "\n";

          // if binary, write a marker to identify the endianness of the file
          if (this->binary())
            {
              const int one = 1;
              std::memcpy(buf, &one, sizeof(int));
              out_stream.write(buf, sizeof(int));
            }

          // the time steps (there is just 1 at the moment)
          if (this->binary())
            {
              double one = 1;
              std::memcpy(buf, &one, sizeof(double));
              out_stream.write(buf, sizeof(double));
            }
          else
            out_stream << "1\n";

          // Loop over the elements and write out the data
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end();

          for ( ; it != end; ++it)
            {
              const Elem * elem = *it;

              // this is quite crappy, but I did not invent that file format!
              for (unsigned int d=0; d<3; d++)  // loop over the dimensions
                {
                  for (unsigned int n=0; n < elem->n_vertices(); n++)   // loop over vertices
                    {
                      const Point & vertex = elem->point(n);
                      if (this->binary())
                        {
                          double tmp = vertex(d);
                          std::memcpy(buf, &tmp, sizeof(double));
                          out_stream.write(reinterpret_cast<char *>(buf), sizeof(double));
                        }
                      else
                        out_stream << vertex(d) << " ";
                    }
                  if (!this->binary())
                    out_stream << "\n";
                }

              // now finally write out the data
              for (unsigned int i=0; i < elem->n_vertices(); i++)   // loop over vertices
                if (this->binary())
                  {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
                    libMesh::out << "WARNING: Gmsh::write_post does not fully support "
                                 << "complex numbers. Will only write the real part of "
                                 << "variable " << varname << std::endl;
#endif
                    double tmp = libmesh_real((*v)[elem->node_id(i)*n_vars + ivar]);
                    std::memcpy(buf, &tmp, sizeof(double));
                    out_stream.write(reinterpret_cast<char *>(buf), sizeof(double));
                  }
                else
                  {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
                    libMesh::out << "WARNING: Gmsh::write_post does not fully support "
                                 << "complex numbers. Will only write the real part of "
                                 << "variable " << varname << std::endl;
#endif
                    out_stream << libmesh_real((*v)[elem->node_id(i)*n_vars + ivar]) << "\n";
                  }
            }
          if (this->binary())
            out_stream << "\n";
          out_stream << "$EndView\n";

        } // end variable loop (writing the views)
    }
}

} // namespace libMesh
