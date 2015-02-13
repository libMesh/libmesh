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


// anonymous namespace to hold local data
namespace
{
using namespace libMesh;

/**
 * Defines mapping from libMesh element types to Gmsh element types.
 */
struct elementDefinition
{
  std::string label;
  std::vector<unsigned int> nodes;
  ElemType type;
  unsigned int exptype;
  unsigned int dim;
  unsigned int nnodes;
};


// maps from a libMesh element type to the proper
// Gmsh elementDefinition.  Placing the data structure
// here in this anonymous namespace gives us the
// benefits of a global variable without the nasty
// side-effects
std::map<ElemType, elementDefinition> eletypes_exp;
std::map<unsigned int, elementDefinition> eletypes_imp;



// ------------------------------------------------------------
// helper function to initialize the eletypes map
void init_eletypes ()
{
  if (eletypes_exp.empty() && eletypes_imp.empty())
    {
      // This should happen only once.  The first time this method
      // is called the eletypes data struture will be empty, and
      // we will fill it.  Any subsequent calls will find an initialized
      // eletypes map and will do nothing.

      // set up the element definitions
      elementDefinition eledef;

      // use "swap trick" from Scott Meyer's "Effective STL" to initialize
      // eledef.nodes vector

      // POINT (only Gmsh)
      {
        eledef.type    = NODEELEM;
        eledef.exptype = 15;
        eledef.dim     = 0;
        eledef.nnodes  = 1;
        eledef.nodes.clear();

        // import only
        eletypes_imp[15] = eledef;
      }

      // EDGE2
      {
        eledef.type    = EDGE2;
        eledef.dim     = 1;
        eledef.nnodes  = 2;
        eledef.exptype = 1;
        eledef.nodes.clear();

        eletypes_exp[EDGE2] = eledef;
        eletypes_imp[1]     = eledef;
      }

      // EDGE3
      {
        eledef.type    = EDGE3;
        eledef.dim     = 1;
        eledef.nnodes  = 3;
        eledef.exptype = 8;
        eledef.nodes.clear();

        eletypes_exp[EDGE3] = eledef;
        eletypes_imp[8]     = eledef;
      }

      // TRI3
      {
        eledef.type    = TRI3;
        eledef.dim     = 2;
        eledef.nnodes  = 3;
        eledef.exptype = 2;
        eledef.nodes.clear();

        eletypes_exp[TRI3] = eledef;
        eletypes_imp[2] = eledef;
      }

      // TRI6
      {
        eledef.type    = TRI6;
        eledef.dim     = 2;
        eledef.nnodes  = 6;
        eledef.exptype = 9;
        eledef.nodes.clear();

        eletypes_exp[TRI6] = eledef;
        eletypes_imp[9]    = eledef;
      }

      // QUAD4
      {
        eledef.type    = QUAD4;
        eledef.dim     = 2;
        eledef.nnodes  = 4;
        eledef.exptype = 3;
        eledef.nodes.clear();

        eletypes_exp[QUAD4] = eledef;
        eletypes_imp[3]     = eledef;
      }

      // QUAD8
      // TODO: what should be done with this on writing?
      {
        eledef.type    = QUAD8;
        eledef.dim     = 2;
        eledef.nnodes  = 8;
        eledef.exptype = 100;
        const unsigned int nodes[] = {1,2,3,4,5,6,7,8};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[QUAD8] = eledef;
        eletypes_imp[10]    = eledef;
      }

      // QUAD9
      {
        eledef.type    = QUAD9;
        eledef.dim     = 2;
        eledef.nnodes  = 9;
        eledef.exptype = 10;
        eledef.nodes.clear();

        eletypes_exp[QUAD9] = eledef;
        eletypes_imp[10]    = eledef;
      }

      // HEX8
      {
        eledef.type    = HEX8;
        eledef.dim     = 3;
        eledef.nnodes  = 8;
        eledef.exptype = 5;
        eledef.nodes.clear();

        eletypes_exp[HEX8] = eledef;
        eletypes_imp[5]    = eledef;
      }

      // HEX20
      // TODO: what should be done with this on writing?
      {
        eledef.type    = HEX20;
        eledef.dim     = 3;
        eledef.nnodes  = 20;
        eledef.exptype = 101;
        const unsigned int nodes[] = {1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,16};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[HEX20] = eledef;
        eletypes_imp[12]    = eledef;
      }

      // HEX27
      {
        eledef.type    = HEX27;
        eledef.dim     = 3;
        eledef.nnodes  = 27;
        eledef.exptype = 12;
        const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,
                                      15,16,19,17,18,20,21,24,22,23,25,26};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[HEX27] = eledef;
        eletypes_imp[12]    = eledef;
      }

      // TET4
      {
        eledef.type    = TET4;
        eledef.dim     = 3;
        eledef.nnodes  = 4;
        eledef.exptype = 4;
        eledef.nodes.clear();

        eletypes_exp[TET4] = eledef;
        eletypes_imp[4]    = eledef;
      }

      // TET10
      {
        eledef.type    = TET10;
        eledef.dim     = 3;
        eledef.nnodes  = 10;
        eledef.exptype = 11;
        const unsigned int nodes[] = {0,1,2,3,4,5,6,7,9,8};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);
        eletypes_exp[TET10] = eledef;
        eletypes_imp[11]    = eledef;
      }

      // PRISM6
      {
        eledef.type    = PRISM6;
        eledef.dim     = 3;
        eledef.nnodes  = 6;
        eledef.exptype = 6;
        eledef.nodes.clear();

        eletypes_exp[PRISM6] = eledef;
        eletypes_imp[6]      = eledef;
      }

      // PRISM15
      // TODO: what should be done with this on writing?
      {
        eledef.type    = PRISM15;
        eledef.dim     = 3;
        eledef.nnodes  = 15;
        eledef.exptype = 103;
        eledef.nodes.clear();

        eletypes_exp[PRISM15] = eledef;
        eletypes_imp[13] = eledef;
      }

      // PRISM18
      {
        eledef.type    = PRISM18;
        eledef.dim     = 3;
        eledef.nnodes  = 18;
        eledef.exptype = 13;
        const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,
                                      12,14,13,15,17,16};
        std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes);

        eletypes_exp[PRISM18] = eledef;
        eletypes_imp[13]      = eledef;
      }

      // PYRAMID5
      {
        eledef.type    = PYRAMID5;
        eledef.dim     = 3;
        eledef.nnodes  = 5;
        eledef.exptype = 7;
        eledef.nodes.clear();

        eletypes_exp[PYRAMID5] = eledef;
        eletypes_imp[7]        = eledef;
      }
    }
}

} // end anonymous namespace


namespace libMesh
{

// ------------------------------------------------------------
// GmshIO  members

GmshIO::GmshIO (const MeshBase& mesh) :
  MeshOutput<MeshBase>(mesh),
  _binary(false)
{
}



GmshIO::GmshIO (MeshBase& mesh) :
  MeshInput<MeshBase>  (mesh),
  MeshOutput<MeshBase> (mesh),
  _binary (false)
{
}



bool & GmshIO::binary ()
{
  return _binary;
}



void GmshIO::read (const std::string& name)
{
  std::ifstream in (name.c_str());
  this->read_mesh (in);
}



void GmshIO::read_mesh(std::istream& in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  libmesh_assert(in.good());

  // initialize the map with element types
  init_eletypes();

  // clear any data in the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();

  // some variables
  int format=0, size=0;
  Real version = 1.0;

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
                  std::map<unsigned int, elementDefinition>::iterator eletypes_it = eletypes_imp.find(type);

                  // Make sure we actually found something
                  if (eletypes_it == eletypes_imp.end())
                    libmesh_error_msg("Element type " << type << " not found!");

                  // Get a reference to the elementDefinition
                  const elementDefinition& eletype = eletypes_it->second;

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
                        Elem* elem = Elem::build(eletype.type).release();
                        elem->set_id(iel);
                        elem = mesh.add_elem(elem);

                        // Make sure that the libmesh element we added has nnodes nodes.
                        if (elem->n_nodes() != nnodes)
                          libmesh_error_msg("Number of nodes for element " \
                                            << id \
                                            << " of type " << eletypes_imp[type].type \
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

              for (unsigned char i=0; i<elem_dimensions_seen.size(); ++i)
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

              if (n_dims_seen > 1)
                {
                  // map from (node ids) -> elem of lower dimensional elements that can provide boundary conditions
                  typedef std::map<std::vector<dof_id_type>, Elem*> provide_container_t;
                  provide_container_t provide_bcs;

                  // 1st loop over active elements - get info about lower-dimensional elements.
                  {
                    MeshBase::element_iterator       it  = mesh.active_elements_begin();
                    const MeshBase::element_iterator end = mesh.active_elements_end();
                    for ( ; it != end; ++it)
                      {
                        Elem* elem = *it;

                        if (elem->dim() < max_elem_dimension_seen)
                          {
                            // Debugging status
                            // libMesh::out << "Processing Elem " << elem->id() << " as a boundary element." << std::endl;

                            // To be pushed into the provide_bcs data structure
                            std::vector<dof_id_type> node_ids(elem->n_nodes());

                            // To be consistent with the previous GmshIO behavior, add all the lower-dimensional elements' nodes to
                            // the Mesh's BoundaryInfo object with the lower-dimensional element's subdomain ID.
                            for (unsigned n=0; n<elem->n_nodes(); n++)
                              {
                                mesh.get_boundary_info().add_node
                                  (elem->node(n), elem->subdomain_id());

                                // And save for our local data structure
                                node_ids[n] = elem->node(n);
                              }

                            // Sort before putting into the map
                            std::sort(node_ids.begin(), node_ids.end());
                            provide_bcs[node_ids] = elem;
                          }
                      }
                  } // end 1st loop over active elements

                  // Debugging: What did we put in the provide_bcs data structure?
                  // {
                  //   provide_container_t::iterator provide_it = provide_bcs.begin();
                  //   provide_container_t::iterator provide_end = provide_bcs.end();
                  //   for ( ; provide_it != provide_end; ++provide_it)
                  //     {
                  //       std::vector<dof_id_type> node_list = (*provide_it).first;
                  //       Elem* elem = (*provide_it).second;
                  //
                  //       libMesh::out << "Elem " << elem->id() << " provides BCs for the face: ";
                  //       for (unsigned i=0; i<node_list.size(); ++i)
                  //         libMesh::out << node_list[i] << " ";
                  //       libMesh::out << std::endl;
                  //     }
                  // }

                  // 2nd loop over active elements - use lower dimensional element data to set BCs for higher dimensional elements
                  {
                    MeshBase::element_iterator       it  = mesh.active_elements_begin();
                    const MeshBase::element_iterator end = mesh.active_elements_end();
                    for ( ; it != end; ++it)
                      {
                        Elem* elem = *it;

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

                            for (unsigned short sn=0;
                                 sn<elem->n_sides(); sn++)
                              {
                                AutoPtr<Elem> side (elem->build_side(sn));

                                // Build up a node_ids vector, which is the key
                                std::vector<dof_id_type> node_ids(side->n_nodes());
                                for (unsigned n=0; n<side->n_nodes(); n++)
                                  node_ids[n] = side->node(n);

                                // Sort the vector before using it as a key
                                std::sort(node_ids.begin(), node_ids.end());

                                // Look for this key in the provide_bcs map
                                provide_container_t::iterator iter = provide_bcs.find(node_ids);

                                if (iter != provide_bcs.end())
                                  {
                                    Elem* lower_dim_elem = (*iter).second;

                                    // libMesh::out << "Elem "
                                    //              << lower_dim_elem->id()
                                    //              << " provides BCs for side "
                                    //              << sn
                                    //              << " of Elem "
                                    //              << elem->id()
                                    //              << std::endl;

                                    // Add boundary information based on the lower-dimensional element's subdomain id.
                                    mesh.get_boundary_info().add_side(elem,
                                                                      sn,
                                                                      cast_int<boundary_id_type>(lower_dim_elem->subdomain_id()));
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
                        Elem* elem = *it;

                        if (elem->dim() < max_elem_dimension_seen)
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



void GmshIO::write (const std::string& name)
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



void GmshIO::write_nodal_data (const std::string& fname,
                               const std::vector<Number>& soln,
                               const std::vector<std::string>& names)
{
  START_LOG("write_nodal_data()", "GmshIO");

  //this->_binary = true;
  if (MeshOutput<MeshBase>::mesh().processor_id() == 0)
    this->write_post  (fname, &soln, &names);

  STOP_LOG("write_nodal_data()", "GmshIO");
}



void GmshIO::write_mesh (std::ostream& out_stream)
{
  // Be sure that the stream is valid.
  libmesh_assert (out_stream.good());

  // initialize the map with element types
  init_eletypes();

  // Get a const reference to the mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  // Note: we are using version 2.0 of the gmsh output format.

  // Write the file header.
  out_stream << "$MeshFormat\n";
  out_stream << "2.0 0 " << sizeof(Real) << '\n';
  out_stream << "$EndMeshFormat\n";

  // write the nodes in (n x y z) format
  out_stream << "$Nodes\n";
  out_stream << mesh.n_nodes() << '\n';

  for (unsigned int v=0; v<mesh.n_nodes(); v++)
    out_stream << mesh.node(v).id()+1 << " "
               << mesh.node(v)(0) << " "
               << mesh.node(v)(1) << " "
               << mesh.node(v)(2) << '\n';
  out_stream << "$EndNodes\n";

  {
    // write the connectivity
    out_stream << "$Elements\n";
    out_stream << mesh.n_active_elem() << '\n';

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    // loop over the elements
    for ( ; it != end; ++it)
      {
        const Elem* elem = *it;

        // Make sure we have a valid entry for
        // the current element type.
        libmesh_assert (eletypes_exp.count(elem->type()));

        // consult the export element table
        const elementDefinition& eletype = eletypes_exp[elem->type()];

        // The element mapper better not require any more nodes
        // than are present in the current element!
        libmesh_assert_less_equal (eletype.nodes.size(), elem->n_nodes());

        // elements ids are 1 based in Gmsh
        out_stream << elem->id()+1 << " ";

        // element type
        out_stream << eletype.exptype;

        // write the number of tags and
        // tag1 (physical entity), tag2 (geometric entity), and tag3 (partition entity)
        out_stream << " 3 "
                   << static_cast<unsigned int>(elem->subdomain_id())
                   << " 1 "
                   << (elem->processor_id()+1) << " ";

        // if there is a node translation table, use it
        if (eletype.nodes.size() > 0)
          for (unsigned int i=0; i < elem->n_nodes(); i++)
            out_stream << elem->node(eletype.nodes[i])+1 << " "; // gmsh is 1-based
        // otherwise keep the same node order
        else
          for (unsigned int i=0; i < elem->n_nodes(); i++)
            out_stream << elem->node(i)+1 << " ";                  // gmsh is 1-based
        out_stream << "\n";
      } // element loop
    out_stream << "$EndElements\n";
  }
}



void GmshIO::write_post (const std::string& fname,
                         const std::vector<Number>* v,
                         const std::vector<std::string>* solution_names)
{

  // Should only do this on processor 0!
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  // Create an output stream
  std::ofstream out_stream(fname.c_str());

  // Make sure it opened correctly
  if (!out_stream.good())
    libmesh_file_error(fname.c_str());

  // initialize the map with element types
  init_eletypes();

  // create a character buffer
  char buf[80];

  // Get a constant reference to the mesh.
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  //  write the data
  if ((solution_names != NULL) && (v != NULL))
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
              const Elem* elem = *it;

              // this is quite crappy, but I did not invent that file format!
              for (unsigned int d=0; d<3; d++)  // loop over the dimensions
                {
                  for (unsigned int n=0; n < elem->n_vertices(); n++)   // loop over vertices
                    {
                      const Point& vertex = elem->point(n);
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
                    double tmp = libmesh_real((*v)[elem->node(i)*n_vars + ivar]);
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
                    out_stream << libmesh_real((*v)[elem->node(i)*n_vars + ivar]) << "\n";
                  }
            }
          if (this->binary())
            out_stream << "\n";
          out_stream << "$EndView\n";

        } // end variable loop (writing the views)
    }
}

} // namespace libMesh
