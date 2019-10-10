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

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/dyna_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/int_range.h"

// C++ includes
#include <fstream>
#include <cstddef>
#include <array>

namespace libMesh
{

// Initialize the static data member
DynaIO::ElementMaps DynaIO::_element_maps = DynaIO::build_element_maps();



// Definition of the static function which constructs the ElementMaps object.
DynaIO::ElementMaps DynaIO::build_element_maps()
{
  // Object to be filled up
  ElementMaps em;

  // em.add_def(ElementDefinition(TRI3, 2, 2, 1)); // node mapping?
  // em.add_def(ElementDefinition(TRI6, 2, 2, 2)); // node mapping?
  // em.add_def(ElementDefinition(TET4,  2, 3, 1));   // node mapping?
  // em.add_def(ElementDefinition(TET10,  2, 3, 2));  // node mapping?
  // em.add_def(ElementDefinition(PRISM6,  3, 3, 1)); // node mapping?
  // em.add_def(ElementDefinition(PRISM18, 3, 3, 2)); // node mapping?

  // Eventually we'll map both tensor-product and non-tensor-product
  // BEXT elements to ours; for now we only support non-tensor-product
  // Bezier elements.
  // for (unsigned int i=0; i != 2; ++i)
  for (unsigned int i=1; i != 2; ++i)
    {
      // We only have one element for whom node orders match...
      em.add_def(ElementDefinition(EDGE2, i, 1, 1));

      em.add_def(ElementDefinition(EDGE3, i, 1, 2,
                                   {0, 2, 1}));
      em.add_def(ElementDefinition(EDGE4, i, 1, 2,
                                   {0, 2, 3, 1}));

      em.add_def(ElementDefinition(QUAD4, i, 2, 1,
                                   {0, 1, 3, 2}));
      em.add_def(ElementDefinition(QUAD9, i, 2, 2,
                                   {0, 4, 1, 7, 8, 5, 3, 6, 2}));

      em.add_def(ElementDefinition(HEX8,  i, 3, 1,
                                   {0, 1, 3, 2, 4, 5, 7, 6}));
      em.add_def(ElementDefinition(HEX27, i, 3, 2,
                                   { 0,  8,  1, 11, 20,  9,  3, 10,  2,
                                    12, 21, 13, 24, 26, 22, 15, 23, 14,
                                     4, 16,  5, 19, 25, 17,  7, 18,  6}));
    }

  return em;
}



DynaIO::ElementDefinition::ElementDefinition
  (ElemType type_in,
   unsigned dyna_type_in,
   unsigned dim_in,
   unsigned p_in) :
  type(type_in),
  dyna_type(dyna_type_in),
  dim(dim_in),
  p(p_in)
{
  const unsigned int n_nodes = Elem::type_to_n_nodes_map[type_in];
  nodes.resize(n_nodes);
  std::iota(nodes.begin(), nodes.end(), 0);
}


DynaIO::ElementDefinition::ElementDefinition
  (ElemType type_in,
   unsigned dyna_type_in,
   unsigned dim_in,
   unsigned p_in,
   std::vector<unsigned int> && nodes_in) :
  type(type_in),
  dyna_type(dyna_type_in),
  dim(dim_in),
  p(p_in),
  nodes(nodes_in)
{}



DynaIO::DynaIO (MeshBase & mesh) :
  MeshInput<MeshBase>  (mesh)
{
}



void DynaIO::read (const std::string & name)
{
  std::ifstream in (name.c_str());
  this->read_mesh (in);
}



void DynaIO::read_mesh(std::istream & in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (MeshInput<MeshBase>::mesh().processor_id(), 0);

  libmesh_assert(in.good());

  // clear any data in the mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();

  // Expect different sections, in this order, perhaps with blank
  // lines and/or comments in between:

  enum FileSection {
    FILE_HEADER,
    PATCH_HEADER,
    NODE_LINES,
    // (repeat NODE_LINES for each node)
    N_ELEM_SUBBLOCKS,
    ELEM_SUBBLOCK_HEADER,
    ELEM_NODES_LINES,
    ELEM_COEF_VEC_IDS,
    // (repeat nodes lines + coef vec ids as necessary)
    // (repeat elem subblock as necessary)
    N_COEF_BLOCKS, // number of coef vec blocks of each type
    N_VECS_PER_BLOCK, // number of coef vecs in each dense block
    COEF_VEC_COMPONENTS,
    //  (repeat coef vec components as necessary)
    //  (repeat coef blocks as necessary)
    //  reserved for sparse block stuff we don't support yet
    END_OF_FILE };

  FileSection section = FILE_HEADER;

  // Values to remember from section to section
  dyna_int_type patch_id, n_nodes, n_elem, n_coef_vec, weight_control_flag;
  dyna_int_type n_elem_blocks;
  dyna_int_type block_elem_type, block_n_elem, block_n_nodes, block_n_coef_vec, block_p,
                block_dim = 1;
  dyna_int_type n_dense_coef_vec_blocks, n_coef_vecs_in_subblock, n_coef_comp;
  unsigned char weight_index = 0;
  const ElementDefinition * current_elem_defn = nullptr;
  Elem * current_elem = nullptr;
  dyna_int_type n_nodes_read = 0,
                n_elem_blocks_read = 0,
                n_elems_read = 0,
                n_elem_nodes_read = 0,
                n_elem_cvids_read = 0,
                n_coef_blocks_read = 0,
                n_coef_comp_read = 0,
                n_coef_vecs_read = 0;

  // For reading the file line by line
  std::string s;

  while (true)
    {
      // Try to read something.  This may set EOF!
      std::getline(in, s);

      if (in)
        {
          // Process s...

          if (s.find("B E X T 2.0") == static_cast<std::string::size_type>(0))
          {
            libmesh_assert_equal_to(section, FILE_HEADER);
            section = PATCH_HEADER;
            continue;
          }

          // Ignore comments
          if (s.find("$#") == static_cast<std::string::size_type>(0))
            continue;

          // Ignore blank lines
          if (s.find_first_not_of(" \t") == std::string::npos)
            continue;

          // Parse each important section
          std::stringstream stream(s);
          switch (section) {
          case PATCH_HEADER:
            stream >> patch_id;
            stream >> n_nodes;
            stream >> n_elem;
            stream >> n_coef_vec;
            stream >> weight_control_flag;
            if (stream.fail())
              libmesh_error_msg("Failure to parse patch header\n");

            if (weight_control_flag)
              {
                weight_index = cast_int<unsigned char>
                  (mesh.add_node_datum<Real>("rational_weight"));
                mesh.set_default_mapping_type(RATIONAL_BERNSTEIN_MAP);
                mesh.set_default_mapping_data(weight_index);
              }

            section = NODE_LINES;
            break;
          case NODE_LINES:
            {
              std::array<dyna_fp_type, 4> xyzw;
              stream >> xyzw[0];
              stream >> xyzw[1];
              stream >> xyzw[2];
              stream >> xyzw[3];

              Point p(xyzw[0], xyzw[1], xyzw[2]);
              Node *n = mesh.add_point(p, n_nodes_read);
              if (weight_control_flag)
                n->set_extra_datum<Real>(weight_index, xyzw[3]);
            }
            ++n_nodes_read;

            if (stream.fail())
              libmesh_error_msg("Failure to parse node line\n");

            if (n_nodes_read >= n_nodes)
              section = N_ELEM_SUBBLOCKS;
            break;
          case N_ELEM_SUBBLOCKS:
            stream >> n_elem_blocks;
            if (stream.fail())
              libmesh_error_msg("Failure to parse n_elem_blocks\n");
            section = ELEM_SUBBLOCK_HEADER;
            break;
          case ELEM_SUBBLOCK_HEADER:
            stream >> block_elem_type;
            stream >> block_n_elem;
            stream >> block_n_nodes;
            stream >> block_n_coef_vec;
            stream >> block_p;

            if (stream.fail())
              libmesh_error_msg("Failure to parse elem block\n");

            block_dim = 1; // All blocks here are at least 1D

            dyna_int_type block_other_p; // Check for isotropic p
            stream >> block_other_p;
            if (!stream.fail())
              {
                block_dim = 2; // Found a second dimension!

                if (block_other_p != block_p)
                  libmesh_not_implemented(); // We don't support p anisotropy

                stream >> block_other_p;
                if (!stream.fail())
                  {
                    block_dim = 3; // Found a third dimension!

                    if (block_other_p != block_p)
                      libmesh_not_implemented();
                  }
              }
            n_elem_blocks_read++;
            n_elems_read = 0;
            section = ELEM_NODES_LINES;
            break;
          case ELEM_NODES_LINES:
            // Start a new Elem with each new line of nodes
            if (n_elem_nodes_read == 0)
              {
                // Consult the import element table to determine which element to build
                auto eletypes_it = _element_maps.in.find(std::make_tuple(block_elem_type, block_dim, block_p));

                // Make sure we actually found something
                if (eletypes_it == _element_maps.in.end())
                  libmesh_error_msg
                    ("Element of type " << block_elem_type <<
                     " dim " << block_dim <<
                     " degree " << block_p << " not found!");

                current_elem_defn = &(eletypes_it->second);
                current_elem = Elem::build(current_elem_defn->type).release();
                libmesh_assert_equal_to(current_elem->n_nodes(), (unsigned int)(block_n_nodes));
                libmesh_assert_equal_to(current_elem->dim(), block_dim);
                n_elem_cvids_read = 0;
              }
            {

              const int end_node_to_read =
                std::min(block_n_nodes, n_elem_nodes_read + max_ints_per_line);
              for (int i = n_elem_nodes_read; i != end_node_to_read; ++i)
                {
                  dyna_int_type node_id;
                  stream >> node_id;
                  node_id--;
                  current_elem->set_node(current_elem_defn->nodes[i]) =
                    mesh.node_ptr(node_id);

                  // Let's assume that our *only* mid-line breaks are
                  // due to the max_ints_per_line limit.  This should be
                  // less flexible but better for error detection.
                  if (stream.fail())
                    libmesh_error_msg("Failure to parse elem nodes\n");
                }
              if (end_node_to_read == block_n_nodes)
                {
                  n_elem_nodes_read = 0;
                  section = ELEM_COEF_VEC_IDS;
                }
              else
                {
                  n_elem_nodes_read = end_node_to_read;
                }
            }
            break;
          case ELEM_COEF_VEC_IDS:
            {
              const int end_cvid_to_read =
                std::min(block_n_nodes, n_elem_cvids_read + max_ints_per_line);
              for (int i = n_elem_cvids_read; i != end_cvid_to_read; ++i)
                {
                  dyna_int_type node_cvid;
                  stream >> node_cvid;
                  node_cvid--;

                  // FIXME - we need to store these somewhere to use for
                  // constraint equations

                  // Let's assume that our *only* mid-line breaks are
                  // due to the max_ints_per_line limit.  This should be
                  // less flexible but better for error detection.
                  if (stream.fail())
                    libmesh_error_msg("Failure to parse elem cvids\n");
                }
              if (end_cvid_to_read == block_n_nodes)
                {
                  current_elem->set_id(n_elems_read);
                  mesh.add_elem(current_elem);
                  n_elems_read++;
                  if (n_elems_read == block_n_elem)
                    section = N_COEF_BLOCKS; // Move on to coefficient vectors
                  else
                    section = ELEM_COEF_VEC_IDS; // Read another elem
                }
              else
                {
                  n_elem_cvids_read = end_cvid_to_read;
                }
            }
            break;
          case N_COEF_BLOCKS:
            {
              stream >> n_dense_coef_vec_blocks;
              dyna_int_type n_sparse_coef_vec_blocks;
              stream >> n_sparse_coef_vec_blocks;

              if (stream.fail())
                libmesh_error_msg("Failure to parse n_*_coef_vec_blocks\n");

              if (n_sparse_coef_vec_blocks != 0)
                libmesh_not_implemented();

              section = N_VECS_PER_BLOCK;
            }
            break;
          case N_VECS_PER_BLOCK:
            stream >> n_coef_vecs_in_subblock;
            stream >> n_coef_comp;

            if (stream.fail())
              libmesh_error_msg("Failure to parse dense coef subblock header\n");

            section = COEF_VEC_COMPONENTS;
            break;
          case COEF_VEC_COMPONENTS:
            // Start a new coefficient line
            if (n_coef_comp_read == 0)
              {
                // FIXME: allocate coef storage
              }
            {

              const int end_coef_to_read =
                std::min(n_coef_comp, n_coef_comp_read + max_fps_per_line);
              for (int i = n_coef_comp_read; i != end_coef_to_read; ++i)
                {
                  dyna_fp_type coef_comp;
                  stream >> coef_comp;

                  // FIXME: store coef

                  // Let's assume that our *only* mid-line breaks are
                  // due to the max_fps_per_line limit.  This should be
                  // less flexible but better for error detection.
                  if (stream.fail())
                    libmesh_error_msg("Failure to parse coefficients\n");
                }
              if (end_coef_to_read == n_coef_comp)
                {
                  n_coef_comp_read = 0;
                  n_coef_vecs_read++;
                  if (n_coef_vecs_read == n_coef_vecs_in_subblock)
                    {
                      n_coef_vecs_read = 0;
                      n_coef_blocks_read++;
                      if (n_coef_blocks_read == n_dense_coef_vec_blocks)
                        section = END_OF_FILE;
                    }
                }
              else
                {
                  n_coef_comp_read = end_coef_to_read;
                }
            }
            break;
          default:
            libmesh_error();
          }

          if (section == END_OF_FILE)
            break;
        } // if (in)
      else if (in.eof())
        libmesh_error_msg("Premature end of file");
      else
        libmesh_error_msg("Input stream failure! Perhaps the file does not exist?");
    }
}


} // namespace libMesh
