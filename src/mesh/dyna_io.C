// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/dyna_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/int_range.h"
#include "libmesh/utility.h"

// TIMPI includes
#include "timpi/parallel_implementation.h"

// gzstream for reading compressed files as a stream
#ifdef LIBMESH_HAVE_GZSTREAM
# include "libmesh/ignore_warnings.h" // shadowing in gzstream.h
# include "gzstream.h"
# include "libmesh/restore_warnings.h"
#endif

// C++ includes
#include <array>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <numeric> // iota

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
   dyna_int_type dyna_type_in,
   dyna_int_type dim_in,
   dyna_int_type p_in) :
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
   dyna_int_type dyna_type_in,
   dyna_int_type dim_in,
   dyna_int_type p_in,
   std::vector<unsigned int> && nodes_in) :
  type(type_in),
  dyna_type(dyna_type_in),
  dim(dim_in),
  p(p_in),
  nodes(nodes_in)
{}



DynaIO::DynaIO (MeshBase & mesh) :
  MeshInput<MeshBase>  (mesh),
  constraint_rows_broadcast (false)
{
}



void DynaIO::read (const std::string & name)
{
  const bool gzipped_file = (name.size() - name.rfind(".gz")  == 3);
  // These will be handled in unzip_file:
  // const bool bzipped_file = (name.size() - name.rfind(".bz2") == 4);
  // const bool xzipped_file = (name.size() - name.rfind(".xz") == 3);

  std::unique_ptr<std::istream> in;

  if (gzipped_file)
    {
#ifdef LIBMESH_HAVE_GZSTREAM
      igzstream * inf = new igzstream;
      libmesh_assert(inf);
      in.reset(inf);
      inf->open(name.c_str(), std::ios::in);
#else
      libmesh_error_msg("ERROR: need gzstream to handle .gz files!!!");
#endif
    }
  else
    {
      std::ifstream * inf = new std::ifstream;
      libmesh_assert(inf);
      in.reset(inf);

      std::string new_name = Utility::unzip_file(name);

      inf->open(new_name.c_str(), std::ios::in);
    }

  libmesh_assert(in.get());

  this->read_mesh (*in);
}



void DynaIO::read_mesh(std::istream & in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (MeshInput<MeshBase>::mesh().processor_id(), 0);

  libmesh_error_msg_if(!in.good(), "Can't read input stream");

  // clear any data in the mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();

  // clear any of our own data
  spline_node_ptrs.clear();
  constraint_rows.clear();
  constraint_rows_broadcast = false;

  // Expect different sections, in this order, perhaps with blank
  // lines and/or comments in between:

  enum FileSection {
    FILE_HEADER,
    PATCH_HEADER,
    NODE_LINES,
    // (repeat for each node)
    N_ELEM_SUBBLOCKS,
    ELEM_SUBBLOCK_HEADER,
    // (repeat for each subblock)
    ELEM_NODES_LINES,
    ELEM_COEF_VEC_IDS,
    // (repeat nodes lines + coef vec ids for each elem, subblock)
    N_COEF_BLOCKS, // number of coef vec blocks of each type
    N_VECS_PER_BLOCK, // number of coef vecs in each dense block
    COEF_VEC_COMPONENTS,
    //  (repeat coef vec components as necessary)
    //  (repeat coef blocks as necessary)
    //
    //  reserved for sparse block stuff we don't support yet
    END_OF_FILE };

  FileSection section = FILE_HEADER;

  // Values to remember from section to section
  dyna_int_type patch_id, n_spline_nodes, n_elem, n_coef_vec, weight_control_flag;
  dyna_int_type n_elem_blocks, n_dense_coef_vec_blocks;
  std::vector<dyna_int_type> // indexed from 0 to n_elem_blocks
    block_elem_type,
    block_n_elem,
    block_n_nodes,    // Number of *spline* nodes constraining elements
    block_n_coef_vec, // Number of coefficient vectors for each elem
    block_p,
    block_dim;
  std::vector<dyna_int_type> // indexed from 0 to n_dense_coef_vec_blocks
    n_coef_vecs_in_subblock, n_coef_comp;
  unsigned char weight_index = 0;
  dyna_int_type n_nodes_read = 0,
                n_elem_blocks_read = 0,
                n_elems_read = 0,
                n_elem_nodes_read = 0,
                n_elem_cvids_read = 0,
                n_coef_headers_read = 0,
                n_coef_blocks_read = 0,
                n_coef_comp_read = 0,
                n_coef_vecs_read = 0;

  // For reading the file line by line
  std::string s;

  // For storing global (spline) weights, until we have enough data to
  // use them for calculating local (Bezier element) nodes
  std::vector<Real> spline_weights;

  // For storing elements' constraint equations:
  // Global node indices (1-based in file, 0-based in memory):
  // elem_global_nodes[block_num][elem_num][local_index] = global_index
  std::vector<std::vector<std::vector<dof_id_type>>> elem_global_nodes;

  // Dense constraint vectors in the Dyna file
  // When first read:
  // dense_constraint_vecs[block_num][vec_in_block][column_num] = coef
  // When used:
  // dense_constraint_vecs[0][vec_num][column_num] = coef
  std::vector<std::vector<std::vector<Real>>> dense_constraint_vecs;

  // Constraint vector indices (1-based in file, 0-based in memory):
  // elem_constraint_rows[block_num][elem_num][row_num] = cv_index
  std::vector<std::vector<std::vector<dof_id_type>>> elem_constraint_rows;

  while (true)
    {
      // Try to read something.  This may set EOF!
      std::getline(in, s);

      if (in)
        {
          // Process s...

          if (s.find("B E X T 2.0") == static_cast<std::string::size_type>(0))
          {
            libmesh_error_msg_if(section != FILE_HEADER,
                                 "Found 'B E X T 2.0' outside file header?");

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
            stream >> n_spline_nodes;
            stream >> n_elem;
            stream >> n_coef_vec;
            stream >> weight_control_flag;
            libmesh_error_msg_if(stream.fail(), "Failure to parse patch header\n");

            spline_node_ptrs.resize(n_spline_nodes);
            spline_weights.resize(n_spline_nodes);

            if (weight_control_flag)
              {
                // If we ever add more nodes that aren't in this file,
                // merge this mesh with a non-spline mesh, etc., 1.0
                // is a good default weight
                const Real default_weight = 1.0;
                weight_index = cast_int<unsigned char>
                  (mesh.add_node_datum<Real>("rational_weight", true,
                                             &default_weight));
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

              if (weight_control_flag)
                spline_weights[n_nodes_read] = xyzw[3];

              // We'll use the spline nodes via NodeElem as the "true"
              // degrees of freedom, to which other Elem degrees of
              // freedom will be tied via constraint equations.
              Node *n = spline_node_ptrs[n_nodes_read] =
                mesh.add_point(Point(xyzw[0], xyzw[1], xyzw[2]));
              Elem * elem = mesh.add_elem(Elem::build(NODEELEM));
              elem->set_node(0) = n;
              elem->subdomain_id() = 1; // Separate id to ease Exodus output
            }
            ++n_nodes_read;

            libmesh_error_msg_if(stream.fail(), "Failure to parse node line\n");

            if (n_nodes_read >= n_spline_nodes)
              section = N_ELEM_SUBBLOCKS;
            break;
          case N_ELEM_SUBBLOCKS:
            stream >> n_elem_blocks;
            libmesh_error_msg_if(stream.fail(), "Failure to parse n_elem_blocks\n");

            block_elem_type.resize(n_elem_blocks);
            block_n_elem.resize(n_elem_blocks);
            block_n_nodes.resize(n_elem_blocks);
            block_n_coef_vec.resize(n_elem_blocks);
            block_p.resize(n_elem_blocks);
            block_dim.resize(n_elem_blocks);

            elem_global_nodes.resize(n_elem_blocks);
            elem_constraint_rows.resize(n_elem_blocks);

            n_elem_blocks_read = 0;
            section = ELEM_SUBBLOCK_HEADER;
            break;
          case ELEM_SUBBLOCK_HEADER:
            stream >> block_elem_type[n_elem_blocks_read];
            stream >> block_n_elem[n_elem_blocks_read];
            stream >> block_n_nodes[n_elem_blocks_read];
            stream >> block_n_coef_vec[n_elem_blocks_read];
            stream >> block_p[n_elem_blocks_read];

            libmesh_error_msg_if(stream.fail(), "Failure to parse elem block\n");

            block_dim[n_elem_blocks_read] = 1; // All blocks here are at least 1D

            dyna_int_type block_other_p; // Check for isotropic p
            stream >> block_other_p;
            if (!stream.fail())
              {
                block_dim[n_elem_blocks_read] = 2; // Found a second dimension!

                if (block_other_p != block_p[n_elem_blocks_read])
                  libmesh_not_implemented(); // We don't support p anisotropy

                stream >> block_other_p;
                if (!stream.fail())
                  {
                    block_dim[n_elem_blocks_read] = 3; // Found a third dimension!

                    if (block_other_p != block_p[n_elem_blocks_read])
                      libmesh_not_implemented();
                  }
              }

            {
              auto & block_global_nodes = elem_global_nodes[n_elem_blocks_read];
              auto & block_constraint_rows = elem_constraint_rows[n_elem_blocks_read];

              block_global_nodes.resize(block_n_elem[n_elem_blocks_read]);
              block_constraint_rows.resize(block_n_elem[n_elem_blocks_read]);

              for (auto e : make_range(block_n_elem[n_elem_blocks_read]))
                {
                  block_global_nodes[e].resize(block_n_nodes[n_elem_blocks_read]);
                  block_constraint_rows[e].resize(block_n_coef_vec[n_elem_blocks_read]);
                }
            }

            n_elem_blocks_read++;
            if (n_elem_blocks_read == n_elem_blocks)
              {
                n_elem_blocks_read = 0;
                n_elems_read = 0;
                section = ELEM_NODES_LINES;
              }
            break;
          case ELEM_NODES_LINES:
            {
              const int end_node_to_read =
                std::min(block_n_nodes[n_elem_blocks_read], n_elem_nodes_read + max_ints_per_line);
              for (; n_elem_nodes_read != end_node_to_read; ++n_elem_nodes_read)
                {
                  dyna_int_type node_id;
                  stream >> node_id;
                  node_id--;
                  elem_global_nodes[n_elem_blocks_read][n_elems_read][n_elem_nodes_read] = node_id;

                  // Let's assume that our *only* mid-line breaks are
                  // due to the max_ints_per_line limit.  This should be
                  // less flexible but better for error detection.
                  libmesh_error_msg_if(stream.fail(), "Failure to parse elem nodes\n");
                }

              if (n_elem_nodes_read == block_n_nodes[n_elem_blocks_read])
                {
                  n_elem_nodes_read = 0;
                  section = ELEM_COEF_VEC_IDS;
                }
            }
            break;
          case ELEM_COEF_VEC_IDS:
            {
              const int end_cvid_to_read =
                std::min(block_n_coef_vec[n_elem_blocks_read], n_elem_cvids_read + max_ints_per_line);
              for (; n_elem_cvids_read != end_cvid_to_read; ++n_elem_cvids_read)
                {
                  dyna_int_type node_cvid;
                  stream >> node_cvid;
                  node_cvid--;

                  elem_constraint_rows[n_elem_blocks_read][n_elems_read][n_elem_cvids_read] = node_cvid;

                  // Let's assume that our *only* mid-line breaks are
                  // due to the max_ints_per_line limit.  This should be
                  // less flexible but better for error detection.
                  libmesh_error_msg_if(stream.fail(), "Failure to parse elem cvids\n");
                }
              if (n_elem_cvids_read == block_n_nodes[n_elem_blocks_read])
                {
                  n_elem_cvids_read = 0;
                  n_elems_read++;
                  section = ELEM_NODES_LINES; // Read another elem, nodes first
                  if (n_elems_read == block_n_elem[n_elem_blocks_read])
                    {
                      n_elems_read = 0;
                      n_elem_blocks_read++;
                      if (n_elem_blocks_read == n_elem_blocks)
                        section = N_COEF_BLOCKS; // Move on to coefficient vectors
                    }
                }
            }
            break;
          case N_COEF_BLOCKS:
            {
              stream >> n_dense_coef_vec_blocks;
              dyna_int_type n_sparse_coef_vec_blocks;
              stream >> n_sparse_coef_vec_blocks;

              libmesh_error_msg_if(stream.fail(), "Failure to parse n_*_coef_vec_blocks\n");

              if (n_sparse_coef_vec_blocks != 0)
                libmesh_not_implemented();

              dense_constraint_vecs.resize(n_dense_coef_vec_blocks);
              n_coef_vecs_in_subblock.resize(n_dense_coef_vec_blocks);
              n_coef_comp.resize(n_dense_coef_vec_blocks);

              section = N_VECS_PER_BLOCK;
            }
            break;
          case N_VECS_PER_BLOCK:
            stream >> n_coef_vecs_in_subblock[n_coef_headers_read];
            stream >> n_coef_comp[n_coef_headers_read];

            libmesh_error_msg_if(stream.fail(), "Failure to parse dense coef subblock header\n");

            dense_constraint_vecs[n_coef_headers_read].resize
              (n_coef_vecs_in_subblock[n_coef_headers_read]);

            for (auto & vec : dense_constraint_vecs[n_coef_headers_read])
              vec.resize(n_coef_comp[n_coef_headers_read]);

            n_coef_headers_read++;
            if (n_coef_headers_read == n_dense_coef_vec_blocks)
              {
                n_coef_headers_read = 0;
                section = COEF_VEC_COMPONENTS;
              }
            break;
          case COEF_VEC_COMPONENTS:
            {
              auto & current_vec =
                dense_constraint_vecs[n_coef_blocks_read][n_coef_vecs_read];

              const int end_coef_to_read =
                std::min(n_coef_comp[n_coef_blocks_read],
                         n_coef_comp_read + max_fps_per_line);
              for (; n_coef_comp_read != end_coef_to_read; ++n_coef_comp_read)
                {
                  dyna_fp_type coef_comp;
                  stream >> coef_comp;

                  current_vec[n_coef_comp_read] = coef_comp;

                  // Let's assume that our *only* mid-line breaks are
                  // due to the max_fps_per_line limit.  This should be
                  // less flexible but better for error detection.
                  libmesh_error_msg_if(stream.fail(), "Failure to parse coefficients\n");
                }
              if (n_coef_comp_read == n_coef_comp[n_coef_blocks_read])
                {
                  n_coef_comp_read = 0;
                  n_coef_vecs_read++;
                  if (n_coef_vecs_read == n_coef_vecs_in_subblock[n_coef_blocks_read])
                    {
                      n_coef_vecs_read = 0;
                      n_coef_blocks_read++;
                      if (n_coef_blocks_read == n_dense_coef_vec_blocks)
                        section = END_OF_FILE;
                    }
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

  // Merge dense_constraint_vecs blocks
  if (n_dense_coef_vec_blocks)
    for (auto coef_vec_block :
         IntRange<dyna_int_type>(1, n_dense_coef_vec_blocks))
      {
        auto & dcv0 = dense_constraint_vecs[0];
        auto & dcvi = dense_constraint_vecs[coef_vec_block];
        dcv0.insert(dcv0.end(),
                    std::make_move_iterator(dcvi.begin()),
                    std::make_move_iterator(dcvi.end()));
      }
  dense_constraint_vecs.resize(1);

  // Constraint matrices:
  // elem_constraint_mat[block_num][elem_num][local_node_index][elem_global_nodes_index] = c
  std::vector<std::vector<std::vector<std::vector<Real>>>> elem_constraint_mat(n_elem_blocks);

  // We need to calculate local nodes on the fly, and we'll be
  // calculating them from constraint matrix columns, and we'll need
  // to make sure that the same node is found each time it's
  // calculated from multiple neighboring elements.
  std::map<std::vector<std::pair<dof_id_type, Real>>, Node *> local_nodes;

  for (auto block_num : make_range(n_elem_blocks))
    {
      elem_constraint_mat[block_num].resize(block_n_elem[block_num]);

      for (auto elem_num :
           make_range(block_n_elem[block_num]))
        {
          // Consult the import element table to determine which element to build
          auto eletypes_it =
            _element_maps.in.find(std::make_tuple(block_elem_type[block_num],
                                                  block_dim[block_num],
                                                  block_p[block_num]));

          // Make sure we actually found something
          libmesh_error_msg_if
            (eletypes_it == _element_maps.in.end(),
             "Element of type " << block_elem_type[block_num] <<
             " dim " << block_dim[block_num] <<
             " degree " << block_p[block_num] << " not found!");

          const ElementDefinition * elem_defn = &(eletypes_it->second);
          auto elem = Elem::build(elem_defn->type);
          libmesh_error_msg_if(elem->dim() != block_dim[block_num],
                               "Elem dim " << elem->dim() <<
                               " != block_dim " << block_dim[block_num]);

          auto & my_constraint_rows = elem_constraint_rows[block_num][elem_num];
          auto & my_global_nodes    = elem_global_nodes[block_num][elem_num];
          auto & my_constraint_mat  = elem_constraint_mat[block_num][elem_num];

          my_constraint_mat.resize(block_n_coef_vec[block_num]);
          for (auto spline_node_index :
               make_range(block_n_coef_vec[block_num]))
            my_constraint_mat[spline_node_index].resize(elem->n_nodes());

          for (auto spline_node_index :
               make_range(block_n_coef_vec[block_num]))
            {
              // Find which coef block this elem's vectors are from
              const dyna_int_type elem_coef_vec_index =
                my_constraint_rows[spline_node_index];

              dyna_int_type coef_block_num = 0;
              dyna_int_type first_block_coef_vec = 0;
              for (; elem_coef_vec_index >= first_block_coef_vec &&
                   coef_block_num != n_dense_coef_vec_blocks; ++coef_block_num)
                {
                  first_block_coef_vec += n_coef_vecs_in_subblock[coef_block_num];
                }

              // Make sure we did find a valid coef block
              libmesh_error_msg_if(coef_block_num == n_dense_coef_vec_blocks &&
                                   first_block_coef_vec <= elem_coef_vec_index,
                                   "Can't find valid constraint coef vector");

              coef_block_num--;

              libmesh_error_msg_if
                (dyna_int_type(elem->n_nodes()) != n_coef_comp[coef_block_num],
                 "Found " << n_coef_comp[coef_block_num] <<
                 " constraint coef vectors for " <<
                 elem->n_nodes() << " nodes");

              for (auto elem_node_index :
                   make_range(elem->n_nodes()))
                my_constraint_mat[spline_node_index][elem_node_index] =
                  dense_constraint_vecs[0][elem_coef_vec_index][elem_node_index];
            }

          for (auto elem_node_index :
               make_range(elem->n_nodes()))
            {
              dof_id_type global_node_idx = DofObject::invalid_id;

              // New finite element node data: dot product of
              // constraint matrix columns with spline node data.
              // Store that column as a key.
              std::vector<std::pair<dof_id_type, Real>> key;

              for (auto spline_node_index :
                   make_range(block_n_coef_vec[block_num]))
                {
                  const dyna_int_type elem_coef_vec_index =
                    my_constraint_rows[spline_node_index];

                  const Real coef =
                    libmesh_vector_at(dense_constraint_vecs[0],
                                      elem_coef_vec_index)[elem_node_index];

                  // Global nodes are supposed to be in sorted order
                  if (global_node_idx != DofObject::invalid_id)
                    libmesh_error_msg_if(my_global_nodes[spline_node_index] <= global_node_idx,
                                         "Found unsorted global node");

                  global_node_idx = my_global_nodes[spline_node_index];

                  if (coef != 0) // Ignore irrelevant spline nodes
                    key.emplace_back(global_node_idx, coef);
                }

              auto local_node_it = local_nodes.find(key);

              if (local_node_it != local_nodes.end())
                elem->set_node(elem_defn->nodes[elem_node_index]) =
                  local_node_it->second;
              else
                {
                  Point p(0);
                  Real w = 0;
                  std::vector<std::pair<dof_id_type, Real>> constraint_row;

                  for (auto spline_node_index :
                       make_range(block_n_coef_vec[block_num]))
                    {
                      const dof_id_type my_node_idx =
                        my_global_nodes[spline_node_index];

                      const dyna_int_type elem_coef_vec_index =
                        my_constraint_rows[spline_node_index];

                      const Node * spline_node =
                          libmesh_vector_at(spline_node_ptrs,
                                            my_node_idx);

                      const Real coef =
                        libmesh_vector_at(dense_constraint_vecs[0],
                                          elem_coef_vec_index)[elem_node_index];
                      p.add_scaled(*spline_node, coef);
                      w += coef * libmesh_vector_at(spline_weights,
                                                    my_node_idx);

                      constraint_row.emplace_back(spline_node->id(), coef);
                    }

                  Node *n = mesh.add_point(p);
                  if (weight_control_flag)
                    n->set_extra_datum<Real>(weight_index, w);
                  local_nodes[key] = n;
                  elem->set_node(elem_defn->nodes[elem_node_index]) = n;

                  constraint_rows[n->id()] = constraint_row;
                }
            }

          mesh.add_elem(std::move(elem));
        }
    }
}


void DynaIO::add_spline_constraints(DofMap & dof_map,
                                    unsigned int sys_num,
                                    unsigned int var_num)
{
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  MeshBase & mesh = this->mesh();

  // We have some strict compatibility requirements here still
  if (mesh.allow_renumbering() ||
      (!mesh.is_replicated() &&
       mesh.allow_remote_element_removal()))
    libmesh_not_implemented();

  // We do mesh reads in serial, and the mesh broadcast doesn't
  // broadcast our internal data, so we may need to do that now
  if (!constraint_rows_broadcast)
    {
      mesh.comm().broadcast(constraint_rows);
      constraint_rows_broadcast = true;
    }

  for (auto & node_row : constraint_rows)
    {
      DofConstraintRow dc_row;
      const Node * node = mesh.query_node_ptr(node_row.first);
      if (!node)
        continue;
      const dof_id_type constrained_id =
        node->dof_number(sys_num, var_num, 0);
      for (auto pr : node_row.second)
        {
          const Node & spline_node = mesh.node_ref(pr.first);
          const dof_id_type spline_dof_id =
            spline_node.dof_number(sys_num, var_num, 0);
          dc_row[spline_dof_id] = pr.second;
        }

      // Don't forbid constraint overwrite, or we're likely to
      // conflict with *any* other constraints.
      dof_map.add_constraint_row(constrained_id, dc_row, false);
    }

  dof_map.process_constraints(mesh);
  dof_map.prepare_send_list();
#else
  libmesh_ignore(dof_map, sys_num, var_num);
  libmesh_not_implemented();
#endif
}


} // namespace libMesh
