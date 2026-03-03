// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/exodusII_io.h"

#include "libmesh/boundary_info.h"
#include "libmesh/dof_map.h"
#include "libmesh/dyna_io.h"  // ElementDefinition for BEX
#include "libmesh/enum_elem_type.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io_helper.h"
#include "libmesh/fpe_disabler.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/parallel.h"
#include "libmesh/system.h"
#include "libmesh/utility.h"

// TIMPI includes
#include "timpi/parallel_sync.h"

// C++ includes
#include <cmath>   // llround
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>

#ifdef LIBMESH_HAVE_EXODUS_API
namespace
{
  using namespace libMesh;

  const std::vector<Real> & bex_constraint_vec(std::size_t i,
                                               const ExodusII_IO_Helper & helper)
  {
    std::size_t vec_offset = 0;
    for (auto block_num : index_range(helper.bex_dense_constraint_vecs))
    {
      const auto & vecblock = helper.bex_dense_constraint_vecs[block_num];
      libmesh_assert_greater_equal(i, vec_offset);
      if (i - vec_offset < vecblock.size())
        return vecblock[i - vec_offset];

      vec_offset += vecblock.size();
    }

    libmesh_error_msg("Requested BEX coefficient vector " << i << " not found");
  }

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  std::vector<Real>
  complex_soln_components (const std::vector<Number> & soln,
                           const unsigned int num_vars,
                           const bool write_complex_abs)
  {
    unsigned int num_values = soln.size();
    unsigned int num_elems = num_values / num_vars;

    // This will contain the real and imaginary parts and the magnitude
    // of the values in soln
    int nco = write_complex_abs ? 3 : 2;
    std::vector<Real> complex_soln(nco * num_values);

    for (unsigned i=0; i<num_vars; ++i)
      {
        for (unsigned int j=0; j<num_elems; ++j)
          {
            Number value = soln[i*num_vars + j];
            complex_soln[nco*i*num_elems + j] = value.real();
          }
        for (unsigned int j=0; j<num_elems; ++j)
          {
            Number value = soln[i*num_vars + j];
            complex_soln[nco*i*num_elems + num_elems + j] = value.imag();
          }
        if (write_complex_abs)
          {
            for (unsigned int j=0; j<num_elems; ++j)
              {
                Number value = soln[i*num_vars + j];
                complex_soln[3*i*num_elems + 2*num_elems + j] = std::abs(value);
              }
          }
      }

    return complex_soln;
  }
#endif // LIBMESH_USE_COMPLEX_NUMBERS

}
#endif

namespace libMesh
{

// ------------------------------------------------------------
// ExodusII_IO class members
ExodusII_IO::ExodusII_IO (MeshBase & mesh,
                          bool single_precision) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase> (mesh,
                        /* is_parallel_format = */ false,
                        /* serial_only_needed_on_proc_0 = */ true),
  ParallelObject(mesh),
#ifdef LIBMESH_HAVE_EXODUS_API
  exio_helper(std::make_unique<ExodusII_IO_Helper>(*this, false, true, single_precision)),
  _timestep(1),
  _verbose(false),
  _append(false),
#endif
  _allow_empty_variables(false),
  _write_complex_abs(true),
  _set_unique_ids_from_maps(false),
  _disc_bex(false)
{
  // if !LIBMESH_HAVE_EXODUS_API, we didn't use this
  libmesh_ignore(single_precision);
}



ExodusII_IO::ExodusII_IO (const MeshBase & mesh,
                          bool single_precision) :
  MeshInput<MeshBase> (),
  MeshOutput<MeshBase> (mesh,
                        /* is_parallel_format = */ false,
                        /* serial_only_needed_on_proc_0 = */ true),
  ParallelObject(mesh),
#ifdef LIBMESH_HAVE_EXODUS_API
  exio_helper(std::make_unique<ExodusII_IO_Helper>(*this, false, true, single_precision)),
  _timestep(1),
  _verbose(false),
  _append(false),
#endif
  _allow_empty_variables(false),
  _write_complex_abs(true),
  _set_unique_ids_from_maps(false),
  _disc_bex(false)
{
  // if !LIBMESH_HAVE_EXODUS_API, we didn't use this
  libmesh_ignore(single_precision);
}



int ExodusII_IO::get_exodus_version()
{
#ifdef LIBMESH_HAVE_EXODUS_API
  return ExodusII_IO_Helper::get_exodus_version();
#else
  return 0;
#endif
}



void ExodusII_IO::set_extra_integer_vars(const std::vector<std::string> & extra_integer_vars)
{
  _extra_integer_vars = extra_integer_vars;
}

void ExodusII_IO::set_output_variables(const std::vector<std::string> & output_variables,
                                       bool allow_empty)
{
  _output_variables = output_variables;
  _allow_empty_variables = allow_empty;
}



void ExodusII_IO::write_discontinuous_exodusII(const std::string & name,
                                               const EquationSystems & es,
                                               const std::set<std::string> * system_names)
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;

  es.build_variable_names  (solution_names, nullptr, system_names);
  es.build_discontinuous_solution_vector (v, system_names,
                                          nullptr, false, /* defaults */
                                          this->get_add_sides());
  this->write_nodal_data_discontinuous(name, v, solution_names);
}


#ifdef LIBMESH_HAVE_EXODUS_API
void ExodusII_IO::write_timestep_discontinuous (const std::string &fname,
                                                const EquationSystems &es,
                                                const int timestep,
                                                const Real time,
                                                const std::set<std::string> * system_names)
{
  _timestep = timestep;
  write_discontinuous_equation_systems (fname,es,system_names);

  if (MeshOutput<MeshBase>::mesh().processor_id())
    return;

  exio_helper->write_timestep(timestep, time);
}

#else
void ExodusII_IO::write_timestep_discontinuous (const std::string & /* fname */,
                                                const EquationSystems & /* es */,
                                                const int /* timestep */,
                                                const Real /* time */,
                                                const std::set<std::string> * /*system_names*/)
{ libmesh_error(); }
#endif


// ------------------------------------------------------------
// When the Exodus API is present...
#ifdef LIBMESH_HAVE_EXODUS_API

ExodusII_IO::~ExodusII_IO ()
{
  exio_helper->close();
}


void ExodusII_IO::read (const std::string & fname)
{
  LOG_SCOPE("read()", "ExodusII_IO");

  // Get a reference to the mesh we are reading
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Add extra integers into the mesh
  std::vector<unsigned int> extra_ids;
  for (auto & name : _extra_integer_vars)
    extra_ids.push_back(mesh.add_elem_integer(name));

  // Clear any existing mesh data
  mesh.clear();

  // Keep track of what kinds of elements this file contains
  elems_of_dimension.clear();
  elems_of_dimension.resize(4, false);

  // Open the exodus file in EX_READ mode
  exio_helper->open(fname.c_str(), /*read_only=*/true);

  // Get header information from exodus file
  exio_helper->read_and_store_header_info();

  // Read the QA records
  exio_helper->read_qa_records();

  // Print header information
  exio_helper->print_header();

  // Read nodes from the exodus file
  exio_helper->read_nodes();

  // Reserve space for the nodes.
  mesh.reserve_nodes(exio_helper->num_nodes);

  // Read the node number map from the Exodus file.  This is
  // required if we want to preserve the numbering of nodes as it
  // exists in the Exodus file.  If the Exodus file does not contain
  // a node_num_map, the identity map is returned by this call.
  exio_helper->read_node_num_map();

  // Read the element number map from the Exodus file.  This is
  // required if we want to preserve the numbering of elements as it
  // exists in the Exodus file.  If the Exodus file does not contain
  // an elem_num_map, the identity map is returned by this call.
  //
  // We now do this before creating nodes, so if we have any spline
  // nodes that need a NodeElem attached we can give them unused elem
  // ids.
  exio_helper->read_elem_num_map();

  // Read any Bezier Extraction coefficient vectors from the file,
  // such as might occur in an IsoGeometric Analysis (IGA) mesh.
  exio_helper->read_bex_cv_blocks();

  // If we have Rational Bezier weights, we'll need to
  // store them.
  unsigned char weight_index = 0;
  const bool weights_exist = !exio_helper->w.empty();

  // If we have Bezier extraction coefficients, we'll need to put
  // NODEELEM elements on spline nodes, since our Rational Bezier
  // elements will be connected to nodes derived from those; nothing
  // else will be directly connected to the spline nodes.
  const bool bex_cv_exist = !exio_helper->bex_dense_constraint_vecs.empty();

  // If the user requests to set Node/Elem unique_ids based on the
  // node/elem_num_maps, then it is very likely that those unique_ids
  // will not be "unique" across the entire set of DofObjects (Nodes
  // and Elems), although they should of course still be unique within
  // the set of Nodes and the set of Elems on an individual basis.  We
  // therefore need to relax some of the uniqueness checking that is
  // done in debug mode by setting the following flag on the Mesh.
  if (_set_unique_ids_from_maps)
    mesh.allow_node_and_elem_unique_id_overlap(true);

  // Even if weights don't exist, we still use RATIONAL_BERNSTEIN for
  // Bezier-Bernstein BEX elements, we just use 1.0 as the weight on
  // every node.
  if (bex_cv_exist)
    {
      const Real default_weight = 1.0;
      weight_index = cast_int<unsigned char>
        (mesh.add_node_datum<Real>("rational_weight", true,
                                   &default_weight));
      mesh.set_default_mapping_type(RATIONAL_BERNSTEIN_MAP);
      mesh.set_default_mapping_data(weight_index);
    }

  std::unordered_map<const Node *, Elem *> spline_nodeelem_ptrs;

  // Loop over the nodes, create Nodes with local processor_id 0.
  for (auto i : make_range(exio_helper->num_nodes))
    {
      // Determine the libmesh node id implied by "i". The
      // get_libmesh_node_id() helper function expects a 1-based
      // Exodus node id, so we construct the "implied" Exodus node id
      // from "i" by adding 1.
      auto libmesh_node_id = exio_helper->get_libmesh_node_id(/*exodus_node_id=*/i+1);

      // Catch the node that was added to the mesh
      Node * added_node = mesh.add_point (Point(exio_helper->x[i], exio_helper->y[i], exio_helper->z[i]), libmesh_node_id);

      // Sanity check: throw an error if the Mesh assigned an ID to
      // the Node which does not match the libmesh_id we just determined.
      libmesh_error_msg_if(added_node->id() != static_cast<unsigned>(libmesh_node_id),
                           "Error!  Mesh assigned node ID "
                           << added_node->id()
                           << " which is different from the (zero-based) Exodus ID "
                           << libmesh_node_id
                           << "!");

      // If the _set_unique_ids_from_maps flag is true, set the
      // unique_id for "node", otherwise do nothing.
      exio_helper->conditionally_set_node_unique_id(mesh, added_node, i);

      // If we have a set of spline weights, these nodes are going to
      // be used as control points for Bezier elements, and we need
      // to attach a NodeElem to each to make sure it doesn't get
      // flagged as an unused node.
      if (weights_exist)
        {
          const auto w = exio_helper->w[i];
          Point & p = *added_node;
          p /= w;  // Exodus Bezier Extraction stores spline nodes in projective space

          added_node->set_extra_datum<Real>(weight_index, exio_helper->w[i]);
        }

      if (bex_cv_exist)
        {
          std::unique_ptr<Elem> elem = Elem::build(NODEELEM);

          // Give the NodeElem ids at the end, so we can match any
          // existing ids in the file for other elements
          //
          // We don't set the unique_id for this NodeElem here even if
          // the user has set the _set_unique_ids_from_maps flag
          // because these NodeElems don't have entries in the
          // elem_num_map. Therefore, we just let the Mesh assign
          // whatever unique_id is "next" as the Elem is added to the
          // Mesh.
          elem->set_id() = exio_helper->end_elem_id() + i;

          elem->set_node(0, added_node);
          Elem * added_elem = mesh.add_elem(std::move(elem));
          spline_nodeelem_ptrs[added_node] = added_elem;
        }
    }

  // This assert is no longer valid if the nodes are not numbered
  // sequentially starting from 1 in the Exodus file.
  // libmesh_assert_equal_to (static_cast<unsigned int>(exio_helper->num_nodes), mesh.n_nodes());

  // Get information about all the element and edge blocks
  exio_helper->read_block_info();

  // Reserve space for the elements.  Account for any NodeElem that
  // have already been attached to spline control nodes.
  mesh.reserve_elem(exio_helper->w.size() + exio_helper->num_elem);

  exio_helper->read_num_time_steps();

  // Read variables for extra integer IDs
  auto extra_integer_values =
    exio_helper->read_extra_integers(_extra_integer_vars);

  // Read in the element connectivity for each block.
  int nelem_last_block = 0;

  // If we're building Bezier elements from spline nodes, we need to
  // calculate those elements' local nodes on the fly, and we'll be
  // calculating them from constraint matrix columns, and we'll need
  // to make sure that the same node is found each time it's
  // calculated from multiple neighboring elements.
  std::map<std::vector<std::pair<dof_id_type, Real>>, Node *> local_nodes;

  // We'll set any spline NodeElem subdomain_id() values to exceed the
  // maximum of subdomain_id() values set via Exodus block ids.
  int max_subdomain_id = std::numeric_limits<int>::min();

  // We've already added all the nodes explicitly specified in the
  // file, but if we have spline nodes we may need to add assembly
  // element nodes based on them.  Add them contiguously so we're
  // compatible with any subsequent code paths (including our
  // ExodusII_IO::write()!) that don't support sparse ids.
  dof_id_type n_nodes = mesh.n_nodes();

  // Loop over all the element blocks
  for (int i=0; i<exio_helper->num_elem_blk; i++)
    {
      // Read the information for block i
      exio_helper->read_elem_in_block (i);
      const int subdomain_id = exio_helper->get_block_id(i);
      max_subdomain_id = std::max(max_subdomain_id, subdomain_id);

      // populate the map of names
      std::string subdomain_name = exio_helper->get_block_name(i);
      if (!subdomain_name.empty())
        mesh.subdomain_name(static_cast<subdomain_id_type>(subdomain_id)) = subdomain_name;

      // Set any relevant node/edge maps for this element
      const std::string type_str (exio_helper->get_elem_type());
      const auto & conv = exio_helper->get_conversion(type_str);

      // Loop over all the faces in this block
      int jmax = nelem_last_block+exio_helper->num_elem_this_blk;
      for (int j=nelem_last_block; j<jmax; j++)
        {
          auto uelem = Elem::build(conv.libmesh_elem_type());

          const int elem_num = j - nelem_last_block;

          // Make sure that Exodus's number of nodes per Elem matches
          // the number of Nodes for this type of Elem. We only check
          // this for the first Elem in each block, since these values
          // are the same for every Elem in the block.
          if (!elem_num)
            libmesh_error_msg_if(exio_helper->num_nodes_per_elem != static_cast<int>(uelem->n_nodes()),
                                 "Error: Exodus file says "
                                 << exio_helper->num_nodes_per_elem
                                 << " nodes per Elem, but Elem type "
                                 << Utility::enum_to_string(uelem->type())
                                 << " has " << uelem->n_nodes() << " nodes.");

          // Assign the current subdomain to this Elem
          uelem->subdomain_id() = static_cast<subdomain_id_type>(subdomain_id);

          // Determine the libmesh elem id implied by "j". The
          // ExodusII_IO_Helper::get_libmesh_elem_id() helper function
          // expects a 1-based Exodus elem id, so we construct the
          // "implied" Exodus elem id from "j" by adding 1.
          auto libmesh_elem_id = exio_helper->get_libmesh_elem_id(/*exodus_elem_id=*/j+1);

          uelem->set_id(libmesh_elem_id);

          // Record that we have seen an element of dimension uelem->dim()
          elems_of_dimension[uelem->dim()] = true;

          // Catch the Elem pointer that the Mesh throws back
          Elem * elem = mesh.add_elem(std::move(uelem));

          // If the _set_unique_ids_from_maps flag is true, set the
          // unique_id for "elem", otherwise do nothing.
          exio_helper->conditionally_set_elem_unique_id(mesh, elem, j);

          // If the Mesh assigned an ID different from the one we
          // tried to give it, we should probably error.
          libmesh_error_msg_if(elem->id() != static_cast<unsigned>(libmesh_elem_id),
                               "Error!  Mesh assigned ID "
                               << elem->id()
                               << " which is different from the (zero-based) Exodus ID "
                               << libmesh_elem_id
                               << "!");

          // Assign extra integer IDs
          for (auto & id : extra_ids)
             elem->set_extra_integer(id, extra_integer_values[id][elem->id()]);

          // Set all the nodes for this element
          //
          // If we don't have any Bezier extraction operators, this
          // is easy: we've already built all our nodes and just need
          // to link to them.
          if (exio_helper->bex_cv_conn.empty())
            {
              for (int k=0; k<exio_helper->num_nodes_per_elem; k++)
                {
                  // Get index into this block's connectivity array
                  int gi = elem_num * exio_helper->num_nodes_per_elem + conv.get_node_map(k);

                  // Get the 1-based Exodus node id from the "connect" array
                  auto exodus_node_id = exio_helper->connect[gi];

                  // Convert this index to a libMesh Node id
                  auto libmesh_node_id = exio_helper->get_libmesh_node_id(exodus_node_id);

                  // Set the node pointer in the Elem
                  elem->set_node(k, mesh.node_ptr(libmesh_node_id));
                }
            }
          else // We have Bezier Extraction data
            {
              auto & constraint_rows = mesh.get_constraint_rows();

              const DynaIO::ElementDefinition & dyna_elem_defn =
                DynaIO::find_elem_definition(elem->type(),
                                             elem->dim(),
                                             int(elem->default_order()));

              std::vector<std::vector<Real>>
                my_constraint_mat(exio_helper->bex_num_elem_cvs);
              for (auto spline_node_index :
                   make_range(exio_helper->bex_num_elem_cvs))
                {
                  my_constraint_mat[spline_node_index].resize(elem->n_nodes());

                  const auto & my_constraint_rows = exio_helper->bex_cv_conn[elem_num];
                  const unsigned long elem_coef_vec_index =
                    my_constraint_rows[spline_node_index] - 1; // Exodus isn't 0-based
                  const auto & my_vec = bex_constraint_vec(elem_coef_vec_index, *exio_helper);
                  for (auto elem_node_index :
                       make_range(elem->n_nodes()))
                    {
                      my_constraint_mat[spline_node_index][elem_node_index] =
                        my_vec[elem_node_index];
                    }

                }

              // The tailing entries in each element's connectivity
              // vector are indices to Exodus constraint coefficient
              // rows.

              // Concatenating these rows gives a matrix with
              // all the constraints for the element nodes: each
              // column of that matrix is the constraint coefficients
              // for the node associated with that column (via the
              // Exodus numbering, not the libMesh numbering).
              const auto & my_constraint_rows = exio_helper->bex_cv_conn[elem_num];

              for (auto elem_node_index :
                   make_range(elem->n_nodes()))
                {
                  // New finite element node data = dot product of
                  // constraint matrix columns with spline node data.
                  // Store each column's non-zero entries, along with
                  // the global spline node indices, as a key to
                  // identify shared finite element nodes.
                  std::vector<std::pair<dof_id_type, Real>> key;

                  for (auto spline_node_index :
                       make_range(exio_helper->bex_num_elem_cvs))
                    {
                      // Pick out a row of the element constraint matrix
                      const unsigned long elem_coef_vec_index =
                        my_constraint_rows[spline_node_index] - 1; // Exodus isn't 0-based

                      auto & coef_vec =
                        bex_constraint_vec(elem_coef_vec_index, *exio_helper);

                      // Get coef from this node's column intersect that row
                      const Real coef =
                        libmesh_vector_at(coef_vec, elem_node_index);

                      // Get the libMesh node corresponding to that row
                      const int gi = elem_num * exio_helper->bex_num_elem_cvs + spline_node_index;

                      // Get the 1-based Exodus node id from the "connect" array
                      auto exodus_node_id = exio_helper->connect[gi];

                      // Convert this index to a libMesh Node id
                      auto libmesh_node_id = exio_helper->get_libmesh_node_id(exodus_node_id);

                      if (coef != 0) // Ignore irrelevant spline nodes
                        key.emplace_back(libmesh_node_id, coef);
                    }

                  // Have we already created this node?  Connect it.
                  if (const auto local_node_it = local_nodes.find(key);
                      local_node_it != local_nodes.end())
                    elem->set_node(dyna_elem_defn.nodes[elem_node_index], local_node_it->second);
                  // Have we not yet created this node?  Construct it,
                  // along with its weight and libMesh constraint row,
                  // then connect it.
                  else
                    {
                      Point p(0);
                      Real w = 0;
                      std::vector<std::pair<std::pair<const Elem *, unsigned int>, Real>> constraint_row;

                      for (auto [libmesh_spline_node_id, coef] : key)
                        {
                          const Node & spline_node = mesh.node_ref(libmesh_spline_node_id);

                          p.add_scaled(spline_node, coef);
                          const Real spline_w = weights_exist ?
                            spline_node.get_extra_datum<Real>(weight_index) : 1;
                          w += coef * spline_w;

                          const Elem * nodeelem =
                            libmesh_map_find(spline_nodeelem_ptrs, &spline_node);
                          constraint_row.emplace_back(std::make_pair(nodeelem, 0), coef);
                        }

                      Node *n = mesh.add_point(p, n_nodes++);
                      if (weights_exist)
                        n->set_extra_datum<Real>(weight_index, w);

                      // If we're building disconnected Bezier
                      // extraction elements then we don't want to
                      // find the new nodes to reuse later; each
                      // finite element node will connect to only one
                      // element.
                      if (!_disc_bex)
                        local_nodes[key] = n;
                      elem->set_node(dyna_elem_defn.nodes[elem_node_index], n);

                      constraint_rows[n] = constraint_row;
                    }
                }
            }
        }

      // running sum of # of elements per block,
      // (should equal total number of elements in the end)
      nelem_last_block += exio_helper->num_elem_this_blk;
    }

  // Now we know enough to fix any spline NodeElem subdomains
  max_subdomain_id++;
  for (auto p : spline_nodeelem_ptrs)
    p.second->subdomain_id() = max_subdomain_id;

  // Read in edge blocks, storing information in the BoundaryInfo object.
  // Edge blocks are treated as BCs.
  exio_helper->read_edge_blocks(mesh);

  // Set the mesh dimension to the largest encountered for an element
  for (unsigned char i=0; i!=4; ++i)
    if (elems_of_dimension[i])
      mesh.set_mesh_dimension(i);

  // Read in sideset information -- this is useful for applying boundary conditions
  {
    // Get basic information about all sidesets
    exio_helper->read_sideset_info();
    int offset=0;
    for (int i=0; i<exio_helper->num_side_sets; i++)
      {
        // Compute new offset
        offset += (i > 0 ? exio_helper->num_sides_per_set[i-1] : 0);
        exio_helper->read_sideset (i, offset);

        std::string sideset_name = exio_helper->get_side_set_name(i);
        if (!sideset_name.empty())
          mesh.get_boundary_info().sideset_name
            (cast_int<boundary_id_type>(exio_helper->get_side_set_id(i)))
            = sideset_name;
      }

    for (auto e : index_range(exio_helper->elem_list))
      {
        // Call helper function to get the libmesh Elem id for the
        // e'th entry in the current elem_list.
        dof_id_type libmesh_elem_id =
          exio_helper->get_libmesh_elem_id(exio_helper->elem_list[e]);

        // Set any relevant node/edge maps for this element
        Elem & elem = mesh.elem_ref(libmesh_elem_id);

        const auto & conv = exio_helper->get_conversion(elem.type());

        // Map the zero-based Exodus side numbering to the libmesh side numbering
        unsigned int raw_side_index = exio_helper->side_list[e]-1;
        std::size_t side_index_offset = conv.get_shellface_index_offset();

        if (raw_side_index < side_index_offset)
          {
            // We assume this is a "shell face"
            int mapped_shellface = raw_side_index;

            // Check for errors
            libmesh_error_msg_if(mapped_shellface < 0 || mapped_shellface >= 2,
                                 "Bad 0-based shellface id: "
                                 << mapped_shellface
                                 << " detected in Exodus file "
                                 << exio_helper->current_filename);

            // Add this (elem,shellface,id) triplet to the BoundaryInfo object.
            mesh.get_boundary_info().add_shellface (libmesh_elem_id,
                                                    cast_int<unsigned short>(mapped_shellface),
                                                    cast_int<boundary_id_type>(exio_helper->id_list[e]));
          }
        else
          {
            unsigned int side_index = static_cast<unsigned int>(raw_side_index - side_index_offset);
            int mapped_side = conv.get_side_map(side_index);

            // Check for errors
            libmesh_error_msg_if(mapped_side == ExodusII_IO_Helper::Conversion::invalid_id,
                                 "Invalid 1-based side id: "
                                 << side_index
                                 << " detected for "
                                 << Utility::enum_to_string(elem.type())
                                 << " in Exodus file "
                                 << exio_helper->current_filename);

            libmesh_error_msg_if(mapped_side < 0 ||
                                 cast_int<unsigned int>(mapped_side) >= elem.n_sides(),
                                 "Bad 0-based side id: "
                                 << mapped_side
                                 << " detected for "
                                 << Utility::enum_to_string(elem.type())
                                 << " in Exodus file "
                                 << exio_helper->current_filename);

            // Add this (elem,side,id) triplet to the BoundaryInfo object.
            mesh.get_boundary_info().add_side (libmesh_elem_id,
                                               cast_int<unsigned short>(mapped_side),
                                               cast_int<boundary_id_type>(exio_helper->id_list[e]));
          }
      } // end for (elem_list)
  } // end read sideset info

  // Read in elemset information and apply to Mesh elements if present
  {
    exio_helper->read_elemset_info();

    // Mimic behavior of sideset case where we store all the set
    // information in a single array with offsets.
    int offset=0;
    for (int i=0; i<exio_helper->num_elem_sets; i++)
      {
        // Compute new offset
        offset += (i > 0 ? exio_helper->num_elems_per_set[i-1] : 0);
        exio_helper->read_elemset (i, offset);

        // TODO: add support for elemset names
        // std::string elemset_name = exio_helper->get_elem_set_name(i);
        // if (!elemset_name.empty())
        //   mesh.get_boundary_info().elemset_name(cast_int<boundary_id_type>(exio_helper->get_elem_set_id(i))) = elemset_name;
      }

    // Debugging: print the concatenated list of elemset ids
    // libMesh::out << "Concatenated list of elemset Elem ids (Exodus numbering):" << std::endl;
    // for (const auto & id : exio_helper->elemset_list)
    //   libMesh::out << id << " ";
    // libMesh::out << std::endl;

    // Next we need to assign the elemset ids to the mesh using the
    // Elem's "extra_integers" support, if we have any.
    if (exio_helper->num_elem_all_elemsets)
      {
        // Build map from Elem -> {elemsets}. This is needed only
        // temporarily to determine a unique set of elemset codes.
        std::map<Elem *, MeshBase::elemset_type> elem_to_elemsets;
        for (auto e : index_range(exio_helper->elemset_list))
          {
           // Call helper function to get the libmesh Elem id for the
           // e'th entry in the current elemset_list.
           dof_id_type libmesh_elem_id =
             exio_helper->get_libmesh_elem_id(exio_helper->elemset_list[e]);

            // Get a pointer to this Elem
            Elem * elem = mesh.elem_ptr(libmesh_elem_id);

            // Debugging:
            // libMesh::out << "Elem " << elem->id() << " is in elemset " << exio_helper->elemset_id_list[e] << std::endl;

            // Store elemset id in the map
            elem_to_elemsets[elem].insert(exio_helper->elemset_id_list[e]);
          }

        // Create a set of unique elemsets
        std::set<MeshBase::elemset_type> unique_elemsets;
        for (const auto & pr : elem_to_elemsets)
          unique_elemsets.insert(pr.second);

        // Debugging: print the unique elemsets
        // libMesh::out << "The set of unique elemsets which exist on the Mesh:" << std::endl;
        // for (const auto & s : unique_elemsets)
        //   {
        //     for (const auto & elemset_id : s)
        //       libMesh::out << elemset_id << " ";
        //     libMesh::out << std::endl;
        //   }

        // Enumerate the unique_elemsets and tell the mesh about them
        dof_id_type code = 0;
        for (const auto & s : unique_elemsets)
          mesh.add_elemset_code(code++, s);

        // Sanity check: make sure that MeshBase::n_elemsets() reports
        // the expected value after calling MeshBase::add_elemset_code()
        // one or more times.
        libmesh_assert_msg(exio_helper->num_elem_sets == cast_int<int>(mesh.n_elemsets()),
                           "Error: mesh.n_elemsets() is " << mesh.n_elemsets()
                           << ", but mesh should have " << exio_helper->num_elem_sets << " elemsets.");

        // Create storage for the extra integer on all Elems. Elems which
        // are not in any set will use the default value of DofObject::invalid_id
        unsigned int elemset_index =
          mesh.add_elem_integer("elemset_code",
                                /*allocate_data=*/true);

        // Store the appropriate extra_integer value on all Elems that need it.
        for (const auto & [elem, s] : elem_to_elemsets)
          elem->set_extra_integer(elemset_index, mesh.get_elemset_code(s));
      }
  } // done reading elemset info

  // Read nodeset info
  {
    // This fills in the following fields of the helper for later use:
    // nodeset_ids
    // num_nodes_per_set
    // num_node_df_per_set
    // node_sets_node_index
    // node_sets_dist_index
    // node_sets_node_list
    // node_sets_dist_fact
    exio_helper->read_all_nodesets();

    for (int nodeset=0; nodeset<exio_helper->num_node_sets; nodeset++)
      {
        boundary_id_type nodeset_id =
          cast_int<boundary_id_type>(exio_helper->nodeset_ids[nodeset]);

        std::string nodeset_name = exio_helper->get_node_set_name(nodeset);
        if (!nodeset_name.empty())
          mesh.get_boundary_info().nodeset_name(nodeset_id) = nodeset_name;

        // Get starting index of node ids for current nodeset.
        unsigned int offset = exio_helper->node_sets_node_index[nodeset];

        for (int i=0; i<exio_helper->num_nodes_per_set[nodeset]; ++i)
          {
            int exodus_node_id = exio_helper->node_sets_node_list[i + offset];
            auto libmesh_node_id = exio_helper->get_libmesh_node_id(exodus_node_id);
            mesh.get_boundary_info().add_node(libmesh_node_id, nodeset_id);
          }
      }
  }

#if LIBMESH_DIM < 3
  libmesh_error_msg_if(mesh.mesh_dimension() > LIBMESH_DIM,
                       "Cannot open dimension "
                       << mesh.mesh_dimension()
                       << " mesh file when configured without "
                       << mesh.mesh_dimension()
                       << "D support.");
#endif
}



ExodusHeaderInfo
ExodusII_IO::read_header (const std::string & fname)
{
  // We will need the Communicator of the Mesh we were created with.
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Eventual return value
  ExodusHeaderInfo header_info;

  // File I/O is done on processor 0, then broadcast to other procs
  if (mesh.processor_id() == 0)
    {
      // Open the exodus file in EX_READ mode
      exio_helper->open(fname.c_str(), /*read_only=*/true);

      // Get header information from exodus file without updating the
      // Helper object's internal data structures.
      header_info = exio_helper->read_header();

      // Close the file, we are now done with it. The goal is to keep the
      // exio_helper object unchanged while calling this function,
      // although it can't quite be marked "const" because we do have to
      // actually open/close the file. This way, it should be possible to
      // use the same ExodusII_IO object to read the headers of multiple
      // different mesh files.
      exio_helper->close();
    }

  // Broadcast header_info to other procs before returning
  header_info.broadcast(mesh.comm());

  // Return the information we read back to the user.
  return header_info;
}



void ExodusII_IO::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;

  // Set the verbose flag in the helper object as well.
  exio_helper->verbose = _verbose;
}



void ExodusII_IO::write_complex_magnitude (bool val)
{
  _write_complex_abs = val;
}

void ExodusII_IO::set_unique_ids_from_maps (bool val)
{
  _set_unique_ids_from_maps = val;

  // Set this flag on the helper object as well. The helper needs to know about this
  // flag, since it sometimes needs to construct libmesh Node ids from nodal connectivity
  // arrays (see e.g. ExodusII_IO_Helper::read_edge_blocks()).
  exio_helper->set_unique_ids_from_maps = val;
}

void ExodusII_IO::use_mesh_dimension_instead_of_spatial_dimension(bool val)
{
  exio_helper->use_mesh_dimension_instead_of_spatial_dimension(val);
}



void ExodusII_IO::write_as_dimension(unsigned dim)
{
  exio_helper->write_as_dimension(dim);
}



void ExodusII_IO::set_coordinate_offset(Point p)
{
  libmesh_warning("This method may be deprecated in the future");
  exio_helper->set_coordinate_offset(p);
}



void ExodusII_IO::append(bool val)
{
  _append = val;
}



void ExodusII_IO::write_added_sides (bool val)
{
  exio_helper->set_add_sides(val);
}



bool ExodusII_IO::get_add_sides ()
{
  return exio_helper->get_add_sides();
}




const std::vector<Real> & ExodusII_IO::get_time_steps()
{
  libmesh_error_msg_if
    (!exio_helper->opened_for_reading,
     "ERROR, ExodusII file must be opened for reading before calling ExodusII_IO::get_time_steps()!");

  exio_helper->read_time_steps();
  return exio_helper->time_steps;
}



int ExodusII_IO::get_num_time_steps()
{
  libmesh_error_msg_if(!exio_helper->opened_for_reading && !exio_helper->opened_for_writing,
                       "ERROR, ExodusII file must be opened for reading or writing before calling ExodusII_IO::get_num_time_steps()!");

  exio_helper->read_num_time_steps();
  return exio_helper->num_time_steps;
}



void ExodusII_IO::copy_nodal_solution(System & system,
                                      std::string system_var_name,
                                      std::string exodus_var_name,
                                      unsigned int timestep)
{
  LOG_SCOPE("copy_nodal_solution()", "ExodusII_IO");

  const unsigned int var_num = system.variable_number(system_var_name);

  const MeshBase & mesh = MeshInput<MeshBase>::mesh();

  libmesh_error_msg_if(mesh.allow_renumbering(),
                       "ERROR, nodal data cannot be loaded if the mesh may be renumbered!");

  // With Exodus files we only open them on processor 0, so that's the
  // where we have to do the data read too.
  if (system.comm().rank() == 0)
    {
      libmesh_error_msg_if(!exio_helper->opened_for_reading,
                           "ERROR, ExodusII file must be opened for reading before copying a nodal solution!");

      exio_helper->read_nodal_var_values(exodus_var_name, timestep);
    }

  auto & node_var_value_map = exio_helper->nodal_var_values;

  const bool serial_on_zero = mesh.is_serial_on_zero();

  // If our mesh isn't serial, then non-root processors need to
  // request the data for their parts of the mesh and insert it
  // themselves.
  if (!serial_on_zero)
    {
      std::unordered_map<processor_id_type, std::vector<dof_id_type>> node_ids_to_request;
      if (this->processor_id() != 0)
        {
          std::vector<dof_id_type> node_ids;
          for (auto & node : mesh.local_node_ptr_range())
            node_ids.push_back(node->id());
          if (!node_ids.empty())
            node_ids_to_request[0] = std::move(node_ids);
        }

      auto value_gather_functor =
        [& node_var_value_map]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         std::vector<Real> & values)
        {
          const std::size_t query_size = ids.size();
          values.resize(query_size);
          for (std::size_t i=0; i != query_size; ++i)
            {
              if (const auto it = node_var_value_map.find(ids[i]);
                  it != node_var_value_map.end())
                {
                  values[i] = it->second;
                  node_var_value_map.erase(it);
                }
              else
                values[i] = std::numeric_limits<Real>::quiet_NaN();
            }
        };

      auto value_action_functor =
        [& node_var_value_map]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         const std::vector<Real> & values)
        {
          const std::size_t query_size = ids.size();
          for (std::size_t i=0; i != query_size; ++i)
            if (!libmesh_isnan(values[i]))
              node_var_value_map[ids[i]] = values[i];
        };

      Real * value_ex = nullptr;
      Parallel::pull_parallel_vector_data
        (system.comm(), node_ids_to_request, value_gather_functor,
         value_action_functor, value_ex);
    }

  // Everybody inserts the data they've received.  If we're
  // serial_on_zero then proc 0 inserts everybody's data and other
  // procs have empty map ranges.
  for (auto p : exio_helper->nodal_var_values)
    {
      dof_id_type i = p.first;
      const Node * node = MeshInput<MeshBase>::mesh().query_node_ptr(i);

      if (node &&
          (serial_on_zero || node->processor_id() == system.processor_id()) &&
          node->n_comp(system.number(), var_num) > 0)
        {
          dof_id_type dof_index = node->dof_number(system.number(), var_num, 0);

          // If the dof_index is local to this processor, set the value
          system.solution->set (dof_index, p.second);
        }
    }

  system.solution->close();
  system.update();
}



void ExodusII_IO::copy_elemental_solution(System & system,
                                          std::string system_var_name,
                                          std::string exodus_var_name,
                                          unsigned int timestep)
{
  LOG_SCOPE("copy_elemental_solution()", "ExodusII_IO");

  const unsigned int var_num = system.variable_number(system_var_name);
  // Assert that variable is an elemental one.
  //
  // NOTE: Currently, this reader is capable of reading only individual components of MONOMIAL_VEC
  //       types, and each must be written out to its own CONSTANT MONOMIAL variable
  libmesh_error_msg_if((system.variable_type(var_num) != FEType(CONSTANT, MONOMIAL))
                       && (system.variable_type(var_num) != FEType(CONSTANT, MONOMIAL_VEC)),
                       "Error! Trying to copy elemental solution into a variable that is not of CONSTANT MONOMIAL nor CONSTANT MONOMIAL_VEC type.");

  const MeshBase & mesh = MeshInput<MeshBase>::mesh();
  const DofMap & dof_map = system.get_dof_map();

  libmesh_error_msg_if(mesh.allow_renumbering(),
                       "ERROR, elemental data cannot be loaded if the mesh may be renumbered!");

  // Map from element ID to elemental variable value.  We need to use
  // a map here rather than a vector (e.g. elem_var_values) since the
  // libmesh element numbering can contain "holes".  This is the case
  // if we are reading elemental var values from an adaptively refined
  // mesh that has not been sequentially renumbered.
  std::map<dof_id_type, Real> elem_var_value_map;

  // With Exodus files we only open them on processor 0, so that's the
  // where we have to do the data read too.
  if (system.comm().rank() == 0)
    {
      libmesh_error_msg_if(!exio_helper->opened_for_reading,
                           "ERROR, ExodusII file must be opened for reading before copying an elemental solution!");

      exio_helper->read_elemental_var_values(exodus_var_name, timestep, elem_var_value_map);
    }

  const bool serial_on_zero = mesh.is_serial_on_zero();

  // If our mesh isn't serial, then non-root processors need to
  // request the data for their parts of the mesh and insert it
  // themselves.
  if (!serial_on_zero)
    {
      std::unordered_map<processor_id_type, std::vector<dof_id_type>> elem_ids_to_request;
      if (this->processor_id() != 0)
        {
          std::vector<dof_id_type> elem_ids;
          for (auto & elem : mesh.active_local_element_ptr_range())
            elem_ids.push_back(elem->id());

          if (!elem_ids.empty())
            elem_ids_to_request[0] = std::move(elem_ids);
        }

      auto value_gather_functor =
        [& elem_var_value_map]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         std::vector<Real> & values)
        {
          const std::size_t query_size = ids.size();
          values.resize(query_size);
          for (std::size_t i=0; i != query_size; ++i)
            {
              if (const auto it = elem_var_value_map.find(ids[i]);
                  it != elem_var_value_map.end())
                {
                  values[i] = it->second;
                  elem_var_value_map.erase(it);
                }
              else
                values[i] = std::numeric_limits<Real>::quiet_NaN();
            }
        };

      auto value_action_functor =
        [& elem_var_value_map]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         const std::vector<Real> & values)
        {
          const std::size_t query_size = ids.size();
          for (std::size_t i=0; i != query_size; ++i)
            if (!libmesh_isnan(values[i]))
              elem_var_value_map[ids[i]] = values[i];
        };

      Real * value_ex = nullptr;
      Parallel::pull_parallel_vector_data
        (system.comm(), elem_ids_to_request, value_gather_functor,
         value_action_functor, value_ex);
    }

  std::map<dof_id_type, Real>::iterator
    it = elem_var_value_map.begin(),
    end = elem_var_value_map.end();

  // Everybody inserts the data they've received.  If we're
  // serial_on_zero then proc 0 inserts everybody's data and other
  // procs have empty map ranges.
  for (; it!=end; ++it)
    {
      const Elem * elem = mesh.query_elem_ptr(it->first);

      if (elem && elem->n_comp(system.number(), var_num) > 0)
        {
          dof_id_type dof_index = elem->dof_number(system.number(), var_num, 0);
          if (serial_on_zero || dof_map.local_index(dof_index ))
            system.solution->set (dof_index, it->second);
        }
    }

  system.solution->close();
  system.update();
}

void ExodusII_IO::copy_scalar_solution(System & system,
                                       std::vector<std::string> system_var_names,
                                       std::vector<std::string> exodus_var_names,
                                       unsigned int timestep)
{
  LOG_SCOPE("copy_scalar_solution()", "ExodusII_IO");

  libmesh_error_msg_if(!exio_helper->opened_for_reading,
                       "ERROR, ExodusII file must be opened for reading before copying a scalar solution!");

  libmesh_error_msg_if(system_var_names.size() != exodus_var_names.size(),
                       "ERROR, the number of system_var_names must match exodus_var_names.");

  std::vector<Real> values_from_exodus;
  read_global_variable(exodus_var_names, timestep, values_from_exodus);

#ifdef LIBMESH_HAVE_MPI
  if (this->n_processors() > 1)
  {
    const Parallel::MessageTag tag = this->comm().get_unique_tag(1);
    if (this->processor_id() == this->n_processors()-1)
      this->comm().receive(0, values_from_exodus, tag);
    if (this->processor_id() == 0)
      this->comm().send(this->n_processors()-1, values_from_exodus, tag);
  }
#endif

  if (system.processor_id() == (system.n_processors()-1))
  {
    const DofMap & dof_map = system.get_dof_map();

    for (auto i : index_range(system_var_names))
    {
      const unsigned int var_num = system.variable_scalar_number(system_var_names[i], 0);

      std::vector<dof_id_type> SCALAR_dofs;
      dof_map.SCALAR_dof_indices(SCALAR_dofs, var_num);

      system.solution->set (SCALAR_dofs[0], values_from_exodus[i]);
    }
  }

  system.solution->close();
  system.update();
}

void ExodusII_IO::read_elemental_variable(std::string elemental_var_name,
                                          unsigned int timestep,
                                          std::map<unsigned int, Real> & unique_id_to_value_map)
{
  LOG_SCOPE("read_elemental_variable()", "ExodusII_IO");

  // Note that this function MUST be called before renumbering
  std::map<dof_id_type, Real> elem_var_value_map;

  exio_helper->read_elemental_var_values(elemental_var_name, timestep, elem_var_value_map);
  for (auto & pr : elem_var_value_map)
    {
      const Elem * elem = MeshInput<MeshBase>::mesh().query_elem_ptr(pr.first);
      unique_id_to_value_map.emplace(elem->top_parent()->unique_id(), pr.second);
    }
}

void ExodusII_IO::read_global_variable(std::vector<std::string> global_var_names,
                                       unsigned int timestep,
                                       std::vector<Real> & global_values)
{
  LOG_SCOPE("read_global_variable()", "ExodusII_IO");

  std::size_t size = global_var_names.size();
  libmesh_error_msg_if(size == 0, "ERROR, empty list of global variables to read from the Exodus file.");

  // read the values for all global variables
  std::vector<Real> values_from_exodus;
  exio_helper->read_var_names(ExodusII_IO_Helper::GLOBAL);
  exio_helper->read_global_values(values_from_exodus, timestep);
  std::vector<std::string> global_var_names_exodus = exio_helper->global_var_names;

  if (values_from_exodus.size() == 0)
    return;   // This will happen in parallel on procs that are not 0

  global_values.clear();
  for (std::size_t i = 0; i != size; ++i)
    {
      // for each global variable in global_var_names, look the corresponding one in global_var_names_from_exodus
      // and fill global_values accordingly
      auto it = find(global_var_names_exodus.begin(), global_var_names_exodus.end(), global_var_names[i]);
      if (it != global_var_names_exodus.end())
        global_values.push_back(values_from_exodus[it - global_var_names_exodus.begin()]);
      else
        libmesh_error_msg("ERROR, Global variable " << global_var_names[i] << \
                          " not found in Exodus file.");
    }

}

void ExodusII_IO::write_element_data (const EquationSystems & es)
{
  LOG_SCOPE("write_element_data()", "ExodusII_IO");

  // Be sure the file has been opened for writing!
  libmesh_error_msg_if(MeshOutput<MeshBase>::mesh().processor_id() == 0 && !exio_helper->opened_for_writing,
                       "ERROR, ExodusII file must be initialized before outputting element variables.");

  // This function currently only works on serialized meshes. We rely
  // on having a reference to a non-const MeshBase object from our
  // MeshInput parent class to construct a MeshSerializer object,
  // similar to what is done in ExodusII_IO::write().  Note that
  // calling ExodusII_IO::write_timestep() followed by
  // ExodusII_IO::write_element_data() when the underlying Mesh is a
  // DistributedMesh will result in an unnecessary additional
  // serialization/re-parallelization step.
  // The "true" specifies that we only need the mesh serialized to processor 0
  MeshSerializer serialize(MeshInput<MeshBase>::mesh(), !MeshOutput<MeshBase>::_is_parallel_format, true);

  // To be (possibly) filled with a filtered list of variable names to output.
  std::vector<std::string> names;

  // If _output_variables is populated, only output the monomials which are
  // also in the _output_variables vector.
  if (_output_variables.size() > 0)
    {
      // Create a list of CONSTANT MONOMIAL variable names
      std::vector<std::string> monomials;
      FEType type(CONSTANT, MONOMIAL);
      es.build_variable_names(monomials, &type);

      // Now concatenate a list of CONSTANT MONOMIAL_VEC variable names
      type = FEType(CONSTANT, MONOMIAL_VEC);
      es.build_variable_names(monomials, &type);

      // Filter that list against the _output_variables list.  Note: if names is still empty after
      // all this filtering, all the monomial variables will be gathered
      for (const auto & var : monomials)
        if (std::find(_output_variables.begin(), _output_variables.end(), var) != _output_variables.end())
          names.push_back(var);
    }

  // If we pass in a list of names to "build_elemental_solution_vector()"
  // it'll filter the variables coming back.
  std::vector<Number> soln;
  es.build_elemental_solution_vector(soln, names);

  // Also, store the list of subdomains on which each variable is active
  std::vector<std::set<subdomain_id_type>> vars_active_subdomains;
  es.get_vars_active_subdomains(names, vars_active_subdomains);

  if (soln.empty()) // If there is nothing to write just return
    return;

  // The data must ultimately be written block by block.  This means that this data
  // must be sorted appropriately.
  if (MeshOutput<MeshBase>::mesh().processor_id())
    return;

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

  std::vector<std::string> complex_names =
    exio_helper->get_complex_names(names, _write_complex_abs);

  std::vector<std::set<subdomain_id_type>>
    complex_vars_active_subdomains =
    exio_helper->get_complex_vars_active_subdomains(vars_active_subdomains,
                                                    _write_complex_abs);
  exio_helper->initialize_element_variables(complex_names, complex_vars_active_subdomains);

  const std::vector<Real> complex_soln =
    complex_soln_components(soln, names.size(), _write_complex_abs);

  exio_helper->write_element_values(mesh, complex_soln, _timestep, complex_vars_active_subdomains);

#else
  exio_helper->initialize_element_variables(names, vars_active_subdomains);
  exio_helper->write_element_values(mesh, soln, _timestep, vars_active_subdomains);
#endif
}



void
ExodusII_IO::write_element_data_from_discontinuous_nodal_data
(const EquationSystems & es,
 const std::set<std::string> * system_names,
 const std::string & var_suffix)
{
  LOG_SCOPE("write_element_data_from_discontinuous_nodal_data()", "ExodusII_IO");

  // Be sure that some other function has already opened the file and prepared it
  // for writing. This is the same behavior as the write_element_data() function
  // which we are trying to mimic.
  libmesh_error_msg_if(MeshOutput<MeshBase>::mesh().processor_id() == 0 && !exio_helper->opened_for_writing,
                       "ERROR, ExodusII file must be initialized before outputting element variables.");

  // This function currently only works on serialized meshes.  The
  // "true" flag specifies that we only need the mesh serialized to
  // processor 0
  MeshSerializer serialize(MeshInput<MeshBase>::mesh(),
                           !MeshOutput<MeshBase>::_is_parallel_format,
                           true);

  // Note: in general we want to respect the contents of
  // _output_variables, only building a solution vector with values
  // from the requested variables. First build a list of all variable
  // names, then throw out ones that aren't in _output_variables, if
  // any.
  std::vector<std::string> var_names;
  es.build_variable_names (var_names, /*fetype=*/nullptr, system_names);

  // Get a subset of all variable names that are CONSTANT,
  // MONOMIALs. We treat those slightly differently since they can
  // truly only have a single value per Elem.
  //
  // Should the same apply here for CONSTANT MONOMIAL_VECs? [CW]
  // That is, get rid of 'const' on 'fe_type' and rerun:
  //    fe_type = FEType(CONSTANT, MONOMIAL_VEC);
  //    es.build_variable_names(monomial_var_names, &fe_type);
  // Then, es.find_variable_numbers() can be used without a type
  // (since we know for sure they're monomials) like:
  //    var_nums = es.find_variable_numbers(monomial_var_names)
  // for which the DOF indices for 'var_nums' have to be resolved
  // manually like in build_elemental_solution_vector()
  std::vector<std::string> monomial_var_names;
  const FEType fe_type(CONSTANT, MONOMIAL);
  es.build_variable_names(monomial_var_names, &fe_type);

  // Remove all names from var_names that are not in _output_variables.
  // Note: This approach avoids errors when the user provides invalid
  // variable names in _output_variables, as the code will not try to
  // write a variable that doesn't exist.
  if (!_output_variables.empty())
    {
      var_names.erase
        (std::remove_if
         (var_names.begin(),
          var_names.end(),
          [this](const std::string & name)
          {return !std::count(_output_variables.begin(),
                              _output_variables.end(),
                              name);}),
         var_names.end());

      // Also filter the monomial variable names.
      monomial_var_names.erase
        (std::remove_if
         (monomial_var_names.begin(),
          monomial_var_names.end(),
          [this](const std::string & name)
          {return !std::count(_output_variables.begin(),
                              _output_variables.end(),
                              name);}),
         monomial_var_names.end());
    }

  // Build a solution vector, limiting the results to the variables in
  // var_names and the Systems in system_names, and only computing values
  // at the vertices.
  std::vector<Number> v;
  es.build_discontinuous_solution_vector
    (v, system_names, &var_names, /*vertices_only=*/true);

  // Get active subdomains for each variable in var_names.
  std::vector<std::set<subdomain_id_type>> vars_active_subdomains;
  es.get_vars_active_subdomains(var_names, vars_active_subdomains);

  // Determine names of variables to write based on the number of
  // nodes/vertices the elements in different subdomains have.
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();
  std::map<subdomain_id_type, unsigned int> subdomain_id_to_vertices_per_elem;
  for (const auto & elem : mesh.active_element_ptr_range())
    {
      // Try to insert key/value pair into the map. If this returns
      // false, check the returned iterator's value to make sure it
      // matches. It shouldn't actually be possible for this to fail
      // (since if the Mesh was like this it would have already
      // failed) but it doesn't hurt to be on the safe side.
      auto pr2 = subdomain_id_to_vertices_per_elem.emplace
        (elem->subdomain_id(), elem->n_vertices());
      libmesh_error_msg_if(!pr2.second && pr2.first->second != elem->n_vertices(),
                           "Elem with different number of vertices found.");
    }

  // Determine "derived" variable names. These names are created by
  // starting with the base variable name and appending the user's
  // variable_suffix (default: "_elem_node_") followed by a node id.
  //
  // Not every derived variable will be active on every subdomain,
  // even if the original variable _is_ active.  Subdomains can have
  // different geometric element types (with differing numbers of
  // nodes), so some of the derived variable names will be inactive on
  // those subdomains.
  //
  // Since we would otherwise generate the same name once per
  // subdomain, we keep the list of names unique as we are creating
  // it. We can't use a std::set for this because we don't want the
  // variables names to be in a different order from the order
  // they were written in the call to: build_discontinuous_solution_vector()
  //
  // The list of derived variable names includes one for each vertex,
  // for higher-order elements we currently only write out vertex
  // values, but this could be changed in the future without too much
  // trouble.
  std::vector<std::string> derived_var_names;

  // Keep track of mapping from derived_name to (orig_name, node_id)
  // pair. We will use this later to determine whether a given
  // variable is active on a given subdomain.
  std::map<std::string, std::pair<std::string, unsigned int>>
    derived_name_to_orig_name_and_node_id;

  for (const auto & pr : subdomain_id_to_vertices_per_elem)
    {
      const subdomain_id_type sbd_id = pr.first;
      const unsigned int vertices_per_elem =
        subdomain_id_to_vertices_per_elem[sbd_id];

      std::ostringstream oss;
      for (unsigned int n=0; n<vertices_per_elem; ++n)
        for (const auto & orig_var_name : var_names)
          {
            oss.str("");
            oss.clear();
            oss << orig_var_name << var_suffix << n;
            std::string derived_name = oss.str();

            // Only add this var name if it's not already in the list.
            if (!std::count(derived_var_names.begin(), derived_var_names.end(), derived_name))
              {
                derived_var_names.push_back(derived_name);
                // Add entry for derived_name -> (orig_name, node_id) mapping.
                derived_name_to_orig_name_and_node_id[derived_name] =
                  std::make_pair(orig_var_name, n);
              }
          }
    }

  // For each derived variable name, determine whether it is active
  // based on how many nodes/vertices the elements in a given subdomain have,
  // and whether they were active on the subdomain to begin with.
  std::vector<std::set<subdomain_id_type>>
    derived_vars_active_subdomains(derived_var_names.size());

  // A new data structure for keeping track of a list of variable names
  // that are in the discontinuous solution vector on each subdomain. Used
  // for indexing. Note: if a variable was inactive at the System level,
  // an entry for it will still be in the discontinuous solution vector,
  // but it will just have a value of zero. On the other hand, when we
  // create the derived variable names some of them are "inactive" on
  // different subdomains in the sense that they don't exist at all, i.e.
  // there is no zero padding for them. We need to be able to distinguish
  // between these two types in order to do the indexing into this vector
  // correctly.
  std::map<subdomain_id_type, std::vector<std::string>>
    subdomain_to_var_names;

  for (auto derived_var_id : index_range(derived_var_names))
    {
      const auto & derived_name = derived_var_names[derived_var_id];
      const auto & [orig_name, node_id] =
        libmesh_map_find (derived_name_to_orig_name_and_node_id,
                  derived_name);

      // For each subdomain, determine whether the current variable
      // should be active on that subdomain.
      for (const auto & pr : subdomain_id_to_vertices_per_elem)
        {
          // Convenience variables for the current subdomain and the
          // number of nodes elements in this subdomain have.
          subdomain_id_type sbd_id = pr.first;
          unsigned int vertices_per_elem_this_sbd =
            subdomain_id_to_vertices_per_elem[sbd_id];

          // Check whether variable orig_name was active on this
          // subdomain to begin with by looking in the
          // vars_active_subdomains container. We assume that the
          // location of orig_name in the var_names vector matches its
          // index in the vars_active_subdomains container.
          auto var_loc = std::find(var_names.begin(), var_names.end(), orig_name);
          libmesh_error_msg_if(var_loc == var_names.end(),
                               "Variable " << orig_name << " somehow not found in var_names array.");
          auto var_id = std::distance(var_names.begin(), var_loc);

          // The derived_var will only be active if this subdomain has
          // enough vertices for that to be the case.
          if (node_id < vertices_per_elem_this_sbd)
            {
              // Regardless of whether the original variable was not active on this subdomain,
              // the discontinuous solution vector will have zero padding for it, and
              // we will need to account for it. Therefore it should still be added to
              // the subdomain_to_var_names data structure!
              subdomain_to_var_names[sbd_id].push_back(derived_name);

              // If the original variable was not active on the
              // current subdomain, it should not be added to the
              // derived_vars_active_subdomains data structure, since
              // it will not be written to the Exodus file.

              // Determine if the original variable was active on the
              // current subdomain.
              bool orig_var_active =
                (vars_active_subdomains[var_id].empty() ||
                 vars_active_subdomains[var_id].count(sbd_id));

              // And only if it was, add it to the
              // derived_vars_active_subdomains data structure.
              if (orig_var_active)
                derived_vars_active_subdomains[derived_var_id].insert(sbd_id);
            }
        } // end loop over subdomain_id_to_vertices_per_elem
    } // end loop over derived_var_names

  // At this point we've built the "true" list of derived names, but
  // if there are any CONSTANT MONOMIALS in this list, we now want to
  // remove all but one copy of them from the derived_var_names list,
  // and rename them in (but not remove them from) the
  // subdomain_to_var_names list, and then update the
  // derived_vars_active_subdomains containers before finally calling
  // the Exodus helper functions.
  for (auto & derived_var_name : derived_var_names)
    {
      // Get the original name associated with this derived name.
      const auto & name_and_id =
        libmesh_map_find (derived_name_to_orig_name_and_node_id,
                  derived_var_name);

      // Convenience variables for the map entry's contents.
      const std::string & orig_name = name_and_id.first;

      // Was the original name a constant monomial?
      if (std::count(monomial_var_names.begin(),
                     monomial_var_names.end(),
                     orig_name))
        {
          // Rename this variable in the subdomain_to_var_names vectors.
          for (auto & pr : subdomain_to_var_names)
            {
              // Reference to ordered list of variable names on this subdomain.
              auto & name_vec = pr.second;

              auto name_vec_it =
                std::find(name_vec.begin(),
                          name_vec.end(),
                          derived_var_name);

              if (name_vec_it != name_vec.end())
                {
                  // Actually rename it back to the orig_name, dropping
                  // the "_elem_corner_" stuff.
                  *name_vec_it = orig_name;
                }
            }

          // Finally, rename the variable in the derived_var_names vector itself.
          derived_var_name = orig_name;
        } // if (monomial)
    } // end loop over derived names

  // Now remove duplicate entries from derived_var_names after the first.
  // Also update the derived_vars_active_subdomains container in a consistent way.
  {
    std::vector<std::string> derived_var_names_edited;
    std::vector<std::set<subdomain_id_type>> derived_vars_active_subdomains_edited;
    std::vector<unsigned int> found_first(monomial_var_names.size());

    for (auto i : index_range(derived_var_names))
      {
        const auto & derived_var_name = derived_var_names[i];
        const auto & active_set = derived_vars_active_subdomains[i];

        // Determine whether we will keep this derived variable name in
        // the final container.
        bool keep = true;
        for (auto j : index_range(monomial_var_names))
          if (derived_var_name == monomial_var_names[j])
            {
              if (!found_first[j])
                found_first[j] = 1;

              else
                keep = false;
            }

        // We also don't keep variables that are not active on any subdomains.
        // Contrary to other uses of the var_active_subdomains container where
        // the empty set means "all" subdomains, here it really means "none".
        if (active_set.empty())
          keep = false;

        if (keep)
          {
            derived_var_names_edited.push_back(derived_var_name);
            derived_vars_active_subdomains_edited.push_back(active_set);
          }
      }

    // We built the filtered ranges, now swap them with the originals.
    derived_var_names.swap(derived_var_names_edited);
    derived_vars_active_subdomains.swap(derived_vars_active_subdomains_edited);
  }

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  // Build complex variable names "r_foo", "i_foo", "a_foo" and the lists of
  // subdomains on which they are active.
  auto complex_var_names =
    exio_helper->get_complex_names(derived_var_names,
                                   _write_complex_abs);
  auto complex_vars_active_subdomains =
    exio_helper->get_complex_vars_active_subdomains(derived_vars_active_subdomains,
                                                    _write_complex_abs);
  auto complex_subdomain_to_var_names =
    exio_helper->get_complex_subdomain_to_var_names(subdomain_to_var_names,
                                                    _write_complex_abs);

  // Make expanded version of vector "v" in which each entry in the
  // original expands to an ("r_", "i_", "a_") triple.
  // "nco" is the number of complex outputs, which depends on whether
  // or not we are writing out the complex magnitudes.
  std::vector<Real> complex_v;
  int nco = _write_complex_abs ? 3 : 2;
  complex_v.reserve(nco * v.size());
  for (const auto & val : v)
    {
      complex_v.push_back(val.real());
      complex_v.push_back(val.imag());
      if (_write_complex_abs)
        complex_v.push_back(std::abs(val));
    }

  // Finally, initialize storage for the variables and write them to file.
  exio_helper->initialize_element_variables
    (complex_var_names, complex_vars_active_subdomains);
  exio_helper->write_element_values_element_major
    (mesh, complex_v, _timestep,
     complex_vars_active_subdomains,
     complex_var_names,
     complex_subdomain_to_var_names);
#else

  // Call function which writes the derived variable names to the
  // Exodus file.
  exio_helper->initialize_element_variables(derived_var_names, derived_vars_active_subdomains);

  // ES::build_discontinuous_solution_vector() creates a vector with
  // an element-major ordering, so call Helper::write_element_values()
  // passing false for the last argument.
  exio_helper->write_element_values_element_major
    (mesh, v, _timestep,
     derived_vars_active_subdomains,
     derived_var_names,
     subdomain_to_var_names);
#endif
}



void ExodusII_IO::write_nodal_data (const std::string & fname,
                                    const std::vector<Number> & soln,
                                    const std::vector<std::string> & names)
{
  LOG_SCOPE("write_nodal_data()", "ExodusII_IO");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  int num_vars = cast_int<int>(names.size());
  dof_id_type num_nodes = mesh.n_nodes();

  // The names of the variables to be output
  std::vector<std::string> output_names;

  if (_allow_empty_variables || !_output_variables.empty())
    output_names = _output_variables;
  else
    output_names = names;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  std::vector<std::string> complex_names =
    exio_helper->get_complex_names(output_names,
                                   _write_complex_abs);

  // Call helper function for opening/initializing data, giving it the
  // complex variable names
  this->write_nodal_data_common(fname, complex_names, /*continuous=*/true);
#else
  // Call helper function for opening/initializing data
  this->write_nodal_data_common(fname, output_names, /*continuous=*/true);
#endif

  if (mesh.processor_id())
    return;

  // This will count the number of variables actually output
  for (int c=0; c<num_vars; c++)
    {
      std::stringstream name_to_find;

      std::vector<std::string>::iterator pos =
        std::find(output_names.begin(), output_names.end(), names[c]);
      if (pos == output_names.end())
        continue;

      unsigned int variable_name_position =
        cast_int<unsigned int>(pos - output_names.begin());

      // Set up temporary vectors to be passed to Exodus to write the
      // nodal values for a single variable at a time.
#ifdef LIBMESH_USE_REAL_NUMBERS
      std::vector<Number> cur_soln;

      // num_nodes is either exactly how much space we will need for
      // each vector, or a safe upper bound for the amount of memory
      // we will require when there are gaps in the numbering.
      cur_soln.reserve(num_nodes);
#else
      std::vector<Real> real_parts;
      std::vector<Real> imag_parts;
      std::vector<Real> magnitudes;
      real_parts.reserve(num_nodes);
      imag_parts.reserve(num_nodes);
      if (_write_complex_abs)
        magnitudes.reserve(num_nodes);
#endif

      // There could be gaps in soln based on node numbering, but in
      // serial the empty numbers are left empty.
      // There could also be offsets in soln based on "fake" nodes
      // inserted on each processor (because NumericVector indices
      // have to be contiguous); the helper keeps track of those.
      // We now copy the proper solution values contiguously into
      // "cur_soln", removing the gaps.
      for (const auto & node : mesh.node_ptr_range())
        {
          const dof_id_type idx =
            (exio_helper->node_id_to_vec_id(node->id()))
            * num_vars + c;
#ifdef LIBMESH_USE_REAL_NUMBERS
          cur_soln.push_back(soln[idx]);
#else
          real_parts.push_back(soln[idx].real());
          imag_parts.push_back(soln[idx].imag());
          if (_write_complex_abs)
            magnitudes.push_back(std::abs(soln[idx]));
#endif
        }

      // If we're adding extra sides, we need to add their data too.
      //
      // Because soln was created from a parallel NumericVector, its
      // numbering was contiguous on each processor; we need to use
      // the same offsets here, and we need to loop through elements
      // from earlier ranks first.
      if (exio_helper->get_add_sides())
        {
          std::vector<std::vector<const Elem *>>
            elems_by_pid(mesh.n_processors());

          for (const auto & elem : mesh.active_element_ptr_range())
            elems_by_pid[elem->processor_id()].push_back(elem);

          for (auto p : index_range(elems_by_pid))
            {
              dof_id_type global_idx =
                exio_helper->added_node_offset_on(p) * num_vars + c;
              for (const Elem * elem : elems_by_pid[p])
                {
                  for (auto s : elem->side_index_range())
                    {
                      if (EquationSystems::redundant_added_side(*elem,s))
                        continue;

                      const std::vector<unsigned int> side_nodes =
                        elem->nodes_on_side(s);

                      for (auto n : index_range(side_nodes))
                        {
                          libmesh_ignore(n);
                          libmesh_assert_less(global_idx, soln.size());
#ifdef LIBMESH_USE_REAL_NUMBERS
                          cur_soln.push_back(soln[global_idx]);
#else
                          real_parts.push_back(soln[global_idx].real());
                          imag_parts.push_back(soln[global_idx].imag());
                          if (_write_complex_abs)
                            magnitudes.push_back(std::abs(soln[global_idx]));
#endif
                          global_idx += num_vars;
                        }
                    }
                }
            }
        }

      // Finally, actually call the Exodus API to write to file.
#ifdef LIBMESH_USE_REAL_NUMBERS
      exio_helper->write_nodal_values(variable_name_position+1, cur_soln, _timestep);
#else
      int nco = _write_complex_abs ? 3 : 2;
      exio_helper->write_nodal_values(nco*variable_name_position+1, real_parts, _timestep);
      exio_helper->write_nodal_values(nco*variable_name_position+2, imag_parts, _timestep);
      if (_write_complex_abs)
        exio_helper->write_nodal_values(3*variable_name_position+3, magnitudes, _timestep);
#endif

    }
}




void ExodusII_IO::write_information_records (const std::vector<std::string> & records)
{
  if (MeshOutput<MeshBase>::mesh().processor_id())
    return;

  libmesh_error_msg_if(!exio_helper->opened_for_writing,
                       "ERROR, ExodusII file must be initialized before outputting information records.");

  exio_helper->write_information_records(records);
}



void ExodusII_IO::write_global_data (const std::vector<Number> & soln,
                                     const std::vector<std::string> & names)
{
  LOG_SCOPE("write_global_data()", "ExodusII_IO");

  if (MeshOutput<MeshBase>::mesh().processor_id())
    return;

  libmesh_error_msg_if(!exio_helper->opened_for_writing,
                       "ERROR, ExodusII file must be initialized before outputting global variables.");

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

  std::vector<std::string> complex_names =
    exio_helper->get_complex_names(names,
                                   _write_complex_abs);

  exio_helper->initialize_global_variables(complex_names);

  const std::vector<Real> complex_soln =
    complex_soln_components(soln, names.size(), _write_complex_abs);

  exio_helper->write_global_values(complex_soln, _timestep);

#else
  exio_helper->initialize_global_variables(names);
  exio_helper->write_global_values(soln, _timestep);
#endif
}



void ExodusII_IO::write_timestep (const std::string & fname,
                                  const EquationSystems & es,
                                  const int timestep,
                                  const Real time,
                                  const std::set<std::string> * system_names)
{
  _timestep = timestep;
  MeshOutput<MeshBase>::write_equation_systems(fname,es,system_names);

  if (MeshOutput<MeshBase>::mesh().processor_id())
    return;

  exio_helper->write_timestep(timestep, time);
}


void ExodusII_IO::write_equation_systems (const std::string & fname,
                                          const EquationSystems & es,
                                          const std::set<std::string> * system_names)
{
  write_timestep(fname, es, 1, 0, system_names);
}


void ExodusII_IO::write_elemsets()
{
  libmesh_error_msg_if(!exio_helper->opened_for_writing,
                       "ERROR, ExodusII file must be opened for writing "
                       "before calling ExodusII_IO::write_elemsets()!");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();
  exio_helper->write_elemsets(mesh);
}

void
ExodusII_IO::
write_sideset_data(int timestep,
                   const std::vector<std::string> & var_names,
                   const std::vector<std::set<boundary_id_type>> & side_ids,
                   const std::vector<std::map<BoundaryInfo::BCTuple, Real>> & bc_vals)
{
  libmesh_error_msg_if(!exio_helper->opened_for_writing,
                       "ERROR, ExodusII file must be opened for writing "
                       "before calling ExodusII_IO::write_sideset_data()!");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();
  exio_helper->write_sideset_data(mesh, timestep, var_names, side_ids, bc_vals);
}



void
ExodusII_IO::
read_sideset_data(int timestep,
                  std::vector<std::string> & var_names,
                  std::vector<std::set<boundary_id_type>> & side_ids,
                  std::vector<std::map<BoundaryInfo::BCTuple, Real>> & bc_vals)
{
  libmesh_error_msg_if(!exio_helper->opened_for_reading,
                       "ERROR, ExodusII file must be opened for reading "
                       "before calling ExodusII_IO::read_sideset_data()!");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();
  exio_helper->read_sideset_data(mesh, timestep, var_names, side_ids, bc_vals);
}



void
ExodusII_IO::
get_sideset_data_indices (std::map<BoundaryInfo::BCTuple, unsigned int> & bc_array_indices)

{
  libmesh_error_msg_if(!exio_helper->opened_for_reading,
                       "ERROR, ExodusII file must be opened for reading "
                       "before calling ExodusII_IO::get_sideset_data_indices()!");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();
  exio_helper->get_sideset_data_indices(mesh, bc_array_indices);
}

void
ExodusII_IO::
get_nodeset_data_indices (std::map<BoundaryInfo::NodeBCTuple, unsigned int> & bc_array_indices)
{
  libmesh_error_msg_if(!exio_helper->opened_for_reading,
                       "ERROR, ExodusII file must be opened for reading "
                       "before calling ExodusII_IO::get_nodeset_data_indices()!");

  exio_helper->get_nodeset_data_indices(bc_array_indices);
}

void
ExodusII_IO::
write_nodeset_data (int timestep,
                    const std::vector<std::string> & var_names,
                    const std::vector<std::set<boundary_id_type>> & node_boundary_ids,
                    const std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> & bc_vals)
{
  libmesh_error_msg_if(!exio_helper->opened_for_writing,
                       "ERROR, ExodusII file must be opened for writing "
                       "before calling ExodusII_IO::write_nodeset_data()!");

  exio_helper->write_nodeset_data(timestep, var_names, node_boundary_ids, bc_vals);
}



void
ExodusII_IO::
read_nodeset_data (int timestep,
                   std::vector<std::string> & var_names,
                   std::vector<std::set<boundary_id_type>> & node_boundary_ids,
                   std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> & bc_vals)
{
  libmesh_error_msg_if(!exio_helper->opened_for_reading,
                       "ERROR, ExodusII file must be opened for reading "
                       "before calling ExodusII_IO::read_nodeset_data()!");

  exio_helper->read_nodeset_data(timestep, var_names, node_boundary_ids, bc_vals);
}

void
ExodusII_IO::
write_elemset_data (int timestep,
                    const std::vector<std::string> & var_names,
                    const std::vector<std::set<elemset_id_type>> & elemset_ids_in,
                    const std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> & elemset_vals)
{
  libmesh_error_msg_if(!exio_helper->opened_for_writing,
                       "ERROR, ExodusII file must be opened for writing "
                       "before calling ExodusII_IO::write_elemset_data()!");

  exio_helper->write_elemset_data(timestep, var_names, elemset_ids_in, elemset_vals);
}



void
ExodusII_IO::
read_elemset_data (int timestep,
                   std::vector<std::string> & var_names,
                   std::vector<std::set<elemset_id_type>> & elemset_ids_in,
                   std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> & elemset_vals)
{
  libmesh_error_msg_if(!exio_helper->opened_for_reading,
                       "ERROR, ExodusII file must be opened for reading "
                       "before calling ExodusII_IO::read_elemset_data()!");

  exio_helper->read_elemset_data(timestep, var_names, elemset_ids_in, elemset_vals);
}

void
ExodusII_IO::get_elemset_data_indices (std::map<std::pair<dof_id_type, elemset_id_type>, unsigned int> & elemset_array_indices)
{
  libmesh_error_msg_if(!exio_helper->opened_for_reading,
                       "ERROR, ExodusII file must be opened for reading "
                       "before calling ExodusII_IO::get_elemset_data_indices()!");

  exio_helper->get_elemset_data_indices(elemset_array_indices);
}


void ExodusII_IO::write (const std::string & fname)
{
  LOG_SCOPE("write()", "ExodusII_IO");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // We may need to gather a DistributedMesh to output it, making that
  // const qualifier in our constructor a dirty lie
  // The "true" specifies that we only need the mesh serialized to processor 0
  MeshSerializer serialize
    (const_cast<MeshBase &>(mesh),
     !MeshOutput<MeshBase>::_is_parallel_format, true);

  libmesh_assert( !exio_helper->opened_for_writing );

  // If the user has set the append flag here, it doesn't really make
  // sense: the intent of this function is to write a Mesh with no
  // data, while "appending" is really intended to add data to an
  // existing file.  If we're verbose, print a message to this effect.
  if (_append && _verbose)
    libmesh_warning("Warning: Appending in ExodusII_IO::write() does not make sense.\n"
                    "Creating a new file instead!");

  exio_helper->create(fname);
  exio_helper->initialize(fname,mesh);
  exio_helper->write_nodal_coordinates(mesh);
  exio_helper->write_elements(mesh);
  exio_helper->write_sidesets(mesh);
  exio_helper->write_nodesets(mesh);
  exio_helper->write_elemsets(mesh);

  if ((mesh.get_boundary_info().n_edge_conds() > 0) && _verbose)
    libmesh_warning("Warning: Mesh contains edge boundary IDs, but these "
                    "are not supported by the ExodusII format.");
}



void ExodusII_IO::write_nodal_data_discontinuous (const std::string & fname,
                                                  const std::vector<Number> & soln,
                                                  const std::vector<std::string> & names)
{
  LOG_SCOPE("write_nodal_data_discontinuous()", "ExodusII_IO");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

  std::vector<std::string> complex_names =
    exio_helper->get_complex_names(names,
                                   _write_complex_abs);

  // Call helper function for opening/initializing data, giving it the
  // complex variable names
  this->write_nodal_data_common(fname, complex_names, /*continuous=*/false);
#else
  // Call helper function for opening/initializing data
  this->write_nodal_data_common(fname, names, /*continuous=*/false);
#endif

  if (mesh.processor_id())
    return;

  int num_vars = cast_int<int>(names.size());
  libmesh_assert_equal_to(soln.size() % num_vars, 0);
  int num_nodes = soln.size() / num_vars;
  libmesh_assert_equal_to(exio_helper->num_nodes, num_nodes);

#ifndef NDEBUG
  if (!this->get_add_sides())
    {
      int num_real_nodes = 0;
      for (const auto & elem : mesh.active_element_ptr_range())
        num_real_nodes += elem->n_nodes();
      libmesh_assert_equal_to(num_real_nodes, num_nodes);
    }
#endif

  for (int c=0; c<num_vars; c++)
    {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      std::vector<Real> real_parts(num_nodes);
      std::vector<Real> imag_parts(num_nodes);
      std::vector<Real> magnitudes;
      if (_write_complex_abs)
        magnitudes.resize(num_nodes);

      // The number of complex outputs depends on whether or not we are
      // writing out the absolute values.
      int nco = _write_complex_abs ? 3 : 2;

      for (int i=0; i<num_nodes; ++i)
        {
          real_parts[i] = soln[i*num_vars + c].real();
          imag_parts[i] = soln[i*num_vars + c].imag();
          if (_write_complex_abs)
            magnitudes[i] = std::abs(soln[i*num_vars + c]);
        }
      exio_helper->write_nodal_values(nco*c+1, real_parts, _timestep);
      exio_helper->write_nodal_values(nco*c+2, imag_parts, _timestep);
      if (_write_complex_abs)
        exio_helper->write_nodal_values(3*c+3, magnitudes, _timestep);
#else
      // Copy out this variable's solution
      std::vector<Number> cur_soln(num_nodes);

      for (int i=0; i<num_nodes; i++)
        cur_soln[i] = soln[i*num_vars + c];

      exio_helper->write_nodal_values(c+1,cur_soln,_timestep);
#endif
    }
}



void ExodusII_IO::write_nodal_data_common(std::string fname,
                                          const std::vector<std::string> & names,
                                          bool continuous)
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // This function can be called multiple times, we only want to open
  // the ExodusII file the first time it's called.
  if (!exio_helper->opened_for_writing)
    {
      // If we're appending, open() the file with read_only=false,
      // otherwise create() it and write the contents of the mesh to
      // it.
      if (_append)
        {
          // We do our writing only from proc 0, to avoid race
          // conditions with Exodus 8
          if (!MeshOutput<MeshBase>::mesh().processor_id())
            {
              exio_helper->open(fname.c_str(), /*read_only=*/false);
              // If we're appending, it's not valid to call exio_helper->initialize()
              // or exio_helper->initialize_nodal_variables(), but we do need to set up
              // certain aspects of the Helper object itself, such as the number of nodes
              // and elements.  We do that by reading the header...
              exio_helper->read_and_store_header_info();

              // ...and reading the block info
              exio_helper->read_block_info();
            }
            // Keep other processors aware of what we've done on root
          else
            {
              exio_helper->opened_for_writing = true;
              exio_helper->current_filename = fname;
            }
        }
      else
        {
          exio_helper->create(fname);

          // But some of our write calls are parallel-only, due to
          // calls to parallel-only getter functions.
          exio_helper->initialize(fname, mesh, !continuous);

          exio_helper->write_nodal_coordinates(mesh, !continuous);
          exio_helper->write_elements(mesh, !continuous);

          exio_helper->write_sidesets(mesh);
          exio_helper->write_nodesets(mesh);
          exio_helper->write_elemsets(mesh);

          exio_helper->initialize_nodal_variables(names);
        }
    }
  else
    {
      // We are already open for writing, so check that the filename
      // passed to this function matches the filename currently in use
      // by the helper.
      libmesh_error_msg_if(fname != exio_helper->current_filename,
                           "Error! This ExodusII_IO object is already associated with file: "
                           << exio_helper->current_filename
                           << ", cannot use it with requested file: "
                           << fname);
    }
}

const std::vector<std::string> & ExodusII_IO::get_nodal_var_names()
{
  exio_helper->read_var_names(ExodusII_IO_Helper::NODAL);
  return exio_helper->nodal_var_names;
}

const std::vector<std::string> & ExodusII_IO::get_elem_var_names()
{
  exio_helper->read_var_names(ExodusII_IO_Helper::ELEMENTAL);
  return exio_helper->elem_var_names;
}

const std::vector<std::string> & ExodusII_IO::get_global_var_names()
{
  exio_helper->read_var_names(ExodusII_IO_Helper::GLOBAL);
  return exio_helper->global_var_names;
}

const std::vector<int> & ExodusII_IO::get_elem_num_map() const
{
  // We could make this function non-const and have it call
  // exio_helper->read_elem_num_map() before returning a reference...
  // but the intention is that this will be called some time after a
  // mesh is read in, in which case it would be doing extra work to
  // read in the elem_num_map twice.
  return exio_helper->elem_num_map;
}

const std::vector<int> & ExodusII_IO::get_node_num_map() const
{
  return exio_helper->node_num_map;
}

ExodusII_IO_Helper & ExodusII_IO::get_exio_helper()
{
  // Provide a warning when accessing the helper object
  // since it is a non-public API and is likely to see
  // future API changes
  libmesh_experimental();

  return *exio_helper;
}


void ExodusII_IO::set_hdf5_writing(bool write_hdf5)
{
  exio_helper->set_hdf5_writing(write_hdf5);
}


void ExodusII_IO::set_max_name_length(unsigned int max_length)
{
  exio_helper->set_max_name_length(max_length);
}


void ExodusII_IO::set_discontinuous_bex(bool disc_bex)
{
  _disc_bex = disc_bex;
}



// LIBMESH_HAVE_EXODUS_API is not defined, declare error() versions of functions...
#else



ExodusII_IO::~ExodusII_IO () = default;



void ExodusII_IO::read (const std::string &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



ExodusHeaderInfo ExodusII_IO::read_header (const std::string &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::verbose (bool)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::use_mesh_dimension_instead_of_spatial_dimension(bool)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::write_as_dimension(unsigned)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::set_coordinate_offset(Point)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::append(bool)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::write_added_sides (bool)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



bool ExodusII_IO::get_add_sides ()
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
  return false;
}



const std::vector<Real> & ExodusII_IO::get_time_steps()
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



int ExodusII_IO::get_num_time_steps()
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}


void ExodusII_IO::copy_nodal_solution(System &,
                                      std::string,
                                      std::string,
                                      unsigned int)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::copy_elemental_solution(System &,
                                          std::string,
                                          std::string,
                                          unsigned int)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::copy_scalar_solution(System &,
                                       std::vector<std::string>,
                                       std::vector<std::string>,
                                       unsigned int)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::write_element_data (const EquationSystems &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void
ExodusII_IO::write_element_data_from_discontinuous_nodal_data
(const EquationSystems &,
 const std::set<std::string> *,
 const std::string & )
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::write_nodal_data (const std::string &,
                                    const std::vector<Number> &,
                                    const std::vector<std::string> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::write_information_records (const std::vector<std::string> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::write_global_data (const std::vector<Number> &,
                                     const std::vector<std::string> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::write_timestep (const std::string &,
                                  const EquationSystems &,
                                  const int,
                                  const Real,
                                  const std::set<std::string> *)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}


void ExodusII_IO::write_equation_systems (const std::string &,
                                          const EquationSystems &,
                                          const std::set<std::string> *)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}


void ExodusII_IO::write_elemsets()
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void
ExodusII_IO::
write_sideset_data (int,
                    const std::vector<std::string> &,
                    const std::vector<std::set<boundary_id_type>> &,
                    const std::vector<std::map<BoundaryInfo::BCTuple, Real>> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void
ExodusII_IO::
read_sideset_data (int,
                   std::vector<std::string> &,
                   std::vector<std::set<boundary_id_type>> &,
                   std::vector<std::map<BoundaryInfo::BCTuple, Real>> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void
ExodusII_IO::
get_sideset_data_indices (std::map<BoundaryInfo::BCTuple, unsigned int> &)

{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void
ExodusII_IO::
write_nodeset_data (int,
                    const std::vector<std::string> &,
                    const std::vector<std::set<boundary_id_type>> &,
                    const std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void
ExodusII_IO::
read_nodeset_data (int,
                   std::vector<std::string> &,
                   std::vector<std::set<boundary_id_type>> &,
                   std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void
ExodusII_IO::
get_nodeset_data_indices (std::map<BoundaryInfo::NodeBCTuple, unsigned int> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void
ExodusII_IO::
write_elemset_data (int,
                    const std::vector<std::string> &,
                    const std::vector<std::set<elemset_id_type>> &,
                    const std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void
ExodusII_IO::
read_elemset_data (int,
                   std::vector<std::string> &,
                   std::vector<std::set<elemset_id_type>> &,
                   std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void
ExodusII_IO::
get_elemset_data_indices (std::map<std::pair<dof_id_type, elemset_id_type>, unsigned int> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void ExodusII_IO::write (const std::string &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::write_nodal_data_discontinuous (const std::string &,
                                                  const std::vector<Number> &,
                                                  const std::vector<std::string> &)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}



void ExodusII_IO::write_nodal_data_common(std::string,
                                          const std::vector<std::string> &,
                                          bool)
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}


const std::vector<std::string> & ExodusII_IO::get_elem_var_names()
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

const std::vector<std::string> & ExodusII_IO::get_nodal_var_names()
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

const std::vector<std::string> & ExodusII_IO::get_global_var_names()
{
  libmesh_error_msg("ERROR, ExodusII API is not defined.");
}

void ExodusII_IO::set_hdf5_writing(bool) {}

#endif // LIBMESH_HAVE_EXODUS_API
} // namespace libMesh
