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


// C++ includes
#include <fstream>
#include <cstring>
#include <sstream>
#include <map>

// Local includes
#include "libmesh/exodusII_io.h"
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_base.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io_helper.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/dof_map.h"

namespace libMesh
{

// ------------------------------------------------------------
// ExodusII_IO class members
ExodusII_IO::ExodusII_IO (MeshBase & mesh,
#ifdef LIBMESH_HAVE_EXODUS_API
                          bool single_precision
#else
                          bool
#endif
                          ) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase> (mesh,
                        /* is_parallel_format = */ false,
                        /* serial_only_needed_on_proc_0 = */ true),
  ParallelObject(mesh),
#ifdef LIBMESH_HAVE_EXODUS_API
  exio_helper(new ExodusII_IO_Helper(*this, false, true, single_precision)),
  _timestep(1),
  _verbose(false),
  _append(false),
#endif
  _allow_empty_variables(false)
{
}


void ExodusII_IO::set_output_variables(const std::vector<std::string> & output_variables,
                                       bool allow_empty)
{
  _output_variables = output_variables;
  _allow_empty_variables = allow_empty;
}



#ifdef LIBMESH_ENABLE_DEPRECATED
void ExodusII_IO::copy_nodal_solution(System & system,
                                      std::string var_name,
                                      unsigned int timestep)
{
  libmesh_deprecated();
  copy_nodal_solution(system, var_name, var_name, timestep);
}
#endif



void ExodusII_IO::write_discontinuous_exodusII(const std::string & name,
                                               const EquationSystems & es,
                                               const std::set<std::string> * system_names)
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;

  es.build_variable_names  (solution_names, nullptr, system_names);
  es.build_discontinuous_solution_vector (v, system_names);
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
  // Get a reference to the mesh we are reading
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Clear any existing mesh data
  mesh.clear();

  // Keep track of what kinds of elements this file contains
  elems_of_dimension.clear();
  elems_of_dimension.resize(4, false);

  // Instantiate the ElementMaps interface
  ExodusII_IO_Helper::ElementMaps em;

  // Open the exodus file in EX_READ mode
  exio_helper->open(fname.c_str(), /*read_only=*/true);

  // Get header information from exodus file
  exio_helper->read_header();

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

  // Loop over the nodes, create Nodes with local processor_id 0.
  for (int i=0; i<exio_helper->num_nodes; i++)
    {
      // Use the node_num_map to get the correct ID for Exodus
      int exodus_id = exio_helper->node_num_map[i];

      // Catch the node that was added to the mesh
      Node * added_node = mesh.add_point (Point(exio_helper->x[i], exio_helper->y[i], exio_helper->z[i]), exodus_id-1);

      // If the Mesh assigned an ID different from what is in the
      // Exodus file, we should probably error.
      if (added_node->id() != static_cast<unsigned>(exodus_id-1))
        libmesh_error_msg("Error!  Mesh assigned node ID "    \
                          << added_node->id()                         \
                          << " which is different from the (zero-based) Exodus ID " \
                          << exodus_id-1                              \
                          << "!");
    }

  // This assert is no longer valid if the nodes are not numbered
  // sequentially starting from 1 in the Exodus file.
  // libmesh_assert_equal_to (static_cast<unsigned int>(exio_helper->num_nodes), mesh.n_nodes());

  // Get information about all the blocks
  exio_helper->read_block_info();

  // Reserve space for the elements
  mesh.reserve_elem(exio_helper->num_elem);

  // Read the element number map from the Exodus file.  This is
  // required if we want to preserve the numbering of elements as it
  // exists in the Exodus file.  If the Exodus file does not contain
  // an elem_num_map, the identity map is returned by this call.
  exio_helper->read_elem_num_map();

  // Read in the element connectivity for each block.
  int nelem_last_block = 0;

  // Loop over all the blocks
  for (int i=0; i<exio_helper->num_elem_blk; i++)
    {
      // Read the information for block i
      exio_helper->read_elem_in_block (i);
      int subdomain_id = exio_helper->get_block_id(i);

      // populate the map of names
      std::string subdomain_name = exio_helper->get_block_name(i);
      if (!subdomain_name.empty())
        mesh.subdomain_name(static_cast<subdomain_id_type>(subdomain_id)) = subdomain_name;

      // Set any relevant node/edge maps for this element
      const std::string type_str (exio_helper->get_elem_type());
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(type_str);

      // Loop over all the faces in this block
      int jmax = nelem_last_block+exio_helper->num_elem_this_blk;
      for (int j=nelem_last_block; j<jmax; j++)
        {
          Elem * elem = Elem::build (conv.get_canonical_type()).release();
          libmesh_assert (elem);
          elem->subdomain_id() = static_cast<subdomain_id_type>(subdomain_id) ;

          // Use the elem_num_map to obtain the ID of this element in the Exodus file
          int exodus_id = exio_helper->elem_num_map[j];

          // Assign this element the same ID it had in the Exodus
          // file, but make it zero-based by subtracting 1.  Note:
          // some day we could use 1-based numbering in libmesh and
          // thus match the Exodus numbering exactly, but at the
          // moment libmesh is zero-based.
          elem->set_id(exodus_id-1);

          // Record that we have seen an element of dimension elem->dim()
          elems_of_dimension[elem->dim()] = true;

          // Catch the Elem pointer that the Mesh throws back
          elem = mesh.add_elem (elem);

          // If the Mesh assigned an ID different from what is in the
          // Exodus file, we should probably error.
          if (elem->id() != static_cast<unsigned>(exodus_id-1))
            libmesh_error_msg("Error!  Mesh assigned ID "       \
                              << elem->id()                             \
                              << " which is different from the (zero-based) Exodus ID " \
                              << exodus_id-1                            \
                              << "!");

          // Set all the nodes for this element
          for (int k=0; k<exio_helper->num_nodes_per_elem; k++)
            {
              // global index
              int gi = (j-nelem_last_block)*exio_helper->num_nodes_per_elem + conv.get_node_map(k);

              // The entries in 'connect' are actually (1-based)
              // indices into the node_num_map, so to get the right
              // node ID we:
              // 1.) Subtract 1 from connect[gi]
              // 2.) Pass it through node_num_map to get the corresponding Exodus ID
              // 3.) Subtract 1 from that, since libmesh node numbering is "zero"-based,
              //     even when the Exodus node numbering doesn't start with 1.
              int libmesh_node_id = exio_helper->node_num_map[exio_helper->connect[gi] - 1] - 1;

              // Set the node pointer in the Elem
              elem->set_node(k) = mesh.node_ptr(libmesh_node_id);
            }
        }

      // running sum of # of elements per block,
      // (should equal total number of elements in the end)
      nelem_last_block += exio_helper->num_elem_this_blk;
    }

  // This assert isn't valid if the Exodus file's numbering doesn't
  // start with 1!  For example, if Exodus's elem_num_map is 21, 22,
  // 23, 24, 25, 26, 27, 28, 29, 30, ... 84, then by the time you are
  // done with the loop above, mesh.n_elem() will report 84 and
  // nelem_last_block will be 64.
  // libmesh_assert_equal_to (static_cast<unsigned>(nelem_last_block), mesh.n_elem());

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
        // The numbers in the Exodus file sidesets should be thought
        // of as (1-based) indices into the elem_num_map array.  So,
        // to get the right element ID we have to:
        // 1.) Subtract 1 from elem_list[e] (to get a zero-based index)
        // 2.) Pass it through elem_num_map (to get the corresponding Exodus ID)
        // 3.) Subtract 1 from that, since libmesh is "zero"-based,
        //     even when the Exodus numbering doesn't start with 1.
        dof_id_type libmesh_elem_id =
          cast_int<dof_id_type>(exio_helper->elem_num_map[exio_helper->elem_list[e] - 1] - 1);

        // Set any relevant node/edge maps for this element
        Elem & elem = mesh.elem_ref(libmesh_elem_id);

        const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(elem.type());

        // Map the zero-based Exodus side numbering to the libmesh side numbering
        unsigned int raw_side_index = exio_helper->side_list[e]-1;
        std::size_t side_index_offset = conv.get_shellface_index_offset();

        if (raw_side_index < side_index_offset)
          {
            // We assume this is a "shell face"
            int mapped_shellface = raw_side_index;

            // Check for errors
            if (mapped_shellface == ExodusII_IO_Helper::Conversion::invalid_id)
              libmesh_error_msg("Invalid 1-based side id: "                 \
                                << mapped_shellface                         \
                                << " detected for "                         \
                                << Utility::enum_to_string(elem.type()));

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
            if (mapped_side == ExodusII_IO_Helper::Conversion::invalid_id)
              libmesh_error_msg("Invalid 1-based side id: "                 \
                                << side_index                               \
                                << " detected for "                         \
                                << Utility::enum_to_string(elem.type()));

            // Add this (elem,side,id) triplet to the BoundaryInfo object.
            mesh.get_boundary_info().add_side (libmesh_elem_id,
                                               cast_int<unsigned short>(mapped_side),
                                               cast_int<boundary_id_type>(exio_helper->id_list[e]));
          }
      }
  }

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
            int exodus_id = exio_helper->node_sets_node_list[i + offset];

            // As before, the entries in 'node_list' are 1-based
            // indices into the node_num_map array, so we have to map
            // them.  See comment above.
            int libmesh_node_id = exio_helper->node_num_map[exodus_id - 1] - 1;
            mesh.get_boundary_info().add_node(cast_int<dof_id_type>(libmesh_node_id),
                                              nodeset_id);
          }
      }
  }

#if LIBMESH_DIM < 3
  if (mesh.mesh_dimension() > LIBMESH_DIM)
    libmesh_error_msg("Cannot open dimension "        \
                      << mesh.mesh_dimension()            \
                      << " mesh file when configured without "        \
                      << mesh.mesh_dimension()                        \
                      << "D support.");
#endif
}



void ExodusII_IO::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;

  // Set the verbose flag in the helper object as well.
  exio_helper->verbose = _verbose;
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



const std::vector<Real> & ExodusII_IO::get_time_steps()
{
  if (!exio_helper->opened_for_reading)
    libmesh_error_msg("ERROR, ExodusII file must be opened for reading before calling ExodusII_IO::get_time_steps()!");

  exio_helper->read_time_steps();
  return exio_helper->time_steps;
}



int ExodusII_IO::get_num_time_steps()
{
  if (!exio_helper->opened_for_reading && !exio_helper->opened_for_writing)
    libmesh_error_msg("ERROR, ExodusII file must be opened for reading or writing before calling ExodusII_IO::get_num_time_steps()!");

  exio_helper->read_num_time_steps();
  return exio_helper->num_time_steps;
}



void ExodusII_IO::copy_nodal_solution(System & system,
                                      std::string system_var_name,
                                      std::string exodus_var_name,
                                      unsigned int timestep)
{
  if (!exio_helper->opened_for_reading)
    libmesh_error_msg("ERROR, ExodusII file must be opened for reading before copying a nodal solution!");

  exio_helper->read_nodal_var_values(exodus_var_name, timestep);

  const unsigned int var_num = system.variable_number(system_var_name);

  for (dof_id_type i=0,
       n_nodal = cast_int<dof_id_type>(exio_helper->nodal_var_values.size());
       i != n_nodal; ++i)
    {
      const Node * node = MeshInput<MeshBase>::mesh().query_node_ptr(i);

      if (node && node->n_comp(system.number(), var_num) > 0)
        {
          dof_id_type dof_index = node->dof_number(system.number(), var_num, 0);

          // If the dof_index is local to this processor, set the value
          if ((dof_index >= system.solution->first_local_index()) && (dof_index < system.solution->last_local_index()))
            system.solution->set (dof_index, exio_helper->nodal_var_values[i]);
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
  if (system.comm().rank() == 0)
    {
      if (!exio_helper->opened_for_reading)
        libmesh_error_msg("ERROR, ExodusII file must be opened for reading before copying an elemental solution!");

      // Map from element ID to elemental variable value.  We need to use
      // a map here rather than a vector (e.g. elem_var_values) since the
      // libmesh element numbering can contain "holes".  This is the case
      // if we are reading elemental var values from an adaptively refined
      // mesh that has not been sequentially renumbered.
      std::map<dof_id_type, Real> elem_var_value_map;
      exio_helper->read_elemental_var_values(exodus_var_name, timestep, elem_var_value_map);

      const unsigned int var_num = system.variable_number(system_var_name);
      if (system.variable_type(var_num) != FEType(CONSTANT, MONOMIAL))
        libmesh_error_msg("Error! Trying to copy elemental solution into a variable that is not of CONSTANT MONOMIAL type.");

      std::map<dof_id_type, Real>::iterator
        it = elem_var_value_map.begin(),
        end = elem_var_value_map.end();

      for (; it!=end; ++it)
        {
          const Elem * elem = MeshInput<MeshBase>::mesh().query_elem_ptr(it->first);

          if (elem && elem->n_comp(system.number(), var_num) > 0)
            {
              dof_id_type dof_index = elem->dof_number(system.number(), var_num, 0);
              system.solution->set (dof_index, it->second);
            }
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
  if (!exio_helper->opened_for_reading)
    libmesh_error_msg("ERROR, ExodusII file must be opened for reading before copying a scalar solution!");

  if (system_var_names.size() != exodus_var_names.size())
    libmesh_error_msg("ERROR, the number of system_var_names must match exodus_var_names.");

  std::vector<Real> values_from_exodus;
  read_global_variable(exodus_var_names, timestep, values_from_exodus);

#ifdef LIBMESH_HAVE_MPI
  if (this->n_processors() > 1)
  {
    const Parallel::MessageTag tag(1);
    if (this->processor_id() == this->n_processors()-1)
      this->comm().receive(0, values_from_exodus, tag);
    if (this->processor_id() == 0)
      this->comm().send(this->n_processors()-1, values_from_exodus, tag);
  }
#endif

  if (system.processor_id() == (system.n_processors()-1))
  {
    const DofMap & dof_map = system.get_dof_map();

    for (std::size_t i = 0; i < system_var_names.size(); i++)
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
  // Note that this function MUST be called before renumbering
  std::map<dof_id_type, Real> elem_var_value_map;

  exio_helper->read_elemental_var_values(elemental_var_name, timestep, elem_var_value_map);
  for (auto & pr : elem_var_value_map)
    {
      const Elem * elem = MeshInput<MeshBase>::mesh().query_elem_ptr(pr.first);
      unique_id_to_value_map.insert(std::make_pair(elem->top_parent()->unique_id(), pr.second));
    }
}

void ExodusII_IO::read_global_variable(std::vector<std::string> global_var_names,
                                       unsigned int timestep,
                                       std::vector<Real> & global_values)
{
  std::size_t size = global_var_names.size();
  if (size == 0)
    libmesh_error_msg("ERROR, empty list of global variables to read from the Exodus file.");

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
  // Be sure the file has been opened for writing!
  if (MeshOutput<MeshBase>::mesh().processor_id() == 0 && !exio_helper->opened_for_writing)
    libmesh_error_msg("ERROR, ExodusII file must be initialized before outputting element variables.");

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
      std::vector<std::string> monomials;
      const FEType type(CONSTANT, MONOMIAL);

      // Create a list of monomial variable names
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

  std::vector<std::string> complex_names = exio_helper->get_complex_names(names);

  std::vector<std::set<subdomain_id_type>> complex_vars_active_subdomains =
    exio_helper->get_complex_vars_active_subdomains(vars_active_subdomains);
  exio_helper->initialize_element_variables(complex_names, complex_vars_active_subdomains);

  unsigned int num_values = soln.size();
  unsigned int num_vars = names.size();
  unsigned int num_elems = num_values / num_vars;

  // This will contain the real and imaginary parts and the magnitude
  // of the values in soln
  std::vector<Real> complex_soln(3*num_values);

  for (unsigned i=0; i<num_vars; ++i)
    {

      for (unsigned int j=0; j<num_elems; ++j)
        {
          Number value = soln[i*num_vars + j];
          complex_soln[3*i*num_elems + j] = value.real();
        }
      for (unsigned int j=0; j<num_elems; ++j)
        {
          Number value = soln[i*num_vars + j];
          complex_soln[3*i*num_elems + num_elems +j] = value.imag();
        }
      for (unsigned int j=0; j<num_elems; ++j)
        {
          Number value = soln[i*num_vars + j];
          complex_soln[3*i*num_elems + 2*num_elems + j] = std::abs(value);
        }
    }

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
  // Be sure that some other function has already opened the file and prepared it
  // for writing. This is the same behavior as the write_element_data() function
  // which we are trying to mimic.
  if (MeshOutput<MeshBase>::mesh().processor_id() == 0 && !exio_helper->opened_for_writing)
    libmesh_error_msg("ERROR, ExodusII file must be initialized before outputting element variables.");

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
      auto pr2 = subdomain_id_to_vertices_per_elem.insert
        (std::make_pair(elem->subdomain_id(), elem->n_vertices()));
      if (!pr2.second && pr2.first->second != elem->n_vertices())
        libmesh_error_msg("Elem with different number of vertices found.");
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
  // that are in the discontinous solution vector on each subdomain. Used
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
      auto it = derived_name_to_orig_name_and_node_id.find(derived_name);
      if (it == derived_name_to_orig_name_and_node_id.end())
        libmesh_error_msg("Unmapped derived variable name: " << derived_name);

      // Convenience variables for the map entry's contents.
      const std::string & orig_name = it->second.first;
      const unsigned int node_id = it->second.second;

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
          if (var_loc == var_names.end())
            libmesh_error_msg("Variable " << orig_name << " somehow not found in var_names array.");
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
      auto it = derived_name_to_orig_name_and_node_id.find(derived_var_name);
      if (it == derived_name_to_orig_name_and_node_id.end())
        libmesh_error_msg("Unmapped derived variable name: " << derived_var_name);

      // Convenience variables for the map entry's contents.
      const std::string & orig_name = it->second.first;

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

  // Call function which writes the derived variable names to the
  // Exodus file.
  exio_helper->initialize_element_variables(derived_var_names, derived_vars_active_subdomains);

  // ES::build_discontinuous_solution_vector() creates a vector with
  // an element-major ordering, so call Helper::write_element_values()
  // passing false for the last argument.
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  // TODO: Make this work when complex variables are enabled by
  // writing the real and complex parts separately. See
  // ExodusII_IO::write_element_data() for details.
  libmesh_not_implemented_msg
    ("write_element_data_from_discontinuous_nodal_data() is not "
     "yet supported when complex variables are enabled.");
#else
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

  std::vector<std::string> complex_names = exio_helper->get_complex_names(names);

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
      magnitudes.reserve(num_nodes);
#endif

      // There could be gaps in "soln", but it will always be in the
      // order of [num_vars * node_id + var_id]. We now copy the
      // proper solution values contiguously into "cur_soln",
      // removing the gaps.
      for (const auto & node : mesh.node_ptr_range())
        {
          dof_id_type idx = node->id()*num_vars + c;
#ifdef LIBMESH_USE_REAL_NUMBERS
          cur_soln.push_back(soln[idx]);
#else
          real_parts.push_back(soln[idx].real());
          imag_parts.push_back(soln[idx].imag());
          magnitudes.push_back(std::abs(soln[idx]));
#endif
        }

      // Finally, actually call the Exodus API to write to file.
#ifdef LIBMESH_USE_REAL_NUMBERS
      exio_helper->write_nodal_values(variable_name_position+1, cur_soln, _timestep);
#else
      exio_helper->write_nodal_values(3*variable_name_position+1, real_parts, _timestep);
      exio_helper->write_nodal_values(3*variable_name_position+2, imag_parts, _timestep);
      exio_helper->write_nodal_values(3*variable_name_position+3, magnitudes, _timestep);
#endif

    }
}




void ExodusII_IO::write_information_records (const std::vector<std::string> & records)
{
  if (MeshOutput<MeshBase>::mesh().processor_id())
    return;

  if (!exio_helper->opened_for_writing)
    libmesh_error_msg("ERROR, ExodusII file must be initialized before outputting information records.");

  exio_helper->write_information_records(records);
}



void ExodusII_IO::write_global_data (const std::vector<Number> & soln,
                                     const std::vector<std::string> & names)
{
  if (MeshOutput<MeshBase>::mesh().processor_id())
    return;

  if (!exio_helper->opened_for_writing)
    libmesh_error_msg("ERROR, ExodusII file must be initialized before outputting global variables.");

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

  std::vector<std::string> complex_names = exio_helper->get_complex_names(names);

  exio_helper->initialize_global_variables(complex_names);

  unsigned int num_values = soln.size();
  unsigned int num_vars = names.size();
  unsigned int num_elems = num_values / num_vars;

  // This will contain the real and imaginary parts and the magnitude
  // of the values in soln
  std::vector<Real> complex_soln(3*num_values);

  for (unsigned i=0; i<num_vars; ++i)
    {

      for (unsigned int j=0; j<num_elems; ++j)
        {
          Number value = soln[i*num_vars + j];
          complex_soln[3*i*num_elems + j] = value.real();
        }
      for (unsigned int j=0; j<num_elems; ++j)
        {
          Number value = soln[i*num_vars + j];
          complex_soln[3*i*num_elems + num_elems +j] = value.imag();
        }
      for (unsigned int j=0; j<num_elems; ++j)
        {
          Number value = soln[i*num_vars + j];
          complex_soln[3*i*num_elems + 2*num_elems + j] = std::abs(value);
        }
    }

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
  write_equation_systems(fname,es,system_names);

  if (MeshOutput<MeshBase>::mesh().processor_id())
    return;

  exio_helper->write_timestep(timestep, time);
}



void ExodusII_IO::write (const std::string & fname)
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // We may need to gather a DistributedMesh to output it, making that
  // const qualifier in our constructor a dirty lie
  // The "true" specifies that we only need the mesh serialized to processor 0
  MeshSerializer serialize(MeshInput<MeshBase>::mesh(), !MeshOutput<MeshBase>::_is_parallel_format, true);

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

  int num_vars = cast_int<int>(names.size());
  int num_nodes = 0;
  for (const auto & elem : mesh.active_element_ptr_range())
    num_nodes += elem->n_nodes();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

  std::vector<std::string> complex_names = exio_helper->get_complex_names(names);

  // Call helper function for opening/initializing data, giving it the
  // complex variable names
  this->write_nodal_data_common(fname, complex_names, /*continuous=*/false);
#else
  // Call helper function for opening/initializing data
  this->write_nodal_data_common(fname, names, /*continuous=*/false);
#endif

  if (mesh.processor_id())
    return;

  for (int c=0; c<num_vars; c++)
    {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      std::vector<Real> real_parts(num_nodes);
      std::vector<Real> imag_parts(num_nodes);
      std::vector<Real> magnitudes(num_nodes);

      for (int i=0; i<num_nodes; ++i)
        {
          real_parts[i] = soln[i*num_vars + c].real();
          imag_parts[i] = soln[i*num_vars + c].imag();
          magnitudes[i] = std::abs(soln[i*num_vars + c]);
        }
      exio_helper->write_nodal_values(3*c+1,real_parts,_timestep);
      exio_helper->write_nodal_values(3*c+2,imag_parts,_timestep);
      exio_helper->write_nodal_values(3*c+3,magnitudes,_timestep);
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
          exio_helper->open(fname.c_str(), /*read_only=*/false);
          // If we're appending, it's not valid to call exio_helper->initialize()
          // or exio_helper->initialize_nodal_variables(), but we do need to set up
          // certain aspects of the Helper object itself, such as the number of nodes
          // and elements.  We do that by reading the header...
          exio_helper->read_header();

          // ...and reading the block info
          exio_helper->read_block_info();
        }
      else
        {
          exio_helper->create(fname);

          exio_helper->initialize(fname, mesh, !continuous);
          exio_helper->write_nodal_coordinates(mesh, !continuous);
          exio_helper->write_elements(mesh, !continuous);

          exio_helper->write_sidesets(mesh);
          exio_helper->write_nodesets(mesh);

          if ((mesh.get_boundary_info().n_edge_conds() > 0) && _verbose)
            libmesh_warning("Warning: Mesh contains edge boundary IDs, but these "
                            "are not supported by the ExodusII format.");

          exio_helper->initialize_nodal_variables(names);
        }
    }
  else
    {
      // We are already open for writing, so check that the filename
      // passed to this function matches the filename currently in use
      // by the helper.
      if (fname != exio_helper->current_filename)
        libmesh_error_msg("Error! This ExodusII_IO object is already associated with file: " \
                          << exio_helper->current_filename              \
                          << ", cannot use it with requested file: "    \
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

ExodusII_IO_Helper & ExodusII_IO::get_exio_helper()
{
  // Provide a warning when accessing the helper object
  // since it is a non-public API and is likely to see
  // future API changes
  libmesh_experimental();

  return *exio_helper;
}


// LIBMESH_HAVE_EXODUS_API is not defined, declare error() versions of functions...
#else



ExodusII_IO::~ExodusII_IO ()
{
}



void ExodusII_IO::read (const std::string &)
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

#endif // LIBMESH_HAVE_EXODUS_API
} // namespace libMesh
