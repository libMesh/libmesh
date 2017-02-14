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
#endif
  _timestep(1),
  _verbose(false),
  _append(false),
  _allow_empty_variables(false)
{
}


void ExodusII_IO::set_output_variables(const std::vector<std::string> & output_variables,
                                       bool allow_empty)
{
  _output_variables = output_variables;
  _allow_empty_variables = allow_empty;
}



void ExodusII_IO::copy_nodal_solution(System & system,
                                      std::string var_name,
                                      unsigned int timestep)
{
  libmesh_deprecated();
  copy_nodal_solution(system, var_name, var_name, timestep);
}



void ExodusII_IO::write_discontinuous_exodusII(const std::string & name,
                                               const EquationSystems & es,
                                               const std::set<std::string> * system_names)
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;

  es.build_variable_names  (solution_names, libmesh_nullptr, system_names);
  es.build_discontinuous_solution_vector (v, system_names);

  this->write_nodal_data_discontinuous(name, v, solution_names);
}




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

#ifdef DEBUG
  this->verbose(true);
#endif

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

    for (std::size_t e=0; e<exio_helper->elem_list.size(); e++)
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
        unsigned int side_index_offset = conv.get_shellface_index_offset();

        if(raw_side_index < side_index_offset)
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
    exio_helper->read_nodeset_info();

    for (int nodeset=0; nodeset<exio_helper->num_node_sets; nodeset++)
      {
        boundary_id_type nodeset_id =
          cast_int<boundary_id_type>(exio_helper->nodeset_ids[nodeset]);

        std::string nodeset_name = exio_helper->get_node_set_name(nodeset);
        if (!nodeset_name.empty())
          mesh.get_boundary_info().nodeset_name(nodeset_id) = nodeset_name;

        exio_helper->read_nodeset(nodeset);

        for (std::size_t node=0; node<exio_helper->node_list.size(); node++)
          {
            // As before, the entries in 'node_list' are 1-based
            // indcies into the node_num_map array, so we have to map
            // them.  See comment above.
            int libmesh_node_id = exio_helper->node_num_map[exio_helper->node_list[node] - 1] - 1;
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
  libmesh_deprecated();
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

  for (std::size_t i=0; i<exio_helper->nodal_var_values.size(); ++i)
    {
      const Node & node = MeshInput<MeshBase>::mesh().node_ref(i);

      if (node.n_comp(system.number(), var_num) > 0)
        {
          dof_id_type dof_index = node.dof_number(system.number(), var_num, 0);

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

      if (!elem)
        libmesh_error_msg("Error! Mesh returned NULL pointer for elem " << it->first);

      if (elem->n_comp(system.number(), var_num) > 0)
        {
          dof_id_type dof_index = elem->dof_number(system.number(), var_num, 0);

          // If the dof_index is local to this processor, set the value
          if ((dof_index >= system.solution->first_local_index()) && (dof_index < system.solution->last_local_index()))
            system.solution->set (dof_index, it->second);
        }
    }

  system.solution->close();
  system.update();
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
      std::vector<std::string>::iterator it = monomials.begin();
      for (; it!=monomials.end(); ++it)
        if (std::find(_output_variables.begin(), _output_variables.end(), *it) != _output_variables.end())
          names.push_back(*it);
    }

  // If we pass in a list of names to "get_solution" it'll filter the variables coming back
  std::vector<Number> soln;
  es.get_solution(soln, names);

  if(soln.empty()) // If there is nothing to write just return
    return;

  // The data must ultimately be written block by block.  This means that this data
  // must be sorted appropriately.
  if(MeshOutput<MeshBase>::mesh().processor_id())
    return;

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

  std::vector<std::string> complex_names = exio_helper->get_complex_names(names);

  exio_helper->initialize_element_variables(complex_names);

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

  exio_helper->write_element_values(mesh, complex_soln, _timestep);

#else
  exio_helper->initialize_element_variables(names);
  exio_helper->write_element_values(mesh, soln, _timestep);
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

  if(_allow_empty_variables || !_output_variables.empty())
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
      {
        MeshBase::const_node_iterator it = mesh.nodes_begin();
        const MeshBase::const_node_iterator end = mesh.nodes_end();
        for (; it != end; ++it)
          {
            const Node * node = *it;
            dof_id_type idx = node->id()*num_vars + c;
#ifdef LIBMESH_USE_REAL_NUMBERS
            cur_soln.push_back(soln[idx]);
#else
            real_parts.push_back(soln[idx].real());
            imag_parts.push_back(soln[idx].imag());
            magnitudes.push_back(std::abs(soln[idx]));
#endif
          }
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
  if(MeshOutput<MeshBase>::mesh().processor_id())
    return;

  if (!exio_helper->opened_for_writing)
    libmesh_error_msg("ERROR, ExodusII file must be initialized before outputting information records.");

  exio_helper->write_information_records(records);
}



void ExodusII_IO::write_global_data (const std::vector<Number> & soln,
                                     const std::vector<std::string> & names)
{
  if(MeshOutput<MeshBase>::mesh().processor_id())
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
                                  const Real time)
{
  _timestep = timestep;
  write_equation_systems(fname,es);

  if(MeshOutput<MeshBase>::mesh().processor_id())
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
    libMesh::out << "Warning: Appending in ExodusII_IO::write() does not make sense.\n"
                 << "Creating a new file instead!"
                 << std::endl;

  exio_helper->create(fname);
  exio_helper->initialize(fname,mesh);
  exio_helper->write_nodal_coordinates(mesh);
  exio_helper->write_elements(mesh);
  exio_helper->write_sidesets(mesh);
  exio_helper->write_nodesets(mesh);

  if( (mesh.get_boundary_info().n_edge_conds() > 0) &&
      _verbose )
    {
      libMesh::out << "Warning: Mesh contains edge boundary IDs, but these "
                   << "are not supported by the ExodusII format."
                   << std::endl;
    }
}



void ExodusII_IO::write_nodal_data_discontinuous (const std::string & fname,
                                                  const std::vector<Number> & soln,
                                                  const std::vector<std::string> & names)
{
  LOG_SCOPE("write_nodal_data_discontinuous()", "ExodusII_IO");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  int num_vars = cast_int<int>(names.size());
  int num_nodes = 0;
  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
  for ( ; it != end; ++it)
    num_nodes += (*it)->n_nodes();

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

          if( (mesh.get_boundary_info().n_edge_conds() > 0) &&
              _verbose )
            {
              libMesh::out << "Warning: Mesh contains edge boundary IDs, but these "
                           << "are not supported by the ExodusII format."
                           << std::endl;
            }

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



void ExodusII_IO::write_element_data (const EquationSystems &)
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
                                  const Real)
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

#endif // LIBMESH_HAVE_EXODUS_API
} // namespace libMesh
