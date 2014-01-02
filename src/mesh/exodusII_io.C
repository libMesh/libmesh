// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{

// ------------------------------------------------------------
// ExodusII_IO class members
ExodusII_IO::ExodusII_IO (MeshBase& mesh) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase> (mesh),
  ParallelObject(mesh),
#ifdef LIBMESH_HAVE_EXODUS_API
  exio_helper(new ExodusII_IO_Helper(*this)),
#endif
  _timestep(1),
  _verbose(false),
  _append(false),
  _allow_empty_variables(false)
{
}


void ExodusII_IO::set_output_variables(const std::vector<std::string> & output_variables, bool allow_empty)
{
  _output_variables = output_variables;
  _allow_empty_variables = allow_empty;
}



void ExodusII_IO::copy_nodal_solution(System& system, std::string var_name, unsigned int timestep)
{
  libmesh_deprecated();
  copy_nodal_solution(system, var_name, var_name, timestep);
}



void ExodusII_IO::write_discontinuous_exodusII(const std::string& name, const EquationSystems& es)
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;

  es.build_variable_names  (solution_names);
  es.build_discontinuous_solution_vector (v);

  this->write_nodal_data_discontinuous(name, v, solution_names);
}




// ------------------------------------------------------------
// When the Exodus API is present...
#ifdef LIBMESH_HAVE_EXODUS_API

ExodusII_IO::~ExodusII_IO ()
{
  exio_helper->close();
  delete exio_helper;
}



void ExodusII_IO::read (const std::string& fname)
{
  // Get a reference to the mesh we are reading
  MeshBase& mesh = MeshInput<MeshBase>::mesh();

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
      Node* added_node = mesh.add_point (Point(exio_helper->x[i], exio_helper->y[i], exio_helper->z[i]), exodus_id-1);

      // If the Mesh assigned an ID different from what is in the
      // Exodus file, we should probably error.
      if (added_node->id() != static_cast<unsigned>(exodus_id-1))
        {
          libMesh::err << "Error!  Mesh assigned node ID "
                       << added_node->id()
                       << " which is different from the (zero-based) Exodus ID "
                       << exodus_id-1
                       << "!"
                       << std::endl;
          libmesh_error();
        }
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
	  Elem* elem = Elem::build (conv.get_canonical_type()).release();
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
            {
              libMesh::err << "Error!  Mesh assigned ID "
                           << elem->id()
                           << " which is different from the (zero-based) Exodus ID "
                           << exodus_id-1
                           << "!"
                           << std::endl;
              libmesh_error();
            }

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
  for (unsigned int i=0; i!=4; ++i)
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
          mesh.boundary_info->sideset_name(exio_helper->get_side_set_id(i)) = sideset_name;
      }

    for (unsigned int e=0; e<exio_helper->elem_list.size(); e++)
      {
        // The numbers in the Exodus file sidesets should be thought
        // of as (1-based) indices into the elem_num_map array.  So,
        // to get the right element ID we have to:
        // 1.) Subtract 1 from elem_list[e] (to get a zero-based index)
        // 2.) Pass it through elem_num_map (to get the corresponding Exodus ID)
        // 3.) Subtract 1 from that, since libmesh is "zero"-based,
        //     even when the Exodus numbering doesn't start with 1.
        int libmesh_elem_id = exio_helper->elem_num_map[exio_helper->elem_list[e] - 1] - 1;

	// Set any relevant node/edge maps for this element
        Elem * elem = mesh.elem(libmesh_elem_id);

	const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(elem->type());

        // Add this (elem,side,id) triplet to the BoundaryInfo object.
	mesh.boundary_info->add_side (libmesh_elem_id,
				      conv.get_side_map(exio_helper->side_list[e]-1),
				      exio_helper->id_list[e]);
      }
  }

  // Read nodeset info
  {
    exio_helper->read_nodeset_info();

    for (int nodeset=0; nodeset<exio_helper->num_node_sets; nodeset++)
      {
        int nodeset_id = exio_helper->nodeset_ids[nodeset];

        std::string nodeset_name = exio_helper->get_node_set_name(nodeset);
        if (!nodeset_name.empty())
          mesh.boundary_info->nodeset_name(nodeset_id) = nodeset_name;

        exio_helper->read_nodeset(nodeset);

        for (unsigned int node=0; node<exio_helper->node_list.size(); node++)
          {
            // As before, the entries in 'node_list' are 1-based
            // indcies into the node_num_map array, so we have to map
            // them.  See comment above.
            int libmesh_node_id = exio_helper->node_num_map[exio_helper->node_list[node] - 1] - 1;
            mesh.boundary_info->add_node(libmesh_node_id, nodeset_id);
          }
      }
  }

#if LIBMESH_DIM < 3
  if (mesh.mesh_dimension() > LIBMESH_DIM)
    {
      libMesh::err << "Cannot open dimension " <<
		      mesh.mesh_dimension() <<
		      " mesh file when configured without " <<
                      mesh.mesh_dimension() << "D support." <<
                      std::endl;
      libmesh_error();
    }
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



void ExodusII_IO::set_coordinate_offset(Point p)
{
  libmesh_deprecated();
  exio_helper->set_coordinate_offset(p);
}



void ExodusII_IO::append(bool val)
{
  _append = val;
}



const std::vector<Real>& ExodusII_IO::get_time_steps()
{
  if (!exio_helper->opened_for_reading)
    {
      libMesh::err << "ERROR, ExodusII file must be opened for reading before calling ExodusII_IO::get_time_steps()!" << std::endl;
      libmesh_error();
    }

  exio_helper->read_time_steps();
  return exio_helper->time_steps;
}



int ExodusII_IO::get_num_time_steps()
{
  if (!exio_helper->opened_for_reading && !exio_helper->opened_for_writing)
    {
      libMesh::err << "ERROR, ExodusII file must be opened for reading or writing before calling ExodusII_IO::get_num_time_steps()!" << std::endl;
      libmesh_error();
    }

  exio_helper->read_num_time_steps();
  return exio_helper->num_time_steps;
}



void ExodusII_IO::copy_nodal_solution(System& system, std::string system_var_name, std::string exodus_var_name, unsigned int timestep)
{
  if (!exio_helper->opened_for_reading)
    {
      libMesh::err << "ERROR, ExodusII file must be opened for reading before copying a nodal solution!" << std::endl;
      libmesh_error();
    }

  exio_helper->read_nodal_var_values(exodus_var_name, timestep);

  const unsigned int var_num = system.variable_number(system_var_name);

  for (unsigned int i=0; i<exio_helper->nodal_var_values.size(); ++i)
    {
      const Node* node = MeshInput<MeshBase>::mesh().query_node_ptr(i);

      if (!node)
        {
          libMesh::err << "Error! Mesh returned NULL pointer for node " << i << std::endl;
          libmesh_error();
        }

      dof_id_type dof_index = node->dof_number(system.number(), var_num, 0);

      // If the dof_index is local to this processor, set the value
      if ((dof_index >= system.solution->first_local_index()) && (dof_index < system.solution->last_local_index()))
	system.solution->set (dof_index, exio_helper->nodal_var_values[i]);
    }

  system.solution->close();
  system.update();
}



void ExodusII_IO::copy_elemental_solution(System& system, std::string system_var_name, std::string exodus_var_name, unsigned int timestep)
{
  if (!exio_helper->opened_for_reading)
    {
      libMesh::err << "ERROR, ExodusII file must be opened for reading before copying an elemental solution!" << std::endl;
      libmesh_error();
    }

  exio_helper->read_elemental_var_values(exodus_var_name, timestep);

  const unsigned int var_num = system.variable_number(system_var_name);
  if (system.variable_type(var_num) != FEType(CONSTANT, MONOMIAL))
  {
    libMesh::err << "Error! Trying to copy elemental solution into a variable that is not of CONSTANT MONOMIAL type. " << std::endl;
    libmesh_error();
  }

  for (unsigned int i=0; i<exio_helper->elem_var_values.size(); ++i)
    {
      const Elem * elem = MeshInput<MeshBase>::mesh().query_elem(i);

      if (!elem)
        {
          libMesh::err << "Error! Mesh returned NULL pointer for elem " << i << std::endl;
          libmesh_error();
        }

      dof_id_type dof_index = elem->dof_number(system.number(), var_num, 0);

      // If the dof_index is local to this processor, set the value
      if ((dof_index >= system.solution->first_local_index()) && (dof_index < system.solution->last_local_index()))
        system.solution->set (dof_index, exio_helper->elem_var_values[i]);
    }

  system.solution->close();
  system.update();
}



void ExodusII_IO::write_element_data (const EquationSystems & es)
{
  // Be sure the file has been opened for writing!
  if (!exio_helper->opened_for_writing)
    {
      libMesh::err << "ERROR, ExodusII file must be initialized "
                   << "before outputting element variables.\n"
                   << std::endl;
      libmesh_error();
    }

  // This function currently only works on SerialMeshes. We rely on
  // having a reference to a non-const MeshBase object from our
  // MeshInput parent class to construct a MeshSerializer object,
  // similar to what is done in ExodusII_IO::write().  Note that
  // calling ExodusII_IO::write_timestep() followed by
  // ExodusII_IO::write_element_data() when the underlying Mesh is a
  // ParallelMesh will result in an unnecessary additional
  // serialization/re-parallelization step.
  MeshSerializer serialize(MeshInput<MeshBase>::mesh(), !MeshOutput<MeshBase>::_is_parallel_format);

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

  // The data must ultimately be written block by block.  This means that this data
  // must be sorted appropriately.

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  exio_helper->initialize_element_variables(names);
  exio_helper->write_element_values(mesh, soln, _timestep);
}



void ExodusII_IO::write_nodal_data (const std::string& fname,
				    const std::vector<Number>& soln,
				    const std::vector<std::string>& names)
{
  START_LOG("write_nodal_data()", "ExodusII_IO");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  int num_vars = libmesh_cast_int<int>(names.size());
  dof_id_type num_nodes = mesh.n_nodes();

  // The names of the variables to be output
  std::vector<std::string> output_names;

  if(_allow_empty_variables || !_output_variables.empty())
    output_names = _output_variables;
  else
    output_names = names;

  // Call helper function for opening/initializing data
  this->write_nodal_data_common(fname, output_names, /*continuous=*/true);

  // This will count the number of variables actually output
  for (int c=0; c<num_vars; c++)
    {
      std::vector<std::string>::iterator pos =
        std::find(output_names.begin(), output_names.end(), names[c]);
      if (pos == output_names.end())
        continue;

      unsigned int variable_name_position =
	libmesh_cast_int<unsigned int>(pos - output_names.begin());

      std::vector<Number> cur_soln(num_nodes);

      // Copy out this variable's solution
      for(dof_id_type i=0; i<num_nodes; i++)
        cur_soln[i] = soln[i*num_vars + c];

      exio_helper->write_nodal_values(variable_name_position+1,cur_soln,_timestep);
    }

  STOP_LOG("write_nodal_data()", "ExodusII_IO");
}




void ExodusII_IO::write_information_records (const std::vector<std::string>& records)
{
  if (!exio_helper->opened_for_writing)
    {
      libMesh::err << "ERROR, ExodusII file must be initialized "
                   << "before outputting information records.\n"
                   << std::endl;
      libmesh_error();
    }

  exio_helper->write_information_records(records);
}



void ExodusII_IO::write_global_data (const std::vector<Number>& soln,
                                     const std::vector<std::string>& names)
{
  if (!exio_helper->opened_for_writing)
    {
      libMesh::err << "ERROR, ExodusII file must be initialized "
                   << "before outputting global variables.\n"
                   << std::endl;
      libmesh_error();
    }

  exio_helper->initialize_global_variables(names);
  exio_helper->write_global_values(soln, _timestep);
}



void ExodusII_IO::write_timestep (const std::string& fname,
				  const EquationSystems& es,
				  const int timestep,
				  const Real time)
{
  _timestep = timestep;
  write_equation_systems(fname,es);
  exio_helper->write_timestep(timestep, time);
}



void ExodusII_IO::write (const std::string& fname)
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // We may need to gather a ParallelMesh to output it, making that
  // const qualifier in our constructor a dirty lie
  MeshSerializer serialize(const_cast<MeshBase&>(mesh), !MeshOutput<MeshBase>::_is_parallel_format);

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

  if( (mesh.boundary_info->n_edge_conds() > 0) &&
       _verbose )
  {
    libMesh::out << "Warning: Mesh contains edge boundary IDs, but these "
                 << "are not supported by the ExodusII format."
                 << std::endl;
  }
}



void ExodusII_IO::write_nodal_data_discontinuous (const std::string& fname,
                                                  const std::vector<Number>& soln,
                                                  const std::vector<std::string>& names)
{
  START_LOG("write_nodal_data_discontinuous()", "ExodusII_IO");

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  int num_vars = libmesh_cast_int<int>(names.size());
  int num_nodes = 0;
  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
  for ( ; it != end; ++it)
    num_nodes += (*it)->n_nodes();

  // Call helper function for opening/initializing data
  this->write_nodal_data_common(fname, names, /*continuous=*/false);

  if (this->processor_id() == 0)
    for (int c=0; c<num_vars; c++)
      {
        // Copy out this variable's solution
        std::vector<Number> cur_soln(num_nodes);

        for(int i=0; i<num_nodes; i++)
          cur_soln[i] = soln[i*num_vars + c];

        exio_helper->write_nodal_values(c+1,cur_soln,_timestep);
      }

  STOP_LOG("write_nodal_data_discontinuous()", "ExodusII_IO");
}



void ExodusII_IO::write_nodal_data_common(std::string fname,
                                          const std::vector<std::string>& names,
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

          if (continuous)
            {
              exio_helper->initialize(fname,mesh);
              exio_helper->write_nodal_coordinates(mesh);
              exio_helper->write_elements(mesh);
            }
          else
            {
              exio_helper->initialize_discontinuous(fname,mesh);
              exio_helper->write_nodal_coordinates_discontinuous(mesh);
              exio_helper->write_elements_discontinuous(mesh);
            }
          exio_helper->write_sidesets(mesh);
          exio_helper->write_nodesets(mesh);

          if( (mesh.boundary_info->n_edge_conds() > 0) &&
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
        {
          libMesh::err << "Error! This ExodusII_IO object is already associated with file: "
                       << exio_helper->current_filename
                       << ", cannot use it with requested file: "
                       << fname
                       << std::endl;
          libmesh_error();
        }
    }
}



// LIBMESH_HAVE_EXODUS_API is not defined, declare error() versions of functions...
#else



ExodusII_IO::~ExodusII_IO ()
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::read (const std::string&)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::verbose (bool)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::use_mesh_dimension_instead_of_spatial_dimension(bool)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::set_coordinate_offset(Point)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::append(bool val)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



const std::vector<Real>& ExodusII_IO::get_time_steps()
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



int ExodusII_IO::get_num_time_steps()
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}


void ExodusII_IO::copy_nodal_solution(System&, std::string, std::string, unsigned int)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::copy_elemental_solution(System&, std::string, std::string, unsigned int)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::write_element_data (const EquationSystems& es)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::write_nodal_data (const std::string&, const std::vector<Number>&, const std::vector<std::string>&)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::write_information_records (const std::vector<std::string>&)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::write_global_data (const std::vector<Number>&, const std::vector<std::string>&)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::write_timestep (const std::string&, const EquationSystems&, const int, const Real)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::write (const std::string&)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



void ExodusII_IO::write_nodal_data_discontinuous (const std::string&, const std::vector<Number>&, const std::vector<std::string>&)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}



  void ExodusII_IO::write_nodal_data_common(std::string,
					    const std::vector<std::string>&,
					    bool)
{
  libMesh::err << "ERROR, ExodusII API is not defined." << std::endl;
  libmesh_error();
}

#endif // LIBMESH_HAVE_EXODUS_API
} // namespace libMesh
