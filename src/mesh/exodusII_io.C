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
  _verbose (false)
{
}



ExodusII_IO::~ExodusII_IO ()
{
#ifndef LIBMESH_HAVE_EXODUS_API

  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << std::endl;
  libmesh_error();

#else

  exio_helper->close();

  delete exio_helper;

#endif
}



void ExodusII_IO::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;

#ifdef LIBMESH_HAVE_EXODUS_API
  // Set the verbose flag in the helper object
  // as well.
  exio_helper->verbose(_verbose);
#endif
}



  void ExodusII_IO::use_mesh_dimension_instead_of_spatial_dimension(bool val)
  {
#ifdef LIBMESH_HAVE_EXODUS_API
    exio_helper->use_mesh_dimension_instead_of_spatial_dimension(val);
#endif
  }



void ExodusII_IO::read (const std::string& fname)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
//  libmesh_assert_equal_to (this->processor_id(), 0);

#ifndef LIBMESH_HAVE_EXODUS_API

  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << "Input file " << fname << " cannot be read"
	        << std::endl;
  libmesh_error();

#else

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


  ExodusII_IO_Helper::ElementMaps em;     // Instantiate the ElementMaps interface

  exio_helper->open(fname.c_str());       // Open the exodus file, if possible
  exio_helper->read_header();             // Get header information from exodus file
  exio_helper->print_header();            // Print header information

  //assertion fails due to inconsistent mesh dimension
//  libmesh_assert_equal_to (static_cast<unsigned int>(exio_helper->get_num_dim()), mesh.mesh_dimension()); // Be sure number of dimensions
                                                                                // is equal to the number of
                                                                                // dimensions in the mesh supplied.

  exio_helper->read_nodes();                        // Read nodes from the exodus file
  mesh.reserve_nodes(exio_helper->get_num_nodes()); // Reserve space for the nodes.

  // Loop over the nodes, create Nodes with local processor_id 0.
  for (int i=0; i<exio_helper->get_num_nodes(); i++)
    mesh.add_point (Point(exio_helper->get_x(i),
			  exio_helper->get_y(i),
			  exio_helper->get_z(i)), i);

  libmesh_assert_equal_to (static_cast<unsigned int>(exio_helper->get_num_nodes()), mesh.n_nodes());

  exio_helper->read_block_info();                 // Get information about all the blocks
  mesh.reserve_elem(exio_helper->get_num_elem()); // Reserve space for the elements


  // Read in the element connectivity for each block.
  int nelem_last_block = 0;

  std::map<int, dof_id_type> exodus_id_to_mesh_id;

  // Loop over all the blocks
  for (int i=0; i<exio_helper->get_num_elem_blk(); i++)
    {
      // Read the information for block i
      exio_helper->read_elem_in_block (i);
      int subdomain_id = exio_helper->get_block_id(i);

      // populate the map of names
      mesh.subdomain_name(static_cast<subdomain_id_type>(subdomain_id)) =
        exio_helper->get_block_name(i);

      // Set any relevant node/edge maps for this element
      const std::string type_str (exio_helper->get_elem_type());
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(type_str);
      //if (_verbose)
      //libMesh::out << "Reading a block of " << type_str << " elements." << std::endl;

      // Loop over all the faces in this block
      int jmax = nelem_last_block+exio_helper->get_num_elem_this_blk();
      for (int j=nelem_last_block; j<jmax; j++)
	{
	  Elem* elem = Elem::build (conv.get_canonical_type()).release();
	  libmesh_assert (elem);
          elem->subdomain_id() = static_cast<subdomain_id_type>(subdomain_id) ;
          //elem->set_id(j);// Don't try to second guess the Element ID setting scheme!

          elems_of_dimension[elem->dim()] = true;

	  elem = mesh.add_elem (elem); // Catch the Elem pointer that the Mesh throws back

          exodus_id_to_mesh_id[j+1] = elem->id();

	  // Set all the nodes for this element
	  for (int k=0; k<exio_helper->get_num_nodes_per_elem(); k++)
	    {
	      int gi = (j-nelem_last_block)*exio_helper->get_num_nodes_per_elem() + conv.get_node_map(k); // global index
	      int node_number   = exio_helper->get_connect(gi);             // Global node number (1-based)
	      elem->set_node(k) = mesh.node_ptr((node_number-1)); // Set node number
	                                                          // Subtract 1 since
		                                                  // exodus is internally 1-based
	    }
	}

      // running sum of # of elements per block,
      // (should equal total number of elements in the end)
      nelem_last_block += exio_helper->get_num_elem_this_blk();
    }
  libmesh_assert_equal_to (static_cast<unsigned int>(nelem_last_block), mesh.n_elem());

   // Set the mesh dimension to the largest encountered for an element
  for (unsigned int i=0; i!=4; ++i)
    if (elems_of_dimension[i])
      mesh.set_mesh_dimension(i);

  // Read in sideset information -- this is useful for applying boundary conditions
  {
    exio_helper->read_sideset_info(); // Get basic information about ALL sidesets
    int offset=0;
    for (int i=0; i<exio_helper->get_num_side_sets(); i++)
      {
	offset += (i > 0 ? exio_helper->get_num_sides_per_set(i-1) : 0); // Compute new offset
	exio_helper->read_sideset (i, offset);

        mesh.boundary_info->sideset_name(exio_helper->get_side_set_id(i)) =
          exio_helper->get_side_set_name(i);
      }

    const std::vector<int>& elem_list = exio_helper->get_elem_list();
    const std::vector<int>& side_list = exio_helper->get_side_list();
    const std::vector<int>& id_list   = exio_helper->get_id_list();

    for (unsigned int e=0; e<elem_list.size(); e++)
      {
	// Set any relevant node/edge maps for this element

        Elem * elem = mesh.elem(exodus_id_to_mesh_id[elem_list[e]]);

	const ExodusII_IO_Helper::Conversion conv =
	  em.assign_conversion(elem->type());

	mesh.boundary_info->add_side (exodus_id_to_mesh_id[elem_list[e]],
				      conv.get_side_map(side_list[e]-1),
				      id_list[e]);
      }
  }

  // Read nodeset info
  {
    exio_helper->read_nodeset_info();

    for (int nodeset=0; nodeset<exio_helper->get_num_node_sets(); nodeset++)
      {
        int nodeset_id = exio_helper->get_nodeset_id(nodeset);

        mesh.boundary_info->nodeset_name(nodeset_id) =
          exio_helper->get_node_set_name(nodeset);

        exio_helper->read_nodeset(nodeset);

        const std::vector<int>& node_list = exio_helper->get_node_list();

        for(unsigned int node=0; node<node_list.size(); node++)
          mesh.boundary_info->add_node(node_list[node]-1, nodeset_id);
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

#endif
}



#ifndef LIBMESH_HAVE_EXODUS_API

const std::vector<Real>& ExodusII_IO::get_time_steps()
{
  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << std::endl;
  libmesh_error();
}

#else

const std::vector<Real>& ExodusII_IO::get_time_steps()
{
  return exio_helper->get_time_steps();
}

#endif




void ExodusII_IO::copy_nodal_solution(System& system, std::string var_name, unsigned int timestep)
{
  libmesh_deprecated();
  copy_nodal_solution(system, var_name, var_name, timestep);
}



#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::copy_nodal_solution(System&, std::string, std::string, unsigned int)
{

  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << std::endl;
  libmesh_error();
}

#else

void ExodusII_IO::copy_nodal_solution(System& system, std::string system_var_name, std::string exodus_var_name, unsigned int timestep)
{
  // FIXME: Do we need to call get_time_steps() at all?
  /*const std::vector<double>& time_steps = */
  exio_helper->get_time_steps();

  const std::vector<Real> & nodal_values = exio_helper->get_nodal_var_values(exodus_var_name, timestep);

  const unsigned int var_num = system.variable_number(system_var_name);

  for (unsigned int i=0; i<nodal_values.size(); ++i)
    {
      const Node* node = MeshInput<MeshBase>::mesh().node_ptr(i);

      if (!node)
        {
          libMesh::err << "Error! Mesh returned NULL pointer for node " << i << std::endl;
          libmesh_error();
        }

      dof_id_type dof_index = node->dof_number(system.number(), var_num, 0);

      // If the dof_index is local to this processor, set the value
      if ((dof_index >= system.solution->first_local_index()) && (dof_index < system.solution->last_local_index()))
	system.solution->set (dof_index, nodal_values[i]);
    }

  system.solution->close();
  system.update();
}

#endif




#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::write_element_data (const EquationSystems& es)
{

  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << std::endl;
  libmesh_error();
}

#else

void ExodusII_IO::write_element_data (const EquationSystems & es)
{
  // The first step is to collect the element data onto this processor.
  // We want the constant monomial data.

  std::vector<Number> soln;
  std::vector<std::string> names;

  // If _output_variables is populated we need to filter the monomials we output
  if (_output_variables.size())
  {
    std::vector<std::string> monomials;
    const FEType type(CONSTANT, MONOMIAL);
    es.build_variable_names(monomials, &type);

    for (std::vector<std::string>::iterator it = monomials.begin(); it != monomials.end(); ++it)
      if (std::find(_output_variables.begin(), _output_variables.end(), *it) != _output_variables.end())
        names.push_back(*it);
  }

  // If we pass in a list of names to "get_solution" it'll filter the variables comming back
  es.get_solution( soln, names );

  // The data must ultimately be written block by block.  This means that this data
  // must be sorted appropriately.

  if (!exio_helper->created())
    {
      libMesh::err << "ERROR, ExodusII file must be initialized "
                   << "before outputting element variables.\n"
                   << std::endl;
      libmesh_error();
    }

  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  exio_helper->initialize_element_variables( mesh, names );
  exio_helper->write_element_values(mesh,soln,_timestep);
}

#endif



#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::write_nodal_data (const std::string& ,
				    const std::vector<Number>& ,
				    const std::vector<std::string>& )
{

  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << std::endl;
  libmesh_error();
}

#else
void ExodusII_IO::write_discontinuous_exodusII (const std::string& name,
				     const EquationSystems& es)
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;

  es.build_variable_names  (solution_names);
  es.build_discontinuous_solution_vector (v);

  this->write_nodal_data_discontinuous(name, v, solution_names);
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

  // FIXME: Will we ever _not_ need to do this?
  // DRG: Yes... when writing multiple timesteps to the same file.
  if (!exio_helper->created())
    {
      exio_helper->create(fname);
      exio_helper->initialize_discontinuous(fname,mesh);
      exio_helper->write_nodal_coordinates_discontinuous(mesh);
      exio_helper->write_elements_discontinuous(mesh);
      exio_helper->write_sidesets(mesh);
      exio_helper->write_nodesets(mesh);
      exio_helper->initialize_nodal_variables(names);
    }

  if (this->processor_id() == 0)
    for (int c=0; c<num_vars; c++)
      {
        //Copy out this variable's solution
        std::vector<Number> cur_soln(num_nodes);

        for(int i=0; i<num_nodes; i++)
          cur_soln[i] = soln[i*num_vars + c];//c*num_nodes+i];

        exio_helper->write_nodal_values(c+1,cur_soln,_timestep);
      }

  STOP_LOG("write_nodal_data_discontinuous()", "ExodusII_IO");
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

  if(_output_variables.size())
    output_names = _output_variables;
  else
    output_names = names;

  // FIXME: Will we ever _not_ need to do this?
  // DRG: Yes... when writing multiple timesteps to the same file.
  if (!exio_helper->created())
    {
      exio_helper->create(fname);
      exio_helper->initialize(fname,mesh);
      exio_helper->write_nodal_coordinates(mesh);
      exio_helper->write_elements(mesh);
      exio_helper->write_sidesets(mesh);
      exio_helper->write_nodesets(mesh);
      exio_helper->initialize_nodal_variables(output_names);
    }

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

      //Copy out this variable's solution
      for(dof_id_type i=0; i<num_nodes; i++)
        cur_soln[i] = soln[i*num_vars + c];//c*num_nodes+i];

      exio_helper->write_nodal_values(variable_name_position+1,cur_soln,_timestep);
    }

  STOP_LOG("write_nodal_data()", "ExodusII_IO");
}

#endif

#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::write_information_records ( const std::vector<std::string>& )
{

  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << std::endl;
  libmesh_error();
}

#else

void ExodusII_IO::write_information_records (const std::vector<std::string>& records)
{
  if (!exio_helper->created())
    {
      libMesh::err << "ERROR, ExodusII file must be initialized "
                   << "before outputting information records.\n"
                   << std::endl;
      libmesh_error();
    }

  exio_helper->write_information_records( records );
}

#endif

#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::write_global_data (const std::vector<Number>& ,
				    const std::vector<std::string>& )
{

  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << std::endl;
  libmesh_error();
}

#else

void ExodusII_IO::write_global_data (const std::vector<Number>& soln,
                                     const std::vector<std::string>& names)
{
  if (!exio_helper->created())
    {
      libMesh::err << "ERROR, ExodusII file must be initialized "
                   << "before outputting global variables.\n"
                   << std::endl;
      libmesh_error();
    }

  exio_helper->initialize_global_variables( names );
  exio_helper->write_global_values( soln, _timestep );
}

#endif


#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::write_timestep (const std::string& ,
				  const EquationSystems& ,
				  const int ,
				  const Real )
{

  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << std::endl;
  libmesh_error();
}

#else

void ExodusII_IO::write_timestep (const std::string& fname,
				  const EquationSystems& es,
				  const int timestep,
				  const Real time)
{
  _timestep=timestep;
  write_equation_systems(fname,es);

  exio_helper->write_timestep(timestep, time);
}

#endif



#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::write (const std::string& )
{
  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
	        << std::endl;
  libmesh_error();
}


#else

void ExodusII_IO::write (const std::string& fname)
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // We may need to gather a ParallelMesh to output it, making that
  // const qualifier in our constructor a dirty lie
  MeshSerializer serialize(const_cast<MeshBase&>(mesh), !MeshOutput<MeshBase>::_is_parallel_format);

  libmesh_assert( !exio_helper->created() );

  exio_helper->create(fname);
  exio_helper->initialize(fname,mesh);
  exio_helper->write_nodal_coordinates(mesh);
  exio_helper->write_elements(mesh);
  exio_helper->write_sidesets(mesh);
  exio_helper->write_nodesets(mesh);
  // Note: the file is closed automatically by the ExodusII_IO destructor.
}

#endif

} // namespace libMesh
