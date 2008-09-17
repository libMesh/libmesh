// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "exodusII_io.h"
#include "boundary_info.h"
#include "mesh_base.h"
#include "enum_elem_type.h"
#include "elem.h"
#include "system.h"
#include "numeric_vector.h"
#include "exodusII_io_helper.h"





// ------------------------------------------------------------
// ExodusII_IO class members
ExodusII_IO::ExodusII_IO (MeshBase& mesh) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase> (mesh),
  _timestep(1),
  _verbose (false)
{
}



ExodusII_IO::~ExodusII_IO ()
{
#ifndef LIBMESH_HAVE_EXODUS_API

  std::cerr <<  "ERROR, ExodusII API is not defined.\n"
	    << std::endl;
  libmesh_error();
    
#else

  exio_helper.close();
  
#endif
}



void ExodusII_IO::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;

#ifdef LIBMESH_HAVE_EXODUS_API
  // Set the verbose flag in the helper object
  // as well.
  exio_helper.verbose(_verbose);
#endif
}



void ExodusII_IO::read (const std::string& fname)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert(libMesh::processor_id() == 0);

#ifndef LIBMESH_HAVE_EXODUS_API

  std::cerr <<  "ERROR, ExodusII API is not defined.\n"
	    << "Input file " << fname << " cannot be read"
	    << std::endl;
  libmesh_error();
    
#else
  
  // Get a reference to the mesh we are reading
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  
  // Clear any existing mesh data
  mesh.clear();
  
  if (mesh.mesh_dimension() == 1) // No support for 1D ExodusII meshes
    libmesh_not_implemented();
  
#ifdef DEBUG
  this->verbose(true);
#endif

  
  ExodusII_IO_Helper::ElementMaps em;     // Instantiate the ElementMaps interface
    
  exio_helper.open(fname.c_str());       // Open the exodus file, if possible
  exio_helper.read_header();             // Get header information from exodus file
  exio_helper.print_header();            // Print header information

  libmesh_assert(static_cast<unsigned int>(exio_helper.get_num_dim()) == mesh.mesh_dimension()); // Be sure number of dimensions
                                                                                // is equal to the number of 
                                                                                // dimensions in the mesh supplied.
  
  exio_helper.read_nodes();                        // Read nodes from the exodus file
  mesh.reserve_nodes(exio_helper.get_num_nodes()); // Reserve space for the nodes.
  
  // Loop over the nodes, create Nodes with local processor_id 0.
  for (int i=0; i<exio_helper.get_num_nodes(); i++)
    mesh.add_point (Point(exio_helper.get_x(i),
			  exio_helper.get_y(i),
			  exio_helper.get_z(i)), i);
  
  libmesh_assert (static_cast<unsigned int>(exio_helper.get_num_nodes()) == mesh.n_nodes());

  exio_helper.read_block_info();                 // Get information about all the blocks
  mesh.reserve_elem(exio_helper.get_num_elem()); // Reserve space for the elements
   

  // Read in the element connectivity for each block.
  int nelem_last_block = 0;

  // Loop over all the blocks
  for (int i=0; i<exio_helper.get_num_elem_blk(); i++)
    {
      // Read the information for block i
      exio_helper.read_elem_in_block (i);
      int subdomain_id = exio_helper.get_block_id(i);

      // Set any relevant node/edge maps for this element
      const std::string type_str (exio_helper.get_elem_type());
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(type_str); 
      //if (_verbose)
      //std::cout << "Reading a block of " << type_str << " elements." << std::endl;
      
      // Loop over all the faces in this block
      int jmax = nelem_last_block+exio_helper.get_num_elem_this_blk();
      for (int j=nelem_last_block; j<jmax; j++)
	{
	  Elem* elem = Elem::build (conv.get_canonical_type()).release();
	  libmesh_assert (elem);
          elem->subdomain_id() = subdomain_id;
          //elem->set_id(j);// Don't try to second guess the Element ID setting scheme!
	  elem = mesh.add_elem (elem); // Catch the Elem pointer that the Mesh throws back
	    
	  // Set all the nodes for this element
	  for (int k=0; k<exio_helper.get_num_nodes_per_elem(); k++)
	    {
	      int gi = (j-nelem_last_block)*exio_helper.get_num_nodes_per_elem() + conv.get_node_map(k); // global index 
	      int node_number   = exio_helper.get_connect(gi);             // Global node number (1-based)
	      elem->set_node(k) = mesh.node_ptr((node_number-1)); // Set node number
	                                                          // Subtract 1 since
		                                                  // exodus is internally 1-based
	    }
	}
      
      // running sum of # of elements per block,
      // (should equal total number of elements in the end)
      nelem_last_block += exio_helper.get_num_elem_this_blk();
    }
  libmesh_assert (static_cast<unsigned int>(nelem_last_block) == mesh.n_elem());
  
  // Read in sideset information -- this is useful for applying boundary conditions
  {
    exio_helper.read_sideset_info(); // Get basic information about ALL sidesets
    int offset=0;
    for (int i=0; i<exio_helper.get_num_side_sets(); i++)
      {
	offset += (i > 0 ? exio_helper.get_num_sides_per_set(i-1) : 0); // Compute new offset
	exio_helper.read_sideset (i, offset);
      }
    
    const std::vector<int>& elem_list = exio_helper.get_elem_list();
    const std::vector<int>& side_list = exio_helper.get_side_list();
    const std::vector<int>& id_list   = exio_helper.get_id_list();

    for (unsigned int e=0; e<elem_list.size(); e++)
      {
	// Set any relevant node/edge maps for this element
	const ExodusII_IO_Helper::Conversion conv =
	  em.assign_conversion(mesh.elem(elem_list[e]-1)->type());
	
	mesh.boundary_info->add_side (elem_list[e]-1,
				      conv.get_side_map(side_list[e]-1),
				      id_list[e]);
      }
  }
#endif
}



#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::copy_nodal_solution(System& , std::string)
{

  std::cerr <<  "ERROR, ExodusII API is not defined.\n"
	    << std::endl;
  libmesh_error();
}

#else

void ExodusII_IO::copy_nodal_solution(System& system, std::string nodal_var_name)
{
  // FIXME: Do we need to call get_time_steps() at all?
  /*const std::vector<double>& time_steps = */
  exio_helper.get_time_steps();

  //For now just read the first timestep (1)
  const std::vector<double> & nodal_values = exio_helper.get_nodal_var_values(nodal_var_name,1);


  //const DofMap & dof_map = system.get_dof_map();

  const unsigned int var_num = system.variable_number(nodal_var_name);

  
  for (unsigned int i=0; i<nodal_values.size(); ++i)
    {
      const unsigned int dof_index = MeshInput<MeshBase>::mesh().node_ptr(i)->dof_number(system.number(),var_num,0);

      // If the dof_index is local to this processor, set the value
      if ((dof_index >= system.solution->first_local_index()) && (dof_index <  system.solution->last_local_index()))
	system.solution->set (dof_index, nodal_values[i]);
    }

  system.update();
}
  
#endif



#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::write_nodal_data (const std::string& ,
				    const std::vector<Number>& ,
				    const std::vector<std::string>& )
{

  std::cerr <<  "ERROR, ExodusII API is not defined.\n"
	    << std::endl;
  libmesh_error();
}

#else

void ExodusII_IO::write_nodal_data (const std::string& fname,
				    const std::vector<Number>& soln,
				    const std::vector<std::string>& names)
{
  if (libMesh::processor_id() == 0)
    {
      const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

      int num_vars = names.size();
      int num_nodes = mesh.n_nodes();
  
      // FIXME: Will we ever _not_ need to do this?
      if (!exio_helper.created())
	{
	  exio_helper.create(fname);
	  exio_helper.initialize(fname,mesh);
	  exio_helper.write_nodal_coordinates(mesh);
	  exio_helper.write_elements(mesh);
	  exio_helper.initialize_nodal_variables(names);
	}
    
      for (int c=0; c<num_vars; c++)
	{
	  std::vector<Number> cur_soln(num_nodes);

	  //Copy out this variable's solution
	  for(int i=0; i<num_nodes; i++)
	    cur_soln[i] = soln[i*num_vars + c];//c*num_nodes+i];
    
	  exio_helper.write_nodal_values(c+1,cur_soln,_timestep);
	}  
    }
}

#endif



#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::write_timestep (const std::string& ,
				  const EquationSystems& ,
				  const int ,
				  const double )
{

  std::cerr <<  "ERROR, ExodusII API is not defined.\n"
	    << std::endl;
  libmesh_error();
}

#else

void ExodusII_IO::write_timestep (const std::string& fname,
				  const EquationSystems& es,
				  const int timestep,
				  const double time)
{
  _timestep=timestep; 
  write_equation_systems(fname,es);

  if (libMesh::processor_id() == 0)
    exio_helper.write_timestep(timestep, time);
}

#endif



#ifndef LIBMESH_HAVE_EXODUS_API

void ExodusII_IO::write (const std::string& )
{
  std::cerr <<  "ERROR, ExodusII API is not defined.\n"
	    << std::endl;
  libmesh_error();
}


#else

void ExodusII_IO::write (const std::string& fname)
{
  if (libMesh::processor_id() == 0)
    {
      const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  
      libmesh_assert( !exio_helper.created() );

      exio_helper.create(fname);
      exio_helper.initialize(fname,mesh);
      exio_helper.write_nodal_coordinates(mesh);
      exio_helper.write_elements(mesh);

      // Note: the file is closed automatically by the ExodusII_IO destructor.
    }
}

#endif
