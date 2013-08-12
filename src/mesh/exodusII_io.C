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
#include "libmesh/nemesis_io_helper.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/parallel_mesh.h"

// Include the ParMETIS header files
namespace Parmetis {
    extern "C" {
#     include "libmesh/ignore_warnings.h"
#     include "parmetis.h"
#     include "libmesh/restore_warnings.h"
    }
}

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

void ExodusII_IO::set_coordinate_offset(Point p)
{
#ifdef LIBMESH_HAVE_EXODUS_API
  libmesh_deprecated();
  exio_helper->set_coordinate_offset(p);
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


    
    
    void ExodusII_IO::read_parallel (const std::string& base_filename)
    {
        // On one processor, Nemesis and ExodusII should be equivalent, so
        // let's cowardly defer to that implementation...
        if (this->n_processors() == 1)
        {
            // We can do this in one line but if the verbose flag was set in this
            // object, it will no longer be set... thus no extra print-outs for serial runs.
            // ExodusII_IO(this->mesh()).read (base_filename); // ambiguous when Nemesis_IO is multiply-inherited
            
            this->read (base_filename);
            return;
        }
        
        START_LOG ("read()","Exodus_IO");
        
        
        libMesh::out << "getting into read" << std::endl;
        
        // This function must be run on all processors at once
        parallel_object_only();
        
        // Open the Exodus file
        this->exio_helper->open(base_filename.c_str()); // just to avoid error from within this class
        ExodusII_IO_Helper ex_io_helper(this->comm(), true, false);
        ex_io_helper.open(base_filename.c_str());
        
        // Get a reference to the mesh.  We need to be specific
        // since Nemesis_IO is multiply-inherited
        // MeshBase& mesh = this->mesh();
        MeshBase& mesh = MeshInput<MeshBase>::mesh();
        
        // Local information: Read the following information from the standard Exodus header
        //  title[0]
        //  num_dim
        //  num_nodes
        //  num_elem
        //  num_elem_blk
        //  num_node_sets
        //  num_side_sets
        ex_io_helper.read_header();
        ex_io_helper.print_header();
        
        libMesh::out << "after header" << std::endl;
        
        ex_io_helper.read_block_info();

        exII::ex_get_elem_block(ex_io_helper.ex_id,
                                ex_io_helper.block_ids[0],
                                &ex_io_helper.elem_type[0],
                                &ex_io_helper.num_elem_this_blk,
                                &ex_io_helper.num_nodes_per_elem,
                                &ex_io_helper.num_attr);
        
        libMesh::out << "after read block info" << std::endl;
        
        //    // Get global information: number of nodes, elems, blocks, nodesets and sidesets
        //    ex_helper.get_init_global();
        
        // the approach is to partition the elements based on a space-filling approach.
        // So, the centroid information of each element is obtained using the nodal information.
        // First, each processor reads in its chunk of the elements
        
        // beginning and end of elem_ids on each processor
        std::vector<int> proc_elems(this->n_processors()+1, 0);
        int n_remaining_elems = ex_io_helper.num_elem,
        n_elems_per_proc = ex_io_helper.num_elem / this->n_processors();
        
        proc_elems[0] = 1; // exodus numbering starts from 1
        for (unsigned int i=0; i<this->n_processors(); i++)
        {
            proc_elems[i+1] = proc_elems[i] + std::min( n_elems_per_proc, n_remaining_elems)+1;
            n_remaining_elems -= (proc_elems[i+1]-proc_elems[i]);
        }
        
        proc_elems[this->n_processors()] += n_remaining_elems; // in case any elements were left unassigned
        
        unsigned int n_local_elems =
        proc_elems[this->processor_id()+1] - proc_elems[this->processor_id()];
        std::vector<char> elem_type(10);
        int n_elem_in_block=0, n_nodes_per_elem=0, n_attr=0;
        
        libMesh::out << "Elem range on proc: " << this->processor_id() << "  "
        << proc_elems[this->processor_id()] << "  " << proc_elems[this->processor_id()+1]
        << " n elem on proc: " << n_local_elems << std::endl;
        
        exII::ex_get_elem_block(ex_io_helper.ex_id, ex_io_helper.block_ids[0], &elem_type[0],
                                &n_elem_in_block,
                                &n_nodes_per_elem,
                                &n_attr);
        
        libMesh::out << "after block elem info" << std::endl;
        
        std::vector<int> elem_conn(n_local_elems * n_nodes_per_elem, 0);
        std::vector<float> elem_xyz(n_local_elems*ex_io_helper.num_dim, 0.);
        
        int err = Nemesis::ne_get_n_elem_conn(ex_io_helper.ex_id, ex_io_helper.block_ids[0], proc_elems[this->processor_id()],
                                              n_local_elems, &elem_conn[0]);
        
        libMesh::out << "after elem conn" << std::endl;
        
        // find the first and last nodes
        dof_id_type first_node=ex_io_helper.num_nodes, last_node=0;
        for (dof_id_type i=0; i<elem_conn.size(); i++)
        {
            if (first_node > elem_conn[i])
                first_node = elem_conn[i];
            if (last_node < elem_conn[i])
                last_node = elem_conn[i];
        }
        
        libMesh::out << "Node range on proc: " << this->processor_id() << "  " << first_node << "  " << last_node << std::endl;
        
        // now read in the detail for each node specified in the connectivity list for each element, and calculate
        // the centroid information
        // hopefully the size of these vectors will not be too huge
        std::vector<Real> node_x(last_node-first_node+1, 0.),
        node_y(last_node-first_node+1, 0.), node_z(last_node-first_node+1, 0.);
        err = Nemesis::ne_get_n_coord(ex_io_helper.ex_id, first_node, last_node-first_node+1,
                                      &node_x[0], &node_y[0], &node_z[0]);
        
        unsigned int node_pos_in_vector;
        for (unsigned int i_elem=0; i_elem<n_local_elems; i_elem++)
        {
            for (unsigned int i_node=0; i_node<n_nodes_per_elem; i_node++)
            {
                node_pos_in_vector = elem_conn[i_elem*n_nodes_per_elem+i_node]-first_node;
                elem_xyz[i_elem*ex_io_helper.num_dim+0] += node_x[node_pos_in_vector];
                elem_xyz[i_elem*ex_io_helper.num_dim+1] += node_y[node_pos_in_vector];
                elem_xyz[i_elem*ex_io_helper.num_dim+2] += node_z[node_pos_in_vector];
            }
            elem_xyz[i_elem*ex_io_helper.num_dim+0] /= n_nodes_per_elem;
            elem_xyz[i_elem*ex_io_helper.num_dim+1] /= n_nodes_per_elem;
            elem_xyz[i_elem*ex_io_helper.num_dim+2] /= n_nodes_per_elem;
        }
        
        // clear unneeded storage
        node_x.clear(); node_y.clear(); node_z.clear();
        
        // now call the partitioner
        std::vector<int> part(n_local_elems, 0);
        MPI_Comm mpi_comm = this->comm().get();
        
        libMesh::out << "getting into Parmetis" << std::endl;
        
        err = Parmetis::ParMETIS_V3_PartGeom(&proc_elems[0], &ex_io_helper.num_dim, &elem_xyz[0],
                                             &part[0], &mpi_comm);
        
        libMesh::out << "after Parmetis" << std::endl;
        
        // clear unneeded storage
        elem_xyz.clear();
        
        // now, this partitioning needs to be communicated to the respective processors
        std::vector<dof_id_type> collected_elems_on_proc, elems_for_comm;

        // add elements from the local processor
        unsigned int n_elems_for_comm = 0;
        elems_for_comm.resize(n_local_elems);
        
        // iterate over the elements locally and create the list of
        // elements that belong to proc dest
        for (unsigned int i=0; i<n_local_elems; i++)
            if (part[i] == this->processor_id())
            {
                elems_for_comm[n_elems_for_comm] = proc_elems[this->processor_id()]+i;
                n_elems_for_comm++;
            }
        collected_elems_on_proc.insert(collected_elems_on_proc.end(),
                                       elems_for_comm.begin(),
                                       elems_for_comm.begin()+n_elems_for_comm);

        
        for (processor_id_type pid=0; pid<mesh.n_processors(); pid++) // pid received data from others
        {
            if (pid == this->processor_id()) // receive from others
            {
                std::cout << "***** Receiving for pid: " << pid << std::endl;
                for (processor_id_type source=0; source<mesh.n_processors(); source++)
                    if (source != this->processor_id()) // don't send to self
                    {
                        elems_for_comm.clear();
                        this->comm().receive(source, elems_for_comm);
                        std::cout << "Received: " << elems_for_comm.size() << " elems from proc: " << source << std::endl;
                        collected_elems_on_proc.insert(collected_elems_on_proc.end(),
                                                       elems_for_comm.begin(), elems_for_comm.end());
                    }
                std::cout << "Total: " << collected_elems_on_proc.size() << " elems on proc: " << pid << std::endl;
            }
            else // send data to pid
            {
                std::cout << "***** Sending from pid: " << this->processor_id() << std::endl;
                elems_for_comm.resize(n_local_elems); // upper limit of the size for this vector

                n_elems_for_comm = 0;
                
                // iterate over the elements locally and create the list of
                // elements that belong to proc dest
                for (unsigned int i=0; i<n_local_elems; i++)
                    if (part[i] == pid)
                    {
                        elems_for_comm[n_elems_for_comm] = proc_elems[this->processor_id()]+i;
                        n_elems_for_comm++;
                    }
                
                std::vector<dof_id_type> elem_send_list;
                elem_send_list.insert(elem_send_list.end(),
                                      elems_for_comm.begin(),
                                      elems_for_comm.begin()+n_elems_for_comm);
                
                std::cout << "Sending: " << elem_send_list.size() << " elems to proc: " << pid << std::endl;
                this->comm().send(pid, elem_send_list);
            }
        }
        
        // clear the unneeded storage
        elems_for_comm.clear(); part.clear();

        // find the range of element ids
        unsigned int first_elem=ex_io_helper.num_elem, last_elem=0;
        for (std::vector<unsigned int>::const_iterator elem_it=collected_elems_on_proc.begin();
             elem_it != collected_elems_on_proc.end(); elem_it++)
        {
            if (*elem_it > last_elem)
                last_elem = *elem_it;
            if (*elem_it < first_elem)
                first_elem = *elem_it;
        }
        
        // now, get the node details for this element range
        n_local_elems = collected_elems_on_proc.size();
        elem_conn.resize((last_elem-first_elem+1) * n_nodes_per_elem); // this stores the connectivity for the entire range

        libMesh::out << "Elem range: " << first_elem << "  --  "  << last_elem << "  : with total elems:  " << n_local_elems << std::endl;
        

        libMesh::out << "Reading updated element connectivity " << std::endl;
        
        err = Nemesis::ne_get_n_elem_conn(ex_io_helper.ex_id, ex_io_helper.block_ids[0], first_elem,
                                          (last_elem-first_elem+1), &elem_conn[0]);

        libMesh::out << "Done reading updated element connectivity: Preparing node process ids" << std::endl;

        // find the node ownership and the range of node ids on this processor.
        // If a node lies on multiple processors, then the smallest processor id would take ownership of the node

        std::map<dof_id_type, processor_id_type> node_processor_id_map;

        unsigned int elem_offset, node_offset;
        for (std::vector<unsigned int>::const_iterator elem_it=collected_elems_on_proc.begin();
             elem_it != collected_elems_on_proc.end(); elem_it++)
        {
            elem_offset = (*elem_it-first_elem)*n_nodes_per_elem; // offset for connectivity data
            
            for (unsigned int j=0; j<n_nodes_per_elem; j++)
                // start by identify each node to be on this processor
                // this will be changed later
                node_processor_id_map[elem_conn[elem_offset+j]] = this->processor_id();
        }
        
        // get the first and last node ids from the map
        first_node = node_processor_id_map.begin()->first;
        last_node = node_processor_id_map.rbegin()->first;

        libMesh::out << "Done preparing node process ids: communicating IDs to processors" << std::endl;

        // now each processor communicates to the higher rank processors about the ownership
        for (processor_id_type pid=0; pid<mesh.n_processors(); pid++) // pid sends data to higher ranked processors
        {
            //prepare the node vector and send it to higher ranked processors
            if (pid == this->processor_id())
            {
                std::vector<dof_id_type> locally_owned_nodes(node_processor_id_map.size());
                dof_id_type index=0;
                for (std::map<dof_id_type, processor_id_type>::const_iterator map_it=node_processor_id_map.begin();
                     map_it != node_processor_id_map.end(); map_it++)
                {
                    if (map_it->second == pid) // pid == this->processor_id() here
                        locally_owned_nodes[index++] = map_it->first;
                }
                
                // send the data to the processor in a vector that is
                // sized for the number of locally owned nodes
                std::vector<dof_id_type> data_to_send;
                data_to_send.insert(data_to_send.end(), locally_owned_nodes.begin(), locally_owned_nodes.end());
                locally_owned_nodes.clear();
                
                // send this to all processors of higher rank
                std::cout << "Sending nodes to processors: n_nodes = " << data_to_send.size() << " : from pid = " <<  pid << std::endl;
                for (processor_id_type dest=pid+1; dest<mesh.n_processors(); dest++)
                    this->comm().send(dest, data_to_send);
            }
            else if (pid < this->processor_id()) // receive from lower ranked processors
            {
                // get the remote node ids
                std::vector<dof_id_type> remote_nodes;
                this->comm().receive(pid, remote_nodes);

                libMesh::out << "Received from : " << pid << " by " << this->processor_id() << " n_nodes : " << remote_nodes.size()  << std::endl;
                
                // now iterate over these nodes and set their processor ids
                std::map<dof_id_type, processor_id_type>::iterator map_it;
                std::map<dof_id_type, processor_id_type>::const_iterator map_end = node_processor_id_map.end();
                for (std::vector<dof_id_type>::const_iterator node_it=remote_nodes.begin();
                     node_it != remote_nodes.end(); node_it++)
                {
                    // check if the node also lies on this processor
                    map_it = node_processor_id_map.find(*node_it);
                    if (map_it != map_end)
                        map_it->second = pid;
                }
            }
        }
                
        libMesh::out << "Node range on proc: " << this->processor_id() << "  " << first_node << "  " << last_node << std::endl;
        
        // now read in the detail for each node specified in the connectivity list for each element, and calculate
        // the centroid information
        // hopefully the size of these vectors will not be too huge
        node_x.resize(last_node-first_node+1, 0.);
        node_y.resize(last_node-first_node+1, 0.);
        node_z.resize(last_node-first_node+1, 0.);
        
        libMesh::out << "Reading node coordinates " << std::endl;

        err = Nemesis::ne_get_n_coord(ex_io_helper.ex_id, first_node, last_node-first_node+1,
                                      &node_x[0], &node_y[0], &node_z[0]);

        libMesh::out << "Done reading node coordinates: Now adding nodes and elements to mesh " << std::endl;

        // add the nodes and elements
        ExodusII_IO_Helper::ElementMaps em;     // Instantiate the ElementMaps interface
        const std::string type_str (ex_io_helper.get_elem_type());
        const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(type_str);

        Node* node_ptr;
        std::map<dof_id_type, processor_id_type>::const_iterator map_it,
        map_end = node_processor_id_map.end();
        // Loop over all the elements in this block
        for (std::vector<unsigned int>::const_iterator elem_it=collected_elems_on_proc.begin();
             elem_it != collected_elems_on_proc.end(); elem_it++)
        {
            elem_offset = (*elem_it-first_elem)*n_nodes_per_elem; // offset for connectivity data

            Elem* elem = Elem::build (conv.get_canonical_type()).release(); // create the element
            elem->processor_id(this->processor_id()); // only add locally
            elem->set_id(*elem_it); // prescribe the elem id
            
            for (unsigned int j=0; j<n_nodes_per_elem; j++)
            {
                node_ptr = mesh.query_node_ptr(elem_conn[elem_offset+j]);
                if (node_ptr == NULL)
                {
                    node_offset = elem_conn[elem_offset+j] - first_node;
                    
                    // local map should certainly have this node id
                    map_it = node_processor_id_map.find(elem_conn[elem_offset+j]);
                    libmesh_assert(map_it != map_end);
                    
                    // use the node id, add processor id identified earlies
                    node_ptr = mesh.add_point(Point(node_x[node_offset],
                                                    node_y[node_offset],
                                                    node_z[node_offset]),
                                              elem_conn[elem_offset+j],
                                              map_it->second);
                }
                
                elem->set_node(j) = node_ptr;
            }

            mesh.add_elem(elem);
        }
        
        // clear unneeded storage
        collected_elems_on_proc.clear(); elem_conn.clear();
        node_x.clear(); node_y.clear(); node_z.clear();
        node_processor_id_map.clear();
        

        unsigned int n_side_sets = 0;
        // Read in sideset information -- this is useful for applying boundary conditions
        {
            ex_io_helper.read_sideset_info(); // Get basic information about ALL sidesets
            int offset=0;
            for (int i=0; i<ex_io_helper.get_num_side_sets(); i++)
            {
                offset += (i > 0 ? ex_io_helper.get_num_sides_per_set(i-1) : 0); // Compute new offset
                ex_io_helper.read_sideset (i, offset);
                
                mesh.boundary_info->sideset_name(ex_io_helper.get_side_set_id(i)) =
                ex_io_helper.get_side_set_name(i);
            }
            
            const std::vector<int>& elem_list = ex_io_helper.get_elem_list();
            const std::vector<int>& side_list = ex_io_helper.get_side_list();
            const std::vector<int>& id_list   = ex_io_helper.get_id_list();
            
            for (unsigned int e=0; e<elem_list.size(); e++)
            {
                // Set any relevant node/edge maps for this element
                
                Elem * elem = mesh.query_elem(elem_list[e]);
                
                if (elem != NULL) // proceed only if this processor contains this elemid
                {
                    
                    const ExodusII_IO_Helper::Conversion conv =
                    em.assign_conversion(elem->type());
                    
                    mesh.boundary_info->add_side (elem_list[e],
                                                  conv.get_side_map(side_list[e]-1),
                                                  id_list[e]);
                    n_side_sets++;
                }
            }
        }
        
        
        std::cout << "Done adding side set: on pid: " << this->processor_id() << " : n_side_sets : " << n_side_sets << std::endl;
        this->comm().sum(n_side_sets);
        std::cout << "Total side sets: " << n_side_sets << std::endl;
        
        // Read nodeset info
        {
            ex_io_helper.read_nodeset_info();
            
            for (int nodeset=0; nodeset<ex_io_helper.get_num_node_sets(); nodeset++)
            {
                int nodeset_id = ex_io_helper.get_nodeset_id(nodeset);
                
                mesh.boundary_info->nodeset_name(nodeset_id) =
                ex_io_helper.get_node_set_name(nodeset);
                
                ex_io_helper.read_nodeset(nodeset);
                
                const std::vector<int>& node_list = ex_io_helper.get_node_list();
                
                for(unsigned int node=0; node<node_list.size(); node++)
                {
                    Node * node_ptr = mesh.query_node_ptr(node_list[node]);
                    if (node_ptr != NULL)
                        mesh.boundary_info->add_node(node_list[node], nodeset_id);
                }
            }
        }

        
        
        libMesh::out << "Done adding elements to mesh: Now preparing for use " << std::endl;
        
        err = exII::ex_close(ex_io_helper.ex_id);

        
        // For ParallelMesh, it seems that _is_serial is true by default.  A hack to
        // make the Mesh think it's parallel might be to call:
        mesh.update_post_partitioning();
        mesh.delete_remote_elements();

        // now prepare for use
        MeshCommunication().gather_neighboring_elements(libmesh_cast_ref<ParallelMesh&>(mesh));
        
        libMesh::out << "Done " << std::endl;

        STOP_LOG ("read()","Exodus_IO");
        
        return;
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

void ExodusII_IO::copy_elemental_solution(System&, std::string, std::string, unsigned int)
{

  libMesh::err <<  "ERROR, ExodusII API is not defined.\n"
            << std::endl;
  libmesh_error();
}

#else

void ExodusII_IO::copy_elemental_solution(System& system, std::string system_var_name, std::string exodus_var_name, unsigned int timestep)
{
  // FIXME: Do we need to call get_time_steps() at all?
  /*const std::vector<double>& time_steps = */
  exio_helper->get_time_steps();

  const std::vector<Real> & elemental_values = exio_helper->get_elemental_var_values(exodus_var_name, timestep);

  const unsigned int var_num = system.variable_number(system_var_name);
  if (system.variable_type(var_num) != FEType(CONSTANT, MONOMIAL))
  {
    libMesh::err << "Error! Trying to copy elemental solution into a variable that is not of CONSTANT MONOMIAL type. " << std::endl;
    libmesh_error();
  }

  for (unsigned int i=0; i<elemental_values.size(); ++i)
    {
      const Elem * elem = MeshInput<MeshBase>::mesh().elem(i);

      if (!elem)
        {
          libMesh::err << "Error! Mesh returned NULL pointer for elem " << i << std::endl;
          libmesh_error();
        }

      dof_id_type dof_index = elem->dof_number(system.number(), var_num, 0);

      // If the dof_index is local to this processor, set the value
      if ((dof_index >= system.solution->first_local_index()) && (dof_index < system.solution->last_local_index()))
    system.solution->set (dof_index, elemental_values[i]);
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
