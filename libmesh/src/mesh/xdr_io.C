// $Id: legacy_xdr_io.C 2560 2007-12-03 17:52:20Z benkirk $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include <iostream>
#include <iomanip>

#include <vector>
#include <string>

// Local includes
#include "xdr_io.h"
#include "legacy_xdr_io.h"
#include "xdr_cxx.h"
#include "enum_xdr_mode.h"
#include "mesh_base.h"
#include "node.h"
#include "elem.h"
#include "boundary_info.h"
#include "parallel.h"
#include "mesh_tools.h"






// ------------------------------------------------------------
// XdrIO static data
const unsigned int XdrIO::io_blksize = 2;



// ------------------------------------------------------------
// XdrIO members
XdrIO::XdrIO (MeshBase& mesh, const bool binary) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _binary             (binary),
  _legacy             (true),
  _version            ("libMesh-0.7.0+"),
  _bc_file_name       ("."),
  _partition_map_file ("n/a"),
  _subdomain_map_file ("n/a"),
  _p_level_file       ("n/a")
{
}



XdrIO::XdrIO (const MeshBase& mesh, const bool binary) :
  MeshOutput<MeshBase>(mesh),
  _binary (binary)
{
}



XdrIO::~XdrIO ()
{
}



void XdrIO::read (const std::string& name)
{
  Xdr io (name, this->binary() ? DECODE : READ);

  // get the version string
  io.data(this->version());
  
  this->legacy() = !(this->version().find("libMesh") < this->version().size());

  std::cout << "version=" << this->version() << std::endl;

  // Check for a legacy version format.
  if (this->legacy())
    {
      io.close();
      LegacyXdrIO(MeshInput<MeshBase>::mesh(), this->binary()).read(name);
      return;
    }

  unsigned int n_elem, n_nodes;
  io.data (n_elem);
  io.data (n_nodes);
  
  io.data (this->boundary_condition_file_name()); std::cout << "bc_file="  << this->boundary_condition_file_name() << std::endl;
  io.data (this->subdomain_map_file_name());      std::cout << "sid_file=" << this->subdomain_map_file_name()      << std::endl;
  io.data (this->partition_map_file_name());      std::cout << "pid_file=" << this->partition_map_file_name()      << std::endl;
  io.data (this->polynomial_level_file_name());   std::cout << "pl_file="  << this->polynomial_level_file_name()  << std::endl;
}



void XdrIO::write (const std::string& name)
{
  if (this->legacy())
    {
      LegacyXdrIO(MeshOutput<MeshBase>::mesh(), this->binary()).write(name);
      return;
    }
  
  Xdr io ((libMesh::processor_id() == 0) ? name : "", this->binary() ? ENCODE : WRITE);

  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  
  unsigned int
    n_elem  = mesh.n_elem(),
    n_nodes = mesh.n_nodes(),
    n_bcs   = mesh.boundary_info->n_boundary_conds();

  // If there are no bcs, don't write anything
  if (!n_bcs)
    this->boundary_condition_file_name() = "n/a";

  // If there are more than one subdomains and the user has not specified an 
  // external file then write the subdomain mapping to the default file "."  
  // Do not assume that we will not write a sudomain file if there is only 1 
  // subdomain, since reading such a mesh would break in the case when all the 
  // elements are on subdomain 100, for example.  In this case, however, the 
  // user will need to explicitly set the subdomain_map_file_name() before 
  // calling the write() method.
  if ((mesh.n_subdomains() > 1) && 
      (this->subdomain_map_file_name() == "n/a"))
    this->subdomain_map_file_name() = ".";

  // write the header
  if (libMesh::processor_id() == 0)
    {
      io.data (this->version());
      
      io.data (n_elem,  "# number of elements");
      io.data (n_nodes, "# number of nodes");
      
      io.data (this->boundary_condition_file_name(), "# boundary condition specification file");
      io.data (this->subdomain_map_file_name(),      "# subdomain id specification file");
      io.data (this->partition_map_file_name(),      "# processor id specification file");
      io.data (this->polynomial_level_file_name(),   "# p-level specification file");      
    }

  // write connectivity
  this->write_serialized_connectivity (io, n_elem);

  // write the nodal locations
  this->write_serialized_nodes (io, n_nodes);
  
  if (libMesh::processor_id() == 0)
    io.data (n_bcs,   "# number of boundary conditions");
}



void XdrIO::write_serialized_connectivity (Xdr &io, const unsigned int n_elem)
{
  assert (io.writing());

  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  assert (n_elem == mesh.n_elem());

  // We will only write active elements and their parents.
  const unsigned int n_active_levels = MeshTools::n_active_levels (mesh);
  std::vector<unsigned int> 
    n_global_elem_at_level(n_active_levels), n_local_elem_at_level(n_active_levels);

  MeshBase::const_element_iterator it  = mesh.local_elements_end(), end=it;

  // Find the number of local and global elements at each level
  for (unsigned int level=0, tot_n_elem=0; level<n_active_levels; level++)
    {
      it  = mesh.local_level_elements_begin(level);
      end = mesh.local_level_elements_end(level);

      n_local_elem_at_level[level] = n_global_elem_at_level[level] = MeshTools::n_elem(it, end);

      Parallel::sum(n_global_elem_at_level[level]);
      tot_n_elem += n_global_elem_at_level[level];
      assert (n_global_elem_at_level[level] <= n_elem);
      assert (tot_n_elem <= n_elem);
    }

  std::vector<unsigned int> 
    xfer_conn, recv_conn, output_buffer, 
    elem_offset(libMesh::n_processors()), 
    xfer_buf_sizes(libMesh::n_processors());
  std::map<unsigned int, unsigned int> parent_id_map, child_id_map;
  unsigned int my_next_elem=0, next_global_elem=0;

  //-------------------------------------------
  // First write the level-0 elements directly.
  it  = mesh.local_level_elements_begin(0);
  end = mesh.local_level_elements_end(0);
  for (; it != end; ++it)
    {
      pack_element (xfer_conn, *it);
      parent_id_map[(*it)->id()] = my_next_elem++;
    }
  xfer_conn.push_back(my_next_elem); // toss in the number of elements transferred.

  unsigned int my_size = xfer_conn.size();
  Parallel::gather (0, my_next_elem, elem_offset);
  Parallel::gather (0, my_size,      xfer_buf_sizes);

  // All elements send their xfer buffers to processor 0.
  Parallel::request request_handle;
  Parallel::isend (0, xfer_conn, request_handle);
  
  // Processor 0 will receive the data and write out the elements. 
  if (libMesh::processor_id() == 0)
    for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
      {
	recv_conn.resize(xfer_buf_sizes[pid]);
	Parallel::recv (pid, recv_conn);	
	// at a minimum, the buffer should contain the number of elements,
	// which could be 0.
	assert (!recv_conn.empty());
	
	const unsigned int n_elem_received = recv_conn.back();
	std::vector<unsigned int>::const_iterator it = recv_conn.begin();

	for (unsigned int elem=0; elem<n_elem_received; elem++, next_global_elem++)
	  {
	    output_buffer.clear();
	    const unsigned int n_nodes = *it; ++it;
	    output_buffer.push_back(*it);     /* type          */ ++it;
	    /*output_buffer.push_back(*it);*/ /* id            */ ++it;
	    /*output_buffer.push_back(*it);*/ /* subdomain id  */ ++it;
	    /*output_buffer.push_back(*it);*/ /* p level       */ ++it;

	    for (unsigned int node=0; node<n_nodes; node++, ++it)
	      output_buffer.push_back(*it);
	    
	    io.data_stream (&output_buffer[0], output_buffer.size(), output_buffer.size());
	  }       	
      }      
#ifdef HAVE_MPI
  MPI_Wait (&request_handle, MPI_STATUS_IGNORE);
#endif

  //--------------------------------------------------------------------
  // Next write the remaining elements indirectly through their parents.
  // This will insure that the children are written in the proper order
  // so they can be reconstructed properly.
  for (unsigned int level=1, my_elem_at_level=0, processor_offset=0; level<n_active_levels; level++)
    {
      xfer_conn.clear();

      it  = mesh.local_level_elements_begin(level-1);
      end = mesh.local_level_elements_end  (level-1);
      
      for (my_elem_at_level=0; it != end; ++it)
	if (!(*it)->active()) // we only want the parents elements at this level, and
	  {                   // there is no direct iterator for this obscure use
	    const Elem *parent = *it;
	    assert (parent_id_map.count(parent->id()));
	    const unsigned int parent_id = parent_id_map[parent->id()];
	    parent_id_map.erase(parent->id());

	    for (unsigned int c=0; c<parent->n_children(); c++, my_next_elem++)
	      {
		const Elem *child = parent->child(c);
		pack_element (xfer_conn, child, parent_id);
		child_id_map[child->id()] = my_elem_at_level++;
	      }
	  }
      xfer_conn.push_back(my_elem_at_level);
      my_size = xfer_conn.size();
      Parallel::gather (0, my_size, xfer_buf_sizes);
      Parallel::isend  (0, xfer_conn, request_handle);
      
      parent_id_map = child_id_map; /**/ child_id_map.clear();
      
      // Processor 0 will receive the data and write the elements.
      if (libMesh::processor_id() == 0)
	for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	  {
	    recv_conn.resize(xfer_buf_sizes[pid]);
	    Parallel::recv (pid, recv_conn);	
	    // at a minimum, the buffer should contain the number of elements,
	    // which could be 0.
	    assert (!recv_conn.empty());
	    const unsigned int n_elem_received = recv_conn.back();
	    std::vector<unsigned int>::const_iterator it = recv_conn.begin();

	    for (unsigned int elem=0; elem<n_elem_received; elem++, next_global_elem++)
	      {
		output_buffer.clear();
		const unsigned int n_nodes = *it; ++it;
		output_buffer.push_back(*it);                   /* type          */ ++it;
		/*output_buffer.push_back(*it);*/               /* id            */ ++it;
		output_buffer.push_back (*it+processor_offset); /* parent id     */ ++it;
		/*output_buffer.push_back(*it);*/               /* subdomain id  */ ++it;
		/*output_buffer.push_back(*it);*/               /* p level       */ ++it;
		
		for (unsigned int node=0; node<n_nodes; node++, ++it)
		  output_buffer.push_back(*it);
	    
		io.data_stream (&output_buffer[0], output_buffer.size(), output_buffer.size());
	      }       	
	    
	    processor_offset += elem_offset[pid];
	  }
#ifdef HAVE_MPI
      MPI_Wait (&request_handle, MPI_STATUS_IGNORE);
#endif
    }
  if (libMesh::processor_id() == 0)
    assert (next_global_elem == n_elem);
  
}



void XdrIO::write_serialized_nodes (Xdr &io, const unsigned int n_nodes)
{
  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  assert (n_nodes == mesh.n_nodes());

  std::vector<unsigned int> xfer_ids;
  std::vector<Real>         xfer_coords, &coords=xfer_coords;

  std::vector<std::vector<unsigned int> > recv_ids   (libMesh::n_processors());;
  std::vector<std::vector<Real> >         recv_coords(libMesh::n_processors());

  unsigned int n_written=0;

  for (unsigned int blk=0, last_node=0; last_node<n_nodes; blk++)
    {
      std::cout << "Writing node block " << blk << std::endl;

      const unsigned int first_node = blk*io_blksize;
                          last_node = std::min((blk+1)*io_blksize, n_nodes);

      // Build up the xfer buffers on each processor
      MeshBase::const_node_iterator
	it  = mesh.local_nodes_begin(),
	end = mesh.local_nodes_end();

      xfer_ids.clear(); xfer_coords.clear();
      
      for (; it!=end; ++it)
	if (((*it)->id() >= first_node) && // node in [first_node, last_node)
	    ((*it)->id() <  last_node))
	  {
	    xfer_ids.push_back((*it)->id());
	    const Point &p = **it;
	    xfer_coords.push_back(p(0));
	    xfer_coords.push_back(p(1));
	    xfer_coords.push_back(p(2));
	  }

      //-------------------------------------
      // Send the xfer buffers to processor 0
      std::vector<unsigned int> ids_size, coords_size;

      const unsigned int my_ids_size = xfer_ids.size();

      // explicitly gather ids_size
      Parallel::gather (0, my_ids_size, ids_size);

      // infer coords_size on processor 0
      if (libMesh::processor_id() == 0)
	{
	  coords_size.reserve(libMesh::n_processors());
	  for (unsigned int p=0; p<ids_size.size(); p++)
	    coords_size.push_back(3*ids_size[p]);
	}			  

      // Note that we will actually send/receive to ourself if we are
      // processor 0, so let's use nonblocking receives.
      std::vector<Parallel::request>
        id_request_handles(libMesh::n_processors()),
        coord_request_handles(libMesh::n_processors());

#ifdef HAVE_MPI
      const unsigned int id_tag=0, coord_tag=1;
      
      // Post the receives -- do this on processor 0 only.
      if (libMesh::processor_id() == 0)
        for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
          {
            recv_ids[pid].resize(ids_size[pid]);
            recv_coords[pid].resize(coords_size[pid]);
            
            Parallel::irecv (pid, recv_ids[pid],    id_request_handles[pid],    id_tag);
            Parallel::irecv (pid, recv_coords[pid], coord_request_handles[pid], coord_tag);
          }
      
      // Send -- do this on all processors.
      Parallel::send(0, xfer_ids,    id_tag);
      Parallel::send(0, xfer_coords, coord_tag);
#else
      // On one processor there's nothing to send
      recv_ids[0]    = xfer_ids;
      recv_coords[0] = xfer_coords;
#endif

      // -------------------------------------------------------
      // Receive the messages and write the output on processor 0.
      if (libMesh::processor_id() == 0)
        {
#ifdef HAVE_MPI
          // Wait for all the receives to complete. We have no
          // need for the statuses since we already know the
          // buffer sizes.
          MPI_Waitall (libMesh::n_processors(),
                       &id_request_handles[0],
                       MPI_STATUSES_IGNORE);
          MPI_Waitall (libMesh::n_processors(),
                       &coord_request_handles[0],
                       MPI_STATUSES_IGNORE);
#endif
          
          // Write the coordinates in this block.
	  unsigned int tot_id_size=0, tot_coord_size=0;
	  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	    {
	      tot_id_size    += recv_ids[pid].size();
	      tot_coord_size += recv_coords[pid].size();
	    }

	  assert (tot_id_size <= std::min(io_blksize, n_nodes));
	  assert (tot_coord_size == 3*tot_id_size);

	  coords.resize (tot_coord_size);
	  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	    for (unsigned int idx=0; idx<recv_ids[pid].size(); idx++)
	      {
		const unsigned int local_idx = recv_ids[pid][idx] - first_node;
		assert ((3*local_idx+2) < coords.size());
		assert ((3*idx+2)      < recv_coords[pid].size());
		
		coords[3*local_idx+0] = recv_coords[pid][3*idx+0];
		coords[3*local_idx+1] = recv_coords[pid][3*idx+1];
		coords[3*local_idx+2] = recv_coords[pid][3*idx+2];

		n_written++;
	      }

	  io.data_stream (coords.empty() ? NULL : &coords[0], coords.size(), 3);
	}
    }
  if (libMesh::processor_id() == 0)
    assert (n_written == n_nodes);
}



void XdrIO::pack_element (std::vector<unsigned int> &conn, const Elem *elem, const unsigned int parent_id) const
{
  assert (elem != NULL);

  conn.push_back(elem->n_nodes());

  conn.push_back (elem->type());
  conn.push_back (elem->id());
  if (parent_id != libMesh::invalid_uint)
    conn.push_back (parent_id);
  conn.push_back (elem->subdomain_id());
  
#ifdef ENABLE_AMR
  conn.push_back (elem->p_level());
#endif
  
  for (unsigned int n=0; n<elem->n_nodes(); n++)
    conn.push_back (elem->node(n));                 
}
