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
#include "partitioner.h"
#include "libmesh_logging.h"


//-----------------------------------------------
// anonymous namespace for implementation details
namespace {
  struct ElemBCData
  {
    unsigned int       elem_id;
    unsigned short int side;
    short int          bc_id;
    
    // Default constructor
    ElemBCData (unsigned int       elem_id_in=0,
		unsigned short int side_in=0,
		short int          bc_id_in=0) :
      elem_id(elem_id_in),
      side(side_in),
      bc_id(bc_id_in)
    {}

    // comparison operator
    bool operator < (const ElemBCData &other) const
    {
      if (this->elem_id == other.elem_id)
	return (this->side < other.side);
      
      return this->elem_id < other.elem_id;
    }  
  };

  // comparison operator
  bool operator < (const unsigned int &other_elem_id,
		   const ElemBCData &elem_bc)
  {
    return other_elem_id < elem_bc.elem_id;
  }

  bool operator < (const ElemBCData &elem_bc,
		   const unsigned int &other_elem_id)
    
  {
    return elem_bc.elem_id < other_elem_id;
  }


  // For some reason SunStudio does not seem to accept the above
  // comparison functions for use in
  // std::equal_range (ElemBCData::iterator, ElemBCData::iterator, unsigned int);
  struct CompareIntElemBCData
  {
    bool operator()(const unsigned int &other_elem_id,
		    const ElemBCData &elem_bc)
    {
      return other_elem_id < elem_bc.elem_id;
    }

    bool operator()(const ElemBCData &elem_bc,
		    const unsigned int &other_elem_id)    
    {
      return elem_bc.elem_id < other_elem_id;
    }
  };
}



// ------------------------------------------------------------
// XdrIO static data
const unsigned int XdrIO::io_blksize = 128000;



// ------------------------------------------------------------
// XdrIO members
XdrIO::XdrIO (MeshBase& mesh, const bool binary) :
  MeshInput<MeshBase> (mesh,/* is_parallel_format = */ true),
  MeshOutput<MeshBase>(mesh,/* is_parallel_format = */ true),
  _binary             (binary),
  _legacy             (false),
  _version            ("libMesh-0.7.0+"),
  _bc_file_name       ("n/a"),
  _partition_map_file ("n/a"),
  _subdomain_map_file ("n/a"),
  _p_level_file       ("n/a")
{
}



XdrIO::XdrIO (const MeshBase& mesh, const bool binary) :
  MeshOutput<MeshBase>(mesh,/* is_parallel_format = */ true),
  _binary (binary)
{
}



XdrIO::~XdrIO ()
{
}



void XdrIO::write (const std::string& name)
{
  if (this->legacy())
    {
      deprecated();
      LegacyXdrIO(MeshOutput<MeshBase>::mesh(), this->binary()).write(name);
      return;
    }
  
  Xdr io ((libMesh::processor_id() == 0) ? name : "", this->binary() ? ENCODE : WRITE);

  START_LOG("write()","XdrIO");
  
  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  
  unsigned int
    n_elem     = mesh.n_elem(),
    n_nodes    = mesh.n_nodes(),
    n_bcs      = mesh.boundary_info->n_boundary_conds(),
    n_p_levels = MeshTools::n_p_levels (mesh); 


  //-------------------------------------------------------------
  // For all the optional files -- the default file name is "n/a".
  // However, the user may specify an optional external file. 
  
  // If there are BCs and the user has not already provided a
  // file name then write to "."
  if (n_bcs &&
      this->boundary_condition_file_name() == "n/a")
    this->boundary_condition_file_name() = ".";
  
  // If there are more than one subdomains and the user has not specified an 
  // external file then write the subdomain mapping to the default file "."  
  if ((mesh.n_subdomains() > 0) && 
      (this->subdomain_map_file_name() == "n/a"))
    this->subdomain_map_file_name() = ".";

  // In general we don't write the partition information.
  
  // If we have p levels and the user has not already provided
  // a file name then write to "."
  if ((n_p_levels > 1) &&
      (this->polynomial_level_file_name() == "n/a"))
    this->polynomial_level_file_name() = ".";

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

  // write the boundary condition information
  this->write_serialized_bcs (io, n_bcs);

  STOP_LOG("write()","XdrIO");
  
  // pause all processes until the write ends -- this will 
  // protect for the pathological case where a write is 
  // followed immediately by a read.  The write must be
  // guaranteed to complete first.
  io.close();
  Parallel::barrier();
}



void XdrIO::write_serialized_connectivity (Xdr &io, const unsigned int n_elem) const
{
  libmesh_assert (io.writing());
  
  const bool
    write_p_level      = ("." == this->polynomial_level_file_name()),
    write_partitioning = ("." == this->partition_map_file_name()),
    write_subdomain_id = ("." == this->subdomain_map_file_name());

  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  libmesh_assert (n_elem == mesh.n_elem());

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
      libmesh_assert (n_global_elem_at_level[level] <= n_elem);
      libmesh_assert (tot_n_elem <= n_elem);
    }

  std::vector<unsigned int> 
    xfer_conn, recv_conn, output_buffer, 
    n_elem_on_proc(libMesh::n_processors()), 
    processor_offsets(libMesh::n_processors()), 
    xfer_buf_sizes(libMesh::n_processors());

  typedef std::map<unsigned int, std::pair<unsigned int, unsigned int> > id_map_type;
  id_map_type parent_id_map, child_id_map;

  unsigned int my_next_elem=0, next_global_elem=0;

  //-------------------------------------------
  // First write the level-0 elements directly.
  it  = mesh.local_level_elements_begin(0);
  end = mesh.local_level_elements_end(0);
  for (; it != end; ++it)
    {
      pack_element (xfer_conn, *it);
      parent_id_map[(*it)->id()] = std::make_pair(libMesh::processor_id(), 
						  my_next_elem++);
    }
  xfer_conn.push_back(my_next_elem); // toss in the number of elements transferred.

  unsigned int my_size = xfer_conn.size();
  Parallel::gather (0, my_next_elem, n_elem_on_proc);
  Parallel::gather (0, my_size,      xfer_buf_sizes);

  processor_offsets[0] = 0;
  for (unsigned int pid=1; pid<libMesh::n_processors(); pid++)
    processor_offsets[pid] = processor_offsets[pid-1] + n_elem_on_proc[pid-1];

  // All processors send their xfer buffers to processor 0.
  // Processor 0 will receive the data and write out the elements. 
  if (libMesh::processor_id() == 0)
    {
      // Write the number of elements at this level.      
      {
	std::string comment = "# n_elem at level 0", legend  = ", [ type ";
	if (write_partitioning)
	  legend += "pid ";
	if (write_subdomain_id)
	  legend += "sid ";
	if (write_p_level)
	  legend += "p_level ";
	legend += "(n0 ... nN-1) ]";
	comment += legend;
	io.data (n_global_elem_at_level[0], comment.c_str());
      }

      for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	{
	  recv_conn.resize(xfer_buf_sizes[pid]);
          if (pid == 0)
	    recv_conn = xfer_conn;
          else
	    Parallel::receive (pid, recv_conn);	

	  // at a minimum, the buffer should contain the number of elements,
	  // which could be 0.
	  libmesh_assert (!recv_conn.empty());
	  
	  const unsigned int n_elem_received = recv_conn.back();
	  std::vector<unsigned int>::const_iterator it = recv_conn.begin();
	  
	  for (unsigned int elem=0; elem<n_elem_received; elem++, next_global_elem++)
	    {
	      output_buffer.clear();
	      const unsigned int n_nodes = *it; ++it;
	      output_buffer.push_back(*it);     /* type       */ ++it;	      
	      /*output_buffer.push_back(*it);*/ /* id         */ ++it;
	      
	      if (write_partitioning)
		output_buffer.push_back(*it); /* processor id */ ++it;

	      if (write_subdomain_id)
		output_buffer.push_back(*it); /* subdomain id */ ++it;

#ifdef LIBMESH_ENABLE_AMR
	      if (write_p_level)
		output_buffer.push_back(*it); /* p level      */ ++it;
#endif
	      for (unsigned int node=0; node<n_nodes; node++, ++it)
		output_buffer.push_back(*it);
	      
	      io.data_stream (&output_buffer[0], output_buffer.size(), output_buffer.size());
	    }       	
	}
    }
  else
    Parallel::send (0, xfer_conn);

#ifdef LIBMESH_ENABLE_AMR  
  //--------------------------------------------------------------------
  // Next write the remaining elements indirectly through their parents.
  // This will insure that the children are written in the proper order
  // so they can be reconstructed properly.
  for (unsigned int level=1, my_n_elem_written_at_level=0; level<n_active_levels; level++)
    {
      xfer_conn.clear();

      it  = mesh.local_level_elements_begin(level-1);
      end = mesh.local_level_elements_end  (level-1);
      
      for (my_n_elem_written_at_level=0; it != end; ++it)
	if (!(*it)->active()) // we only want the parents elements at this level, and
	  {                   // there is no direct iterator for this obscure use
	    const Elem *parent = *it;
	    id_map_type::iterator pos = parent_id_map.find(parent->id());
	    libmesh_assert (pos != parent_id_map.end());
	    const unsigned int parent_pid = pos->second.first;
	    const unsigned int parent_id  = pos->second.second;
	    parent_id_map.erase(pos);

	    for (unsigned int c=0; c<parent->n_children(); c++, my_next_elem++)
	      {
		const Elem *child = parent->child(c);
		pack_element (xfer_conn, child, parent_id, parent_pid);

		// this aproach introduces the possibility that we write
		// non-local elements.  These elements may well be parents
		// at the next step
		child_id_map[child->id()] = std::make_pair (child->processor_id(), 
							    my_n_elem_written_at_level++);
	      }
	  }
      xfer_conn.push_back(my_n_elem_written_at_level);
      my_size = xfer_conn.size();
      Parallel::gather (0, my_size,   xfer_buf_sizes);
      
      // Processor 0 will receive the data and write the elements.
      if (libMesh::processor_id() == 0)
	{
	  // Write the number of elements at this level.      
	  {
	    char buf[80];
	    std::sprintf(buf, "# n_elem at level %d", level);
	    std::string comment(buf), legend  = ", [ type parent ";

	    if (write_partitioning)
	      legend += "pid ";
	    if (write_subdomain_id)
	      legend += "sid ";
	    if (write_p_level)
	      legend += "p_level ";
	    legend += "(n0 ... nN-1) ]";
	    comment += legend;
	    io.data (n_global_elem_at_level[level], comment.c_str());
	  }
	  
	  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	    {
	      recv_conn.resize(xfer_buf_sizes[pid]);
              if (pid == 0)
	        recv_conn = xfer_conn;
              else
	        Parallel::receive (pid, recv_conn);	

	      // at a minimum, the buffer should contain the number of elements,
	      // which could be 0.
	      libmesh_assert (!recv_conn.empty());
	      const unsigned int n_elem_received = recv_conn.back();
	      std::vector<unsigned int>::const_iterator it = recv_conn.begin();
	      
	      for (unsigned int elem=0; elem<n_elem_received; elem++, next_global_elem++)
		{
		  output_buffer.clear();
		  const unsigned int n_nodes = *it; ++it;
		  output_buffer.push_back(*it);                   /* type          */ ++it;
		  /*output_buffer.push_back(*it);*/               /* id            */ ++it;
		  
		  const unsigned int parent_local_id = *it; ++it; 
		  const unsigned int parent_pid      = *it; ++it; 
	      	      
		  output_buffer.push_back (parent_local_id+processor_offsets[parent_pid]);

		  if (write_partitioning)
		    output_buffer.push_back(*it); /* processor id */ ++it;

		  if (write_subdomain_id)
		    output_buffer.push_back(*it); /* subdomain id */ ++it;
		  
		  if (write_p_level)
		    output_buffer.push_back(*it); /* p level       */ ++it;
		  
		  for (unsigned int node=0; node<n_nodes; node++, ++it)
		    output_buffer.push_back(*it);
		  
		  io.data_stream (&output_buffer[0], output_buffer.size(), output_buffer.size());
		}       		    
	    }
	}  
      else
        Parallel::send  (0, xfer_conn);
  
      // update the processor_offsets
      processor_offsets[0] = processor_offsets.back() + n_elem_on_proc.back();	  
      Parallel::gather (0, my_n_elem_written_at_level, n_elem_on_proc);
      for (unsigned int pid=1; pid<libMesh::n_processors(); pid++)
	processor_offsets[pid] = processor_offsets[pid-1] + n_elem_on_proc[pid-1];	    

      // Now, at the next level we will again iterate over local parents.  However, 
      // those parents may have been written by other processors (at this step), 
      // so we need to gather them into our *_id_maps.
      {
	std::vector<std::vector<unsigned int> > requested_ids(libMesh::n_processors());
	std::vector<unsigned int> request_to_fill;	

	it  = mesh.local_level_elements_begin(level);
	end = mesh.local_level_elements_end(level);

	for (; it!=end; ++it)
	  if (!child_id_map.count((*it)->id()))
	    {
	      libmesh_assert ((*it)->parent()->processor_id() != libMesh::processor_id());
	      requested_ids[(*it)->parent()->processor_id()].push_back((*it)->id());
	    }

	// Next set the child_ids 
	for (unsigned int p=1; p != libMesh::n_processors(); ++p)
	  {
	    // Trade my requests with processor procup and procdown
	    unsigned int procup = (libMesh::processor_id() + p) %
                                   libMesh::n_processors();
	    unsigned int procdown = (libMesh::n_processors() +
                                     libMesh::processor_id() - p) %
	                             libMesh::n_processors();

	    Parallel::send_receive(procup, requested_ids[procup],
				   procdown, request_to_fill);

	    // Fill those requests by overwriting the requested ids
	    for (unsigned int i=0; i<request_to_fill.size(); i++)	   
	      {
		libmesh_assert (child_id_map.count(request_to_fill[i]));
		libmesh_assert (child_id_map[request_to_fill[i]].first == procdown);

		request_to_fill[i] = child_id_map[request_to_fill[i]].second;
	      }

	    // Trade back the results
	    std::vector<unsigned int> filled_request;
	    Parallel::send_receive(procdown, request_to_fill,
				   procup, filled_request);

	    libmesh_assert (filled_request.size() == requested_ids[procup].size());
	    
	    for (unsigned int i=0; i<filled_request.size(); i++)
	      child_id_map[requested_ids[procup][i]] = 
		std::make_pair (procup,
				filled_request[i]);
	  }
	// overwrite the parent_id_map with the child_id_map, but
	// use std::map::swap() for efficiency.
	parent_id_map.swap(child_id_map); /**/ child_id_map.clear();
      }
    }
#endif // LIBMESH_ENABLE_AMR
  if (libMesh::processor_id() == 0)
    libmesh_assert (next_global_elem == n_elem);
  
}



void XdrIO::write_serialized_nodes (Xdr &io, const unsigned int n_nodes) const
{
  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  libmesh_assert (n_nodes == mesh.n_nodes());

  std::vector<unsigned int> xfer_ids;
  std::vector<Real>         xfer_coords, &coords=xfer_coords;

  std::vector<std::vector<unsigned int> > recv_ids   (libMesh::n_processors());;
  std::vector<std::vector<Real> >         recv_coords(libMesh::n_processors());

  unsigned int n_written=0;

  for (unsigned int blk=0, last_node=0; last_node<n_nodes; blk++)
    {
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

      // We will have lots of simultaneous receives if we are
      // processor 0, so let's use nonblocking receives.
      std::vector<Parallel::request>
        id_request_handles(libMesh::n_processors()-1),
        coord_request_handles(libMesh::n_processors()-1);

      const unsigned int id_tag=0, coord_tag=1;
      
      // Post the receives -- do this on processor 0 only.
      if (libMesh::processor_id() == 0)
        {
          for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
            {
              recv_ids[pid].resize(ids_size[pid]);
              recv_coords[pid].resize(coords_size[pid]);
            
              if (pid == 0)
                {
                  recv_ids[0] = xfer_ids;
                  recv_coords[0] = xfer_coords;
                }
              else
                {
	          Parallel::nonblocking_receive (pid, recv_ids[pid],
					         id_request_handles[pid-1],
                                                 id_tag);
	          Parallel::nonblocking_receive (pid, recv_coords[pid],
					         coord_request_handles[pid-1],
                                                 coord_tag);
                }
            }
        }
      else
        {
          // Send -- do this on all other processors.
          Parallel::send(0, xfer_ids,    id_tag);
          Parallel::send(0, xfer_coords, coord_tag);
        }

      // -------------------------------------------------------
      // Receive the messages and write the output on processor 0.
      if (libMesh::processor_id() == 0)
        {
          // Wait for all the receives to complete. We have no
          // need for the statuses since we already know the
          // buffer sizes.
	  Parallel::wait (id_request_handles);
	  Parallel::wait (coord_request_handles);
          
          // Write the coordinates in this block.
	  unsigned int tot_id_size=0, tot_coord_size=0;
	  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	    {
	      tot_id_size    += recv_ids[pid].size();
	      tot_coord_size += recv_coords[pid].size();
	    }

	  libmesh_assert (tot_id_size <= std::min(io_blksize, n_nodes));
	  libmesh_assert (tot_coord_size == 3*tot_id_size);

	  coords.resize (tot_coord_size);
	  for (unsigned int pid=0; pid<libMesh::n_processors(); pid++)
	    for (unsigned int idx=0; idx<recv_ids[pid].size(); idx++)
	      {
		const unsigned int local_idx = recv_ids[pid][idx] - first_node;
		libmesh_assert ((3*local_idx+2) < coords.size());
		libmesh_assert ((3*idx+2)      < recv_coords[pid].size());
		
		coords[3*local_idx+0] = recv_coords[pid][3*idx+0];
		coords[3*local_idx+1] = recv_coords[pid][3*idx+1];
		coords[3*local_idx+2] = recv_coords[pid][3*idx+2];

		n_written++;
	      }

	  io.data_stream (coords.empty() ? NULL : &coords[0], coords.size(), 3);
	}
    }
  if (libMesh::processor_id() == 0)
    libmesh_assert (n_written == n_nodes);
}



void XdrIO::write_serialized_bcs (Xdr &io, const unsigned int n_bcs) const
{
  libmesh_assert (io.writing());
  
  if (!n_bcs) return;
  
  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  
  // and our boundary info object
  const BoundaryInfo &boundary_info = *mesh.boundary_info;

  unsigned int n_bcs_out = n_bcs;
  if (libMesh::processor_id() == 0)
    io.data (n_bcs_out, "# number of boundary conditions");
  n_bcs_out = 0;

  std::vector<int> xfer_bcs, recv_bcs;
  std::vector<unsigned int> bc_sizes(libMesh::n_processors());;
  
  // Boundary conditions are only specified for level-0 elements
  MeshBase::const_element_iterator
    it  = mesh.local_level_elements_begin(0),
    end = mesh.local_level_elements_end(0);

  unsigned int n_local_level_0_elem=0;
  for (; it!=end; ++it, n_local_level_0_elem++)
    {
      const Elem *elem = *it;

      for (unsigned int s=0; s<elem->n_sides(); s++)
	if (elem->neighbor(s) == NULL)
	  {
	    const short int bc_id = 
	      boundary_info.boundary_id (elem, s);

	    if (bc_id != BoundaryInfo::invalid_id)
	      {
		xfer_bcs.push_back (n_local_level_0_elem); 
		xfer_bcs.push_back (s) ;
		xfer_bcs.push_back (bc_id);
	      }
	  }
    }

  xfer_bcs.push_back(n_local_level_0_elem);
  unsigned int my_size = xfer_bcs.size();
  Parallel::gather (0, my_size, bc_sizes);
  
  // All processors send their xfer buffers to processor 0
  // Processor 0 will receive all buffers and write out the bcs
  if (libMesh::processor_id() == 0)
    {
      for (unsigned int pid=0, elem_offset=0; pid<libMesh::n_processors(); pid++)
	{
	  recv_bcs.resize(bc_sizes[pid]);
          if (pid == 0)
	    recv_bcs = xfer_bcs;
          else
	    Parallel::receive (pid, recv_bcs);

	  const unsigned int n_local_level_0_elem 
	    = recv_bcs.back(); recv_bcs.pop_back();
	  
	  for (unsigned int idx=0; idx<recv_bcs.size(); idx += 3, n_bcs_out++)
	    recv_bcs[idx+0] += elem_offset;

	  io.data_stream (recv_bcs.empty() ? NULL : &recv_bcs[0], recv_bcs.size(), 3);
	  elem_offset += n_local_level_0_elem;
	}    
      libmesh_assert (n_bcs == n_bcs_out);
    }
  else
    Parallel::send (0, xfer_bcs);
} 



void XdrIO::read (const std::string& name)
{
  // Only open the file on processor 0 -- this is especially important because
  // there may be an underlying bzip/bunzip going on, and multiple simultaneous
  // calls will produce a race condition.
  Xdr io (libMesh::processor_id() == 0 ? name : "", this->binary() ? DECODE : READ);
   
  // convenient reference to our mesh
  MeshBase &mesh = MeshInput<MeshBase>::mesh();

  // get the version string.   
  if (libMesh::processor_id() == 0)    
    io.data (this->version());
  Parallel::broadcast (this->version());  
  
  // note that for "legacy" files the first entry is an
  // integer -- not a string at all.
  this->legacy() = !(this->version().find("libMesh") < this->version().size());

  // Check for a legacy version format.
  if (this->legacy())
    {
      io.close();
      LegacyXdrIO(mesh, this->binary()).read(name);
      return;
    }

  START_LOG("read()","XdrIO");
  
  unsigned int n_elem, n_nodes;
  if (libMesh::processor_id() == 0)
    {
      io.data (n_elem);
      io.data (n_nodes);
      io.data (this->boundary_condition_file_name()); // std::cout << "bc_file="  << this->boundary_condition_file_name() << std::endl;
      io.data (this->subdomain_map_file_name());      // std::cout << "sid_file=" << this->subdomain_map_file_name()      << std::endl;
      io.data (this->partition_map_file_name());      // std::cout << "pid_file=" << this->partition_map_file_name()      << std::endl;
      io.data (this->polynomial_level_file_name());   // std::cout << "pl_file="  << this->polynomial_level_file_name()   << std::endl;
    }

  //TODO:[BSK] a little extra effort here could change this to two broadcasts...
  Parallel::broadcast (n_elem);
  Parallel::broadcast (n_nodes);
  Parallel::broadcast (this->boundary_condition_file_name());
  Parallel::broadcast (this->subdomain_map_file_name());
  Parallel::broadcast (this->partition_map_file_name());
  Parallel::broadcast (this->polynomial_level_file_name());
  
  // Tell the mesh how many nodes/elements to expect. Depending on the mesh type,
  // this may allow for efficient adding of nodes/elements.
  mesh.reserve_elem(n_elem);
  mesh.reserve_nodes(n_nodes);

  // read connectivity
  this->read_serialized_connectivity (io, n_elem);

  // read the nodal locations
  this->read_serialized_nodes (io, n_nodes);  

  // read the boundary conditions
  this->read_serialized_bcs (io);

  STOP_LOG("read()","XdrIO");
   
  // set the node processor ids
  Partitioner::set_node_processor_ids(mesh);
}



void XdrIO::read_serialized_connectivity (Xdr &io, const unsigned int n_elem)
{
  libmesh_assert (io.reading());

  if (!n_elem) return;
  
  const bool
    read_p_level      = ("." == this->polynomial_level_file_name()),
    read_partitioning = ("." == this->partition_map_file_name()),
    read_subdomain_id = ("." == this->subdomain_map_file_name());
  
  // convenient reference to our mesh
  MeshBase &mesh = MeshInput<MeshBase>::mesh();

  std::vector<unsigned int> conn, input_buffer(100 /* oversized ! */);

  int level=-1;
  
  for (unsigned int blk=0, first_elem=0, last_elem=0, n_elem_at_level=0, n_processed_at_level=0;
       last_elem<n_elem; blk++)
    {
      first_elem = blk*io_blksize;
      last_elem  = std::min((blk+1)*io_blksize, n_elem);

      conn.clear();

      if (libMesh::processor_id() == 0)
	for (unsigned int e=first_elem; e<last_elem; e++, n_processed_at_level++)
	  {
	    if (n_processed_at_level == n_elem_at_level)
	      {
		// get the number of elements to read at this level
		io.data (n_elem_at_level);
		n_processed_at_level = 0;
		level++;
	      }
	    
	    // get the element type, 
	    io.data_stream (&input_buffer[0], 1);

	    // maybe the parent
	    if (level)
	      io.data_stream (&input_buffer[1], 1);
	    else
	      input_buffer[1] = libMesh::invalid_uint;

	    // maybe the processor id
	    if (read_partitioning)
	      io.data_stream (&input_buffer[2], 1);
	    else
	      input_buffer[2] = 0;

	    // maybe the subdomain id
	    if (read_subdomain_id)
	      io.data_stream (&input_buffer[3], 1);
	    else 
	      input_buffer[3] = 0;

	    // maybe the p level
	    if (read_p_level)
	      io.data_stream (&input_buffer[4], 1);
	    else
	      input_buffer[4] = 0;
	    
	    // and all the nodes
	    libmesh_assert (5+Elem::type_to_n_nodes_map[input_buffer[0]] < input_buffer.size());
	    io.data_stream (&input_buffer[5], Elem::type_to_n_nodes_map[input_buffer[0]]);
	    conn.insert (conn.end(),
			 input_buffer.begin(),
			 input_buffer.begin() + 5 + Elem::type_to_n_nodes_map[input_buffer[0]]);
	  }
      
      unsigned int conn_size = conn.size();
      Parallel::broadcast(conn_size);	    
      conn.resize (conn_size);
      Parallel::broadcast (conn);

      // All processors now have the connectivity for this block.
      std::vector<unsigned int>::const_iterator it = conn.begin();      
      for (unsigned int e=first_elem; e<last_elem; e++)
	{
	  const ElemType elem_type        = static_cast<ElemType>(*it); ++it;
	  const unsigned int parent_id    = *it; ++it;
	  const unsigned int processor_id = *it; ++it;
	  const unsigned int subdomain_id = *it; ++it;
#ifdef LIBMESH_ENABLE_AMR
	  const unsigned int p_level      = *it;
#endif
	  ++it;

	  Elem *parent = (parent_id == libMesh::invalid_uint) ? NULL : mesh.elem(parent_id);

	  Elem *elem = Elem::build (elem_type, parent).release();
	  elem->set_id() = e;
	  elem->processor_id() = processor_id;
	  elem->subdomain_id() = subdomain_id;
#ifdef LIBMESH_ENABLE_AMR
	  elem->hack_p_level(p_level);

	  if (parent)
	    {
	      parent->add_child(elem);
	      parent->set_refinement_flag (Elem::INACTIVE);
	      elem->set_refinement_flag   (Elem::JUST_REFINED);
	    }
#endif

	  for (unsigned int n=0; n<elem->n_nodes(); n++, ++it)
	    {
	      const unsigned int global_node_number = *it;

	      elem->set_node(n) = 
		mesh.add_point (Point(), global_node_number);
	    }
	  
	  mesh.add_elem(elem);
	}
    }
}



void XdrIO::read_serialized_nodes (Xdr &io, const unsigned int n_nodes)
{
  libmesh_assert (io.reading());
  
  // convenient reference to our mesh
  MeshBase &mesh = MeshInput<MeshBase>::mesh();
  
  if (!mesh.n_nodes()) return;

  // At this point the elements have been read from file and placeholder nodes
  // have been assigned.  These nodes, however, do not have the proper (x,y,z)
  // locations.  This method will read all the nodes from disk, and each processor
  // can then grab the individual values it needs.

  // build up a list of the nodes contained in our local mesh.  These are the nodes
  // stored on the local processor whose (x,y,z) values need to be corrected.
  std::vector<unsigned int> needed_nodes; needed_nodes.reserve (mesh.n_nodes());
  {
    MeshBase::node_iterator
      it  = mesh.nodes_begin(),
      end = mesh.nodes_end();

    for (; it!=end; ++it)
      needed_nodes.push_back((*it)->id());

    std::sort (needed_nodes.begin(), needed_nodes.end());

    // We should not have any duplicate node->id()s
    libmesh_assert (std::unique(needed_nodes.begin(), needed_nodes.end()) == needed_nodes.end());
  }

  // Get the nodes in blocks.
  std::vector<Real> coords;
  std::pair<std::vector<unsigned int>::iterator,
            std::vector<unsigned int>::iterator> pos;
  pos.first = needed_nodes.begin();

  for (unsigned int blk=0, first_node=0, last_node=0; last_node<n_nodes; blk++)
    {
      first_node = blk*io_blksize;
      last_node  = std::min((blk+1)*io_blksize, n_nodes);
      
      coords.resize(3*(last_node - first_node));
	
      if (libMesh::processor_id() == 0)
	io.data_stream (coords.empty() ? NULL : &coords[0], coords.size());

      // For large numbers of processors the majority of processors at any given
      // block may not actually need these data.  It may be worth profiling this,
      // although it is expected that disk IO will be the bottleneck
      Parallel::broadcast (coords);
      
      for (unsigned int n=first_node, idx=0; n<last_node; n++, idx+=3)
	{
	  // first see if we need this node.  use pos.first as a smart lower
	  // bound, this will ensure that the size of the searched range
	  // decreases as we match nodes.
	  pos = std::equal_range (pos.first, needed_nodes.end(), n);

	  if (pos.first != pos.second) // we need this node.
	    {
	      libmesh_assert (*pos.first == n);
	      mesh.node(n) = 
		Point (coords[idx+0],
		       coords[idx+1],
		       coords[idx+2]);

	    }
	}	  
    }    
}



void XdrIO::read_serialized_bcs (Xdr &io)
{
  if (this->boundary_condition_file_name() == "n/a") return;

  libmesh_assert (io.reading());
  
  // convenient reference to our mesh
  MeshBase &mesh = MeshInput<MeshBase>::mesh();
  
  // and our boundary info object
  BoundaryInfo &boundary_info = *mesh.boundary_info;

  std::vector<ElemBCData> elem_bc_data;
  std::vector<int> input_buffer;

  unsigned int n_bcs=0;
  if (libMesh::processor_id() == 0) 
    io.data (n_bcs);
  Parallel::broadcast (n_bcs);

  for (unsigned int blk=0, first_bc=0, last_bc=0; last_bc<n_bcs; blk++)
    {
      first_bc = blk*io_blksize;
      last_bc  = std::min((blk+1)*io_blksize, n_bcs);

      input_buffer.resize (3*(last_bc - first_bc));

      if (libMesh::processor_id() == 0)
	io.data_stream (input_buffer.empty() ? NULL : &input_buffer[0], input_buffer.size());
      
      Parallel::broadcast (input_buffer);
      elem_bc_data.clear(); /**/ elem_bc_data.reserve (input_buffer.size()/3);

      // convert the input_buffer to ElemBCData to facilitate searching
      for (unsigned int idx=0; idx<input_buffer.size(); idx+=3)
	elem_bc_data.push_back (ElemBCData(input_buffer[idx+0],
					   input_buffer[idx+1],
					   input_buffer[idx+2]));
      input_buffer.clear();
      // note that while the files *we* write should already be sorted by
      // element id this is not necessarily guaranteed.
      std::sort (elem_bc_data.begin(), elem_bc_data.end());

      MeshBase::const_element_iterator
	it  = mesh.level_elements_begin(0),
	end = mesh.level_elements_end(0);

      // Look for BCs in this block for all the level-0 elements we have 
      // (not just local ones).  Do this by finding all the entries
      // in elem_bc_data whose elem_id match the ID of the current element.
      // We cannot rely on NULL neighbors at this point since the neighbor
      // data structure has not been initialized.
      for (std::pair<std::vector<ElemBCData>::iterator,
	             std::vector<ElemBCData>::iterator> pos; it!=end; ++it)
#if defined(__SUNPRO_CC) || defined(__PGI)
	for (pos = std::equal_range (elem_bc_data.begin(), elem_bc_data.end(), (*it)->id(), CompareIntElemBCData());
#else
	for (pos = std::equal_range (elem_bc_data.begin(), elem_bc_data.end(), (*it)->id());
#endif
	     pos.first != pos.second; ++pos.first)
	  {
	    libmesh_assert (pos.first->elem_id == (*it)->id());
	    libmesh_assert (pos.first->side < (*it)->n_sides());

	    boundary_info.add_side (*it, pos.first->side, pos.first->bc_id);
	  }	
    }
}



void XdrIO::pack_element (std::vector<unsigned int> &conn, const Elem *elem, 
			  const unsigned int parent_id, const unsigned int parent_pid) const
{
  libmesh_assert (elem != NULL);
  libmesh_assert (elem->n_nodes() == Elem::type_to_n_nodes_map[elem->type()]);
  
  conn.push_back(elem->n_nodes());

  conn.push_back (elem->type());
  conn.push_back (elem->id());
  
  if (parent_id != libMesh::invalid_uint)
    {
      conn.push_back (parent_id);
      libmesh_assert (parent_pid != libMesh::invalid_uint);
      conn.push_back (parent_pid);
    }
  
  conn.push_back (elem->processor_id());
  conn.push_back (elem->subdomain_id());
  
#ifdef LIBMESH_ENABLE_AMR
  conn.push_back (elem->p_level());
#endif
  
  for (unsigned int n=0; n<elem->n_nodes(); n++)
    conn.push_back (elem->node(n));                 
}
