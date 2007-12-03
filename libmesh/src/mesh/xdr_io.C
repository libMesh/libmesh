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
#include "boundary_info.h"
#include "parallel.h"


// ------------------------------------------------------------
// XdrIO static data
const unsigned int XdrIO::io_blksize = 2;



// ------------------------------------------------------------
// XdrIO members
XdrIO::XdrIO (MeshBase& mesh, const bool binary) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _binary (binary)
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




bool & XdrIO::binary ()
{
  return _binary;
}



bool XdrIO::binary () const
{
  return _binary;
}



void XdrIO::read (const std::string& name)
{
  Xdr io (name, this->binary() ? DECODE : READ);

  // get the version string
  std::string version;
  io.data(version);
  
  std::cout << "version=" << version << std::endl;

  // Check for a legacy version format.
  if (!(version.find("libMesh") < version.size()))
    {
      io.close();
      LegacyXdrIO(MeshInput<MeshBase>::mesh(), this->binary()).read(name);
      return;
    }


  unsigned int n_elem, n_nodes;
  io.data (n_elem);
  io.data (n_nodes);
  
  std::cout << "n_elem=" << n_elem << ", n_nodes=" << n_nodes << std::endl;
  
  std::string bc_file;
  std::string sid_file;
  std::string pid_file;
  std::string pl_file;

  io.data (bc_file);  std::cout << "bc_file="  << bc_file  << std::endl;
  io.data (sid_file); std::cout << "sid_file=" << sid_file << std::endl;
  io.data (pid_file); std::cout << "pid_file=" << pid_file << std::endl;
  io.data (pl_file);  std::cout << "pl_file="  << pl_file  << std::endl;
}



void XdrIO::write (const std::string& name)
{
  if (true)
    {
      LegacyXdrIO(MeshOutput<MeshBase>::mesh(), this->binary()).write(name);
      return;
    }
  
  Xdr io ((libMesh::processor_id() == 0) ? name : "", this->binary() ? ENCODE : WRITE);

  // set the version string
  std::string version("libMesh-0.7.0+");

  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  
  unsigned int
    n_elem  = mesh.n_elem(),
    n_nodes = mesh.n_nodes(),
    n_bcs   = mesh.boundary_info->n_boundary_conds();

  std::string bc_file(n_bcs ? "." : "n/a");
  std::string sid_file("n/a");
  std::string pid_file("n/a");
  std::string pl_file("n/a");

  
  // write the header
  if (libMesh::processor_id() == 0)
    {
      io.data(version);
      
      io.data (n_elem,  "# number of elements");
      io.data (n_nodes, "# number of nodes");
      
      io.data (bc_file, "# boundary condition specification file");
      io.data (sid_file,"# subdomain id specification file");
      io.data (pid_file,"# processor id specification file");
      io.data (pl_file, "# p-level specification file");
      
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


  unsigned int last_node=0, n_written=0;

  for (unsigned int blk=0; last_node<n_nodes; blk++)
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
