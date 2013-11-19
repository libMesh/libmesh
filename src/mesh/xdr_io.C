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
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <string>

// Local includes
#include "libmesh/xdr_io.h"
#include "libmesh/legacy_xdr_io.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/parallel.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/partitioner.h"
#include "libmesh/libmesh_logging.h"

namespace libMesh
{

//-----------------------------------------------
// anonymous namespace for implementation details
namespace {
  struct DofBCData
  {
    unsigned int       dof_id;
    unsigned short int side;
    boundary_id_type   bc_id;

    // Default constructor
    DofBCData (unsigned int       dof_id_in=0,
               unsigned short int side_in=0,
               boundary_id_type   bc_id_in=0) :
      dof_id(dof_id_in),
      side(side_in),
      bc_id(bc_id_in)
    {}

    // comparison operator
    bool operator < (const DofBCData &other) const
    {
      if (this->dof_id == other.dof_id)
	return (this->side < other.side);

      return this->dof_id < other.dof_id;
    }
  };

  // comparison operator
  bool operator < (const unsigned int &other_dof_id,
		   const DofBCData &dof_bc)
  {
    return other_dof_id < dof_bc.dof_id;
  }

  bool operator < (const DofBCData &dof_bc,
		   const unsigned int &other_dof_id)

  {
    return dof_bc.dof_id < other_dof_id;
  }


  // For some reason SunStudio does not seem to accept the above
  // comparison functions for use in
  // std::equal_range (ElemBCData::iterator, ElemBCData::iterator, unsigned int);
#if defined(__SUNPRO_CC) || defined(__PGI)
  struct CompareIntDofBCData
  {
    bool operator()(const unsigned int &other_dof_id,
		    const DofBCData &dof_bc)
    {
      return other_dof_id < dof_bc.dof_id;
    }

    bool operator()(const DofBCData &dof_bc,
		    const unsigned int &other_dof_id)
    {
      return dof_bc.dof_id < other_dof_id;
    }
  };
#endif
}



// ------------------------------------------------------------
// XdrIO static data
const std::size_t XdrIO::io_blksize = 128000;



// ------------------------------------------------------------
// XdrIO members
XdrIO::XdrIO (MeshBase& mesh, const bool binary_in) :
  MeshInput<MeshBase> (mesh,/* is_parallel_format = */ true),
  MeshOutput<MeshBase>(mesh,/* is_parallel_format = */ true),
  ParallelObject      (mesh),
  _binary             (binary_in),
  _legacy             (false),
  _write_serial       (false),
  _write_parallel     (false),
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  _write_unique_id    (true),
#else
  _write_unique_id    (false),
#endif
  _field_width        (4),   // In 0.7.0, all fields are 4 bytes, in 0.9.2 they can vary
  _version            ("libMesh-0.9.2+"),
  _bc_file_name       ("n/a"),
  _partition_map_file ("n/a"),
  _subdomain_map_file ("n/a"),
  _p_level_file       ("n/a")
{
}



XdrIO::XdrIO (const MeshBase& mesh, const bool binary_in) :
  MeshOutput<MeshBase>(mesh,/* is_parallel_format = */ true),
  ParallelObject      (mesh),
  _binary (binary_in)
{
}



XdrIO::~XdrIO ()
{
}



void XdrIO::write (const std::string& name)
{
  if (this->legacy())
    {
      libmesh_deprecated();

      // We don't support writing parallel files in the legacy format
      libmesh_assert(!this->_write_parallel);

      LegacyXdrIO(MeshOutput<MeshBase>::mesh(), this->binary()).write(name);
      return;
    }

  Xdr io ((this->processor_id() == 0) ? name : "", this->binary() ? ENCODE : WRITE);

  START_LOG("write()","XdrIO");

  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();

  header_id_type
    n_elem     = mesh.n_elem(),
    n_nodes    = mesh.n_nodes();

  libmesh_assert(n_elem == mesh.n_elem());
  libmesh_assert(n_nodes == mesh.n_nodes());

  std::size_t
    n_bcs      = mesh.boundary_info->n_boundary_conds();
  std::size_t
    n_nodesets = mesh.boundary_info->n_nodeset_conds();
  unsigned int
    n_p_levels = MeshTools::n_p_levels (mesh);

  bool write_parallel_files = this->write_parallel();

  //-------------------------------------------------------------
  // For all the optional files -- the default file name is "n/a".
  // However, the user may specify an optional external file.

  // If there are BCs and the user has not already provided a
  // file name then write to "."
  if ((n_bcs || n_nodesets) &&
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
  if (this->processor_id() == 0)
    {
      std::string full_ver = this->version() + (write_parallel_files ?  " parallel" : "");
      io.data (full_ver);

      io.data (n_elem,  "# number of elements");
      io.data (n_nodes, "# number of nodes");

      io.data (this->boundary_condition_file_name(), "# boundary condition specification file");
      io.data (this->subdomain_map_file_name(),      "# subdomain id specification file");
      io.data (this->partition_map_file_name(),      "# processor id specification file");
      io.data (this->polynomial_level_file_name(),   "# p-level specification file");

      // Version 0.9.2+ introduces sizes for each type
      header_id_type write_size = sizeof(xdr_id_type), zero_size = 0;

      const bool
        write_p_level      = ("." == this->polynomial_level_file_name()),
        write_partitioning = ("." == this->partition_map_file_name()),
        write_subdomain_id = ("." == this->subdomain_map_file_name()),
        write_bcs          = ("." == this->boundary_condition_file_name());

      io.data (write_size, "# type size");
      io.data (_write_unique_id   ? write_size : zero_size, "# uid size");
      io.data (write_partitioning ? write_size : zero_size, "# pid size");
      io.data (write_subdomain_id ? write_size : zero_size, "# sid size");
      io.data (write_p_level      ? write_size : zero_size, "# p-level size");
      // Boundary Condition sizes
      io.data (write_bcs          ? write_size : zero_size, "# eid size");   // elem id
      io.data (write_bcs          ? write_size : zero_size, "# side size "); // side number
      io.data (write_bcs          ? write_size : zero_size, "# bid size");   // boundary id
    }

  if (write_parallel_files)
    {
      // Parallel xdr mesh files aren't implemented yet; until they
      // are we'll just warn the user and write a serial file.
      libMesh::out << "Warning!  Parallel xda/xdr is not yet implemented.\n";
      libMesh::out << "Writing a serialized file instead." << std::endl;

      // write subdomain names
      this->write_serialized_subdomain_names(io);

      // write connectivity
      this->write_serialized_connectivity (io, n_elem);

      // write the nodal locations
      this->write_serialized_nodes (io, n_nodes);

      // write the boundary condition information
      this->write_serialized_bcs (io, n_bcs);

      // write the nodeset information
      this->write_serialized_nodesets (io, n_nodesets);
    }
  else
    {
      // write subdomain names
      this->write_serialized_subdomain_names(io);

      // write connectivity
      this->write_serialized_connectivity (io, n_elem);

      // write the nodal locations
      this->write_serialized_nodes (io, n_nodes);

      // write the boundary condition information
      this->write_serialized_bcs (io, n_bcs);

      // write the nodeset information
      this->write_serialized_nodesets (io, n_nodesets);
    }

  if(mesh.boundary_info->n_edge_conds() > 0)
  {
    libMesh::out << "Warning: Mesh contains edge boundary IDs, but these "
                 << "are not supported by the XDR format."
                 << std::endl;
  }

  STOP_LOG("write()","XdrIO");

  // pause all processes until the writing ends -- this will
  // protect for the pathological case where a write is
  // followed immediately by a read.  The write must be
  // guaranteed to complete first.
  io.close();
  this->comm().barrier();
}



void XdrIO::write_serialized_subdomain_names(Xdr &io) const
{
  if (this->processor_id() == 0)
  {
    const MeshBase &mesh = MeshOutput<MeshBase>::mesh();

    const std::map<subdomain_id_type, std::string> & subdomain_map = mesh.get_subdomain_name_map();

    std::vector<header_id_type> subdomain_ids;   subdomain_ids.reserve(subdomain_map.size());
    std::vector<std::string>  subdomain_names; subdomain_names.reserve(subdomain_map.size());

    // We need to loop over the map and make sure that there aren't any invalid entries.  Since we
    // return writable references in mesh_base, it's possible for the user to leave some entity names
    // blank.  We can't write those to the XDA file.
    header_id_type n_subdomain_names = 0;
    std::map<subdomain_id_type, std::string>::const_iterator it_end = subdomain_map.end();
    for (std::map<subdomain_id_type, std::string>::const_iterator it = subdomain_map.begin(); it != it_end; ++it)
    {
      if (!it->second.empty())
      {
        n_subdomain_names++;
        subdomain_ids.push_back(it->first);
        subdomain_names.push_back(it->second);
      }
    }

    io.data(n_subdomain_names, "# subdomain id to name map");
    // Write out the ids and names in two vectors
    if (n_subdomain_names)
    {
      io.data(subdomain_ids);
      io.data(subdomain_names);
    }
  }
}



void XdrIO::write_serialized_connectivity (Xdr &io, const dof_id_type libmesh_dbg_var(n_elem)) const
{
  libmesh_assert (io.writing());

  const bool
    write_p_level      = ("." == this->polynomial_level_file_name()),
    write_partitioning = ("." == this->partition_map_file_name()),
    write_subdomain_id = ("." == this->subdomain_map_file_name());

  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  libmesh_assert_equal_to (n_elem, mesh.n_elem());

  // We will only write active elements and their parents.
  const unsigned int n_active_levels = MeshTools::n_active_levels (mesh);
  std::vector<xdr_id_type>
    n_global_elem_at_level(n_active_levels), n_local_elem_at_level(n_active_levels);

  MeshBase::const_element_iterator it  = mesh.local_elements_end(), end=it;

  // Find the number of local and global elements at each level
#ifndef NDEBUG
  dof_id_type tot_n_elem = 0;
#endif
  for (unsigned int level=0; level<n_active_levels; level++)
    {
      it  = mesh.local_level_elements_begin(level);
      end = mesh.local_level_elements_end(level);

      n_local_elem_at_level[level] = n_global_elem_at_level[level] = MeshTools::n_elem(it, end);

      this->comm().sum(n_global_elem_at_level[level]);
#ifndef NDEBUG
      tot_n_elem += n_global_elem_at_level[level];
#endif
      libmesh_assert_less_equal (n_global_elem_at_level[level], n_elem);
      libmesh_assert_less_equal (tot_n_elem, n_elem);
    }

  std::vector<xdr_id_type>
    xfer_conn, recv_conn;
  std::vector<dof_id_type>
    n_elem_on_proc(this->n_processors()), processor_offsets(this->n_processors());
  std::vector<xdr_id_type> output_buffer;
  std::vector<std::size_t>
    xfer_buf_sizes(this->n_processors());

  typedef std::map<dof_id_type, std::pair<processor_id_type, dof_id_type> > id_map_type;
  id_map_type parent_id_map, child_id_map;

  dof_id_type my_next_elem=0, next_global_elem=0;

  //-------------------------------------------
  // First write the level-0 elements directly.
  it  = mesh.local_level_elements_begin(0);
  end = mesh.local_level_elements_end(0);
  for (; it != end; ++it)
    {
      pack_element (xfer_conn, *it);
      parent_id_map[(*it)->id()] = std::make_pair(this->processor_id(),
						  my_next_elem++);
    }
  xfer_conn.push_back(my_next_elem); // toss in the number of elements transferred.

  std::size_t my_size = xfer_conn.size();
  this->comm().gather (0, my_next_elem, n_elem_on_proc);
  this->comm().gather (0, my_size,      xfer_buf_sizes);

  processor_offsets[0] = 0;
  for (unsigned int pid=1; pid<this->n_processors(); pid++)
    processor_offsets[pid] = processor_offsets[pid-1] + n_elem_on_proc[pid-1];

  // All processors send their xfer buffers to processor 0.
  // Processor 0 will receive the data and write out the elements.
  if (this->processor_id() == 0)
    {
      // Write the number of elements at this level.
      {
	std::string comment = "# n_elem at level 0", legend  = ", [ type ";
        if (_write_unique_id)
          legend += "uid ";
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

      for (unsigned int pid=0; pid<this->n_processors(); pid++)
	{
	  recv_conn.resize(xfer_buf_sizes[pid]);
          if (pid == 0)
	    recv_conn = xfer_conn;
          else
	    this->comm().receive (pid, recv_conn);

	  // at a minimum, the buffer should contain the number of elements,
	  // which could be 0.
	  libmesh_assert (!recv_conn.empty());

          {
	  const xdr_id_type n_elem_received = recv_conn.back();
	  std::vector<xdr_id_type>::const_iterator recv_conn_iter = recv_conn.begin();

	  for (xdr_id_type elem=0; elem<n_elem_received; elem++, next_global_elem++)
	    {
	      output_buffer.clear();
	      const xdr_id_type n_nodes = *recv_conn_iter; ++recv_conn_iter;
	      output_buffer.push_back(*recv_conn_iter);     /* type       */ ++recv_conn_iter;

              if (_write_unique_id)
                output_buffer.push_back(*recv_conn_iter);   /* unique_id  */ ++recv_conn_iter;

	      if (write_partitioning)
		output_buffer.push_back(*recv_conn_iter); /* processor id */ ++recv_conn_iter;

	      if (write_subdomain_id)
		output_buffer.push_back(*recv_conn_iter); /* subdomain id */ ++recv_conn_iter;

#ifdef LIBMESH_ENABLE_AMR
	      if (write_p_level)
		output_buffer.push_back(*recv_conn_iter); /* p level      */ ++recv_conn_iter;
#endif
	      for (dof_id_type node=0; node<n_nodes; node++, ++recv_conn_iter)
		output_buffer.push_back(*recv_conn_iter);

	      io.data_stream (&output_buffer[0], output_buffer.size(), output_buffer.size());
	    }
          }
	}
    }
  else
    this->comm().send (0, xfer_conn);

#ifdef LIBMESH_ENABLE_AMR
  //--------------------------------------------------------------------
  // Next write the remaining elements indirectly through their parents.
  // This will insure that the children are written in the proper order
  // so they can be reconstructed properly.
  for (unsigned int level=1; level<n_active_levels; level++)
    {
      xfer_conn.clear();

      it  = mesh.local_level_elements_begin(level-1);
      end = mesh.local_level_elements_end  (level-1);

      dof_id_type my_n_elem_written_at_level = 0;
      for (; it != end; ++it)
	if (!(*it)->active()) // we only want the parents elements at this level, and
	  {                   // there is no direct iterator for this obscure use
	    const Elem *parent = *it;
	    id_map_type::iterator pos = parent_id_map.find(parent->id());
	    libmesh_assert (pos != parent_id_map.end());
	    const processor_id_type parent_pid = pos->second.first;
	    const dof_id_type parent_id  = pos->second.second;
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
      this->comm().gather (0, my_size,   xfer_buf_sizes);

      // Processor 0 will receive the data and write the elements.
      if (this->processor_id() == 0)
	{
	  // Write the number of elements at this level.
	  {
	    char buf[80];
	    std::sprintf(buf, "# n_elem at level %d", level);
	    std::string comment(buf), legend  = ", [ type ";

            if (_write_unique_id)
              legend += "uid ";
            legend += "parent ";
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

	  for (unsigned int pid=0; pid<this->n_processors(); pid++)
	    {
	      recv_conn.resize(xfer_buf_sizes[pid]);
              if (pid == 0)
	        recv_conn = xfer_conn;
              else
	        this->comm().receive (pid, recv_conn);

	      // at a minimum, the buffer should contain the number of elements,
	      // which could be 0.
	      libmesh_assert (!recv_conn.empty());

              {
	      const xdr_id_type n_elem_received = recv_conn.back();
	      std::vector<xdr_id_type>::const_iterator recv_conn_iter = recv_conn.begin();

	      for (xdr_id_type elem=0; elem<n_elem_received; elem++, next_global_elem++)
		{
		  output_buffer.clear();
		  const xdr_id_type n_nodes = *recv_conn_iter; ++recv_conn_iter;
		  output_buffer.push_back(*recv_conn_iter);                   /* type          */ ++recv_conn_iter;
                  if (_write_unique_id)
                    output_buffer.push_back(*recv_conn_iter);                 /* unique_id     */ ++recv_conn_iter;

                  const xdr_id_type parent_local_id = *recv_conn_iter; ++recv_conn_iter;
		  const xdr_id_type parent_pid      = *recv_conn_iter; ++recv_conn_iter;

		  output_buffer.push_back (parent_local_id+processor_offsets[parent_pid]);

		  if (write_partitioning)
		    output_buffer.push_back(*recv_conn_iter); /* processor id */ ++recv_conn_iter;

		  if (write_subdomain_id)
		    output_buffer.push_back(*recv_conn_iter); /* subdomain id */ ++recv_conn_iter;

		  if (write_p_level)
		    output_buffer.push_back(*recv_conn_iter); /* p level       */ ++recv_conn_iter;

		  for (xdr_id_type node=0; node<n_nodes; node++, ++recv_conn_iter)
		    output_buffer.push_back(*recv_conn_iter);

		  io.data_stream (&output_buffer[0], output_buffer.size(), output_buffer.size());
		}
              }
	    }
	}
      else
        this->comm().send  (0, xfer_conn);

      // update the processor_offsets
      processor_offsets[0] = processor_offsets.back() + n_elem_on_proc.back();
      this->comm().gather (0, my_n_elem_written_at_level, n_elem_on_proc);
      for (unsigned int pid=1; pid<this->n_processors(); pid++)
	processor_offsets[pid] = processor_offsets[pid-1] + n_elem_on_proc[pid-1];

      // Now, at the next level we will again iterate over local parents.  However,
      // those parents may have been written by other processors (at this step),
      // so we need to gather them into our *_id_maps.
      {
	std::vector<std::vector<dof_id_type> > requested_ids(this->n_processors());
	std::vector<dof_id_type> request_to_fill;

	it  = mesh.local_level_elements_begin(level);
	end = mesh.local_level_elements_end(level);

	for (; it!=end; ++it)
	  if (!child_id_map.count((*it)->id()))
	    {
	      libmesh_assert_not_equal_to ((*it)->parent()->processor_id(), this->processor_id());
	      requested_ids[(*it)->parent()->processor_id()].push_back((*it)->id());
	    }

	// Next set the child_ids
	for (unsigned int p=1; p != this->n_processors(); ++p)
	  {
	    // Trade my requests with processor procup and procdown
	    unsigned int procup = (this->processor_id() + p) %
                                   this->n_processors();
	    unsigned int procdown = (this->n_processors() +
                                     this->processor_id() - p) %
	                             this->n_processors();

	    this->comm().send_receive(procup, requested_ids[procup],
					      procdown, request_to_fill);

	    // Fill those requests by overwriting the requested ids
	    for (std::size_t i=0; i<request_to_fill.size(); i++)
	      {
		libmesh_assert (child_id_map.count(request_to_fill[i]));
		libmesh_assert_equal_to (child_id_map[request_to_fill[i]].first, procdown);

		request_to_fill[i] = child_id_map[request_to_fill[i]].second;
	      }

	    // Trade back the results
	    std::vector<dof_id_type> filled_request;
	    this->comm().send_receive(procdown, request_to_fill,
					      procup, filled_request);

	    libmesh_assert_equal_to (filled_request.size(), requested_ids[procup].size());

	    for (std::size_t i=0; i<filled_request.size(); i++)
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
  if (this->processor_id() == 0)
    libmesh_assert_equal_to (next_global_elem, n_elem);

}



void XdrIO::write_serialized_nodes (Xdr &io, const dof_id_type n_nodes) const
{
  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  libmesh_assert_equal_to (n_nodes, mesh.n_nodes());

  std::vector<dof_id_type> xfer_ids;
  std::vector<Real>         xfer_coords, &coords=xfer_coords;

  std::vector<std::vector<dof_id_type> > recv_ids   (this->n_processors());;
  std::vector<std::vector<Real> >         recv_coords(this->n_processors());

  std::size_t n_written=0;

  for (std::size_t blk=0, last_node=0; last_node<n_nodes; blk++)
    {
      const std::size_t first_node = blk*io_blksize;
      last_node = std::min((blk+1)*io_blksize, std::size_t(n_nodes));

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
#if LIBMESH_DIM > 1
	    xfer_coords.push_back(p(1));
#endif
#if LIBMESH_DIM > 2
	    xfer_coords.push_back(p(2));
#endif
	  }

      //-------------------------------------
      // Send the xfer buffers to processor 0
      std::vector<std::size_t> ids_size, coords_size;

      const std::size_t my_ids_size = xfer_ids.size();

      // explicitly gather ids_size
      this->comm().gather (0, my_ids_size, ids_size);

      // infer coords_size on processor 0
      if (this->processor_id() == 0)
	{
	  coords_size.reserve(this->n_processors());
	  for (std::size_t p=0; p<ids_size.size(); p++)
	    coords_size.push_back(LIBMESH_DIM*ids_size[p]);
	}

      // We will have lots of simultaneous receives if we are
      // processor 0, so let's use nonblocking receives.
      std::vector<Parallel::Request>
        id_request_handles(this->n_processors()-1),
        coord_request_handles(this->n_processors()-1);

      Parallel::MessageTag
        id_tag    = mesh.comm().get_unique_tag(1234),
        coord_tag = mesh.comm().get_unique_tag(1235);

      // Post the receives -- do this on processor 0 only.
      if (this->processor_id() == 0)
        {
          for (unsigned int pid=0; pid<this->n_processors(); pid++)
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
	          this->comm().receive (pid, recv_ids[pid],
						id_request_handles[pid-1],
						id_tag);
	          this->comm().receive (pid, recv_coords[pid],
						coord_request_handles[pid-1],
						coord_tag);
                }
            }
        }
      else
        {
          // Send -- do this on all other processors.
          this->comm().send(0, xfer_ids,    id_tag);
          this->comm().send(0, xfer_coords, coord_tag);
        }

      // -------------------------------------------------------
      // Receive the messages and write the output on processor 0.
      if (this->processor_id() == 0)
        {
          // Wait for all the receives to complete. We have no
          // need for the statuses since we already know the
          // buffer sizes.
	  Parallel::wait (id_request_handles);
	  Parallel::wait (coord_request_handles);

          // Write the coordinates in this block.
	  std::size_t tot_id_size=0, tot_coord_size=0;
	  for (unsigned int pid=0; pid<this->n_processors(); pid++)
	    {
	      tot_id_size    += recv_ids[pid].size();
	      tot_coord_size += recv_coords[pid].size();
	    }

	  libmesh_assert_less_equal
	    (tot_id_size, std::min(io_blksize, std::size_t(n_nodes)));
	  libmesh_assert_equal_to (tot_coord_size, LIBMESH_DIM*tot_id_size);

	  coords.resize (3*tot_id_size);
	  for (unsigned int pid=0; pid<this->n_processors(); pid++)
	    for (std::size_t idx=0; idx<recv_ids[pid].size(); idx++)
	      {
		const std::size_t local_idx = recv_ids[pid][idx] - first_node;
		libmesh_assert_less ((3*local_idx+2), coords.size());
		libmesh_assert_less ((LIBMESH_DIM*idx+LIBMESH_DIM-1), recv_coords[pid].size());

		coords[3*local_idx+0] = recv_coords[pid][LIBMESH_DIM*idx+0];
#if LIBMESH_DIM > 1
		coords[3*local_idx+1] = recv_coords[pid][LIBMESH_DIM*idx+1];
#else
		coords[3*local_idx+1] = 0.;
#endif
#if LIBMESH_DIM > 2
		coords[3*local_idx+2] = recv_coords[pid][LIBMESH_DIM*idx+2];
#else
		coords[3*local_idx+2] = 0.;
#endif

		n_written++;
	      }

	  io.data_stream (coords.empty() ? NULL : &coords[0], coords.size(), 3);
	}
    }
  if (this->processor_id() == 0)
    libmesh_assert_equal_to (n_written, n_nodes);
}



void XdrIO::write_serialized_bcs (Xdr &io, const std::size_t n_bcs) const
{
  libmesh_assert (io.writing());

  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();

  // and our boundary info object
  const BoundaryInfo &boundary_info = *mesh.boundary_info;

  // Version 0.9.2+ introduces entity names
  write_serialized_bc_names(io, boundary_info, true);  // sideset names

  header_id_type n_bcs_out = n_bcs;
  if (this->processor_id() == 0)
    io.data (n_bcs_out, "# number of boundary conditions");
  n_bcs_out = 0;

  if (!n_bcs) return;

  std::vector<xdr_id_type> xfer_bcs, recv_bcs;
  std::vector<std::size_t> bc_sizes(this->n_processors());

  // Boundary conditions are only specified for level-0 elements
  MeshBase::const_element_iterator
    it  = mesh.local_level_elements_begin(0),
    end = mesh.local_level_elements_end(0);

  dof_id_type n_local_level_0_elem=0;
  for (; it!=end; ++it, n_local_level_0_elem++)
    {
      const Elem *elem = *it;

      for (unsigned int s=0; s<elem->n_sides(); s++)
// We're supporting boundary ids on internal sides now
//	if (elem->neighbor(s) == NULL)
	  {
	    const std::vector<boundary_id_type>& bc_ids =
	      boundary_info.boundary_ids (elem, s);
	    for (std::vector<boundary_id_type>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
	      {
		const boundary_id_type bc_id = *id_it;
		if (bc_id != BoundaryInfo::invalid_id)
		  {
		    xfer_bcs.push_back (n_local_level_0_elem);
		    xfer_bcs.push_back (s) ;
		    xfer_bcs.push_back (bc_id);
		  }
	      }
	  }
    }

  xfer_bcs.push_back(n_local_level_0_elem);
  std::size_t my_size = xfer_bcs.size();
  this->comm().gather (0, my_size, bc_sizes);

  // All processors send their xfer buffers to processor 0
  // Processor 0 will receive all buffers and write out the bcs
  if (this->processor_id() == 0)
    {
      dof_id_type elem_offset = 0;
      for (unsigned int pid=0; pid<this->n_processors(); pid++)
	{
	  recv_bcs.resize(bc_sizes[pid]);
          if (pid == 0)
	    recv_bcs = xfer_bcs;
          else
	    this->comm().receive (pid, recv_bcs);

	  const dof_id_type my_n_local_level_0_elem
	    = recv_bcs.back(); recv_bcs.pop_back();

	  for (std::size_t idx=0; idx<recv_bcs.size(); idx += 3, n_bcs_out++)
	    recv_bcs[idx+0] += elem_offset;

	  io.data_stream (recv_bcs.empty() ? NULL : &recv_bcs[0], recv_bcs.size(), 3);
	  elem_offset += my_n_local_level_0_elem;
	}
      libmesh_assert_equal_to (n_bcs, n_bcs_out);
    }
  else
    this->comm().send (0, xfer_bcs);
}



void XdrIO::write_serialized_nodesets (Xdr &io, const std::size_t n_nodesets) const
{
  libmesh_assert (io.writing());

  // convenient reference to our mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();

  // and our boundary info object
  const BoundaryInfo &boundary_info = *mesh.boundary_info;

  // Version 0.9.2+ introduces entity names
  write_serialized_bc_names(io, boundary_info, false);  // nodeset names

  header_id_type n_nodesets_out = n_nodesets;
  if (this->processor_id() == 0)
    io.data (n_nodesets_out, "# number of nodesets");
  n_nodesets_out = 0;

  if (!n_nodesets) return;

  std::vector<xdr_id_type> xfer_bcs, recv_bcs;
  std::vector<std::size_t> bc_sizes(this->n_processors());

  MeshBase::const_node_iterator
    it  = mesh.local_nodes_begin(),
    end = mesh.local_nodes_end();

  dof_id_type n_node=0;
  for (; it!=end; ++it)
    {
      const Node *node = *it;
      const std::vector<boundary_id_type>& nodeset_ids =
        boundary_info.boundary_ids (node);
      for (std::vector<boundary_id_type>::const_iterator id_it=nodeset_ids.begin(); id_it!=nodeset_ids.end(); ++id_it)
        {
          const boundary_id_type bc_id = *id_it;
          if (bc_id != BoundaryInfo::invalid_id)
            {
              xfer_bcs.push_back ((*it)->id());
              xfer_bcs.push_back (bc_id);
            }
        }
    }

  xfer_bcs.push_back(n_node);
  std::size_t my_size = xfer_bcs.size();
  this->comm().gather (0, my_size, bc_sizes);

  // All processors send their xfer buffers to processor 0
  // Processor 0 will receive all buffers and write out the bcs
  if (this->processor_id() == 0)
    {
      dof_id_type node_offset = 0;
      for (unsigned int pid=0; pid<this->n_processors(); pid++)
	{
	  recv_bcs.resize(bc_sizes[pid]);
          if (pid == 0)
	    recv_bcs = xfer_bcs;
          else
	    this->comm().receive (pid, recv_bcs);

	  const dof_id_type my_n_node
	    = recv_bcs.back(); recv_bcs.pop_back();

	  for (std::size_t idx=0; idx<recv_bcs.size(); idx += 2, n_nodesets_out++)
	    recv_bcs[idx+0] += node_offset;

	  io.data_stream (recv_bcs.empty() ? NULL : &recv_bcs[0], recv_bcs.size(), 2);
	  node_offset += my_n_node;
	}
      libmesh_assert_equal_to (n_nodesets, n_nodesets_out);
    }
  else
    this->comm().send (0, xfer_bcs);
}



void XdrIO::write_serialized_bc_names (Xdr &io, const BoundaryInfo & info, bool is_sideset) const
{
  if (this->processor_id() == 0)
  {
    const std::map<boundary_id_type, std::string> & boundary_map = is_sideset ?
      info.get_sideset_name_map() : info.get_nodeset_name_map();

    std::vector<header_id_type> boundary_ids;   boundary_ids.reserve(boundary_map.size());
    std::vector<std::string>  boundary_names; boundary_names.reserve(boundary_map.size());

    // We need to loop over the map and make sure that there aren't any invalid entries.  Since we
    // return writable references in boundary_info, it's possible for the user to leave some entity names
    // blank.  We can't write those to the XDA file.
    header_id_type n_boundary_names = 0;
    std::map<boundary_id_type, std::string>::const_iterator it_end = boundary_map.end();
    for (std::map<boundary_id_type, std::string>::const_iterator it = boundary_map.begin(); it != it_end; ++it)
    {
      if (!it->second.empty())
      {
        n_boundary_names++;
        boundary_ids.push_back(it->first);
        boundary_names.push_back(it->second);
      }
    }

    if (is_sideset)
      io.data(n_boundary_names, "# sideset id to name map");
    else
      io.data(n_boundary_names, "# nodeset id to name map");
    // Write out the ids and names in two vectors
    if (n_boundary_names)
    {
      io.data(boundary_ids);
      io.data(boundary_names);
    }
  }
}



void XdrIO::read (const std::string& name)
{
  // Only open the file on processor 0 -- this is especially important because
  // there may be an underlying bzip/bunzip going on, and multiple simultaneous
  // calls will produce a race condition.
  Xdr io (this->processor_id() == 0 ? name : "", this->binary() ? DECODE : READ);

  // convenient reference to our mesh
  MeshBase &mesh = MeshInput<MeshBase>::mesh();

  // get the version string.
  if (this->processor_id() == 0)
    io.data (this->version());
  this->comm().broadcast (this->version());

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

  std::vector<header_id_type> meta_data(10, sizeof(xdr_id_type));
  // type_size, uid_size, pid_size, sid_size, p_level_size, eid_size, side_size, bid_size;
  header_id_type pos=0;

  const bool is_version_0_9_2 = this->version().find("0.9.2") != std::string::npos ? true : false;

  if (this->processor_id() == 0)
    {
      io.data (meta_data[pos++]);
      io.data (meta_data[pos++]);
      io.data (this->boundary_condition_file_name()); // libMesh::out << "bc_file="  << this->boundary_condition_file_name() << std::endl;
      io.data (this->subdomain_map_file_name());      // libMesh::out << "sid_file=" << this->subdomain_map_file_name()      << std::endl;
      io.data (this->partition_map_file_name());      // libMesh::out << "pid_file=" << this->partition_map_file_name()      << std::endl;
      io.data (this->polynomial_level_file_name());   // libMesh::out << "pl_file="  << this->polynomial_level_file_name()   << std::endl;

      if (is_version_0_9_2)
      {
        io.data (meta_data[pos++], "# type size");
        io.data (meta_data[pos++], "# uid size");
        io.data (meta_data[pos++], "# pid size");
        io.data (meta_data[pos++], "# sid size");
        io.data (meta_data[pos++], "# p-level size");
        // Boundary Condition sizes
        io.data (meta_data[pos++], "# eid size");   // elem id
        io.data (meta_data[pos++], "# side size "); // side number
        io.data (meta_data[pos++], "# bid size");   // boundary id
      }
    }

  // broadcast the n_elems, n_nodes, and size information
  this->comm().broadcast (meta_data);

  this->comm().broadcast (this->boundary_condition_file_name());
  this->comm().broadcast (this->subdomain_map_file_name());
  this->comm().broadcast (this->partition_map_file_name());
  this->comm().broadcast (this->polynomial_level_file_name());

  // Tell the mesh how many nodes/elements to expect. Depending on the mesh type,
  // this may allow for efficient adding of nodes/elements.
  header_id_type n_elem = meta_data[0];
  header_id_type n_nodes = meta_data[1];

  mesh.reserve_elem(n_elem);
  mesh.reserve_nodes(n_nodes);

  /**
   * We are future proofing the layout of this file by adding in size information for all stored types.
   * TODO: All types are stored as the same size. Use the size information to pack things efficiently.
   * For now we will assume that "type size" is how the entire file will be encoded.
   */
  if (is_version_0_9_2)
    _field_width = meta_data[2];

  if (_field_width == 4)
    {
      uint32_t type_size = 0;

      // read subdomain names
      this->read_serialized_subdomain_names(io);

      // read connectivity
      this->read_serialized_connectivity (io, n_elem, meta_data, type_size);

      // read the nodal locations
      this->read_serialized_nodes (io, n_nodes);

      // read the boundary conditions
      this->read_serialized_bcs (io, type_size);

      if (is_version_0_9_2)
        // read the nodesets
        this->read_serialized_nodesets (io, type_size);
    }
  else if (_field_width == 8)
    {
      uint64_t type_size = 0;

      // read subdomain names
      this->read_serialized_subdomain_names(io);

      // read connectivity
      this->read_serialized_connectivity (io, n_elem, meta_data, type_size);

      // read the nodal locations
      this->read_serialized_nodes (io, n_nodes);

      // read the boundary conditions
      this->read_serialized_bcs (io, type_size);

      if (is_version_0_9_2)
        // read the nodesets
        this->read_serialized_nodesets (io, type_size);
    }


  STOP_LOG("read()","XdrIO");

  // set the node processor ids
  Partitioner::set_node_processor_ids(mesh);
}



void XdrIO::read_serialized_subdomain_names(Xdr &io)
{
  const bool read_entity_info = this->version().find("0.9.2") != std::string::npos ? true : false;
  if (read_entity_info)
  {
    MeshBase &mesh = MeshInput<MeshBase>::mesh();

    unsigned int n_subdomain_names = 0;
    std::vector<header_id_type> subdomain_ids;
    std::vector<std::string>  subdomain_names;

    // Read the sideset names
    if (this->processor_id() == 0)
    {
      io.data(n_subdomain_names);

      subdomain_ids.resize(n_subdomain_names);
      subdomain_names.resize(n_subdomain_names);

      if (n_subdomain_names)
      {
        io.data(subdomain_ids);
        io.data(subdomain_names);
      }
    }

    // Broadcast the subdomain names to all processors
    this->comm().broadcast(n_subdomain_names);
    if (n_subdomain_names == 0)
      return;

    subdomain_ids.resize(n_subdomain_names);
    subdomain_names.resize(n_subdomain_names);
    this->comm().broadcast(subdomain_ids);
    this->comm().broadcast(subdomain_names);

    // Reassemble the named subdomain information
    std::map<subdomain_id_type, std::string> & subdomain_map = mesh.set_subdomain_name_map();

    for (unsigned int i=0; i<n_subdomain_names; ++i)
      subdomain_map.insert(std::make_pair(subdomain_ids[i], subdomain_names[i]));
  }
}


template <typename T>
void XdrIO::read_serialized_connectivity (Xdr &io, const dof_id_type n_elem, std::vector<header_id_type> & sizes, T)
{
  libmesh_assert (io.reading());

  if (!n_elem) return;

  const bool
    read_p_level      = ("." == this->polynomial_level_file_name()),
    read_partitioning = ("." == this->partition_map_file_name()),
    read_subdomain_id = ("." == this->subdomain_map_file_name());

  // convenient reference to our mesh
  MeshBase &mesh = MeshInput<MeshBase>::mesh();

  // Keep track of what kinds of elements this file contains
  elems_of_dimension.clear();
  elems_of_dimension.resize(4, false);

  std::vector<T> conn, input_buffer(100 /* oversized ! */);

  int level=-1;

  // Version 0.9.2+ introduces unique ids
  const size_t unique_id_size_index = 3;
  const bool read_unique_id = this->version().find("0.9.2") != std::string::npos &&
    sizes[unique_id_size_index] ? true : false;

  T n_elem_at_level=0, n_processed_at_level=0;
  for (std::size_t blk=0, first_elem=0, last_elem=0;
       last_elem<n_elem; blk++)
    {
      first_elem = blk*io_blksize;
      last_elem  = std::min((blk+1)*io_blksize, std::size_t(n_elem));

      conn.clear();

      if (this->processor_id() == 0)
	for (dof_id_type e=first_elem; e<last_elem; e++, n_processed_at_level++)
	  {
	    if (n_processed_at_level == n_elem_at_level)
	      {
		// get the number of elements to read at this level
		io.data (n_elem_at_level);
		n_processed_at_level = 0;
		level++;
	      }

            unsigned int pos = 0;
	    // get the element type,
	    io.data_stream (&input_buffer[pos++], 1);

            if (read_unique_id)
              io.data_stream (&input_buffer[pos++], 1);
            // Older versions won't have this field at all (no increment on pos)

	    // maybe the parent
	    if (level)
	      io.data_stream (&input_buffer[pos++], 1);
	    else
	      // We can't always fit DofObject::invalid_id in an
	      // xdr_id_type
	      input_buffer[pos++] = static_cast<T>(-1);

	    // maybe the processor id
	    if (read_partitioning)
	      io.data_stream (&input_buffer[pos++], 1);
	    else
	      input_buffer[pos++] = 0;

	    // maybe the subdomain id
	    if (read_subdomain_id)
	      io.data_stream (&input_buffer[pos++], 1);
	    else
	      input_buffer[pos++] = 0;

	    // maybe the p level
	    if (read_p_level)
	      io.data_stream (&input_buffer[pos++], 1);
	    else
	      input_buffer[pos++] = 0;

	    // and all the nodes
	    libmesh_assert_less (pos+Elem::type_to_n_nodes_map[input_buffer[0]], input_buffer.size());
	    io.data_stream (&input_buffer[pos], Elem::type_to_n_nodes_map[input_buffer[0]]);
	    conn.insert (conn.end(),
			 input_buffer.begin(),
			 input_buffer.begin() + pos + Elem::type_to_n_nodes_map[input_buffer[0]]);
	  }

      std::size_t conn_size = conn.size();
      this->comm().broadcast(conn_size);
      conn.resize (conn_size);
      this->comm().broadcast (conn);

      // All processors now have the connectivity for this block.
      typename std::vector<T>::const_iterator it = conn.begin();
      for (dof_id_type e=first_elem; e<last_elem; e++)
	{
	  const ElemType elem_type        = static_cast<ElemType>(*it); ++it;
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          unique_id_type unique_id = 0;
#endif
          if (read_unique_id)
          {
#ifdef LIBMESH_ENABLE_UNIQUE_ID
            unique_id  = libmesh_cast_int<unique_id_type>(*it);
#endif
            ++it;
          }
          const dof_id_type parent_id =
            (*it == static_cast<T>(-1)) ?
              DofObject::invalid_id :
              libmesh_cast_int<dof_id_type>(*it);
          ++it;
	  const processor_id_type processor_id =
            libmesh_cast_int<processor_id_type>(*it);
          ++it;
	  const subdomain_id_type subdomain_id =
            libmesh_cast_int<subdomain_id_type>(*it);
          ++it;
#ifdef LIBMESH_ENABLE_AMR
	  const unsigned int p_level =
            libmesh_cast_int<unsigned int>(*it);
#endif
	  ++it;

	  Elem *parent =
            (parent_id == DofObject::invalid_id) ? NULL : mesh.elem(parent_id);

	  Elem *elem = Elem::build (elem_type, parent).release();

	  elem->set_id() = e;
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          elem->set_unique_id() = unique_id;
#endif
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
	      const dof_id_type global_node_number = *it;

	      elem->set_node(n) =
		mesh.add_point (Point(), global_node_number);
	    }

          elems_of_dimension[elem->dim()] = true;
	  mesh.add_elem(elem);
	}
    }

  // Set the mesh dimension to the largest encountered for an element
  for (unsigned int i=0; i!=4; ++i)
    if (elems_of_dimension[i])
      mesh.set_mesh_dimension(i);

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



void XdrIO::read_serialized_nodes (Xdr &io, const dof_id_type n_nodes)
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
  std::vector<dof_id_type> needed_nodes; needed_nodes.reserve (mesh.n_nodes());
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
  std::pair<std::vector<dof_id_type>::iterator,
            std::vector<dof_id_type>::iterator> pos;
  pos.first = needed_nodes.begin();

  for (std::size_t blk=0, first_node=0, last_node=0; last_node<n_nodes; blk++)
    {
      first_node = blk*io_blksize;
      last_node  = std::min((blk+1)*io_blksize, std::size_t(n_nodes));

      coords.resize(3*(last_node - first_node));

      if (this->processor_id() == 0)
	io.data_stream (coords.empty() ? NULL : &coords[0], coords.size());

      // For large numbers of processors the majority of processors at any given
      // block may not actually need these data.  It may be worth profiling this,
      // although it is expected that disk IO will be the bottleneck
      this->comm().broadcast (coords);

      for (std::size_t n=first_node, idx=0; n<last_node; n++, idx+=3)
	{
	  // first see if we need this node.  use pos.first as a smart lower
	  // bound, this will ensure that the size of the searched range
	  // decreases as we match nodes.
	  pos = std::equal_range (pos.first, needed_nodes.end(), n);

	  if (pos.first != pos.second) // we need this node.
	    {
	      libmesh_assert_equal_to (*pos.first, n);
	      mesh.node(n) =
		Point (coords[idx+0],
		       coords[idx+1],
		       coords[idx+2]);

	    }
	}
    }
}



template <typename T>
void XdrIO::read_serialized_bcs (Xdr &io, T)
{
  if (this->boundary_condition_file_name() == "n/a") return;

  libmesh_assert (io.reading());

  // convenient reference to our mesh
  MeshBase &mesh = MeshInput<MeshBase>::mesh();

  // and our boundary info object
  BoundaryInfo &boundary_info = *mesh.boundary_info;

  // Version 0.9.2+ introduces unique ids
  read_serialized_bc_names(io, boundary_info, true);  // sideset names

  std::vector<DofBCData> dof_bc_data;
  std::vector<T> input_buffer;

  header_id_type n_bcs=0;
  if (this->processor_id() == 0)
    io.data (n_bcs);
  this->comm().broadcast (n_bcs);

  for (std::size_t blk=0, first_bc=0, last_bc=0; last_bc<n_bcs; blk++)
    {
      first_bc = blk*io_blksize;
      last_bc  = std::min((blk+1)*io_blksize, std::size_t(n_bcs));

      input_buffer.resize (3*(last_bc - first_bc));

      if (this->processor_id() == 0)
	io.data_stream (input_buffer.empty() ? NULL : &input_buffer[0], input_buffer.size());

      this->comm().broadcast (input_buffer);
      dof_bc_data.clear(); /**/ dof_bc_data.reserve (input_buffer.size()/3);

      // convert the input_buffer to DofBCData to facilitate searching
      for (std::size_t idx=0; idx<input_buffer.size(); idx+=3)
	dof_bc_data.push_back (DofBCData(input_buffer[idx+0],
                                          input_buffer[idx+1],
                                          input_buffer[idx+2]));
      input_buffer.clear();
      // note that while the files *we* write should already be sorted by
      // element id this is not necessarily guaranteed.
      std::sort (dof_bc_data.begin(), dof_bc_data.end());

      MeshBase::const_element_iterator
	it  = mesh.level_elements_begin(0),
	end = mesh.level_elements_end(0);

      // Look for BCs in this block for all the level-0 elements we have
      // (not just local ones).  Do this by finding all the entries
      // in dof_bc_data whose elem_id match the ID of the current element.
      // We cannot rely on NULL neighbors at this point since the neighbor
      // data structure has not been initialized.
      for (std::pair<std::vector<DofBCData>::iterator,
	             std::vector<DofBCData>::iterator> pos; it!=end; ++it)
#if defined(__SUNPRO_CC) || defined(__PGI)
	for (pos = std::equal_range (dof_bc_data.begin(), dof_bc_data.end(), (*it)->id(), CompareIntDofBCData());
	     pos.first != pos.second; ++pos.first)
#else
	for (pos = std::equal_range (dof_bc_data.begin(), dof_bc_data.end(), (*it)->id());
	     pos.first != pos.second; ++pos.first)
#endif
	  {
	    libmesh_assert_equal_to (pos.first->dof_id, (*it)->id());
	    libmesh_assert_less (pos.first->side, (*it)->n_sides());

	    boundary_info.add_side (*it, pos.first->side, pos.first->bc_id);
	  }
    }
}



template <typename T>
void XdrIO::read_serialized_nodesets (Xdr &io, T)
{
  if (this->boundary_condition_file_name() == "n/a") return;

  libmesh_assert (io.reading());

  // convenient reference to our mesh
  MeshBase &mesh = MeshInput<MeshBase>::mesh();

  // and our boundary info object
  BoundaryInfo &boundary_info = *mesh.boundary_info;

  // Version 0.9.2+ introduces unique ids
  read_serialized_bc_names(io, boundary_info, false); // nodeset names

  // TODO: Make a data object that works with both the element and nodal bcs
  std::vector<DofBCData> node_bc_data;
  std::vector<T> input_buffer;

  header_id_type n_nodesets=0;
  if (this->processor_id() == 0)
    io.data (n_nodesets);
  this->comm().broadcast (n_nodesets);

  for (std::size_t blk=0, first_bc=0, last_bc=0; last_bc<n_nodesets; blk++)
    {
      first_bc = blk*io_blksize;
      last_bc  = std::min((blk+1)*io_blksize, std::size_t(n_nodesets));

      input_buffer.resize (2*(last_bc - first_bc));

      if (this->processor_id() == 0)
	io.data_stream (input_buffer.empty() ? NULL : &input_buffer[0], input_buffer.size());

      this->comm().broadcast (input_buffer);
      node_bc_data.clear(); /**/ node_bc_data.reserve (input_buffer.size()/2);

      // convert the input_buffer to DofBCData to facilitate searching
      for (std::size_t idx=0; idx<input_buffer.size(); idx+=2)
	node_bc_data.push_back (DofBCData(input_buffer[idx+0],
					   0,
					   input_buffer[idx+1]));
      input_buffer.clear();
      // note that while the files *we* write should already be sorted by
      // node id this is not necessarily guaranteed.
      std::sort (node_bc_data.begin(), node_bc_data.end());

      MeshBase::const_node_iterator
	it  = mesh.nodes_begin(),
	end = mesh.nodes_end();

      // Look for BCs in this block for all nodes we have
      // (not just local ones).  Do this by finding all the entries
      // in node_bc_data whose dof_id(node_id)  match the ID of the current node.
      for (std::pair<std::vector<DofBCData>::iterator,
	             std::vector<DofBCData>::iterator> pos; it!=end; ++it)
#if defined(__SUNPRO_CC) || defined(__PGI)
	for (pos = std::equal_range (node_bc_data.begin(), node_bc_data.end(), (*it)->id(), CompareIntDofBCData());
	     pos.first != pos.second; ++pos.first)
#else
	for (pos = std::equal_range (node_bc_data.begin(), node_bc_data.end(), (*it)->id());
	     pos.first != pos.second; ++pos.first)
#endif
	  {
            // Note: dof_id from ElmeBCData is being used to hold node_id here
	    libmesh_assert_equal_to (pos.first->dof_id, (*it)->id());

	    boundary_info.add_node (*it, pos.first->bc_id);
	  }
    }
}



void XdrIO::read_serialized_bc_names(Xdr &io, BoundaryInfo & info, bool is_sideset)
{
  const bool read_entity_info = this->version().find("0.9.2") != std::string::npos ? true : false;
  if (read_entity_info)
  {
    header_id_type n_boundary_names = 0;
    std::vector<header_id_type> boundary_ids;
    std::vector<std::string>  boundary_names;

    // Read the sideset names
    if (this->processor_id() == 0)
    {
      io.data(n_boundary_names);

      boundary_ids.resize(n_boundary_names);
      boundary_names.resize(n_boundary_names);

      if (n_boundary_names)
      {
        io.data(boundary_ids);
        io.data(boundary_names);
      }
    }

    // Broadcast the boundary names to all processors
    this->comm().broadcast(n_boundary_names);
    if (n_boundary_names == 0)
      return;

    boundary_ids.resize(n_boundary_names);
    boundary_names.resize(n_boundary_names);
    this->comm().broadcast(boundary_ids);
    this->comm().broadcast(boundary_names);

    // Reassemble the named boundary information
    std::map<boundary_id_type, std::string> & boundary_map = is_sideset ?
      info.set_sideset_name_map() : info.set_nodeset_name_map();

    for (unsigned int i=0; i<n_boundary_names; ++i)
      boundary_map.insert(std::make_pair(boundary_ids[i], boundary_names[i]));
  }
}



void XdrIO::pack_element (std::vector<xdr_id_type> &conn, const Elem *elem,
			  const dof_id_type parent_id, const dof_id_type parent_pid) const
{
  libmesh_assert(elem);
  libmesh_assert_equal_to (elem->n_nodes(), Elem::type_to_n_nodes_map[elem->type()]);

  conn.push_back(elem->n_nodes());

  conn.push_back (elem->type());

  // In version 0.7.0+ "id" is stored but it not used.  In version 0.9.2+
  // we will store unique_id instead, therefore there is no need to
  // check for the older version when writing the unique_id.
  conn.push_back (elem->unique_id());

  if (parent_id != DofObject::invalid_id)
    {
      conn.push_back (parent_id);
      libmesh_assert_not_equal_to (parent_pid, DofObject::invalid_id);
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

} // namespace libMesh
