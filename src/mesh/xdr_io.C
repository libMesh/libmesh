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
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <string>

// Local includes
#include "libmesh/xdr_io.h"
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
  dof_id_type        dof_id;
  unsigned short int side;
  boundary_id_type   bc_id;

  // Default constructor
  DofBCData (dof_id_type        dof_id_in=0,
             unsigned short int side_in=0,
             boundary_id_type   bc_id_in=0) :
    dof_id(dof_id_in),
    side(side_in),
    bc_id(bc_id_in)
  {}

  // comparison operator
  bool operator < (const DofBCData & other) const
  {
    if (this->dof_id == other.dof_id)
      return (this->side < other.side);

    return this->dof_id < other.dof_id;
  }
};

// comparison operator
bool operator < (const unsigned int & other_dof_id,
                 const DofBCData & dof_bc)
{
  return other_dof_id < dof_bc.dof_id;
}

bool operator < (const DofBCData & dof_bc,
                 const unsigned int & other_dof_id)

{
  return dof_bc.dof_id < other_dof_id;
}


// For some reason SunStudio does not seem to accept the above
// comparison functions for use in
// std::equal_range (ElemBCData::iterator, ElemBCData::iterator, unsigned int);
#if defined(__SUNPRO_CC) || defined(__PGI)
struct CompareIntDofBCData
{
  bool operator()(const unsigned int & other_dof_id,
                  const DofBCData & dof_bc)
  {
    return other_dof_id < dof_bc.dof_id;
  }

  bool operator()(const DofBCData & dof_bc,
                  const unsigned int & other_dof_id)
  {
    return dof_bc.dof_id < other_dof_id;
  }
};
#endif

template <class T, class U>
struct libmesh_type_is_same {
  static const bool value = false;
};

template <class T>
struct libmesh_type_is_same<T, T> {
  static const bool value = true;
};

}



// ------------------------------------------------------------
// XdrIO static data
const std::size_t XdrIO::io_blksize = 128000;



// ------------------------------------------------------------
// XdrIO members
XdrIO::XdrIO (MeshBase & mesh, const bool binary_in) :
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
  _field_width        (4),   // In 0.7.0, all fields are 4 bytes, in 0.9.2+ they can vary
  _version            ("libMesh-1.3.0"),
  _bc_file_name       ("n/a"),
  _partition_map_file ("n/a"),
  _subdomain_map_file ("n/a"),
  _p_level_file       ("n/a")
{
}



XdrIO::XdrIO (const MeshBase & mesh, const bool binary_in) :
  MeshOutput<MeshBase>(mesh,/* is_parallel_format = */ true),
  ParallelObject      (mesh),
  _binary (binary_in)
{
}



XdrIO::~XdrIO ()
{
}



void XdrIO::write (const std::string & name)
{
  if (this->legacy())
    libmesh_error_msg("We don't support writing parallel files in the legacy format.");

  Xdr io ((this->processor_id() == 0) ? name : "", this->binary() ? ENCODE : WRITE);

  START_LOG("write()","XdrIO");

  // convenient reference to our mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  new_header_id_type n_elem  = mesh.n_elem();
  new_header_id_type n_nodes = mesh.n_nodes();

  libmesh_assert(n_elem == mesh.n_elem());
  libmesh_assert(n_nodes == mesh.n_nodes());

  new_header_id_type n_side_bcs = mesh.get_boundary_info().n_boundary_conds();
  new_header_id_type n_edge_bcs = mesh.get_boundary_info().n_edge_conds();
  new_header_id_type n_shellface_bcs = mesh.get_boundary_info().n_shellface_conds();
  new_header_id_type n_nodesets = mesh.get_boundary_info().n_nodeset_conds();
  unsigned int n_p_levels = MeshTools::n_p_levels (mesh);

  bool write_parallel_files = this->write_parallel();

  //-------------------------------------------------------------
  // For all the optional files -- the default file name is "n/a".
  // However, the user may specify an optional external file.

  // If there are BCs and the user has not already provided a
  // file name then write to "."
  if ((n_side_bcs || n_edge_bcs || n_shellface_bcs || n_nodesets) &&
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
      new_header_id_type write_size = sizeof(xdr_id_type), zero_size = 0;

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
      io.data (write_bcs          ? write_size : zero_size, "# side size");  // side number
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

      // write the side boundary condition information
      this->write_serialized_side_bcs (io, n_side_bcs);

      // write the nodeset information
      this->write_serialized_nodesets (io, n_nodesets);

      // write the edge boundary condition information
      this->write_serialized_edge_bcs (io, n_edge_bcs);

      // write the "shell face" boundary condition information
      this->write_serialized_shellface_bcs (io, n_shellface_bcs);
    }
  else
    {
      // write subdomain names
      this->write_serialized_subdomain_names(io);

      // write connectivity
      this->write_serialized_connectivity (io, n_elem);

      // write the nodal locations
      this->write_serialized_nodes (io, n_nodes);

      // write the side boundary condition information
      this->write_serialized_side_bcs (io, n_side_bcs);

      // write the nodeset information
      this->write_serialized_nodesets (io, n_nodesets);

      // write the edge boundary condition information
      this->write_serialized_edge_bcs (io, n_edge_bcs);

      // write the "shell face" boundary condition information
      this->write_serialized_shellface_bcs (io, n_shellface_bcs);
    }

  STOP_LOG("write()","XdrIO");

  // pause all processes until the writing ends -- this will
  // protect for the pathological case where a write is
  // followed immediately by a read.  The write must be
  // guaranteed to complete first.
  io.close();
  this->comm().barrier();
}



void XdrIO::write_serialized_subdomain_names(Xdr & io) const
{
  if (this->processor_id() == 0)
    {
      const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

      const std::map<subdomain_id_type, std::string> & subdomain_map = mesh.get_subdomain_name_map();

      std::vector<new_header_id_type> subdomain_ids;
      subdomain_ids.reserve(subdomain_map.size());

      std::vector<std::string> subdomain_names;
      subdomain_names.reserve(subdomain_map.size());

      // We need to loop over the map and make sure that there aren't any invalid entries.  Since we
      // return writable references in mesh_base, it's possible for the user to leave some entity names
      // blank.  We can't write those to the XDA file.
      new_header_id_type n_subdomain_names = 0;
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



void XdrIO::write_serialized_connectivity (Xdr & io, const dof_id_type libmesh_dbg_var(n_elem)) const
{
  libmesh_assert (io.writing());

  const bool
    write_p_level      = ("." == this->polynomial_level_file_name()),
    write_partitioning = ("." == this->partition_map_file_name()),
    write_subdomain_id = ("." == this->subdomain_map_file_name());

  // convenient reference to our mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();
  libmesh_assert_equal_to (n_elem, mesh.n_elem());

  // We will only write active elements and their parents.
  const unsigned int n_active_levels = MeshTools::n_active_levels (mesh);
  std::vector<xdr_id_type> n_global_elem_at_level(n_active_levels);

  MeshBase::const_element_iterator it  = mesh.local_elements_end(), end=it;

  // Find the number of local and global elements at each level
#ifndef NDEBUG
  xdr_id_type tot_n_elem = 0;
#endif
  for (unsigned int level=0; level<n_active_levels; level++)
    {
      it  = mesh.local_level_elements_begin(level);
      end = mesh.local_level_elements_end(level);

      n_global_elem_at_level[level] = MeshTools::n_elem(it, end);

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

#ifdef LIBMESH_ENABLE_AMR
  typedef std::map<dof_id_type, std::pair<processor_id_type, dof_id_type>> id_map_type;
  id_map_type parent_id_map, child_id_map;
#endif

  dof_id_type my_next_elem=0, next_global_elem=0;

  //-------------------------------------------
  // First write the level-0 elements directly.
  it  = mesh.local_level_elements_begin(0);
  end = mesh.local_level_elements_end(0);
  for (; it != end; ++it, ++my_next_elem)
    {
      pack_element (xfer_conn, *it);
#ifdef LIBMESH_ENABLE_AMR
      parent_id_map[(*it)->id()] = std::make_pair(this->processor_id(),
                                                  my_next_elem);
#endif
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

                // n. nodes
                const xdr_id_type n_nodes = *recv_conn_iter;
                ++recv_conn_iter;

                // type
                output_buffer.push_back(*recv_conn_iter);
                ++recv_conn_iter;

                // unique_id
                if (_write_unique_id)
                  output_buffer.push_back(*recv_conn_iter);
                ++recv_conn_iter;

                // processor id
                if (write_partitioning)
                  output_buffer.push_back(*recv_conn_iter);
                ++recv_conn_iter;

                // subdomain id
                if (write_subdomain_id)
                  output_buffer.push_back(*recv_conn_iter);
                ++recv_conn_iter;

#ifdef LIBMESH_ENABLE_AMR
                // p level
                if (write_p_level)
                  output_buffer.push_back(*recv_conn_iter);
                ++recv_conn_iter;
#endif
                for (dof_id_type node=0; node<n_nodes; node++, ++recv_conn_iter)
                  output_buffer.push_back(*recv_conn_iter);

                io.data_stream
                  (&output_buffer[0],
                   cast_int<unsigned int>(output_buffer.size()),
                   cast_int<unsigned int>(output_buffer.size()));
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
            const Elem * parent = *it;
            id_map_type::iterator pos = parent_id_map.find(parent->id());
            libmesh_assert (pos != parent_id_map.end());
            const processor_id_type parent_pid = pos->second.first;
            const dof_id_type parent_id  = pos->second.second;
            parent_id_map.erase(pos);

            for (auto & child : parent->child_ref_range())
              {
                pack_element (xfer_conn, &child, parent_id, parent_pid);

                // this aproach introduces the possibility that we write
                // non-local elements.  These elements may well be parents
                // at the next step
                child_id_map[child.id()] = std::make_pair (child.processor_id(),
                                                           my_n_elem_written_at_level++);
                my_next_elem++;
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
            std::sprintf(buf, "# n_elem at level %u", level);
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

                    // n. nodes
                    const xdr_id_type n_nodes = *recv_conn_iter;
                    ++recv_conn_iter;

                    // type
                    output_buffer.push_back(*recv_conn_iter);
                    ++recv_conn_iter;

                    // unique_id
                    if (_write_unique_id)
                      output_buffer.push_back(*recv_conn_iter);
                    ++recv_conn_iter;

                    // parent local id
                    const xdr_id_type parent_local_id = *recv_conn_iter;
                    ++recv_conn_iter;

                    // parent processor id
                    const xdr_id_type parent_pid = *recv_conn_iter;
                    ++recv_conn_iter;

                    output_buffer.push_back (parent_local_id+processor_offsets[parent_pid]);

                    // processor id
                    if (write_partitioning)
                      output_buffer.push_back(*recv_conn_iter);
                    ++recv_conn_iter;

                    // subdomain id
                    if (write_subdomain_id)
                      output_buffer.push_back(*recv_conn_iter);
                    ++recv_conn_iter;

                    // p level
                    if (write_p_level)
                      output_buffer.push_back(*recv_conn_iter);
                    ++recv_conn_iter;

                    for (xdr_id_type node=0; node<n_nodes; node++, ++recv_conn_iter)
                      output_buffer.push_back(*recv_conn_iter);

                    io.data_stream
                      (&output_buffer[0],
                       cast_int<unsigned int>(output_buffer.size()),
                       cast_int<unsigned int>(output_buffer.size()));
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
        std::vector<std::vector<dof_id_type>> requested_ids(this->n_processors());
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
        parent_id_map.swap(child_id_map);
        child_id_map.clear();
      }
    }
#endif // LIBMESH_ENABLE_AMR
  if (this->processor_id() == 0)
    libmesh_assert_equal_to (next_global_elem, n_elem);

}



void XdrIO::write_serialized_nodes (Xdr & io, const dof_id_type n_nodes) const
{
  // convenient reference to our mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();
  libmesh_assert_equal_to (n_nodes, mesh.n_nodes());

  std::vector<dof_id_type> xfer_ids;
  std::vector<Real> xfer_coords;
  std::vector<Real> & coords=xfer_coords;

  std::vector<std::vector<dof_id_type>> recv_ids   (this->n_processors());
  std::vector<std::vector<Real>>         recv_coords(this->n_processors());

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  std::vector<xdr_id_type> xfer_unique_ids;
  std::vector<xdr_id_type> & unique_ids=xfer_unique_ids;
  std::vector<std::vector<xdr_id_type>> recv_unique_ids (this->n_processors());
#endif // LIBMESH_ENABLE_UNIQUE_ID

  std::size_t n_written=0;

  for (std::size_t blk=0, last_node=0; last_node<n_nodes; blk++)
    {
      const std::size_t first_node = blk*io_blksize;
      last_node = std::min((blk+1)*io_blksize, std::size_t(n_nodes));

      // Build up the xfer buffers on each processor
      xfer_ids.clear();
      xfer_coords.clear();
#ifdef LIBMESH_ENABLE_UNIQUE_ID
      xfer_unique_ids.clear();
#endif // LIBMESH_ENABLE_UNIQUE_ID

      for (const auto & node : mesh.local_node_ptr_range())
        if ((node->id() >= first_node) && // node in [first_node, last_node)
            (node->id() <  last_node))
          {
            xfer_ids.push_back(node->id());
#ifdef LIBMESH_ENABLE_UNIQUE_ID
            xfer_unique_ids.push_back(node->unique_id());
#endif // LIBMESH_ENABLE_UNIQUE_ID
            xfer_coords.push_back((*node)(0));
#if LIBMESH_DIM > 1
            xfer_coords.push_back((*node)(1));
#endif
#if LIBMESH_DIM > 2
            xfer_coords.push_back((*node)(2));
#endif
          }

      //-------------------------------------
      // Send the xfer buffers to processor 0
      std::vector<std::size_t> ids_size;

      const std::size_t my_ids_size = xfer_ids.size();

      // explicitly gather ids_size
      this->comm().gather (0, my_ids_size, ids_size);

      // We will have lots of simultaneous receives if we are
      // processor 0, so let's use nonblocking receives.
      std::vector<Parallel::Request>
#ifdef LIBMESH_ENABLE_UNIQUE_ID
        unique_id_request_handles(this->n_processors()-1),
#endif // LIBMESH_ENABLE_UNIQUE_ID
        id_request_handles(this->n_processors()-1),
        coord_request_handles(this->n_processors()-1);

      Parallel::MessageTag
#ifdef LIBMESH_ENABLE_UNIQUE_ID
        unique_id_tag = mesh.comm().get_unique_tag(1233),
#endif // LIBMESH_ENABLE_UNIQUE_ID
        id_tag    = mesh.comm().get_unique_tag(1234),
        coord_tag = mesh.comm().get_unique_tag(1235);

      // Post the receives -- do this on processor 0 only.
      if (this->processor_id() == 0)
        {
          for (unsigned int pid=0; pid<this->n_processors(); pid++)
            {
              recv_ids[pid].resize(ids_size[pid]);
              recv_coords[pid].resize(ids_size[pid]*LIBMESH_DIM);
#ifdef LIBMESH_ENABLE_UNIQUE_ID
              recv_unique_ids[pid].resize(ids_size[pid]);
#endif // LIBMESH_ENABLE_UNIQUE_ID

              if (pid == 0)
                {
                  recv_ids[0] = xfer_ids;
                  recv_coords[0] = xfer_coords;
#ifdef LIBMESH_ENABLE_UNIQUE_ID
                  recv_unique_ids[0] = xfer_unique_ids;
#endif // LIBMESH_ENABLE_UNIQUE_ID
                }
              else
                {
                  this->comm().receive (pid, recv_ids[pid],
                                        id_request_handles[pid-1],
                                        id_tag);
                  this->comm().receive (pid, recv_coords[pid],
                                        coord_request_handles[pid-1],
                                        coord_tag);
#ifdef LIBMESH_ENABLE_UNIQUE_ID
                  this->comm().receive (pid, recv_unique_ids[pid],
                                        unique_id_request_handles[pid-1],
                                        unique_id_tag);
#endif // LIBMESH_ENABLE_UNIQUE_ID
                }
            }
        }
      else
        {
          // Send -- do this on all other processors.
          this->comm().send(0, xfer_ids,    id_tag);
          this->comm().send(0, xfer_coords, coord_tag);
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          this->comm().send(0, xfer_unique_ids, unique_id_tag);
#endif // LIBMESH_ENABLE_UNIQUE_ID
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
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          Parallel::wait (unique_id_request_handles);
#endif // LIBMESH_ENABLE_UNIQUE_ID

          // Write the coordinates in this block.
          std::size_t tot_id_size=0;

          for (unsigned int pid=0; pid<this->n_processors(); pid++)
            {
              tot_id_size    += recv_ids[pid].size();
              libmesh_assert_equal_to(recv_coords[pid].size(),
                                      recv_ids[pid].size()*LIBMESH_DIM);
#ifdef LIBMESH_ENABLE_UNIQUE_ID
              libmesh_assert_equal_to
                (recv_ids[pid].size(), recv_unique_ids[pid].size());
#endif // LIBMESH_ENABLE_UNIQUE_ID
            }

          libmesh_assert_less_equal
            (tot_id_size, std::min(io_blksize, std::size_t(n_nodes)));

          coords.resize (3*tot_id_size);
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          unique_ids.resize(tot_id_size);
#endif

          for (unsigned int pid=0; pid<this->n_processors(); pid++)
            for (std::size_t idx=0; idx<recv_ids[pid].size(); idx++)
              {
                const std::size_t local_idx = recv_ids[pid][idx] - first_node;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
                libmesh_assert_less (local_idx, unique_ids.size());

                unique_ids[local_idx] = recv_unique_ids[pid][idx];
#endif

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

          io.data_stream (coords.empty() ? libmesh_nullptr : &coords[0],
                          cast_int<unsigned int>(coords.size()), 3);
        }
    }

  if (this->processor_id() == 0)
    libmesh_assert_equal_to (n_written, n_nodes);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  // XDR unsigned char doesn't work as anticipated
  unsigned short write_unique_ids = 1;
#else
  unsigned short write_unique_ids = 0;
#endif
  if (this->processor_id() == 0)
    io.data (write_unique_ids, "# presence of unique ids");

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  n_written = 0;

  for (std::size_t blk=0, last_node=0; last_node<n_nodes; blk++)
    {
      const std::size_t first_node = blk*io_blksize;
      last_node = std::min((blk+1)*io_blksize, std::size_t(n_nodes));

      // Build up the xfer buffers on each processor
      xfer_ids.clear();
      xfer_unique_ids.clear();

      for (const auto & node : mesh.local_node_ptr_range())
        if ((node->id() >= first_node) && // node in [first_node, last_node)
            (node->id() <  last_node))
          {
            xfer_ids.push_back(node->id());
            xfer_unique_ids.push_back(node->unique_id());
          }

      //-------------------------------------
      // Send the xfer buffers to processor 0
      std::vector<std::size_t> ids_size;

      const std::size_t my_ids_size = xfer_ids.size();

      // explicitly gather ids_size
      this->comm().gather (0, my_ids_size, ids_size);

      // We will have lots of simultaneous receives if we are
      // processor 0, so let's use nonblocking receives.
      std::vector<Parallel::Request>
        unique_id_request_handles(this->n_processors()-1),
        id_request_handles(this->n_processors()-1);

      Parallel::MessageTag
        unique_id_tag = mesh.comm().get_unique_tag(1236),
        id_tag    = mesh.comm().get_unique_tag(1237);

      // Post the receives -- do this on processor 0 only.
      if (this->processor_id() == 0)
        {
          for (unsigned int pid=0; pid<this->n_processors(); pid++)
            {
              recv_ids[pid].resize(ids_size[pid]);
              recv_unique_ids[pid].resize(ids_size[pid]);

              if (pid == 0)
                {
                  recv_ids[0] = xfer_ids;
                  recv_unique_ids[0] = xfer_unique_ids;
                }
              else
                {
                  this->comm().receive (pid, recv_ids[pid],
                                        id_request_handles[pid-1],
                                        id_tag);
                  this->comm().receive (pid, recv_unique_ids[pid],
                                        unique_id_request_handles[pid-1],
                                        unique_id_tag);
                }
            }
        }
      else
        {
          // Send -- do this on all other processors.
          this->comm().send(0, xfer_ids,    id_tag);
          this->comm().send(0, xfer_unique_ids, unique_id_tag);
        }

      // -------------------------------------------------------
      // Receive the messages and write the output on processor 0.
      if (this->processor_id() == 0)
        {
          // Wait for all the receives to complete. We have no
          // need for the statuses since we already know the
          // buffer sizes.
          Parallel::wait (id_request_handles);
          Parallel::wait (unique_id_request_handles);

          // Write the unique ids in this block.
          std::size_t tot_id_size=0;

          for (unsigned int pid=0; pid<this->n_processors(); pid++)
            {
              tot_id_size    += recv_ids[pid].size();
              libmesh_assert_equal_to
                (recv_ids[pid].size(), recv_unique_ids[pid].size());
            }

          libmesh_assert_less_equal
            (tot_id_size, std::min(io_blksize, std::size_t(n_nodes)));

          unique_ids.resize(tot_id_size);

          for (unsigned int pid=0; pid<this->n_processors(); pid++)
            for (std::size_t idx=0; idx<recv_ids[pid].size(); idx++)
              {
                const std::size_t local_idx = recv_ids[pid][idx] - first_node;

                libmesh_assert_less (local_idx, unique_ids.size());

                unique_ids[local_idx] = recv_unique_ids[pid][idx];

                n_written++;
              }

          io.data_stream (unique_ids.empty() ? libmesh_nullptr : &unique_ids[0],
                          cast_int<unsigned int>(unique_ids.size()), 1);
        }
    }

  if (this->processor_id() == 0)
    libmesh_assert_equal_to (n_written, n_nodes);

#endif // LIBMESH_ENABLE_UNIQUE_ID
}



void XdrIO::write_serialized_bcs_helper (Xdr & io, const new_header_id_type n_bcs, const std::string bc_type) const
{
  libmesh_assert (io.writing());

  // convenient reference to our mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // and our boundary info object
  const BoundaryInfo & boundary_info = mesh.get_boundary_info();

  // Version 0.9.2+ introduces entity names
  write_serialized_bc_names(io, boundary_info, true);  // sideset names

  new_header_id_type n_bcs_out = n_bcs;
  if (this->processor_id() == 0)
    {
      std::stringstream comment_string;
      comment_string << "# number of " << bc_type << " boundary conditions";
      io.data (n_bcs_out, comment_string.str().c_str());
    }
  n_bcs_out = 0;

  if (!n_bcs) return;

  std::vector<xdr_id_type> xfer_bcs, recv_bcs;
  std::vector<std::size_t> bc_sizes(this->n_processors());

  // Boundary conditions are only specified for level-0 elements
  MeshBase::const_element_iterator
    it  = mesh.local_level_elements_begin(0),
    end = mesh.local_level_elements_end(0);

  // Container to catch boundary IDs handed back by BoundaryInfo
  std::vector<boundary_id_type> bc_ids;

  dof_id_type n_local_level_0_elem=0;
  for (; it!=end; ++it, n_local_level_0_elem++)
    {
      const Elem * elem = *it;

      if (bc_type == "side")
        {
          for (auto s : elem->side_index_range())
            {
              boundary_info.boundary_ids (elem, s, bc_ids);
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
      else if (bc_type == "edge")
        {
          for (auto e : elem->edge_index_range())
            {
              boundary_info.edge_boundary_ids (elem, e, bc_ids);
              for (std::vector<boundary_id_type>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
                {
                  const boundary_id_type bc_id = *id_it;
                  if (bc_id != BoundaryInfo::invalid_id)
                    {
                      xfer_bcs.push_back (n_local_level_0_elem);
                      xfer_bcs.push_back (e) ;
                      xfer_bcs.push_back (bc_id);
                    }
                }
            }
        }
      else if (bc_type == "shellface")
        {
          for (unsigned short sf=0; sf<2; sf++)
            {
              boundary_info.shellface_boundary_ids (elem, sf, bc_ids);
              for (std::vector<boundary_id_type>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
                {
                  const boundary_id_type bc_id = *id_it;
                  if (bc_id != BoundaryInfo::invalid_id)
                    {
                      xfer_bcs.push_back (n_local_level_0_elem);
                      xfer_bcs.push_back (sf) ;
                      xfer_bcs.push_back (bc_id);
                    }
                }
            }
        }
      else
        {
          libmesh_error_msg("bc_type not recognized: " + bc_type);
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
            = cast_int<dof_id_type>(recv_bcs.back());
          recv_bcs.pop_back();

          for (std::size_t idx=0; idx<recv_bcs.size(); idx += 3, n_bcs_out++)
            recv_bcs[idx+0] += elem_offset;

          io.data_stream (recv_bcs.empty() ? libmesh_nullptr : &recv_bcs[0],
                          cast_int<unsigned int>(recv_bcs.size()), 3);
          elem_offset += my_n_local_level_0_elem;
        }
      libmesh_assert_equal_to (n_bcs, n_bcs_out);
    }
  else
    this->comm().send (0, xfer_bcs);
}



void XdrIO::write_serialized_side_bcs (Xdr & io, const new_header_id_type n_side_bcs) const
{
  write_serialized_bcs_helper(io, n_side_bcs, "side");
}



void XdrIO::write_serialized_edge_bcs (Xdr & io, const new_header_id_type n_edge_bcs) const
{
  write_serialized_bcs_helper(io, n_edge_bcs, "edge");
}



void XdrIO::write_serialized_shellface_bcs (Xdr & io, const new_header_id_type n_shellface_bcs) const
{
  write_serialized_bcs_helper(io, n_shellface_bcs, "shellface");
}



void XdrIO::write_serialized_nodesets (Xdr & io, const new_header_id_type n_nodesets) const
{
  libmesh_assert (io.writing());

  // convenient reference to our mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // and our boundary info object
  const BoundaryInfo & boundary_info = mesh.get_boundary_info();

  // Version 0.9.2+ introduces entity names
  write_serialized_bc_names(io, boundary_info, false);  // nodeset names

  new_header_id_type n_nodesets_out = n_nodesets;
  if (this->processor_id() == 0)
    io.data (n_nodesets_out, "# number of nodesets");
  n_nodesets_out = 0;

  if (!n_nodesets) return;

  std::vector<xdr_id_type> xfer_bcs, recv_bcs;
  std::vector<std::size_t> bc_sizes(this->n_processors());

  // Container to catch boundary IDs handed back by BoundaryInfo
  std::vector<boundary_id_type> nodeset_ids;

  dof_id_type n_node=0;
  for (const auto & node : mesh.local_node_ptr_range())
    {
      boundary_info.boundary_ids (node, nodeset_ids);
      for (std::vector<boundary_id_type>::const_iterator id_it=nodeset_ids.begin(); id_it!=nodeset_ids.end(); ++id_it)
        {
          const boundary_id_type bc_id = *id_it;
          if (bc_id != BoundaryInfo::invalid_id)
            {
              xfer_bcs.push_back (node->id());
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

          const dof_id_type my_n_node =
            cast_int<dof_id_type>(recv_bcs.back());
          recv_bcs.pop_back();

          for (std::size_t idx=0; idx<recv_bcs.size(); idx += 2, n_nodesets_out++)
            recv_bcs[idx+0] += node_offset;

          io.data_stream (recv_bcs.empty() ? libmesh_nullptr : &recv_bcs[0],
                          cast_int<unsigned int>(recv_bcs.size()), 2);
          node_offset += my_n_node;
        }
      libmesh_assert_equal_to (n_nodesets, n_nodesets_out);
    }
  else
    this->comm().send (0, xfer_bcs);
}



void XdrIO::write_serialized_bc_names (Xdr & io, const BoundaryInfo & info, bool is_sideset) const
{
  if (this->processor_id() == 0)
    {
      const std::map<boundary_id_type, std::string> & boundary_map = is_sideset ?
        info.get_sideset_name_map() : info.get_nodeset_name_map();

      std::vector<new_header_id_type> boundary_ids;
      boundary_ids.reserve(boundary_map.size());

      std::vector<std::string>  boundary_names;
      boundary_names.reserve(boundary_map.size());

      // We need to loop over the map and make sure that there aren't any invalid entries.  Since we
      // return writable references in boundary_info, it's possible for the user to leave some entity names
      // blank.  We can't write those to the XDA file.
      new_header_id_type n_boundary_names = 0;
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



void XdrIO::read (const std::string & name)
{
  LOG_SCOPE("read()","XdrIO");

  // Only open the file on processor 0 -- this is especially important because
  // there may be an underlying bzip/bunzip going on, and multiple simultaneous
  // calls will produce a race condition.
  Xdr io (this->processor_id() == 0 ? name : "", this->binary() ? DECODE : READ);

  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // get the version string.
  if (this->processor_id() == 0)
    io.data (this->version());
  this->comm().broadcast (this->version());

  // note that for "legacy" files the first entry is an
  // integer -- not a string at all.
  this->legacy() = !(this->version().find("libMesh") < this->version().size());

  // Check for a legacy version format.
  if (this->legacy())
    libmesh_error_msg("We no longer support reading files in the legacy format.");

  // Read headers with the old id type if they're pre-1.3.0, or with
  // the new id type if they're post-1.3.0
  std::vector<new_header_id_type> meta_data(10, sizeof(xdr_id_type));
  if (this->version_at_least_1_3_0())
    {
      this->read_header(io, meta_data);
    }
  else
    {
      std::vector<old_header_id_type> old_data(10, sizeof(xdr_id_type));

      this->read_header(io, old_data);

      meta_data.assign(old_data.begin(), old_data.end());
    }

  const new_header_id_type & n_elem = meta_data[0];
  const new_header_id_type & n_nodes = meta_data[1];

  /**
   * We are future proofing the layout of this file by adding in size information for all stored types.
   * TODO: All types are stored as the same size. Use the size information to pack things efficiently.
   * For now we will assume that "type size" is how the entire file will be encoded.
   */
  if (version_at_least_0_9_2())
    _field_width = meta_data[2];

  // On systems where uint64_t==unsigned long, we were previously
  // writing 64-bit unsigned integers via xdr_u_long(), a function
  // which is literally less suited for that task than abort() would
  // have been, because at least abort() would have *known* it
  // couldn't write rather than truncating writes to 32 bits.
  //
  // If we have files with version < 1.3.0, then we'll continue to use
  // 32 bit field width, regardless of whether the file thinks we
  // should, whenever we're on a system where the problem would have
  // occurred.
  if ((_field_width == 4) ||
      (!version_at_least_1_3_0() &&
       libmesh_type_is_same<uint64_t, unsigned long>::value))
    {
      uint32_t type_size = 0;

      // read subdomain names
      this->read_serialized_subdomain_names(io);

      // read connectivity
      this->read_serialized_connectivity (io, n_elem, meta_data, type_size);

      // read the nodal locations
      this->read_serialized_nodes (io, n_nodes);

      // read the side boundary conditions
      this->read_serialized_side_bcs (io, type_size);

      if (version_at_least_0_9_2())
        // read the nodesets
        this->read_serialized_nodesets (io, type_size);

      if (version_at_least_1_1_0())
        {
          // read the edge boundary conditions
          this->read_serialized_edge_bcs (io, type_size);

          // read the "shell face" boundary conditions
          this->read_serialized_shellface_bcs (io, type_size);
        }
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
      this->read_serialized_side_bcs (io, type_size);

      if (version_at_least_0_9_2())
        // read the nodesets
        this->read_serialized_nodesets (io, type_size);

      if (version_at_least_1_1_0())
        {
          // read the edge boundary conditions
          this->read_serialized_edge_bcs (io, type_size);

          // read the "shell face" boundary conditions
          this->read_serialized_shellface_bcs (io, type_size);
        }
    }

  // set the node processor ids
  Partitioner::set_node_processor_ids(mesh);
}



template <typename T>
void XdrIO::read_header (Xdr & io, std::vector<T> & meta_data)
{
  LOG_SCOPE("read_header()","XdrIO");

  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  if (this->processor_id() == 0)
    {
      unsigned int pos=0;

      io.data (meta_data[pos++]);
      io.data (meta_data[pos++]);
      io.data (this->boundary_condition_file_name()); // libMesh::out << "bc_file="  << this->boundary_condition_file_name() << std::endl;
      io.data (this->subdomain_map_file_name());      // libMesh::out << "sid_file=" << this->subdomain_map_file_name()      << std::endl;
      io.data (this->partition_map_file_name());      // libMesh::out << "pid_file=" << this->partition_map_file_name()      << std::endl;
      io.data (this->polynomial_level_file_name());   // libMesh::out << "pl_file="  << this->polynomial_level_file_name()   << std::endl;

      if (version_at_least_0_9_2())
        {
          io.data (meta_data[pos++], "# type size");
          io.data (meta_data[pos++], "# uid size");
          io.data (meta_data[pos++], "# pid size");
          io.data (meta_data[pos++], "# sid size");
          io.data (meta_data[pos++], "# p-level size");
          // Boundary Condition sizes
          io.data (meta_data[pos++], "# eid size");   // elem id
          io.data (meta_data[pos++], "# side size");  // side number
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
  const T & n_elem = meta_data[0];
  const T & n_nodes = meta_data[1];

  mesh.reserve_elem(n_elem);
  mesh.reserve_nodes(n_nodes);

  // Our mesh is pre-partitioned as it's created
  this->set_n_partitions(this->n_processors());

  /**
   * We are future proofing the layout of this file by adding in size information for all stored types.
   * TODO: All types are stored as the same size. Use the size information to pack things efficiently.
   * For now we will assume that "type size" is how the entire file will be encoded.
   */
  if (version_at_least_0_9_2())
    _field_width = meta_data[2];
}



void XdrIO::read_serialized_subdomain_names(Xdr & io)
{
  const bool read_entity_info = version_at_least_0_9_2();
  const bool use_new_header_type (this->version_at_least_1_3_0());
  if (read_entity_info)
    {
      MeshBase & mesh = MeshInput<MeshBase>::mesh();

      new_header_id_type n_subdomain_names = 0;
      std::vector<new_header_id_type> subdomain_ids;
      std::vector<std::string>  subdomain_names;

      // Read the sideset names
      if (this->processor_id() == 0)
        {
          if (use_new_header_type)
            io.data(n_subdomain_names);
          else
            {
              old_header_id_type temp;
              io.data(temp);
              n_subdomain_names = temp;
            }

          subdomain_ids.resize(n_subdomain_names);
          subdomain_names.resize(n_subdomain_names);

          if (n_subdomain_names)
            {
              if (use_new_header_type)
                io.data(subdomain_ids);
              else
                {
                  std::vector<old_header_id_type> temp;
                  io.data(temp);
                  subdomain_ids.assign(temp.begin(), temp.end());
                }

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
void XdrIO::read_serialized_connectivity (Xdr & io, const dof_id_type n_elem, std::vector<new_header_id_type> & sizes, T)
{
  libmesh_assert (io.reading());

  if (!n_elem) return;

  const bool
    read_p_level      = ("." == this->polynomial_level_file_name()),
    read_partitioning = ("." == this->partition_map_file_name()),
    read_subdomain_id = ("." == this->subdomain_map_file_name());

  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Keep track of what kinds of elements this file contains
  elems_of_dimension.clear();
  elems_of_dimension.resize(4, false);

  std::vector<T> conn, input_buffer(100 /* oversized ! */);

  int level=-1;

  // Version 0.9.2+ introduces unique ids
  const size_t unique_id_size_index = 3;

  const bool read_unique_id =
    (version_at_least_0_9_2()) &&
    sizes[unique_id_size_index];

  T n_elem_at_level=0, n_processed_at_level=0;
  for (dof_id_type blk=0, first_elem=0, last_elem=0;
       last_elem<n_elem; blk++)
    {
      first_elem = cast_int<dof_id_type>(blk*io_blksize);
      last_elem  = cast_int<dof_id_type>(std::min(cast_int<std::size_t>((blk+1)*io_blksize),
                                                  cast_int<std::size_t>(n_elem)));

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
          // We are on all processors here, so we can easily assign
          // consistent unique ids if the file doesn't specify them
          // later.
          unique_id_type unique_id = e;
#endif
          if (read_unique_id)
            {
#ifdef LIBMESH_ENABLE_UNIQUE_ID
              unique_id  = cast_int<unique_id_type>(*it);
#endif
              ++it;
            }
          const dof_id_type parent_id =
            (*it == static_cast<T>(-1)) ?
            DofObject::invalid_id :
            cast_int<dof_id_type>(*it);
          ++it;
          const processor_id_type proc_id =
            cast_int<processor_id_type>(*it);
          ++it;
          const subdomain_id_type subdomain_id =
            cast_int<subdomain_id_type>(*it);
          ++it;
#ifdef LIBMESH_ENABLE_AMR
          const unsigned int p_level =
            cast_int<unsigned int>(*it);
#endif
          ++it;

          Elem * parent = (parent_id == DofObject::invalid_id) ?
            libmesh_nullptr : mesh.elem_ptr(parent_id);

          Elem * elem = Elem::build (elem_type, parent).release();

          elem->set_id() = e;
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          elem->set_unique_id() = unique_id;
#endif
          elem->processor_id() = proc_id;
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

          for (unsigned int n=0, n_n = elem->n_nodes(); n != n_n;
               n++, ++it)
            {
              const dof_id_type global_node_number =
                cast_int<dof_id_type>(*it);

              elem->set_node(n) =
                mesh.add_point (Point(), global_node_number);
            }

          elems_of_dimension[elem->dim()] = true;
          mesh.add_elem(elem);
        }
    }

  // Set the mesh dimension to the largest encountered for an element
  for (unsigned char i=0; i!=4; ++i)
    if (elems_of_dimension[i])
      mesh.set_mesh_dimension(i);

#if LIBMESH_DIM < 3
  if (mesh.mesh_dimension() > LIBMESH_DIM)
    libmesh_error_msg("Cannot open dimension "              \
                      << mesh.mesh_dimension()                          \
                      << " mesh file when configured without "          \
                      << mesh.mesh_dimension()                          \
                      << "D support.");
#endif
}



void XdrIO::read_serialized_nodes (Xdr & io, const dof_id_type n_nodes)
{
  libmesh_assert (io.reading());

  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  if (!mesh.n_nodes()) return;

  // At this point the elements have been read from file and placeholder nodes
  // have been assigned.  These nodes, however, do not have the proper (x,y,z)
  // locations or unique_id values.  This method will read all the
  // nodes from disk, and each processor can then grab the individual
  // values it needs.

  // If the file includes unique ids for nodes (as indicated by a
  // flag in 0.9.6+ files), those will be read next.

  // build up a list of the nodes contained in our local mesh.  These are the nodes
  // stored on the local processor whose (x,y,z) and unique_id values
  // need to be corrected.
  std::vector<dof_id_type> needed_nodes; needed_nodes.reserve (mesh.n_nodes());
  {
    for (auto & node : mesh.node_ptr_range())
      needed_nodes.push_back(node->id());

    std::sort (needed_nodes.begin(), needed_nodes.end());

    // We should not have any duplicate node->id()s
    libmesh_assert (std::unique(needed_nodes.begin(), needed_nodes.end()) == needed_nodes.end());
  }

  // Get the nodes in blocks.
  std::vector<Real> coords;
  std::pair<std::vector<dof_id_type>::iterator,
            std::vector<dof_id_type>::iterator> pos;
  pos.first = needed_nodes.begin();

  // Broadcast node coordinates
  for (std::size_t blk=0, first_node=0, last_node=0; last_node<n_nodes; blk++)
    {
      first_node = blk*io_blksize;
      last_node  = std::min((blk+1)*io_blksize, std::size_t(n_nodes));

      coords.resize(3*(last_node - first_node));

      if (this->processor_id() == 0)
        io.data_stream (coords.empty() ? libmesh_nullptr : &coords[0],
                        cast_int<unsigned int>(coords.size()));

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
              mesh.node_ref(cast_int<dof_id_type>(n)) =
                Point (coords[idx+0],
                       coords[idx+1],
                       coords[idx+2]);

            }
        }
    }

  if (version_at_least_0_9_6())
    {
      // Check for node unique ids
      unsigned short read_unique_ids;

      if (this->processor_id() == 0)
        io.data (read_unique_ids);

      this->comm().broadcast (read_unique_ids);

      // If no unique ids are in the file, we're done.
      if (!read_unique_ids)
        return;

      std::vector<uint32_t> unique_32;
      std::vector<uint64_t> unique_64;

      // We're starting over from node 0 again
      pos.first = needed_nodes.begin();

      for (std::size_t blk=0, first_node=0, last_node=0; last_node<n_nodes; blk++)
        {
          first_node = blk*io_blksize;
          last_node  = std::min((blk+1)*io_blksize, std::size_t(n_nodes));

          libmesh_assert((_field_width == 8) || (_field_width == 4));

          if (_field_width == 8)
            unique_64.resize(last_node - first_node);
          else
            unique_32.resize(last_node - first_node);

          if (this->processor_id() == 0)
            {
              if (_field_width == 8)
                io.data_stream (unique_64.empty() ? libmesh_nullptr : &unique_64[0],
                                cast_int<unsigned int>(unique_64.size()));
              else
                io.data_stream (unique_32.empty() ? libmesh_nullptr : &unique_32[0],
                                cast_int<unsigned int>(unique_32.size()));
            }

#ifdef LIBMESH_ENABLE_UNIQUE_ID
          if (_field_width == 8)
            this->comm().broadcast (unique_64);
          else
            this->comm().broadcast (unique_32);

          for (std::size_t n=first_node, idx=0; n<last_node; n++, idx++)
            {
              // first see if we need this node.  use pos.first as a smart lower
              // bound, this will ensure that the size of the searched range
              // decreases as we match nodes.
              pos = std::equal_range (pos.first, needed_nodes.end(), n);

              if (pos.first != pos.second) // we need this node.
                {
                  libmesh_assert_equal_to (*pos.first, n);
                  if (_field_width == 8)
                    mesh.node_ref(cast_int<dof_id_type>(n)).set_unique_id()
                      = unique_64[idx];
                  else
                    mesh.node_ref(cast_int<dof_id_type>(n)).set_unique_id()
                      = unique_32[idx];
                }
            }
#endif // LIBMESH_ENABLE_UNIQUE_ID
        }
    }
}



template <typename T>
void XdrIO::read_serialized_bcs_helper (Xdr & io, T, const std::string bc_type)
{
  if (this->boundary_condition_file_name() == "n/a") return;

  libmesh_assert (io.reading());

  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // and our boundary info object
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  // Version 0.9.2+ introduces unique ids
  read_serialized_bc_names(io, boundary_info, true);  // sideset names

  std::vector<DofBCData> dof_bc_data;
  std::vector<T> input_buffer;

  new_header_id_type n_bcs=0;
  if (this->processor_id() == 0)
    {
      if (this->version_at_least_1_3_0())
        io.data (n_bcs);
      else
        {
          old_header_id_type temp;
          io.data (temp);
          n_bcs = temp;
        }
    }
  this->comm().broadcast (n_bcs);

  for (std::size_t blk=0, first_bc=0, last_bc=0; last_bc<n_bcs; blk++)
    {
      first_bc = blk*io_blksize;
      last_bc  = std::min((blk+1)*io_blksize, std::size_t(n_bcs));

      input_buffer.resize (3*(last_bc - first_bc));

      if (this->processor_id() == 0)
        io.data_stream (input_buffer.empty() ? libmesh_nullptr : &input_buffer[0],
                        cast_int<unsigned int>(input_buffer.size()));

      this->comm().broadcast (input_buffer);
      dof_bc_data.clear();
      dof_bc_data.reserve (input_buffer.size()/3);

      // convert the input_buffer to DofBCData to facilitate searching
      for (std::size_t idx=0; idx<input_buffer.size(); idx+=3)
        dof_bc_data.push_back
          (DofBCData(cast_int<dof_id_type>(input_buffer[idx+0]),
                     cast_int<unsigned short>(input_buffer[idx+1]),
                     cast_int<boundary_id_type>(input_buffer[idx+2])));
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
      // We cannot rely on libmesh_nullptr neighbors at this point since the neighbor
      // data structure has not been initialized.
      for (; it!=end; ++it)
        {
          std::pair<std::vector<DofBCData>::iterator,
                    std::vector<DofBCData>::iterator> bounds =
            std::equal_range (dof_bc_data.begin(),
                              dof_bc_data.end(),
                              (*it)->id()
#if defined(__SUNPRO_CC) || defined(__PGI)
                              , CompareIntDofBCData()
#endif
                              );

          for (const auto & data : as_range(bounds))
            {
              libmesh_assert_equal_to (data.dof_id, (*it)->id());

              if (bc_type == "side")
                {
                  libmesh_assert_less (data.side, (*it)->n_sides());
                  boundary_info.add_side (*it, data.side, data.bc_id);
                }
              else if (bc_type == "edge")
                {
                  libmesh_assert_less (data.side, (*it)->n_edges());
                  boundary_info.add_edge (*it, data.side, data.bc_id);
                }
              else if (bc_type == "shellface")
                {
                  // Shell face IDs can only be 0 or 1.
                  libmesh_assert_less(data.side, 2);

                  boundary_info.add_shellface (*it, data.side, data.bc_id);
                }
              else
                {
                  libmesh_error_msg("bc_type not recognized: " + bc_type);
                }
            }
        }
    }
}



template <typename T>
void XdrIO::read_serialized_side_bcs (Xdr & io, T value)
{
  read_serialized_bcs_helper(io, value, "side");
}



template <typename T>
void XdrIO::read_serialized_edge_bcs (Xdr & io, T value)
{
  read_serialized_bcs_helper(io, value, "edge");
}



template <typename T>
void XdrIO::read_serialized_shellface_bcs (Xdr & io, T value)
{
  read_serialized_bcs_helper(io, value, "shellface");
}



template <typename T>
void XdrIO::read_serialized_nodesets (Xdr & io, T)
{
  if (this->boundary_condition_file_name() == "n/a") return;

  libmesh_assert (io.reading());

  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // and our boundary info object
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  // Version 0.9.2+ introduces unique ids
  read_serialized_bc_names(io, boundary_info, false); // nodeset names

  // TODO: Make a data object that works with both the element and nodal bcs
  std::vector<DofBCData> node_bc_data;
  std::vector<T> input_buffer;

  new_header_id_type n_nodesets=0;
  if (this->processor_id() == 0)
    {
      if (this->version_at_least_1_3_0())
        io.data (n_nodesets);
      else
        {
          old_header_id_type temp;
          io.data (temp);
          n_nodesets = temp;
        }
    }
  this->comm().broadcast (n_nodesets);

  for (std::size_t blk=0, first_bc=0, last_bc=0; last_bc<n_nodesets; blk++)
    {
      first_bc = blk*io_blksize;
      last_bc  = std::min((blk+1)*io_blksize, std::size_t(n_nodesets));

      input_buffer.resize (2*(last_bc - first_bc));

      if (this->processor_id() == 0)
        io.data_stream (input_buffer.empty() ? libmesh_nullptr : &input_buffer[0],
                        cast_int<unsigned int>(input_buffer.size()));

      this->comm().broadcast (input_buffer);
      node_bc_data.clear();
      node_bc_data.reserve (input_buffer.size()/2);

      // convert the input_buffer to DofBCData to facilitate searching
      for (std::size_t idx=0; idx<input_buffer.size(); idx+=2)
        node_bc_data.push_back
          (DofBCData(cast_int<dof_id_type>(input_buffer[idx+0]),
                     0,
                     cast_int<boundary_id_type>(input_buffer[idx+1])));
      input_buffer.clear();
      // note that while the files *we* write should already be sorted by
      // node id this is not necessarily guaranteed.
      std::sort (node_bc_data.begin(), node_bc_data.end());

      // Look for BCs in this block for all nodes we have
      // (not just local ones).  Do this by finding all the entries
      // in node_bc_data whose dof_id(node_id)  match the ID of the current node.
      for (auto & node : mesh.node_ptr_range())
        {
          std::pair<std::vector<DofBCData>::iterator,
                    std::vector<DofBCData>::iterator> bounds =
            std::equal_range (node_bc_data.begin(),
                              node_bc_data.end(),
                              node->id()
#if defined(__SUNPRO_CC) || defined(__PGI)
                              , CompareIntDofBCData()
#endif
                              );

          for (const auto & data : as_range(bounds))
            {
              // Note: dof_id from ElmeBCData is being used to hold node_id here
              libmesh_assert_equal_to (data.dof_id, node->id());

              boundary_info.add_node (node, data.bc_id);
            }
        }
    }
}



void XdrIO::read_serialized_bc_names(Xdr & io, BoundaryInfo & info, bool is_sideset)
{
  const bool read_entity_info = version_at_least_0_9_2();
  const bool use_new_header_type (this->version_at_least_1_3_0());
  if (read_entity_info)
    {
      new_header_id_type n_boundary_names = 0;
      std::vector<new_header_id_type> boundary_ids;
      std::vector<std::string>  boundary_names;

      // Read the sideset names
      if (this->processor_id() == 0)
        {
          if (use_new_header_type)
            io.data(n_boundary_names);
          else
            {
              old_header_id_type temp;
              io.data(temp);
              n_boundary_names = temp;
            }

          boundary_names.resize(n_boundary_names);

          if (n_boundary_names)
            {
              if (use_new_header_type)
                io.data(boundary_ids);
              else
                {
                  std::vector<old_header_id_type> temp(n_boundary_names);
                  io.data(temp);
                  boundary_ids.assign(temp.begin(), temp.end());
                }
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
        boundary_map.insert(std::make_pair(cast_int<boundary_id_type>(boundary_ids[i]), boundary_names[i]));
    }
}



void XdrIO::pack_element (std::vector<xdr_id_type> & conn, const Elem * elem,
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

  for (auto n : elem->node_index_range())
    conn.push_back (elem->node_id(n));
}

bool XdrIO::version_at_least_0_9_2() const
{
  return
    (this->version().find("0.9.2") != std::string::npos) ||
    (this->version().find("0.9.6") != std::string::npos) ||
    (this->version().find("1.1.0") != std::string::npos) ||
    (this->version().find("1.3.0") != std::string::npos);
}

bool XdrIO::version_at_least_0_9_6() const
{
  return
    (this->version().find("0.9.6") != std::string::npos) ||
    (this->version().find("1.1.0") != std::string::npos) ||
    (this->version().find("1.3.0") != std::string::npos);
}

bool XdrIO::version_at_least_1_1_0() const
{
  return
    (this->version().find("1.1.0") != std::string::npos) ||
    (this->version().find("1.3.0") != std::string::npos);
}

bool XdrIO::version_at_least_1_3_0() const
{
  return
    (this->version().find("1.3.0") != std::string::npos);
}

} // namespace libMesh
