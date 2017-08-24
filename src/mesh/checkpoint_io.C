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

#include "libmesh/checkpoint_io.h"

// C++ includes
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include <sstream> // for ostringstream

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/node.h"
#include "libmesh/parallel.h"
#include "libmesh/partitioner.h"
#include "libmesh/remote_elem.h"
#include "libmesh/xdr_io.h"
#include "libmesh/xdr_cxx.h"

namespace libMesh
{

// ------------------------------------------------------------
// CheckpointIO members
CheckpointIO::CheckpointIO (MeshBase & mesh, const bool binary_in) :
  MeshInput<MeshBase> (mesh,/* is_parallel_format = */ true),
  MeshOutput<MeshBase>(mesh,/* is_parallel_format = */ true),
  ParallelObject      (mesh),
  _binary             (binary_in),
  _parallel           (false),
  _version            ("checkpoint-1.1"),
  _my_processor_ids   (1, processor_id()),
  _my_n_processors    (n_processors())
{
}



CheckpointIO::CheckpointIO (const MeshBase & mesh, const bool binary_in) :
  MeshOutput<MeshBase>(mesh,/* is_parallel_format = */ true),
  ParallelObject      (mesh),
  _binary             (binary_in),
  _parallel           (false),
  _my_processor_ids   (1, processor_id()),
  _my_n_processors    (n_processors())
{
}



CheckpointIO::~CheckpointIO ()
{
}




void CheckpointIO::write (const std::string & name)
{
  LOG_SCOPE("write()", "CheckpointIO");

  // convenient reference to our mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // FIXME: For backwards compatibility, we'll assume for now that we
  // only want to write distributed meshes in parallel.  Later we can
  // do a gather_to_zero() and support that case too.
  _parallel = _parallel || !mesh.is_serial();

  // We'll write a header file from processor 0 to make it easier to do unambiguous
  // restarts later:
  if (this->processor_id() == 0)
    {
      Xdr io (name, this->binary() ? ENCODE : WRITE);

      // write the version
      io.data(_version, "# version");

      // Write out the max mesh dimension for backwards compatibility
      // with code that sets it independently of element dimensions
      {
        unsigned int mesh_dimension = mesh.mesh_dimension();
        io.data(mesh_dimension, "# dimensions");
      }

      // Write out whether or not this is serial output
      {
        unsigned int parallel = _parallel;
        io.data(parallel, "# parallel");
      }

      // If we're writing out a parallel mesh then we need to write the number of processors
      // so we can check it upon reading the file
      if (_parallel)
        {
          processor_id_type n_procs = _my_n_processors;
          io.data(n_procs, "# n_procs");
        }

      // write subdomain names
      this->write_subdomain_names(io);

      // write boundary id names
      const BoundaryInfo & boundary_info = mesh.get_boundary_info();
      write_bc_names(io, boundary_info, true);  // sideset names
      write_bc_names(io, boundary_info, false); // nodeset names
    }

  // If this is a serial mesh written to a serial file then we're only
  // going to write local data from processor 0.  If this is a mesh being
  // written in parallel then we're going to write from every
  // processor.
  std::vector<processor_id_type> ids_to_write;

  if (_parallel)
    {
      ids_to_write = _my_processor_ids;
    }
  else if (mesh.is_serial())
    {
      if (mesh.processor_id() == 0)
        {
          // placeholder
          ids_to_write.push_back(0);
        }
    }
  else
    {
      libmesh_error_msg("Cannot write serial checkpoint from distributed mesh");
    }

  for (std::vector<processor_id_type>::const_iterator
         id_it = ids_to_write.begin(), id_end = ids_to_write.end();
         id_it != id_end; ++id_it)
    {
      const processor_id_type my_pid = *id_it;

      std::ostringstream file_name_stream;

      file_name_stream << name << "-" << (_parallel ? _my_n_processors : 1) << "-" << my_pid;

      Xdr io (file_name_stream.str(), this->binary() ? ENCODE : WRITE);

      std::set<const Elem *, CompareElemIdsByLevel> elements;

      // For serial files or for already-distributed meshs, we write
      // everything we can see.
      if (!_parallel || !mesh.is_serial())
        elements.insert(mesh.elements_begin(), mesh.elements_end());
      // For parallel files written from serial meshes we write what
      // we'd be required to keep if we were to be deleting remote
      // elements.  This allows us to write proper parallel files even
      // from a ReplicateMesh.
      //
      // WARNING: If we have a DistributedMesh which used
      // "add_extra_ghost_elem" rather than ghosting functors to
      // preserve elements and which is *also* currently serialized
      // then we're not preserving those elements here.  As a quick
      // workaround user code should delete_remote_elements() before
      // writing the checkpoint; as a long term workaround user code
      // should use ghosting functors instead of extra_ghost_elem
      // lists.
      else
        {
          query_ghosting_functors(mesh, my_pid, false, elements);
          query_ghosting_functors(mesh, DofObject::invalid_processor_id, false, elements);
          connect_children(mesh, my_pid, elements);
          connect_children(mesh, DofObject::invalid_processor_id, elements);
          connect_families(elements);
        }

      std::set<const Node *> connected_nodes;
      reconnect_nodes(elements, connected_nodes);

      // write the nodal locations
      this->write_nodes (io, connected_nodes);

      // write connectivity
      this->write_connectivity (io, elements);

      // write remote_elem connectivity
      this->write_remote_elem (io, elements);

      // write the boundary condition information
      this->write_bcs (io, elements);

      // write the nodeset information
      this->write_nodesets (io, connected_nodes);

      // close it up
      io.close();
    }

  // this->comm().barrier();
}

void CheckpointIO::write_subdomain_names(Xdr & io) const
{
  {
    const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

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



void CheckpointIO::write_nodes (Xdr & io,
                                const std::set<const Node *> & nodeset) const
{
  dof_id_type n_nodes_here = nodeset.size();

  io.data(n_nodes_here, "# n_nodes on proc");

  // Will hold the node id and pid
  std::vector<largest_id_type> id_pid(2);

  // For the coordinates
  std::vector<Real> coords(LIBMESH_DIM);

  for (std::set<const Node *>::iterator it = nodeset.begin(),
         end = nodeset.end(); it != end; ++it)
    {
      const Node & node = **it;

      id_pid[0] = node.id();
      id_pid[1] = node.processor_id();

      io.data_stream(&id_pid[0], 2, 2);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      largest_id_type unique_id = node.unique_id();

      io.data(unique_id, "# unique id");
#endif

      coords[0] = node(0);

#if LIBMESH_DIM > 1
      coords[1] = node(1);
#endif

#if LIBMESH_DIM > 2
      coords[2] = node(2);
#endif

      io.data_stream(&coords[0], LIBMESH_DIM, 3);
    }
}



void CheckpointIO::write_connectivity (Xdr & io,
                                       const std::set<const Elem *, CompareElemIdsByLevel> & elements) const
{
  libmesh_assert (io.writing());

  // Put these out here to reduce memory churn
  // id type pid subdomain_id parent_id
  std::vector<largest_id_type> elem_data(6);
  std::vector<largest_id_type> conn_data;

  dof_id_type n_elems_here = elements.size();

  io.data(n_elems_here, "# number of elements");

  for (std::set<const Elem *, CompareElemIdsByLevel>::const_iterator it = elements.begin(),
         end = elements.end(); it != end; ++it)
    {
      const Elem & elem = **it;

      unsigned int n_nodes = elem.n_nodes();

      elem_data[0] = elem.id();
      elem_data[1] = elem.type();
      elem_data[2] = elem.processor_id();
      elem_data[3] = elem.subdomain_id();

#ifdef LIBMESH_ENABLE_AMR
      if (elem.parent() != libmesh_nullptr)
        {
          elem_data[4] = elem.parent()->id();
          elem_data[5] = elem.parent()->which_child_am_i(&elem);
        }
      else
#endif
        {
          elem_data[4] = DofObject::invalid_processor_id;
          elem_data[5] = DofObject::invalid_processor_id;
        }

      conn_data.resize(n_nodes);

      for (unsigned int i=0; i<n_nodes; i++)
        conn_data[i] = elem.node_id(i);

      io.data_stream(&elem_data[0],
                     cast_int<unsigned int>(elem_data.size()),
                     cast_int<unsigned int>(elem_data.size()));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      largest_id_type unique_id = elem.unique_id();

      io.data(unique_id, "# unique id");
#endif

#ifdef LIBMESH_ENABLE_AMR
      unsigned int p_level = elem.p_level();
      io.data(p_level, "# p_level");

      unsigned short rflag = elem.refinement_flag();
      io.data(rflag, "# rflag");

      unsigned short pflag = elem.p_refinement_flag();
      io.data(pflag, "# pflag");
#endif
      io.data_stream(&conn_data[0],
                     cast_int<unsigned int>(conn_data.size()),
                     cast_int<unsigned int>(conn_data.size()));
    }
}


void CheckpointIO::write_remote_elem (Xdr & io,
                                      const std::set<const Elem *, CompareElemIdsByLevel> & elements) const
{
  libmesh_assert (io.writing());

  // Find the remote_elem neighbor and child links
  std::vector<dof_id_type> elem_ids, parent_ids;
  std::vector<unsigned short int> elem_sides, child_numbers;

  for (std::set<const Elem *, CompareElemIdsByLevel>::const_iterator it = elements.begin(),
         end = elements.end(); it != end; ++it)
    {
      const Elem & elem = **it;

      for (unsigned int n=0; n != elem.n_neighbors(); ++n)
        {
          const Elem * neigh = elem.neighbor_ptr(n);
          if (neigh == remote_elem ||
              (neigh && !elements.count(neigh)))
            {
              elem_ids.push_back(elem.id());
              elem_sides.push_back(n);
            }
        }

#ifdef LIBMESH_ENABLE_AMR
      if (elem.has_children())
        for (unsigned int c = 0; c != elem.n_children(); ++c)
        {
          const Elem * child = elem.child_ptr(c);
          if (child == remote_elem ||
              (child && !elements.count(child)))
            {
              parent_ids.push_back(elem.id());
              child_numbers.push_back(c);
            }
        }
#endif
    }

  io.data(elem_ids, "# remote neighbor elem_ids");
  io.data(elem_sides, "# remote neighbor elem_sides");
  io.data(parent_ids, "# remote child parent_ids");
  io.data(child_numbers, "# remote child_numbers");
}



void CheckpointIO::write_bcs (Xdr & io,
                              const std::set<const Elem *, CompareElemIdsByLevel> & elements) const
{
  libmesh_assert (io.writing());

  // convenient reference to our mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // and our boundary info object
  const BoundaryInfo & boundary_info = mesh.get_boundary_info();

  std::vector<dof_id_type> full_element_id_list;
  std::vector<unsigned short int> full_side_list;
  std::vector<boundary_id_type> full_bc_id_list;

  boundary_info.build_side_list(full_element_id_list, full_side_list, full_bc_id_list);

  std::size_t bc_size = full_element_id_list.size();
  libmesh_assert_equal_to(bc_size, full_side_list.size());
  libmesh_assert_equal_to(bc_size, full_bc_id_list.size());

  std::vector<dof_id_type> element_id_list;
  std::vector<unsigned short int> side_list;
  std::vector<boundary_id_type> bc_id_list;

  element_id_list.reserve(bc_size);
  side_list.reserve(bc_size);
  bc_id_list.reserve(bc_size);

  for (std::size_t i = 0; i != bc_size; ++i)
    if (elements.count(mesh.elem_ptr(full_element_id_list[i])))
      {
        element_id_list.push_back(full_element_id_list[i]);
        side_list.push_back(full_side_list[i]);
        bc_id_list.push_back(full_bc_id_list[i]);
      }


  io.data(element_id_list, "# element ids for bcs");
  io.data(side_list, "# sides of elements for bcs");
  io.data(bc_id_list, "# bc ids");
}



void CheckpointIO::write_nodesets (Xdr & io,
                                   const std::set<const Node *> & nodeset) const
{
  libmesh_assert (io.writing());

  // convenient reference to our mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // and our boundary info object
  const BoundaryInfo & boundary_info = mesh.get_boundary_info();

  std::vector<dof_id_type> full_node_id_list;
  std::vector<boundary_id_type> full_bc_id_list;

  boundary_info.build_node_list(full_node_id_list, full_bc_id_list);

  std::size_t nodeset_size = full_node_id_list.size();
  libmesh_assert_equal_to(nodeset_size, full_bc_id_list.size());

  std::vector<dof_id_type> node_id_list;
  std::vector<boundary_id_type> bc_id_list;

  node_id_list.reserve(nodeset_size);
  bc_id_list.reserve(nodeset_size);

  for (std::size_t i = 0; i != nodeset_size; ++i)
    if (nodeset.count(mesh.node_ptr(full_node_id_list[i])))
      {
        node_id_list.push_back(full_node_id_list[i]);
        bc_id_list.push_back(full_bc_id_list[i]);
      }

  io.data(node_id_list, "# node id list");
  io.data(bc_id_list, "# nodeset bc id list");
}



void CheckpointIO::write_bc_names (Xdr & io, const BoundaryInfo & info, bool is_sideset) const
{
  const std::map<boundary_id_type, std::string> & boundary_map = is_sideset ?
    info.get_sideset_name_map() : info.get_nodeset_name_map();

  std::vector<boundary_id_type> boundary_ids;   boundary_ids.reserve(boundary_map.size());
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



void CheckpointIO::read (const std::string & name)
{
  LOG_SCOPE("read()","CheckpointIO");

  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  libmesh_assert(!mesh.n_elem());

  // Will this be a parallel input file?  With how many processors?  Stay tuned!
  unsigned int input_parallel;
  processor_id_type input_n_procs;

  // We'll write a header file from processor 0 and broadcast.
  if (this->processor_id() == 0)
    {
      {
        std::ifstream in (name.c_str());

        if (!in.good())
          libmesh_error_msg("ERROR: cannot locate header file:\n\t" << name);
      }

      Xdr io (name, this->binary() ? DECODE : READ);

      // read the version, but don't care about it
      std::string input_version;
      io.data(input_version);

      // read the dimension
      io.data (_mesh_dimension);
      this->comm().broadcast(_mesh_dimension);
      mesh.set_mesh_dimension(_mesh_dimension);

      // Read whether or not this is a parallel file
      io.data(input_parallel);
      this->comm().broadcast(input_parallel);

      // With how many processors?
      if (input_parallel)
        {
          io.data(input_n_procs);
          this->comm().broadcast(input_n_procs);
        }
      else
        input_n_procs = 1;

      // read and broadcast subdomain names
      this->read_subdomain_names(io);

      std::map<subdomain_id_type, std::string> & subdomain_map =
        mesh.set_subdomain_name_map();

      this->comm().broadcast(subdomain_map);

      // read and broadcast boundary names
      BoundaryInfo & boundary_info = mesh.get_boundary_info();

      read_bc_names(io, boundary_info, true);  // sideset names
      this->comm().broadcast(boundary_info.set_sideset_name_map());

      read_bc_names(io, boundary_info, false); // nodeset names
      this->comm().broadcast(boundary_info.set_nodeset_name_map());
    }
  // We'll receive the header broadcast everywhere else.
  else
    {
      this->comm().broadcast(_mesh_dimension);
      mesh.set_mesh_dimension(_mesh_dimension);

      this->comm().broadcast(input_parallel);
      if (input_parallel)
        this->comm().broadcast(input_n_procs);
      else
        input_n_procs = 1;

      std::map<subdomain_id_type, std::string> & subdomain_map =
        mesh.set_subdomain_name_map();
      this->comm().broadcast(subdomain_map);

      BoundaryInfo & boundary_info = mesh.get_boundary_info();
      this->comm().broadcast(boundary_info.set_sideset_name_map());
      this->comm().broadcast(boundary_info.set_nodeset_name_map());
    }


  // If this is a serial read then we're going to only read the mesh
  // on processor 0, then broadcast it
  if ((input_parallel && !mesh.is_replicated()) || mesh.processor_id() == 0)
    {
      // If we're trying to read a parallel checkpoint file on a
      // replicated mesh, we'll read every file on processor 0 so we
      // can broadcast it later.  If we're on a distributed mesh then
      // we'll read every id to it's own processor and we'll "wrap
      // around" with any ids that exceed our processor count.
      const processor_id_type begin_proc_id =
        (input_parallel && !mesh.is_replicated()) ?
        mesh.processor_id() : 0;
      const processor_id_type stride =
        (input_parallel && !mesh.is_replicated()) ?
        mesh.n_processors() : 1;

      for (processor_id_type proc_id = begin_proc_id; proc_id < input_n_procs; proc_id += stride)
        {
          std::ostringstream file_name_stream;

          file_name_stream << name;

          file_name_stream << "-" << (input_parallel ? input_n_procs : 1) << "-" << proc_id;

          {
            std::ifstream in (file_name_stream.str().c_str());

            if (!in.good())
              libmesh_error_msg("ERROR: cannot locate specified file:\n\t" << file_name_stream.str());
          }

          Xdr io (file_name_stream.str(), this->binary() ? DECODE : READ);

          // read the nodal locations
          this->read_nodes (io);

          // read connectivity
          this->read_connectivity (io);

          // Do we expect all our files' remote_elem entries to really
          // be remote?  Only if we're not reading multiple input
          // files on the same processor.
          const bool expect_all_remote =
            (input_n_procs <= mesh.n_processors() &&
             !mesh.is_replicated());

          // read remote_elem connectivity
          this->read_remote_elem (io, expect_all_remote);

          // read the boundary conditions
          this->read_bcs (io);

          // read the nodesets
          this->read_nodesets (io);

          io.close();
        }
    }

  // If the mesh was only read on processor 0 then we need to broadcast it
  if (mesh.is_replicated())
    MeshCommunication().broadcast(mesh);
  // If the mesh is really distributed then we need to make sure it
  // knows that
  else if (mesh.n_processors() > 1)
    mesh.set_distributed();
}



void CheckpointIO::read_subdomain_names(Xdr & io)
{
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  std::map<subdomain_id_type, std::string> & subdomain_map =
    mesh.set_subdomain_name_map();

  std::vector<header_id_type> subdomain_ids;
  subdomain_ids.reserve(subdomain_map.size());

  std::vector<std::string>  subdomain_names;
  subdomain_names.reserve(subdomain_map.size());

  header_id_type n_subdomain_names = 0;
  io.data(n_subdomain_names, "# subdomain id to name map");

  if (n_subdomain_names)
    {
      io.data(subdomain_ids);
      io.data(subdomain_names);

      for (std::size_t i=0; i<subdomain_ids.size(); i++)
        subdomain_map[cast_int<subdomain_id_type>(subdomain_ids[i])] =
          subdomain_names[i];
    }
}



void CheckpointIO::read_nodes (Xdr & io)
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  unsigned int n_nodes_here;
  io.data(n_nodes_here, "# n_nodes on proc");

  // Will hold the node id and pid
  std::vector<largest_id_type> id_pid(2);

  // For the coordinates
  std::vector<Real> coords(LIBMESH_DIM);

  for (unsigned int i=0; i<n_nodes_here; i++)
    {
      io.data_stream(&id_pid[0], 2, 2);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      largest_id_type unique_id = 0;
      io.data(unique_id, "# unique id");
#endif

      io.data_stream(&coords[0], LIBMESH_DIM, LIBMESH_DIM);

      Point p;
      p(0) = coords[0];

#if LIBMESH_DIM > 1
      p(1) = coords[1];
#endif

#if LIBMESH_DIM > 2
      p(2) = coords[2];
#endif

      const dof_id_type id = cast_int<dof_id_type>(id_pid[0]);

      // "Wrap around" if we see more processors than we're using.
      processor_id_type pid =
        cast_int<processor_id_type>(id_pid[1] % mesh.n_processors());

      // If we already have this node (e.g. from another file, when
      // reading multiple distributed CheckpointIO files into a
      // ReplicatedMesh) then we don't want to add it again (because
      // ReplicatedMesh can't handle that) but we do want to assert
      // consistency between what we're reading and what we have.
      const Node * old_node = mesh.query_node_ptr(id);

      if (old_node)
        {
          libmesh_assert_equal_to(pid, old_node->processor_id());
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          libmesh_assert_equal_to(unique_id, old_node->unique_id());
#endif
        }
      else
        {
#ifdef LIBMESH_ENABLE_UNIQUE_ID
          Node * node =
#endif
            mesh.add_point(p, id, pid);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
          node->set_unique_id() = unique_id;
#endif
        }
    }
}



void CheckpointIO::read_connectivity (Xdr & io)
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  dof_id_type n_elems_here;
  io.data(n_elems_here);

  // Keep track of the highest dimensional element we've added to the mesh
  unsigned int highest_elem_dim = 1;

  for (unsigned int i=0; i<n_elems_here; i++)
    {
      // id type pid subdomain_id parent_id
      std::vector<largest_id_type> elem_data(6);
      io.data_stream
        (&elem_data[0], cast_int<unsigned int>(elem_data.size()),
         cast_int<unsigned int>(elem_data.size()));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      largest_id_type unique_id = 0;
      io.data(unique_id, "# unique id");
#endif

#ifdef LIBMESH_ENABLE_AMR
      unsigned int p_level = 0;
      io.data(p_level, "# p_level");

      unsigned short rflag, pflag;
      io.data(rflag, "# rflag");
      io.data(pflag, "# pflag");
#endif

      unsigned int n_nodes = Elem::type_to_n_nodes_map[elem_data[1]];

      // Snag the node ids this element was connected to
      std::vector<largest_id_type> conn_data(n_nodes);
      io.data_stream
        (&conn_data[0], cast_int<unsigned int>(conn_data.size()),
         cast_int<unsigned int>(conn_data.size()));

      const dof_id_type id                 =
        cast_int<dof_id_type>      (elem_data[0]);
      const ElemType elem_type             =
        static_cast<ElemType>      (elem_data[1]);
      const processor_id_type proc_id      =
        cast_int<processor_id_type>
        (elem_data[2] % mesh.n_processors());
      const subdomain_id_type subdomain_id =
        cast_int<subdomain_id_type>(elem_data[3]);
      const dof_id_type parent_id          =
        cast_int<dof_id_type>      (elem_data[4]);
      const unsigned short int child_num   =
        cast_int<dof_id_type>      (elem_data[5]);

      Elem * parent =
        (parent_id == DofObject::invalid_processor_id) ?
        libmesh_nullptr : mesh.elem_ptr(parent_id);

      if (!parent)
        libmesh_assert_equal_to
          (child_num, DofObject::invalid_processor_id);

      Elem * old_elem = mesh.query_elem_ptr(id);

      // If we already have this element (e.g. from another file,
      // when reading multiple distributed CheckpointIO files into
      // a ReplicatedMesh) then we don't want to add it again
      // (because ReplicatedMesh can't handle that) but we do want
      // to assert consistency between what we're reading and what
      // we have.
      if (old_elem)
        {
          libmesh_assert_equal_to(elem_type, old_elem->type());
          libmesh_assert_equal_to(proc_id, old_elem->processor_id());
          libmesh_assert_equal_to(subdomain_id, old_elem->subdomain_id());
          if (parent)
            libmesh_assert_equal_to(parent, old_elem->parent());
          else
            libmesh_assert(!old_elem->parent());

          libmesh_assert_equal_to(old_elem->n_nodes(), conn_data.size());

          for (std::size_t n=0; n != conn_data.size(); ++n)
            libmesh_assert_equal_to
              (old_elem->node_id(n),
               cast_int<dof_id_type>(conn_data[n]));
        }
      else
        {
          // Create the element
          Elem * elem = Elem::build(elem_type, parent).release();

#ifdef LIBMESH_ENABLE_UNIQUE_ID
          elem->set_unique_id() = unique_id;
#endif

          if (elem->dim() > highest_elem_dim)
            highest_elem_dim = elem->dim();

          elem->set_id()       = id;
          elem->processor_id() = proc_id;
          elem->subdomain_id() = subdomain_id;

#ifdef LIBMESH_ENABLE_AMR
          elem->hack_p_level(p_level);

          elem->set_refinement_flag  (cast_int<Elem::RefinementState>(rflag));
          elem->set_p_refinement_flag(cast_int<Elem::RefinementState>(pflag));

          // Set parent connections
          if (parent)
            {
              // We must specify a child_num, because we will have
              // skipped adding any preceding remote_elem children
              parent->add_child(elem, child_num);
            }
#endif

          libmesh_assert(elem->n_nodes() == conn_data.size());

          // Connect all the nodes to this element
          for (std::size_t n=0; n<conn_data.size(); n++)
            elem->set_node(n) =
              mesh.node_ptr(cast_int<dof_id_type>(conn_data[n]));

          mesh.add_elem(elem);
        }
    }

  mesh.set_mesh_dimension(cast_int<unsigned char>(highest_elem_dim));
}


void CheckpointIO::read_remote_elem (Xdr & io, bool libmesh_dbg_var(expect_all_remote))
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Find the remote_elem neighbor links
  std::vector<dof_id_type> elem_ids;
  std::vector<unsigned short int> elem_sides;

  io.data(elem_ids, "# remote neighbor elem_ids");
  io.data(elem_sides, "# remote neighbor elem_sides");

  libmesh_assert_equal_to(elem_ids.size(), elem_sides.size());

  for (std::size_t i=0; i != elem_ids.size(); ++i)
    {
      Elem & elem = mesh.elem_ref(elem_ids[i]);
      if (!elem.neighbor(elem_sides[i]))
        elem.set_neighbor(elem_sides[i],
                          const_cast<RemoteElem *>(remote_elem));
      else
        libmesh_assert(!expect_all_remote);
    }

  // Find the remote_elem children links
  std::vector<dof_id_type> parent_ids;
  std::vector<unsigned short int> child_numbers;

  io.data(parent_ids, "# remote child parent_ids");
  io.data(child_numbers, "# remote child_numbers");

#ifdef LIBMESH_ENABLE_AMR
  for (std::size_t i=0; i != parent_ids.size(); ++i)
    {
      Elem & elem = mesh.elem_ref(parent_ids[i]);

      // We'd like to assert that no child pointer already exists to
      // be overwritten by remote_elem, but Elem doesn't actually have
      // an API that will return a child pointer without asserting
      // that it isn't NULL
      //
      const Elem * child = elem.raw_child_ptr(child_numbers[i]);

      if (!child)
        elem.add_child(const_cast<RemoteElem *>(remote_elem),
                       child_numbers[i]);
      else
        libmesh_assert(!expect_all_remote);
    }
#endif
}



void CheckpointIO::read_bcs (Xdr & io)
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // and our boundary info object
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  std::vector<dof_id_type> element_id_list;
  std::vector<unsigned short int> side_list;
  std::vector<boundary_id_type> bc_id_list;

  io.data(element_id_list, "# element ids for bcs");
  io.data(side_list, "# sides of elements for bcs");
  io.data(bc_id_list, "# bc ids");

  for (std::size_t i=0; i<element_id_list.size(); i++)
    boundary_info.add_side(element_id_list[i], side_list[i], bc_id_list[i]);
}



void CheckpointIO::read_nodesets (Xdr & io)
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // and our boundary info object
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  std::vector<dof_id_type> node_id_list;
  std::vector<boundary_id_type> bc_id_list;

  io.data(node_id_list, "# node id list");
  io.data(bc_id_list, "# nodeset bc id list");

  for (std::size_t i=0; i<node_id_list.size(); i++)
    boundary_info.add_node(node_id_list[i], bc_id_list[i]);
}



void CheckpointIO::read_bc_names(Xdr & io, BoundaryInfo & info, bool is_sideset)
{
  std::map<boundary_id_type, std::string> & boundary_map = is_sideset ?
    info.set_sideset_name_map() : info.set_nodeset_name_map();

  std::vector<boundary_id_type> boundary_ids;
  std::vector<std::string>  boundary_names;

  header_id_type n_boundary_names = 0;

  if (is_sideset)
    io.data(n_boundary_names, "# sideset id to name map");
  else
    io.data(n_boundary_names, "# nodeset id to name map");

  if (n_boundary_names)
    {
      io.data(boundary_ids);
      io.data(boundary_names);
    }

  // Add them back into the map
  for (std::size_t i=0; i<boundary_ids.size(); i++)
    boundary_map[boundary_ids[i]] = boundary_names[i];
}


unsigned int CheckpointIO::n_active_levels_in(MeshBase::const_element_iterator begin,
                                              MeshBase::const_element_iterator end) const
{
  unsigned int max_level = 0;

  for (MeshBase::const_element_iterator it = begin;
         it != end; ++it)
    max_level = std::max((*it)->level(), max_level);

  return max_level + 1;
}

} // namespace libMesh
