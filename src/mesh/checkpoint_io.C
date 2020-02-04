// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/metis_partitioner.h"
#include "libmesh/remote_elem.h"
#include "libmesh/xdr_io.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/utility.h"
#include "libmesh/int_range.h"

// C++ includes
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <unistd.h>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream> // for ostringstream
#include <unordered_map>
#include <unordered_set>

namespace
{
// chunking computes the number of chunks and first-chunk-offset when splitting a mesh
// into nsplits pieces using size procs for the given MPI rank.  The number of chunks and offset
// are stored in nchunks and first_chunk respectively.
void chunking(libMesh::processor_id_type size, libMesh::processor_id_type rank, libMesh::processor_id_type nsplits,
              libMesh::processor_id_type & nchunks, libMesh::processor_id_type & first_chunk)
{
  if (nsplits % size == 0) // the chunks divide evenly over the processors
    {
      nchunks = nsplits / size;
      first_chunk = libMesh::cast_int<libMesh::processor_id_type>(nchunks * rank);
      return;
    }

  libMesh::processor_id_type nextra = nsplits % size;
  if (rank < nextra) // leftover chunks cause an extra chunk to be added to this processor
    {
      nchunks = libMesh::cast_int<libMesh::processor_id_type>(nsplits / size + 1);
      first_chunk = libMesh::cast_int<libMesh::processor_id_type>(nchunks * rank);
    }
  else // no extra chunks, but first chunk is offset by extras on earlier ranks
    {
      nchunks = nsplits / size;
      // account for the case where nchunks is zero where we want max int
      first_chunk = libMesh::cast_int<libMesh::processor_id_type>
        (std::max((int)((nchunks + 1) * (nsplits % size) + nchunks * (rank - nsplits % size)),
         (1 - (int)nchunks) * std::numeric_limits<int>::max()));
    }
}

std::string extension(const std::string & s)
{
  auto pos = s.rfind(".");
  if (pos == std::string::npos)
    return "";
  return s.substr(pos, s.size() - pos);
}

std::string split_dir(const std::string & input_name, libMesh::processor_id_type n_procs)
{
  return input_name + "/" + std::to_string(n_procs);
}


std::string header_file(const std::string & input_name, libMesh::processor_id_type n_procs)
{
  return split_dir(input_name, n_procs) + "/header" + extension(input_name);
}

std::string
split_file(const std::string & input_name,
                        libMesh::processor_id_type n_procs,
                        libMesh::processor_id_type proc_id)
{
  return split_dir(input_name, n_procs) + "/split-" + std::to_string(n_procs) + "-" +
         std::to_string(proc_id) + extension(input_name);
}

void make_dir(const std::string & input_name, libMesh::processor_id_type n_procs)
{
  auto ret = libMesh::Utility::mkdir(input_name.c_str());
  // error only if we failed to create dir - don't care if it was already there
  if (ret != 0 && ret != -1)
    libmesh_error_msg(
        "Failed to create mesh split directory '" << input_name << "': " << std::strerror(ret));

  auto dir_name = split_dir(input_name, n_procs);
  ret = libMesh::Utility::mkdir(dir_name.c_str());
  if (ret == -1)
    libmesh_warning("In CheckpointIO::write, directory '"
                    << dir_name << "' already exists, overwriting contents.");
  else if (ret != 0)
    libmesh_error_msg(
        "Failed to create mesh split directory '" << dir_name << "': " << std::strerror(ret));
}

} // namespace

namespace libMesh
{

std::unique_ptr<CheckpointIO> split_mesh(MeshBase & mesh, processor_id_type nsplits)
{
  // There is currently an issue with DofObjects not being properly
  // reset if the mesh is not first repartitioned onto 1 processor
  // *before* being repartitioned onto the desired number of
  // processors. So, this is a workaround, but not a particularly
  // onerous one.
  mesh.partition(1);
  mesh.partition(nsplits);

  processor_id_type my_num_chunks = 0;
  processor_id_type my_first_chunk = 0;
  chunking(mesh.comm().size(), mesh.comm().rank(), nsplits, my_num_chunks, my_first_chunk);

  auto cpr = libmesh_make_unique<CheckpointIO>(mesh);
  cpr->current_processor_ids().clear();
  for (processor_id_type i = my_first_chunk; i < my_first_chunk + my_num_chunks; i++)
    cpr->current_processor_ids().push_back(i);
  cpr->current_n_processors() = nsplits;
  cpr->parallel() = true;
  return cpr;
}


// ------------------------------------------------------------
// CheckpointIO members
CheckpointIO::CheckpointIO (MeshBase & mesh, const bool binary_in) :
  MeshInput<MeshBase> (mesh,/* is_parallel_format = */ true),
  MeshOutput<MeshBase>(mesh,/* is_parallel_format = */ true),
  ParallelObject      (mesh),
  _binary             (binary_in),
  _parallel           (false),
  _version            ("checkpoint-1.5"),
  _my_processor_ids   (1, processor_id()),
  _my_n_processors    (mesh.is_replicated() ? 1 : n_processors())
{
}

CheckpointIO::CheckpointIO (const MeshBase & mesh, const bool binary_in) :
  MeshOutput<MeshBase>(mesh,/* is_parallel_format = */ true),
  ParallelObject      (mesh),
  _binary             (binary_in),
  _parallel           (false),
  _my_processor_ids   (1, processor_id()),
  _my_n_processors    (mesh.is_replicated() ? 1 : n_processors())
{
}

CheckpointIO::~CheckpointIO ()
{
}

processor_id_type CheckpointIO::select_split_config(const std::string & input_name, header_id_type & data_size)
{
  std::string header_name;

  // We'll read a header file from processor 0 and broadcast.
  if (this->processor_id() == 0)
    {
      header_name = header_file(input_name, _my_n_processors);

      {
        // look for header+splits with nprocs equal to _my_n_processors
        std::ifstream in (header_name.c_str());
        if (!in.good())
          {
            // otherwise fall back to a serial/single-split mesh
            auto orig_header_name = header_name;
            header_name = header_file(input_name, 1);
            std::ifstream in2 (header_name.c_str());
            if (!in2.good())
              {
                libmesh_error_msg("ERROR: Neither one of the following files can be located:\n\t'"
                                  << orig_header_name << "' nor\n\t'" << input_name << "'\n"
                                  << "If you are running a parallel job, double check that you've "
                                  << "created a split for " << _my_n_processors << " ranks.\n"
                                  << "Note: One of paths above may refer to a valid directory on your "
                                  << "system, however we are attempting to read a valid header file.");
              }
          }
      }

      Xdr io (header_name, this->binary() ? DECODE : READ);

      // read the version, but don't care about it
      std::string input_version;
      io.data(input_version);

      // read the data type
      io.data (data_size);
    }

  this->comm().broadcast(data_size);
  this->comm().broadcast(header_name);

  // How many per-processor files are here?
  largest_id_type input_n_procs;

  switch (data_size) {
  case 2:
    input_n_procs = this->read_header<uint16_t>(header_name);
    break;
  case 4:
    input_n_procs = this->read_header<uint32_t>(header_name);
    break;
  case 8:
    input_n_procs = this->read_header<uint64_t>(header_name);
    break;
  default:
    libmesh_error();
  }

  if (!input_n_procs)
    input_n_procs = 1;
  return cast_int<processor_id_type>(input_n_procs);
}

void CheckpointIO::cleanup(const std::string & input_name, processor_id_type n_procs)
{
  auto header = header_file(input_name, n_procs);
  auto ret = std::remove(header.c_str());
  if (ret != 0)
    libmesh_warning("Failed to clean up checkpoint header '" << header << "': " << std::strerror(ret));

  for (processor_id_type i = 0; i < n_procs; i++)
    {
      auto split = split_file(input_name, n_procs, i);
      ret = std::remove(split.c_str());
      if (ret != 0)
        libmesh_warning("Failed to clean up checkpoint split file '" << split << "': " << std::strerror(ret));
    }

  auto dir = split_dir(input_name, n_procs);
  ret = rmdir(dir.c_str());
  if (ret != 0)
    libmesh_warning("Failed to clean up checkpoint split dir '" << dir << "': " << std::strerror(ret));

  // We expect that this may fail if there are other split configurations still present in this
  // directory - so don't bother to check/warn for failure.
  rmdir(input_name.c_str());
}


bool CheckpointIO::version_at_least_1_5() const
{
  return (this->version().find("1.5") != std::string::npos);
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

  processor_id_type use_n_procs = 1;
  if (_parallel)
    use_n_procs = _my_n_processors;

  std::string header_file_name = header_file(name, use_n_procs);
  make_dir(name, use_n_procs);

  // We'll write a header file from processor 0 to make it easier to do unambiguous
  // restarts later:
  if (this->processor_id() == 0)
    {
      Xdr io (header_file_name, this->binary() ? ENCODE : WRITE);

      // write the version
      io.data(_version, "# version");

      // write what kind of data type we're using
      header_id_type data_size = sizeof(largest_id_type);
      io.data(data_size, "# integer size");

      // Write out the max mesh dimension for backwards compatibility
      // with code that sets it independently of element dimensions
      {
        uint16_t mesh_dimension = cast_int<uint16_t>(mesh.mesh_dimension());
        io.data(mesh_dimension, "# dimensions");
      }

      // Write out whether or not this is serial output
      {
        uint16_t parallel = _parallel;
        io.data(parallel, "# parallel");
      }

      // If we're writing out a parallel mesh then we need to write the number of processors
      // so we can check it upon reading the file
      if (_parallel)
        {
          largest_id_type n_procs = _my_n_processors;
          io.data(n_procs, "# n_procs");
        }

      // write subdomain names
      this->write_subdomain_names(io);

      // write boundary id names
      const BoundaryInfo & boundary_info = mesh.get_boundary_info();
      write_bc_names(io, boundary_info, true);  // sideset names
      write_bc_names(io, boundary_info, false); // nodeset names

      // write extra integer names
      const bool write_extra_integers = this->version_at_least_1_5();

      if (write_extra_integers)
        {
          largest_id_type n_node_integers = mesh.n_node_integers();
          io.data(n_node_integers, "# n_extra_integers per node");

          std::vector<std::string> node_integer_names;
          for (unsigned int i=0; i != n_node_integers; ++i)
            node_integer_names.push_back(mesh.get_node_integer_name(i));

          io.data(node_integer_names);

          largest_id_type n_elem_integers = mesh.n_elem_integers();
          io.data(n_elem_integers, "# n_extra_integers per elem");

          std::vector<std::string> elem_integer_names;
          for (unsigned int i=0; i != n_elem_integers; ++i)
            elem_integer_names.push_back(mesh.get_elem_integer_name(i));

          io.data(elem_integer_names);
        }


    }

  // If this is a serial mesh written to a serial file then we're only
  // going to write local data from processor 0.  If this is a mesh being
  // written in parallel then we're going to write from every
  // processor.
  std::vector<processor_id_type> ids_to_write;

  // We're going to sort elements by pid in one pass, to avoid sending
  // predicated iterators through the whole mesh N_p times
  std::unordered_map<processor_id_type, std::vector<Elem *>> elements_on_pid;

  if (_parallel)
    {
      ids_to_write = _my_processor_ids;
      for (processor_id_type p : ids_to_write)
        elements_on_pid[p].clear();
      auto eop_end = elements_on_pid.end();
      for (auto & elem : mesh.element_ptr_range())
        {
          const processor_id_type p = elem->processor_id();
          auto eop_it = elements_on_pid.find(p);
          if (eop_it != eop_end)
            eop_it->second.push_back(elem);
        }
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

  // Call build_side_list() and build_node_list() just *once* to avoid
  // redundant expensive sorts during mesh splitting.
  const BoundaryInfo & boundary_info = mesh.get_boundary_info();
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
    bc_triples = boundary_info.build_side_list();
  std::vector<std::tuple<dof_id_type, boundary_id_type>>
    bc_tuples = boundary_info.build_node_list();

  for (const auto & my_pid : ids_to_write)
    {
      auto file_name = split_file(name, use_n_procs, my_pid);
      Xdr io (file_name, this->binary() ? ENCODE : WRITE);

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
          for (processor_id_type p : {my_pid, DofObject::invalid_processor_id})
            {
              const auto elements_vec_it = elements_on_pid.find(p);
              if (elements_vec_it != elements_on_pid.end())
                {
                  const auto & p_elements = elements_vec_it->second;
                  Elem * const * elempp = p_elements.data();
                  Elem * const * elemend = elempp + p_elements.size();

                  const MeshBase::const_element_iterator
                    pid_elements_begin = MeshBase::const_element_iterator
                      (elempp, elemend, Predicates::NotNull<Elem * const *>()),
                    pid_elements_end = MeshBase::const_element_iterator
                      (elemend, elemend, Predicates::NotNull<Elem * const *>()),
                    active_pid_elements_begin = MeshBase::const_element_iterator
                      (elempp, elemend, Predicates::Active<Elem * const *>()),
                    active_pid_elements_end = MeshBase::const_element_iterator
                      (elemend, elemend, Predicates::Active<Elem * const *>());

                  query_ghosting_functors
                    (mesh, p, active_pid_elements_begin,
                     active_pid_elements_end, elements);
                  connect_children(mesh, pid_elements_begin,
                                   pid_elements_end, elements);
                }
              connect_families(elements);
            }
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
      this->write_bcs (io, elements, bc_triples);

      // write the nodeset information
      this->write_nodesets (io, connected_nodes, bc_tuples);

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

    std::vector<largest_id_type> subdomain_ids;   subdomain_ids.reserve(subdomain_map.size());
    std::vector<std::string>  subdomain_names; subdomain_names.reserve(subdomain_map.size());

    // We need to loop over the map and make sure that there aren't any invalid entries.  Since we
    // return writable references in mesh_base, it's possible for the user to leave some entity names
    // blank.  We can't write those to the XDA file.
    largest_id_type n_subdomain_names = 0;
    for (const auto & pr : subdomain_map)
      if (!pr.second.empty())
        {
          n_subdomain_names++;
          subdomain_ids.push_back(pr.first);
          subdomain_names.push_back(pr.second);
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
  largest_id_type n_nodes_here = nodeset.size();

  io.data(n_nodes_here, "# n_nodes on proc");

  const bool write_extra_integers = this->version_at_least_1_5();
  const unsigned int n_extra_integers =
    write_extra_integers ? MeshOutput<MeshBase>::mesh().n_node_integers() : 0;

  // Will hold the node id and pid and extra integers
  std::vector<largest_id_type> id_pid(2 + n_extra_integers);

  // For the coordinates
  std::vector<Real> coords(LIBMESH_DIM);

  for (const auto & node : nodeset)
    {
      id_pid[0] = node->id();
      id_pid[1] = node->processor_id();

      libmesh_assert_equal_to(n_extra_integers, node->n_extra_integers());
      for (unsigned int i=0; i != n_extra_integers; ++i)
        id_pid[2+i] = node->get_extra_integer(i);

      io.data_stream(id_pid.data(), 2 + n_extra_integers, 2 + n_extra_integers);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      largest_id_type unique_id = node->unique_id();

      io.data(unique_id, "# unique id");
#endif

      coords[0] = (*node)(0);

#if LIBMESH_DIM > 1
      coords[1] = (*node)(1);
#endif

#if LIBMESH_DIM > 2
      coords[2] = (*node)(2);
#endif

      io.data_stream(coords.data(), LIBMESH_DIM, 3);
    }
}



void CheckpointIO::write_connectivity (Xdr & io,
                                       const std::set<const Elem *, CompareElemIdsByLevel> & elements) const
{
  libmesh_assert (io.writing());

  const bool write_extra_integers = this->version_at_least_1_5();
  const unsigned int n_extra_integers =
    write_extra_integers ? MeshOutput<MeshBase>::mesh().n_elem_integers() : 0;

  // Put these out here to reduce memory churn
  // id type pid subdomain_id parent_id extra_integer_0 ...
  std::vector<largest_id_type> elem_data(6 + n_extra_integers);
  std::vector<largest_id_type> conn_data;

  largest_id_type n_elems_here = elements.size();

  io.data(n_elems_here, "# number of elements");

  for (const auto & elem : elements)
    {
      unsigned int n_nodes = elem->n_nodes();

      elem_data[0] = elem->id();
      elem_data[1] = elem->type();
      elem_data[2] = elem->processor_id();
      elem_data[3] = elem->subdomain_id();

#ifdef LIBMESH_ENABLE_AMR
      if (elem->parent() != nullptr)
        {
          elem_data[4] = elem->parent()->id();
          elem_data[5] = elem->parent()->which_child_am_i(elem);
        }
      else
#endif
        {
          elem_data[4] = static_cast<largest_id_type>(-1);
          elem_data[5] = static_cast<largest_id_type>(-1);
        }

      for (unsigned int i=0; i != n_extra_integers; ++i)
        elem_data[6+i] = elem->get_extra_integer(i);

      conn_data.resize(n_nodes);

      for (unsigned int i=0; i<n_nodes; i++)
        conn_data[i] = elem->node_id(i);

      io.data_stream(elem_data.data(),
                     cast_int<unsigned int>(elem_data.size()),
                     cast_int<unsigned int>(elem_data.size()));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      largest_id_type unique_id = elem->unique_id();

      io.data(unique_id, "# unique id");
#endif

#ifdef LIBMESH_ENABLE_AMR
      uint16_t p_level = cast_int<uint16_t>(elem->p_level());
      io.data(p_level, "# p_level");

      uint16_t rflag = elem->refinement_flag();
      io.data(rflag, "# rflag");

      uint16_t pflag = elem->p_refinement_flag();
      io.data(pflag, "# pflag");
#endif
      io.data_stream(conn_data.data(),
                     cast_int<unsigned int>(conn_data.size()),
                     cast_int<unsigned int>(conn_data.size()));
    }
}


void CheckpointIO::write_remote_elem (Xdr & io,
                                      const std::set<const Elem *, CompareElemIdsByLevel> & elements) const
{
  libmesh_assert (io.writing());

  // Find the remote_elem neighbor and child links
  std::vector<largest_id_type> elem_ids, parent_ids;
  std::vector<uint16_t> elem_sides, child_numbers;

  for (const auto & elem : elements)
    {
      for (auto n : elem->side_index_range())
        {
          const Elem * neigh = elem->neighbor_ptr(n);
          if (neigh == RemoteElem::get_instance() ||
              (neigh && !elements.count(neigh)))
            {
              elem_ids.push_back(elem->id());
              elem_sides.push_back(n);
            }
        }

#ifdef LIBMESH_ENABLE_AMR
      if (elem->has_children())
        {
          for (unsigned short c = 0,
               nc = cast_int<unsigned short>(elem->n_children());
               c != nc; ++c)
            {
              const Elem * child = elem->child_ptr(c);
              if (child == RemoteElem::get_instance() ||
                  (child && !elements.count(child)))
                {
                  parent_ids.push_back(elem->id());
                  child_numbers.push_back(c);
                }
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
                              const std::set<const Elem *, CompareElemIdsByLevel> & elements,
                              const std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> & bc_triples) const
{
  libmesh_assert (io.writing());

  // Build a list of (elem, side, bc) tuples.
  std::size_t bc_size = bc_triples.size();

  std::vector<largest_id_type> element_id_list;
  std::vector<uint16_t> side_list;
  std::vector<largest_id_type> bc_id_list;

  element_id_list.reserve(bc_size);
  side_list.reserve(bc_size);
  bc_id_list.reserve(bc_size);

  std::unordered_set<dof_id_type> elems;
  for (auto & e : elements)
    elems.insert(e->id());

  for (const auto & t : bc_triples)
    if (elems.count(std::get<0>(t)))
      {
        element_id_list.push_back(std::get<0>(t));
        side_list.push_back(std::get<1>(t));
        bc_id_list.push_back(std::get<2>(t));
      }


  io.data(element_id_list, "# element ids for bcs");
  io.data(side_list, "# sides of elements for bcs");
  io.data(bc_id_list, "# bc ids");
}



void CheckpointIO::write_nodesets (Xdr & io,
                                   const std::set<const Node *> & nodeset,
                                   const std::vector<std::tuple<dof_id_type, boundary_id_type>> & bc_tuples) const
{
  libmesh_assert (io.writing());

  // convenient reference to our mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // Build a list of (node, bc) tuples
  std::size_t nodeset_size = bc_tuples.size();

  std::vector<largest_id_type> node_id_list;
  std::vector<largest_id_type> bc_id_list;

  node_id_list.reserve(nodeset_size);
  bc_id_list.reserve(nodeset_size);

  for (const auto & t : bc_tuples)
    if (nodeset.count(mesh.node_ptr(std::get<0>(t))))
      {
        node_id_list.push_back(std::get<0>(t));
        bc_id_list.push_back(std::get<1>(t));
      }

  io.data(node_id_list, "# node id list");
  io.data(bc_id_list, "# nodeset bc id list");
}



void CheckpointIO::write_bc_names (Xdr & io, const BoundaryInfo & info, bool is_sideset) const
{
  const std::map<boundary_id_type, std::string> & boundary_map = is_sideset ?
    info.get_sideset_name_map() : info.get_nodeset_name_map();

  std::vector<largest_id_type> boundary_ids;   boundary_ids.reserve(boundary_map.size());
  std::vector<std::string>  boundary_names; boundary_names.reserve(boundary_map.size());

  // We need to loop over the map and make sure that there aren't any invalid entries.  Since we
  // return writable references in boundary_info, it's possible for the user to leave some entity names
  // blank.  We can't write those to the XDA file.
  largest_id_type n_boundary_names = 0;
  for (const auto & pr : boundary_map)
    if (!pr.second.empty())
      {
        n_boundary_names++;
        boundary_ids.push_back(pr.first);
        boundary_names.push_back(pr.second);
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

void CheckpointIO::read (const std::string & input_name)
{
  LOG_SCOPE("read()","CheckpointIO");

  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  libmesh_assert(!mesh.n_elem());

  header_id_type data_size;
  processor_id_type input_n_procs = select_split_config(input_name, data_size);
  auto header_name = header_file(input_name, input_n_procs);
  bool input_parallel = input_n_procs > 0;

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

      for (processor_id_type proc_id = begin_proc_id; proc_id < input_n_procs;
           proc_id = cast_int<processor_id_type>(proc_id + stride))
        {
          auto file_name = split_file(input_name, input_n_procs, proc_id);

          {
            std::ifstream in (file_name.c_str());

            if (!in.good())
              libmesh_error_msg("ERROR: cannot locate specified file:\n\t" << file_name);
          }

          // Do we expect all our files' remote_elem entries to really
          // be remote?  Only if we're not reading multiple input
          // files on the same processor.
          const bool expect_all_remote =
            (input_n_procs <= mesh.n_processors() &&
             !mesh.is_replicated());

          Xdr io (file_name, this->binary() ? DECODE : READ);

          switch (data_size) {
          case 2:
            this->read_subfile<uint16_t>(io, expect_all_remote);
            break;
          case 4:
            this->read_subfile<uint32_t>(io, expect_all_remote);
            break;
          case 8:
            this->read_subfile<uint64_t>(io, expect_all_remote);
            break;
          default:
            libmesh_error();
          }

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



template <typename file_id_type>
file_id_type CheckpointIO::read_header (const std::string & name)
{
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Hack for codes which don't look at all elem dimensions
  uint16_t mesh_dimension;

  // Will this be a parallel input file?  With how many processors?  Stay tuned!
  uint16_t input_parallel;
  file_id_type input_n_procs;

  std::vector<std::string> node_integer_names, elem_integer_names;

  // We'll write a header file from processor 0 and broadcast.
  if (this->processor_id() == 0)
    {
      Xdr io (name, this->binary() ? DECODE : READ);

      // read the version, but don't care about it
      std::string input_version;
      io.data(input_version);

      // read the data type, don't care about it this time
      header_id_type data_size;
      io.data (data_size);

      // read the dimension
      io.data (mesh_dimension);

      // Read whether or not this is a parallel file
      io.data(input_parallel);

      // With how many processors?
      if (input_parallel)
        io.data(input_n_procs);

      // read subdomain names
      this->read_subdomain_names<file_id_type>(io);

      // read boundary names
      BoundaryInfo & boundary_info = mesh.get_boundary_info();

      this->read_bc_names<file_id_type>(io, boundary_info, true);  // sideset names
      this->read_bc_names<file_id_type>(io, boundary_info, false); // nodeset names

      // read extra integer names?
      std::swap(input_version, this->version());
      const bool read_extra_integers = this->version_at_least_1_5();
      std::swap(input_version, this->version());

      if (read_extra_integers)
        this->read_integers_names<file_id_type>
          (io, node_integer_names, elem_integer_names);
    }

  // broadcast data from processor 0, set values everywhere
  this->comm().broadcast(mesh_dimension);
  mesh.set_mesh_dimension(cast_int<unsigned char>(mesh_dimension));

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

  this->comm().broadcast(node_integer_names);
  this->comm().broadcast(elem_integer_names);

  for (auto & int_name : node_integer_names)
    mesh.add_node_integer(int_name);

  for (auto & int_name : elem_integer_names)
    mesh.add_elem_integer(int_name);

  return input_parallel ? input_n_procs : 0;
}



template <typename file_id_type>
void CheckpointIO::read_subfile (Xdr & io, bool expect_all_remote)
{
  // read the nodal locations
  this->read_nodes<file_id_type> (io);

  // read connectivity
  this->read_connectivity<file_id_type> (io);

  // read remote_elem connectivity
  this->read_remote_elem<file_id_type> (io, expect_all_remote);

  // read the boundary conditions
  this->read_bcs<file_id_type> (io);

  // read the nodesets
  this->read_nodesets<file_id_type> (io);
}



template <typename file_id_type>
void CheckpointIO::read_subdomain_names(Xdr & io)
{
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  std::map<subdomain_id_type, std::string> & subdomain_map =
    mesh.set_subdomain_name_map();

  std::vector<file_id_type> subdomain_ids;
  subdomain_ids.reserve(subdomain_map.size());

  std::vector<std::string>  subdomain_names;
  subdomain_names.reserve(subdomain_map.size());

  file_id_type n_subdomain_names = 0;
  io.data(n_subdomain_names, "# subdomain id to name map");

  if (n_subdomain_names)
    {
      io.data(subdomain_ids);
      io.data(subdomain_names);

      for (auto i : index_range(subdomain_ids))
        subdomain_map[cast_int<subdomain_id_type>(subdomain_ids[i])] =
          subdomain_names[i];
    }
}



template <typename file_id_type>
void CheckpointIO::read_nodes (Xdr & io)
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  file_id_type n_nodes_here;
  io.data(n_nodes_here, "# n_nodes on proc");

  const bool read_extra_integers = this->version_at_least_1_5();

  const unsigned int n_extra_integers =
    read_extra_integers ? mesh.n_node_integers() : 0;

  // Will hold the node id and pid and extra integers
  std::vector<file_id_type> id_pid(2 + n_extra_integers);

  // For the coordinates
  std::vector<Real> coords(LIBMESH_DIM);

  for (unsigned int i=0; i<n_nodes_here; i++)
    {
      io.data_stream(id_pid.data(), 2 + n_extra_integers, 2 + n_extra_integers);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      file_id_type unique_id = 0;
      io.data(unique_id, "# unique id");
#endif

      io.data_stream(coords.data(), LIBMESH_DIM, LIBMESH_DIM);

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

          libmesh_assert_equal_to(n_extra_integers, old_node->n_extra_integers());
#ifndef NDEBUG
          for (unsigned int ei=0; ei != n_extra_integers; ++ei)
            {
              const dof_id_type extra_int = cast_int<dof_id_type>(id_pid[2+ei]);
              libmesh_assert_equal_to(extra_int, old_node->get_extra_integer(ei));
            }
#endif

#ifdef LIBMESH_ENABLE_UNIQUE_ID
          libmesh_assert_equal_to(unique_id, old_node->unique_id());
#endif
        }
      else
        {
          Node * node =
            mesh.add_point(p, id, pid);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
          node->set_unique_id() = unique_id;
#endif

          libmesh_assert_equal_to(n_extra_integers, node->n_extra_integers());

          for (unsigned int ei=0; ei != n_extra_integers; ++ei)
            {
              const dof_id_type extra_int = cast_int<dof_id_type>(id_pid[2+ei]);
              node->set_extra_integer(ei, extra_int);
            }
        }
    }
}



template <typename file_id_type>
void CheckpointIO::read_connectivity (Xdr & io)
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  const bool read_extra_integers = this->version_at_least_1_5();

  const unsigned int n_extra_integers =
    read_extra_integers ? mesh.n_elem_integers() : 0;

  file_id_type n_elems_here;
  io.data(n_elems_here);

  // Keep track of the highest dimensional element we've added to the mesh
  unsigned int highest_elem_dim = 1;

  // RHS: Originally we used invalid_processor_id as a "no parent" tag
  // number, because I'm an idiot.  Let's try to support broken files
  // as much as possible.
  bool file_is_broken = false;

  for (unsigned int i=0; i<n_elems_here; i++)
    {
      // id type pid subdomain_id parent_id
      std::vector<file_id_type> elem_data(6 + n_extra_integers);
      io.data_stream
        (elem_data.data(), cast_int<unsigned int>(elem_data.size()),
         cast_int<unsigned int>(elem_data.size()));

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      file_id_type unique_id = 0;
      io.data(unique_id, "# unique id");
#endif

#ifdef LIBMESH_ENABLE_AMR
      uint16_t p_level = 0;
      io.data(p_level, "# p_level");

      uint16_t rflag, pflag;
      io.data(rflag, "# rflag");
      io.data(pflag, "# pflag");
#endif

      unsigned int n_nodes = Elem::type_to_n_nodes_map[elem_data[1]];

      // Snag the node ids this element was connected to
      std::vector<file_id_type> conn_data(n_nodes);
      io.data_stream
        (conn_data.data(), cast_int<unsigned int>(conn_data.size()),
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

      // Old broken files used processsor_id_type(-1)...
      // But we *know* our first element will be level 0
      if (i == 0 && elem_data[4] == 65535)
        file_is_broken = true;

      // On a broken file we can't tell whether a parent of 65535 is a
      // null parent or an actual parent of 65535.  Assuming the
      // former will cause less breakage.
      Elem * parent =
        (elem_data[4] == static_cast<largest_id_type>(-1) ||
         (file_is_broken && elem_data[4] == 65535)) ?
        nullptr : mesh.elem_ptr(cast_int<dof_id_type>(elem_data[4]));

      const unsigned short int child_num   =
        (elem_data[5] == static_cast<largest_id_type>(-1) ||
         (file_is_broken && elem_data[5] == 65535)) ?
        static_cast<unsigned short>(-1) :
        cast_int<unsigned short>(elem_data[5]);

      if (!parent)
        libmesh_assert_equal_to
          (child_num, static_cast<unsigned short>(-1));

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

          libmesh_assert_equal_to(n_extra_integers, old_elem->n_extra_integers());
#ifndef NDEBUG
          for (unsigned int ei=0; ei != n_extra_integers; ++ei)
            {
              const dof_id_type extra_int = cast_int<dof_id_type>(elem_data[6+ei]);
              libmesh_assert_equal_to(extra_int, old_elem->get_extra_integer(ei));
            }
#endif

          libmesh_assert_equal_to(old_elem->n_nodes(), conn_data.size());

          for (unsigned int n=0,
               n_conn = cast_int<unsigned int>(conn_data.size());
               n != n_conn; n++)
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
#else
          libmesh_ignore(child_num);
#endif

          libmesh_assert(elem->n_nodes() == conn_data.size());

          // Connect all the nodes to this element
          for (unsigned int n=0,
               n_conn = cast_int<unsigned int>(conn_data.size());
               n != n_conn; n++)
            elem->set_node(n) =
              mesh.node_ptr(cast_int<dof_id_type>(conn_data[n]));

          mesh.add_elem(elem);

          libmesh_assert_equal_to(n_extra_integers, elem->n_extra_integers());
          for (unsigned int ei=0; ei != n_extra_integers; ++ei)
            {
              const dof_id_type extra_int = cast_int<dof_id_type>(elem_data[6+ei]);
              elem->set_extra_integer(ei, extra_int);
            }
        }
    }

  mesh.set_mesh_dimension(cast_int<unsigned char>(highest_elem_dim));
}


template <typename file_id_type>
void CheckpointIO::read_remote_elem (Xdr & io, bool libmesh_dbg_var(expect_all_remote))
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Find the remote_elem neighbor links
  std::vector<file_id_type> elem_ids;
  std::vector<uint16_t> elem_sides;

  io.data(elem_ids, "# remote neighbor elem_ids");
  io.data(elem_sides, "# remote neighbor elem_sides");

  libmesh_assert_equal_to(elem_ids.size(), elem_sides.size());

  for (auto i : index_range(elem_ids))
    {
      Elem & elem = mesh.elem_ref(cast_int<dof_id_type>(elem_ids[i]));
      if (!elem.neighbor_ptr(elem_sides[i]))
        elem.set_neighbor(elem_sides[i],
                          const_cast<RemoteElem *>(RemoteElem::get_instance()));
      else
        libmesh_assert(!expect_all_remote);
    }

  // Find the remote_elem children links
  std::vector<file_id_type> parent_ids;
  std::vector<uint16_t> child_numbers;

  io.data(parent_ids, "# remote child parent_ids");
  io.data(child_numbers, "# remote child_numbers");

#ifdef LIBMESH_ENABLE_AMR
  for (auto i : index_range(parent_ids))
    {
      Elem & elem = mesh.elem_ref(cast_int<dof_id_type>(parent_ids[i]));

      // We'd like to assert that no child pointer already exists to
      // be overwritten by remote_elem, but Elem doesn't actually have
      // an API that will return a child pointer without asserting
      // that it isn't nullptr.
      const Elem * child = elem.raw_child_ptr(child_numbers[i]);

      if (!child)
        elem.add_child(const_cast<RemoteElem *>(RemoteElem::get_instance()),
                       child_numbers[i]);
      else
        libmesh_assert(!expect_all_remote);
    }
#endif
}



template <typename file_id_type>
void CheckpointIO::read_bcs (Xdr & io)
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // and our boundary info object
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  std::vector<file_id_type> element_id_list;
  std::vector<uint16_t> side_list;
  std::vector<file_id_type> bc_id_list;

  io.data(element_id_list, "# element ids for bcs");
  io.data(side_list, "# sides of elements for bcs");
  io.data(bc_id_list, "# bc ids");

  for (auto i : index_range(element_id_list))
    boundary_info.add_side
      (cast_int<dof_id_type>(element_id_list[i]), side_list[i],
       cast_int<boundary_id_type>(bc_id_list[i]));
}



template <typename file_id_type>
void CheckpointIO::read_nodesets (Xdr & io)
{
  // convenient reference to our mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // and our boundary info object
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  std::vector<file_id_type> node_id_list;
  std::vector<file_id_type> bc_id_list;

  io.data(node_id_list, "# node id list");
  io.data(bc_id_list, "# nodeset bc id list");

  for (auto i : index_range(node_id_list))
    boundary_info.add_node
      (cast_int<dof_id_type>(node_id_list[i]),
       cast_int<boundary_id_type>(bc_id_list[i]));
}



template <typename file_id_type>
void CheckpointIO::read_bc_names(Xdr & io, BoundaryInfo & info, bool is_sideset)
{
  std::map<boundary_id_type, std::string> & boundary_map = is_sideset ?
    info.set_sideset_name_map() : info.set_nodeset_name_map();

  std::vector<file_id_type> boundary_ids;
  std::vector<std::string>  boundary_names;

  file_id_type n_boundary_names = 0;

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
  for (auto i : index_range(boundary_ids))
    boundary_map[cast_int<boundary_id_type>(boundary_ids[i])] =
      boundary_names[i];
}


template <typename file_id_type>
void CheckpointIO::read_integers_names
  (Xdr & io,
   std::vector<std::string> & node_integer_names,
   std::vector<std::string> & elem_integer_names)
{
  file_id_type n_node_integers, n_elem_integers;

  io.data(n_node_integers, "# n_extra_integers per node");
  io.data(node_integer_names);
  io.data(n_elem_integers, "# n_extra_integers per elem");
  io.data(elem_integer_names);
}


unsigned int CheckpointIO::n_active_levels_in(MeshBase::const_element_iterator begin,
                                              MeshBase::const_element_iterator end) const
{
  unsigned int max_level = 0;

  for (const auto & elem : as_range(begin, end))
    max_level = std::max(elem->level(), max_level);

  return max_level + 1;
}

} // namespace libMesh
