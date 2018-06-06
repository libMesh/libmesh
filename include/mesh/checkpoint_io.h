// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CHECKPOINT_IO_H
#define LIBMESH_CHECKPOINT_IO_H


// Local includes
#include "libmesh/compare_elems_by_level.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"
#include "libmesh/parallel_object.h"

// C++ includes
#include <string>
#include <vector>

namespace libMesh
{
// Forward declarations
class Xdr;
class CheckpointIO;

/**
 * split_mesh takes the given initialized/opened mesh and partitions it into nsplits pieces or
 * chunks.  It returns a CheckpointIO object that can be used to write the mesh chunks into
 * individual files (e.g. by calling checkpoint_obj.write(out_file_name)) - the number of files is
 * equal to the number of chunks.  This function supports MPI parallelism and can be used with
 * several MPI procs to speed up splitting.
 */
std::unique_ptr<CheckpointIO> split_mesh(MeshBase & mesh, unsigned int nsplits);

/**
 * The CheckpointIO class can be used to write simplified restart
 * files that can be used to restart simulations that have
 * crashed.
 *
 * \author Benjamin Kirk
 * \author John Peterson
 * \author Derek Gaston
 * \author Roy Stogner
 * \date 2017
 */
class CheckpointIO : public MeshInput<MeshBase>,
                     public MeshOutput<MeshBase>,
                     public ParallelObject
{
public:
  // The size used for encoding all id types in this file
  typedef largest_id_type xdr_id_type;

  // The size type used to read header sizes (meta data information)
  typedef uint64_t header_id_type;

  /**
   * Constructor.  Takes a writable reference to a mesh object.
   * This is the constructor required to read a mesh.
   * The optional parameter \p binary can be used to switch
   * between ASCII (\p false, the default) or binary (\p true)
   * files.
   */
  explicit
  CheckpointIO (MeshBase &, const bool=false);

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   * The optional parameter \p binary can be used to switch
   * between ASCII (\p false, the default) or binary (\p true)
   * files.
   */
  explicit
  CheckpointIO (const MeshBase &, const bool=false);

  /**
   * Destructor.
   */
  virtual ~CheckpointIO ();

  /**
   * This method implements reading a mesh from a specified file.  If the mesh has been split for
   * running on several processors, input_name should simply be the name of the mesh split
   * directory without the "-split[n]" suffix.  The number of splits will be determined
   * automatically by the number of processes being used for the mesh at the time of reading.
   */
  virtual void read (const std::string & input_name) override;

  /**
   * This method implements writing a mesh to a specified file.  If the mesh has been split
   * for running on several processors, this will create a subdirectory named
   * "[name]-split[n]" where name is the given name argument and n is the number of
   * processors the mesh is split for running on.  For example:
   *
   *     unsigned int n_splits = 42;
   *     std::unique_ptr<CheckpointIO> cp = split_mesh(my_mesh, n_splits);
   *     // ...
   *     cp->write("foo.cpr");
   *
   * would create a directory named "foo.cpr-split42".
   */
  virtual void write (const std::string & name) override;

  /**
   * Used to remove a checkpoint directory and its corresponding files.  This effectively undoes
   * all the work done be calls to write(...).  For example, if a checkpoint configuration was
   * written via:
   *
   *     unsigned int n_splits = 42;
   *     std::unique_ptr<CheckpointIO> cp = split_mesh(my_mesh, n_splits);
   *     // ...
   *     cp->write("foo.cpr");
   *
   * then you could remove all the corresponding created files/dirs by:
   *
   *     CheckpointIO::cleanup(my_mesh, n_splits);
   *
   * Or for cases where the split configuration was determined automatically (e.g. via number of
   * running procs with distributed/parallel mesh), then you could:
   *
   *     CheckpointIO::cleanup(your_mesh, your_mesh.comm().size());
   *
   * Other remaining checkpoint split configurations for the mesh are left unmodified.
   */
  static void cleanup(const std::string & input_name, processor_id_type n_procs);

  /**
   * Get/Set the flag indicating if we should read/write binary.
   */
  bool   binary() const { return _binary; }
  bool & binary()       { return _binary; }

  /**
   * Get/Set the flag indicating if we should read/write binary.
   */
  bool   parallel() const { return _parallel; }
  bool & parallel()       { return _parallel; }

  /**
   * Get/Set the version string.
   */
  const std::string & version () const { return _version; }
  std::string &       version ()       { return _version; }

  /**
   * Get/Set the processor id or processor ids to use.
   *
   * The default processor_id to use is the processor_id() of the
   * mesh.
   *
   * This is used for m->n parallel checkpoint file writing:
   * You can force CheckpointIO to write out different partitions of a
   * mesh by setting which partitions to write from each processor here.
   */
  const std::vector<processor_id_type> & current_processor_ids() const { return _my_processor_ids; }
  std::vector<processor_id_type> & current_processor_ids() { return _my_processor_ids; }

  /**
   * Get/Set the n_processors to use
   *
   * The default n_processors to use is the n_processors() of the
   * mesh.
   *
   * This is used for m->n parallel checkpoint file writing:
   * You can force CheckpointIO to view the world as if it contains this number of
   * processors by setting it here
   */
  const processor_id_type & current_n_processors() const { return _my_n_processors; }
  processor_id_type & current_n_processors() { return _my_n_processors; }

private:
  //---------------------------------------------------------------------------
  // Write Implementation

  /**
   * Write subdomain name information
   */
  void write_subdomain_names(Xdr & io) const;

  /**
   * Write the connectivity for part of a mesh
   */
  void write_connectivity (Xdr & io,
                           const std::set<const Elem *, CompareElemIdsByLevel> & elements) const;

  /**
   * Write the remote_elem neighbor and child links for part of a mesh
   */
  void write_remote_elem (Xdr & io,
                          const std::set<const Elem *, CompareElemIdsByLevel> & elements) const;

  /**
   * Write the nodal locations for part of a mesh
   */
  void write_nodes (Xdr & io,
                    const std::set<const Node *> & nodeset) const;

  /**
   * Write the side boundary conditions for part of a mesh
   */
  void write_bcs (Xdr & io,
                  const std::set<const Elem *, CompareElemIdsByLevel> & elements) const;

  /**
   * Write the nodal boundary conditions for part of a mesh
   */
  void write_nodesets (Xdr & io,
                       const std::set<const Node *> & nodeset) const;

  /**
   * Write boundary names information (sideset and nodeset)
   */
  void write_bc_names (Xdr & io, const BoundaryInfo & info, bool is_sideset) const;


  //---------------------------------------------------------------------------
  // Read Implementation

  /**
   * Read header data on processor 0, then broadcast.  Returns the
   * number of processors for which parallel non-header files have
   * been written, or 0 if files were written in serial.
   */
  template <typename file_id_type>
  file_id_type read_header(const std::string & name);

  /**
   * Read a non-header file
   */
  template <typename file_id_type>
  void read_subfile(Xdr & io, bool expect_all_remote);

  /**
   * Read subdomain name information
   */
  template <typename file_id_type>
  void read_subdomain_names(Xdr & io);

  /**
   * Read the connectivity for a parallel, distributed mesh
   */
  template <typename file_id_type>
  void read_connectivity (Xdr & io);

  /**
   * Read the remote_elem neighbor and child links for a parallel, distributed mesh
   *
   * If we expect all these remote_elem links to truly be remote,
   * because we aren't doing an N -> M restart with M < N, then we set
   * \p expect_all_remote to true and test more assertions.
   */
  template <typename file_id_type>
  void read_remote_elem (Xdr & io, bool expect_all_remote);

  /**
   * Read the nodal locations for a parallel, distributed mesh
   */
  template <typename file_id_type>
  void read_nodes (Xdr & io);

  /**
   * Read the boundary conditions for a parallel, distributed mesh
   */
  template <typename file_id_type>
  void read_bcs (Xdr & io);

  /**
   * Read the nodeset conditions for a parallel, distributed mesh
   */
  template <typename file_id_type>
  void read_nodesets (Xdr & io);

  /**
   * Read boundary names information (sideset and nodeset)
   */
  template <typename file_id_type>
  void read_bc_names(Xdr & io, BoundaryInfo & info, bool is_sideset);

  /**
   * \returns The number of levels of refinement in the active mesh on
   * this processor.
   *
   * \note This includes _all_ elements on this processor even those
   * not owned by this processor!  Implemented by looping over all the
   * active elements and finding the maximum level.
   */
  unsigned int n_active_levels_in(MeshBase::const_element_iterator begin,
                                  MeshBase::const_element_iterator end) const;

  processor_id_type select_split_config(const std::string & input_name, header_id_type & data_size);

  bool _binary;
  bool _parallel;
  std::string _version;

  // The processor ids to write
  std::vector<processor_id_type> _my_processor_ids;

  // The largest processor id to write
  processor_id_type _my_n_processors;
};


} // namespace libMesh

#endif // LIBMESH_CHECKPOINT_IO_H
