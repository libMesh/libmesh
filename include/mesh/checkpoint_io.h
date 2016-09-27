// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh.h"
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

/**
 * The CheckpointIO class can be used to write simplified restart
 * files that can be used to restart simulations that have
 * crashed. Only N-to-N (procs) restart is supported with CheckpointIO
 * files.
 *
 * \author Benjamin Kirk
 * \author John Peterson
 * \author Derek Gaston
 * \date 2013
 */
class CheckpointIO : public MeshInput<MeshBase>,
                     public MeshOutput<MeshBase>,
                     public ParallelObject
{
public:
  // The size used for encoding all id types in this file
  typedef largest_id_type xdr_id_type;

  // The size type used to read header sizes (meta data information)
  typedef uint32_t header_id_type;

  /**
   * Constructor.  Takes a writeable reference to a mesh object.
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
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string &) libmesh_override;

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string &) libmesh_override;

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
   * Get/Set the processor_id to use
   *
   * This is used for m->n parallel checkpoint file writing:
   * You can force CheckpointIO to view the world as if it is on a particular
   * processor_id by setting it here
   */
  const processor_id_type & current_processor_id() const { return _my_processor_id; }
  processor_id_type & current_processor_id() { return _my_processor_id; }

  /**
   * Get/Set the n_processors to use
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
   * Build up the elem list
   */
  void build_elem_list();

  /**
   * Build up the node list
   */
  void build_node_list();

  /**
   * Write subdomain name information - NEW in 0.9.2 format
   */
  void write_subdomain_names(Xdr & io) const;

  /**
   * Write the connectivity for a parallel, distributed mesh
   */
  void write_connectivity (Xdr & io) const;

  /**
   * Write the nodal locations for a parallel, distributed mesh
   */
  void write_nodes (Xdr & io) const;

  /**
   * Write the boundary conditions for a parallel, distributed mesh
   */
  void write_bcs (Xdr & io) const;

  /**
   * Write the boundary conditions for a parallel, distributed mesh
   */
  void write_nodesets (Xdr & io) const;

  /**
   * Write boundary names information (sideset and nodeset) - NEW in 0.9.2 format
   */
  void write_bc_names (Xdr & io, const BoundaryInfo & info, bool is_sideset) const;


  //---------------------------------------------------------------------------
  // Read Implementation
  /**
   * Read subdomain name information - NEW in 0.9.2 format
   */
  void read_subdomain_names(Xdr & io);

  /**
   * Read the connectivity for a parallel, distributed mesh
   */
  void read_connectivity (Xdr & io);

  /**
   * Read the nodal locations for a parallel, distributed mesh
   */
  void read_nodes (Xdr & io);

  /**
   * Read the boundary conditions for a parallel, distributed mesh
   */
  void read_bcs (Xdr & io);

  /**
   * Read the nodeset conditions for a parallel, distributed mesh
   */
  void read_nodesets (Xdr & io);

  /**
   * Read boundary names information (sideset and nodeset) - NEW in 0.9.2 format
   */
  void read_bc_names(Xdr & io, BoundaryInfo & info, bool is_sideset);

  /**
   * Return the number of levels of refinement in the active mesh on this processor.
   * NOTE: This includes _all_ elements on this processor even those not owned by this processor!
   * Implemented by looping over all the active elements and finding
   * the maximum level.
   */
  unsigned int n_active_levels_on_processor(const MeshBase & mesh) const;

  bool _binary;
  bool _parallel;
  std::string _version;
  unsigned int _mesh_dimension;

  /// These are sets of IDs to make the lookup for boundary conditions simpler
  std::set<largest_id_type> _local_elements;
  std::set<largest_id_type> _nodes_connected_to_local_elements;

  /// The processor_id to use
  processor_id_type _my_processor_id;

  /// The number of processors to use
  processor_id_type _my_n_processors;
};


} // namespace libMesh

#endif // LIBMESH_CHECKPOINT_IO_H
