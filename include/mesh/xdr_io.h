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



#ifndef LIBMESH_XDR_IO_H
#define LIBMESH_XDR_IO_H


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
class MeshBase;
class MeshData;
class Xdr;
class Elem;


/**
 *
 * @author Benjamin Kirk, John Peterson, 2004.
 */

// ------------------------------------------------------------
// XdrIO class definition
class XdrIO : public MeshInput<MeshBase>,
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
  XdrIO (MeshBase&,       const bool=false);

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   * The optional parameter \p binary can be used to switch
   * between ASCII (\p false, the default) or binary (\p true)
   * files.
   */
  explicit
  XdrIO (const MeshBase&, const bool=false);

  /**
   * Destructor.
   */
  virtual ~XdrIO ();

  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string&);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string&);

  /**
   * Get/Set the flag indicating if we should read/write binary.
   */
  bool   binary() const { return _binary; }
  bool & binary()       { return _binary; }

  /**
   * Get/Set the flag indicating if we should read/write legacy.
   */
  bool   legacy() const { return _legacy; }
  bool & legacy()       { return _legacy; }

  /**
   * Report whether we should write parallel files.
   */
  bool write_parallel() const;

  /**
   * Insist that we should/shouldn't write parallel files.
   */
  void set_write_parallel      (bool do_parallel = true);

  /**
   * Insist that we should write parallel files if and only if the
   * mesh is an already distributed ParallelMesh.
   */
  void set_auto_parallel ();

  /**
   * Get/Set the version string.  Vailid version strings:
     \verbatim

     "libMesh-0.7.0+"
     "libMesh-0.7.0+ parallel"

     \endverbatim
     If "libMesh" is not detected in the version string the
     \p LegacyXdrIO class will be used to read older
     (pre version 0.7.0) mesh files.
   */
  const std::string & version () const { return _version; }
  std::string &       version ()       { return _version; }

  /**
   * Get/Set the boundary condition file name.
   */
  const std::string & boundary_condition_file_name() const { return _bc_file_name; }
  std::string &       boundary_condition_file_name()       { return _bc_file_name; }

  /**
   * Get/Set the partitioning file name.
   */
  const std::string & partition_map_file_name() const { return _partition_map_file; }
  std::string &       partition_map_file_name()       { return _partition_map_file; }

  /**
   * Get/Set the subdomain file name.
   */
  const std::string & subdomain_map_file_name() const { return _subdomain_map_file; }
  std::string &       subdomain_map_file_name()       { return _subdomain_map_file; }

  /**
   * Get/Set the polynomial degree file name.
   */
  const std::string & polynomial_level_file_name() const { return _p_level_file; }
  std::string &       polynomial_level_file_name()       { return _p_level_file; }


 private:


  //---------------------------------------------------------------------------
  // Write Implementation
  /**
   * Write subdomain name information - NEW in 0.9.2 format
   */
  void write_serialized_subdomain_names(Xdr &io) const;

  /**
   * Write the connectivity for a parallel, distributed mesh
   */
  void write_serialized_connectivity (Xdr &io, const dof_id_type n_elem) const;

  /**
   * Write the nodal locations for a parallel, distributed mesh
   */
  void write_serialized_nodes (Xdr &io, const dof_id_type n_nodes) const;

  /**
   * Write the boundary conditions for a parallel, distributed mesh
   */
  void write_serialized_bcs (Xdr &io, const std::size_t n_bcs) const;

  /**
   * Write the boundary conditions for a parallel, distributed mesh
   */
  void write_serialized_nodesets (Xdr &io, const std::size_t n_nodesets) const;

  /**
   * Write boundary names information (sideset and nodeset) - NEW in 0.9.2 format
   */
  void write_serialized_bc_names (Xdr &io, const BoundaryInfo & info, bool is_sideset) const;


  //---------------------------------------------------------------------------
  // Read Implementation
  /**
   * Read subdomain name information - NEW in 0.9.2 format
   */
  void read_serialized_subdomain_names(Xdr &io);

  /**
   * Read the connectivity for a parallel, distributed mesh
   */
  template <typename T>
  void read_serialized_connectivity (Xdr &io, const dof_id_type n_elem, std::vector<header_id_type> & sizes, T);

  /**
   * Read the nodal locations for a parallel, distributed mesh
   */
  void read_serialized_nodes (Xdr &io, const dof_id_type n_nodes);

  /**
   * Read the boundary conditions for a parallel, distributed mesh
   * @return the number of bcs read
   */
  template <typename T>
  void read_serialized_bcs (Xdr &io, T);

  /**
   * Read the nodeset conditions for a parallel, distributed mesh
   * @return the number of nodesets read
   */
  template <typename T>
  void read_serialized_nodesets (Xdr &io, T);

  /**
   * Read boundary names information (sideset and nodeset) - NEW in 0.9.2 format
   */
  void read_serialized_bc_names(Xdr &io, BoundaryInfo & info, bool is_sideset);

  //-------------------------------------------------------------------------
  /**
   * Pack an element into a transfer buffer for parallel communication.
   */
  void pack_element (std::vector<xdr_id_type> &conn,
		     const Elem *elem,
		     const dof_id_type parent_id  = DofObject::invalid_id,
		     const dof_id_type parent_pid = DofObject::invalid_id) const;

  bool _binary;
  bool _legacy;
  bool _write_serial;
  bool _write_parallel;
  bool _write_unique_id;
  header_id_type _field_width;
  std::string _version;
  std::string _bc_file_name;
  std::string _partition_map_file;
  std::string _subdomain_map_file;
  std::string _p_level_file;

  /**
   * Define the block size to use for chunked IO.
   */
  static const std::size_t io_blksize;
};


// ------------------------------------------------------------
// XdrIO inline members

inline
bool XdrIO::write_parallel() const
{
  // We can't insist on both serial and parallel
  libmesh_assert (!this->_write_serial || !this->_write_parallel);

  // If we insisted on serial, do that
  if (this->_write_serial)
    return false;

  // If we insisted on parallel, do that
  if (this->_write_parallel)
    return true;

  // If we're doing things automatically, check the mesh
  const MeshBase &mesh = MeshOutput<MeshBase>::mesh();
  return !mesh.is_serial();
}



inline
void XdrIO::set_write_parallel (bool do_parallel)
{
  this->_write_parallel = do_parallel;

  this->_write_serial = !do_parallel;
}



inline
void XdrIO::set_auto_parallel ()
{
  this->_write_serial   = false;
  this->_write_parallel = false;
}


} // namespace libMesh



#endif // LIBMESH_XDR_IO_H
