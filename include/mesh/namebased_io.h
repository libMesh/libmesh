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



#ifndef LIBMESH_NAMEBASED_IO_H
#define LIBMESH_NAMEBASED_IO_H

// Local includes
#include "libmesh/mesh_output.h"
#include "libmesh/mesh_input.h"

namespace libMesh
{

// Forward declarations
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;

/**
 * This class supports simple reads and writes in any
 * libMesh-supported format, by dispatching to one of the other I/O
 * classes based on filename.
 *
 * Other I/O classes may have more advanced features that are not
 * accessible via this interface.
 *
 * \author Roy H. Stogner
 * \date 2015
 */
class NameBasedIO : public MeshInput<MeshBase>,
                    public MeshOutput<MeshBase>
{
public:

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  explicit
  NameBasedIO (const MeshBase &);

  /**
   * Constructor.  Takes a writable reference to a mesh object.
   * This constructor is required to let us read in a mesh.
   */
  explicit
  NameBasedIO (MeshBase &);

  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string & mesh_file) override;

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string & mesh_file) override;

  /**
   * This method implements writing a mesh with data to a specified file
   * where the data is taken from the \p EquationSystems object.
   *
   * We override the default MeshOutput::write_equation_systems
   * because it only outputs nodal data by default, whereas we want to
   * output a proper restart file if the requested filename is an XDA
   * or XDR type.
   */
  virtual void write_equation_systems (const std::string & filename,
                                       const EquationSystems & es,
                                       const std::set<std::string> * system_names=nullptr) override;

  /**
   * Bring in base class functionality for name resolution and to
   * avoid warnings about hidden overloaded virtual functions.
   */
  using MeshOutput<MeshBase>::write_nodal_data;

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string &,
                                 const std::vector<Number> &,
                                 const std::vector<std::string> &) override;

  // Certain mesh formats can support parallel I/O, including the
  // "new" Xdr format and the Nemesis format.
  bool is_parallel_file_format (const std::string & filename);
};



// ------------------------------------------------------------
// NameBasedIO inline members
inline
NameBasedIO::NameBasedIO (const MeshBase & mesh) :
  MeshOutput<MeshBase>    (mesh)
{
}

inline
NameBasedIO::NameBasedIO (MeshBase & mesh) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh)
{
}

inline
bool
NameBasedIO::is_parallel_file_format (const std::string & name)
{
  return ((name.rfind(".xda") < name.size()) ||
          (name.rfind(".xdr") < name.size()) ||
          (name.rfind(".nem") + 4 == name.size()) ||
          (name.rfind(".n") + 2 == name.size()) ||
          (name.rfind(".cp") < name.size())
          );
}


} // namespace libMesh


#endif // LIBMESH_NAMEBASED_IO_H
