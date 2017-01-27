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



#ifndef LIBMESH_MESH_OUTPUT_H
#define LIBMESH_MESH_OUTPUT_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_serializer.h"

// C++ includes
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace libMesh
{

// Forward declares
class EquationSystems;
template <typename T> class NumericVector;


/**
 * This class defines an abstract interface for \p Mesh output.
 * Specific classes derived from this class actually implement
 * writing various mesh formats.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
template <class MT>
class MeshOutput
{
protected:

  /**
   * Default constructor. Will set the _obj to NULL, effectively
   * rendering this object useless.
   */
  explicit
  MeshOutput (const bool is_parallel_format = false, const bool serial_only_needed_on_proc_0 = false);

  /**
   * Constructor.  Takes a reference to a constant object.
   * This constructor will only allow us to write the object.
   */
  explicit
  MeshOutput (const MT &, const bool is_parallel_format = false, const bool serial_only_needed_on_proc_0 = false);


public:

  /**
   * Destructor.
   */
  virtual ~MeshOutput ();

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string &) = 0;

  /**
   * This method implements writing a mesh with data to a specified file
   * where the data is taken from the \p EquationSystems object.
   */
  virtual void write_equation_systems (const std::string &,
                                       const EquationSystems &,
                                       const std::set<std::string> * system_names=libmesh_nullptr);

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string &,
                                 const std::vector<Number> &,
                                 const std::vector<std::string> &)
  { libmesh_not_implemented(); }

  /**
   * This method should be overridden by "parallel" output formats for
   * writing nodal data.  Instead of getting a localized copy of the
   * nodal solution vector, it is passed a NumericVector of
   * type=PARALLEL which is in node-major order i.e.
   * (u0,v0,w0, u1,v1,w1, u2,v2,w2, u3,v3,w3, ...)
   * and contains n_nodes*n_vars total entries.  Then, it is up to the
   * individual I/O class to extract the required solution values from
   * this vector and write them in parallel.
   *
   * If not implemented, localizes the parallel vector into a std::vector
   * and calls the other version of this function.
   */
  virtual void write_nodal_data (const std::string &,
                                 const NumericVector<Number> &,
                                 const std::vector<std::string> &);

  /**
   * Return/set the precision to use when writing ASCII files.
   *
   * By default we use numeric_limits<Real>::digits10 + 2, which
   * should be enough to write out to ASCII and get the exact same
   * Real back when reading in.
   */
  unsigned int & ascii_precision ();

protected:


  /**
   * Returns the object as a read-only reference.
   */
  const MT & mesh() const;


  /**
   * Flag specifying whether this format is parallel-capable.
   * If this is false (default) I/O is only permitted when the mesh
   * has been serialized.
   */
  const bool _is_parallel_format;

  /**
   * Flag specifying whether this format can be written by only
   * serializing the mesh to processor zero
   *
   * If this is false (default) the mesh will be serialized to
   * all processors
   */
  const bool _serial_only_needed_on_proc_0;

private:


  /**
   *  A pointer to a constant object.
   * This allows us to write the object to file.
   */
  const MT * const _obj;

  /**
   * Precision to use when writing ASCII files.
   */
  unsigned int _ascii_precision;
};






// ------------------------------------------------------------
// MeshOutput inline members
template <class MT>
inline
MeshOutput<MT>::MeshOutput (const bool is_parallel_format, const bool serial_only_needed_on_proc_0) :
  _is_parallel_format(is_parallel_format),
  _serial_only_needed_on_proc_0(serial_only_needed_on_proc_0),
  _obj(libmesh_nullptr),
  _ascii_precision (std::numeric_limits<Real>::digits10 + 2)
{}



template <class MT>
inline
MeshOutput<MT>::MeshOutput (const MT & obj, const bool is_parallel_format, const bool serial_only_needed_on_proc_0) :
  _is_parallel_format(is_parallel_format),
  _serial_only_needed_on_proc_0(serial_only_needed_on_proc_0),
  _obj (&obj),
  _ascii_precision (std::numeric_limits<Real>::digits10 + 2)
{
  if (!_is_parallel_format && !this->mesh().is_serial())
    {
      if (this->mesh().processor_id() == 0)
        {
          libmesh_do_once(libMesh::out <<
                          "Warning:  This MeshOutput subclass only supports meshes which have been serialized!"
                          << std::endl;);
        }
    }
}



template <class MT>
inline
MeshOutput<MT>::~MeshOutput ()
{
}



template <class MT>
inline
const MT & MeshOutput<MT>::mesh () const
{
  libmesh_assert(_obj);
  return *_obj;
}



template <class MT>
inline
unsigned int & MeshOutput<MT>::ascii_precision ()
{
  return _ascii_precision;
}


} // namespace libMesh


#endif // LIBMESH_MESH_OUTPUT_H
