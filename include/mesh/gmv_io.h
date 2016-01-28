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



#ifndef LIBMESH_GMV_IO_H
#define LIBMESH_GMV_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_output.h"
#include "libmesh/mesh_input.h"
#include "libmesh/enum_elem_type.h"

// C++ includes
#include <cstddef>
#include <cstring>  // for memcpy
#include <map>

namespace libMesh
{

// Forward declarations
class MeshBase;

/**
 * This class implements writing meshes in the GMV format.
 * For a full description of the GMV format and to obtain the
 * GMV software see
 * <a href="http://laws.lanl.gov/XCM/gmv/GMVHome.html">the GMV home page</a>
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
class GMVIO : public MeshInput<MeshBase>,
              public MeshOutput<MeshBase>
{
public:

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  explicit
  GMVIO (const MeshBase &);

  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This constructor is required to let us read in a mesh.
   */
  explicit
  GMVIO (MeshBase &);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string &) libmesh_override;

  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string & mesh_file) libmesh_override;

  //   /**
  //    * This method implements reading a mesh from a specified file.
  //    */
  //   virtual void read (const std::string & mesh_file)
  //   { this->read_mesh_and_nodal_data(mesh_file, libmesh_nullptr); }

  //   /**
  //    * Extension of the MeshInput::read() routine which
  //    * also takes an optional EquationSystems pointer and
  //    * tries to read field variables from the GMV file
  //    * into the EquationSystems object.
  //    */
  //   virtual void read_mesh_and_nodal_data (const std::string & ,
  //  EquationSystems * es=libmesh_nullptr);

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string &,
                                 const std::vector<Number> &,
                                 const std::vector<std::string> &) libmesh_override;

  /**
   * Flag indicating whether or not to write a binary file.  While binary
   * files may end up being smaller than equivalent ASCII files, they will
   * almost certainly take longer to write.  The reason for this is that
   * the ostream::write() function which is used to write "binary" data to
   * streams, only takes a pointer to char as its first argument.  This means
   * if you want to write anything other than a buffer of chars, you first
   * have to use a strange memcpy hack to get the data into the desired format.
   * See the templated to_binary_stream() function below.
   */
  bool & binary ();

  /**
   * Flag indicating whether or not to write the mesh
   * as discontinuous cell patches
   */
  bool & discontinuous();

  /**
   * Flag indicating whether or not to write the partitioning
   * information for the mesh.
   */
  bool & partitioning();

  /**
   * Flag to write element subdomain_id's as GMV "materials" instead
   * of element processor_id's.  Allows you to generate exploded views
   * on user-defined subdomains, potentially creating a pretty picture.
   */
  bool & write_subdomain_id_as_material();

  /**
   * Flag indicating whether or not to subdivide second order
   * elements
   */
  bool & subdivide_second_order();

  /**
   * Flag indicating whether or not to write p level
   * information for p refined meshes
   */
  bool & p_levels();

  /**
   * Writes a GMV file with discontinuous data
   */
  void write_discontinuous_gmv (const std::string & name,
                                const EquationSystems & es,
                                const bool write_partitioning,
                                const std::set<std::string> * system_names=libmesh_nullptr) const;


  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write an ASCII file.   This is the new implementation
   * (without subcells).
   */
  void write_ascii_new_impl (const std::string &,
                             const std::vector<Number> * = libmesh_nullptr,
                             const std::vector<std::string> * = libmesh_nullptr);

  /**
   * Takes a vector of cell-centered data to be plotted.
   * You must ensure that for every active element e,
   * v[e->id()] is a valid number.  You can add an arbitrary
   * number of different cell-centered data sets by calling
   * this function multiple times.
   *
   * .) GMV does not like spaces in the cell_centered_data_name
   * .) No matter what order you add cell-centered data, it will be
   *    output alphabetically.
   */
  void add_cell_centered_data (const std::string &       cell_centered_data_name,
                               const std::vector<Real> * cell_centered_data_vals);

  /**
   * If we read in a nodal solution while reading in a mesh, we can attempt
   * to copy that nodal solution into an EquationSystems object.
   */
  void copy_nodal_solution(EquationSystems & es);

private:

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write an ASCII file.  This is the old implementation
   * (using subcells) which was the default in libMesh-0.4.3-rc2.
   */
  void write_ascii_old_impl (const std::string &,
                             const std::vector<Number> * = libmesh_nullptr,
                             const std::vector<std::string> * = libmesh_nullptr);

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.
   */
  void write_binary (const std::string &,
                     const std::vector<Number> * = libmesh_nullptr,
                     const std::vector<std::string> * = libmesh_nullptr);

  /**
   * Helper function for writing unsigned ints to an ostream in binary format.
   * Implemented via memcpy as suggested in the standard.
   */
  template <typename T>
  void to_binary_stream(std::ostream & out,
                        const T i);

  /**
   * Flag to write binary data.
   */
  bool _binary;

  /**
   * Flag to write the mesh as discontinuous patches.
   */
  bool _discontinuous;

  /**
   * Flag to write the mesh partitioning.
   */
  bool _partitioning;

  /**
   * Flag to write element subdomain_id's as GMV "materials" instead
   * of element processor_id's.
   */
  bool _write_subdomain_id_as_material;

  /**
   * Flag to subdivide second order elements
   */
  bool _subdivide_second_order;

  /**
   * Flag to write the mesh p refinement levels.
   */
  bool _p_levels;

  /**
   * Storage for arbitrary cell-centered data.  Ex: You can use this
   * to plot the effectivity index for a given cell.  The map is
   * between the string representing the variable name and a pointer
   * to a vector containing the data.
   */
  std::map<std::string, const std::vector<Real> * > _cell_centered_data;

  /**
   * Helper functions for reading nodes/cells from a GMV file
   */
  void _read_nodes();
  unsigned int _next_elem_id;
  void _read_one_cell();
  ElemType _gmv_elem_to_libmesh_elem(const char * elemname);
  void _read_materials();
  void _read_var();
  std::map<std::string, std::vector<Number> > _nodal_data;
};



// ------------------------------------------------------------
// GMVIO inline members
inline
GMVIO::GMVIO (const MeshBase & mesh) :
  MeshOutput<MeshBase>    (mesh),
  _binary                 (false),
  _discontinuous          (false),
  _partitioning           (true),
  _write_subdomain_id_as_material (false),
  _subdivide_second_order (true),
  _p_levels               (true),
  _next_elem_id           (0)
{
}

inline
GMVIO::GMVIO (MeshBase & mesh) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _binary (false),
  _discontinuous          (false),
  _partitioning           (true),
  _write_subdomain_id_as_material (false),
  _subdivide_second_order (true),
  _p_levels               (true),
  _next_elem_id           (0)
{
}



inline
bool & GMVIO::binary ()
{
  return _binary;
}



inline
bool & GMVIO::discontinuous ()
{
  return _discontinuous;
}



inline
bool & GMVIO::partitioning ()
{
  return _partitioning;
}


inline
bool & GMVIO::write_subdomain_id_as_material ()
{
  return _write_subdomain_id_as_material;
}



inline
bool & GMVIO::subdivide_second_order ()
{
  return _subdivide_second_order;
}



inline
bool & GMVIO::p_levels()
{
  return _p_levels;
}



template <typename T>
void GMVIO::to_binary_stream(std::ostream & out_str,
                             const T i)
{
  static char buf[sizeof(T)];
  memcpy(buf, &i, sizeof(T));
  out_str.write(buf, sizeof(T));
}

} // namespace libMesh


#endif // LIBMESH_GMV_IO_H
