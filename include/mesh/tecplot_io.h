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



#ifndef LIBMESH_TECPLOT_IO_H
#define LIBMESH_TECPLOT_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_output.h"

// C++ Includes
#include <cstddef>
#include <set>

namespace libMesh
{

// Forward declarations
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;

/**
 * This class implements writing meshes in the Tecplot format.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
class TecplotIO : public MeshOutput<MeshBase>
{
public:

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   * The optional parameter \p binary can be used to switch
   * between ASCII (\p false, the default) or binary (\p true)
   * output files.
   */
  explicit
  TecplotIO (const MeshBase &, const bool binary=false,
             const double time=0., const int strand_offset=0);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string &) override;

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

  /**
   * Flag indicating whether or not to write a binary file
   * (if the tecio.a library was found by \p configure).
   */
  bool & binary ();

  /**
   * Solution time for transient data.
   * Written to newer binary formats that are time-aware.
   */
  double & time ();

  /**
   * Strand offset for this file.  Each mesh block will
   * be written to (strand_id=block_id+1+strand_offset).
   * Written to newer binary formats that are time-aware,
   * defaults to 0.
   */
  int & strand_offset ();

  /**
   *  The zone title to write.
   */
  std::string & zone_title ();

  /**
   * Set to true to write multiple solutions to a single file (ASCII
   * only).  Tecplot will read multiple zones in a single file, but
   * currently you have to repeat the mesh information each time.
   */
  bool & ascii_append ();

private:

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write an ASCII file.
   */
  void write_ascii (const std::string &,
                    const std::vector<Number> * = nullptr,
                    const std::vector<std::string> * = nullptr);

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write a binary file if the tecio.a library was
   * found at compile time, otherwise a warning message will be printed and
   * an ASCII file will be created.
   */
  void write_binary (const std::string &,
                     const std::vector<Number> * = nullptr,
                     const std::vector<std::string> * = nullptr);

  /**
   * Determines the logical spatial dimension of the elements in the
   * Mesh.  Ex: A 1D edge element living in 3D is a logically
   * one-dimensional element as far as Tecplot is concerned.  Throws
   * an error if mixed-dimension element types are found, since I'm
   * not sure how to handle that case currently.
   */
  unsigned elem_dimension();

  //---------------------------------------------------------------------------
  // local data

  /**
   * Flag to write binary data.
   */
  bool _binary;

  /**
   * Solution time.
   */
  double _time;

  /**
   * Offset for Tecplot's STRANDID.
   */
  int _strand_offset;

  /**
   * The zone title to write.
   */
  std::string _zone_title;

  /**
   * If true, when writing in ASCII format, open the file in
   * std::ofstream::app mode.
   */
  bool _ascii_append;

  /**
   * The subdomains in the mesh.
   */
  std::set<subdomain_id_type> _subdomain_ids;
};

} // namespace libMesh


#endif // LIBMESH_TECPLOT_IO_H
