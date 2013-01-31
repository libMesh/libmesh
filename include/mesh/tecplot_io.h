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
class MeshBase;



/**
 * This class implements writing meshes in the Tecplot format.
 *
 * @author Benjamin S. Kirk, 2004
 */

// ------------------------------------------------------------
// TecplotIO class definition
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
  TecplotIO (const MeshBase&, const bool binary=false,
	     const double time=0., const int strand_offset=0);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string& );

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string&,
				 const std::vector<Number>&,
				 const std::vector<std::string>&);

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
  
 private:

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write an ASCII file.
   */
  void write_ascii (const std::string&,
		    const std::vector<Number>* = NULL,
		    const std::vector<std::string>* = NULL);

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write a binary file if the tecio.a library was
   * found at compile time, otherwise a warning message will be printed and
   * an ASCII file will be created.
   */
  void write_binary (const std::string&,
		     const std::vector<Number>* = NULL,
		     const std::vector<std::string>* = NULL);

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
   * The subdomains in the mesh.
   */
  std::set<subdomain_id_type> _subdomain_ids;
};



// ------------------------------------------------------------
// TecplotIO inline members
inline
bool & TecplotIO::binary ()
{
  return _binary;
}



inline
double & TecplotIO::time ()
{
  return _time;
}



inline
int & TecplotIO::strand_offset ()
{
  return _strand_offset;
}



inline
std::string & TecplotIO::zone_title ()
{
  return _zone_title;
}


} // namespace libMesh


#endif // LIBMESH_TECPLOT_IO_H
