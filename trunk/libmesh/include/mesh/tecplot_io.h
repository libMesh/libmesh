// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __tecplot_io_h__
#define __tecplot_io_h__

// Local includes
#include "libmesh_common.h"
#include "mesh_output.h"

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
  TecplotIO (const MeshBase&, const bool binary=false);
  
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
   * Return/set the precision to use when writing ASCII files.
   */
  unsigned int & ascii_precision ();
  
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
   * Precision to use when writing ACII files.
   */
  unsigned int _ascii_precision;
};



// ------------------------------------------------------------
// TecplotIO inline members
inline
TecplotIO::TecplotIO (const MeshBase& mesh, const bool binary) :
  MeshOutput<MeshBase> (mesh),
  _binary (binary),
  _ascii_precision (6)
{
}



inline
bool & TecplotIO::binary ()
{
  return _binary;
}



inline
unsigned int & TecplotIO::ascii_precision ()
{
  return _ascii_precision;
}



#endif // #define __tecplot_io_h__
