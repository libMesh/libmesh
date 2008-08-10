// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __exodusII_io_h__
#define __exodusiI_io_h__


// C++ inludes

// Local includes
#include "mesh_input.h"
#include "mesh_output.h"

// Forward declarations
class MeshBase;
class System;
class ExodusII;

/**
 * The \p ExodusII_IO class implements reading meshes in the
 * \p ExodusII file format from Sandia National Labs.  Due to
 * licensing restrictions the external libraries required to
 * read such files cannot be redistributed with \p libMesh.
 * However, if you have these libraries, this class will allow
 * you to use the API to read ExodusII, or \p .exd/.e, files.
 *
 * @author Benjamin Kirk, John Peterson, 2004.
 */

// ------------------------------------------------------------
// ExodusII_IO class definition
class ExodusII_IO : public MeshInput<MeshBase>,
		    public MeshOutput<MeshBase>
{

 public:
  
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  ExodusII_IO (MeshBase& mesh);
  
  /**
   * Destructor.
   */
  virtual ~ExodusII_IO ();
  
  /**
   * This method implements reading a mesh from a specified file.
   * Open the file named \p name and read the mesh in Sandia National Lab's
   * ExodusII format. This is the method to use for reading in meshes generated
   * by cubit.  Works in 2D for \p TRIs, \p TRI6s, \p QUAD s, and \p QUAD9s.
   * Works in 3D for \p TET4s, \p TET10s, \p HEX8s, and \p HEX27s.
   */
  virtual void read (const std::string& name);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string& )
  {
    std::cerr << "The ExodusII_IO::write() method is not yet implemented." << std::endl;
    libmesh_error();
  }


  /**
   * Set the flag indicationg if we should be verbose.
   */
  bool & verbose ();

  /**
   * If we read in a nodal solution while reading in a mesh, we can attempt
   * to copy that nodal solution into an EquationSystems object.
   */
  void copy_nodal_solution(System& es, std::string nodal_var_name);

  /**
   * Write out a nodal solution.
   */
  void write_nodal_data (const std::string&,
			 const std::vector<Number>&,
			 const std::vector<std::string>&);

  /**
   * Writes out the solution at a specific timestep.
   * @param timestep The timestep to write out, should be _1_ indexed.
   */
  void write_timestep (const std::string& fname,
		       const EquationSystems& es,
		       const int timestep,
		       const double time);
  
 private:
  ExodusII * ex_ptr;
  int _timestep;

//-------------------------------------------------------------
  // local data

  /**
   * should be be verbose?
   */
  bool _verbose;
};

inline
bool & ExodusII_IO::verbose ()
{
  return _verbose;
}


#endif // #define __exodusII_io_h__
