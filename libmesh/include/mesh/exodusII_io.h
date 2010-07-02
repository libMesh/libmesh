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



#ifndef __exodusII_io_h__
#define __exodusII_io_h__


// C++ inludes

// Local includes
#include "libmesh_common.h"
#include "mesh_input.h"
#include "mesh_output.h"
#include "exodusII_io_helper.h"

namespace libMesh
{

// Forward declarations
class EquationSystems;
class MeshBase;
class System;

/**
 * The \p ExodusII_IO class implements reading meshes in the
 * \p ExodusII file format from Sandia National Labs.  By
 * default, LibMesh expects ExodusII files to have a ".exd"
 * or ".e" file extension.
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
  virtual void write (const std::string& fname);


  /**
   * Set the flag indicationg if we should be verbose.
   */
  void verbose (bool set_verbosity);

  /**
   * If we read in a nodal solution while reading in a mesh, we can attempt
   * to copy that nodal solution into an EquationSystems object.
   */
  void copy_nodal_solution(System& es, std::string nodal_var_name, unsigned int timestep=1);
  
  /**
   * Writes a exodusII file with discontinuous data
   */ 
  void write_discontinuous_exodusII (const std::string& name, 
				const EquationSystems& es);
  

  /**
   * Write out a nodal solution.
   */
  void write_nodal_data (const std::string&,
			 const std::vector<Number>&,
			 const std::vector<std::string>&);

  /**
   * Write out a discontinuous nodal solution.
   */
  void write_nodal_data_discontinuous (const std::string&,
			 const std::vector<Number>&,
			 const std::vector<std::string>&);

  /**
   * Writes out the solution at a specific timestep.
   * @param timestep The timestep to write out, should be _1_ indexed.
   */
  void write_timestep (const std::string& fname,
		       const EquationSystems& es,
		       const int timestep,
		       const Real time);
  
 private:
  /**
   * Only attempt to instantiate an ExodusII helper class
   * if the Exodus API is defined.  This class will have no
   * functionality when LIBMESH_HAVE_EXODUS_API is not defined.
   */
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO_Helper exio_helper;
#endif

  /**
   *
   */
  int _timestep;

  /**
   * should we be verbose?
   */
  bool _verbose;
};


} // namespace libMesh


#endif // #define __exodusII_io_h__
