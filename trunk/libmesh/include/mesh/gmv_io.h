// $Id: gmv_io.h,v 1.3 2004-06-15 20:17:20 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __gmv_io_h__
#define __gmv_io_h__

// Local includes
#include "libmesh_common.h"
#include "mesh_io.h"

// Forward declarations
class MeshBase;



/**
 * This class implements writing meshes in the GMV format.
 * For a full description of the GMV format and to obtain the
 * GMV software see
 * <a href="http://laws.lanl.gov/XCM/gmv/GMVHome.html">the GMV home page</a>
 *
 * @author Benjamin S. Kirk, 2004
 */

// ------------------------------------------------------------
// GMVIO class definition
class GMVIO : public MeshIO<MeshBase>
{
 public:

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  GMVIO (const MeshBase&);
  
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
   * Writes a GMV file with discontinuous data
   */ 
  void write_discontinuous_gmv (const std::string& name, 
				const EquationSystems& es,
				const bool write_partitioning) const;
  

private:

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write an ASCII file.
   */
  virtual void write_ascii (const std::string&,
			    const std::vector<Number>* = NULL,
			    const std::vector<std::string>* = NULL);

//   /**
//    * This method implements writing a mesh with nodal data to a
//    * specified file where the nodal data and variable names are optionally
//    * provided.  This will write an ASCII file as discontinuous patches
//    */
//   virtual void write_discontinuous_ascii (const std::string&,
// 					  const std::vector<Number>* = NULL,
// 					  const std::vector<std::string>* = NULL);

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are optionally
   * provided.  This will write a binary file if the tecio.a library was
   * found at compile time, otherwise a warning message will be printed and
   * an ASCII file will be created.
   */
  virtual void write_binary (const std::string&,
			     const std::vector<Number>* = NULL,
			     const std::vector<std::string>* = NULL);

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
};



// ------------------------------------------------------------
// GMVIO inline members
inline
GMVIO::GMVIO (const MeshBase& mesh) :
  MeshIO<MeshBase> (mesh),
  _binary        (false),
  _discontinuous (false),  
  _partitioning  (true)
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



#endif // #define __gmv_io_h__
