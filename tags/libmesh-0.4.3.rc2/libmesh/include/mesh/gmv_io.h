// $Id: gmv_io.h,v 1.4 2004-08-05 20:20:09 jwpeterson Exp $

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

// C++ includes
#include <string.h>  // for memcpy

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
   * provided.  
   */
  virtual void write_binary (const std::string&,
			     const std::vector<Number>* = NULL,
			     const std::vector<std::string>* = NULL);

  /**
   * Helper function for writing unsigned ints to an ostream in binary format.
   * Implemented via memcpy as suggested in the standard.
   */
  template <typename T>
  void to_binary_stream(std::ostream& out,
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


template <typename T>
void GMVIO::to_binary_stream(std::ostream& out,
			     const T i)
{
  static char buf[sizeof(T)];
  memcpy(buf, &i, sizeof(T));
  out.write(buf, sizeof(T));
}


#endif // #define __gmv_io_h__
