// $Id: mesh.h,v 1.12 2005-02-22 22:17:33 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __mesh_h__
#define __mesh_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "mesh_base.h"
//#include "elem_iterators.h"



/**
 * The \p Mesh class is derived from the \p MeshBase class.
 * The user will typically want to instantiate and use the
 * Mesh class in her applications.
 * In order to use the adaptive mesh refinment capabilities
 * of the library, first instantiate a MeshRefinement object
 * with a reference to this class.  Then call the appropriate
 * refinement functions from that object.  To interact with the 
 * boundary, instantiate a BoundaryMesh with a reference to
 * this class, and then use that object's functionality.
*/

// ------------------------------------------------------------
// Mesh class definition
class Mesh : public MeshBase
{
 public:

  /**
   * Constructor.  Requires the dimension and optionally
   * a processor id.  Note that \p proc_id should always
   * be provided for multiprocessor applications.
   */
  Mesh (unsigned int d);
  
  /**
   * Destructor.
   */
  ~Mesh();

//   /**
//    * Builds a \f$ nx \times ny \times nz \f$ (elements) cube.
//    * Defaults to a unit cube (or line in 1D, square in 2D),
//    * but the dimensions can be specified through the optional
//    * arguments.
//    */  
//   void build_cube (const unsigned int nx=0,
// 		   const unsigned int ny=0,
// 		   const unsigned int nz=0,
// 		   const Real xmin=0., const Real xmax=1.,
// 		   const Real ymin=0., const Real ymax=1.,
// 		   const Real zmin=0., const Real zmax=1.,
// 		   const ElemType type=INVALID_ELEM);

//   /**
//    * A specialized \p build_cube() for 2D meshes.
//    */
//   void build_square (const unsigned int nx,
// 		     const unsigned int ny,
// 		     const Real xmin=0., const Real xmax=1.,
// 		     const Real ymin=0., const Real ymax=1.,
// 		     const ElemType type=INVALID_ELEM);

//   /**
//    * Meshes a spherical or mapped-spherical domain.
//    */
//   void build_sphere (const Real rad=1,
// 		     const unsigned int nr=2,
// 		     const ElemType type=INVALID_ELEM);

  /**
   * Reads the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.  This is now the only
   * way to read a mesh.  The \p Mesh then initializes its data
   * structures and is ready for use.
   *
   * In order to read the UNV and TetGen file types, you must
   * also pass a separate pointer to the MeshData object you will
   * use with this mesh, since these read methods expect it.
   */
  void read (const std::string& name,
	     MeshData* mesh_data=NULL);
  /**
   * Write the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.
   *
   * In order to write the UNV and TetGen file types, you must
   * also pass a separate pointer to the MeshData object you have been
   * using with this mesh, since these write methods expect it.
   */
  void write (const std::string& name,
	      MeshData* mesh_data=NULL);
  
  /**
   * Write to the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension. Also writes data.
   */
  void write (const std::string& name,
	      const std::vector<Number>& values,
	      const std::vector<std::string>& variable_names);

  /**
   * Clear all internal data.
   */
  void clear();

  /**
   * Converts a (conforming, non-refined) mesh with linear 
   * elements into a mesh with second-order elements.  For 
   * example, a mesh consisting of \p Tet4 will be converted
   * to a mesh with \p Tet10 etc.  Note that for some elements
   * like \p Hex8 there exist @e two higher order equivalents,
   * \p Hex20 and \p Hex27.  When \p full_ordered is \p true
   * (default), then \p Hex27 is built.  Otherwise, \p Hex20
   * is built.  The same holds obviously for \p Quad4, \p Prism6
   * ...
   */
  void all_second_order (const bool full_ordered=true);
  
  /**
   * Generates a new mesh containing all the elements which
   * are assigned to processor \p pid.  This mesh is written
   * to the pid_mesh reference which you must create and pass
   * to the function.
   */
  void create_pid_mesh (Mesh& pid_mesh,
			const unsigned int pid) const;
  
  /**
   * Constructs a mesh called "new_mesh" from the current mesh by
   * iterating over the elements between it and it_end and adding
   * them to the new mesh.
   */
  void create_submesh (Mesh& new_mesh,
		       const_element_iterator& it,
		       const const_element_iterator& it_end) const;


  
 private:
};




#endif
