// $Id: mesh.h,v 1.25 2003-10-01 04:02:52 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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

  /**
   * Builds a \f$ nx \times ny \times nz \f$ (elements) cube.
   * Defaults to a unit cube (or line in 1D, square in 2D),
   * but the dimensions can be specified through the optional
   * arguments.
   */  
  void build_cube (const unsigned int nx=0,
		   const unsigned int ny=0,
		   const unsigned int nz=0,
		   const Real xmin=0., const Real xmax=1.,
		   const Real ymin=0., const Real ymax=1.,
		   const Real zmin=0., const Real zmax=1.,
		   const ElemType type=INVALID_ELEM);

  /**
   * A specialized \p build_cube() for 2D meshes.
   */
  void build_square (const unsigned int nx,
		     const unsigned int ny,
		     const Real xmin=0., const Real xmax=1.,
		     const Real ymin=0., const Real ymax=1.,
		     const ElemType type=INVALID_ELEM);

  /**
   * Meshes a spherical or mapped-spherical domain.
   */
  void build_sphere (const Real rad=1,
		     const unsigned int nr=2,
		     const ElemType type=INVALID_ELEM);

  /**
   * Reads the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.  This is now the only
   * way to read a mesh.  It allows the library to so some initialization
   * after the mesh is read (particularly for for parallel data layout).
   * This simplifies life by centralizing this code.
   */
  void read (const std::string& name,
	     const bool do_prepare_for_use = true);

  
  /**
   * Read solutions in mgflo's XDR format.
   * Should be compatible with the MGF
   * solution file format. This method expects
   * an ASCII file.  What is MGF?  It was a microgravity
   * fluid physics code developed under a NASA ESS Grand
   * Challenge Grant.  This method exists solely for
   * backwards compatiblity with MGF and could be
   * deprecated at any time.
   */
  void read_xdr_soln (const std::string& name,
		      std::vector<Number>& soln,
		      std::vector<std::string>& var_names);

  /**
   * Same, but expects a true XDR-Encoded binary file.
   */
  void read_xdr_soln_binary (const std::string& name,
			     std::vector<Number>& soln,
			     std::vector<std::string>& var_names);
  
  /**
   * Write the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.
   */
  void write (const std::string& name);
  
  /**
   * Write to the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension. Also writes data.
   */
  void write (const std::string& name,
	      std::vector<Number>& values,
	      std::vector<std::string>& variable_names);
  
  /**
   * Write meshes in DIVA's ASCII
   * format for visualization.
   */
  void write_diva (const std::string& name);
    
  /**
   * Write meshes in mgflo's XDR format.
   * Note: MGF won't be able to read these
   * meshes in general since they will be
   * hybrid meshes.
   */
  void write_xdr (const std::string& name);

  /**
   * Write meshes in mgflo's binary XDR format.
   * Note: MGF won't be able to read these
   * meshes in general since they will be
   * hybrid meshes.
   */
  void write_xdr_binary (const std::string& name);

  /**
   * Write solutions in mgflo's XDR format.
   * Should be compatible with the mgf
   * solution file format.  Writes an ASCII file.
   */
  void write_xdr_soln (const std::string& name,
		       std::vector<Number>& soln,
		       std::vector<std::string>& var_names);
  
  /**
   * Same, but writes an XDR-Encoded binary file.
   */
  void write_xdr_soln_binary (const std::string& name,
			      std::vector<Number>& soln,
			      std::vector<std::string>& var_names);

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
   * document me!
   */
  void create_submesh (Mesh& new_mesh,
		       const_elem_iterator& it,
		       const const_elem_iterator& it_end) const;
 protected:

  /**
   * Open the file named \p name and read the mesh in Sandia National Lab's
   * ExodusII format. This is the method to use for reading in meshes generated
   * by cubit.  Works in 2D for \p TRIs, \p TRI6s, \p QUAD s, and \p QUAD9s.
   * Works in 3D for \p TET4s, \p TET10s, \p HEX8s, and \p HEX27s.  This only
   * works for \p Mesh since it needs access to the \p boundary_info
   * structure.
   */
  void read_exd (const std::string& name);

  /**
   * Read meshes in mgflo's XDR format.
   * Should be compatible with the mgf
   * mesh file formats.  This method actualy
   * expects an ASCII-file.
   */
  void read_xdr (const std::string& name);

  /**
   * Same, but expects a true XDR-Encoded binary file.
   */
  void read_xdr_binary (const std::string& name);

  /**
   * Actual implementation of writing meshes in DIVA's ASCII
   * format.
   */
  void write_diva (std::ostream& out);

  
 private:
};




#endif
