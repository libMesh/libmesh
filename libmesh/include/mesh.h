// $Id: mesh.h,v 1.3 2003-01-20 17:06:11 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_refinement.h"
#include "boundary_info.h"


/**
 * The \p Mesh class is identical to the \p MeshBase class except
 * it contains information about the boundary of the domain and
 * can be refined since it contains a \p MeshRefinement object.
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
  Mesh(unsigned int d,
       unsigned int proc_id=0);

  /**
   * Destructor.
   */
  ~Mesh();
  
  /**
   * This class holds the boundary information.  It can store nodes, edges,
   * and faces with a corresponding id that facilitates setting boundary
   * conditions.
   */
  BoundaryInfo boundary_info;
  
#ifdef ENABLE_AMR

  /**
   * Class that handles adaptive mesh refinement implementation.
   */  
  MeshRefinement mesh_refinement;
  
#endif
  
  /**
   * Builds a \f$ nx \times ny \times nz \f$ (elements) cube.
   * Defaults to a unit cube (or line in 1D, square in 2D),
   * but the dimensions can be specified through the optional
   * arguments.
   */  
  void build_cube (const unsigned int nx=0,
		   const unsigned int ny=0,
		   const unsigned int nz=0,
		   const real xmin=0., const real xmax=1.,
		   const real ymin=0., const real ymax=1.,
		   const real zmin=0., const real zmax=1.,
		   const ElemType type=INVALID_ELEM);


  /**
   * Meshes a spherical or mapped-spherical domain.
   */
  void build_sphere (const real rad=1,
		     const unsigned int nr=2,
		     const ElemType type=INVALID_ELEM);
  
  /**
   * Reads the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.
   */
  void read(const std::string& name);

  
  /**
   * Open the file named \p name and read the mesh in Sandia National Lab's
   * ExodusII format. This is the method to use for reading in meshes generated
   * by cubit.  Works in 2D for \p TRIs, \p TRI6s, \p QUAD s, and \p QUAD9s.
   * Works in 3D for \p TET4s, \p TET10s, \p HEX8s, and \p HEX27s.  This only
   * works for \p Mesh since it needs access to the \p boundary_info
   * structure.
   */
  void read_exd(const std::string& name);

  /**
   * Read meshes in mgflo's XDR format.
   * Should be compatible with the mgf
   * mesh file formats.  This method actualy
   * expects an ASCII-file.
   */
  void read_xdr(const std::string& name);

  /**
   * Same, but expects a true XDR-Encoded binary file.
   */
  void read_xdr_binary(const std::string& name);

  /**
   * Read solutions in mgflo's XDR format.
   * Should be compatible with the mgf
   * solution file format. This method expects
   * an ASCII file.
   */
  void read_xdr_soln(const std::string& name,
		     std::vector<number>& soln,
		     std::vector<std::string>& var_names);

  /**
   * Same, but expects a true XDR-Encoded binary file.
   */
  void read_xdr_soln_binary(const std::string& name,
			    std::vector<number>& soln,
			    std::vector<std::string>& var_names);
  
  /** 
   * Read mesh from the file specified by \p name in Universal (unv) format.  
   */
  void read_unv(const std::string& name);

  /**
   * Write the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.
   */
  void write(const std::string& name);
  
  /**
   * Write to the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension. Also writes data.
   */
  void write(const std::string& name,
	     std::vector<number>& values,
	     std::vector<std::string>& variable_names);
  
  /**
   * Write meshes in DIVA's ASCII
   * format for visualization.
   */
  void write_diva(const std::string& name);
    
  /**
   * Write meshes in mgflo's XDR format.
   * Note: MGF won't be able to read these
   * meshes in general since they will be
   * hybrid meshes.
   */
  void write_xdr(const std::string& name);

  /**
   * Write meshes in mgflo's binary XDR format.
   * Note: MGF won't be able to read these
   * meshes in general since they will be
   * hybrid meshes.
   */
  void write_xdr_binary(const std::string& name);

  /**
   * Write solutions in mgflo's XDR format.
   * Should be compatible with the mgf
   * solution file format.  Writes an ASCII file.
   */
  void write_xdr_soln(const std::string& name,
		      std::vector<number>& soln,
		      std::vector<std::string>& var_names);

  /**
   * Same, but writes an XDR-Encoded binary file.
   */
  void write_xdr_soln_binary(const std::string& name,
			     std::vector<number>& soln,
			     std::vector<std::string>& var_names);

  /**
   * Clear all internal data.
   */
  void clear();


 protected:

  /** 
   * Actual implementation of reading a mesh in Universal (unv) format from
   * a stream.
   */
  void read_unv(std::istream& in);

  /**
   * Actual implementation of writing  meshes in DIVA's ASCII
   * format.
   */
  void write_diva(std::ostream& out);

  
 private:

  
#ifdef ENABLE_AMR

  /**
   * After coarsening a mesh it is possible that
   * there are some voids in the \p nodes and
   * \p elements vectors that need to be cleaned
   * up.  This functions does that.  The indices
   * of the entities to be removed are contained
   * in the two input sets.
   */
  void trim_unused_elements(std::set<unsigned int>& unused_elements);

  friend class MeshRefinement;
  
#endif
  
};




#endif
