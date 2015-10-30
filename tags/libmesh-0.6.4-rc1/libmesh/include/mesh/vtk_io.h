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


#ifndef __vtk_io_h__
#define __vtk_io_h__

// C++ includes
#include <map>

// Local includes
#include "libmesh_common.h"
#include "mesh_input.h"
#include "mesh_output.h"


// Forward declarations
class MeshBase;
class MeshData;

class vtkUnstructuredGrid;
class vtkPoints;
class vtkCellArray;

/**
 * This class implements reading and writing meshes in the VTK format.
 * Format description: 
 * cf. <a href="http://www.vtk.org/">VTK home page</a>.
 *
 * This class will not have any functionality unless VTK is detected
 * during configure and hence LIBMESH_HAVE_VTK is defined.
 *
 * @author Wout Ruijter, 2007
 * (Checked in to LibMesh by J.W. Peterson)
 */

// ------------------------------------------------------------
// VTKIO class definition
class VTKIO : public MeshInput<MeshBase>,
	      public MeshOutput<MeshBase>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  VTKIO (MeshBase& mesh, MeshData* mesh_data=NULL);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  VTKIO (const MeshBase& mesh, MeshData* mesh_data=NULL);

 /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
//  virtual void write_nodal_data (const std::string&,
//             const std::vector<Number>&,
//             const std::vector<std::string>&);

  /**
   * Overloads writing equation systems, this is done because when overloading
   * write_nodal_data there would be no way to export cell centered data
   */

  virtual void write_equation_systems(const std::string& fname, const EquationSystems& es); 

  /**
   * This method implements reading a mesh from a specified file
   * in VTK format.
   */
  virtual void read (const std::string& );

  /**
   * This method implements writing a mesh to a specified ".poly"   file.
   * ".poly" files defines so called Piecewise Linear Complex   (PLC).
  */
  virtual void write (const std::string& );  

  /** 
   * Get a pointer to the VTK datastructure
   */
   vtkUnstructuredGrid* get_vtk_grid(){return _vtk_grid;}
	
private:

  /**
   * write the nodes from the mesh into a vtkUnstructuredGrid
   */
  vtkPoints* nodes_to_vtk(const MeshBase& mesh);
//  , vtkUnstructuredGrid*& grid);

  /** 
   * write the cells from the mesh into a vtkUnstructuredGrid
   */
  vtkCellArray* cells_to_vtk(const MeshBase& mesh, int*& types);
//  , vtkUnstructuredGrid*& grid);

  /**
   * write the solution to a vtkUnstructuredGrid
   */
  void solution_to_vtk(const EquationSystems& es,vtkUnstructuredGrid*& grid);

  /**
   * write the system vectors to vtk 
   */
  void system_vectors_to_vtk(const EquationSystems& es,vtkUnstructuredGrid*& grid);

   /**
    * pointer to the VTK grid
    */
   vtkUnstructuredGrid* _vtk_grid;
  
  /**
   * A pointer to the MeshData object you would like to use.
   * with this VTKIO object.  Can be NULL.
   */
  MeshData* _mesh_data;
};



// ------------------------------------------------------------
// VTKIO inline members
inline
VTKIO::VTKIO (MeshBase& mesh, MeshData* mesh_data) :
	MeshInput<MeshBase> (mesh),
	MeshOutput<MeshBase>(mesh),
	_mesh_data(mesh_data)
{
  _vtk_grid = NULL;	
  untested();
}



inline
VTKIO::VTKIO (const MeshBase& mesh, MeshData* mesh_data) :
	MeshOutput<MeshBase>(mesh),
	_mesh_data(mesh_data)
{
  _vtk_grid = NULL;	
  untested();
}



#endif // #define __vtk_io_h__
