// $Id: vtk_io.C,v 1.3 2007-09-17 19:27:47 woutruijter Exp $

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


// C++ includes
#include <fstream>

// Local includes
#include "vtk_io.h"
#include "mesh_base.h"
#include "cell_tet4.h"
#include "cell_tet10.h"
#include "cell_prism6.h"
#include "cell_pyramid5.h"
#include "cell_hex8.h"

//#include "cell_inf.h"
//#include "cell_inf_prism12.h"
//#include "cell_prism6.h"
//#include "cell_tet4.h"

//#include "cell_hex20.h"
//#include "cell_inf_hex16.h"
//#include "cell_inf_prism6.h"
//#include "cell_prism.h"
//#include "cell_tet.h"

//#include "cell_hex27.h"
//#include "cell_inf_hex18.h"
//#include "cell_inf_prism.h"
//#include "cell_pyramid5.h"

//#include "cell_hex8.h"
//#include "cell_inf_hex8.h"
//#include "cell_prism15.h"
//#include "cell_pyramid.h"

//#include "cell_hex.h"
//#include "cell_inf_hex.h"
//#include "cell_prism18.h"
//#include "cell_tet10.h"



#include "mesh_data.h"

#ifdef HAVE_VTK

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGenericGeometryFilter.h"
#include "vtkCellArray.h"
#include "vtkTetra.h"
#include "vtkTriangle.h"
#include "vtkWedge.h"
#include "vtkPyramid.h"
#include "vtkQuad.h"
#include "vtkHexahedron.h"
#include "vtkQuadraticTriangle.h"
#include "vtkQuadraticQuad.h"
#include "vtkQuadraticHexahedron.h"
#include "vtkQuadraticWedge.h"
#include "vtkQuadraticTetra.h"

#endif //HAVE_VTK



// ------------------------------------------------------------
// vtkIO class members
//
void VTKIO::read (const std::string& name)
{
#ifndef HAVE_VTK
  std::cerr << "Cannot read VTK file: " << name
	    << "\nYou must have VTK installed and correctly configured to read VTK meshes."
	    << std::endl;
  error();

#else
  //std::cout<<"read "<<name <<std::endl;  
  vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
  reader->SetFileName( name.c_str() );
  //std::cout<<"force read"<<std::endl;  
  // Force reading
  reader->Update();

  // read in the grid   
  vtkUnstructuredGrid *grid = reader->GetOutput();
  grid->Update();

  // Get a reference to the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  // Clear out any pre-existing data from the Mesh
  mesh.clear();
  
  // read in the nodes
  vtkPoints *points = grid->GetPoints();
	
  // always numbered nicely??, so we can loop like this
  for (unsigned int i=0; i < static_cast<unsigned int>(points->GetNumberOfPoints()); ++i)
    {
      // add to the id map
      // and add the actual point
      double * pnt = points->GetPoint(i);
      Point xyz(pnt[0],pnt[1],pnt[2]);
      Node* newnode = mesh.add_point(xyz);
	    
      // Add node to the nodes vector &
      // tell the MeshData object the foreign node id.
      if (this->_mesh_data != NULL)
	this->_mesh_data->add_foreign_node_id (newnode, i);
    }
	
  for (unsigned int i=0; i < static_cast<unsigned int>(grid->GetNumberOfCells()); ++i)
    {
      vtkCell* cell = grid->GetCell(i);
      Elem* elem = NULL;  // Initialize to avoid compiler warning
	    
      switch(cell->GetCellType())
	{
	case VTK_TETRA:
	  elem = new Tet4();
	  break;
		
	case VTK_WEDGE:
	  elem = new Prism6();
	  break;
		
	case VTK_HEXAHEDRON:
	  elem = new Hex8();
	  break;
		
	case VTK_PYRAMID:	
	  elem = new Pyramid5();					
	  break;
		
	default:
	  std::cerr << "element type not implemented in vtkinterface " << cell->GetCellType() << std::endl;
	  error();

	}
	    
      for(unsigned int j=0;j<elem->n_nodes();++j)
	{
	  elem->set_node(j) = mesh.node_ptr(cell->GetPointId(j));
	} 
      mesh.add_elem(elem);
    } // end loop over VTK cells
  
#endif // HAVE_VTK
}



// // a copy from the VTK docs
//    VTKCellType {
//  VTK_EMPTY_CELL = 0,
//  VTK_VERTEX = 1,
//  VTK_POLY_VERTEX = 2,
//  VTK_LINE = 3,
//  VTK_POLY_LINE = 4,
//  VTK_TRIANGLE = 5,
//  VTK_TRIANGLE_STRIP = 6,
//  VTK_POLYGON =7,
//  VTK_PIXEL = 8,
//  VTK_QUAD = 9,
//  VTK_TETRA = 10,
//  VTK_VOXEL = 11,
//  VTK_HEXAHEDRON = 12,
//  VTK_WEDGE = 13,
//  VTK_PYRAMID = 14,
//  VTK_PENTAGONAL_PRISM= 15,
//  VTK_HEXAGONAL_PRISM = 16,
//  VTK_QUADRATIC_EDGE = 21,
//  VTK_QUADRATIC_TRIANGLE =22,
//  VTK_QUADRATIC_QUAD = 23,
//  VTK_QUADRATIC_TETRA = 24,
//  VTK_QUADRATIC_HEXAHEDRON = 25,
//  VTK_QUADRATIC_WEDGE= 26,
//  VTK_QUADRATIC_PYRAMID = 27,
//  VTK_BIQUADRATIC_QUAD = 28,
//  VTK_TRIQUADRATIC_HEXAHEDRON = 29,
//  VTK_QUADRATIC_LINEAR_QUAD = 30,
//  VTK_QUADRATIC_LINEAR_WEDGE = 31,
//  VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32,
//  VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON =33,
//  VTK_CONVEX_POINT_SET = 41,
//  VTK_PARAMETRIC_CURVE = 51,
//  VTK_PARAMETRIC_SURFACE = 52,
//  VTK_PARAMETRIC_TRI_SURFACE = 53,
//  VTK_PARAMETRIC_QUAD_SURFACE = 54,
//  VTK_PARAMETRIC_TETRA_REGION = 55,
//  VTK_PARAMETRIC_HEX_REGION = 56,
//  VTK_HIGHER_ORDER_EDGE = 60,
//  VTK_HIGHER_ORDER_TRIANGLE = 61,
//  VTK_HIGHER_ORDER_QUAD = 62,
//  VTK_HIGHER_ORDER_POLYGON = 63,
//  VTK_HIGHER_ORDER_TETRAHEDRON = 64,
//  VTK_HIGHER_ORDER_WEDGE = 65,
//  VTK_HIGHER_ORDER_PYRAMID = 66,
//  VTK_HIGHER_ORDER_HEXAHEDRON = 67,
//  VTK_NUMBER_OF_CELL_TYPES
// }



/**
 * This method implements writing to a .vtu (VTK Unstructured Grid) file. 
 * This is one of the new style XML dataformats, binary output is used to keep
 * the file size down. 
 */
void VTKIO::write (const std::string& name)
{	
#ifndef HAVE_VTK
  std::cerr << "Cannot write VTK file: " << name
	    << "\nYou must have VTK installed and correctly configured to write VTK meshes."
	    << std::endl;
  error();

#else
  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  vtkPoints* points = vtkPoints::New();
  vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
  vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
  MeshBase::const_node_iterator nd = mesh.active_nodes_begin();
  MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();
  // write out nodal data
  for (;nd!=nd_end;++nd)
    {
      float tuple[3]; //, dsp[3];
      Node* node = (*nd);

      for(unsigned int i=0;i<3;++i)
	{
	  tuple[i] = (*node)(i);
	}
      points->InsertPoint(node->id(),tuple);
    }
	
  // write out element data 
  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
  for ( ; it != end; ++it)
    {
      Elem *elem  = (*it);
      vtkIdType celltype = VTK_EMPTY_CELL; // initialize to something to avoid compiler warning
		
      switch(elem->type())
	{
	case EDGE2:				
	  celltype = VTK_LINE;
	  break;
	case EDGE3:      
	  celltype = VTK_QUADRATIC_EDGE;
	  break;// 1
	case TRI3:       
	  celltype = VTK_TRIANGLE;
	  break;// 3
	case TRI6:       
	  celltype = VTK_QUADRATIC_TRIANGLE;
	  break;// 4
	case QUAD4:      
	  celltype = VTK_QUAD;
	  break;// 5
	case QUAD8:      
	  celltype = VTK_QUADRATIC_QUAD;
	  break;// 6
	case TET4:      
	  celltype = VTK_TETRA;
	  break;// 8
	case TET10:      
	  celltype = VTK_QUADRATIC_TETRA;
	  break;// 9
	case HEX8:    
	  celltype = VTK_HEXAHEDRON;
	  break;// 10
	case HEX20:      
	  celltype = VTK_QUADRATIC_HEXAHEDRON;
	  break;// 12
	case PRISM6:     
	  celltype = VTK_WEDGE;
	  break;// 13
	case PRISM15:   
	  celltype = VTK_HIGHER_ORDER_WEDGE;
	  break;// 14
	case PRISM18:    
	  break;// 15
	case PYRAMID5:
	  celltype = VTK_PYRAMID;
	  break;// 16
	case EDGE4:      
	case QUAD9:      
	case HEX27:      
	case INFEDGE2:   
	case INFQUAD4:   
	case INFQUAD6:   
	case INFHEX8:    
	case INFHEX16:   
	case INFHEX18:   
	case INFPRISM6:  
	case INFPRISM12: 
	case NODEELEM:   
	case INVALID_ELEM:
	default:
	  {
		  std::cerr<<"element type "<<elem->type()<<" not implemented"<<std::endl;
		  error();
	  }
	}

	  vtkIdList *pts = vtkIdList::New();
	  pts->SetNumberOfIds(elem->n_nodes());
	  // get the connectivity for this element
	  std::vector<unsigned int> conn;
	  elem->connectivity(0,VTK,conn);
	  for(unsigned int i=0;i<elem->n_nodes();++i)
	  {
		  pts->SetId(i,conn[i]);
	  } 
	  grid->InsertNextCell(celltype,pts);
	} // end loop over active elements

  grid->SetPoints(points);
  writer->SetInput(grid);
  writer->SetFileName(name.c_str());
  writer->Write();
#endif // HAVE_VTK
}
