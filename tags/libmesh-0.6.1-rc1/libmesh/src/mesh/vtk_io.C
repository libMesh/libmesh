// $Id: vtk_io.C,v 1.6 2007-10-08 21:26:06 woutruijter Exp $

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
#include "mesh.h"
#include "equation_systems.h"
#include "cell_tet4.h"
#include "cell_tet10.h"
#include "cell_prism6.h"
#include "cell_pyramid5.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "petsc_vector.h"
//#include "cell_inf.h"
//#include "cell_inf_prism12.h"
//#include "cell_prism6.h"
//#include "cell_tet4.h"

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
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGenericGeometryFilter.h"
#include "vtkCellArray.h"
#include "vtkConfigure.h"
/*
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
#include "vtkBiQuadraticQuad.h"
*/
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#endif //HAVE_VTK

//static unsigned int vtk_tet4_mapping[4]= {0,1,2,3};
//static unsigned int vtk_pyramid5_mapping[5]= {0,3,2,1,4};
//static unsigned int vtk_wedge6_mapping[6]= {0,2,1,3,5,4};
//static unsigned int vtk_hex8_mapping[8]= {0,1,2,3,4,5,6,7};

#ifdef HAVE_VTK // private functions
void VTKIO::nodes_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid*& grid){
  vtkPoints* points = vtkPoints::New();
// FIXME change to local iterators when I've figured out how to write in parallel
//  MeshBase::const_node_iterator nd = mesh.local_nodes_begin();
//  MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();
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
  grid->SetPoints(points);
}

void VTKIO::cells_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid*& grid){
//  MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
//  const MeshBase::const_element_iterator end = mesh.active_local_elements_end(); 
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
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
	case QUAD9:      
	  celltype = VTK_BIQUADRATIC_QUAD;
	  break;
#else
	case QUAD9:
#endif 
	case EDGE4:      
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
	  for(unsigned int i=0;i<conn.size();++i)
	  {
		  pts->SetId(i,conn[i]);
	  } 
	  grid->InsertNextCell(celltype,pts);
	} // end loop over active elements
}

/*
 * single processor implementation, this can be done in parallel, but requires a
 * bit more thinking on how to deal with overlapping local solution vectors
 */
void VTKIO::solution_to_vtk(const EquationSystems& es,vtkUnstructuredGrid*& grid){
	const MeshBase& mesh = (MeshBase&)es.get_mesh();
	const unsigned int n_nodes = es.get_mesh().n_nodes();
	const unsigned int n_vars = es.n_vars();
	std::vector<Real> soln;
	es.build_solution_vector(soln);

//   if(libMesh::processor_id()==0){
		// write the solutions belonging to the system
		for(unsigned int i=0;i<es.n_systems();++i){ // for all systems, regardless of whether they are active or not
			for(unsigned int j=0;j<n_vars;++j){
				const System& sys = es.get_system(i);
				const std::string& varname = sys.variable_name(j); 
				vtkFloatArray *data = vtkFloatArray::New(); 
				data->SetName(varname.c_str());
				data->SetNumberOfValues(n_nodes);
				MeshBase::const_node_iterator it = mesh.nodes_begin();
				const MeshBase::const_node_iterator n_end = mesh.nodes_end();
				for(;it!=n_end;++it){
					if((*it)->n_comp(i,j)>0){
						const unsigned int nd_id = (*it)->id();
//                  const unsigned int dof_nr = (*it)->dof_number(i,j,0);
//                  data->InsertValue((*it)->id(),sys.current_solution(count++));
						data->InsertValue(nd_id,soln[nd_id*n_vars+j]);
//               data->InsertValue((*it)->id(),sys.current_solution(dof_nr));
					}else{
						data->InsertValue((*it)->id(),0);
					}
				} 
				grid->GetPointData()->AddArray(data);
			} 
		} 
//   }
}
/*
 * FIXME now this is known to write nonsense on AMR meshes
 */
void VTKIO::system_vectors_to_vtk(const EquationSystems& es,vtkUnstructuredGrid*& grid){
	// write out the vectors added to the systems
	const MeshBase& mesh = (MeshBase&)es.get_mesh();
	const unsigned int n_nodes = mesh.n_nodes();
	for(unsigned int i=0;i<es.n_systems();++i){ // for all systems, regardless of whether they are active or not
		const System& sys = es.get_system(i);
		System::const_vectors_iterator v_end = sys.vectors_end();
		System::const_vectors_iterator it = sys.vectors_begin();
		for(;it!= v_end;++it){ // for all vectors on this system
			vtkFloatArray *data = vtkFloatArray::New(); 
			data->SetName(it->first.c_str());
			std::vector<Real> values; 	
			it->second->localize(values);
			data->SetNumberOfValues(n_nodes);
//         MeshBase::const_node_iterator it = mesh.active_nodes_begin();
//         const MeshBase::const_node_iterator n_end = mesh.nodes_end();
//         for(unsigned int count=0;it!=n_end;++it,++count){			
			for(unsigned int j=0;j<values.size();++j){
				data->InsertValue(j,values[j]);
//            it++;
			} 
			grid->GetPointData()->AddArray(data);		
		} 
	} 
}

/*
// write out mesh data to the VTK file, this might come in handy to display
// boundary conditions and material data
inline void meshdata_to_vtk(const MeshData& meshdata,
										vtkUnstructuredGrid* grid){
	vtkPointData* pointdata = vtkPointData::New();
	
	const unsigned int n_vn = meshdata.n_val_per_node();
	const unsigned int n_dat = meshdata.n_node_data();

	pointdata->SetNumberOfTuples(n_dat);
}*/

#endif


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
//  vtkUnstructuredGrid *grid = reader->GetOutput();
  _vtk_grid = reader->GetOutput();
  _vtk_grid->Update();

  // Get a reference to the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  // Clear out any pre-existing data from the Mesh
  mesh.clear();
  
  // read in the nodes
  vtkPoints *points = _vtk_grid->GetPoints();
	
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
	
  for (unsigned int i=0; i < static_cast<unsigned int>(_vtk_grid->GetNumberOfCells()); ++i)
    {
      vtkCell* cell = _vtk_grid->GetCell(i);
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
	case VTK_QUADRATIC_HEXAHEDRON:
  	  elem = new Hex20();
	  break;
	case VTK_QUADRATIC_TETRA:
	  elem = new Tet10();
	  break;
	default:
	  std::cerr << "element type not implemented in vtkinterface " << cell->GetCellType() << std::endl;
	  error();
	}
  // get the straightforward numbering from the VTK cells
  for(unsigned int j=0;j<elem->n_nodes();++j){
	  elem->set_node(j) = mesh.node_ptr(cell->GetPointId(j));
  } 
  // then get the connectivity 
  std::vector<unsigned int> conn;
  elem->connectivity(0,VTK,conn);
  // then reshuffle the nodes according to the connectivity, this
  // two-time-assign would evade the definition of the vtk_mapping
  for(unsigned int j=0;j<conn.size();++j){
	  elem->set_node(j) = mesh.node_ptr(conn[j]);
  } 
  mesh.add_elem(elem);
  } // end loop over VTK cells
  
#endif // HAVE_VTK
}

/*
 * FIXME this operates on the mesh it "gets" from the ES only, this would
 * prevent passing in a mesh that doesn't belong to the ES
 */
void VTKIO::write_equation_systems(const std::string& fname, const EquationSystems& es){
#ifndef HAVE_VTK
	std::cerr << "Cannot write VTK file: " << fname
	    << "\nYou must have VTK installed and correctly configured to read VTK meshes."
	    << std::endl;
	error();
#else
	/*
	 * we only use Unstructured grids
	 */
	_vtk_grid = vtkUnstructuredGrid::New();
	vtkXMLPUnstructuredGridWriter* writer= vtkXMLPUnstructuredGridWriter::New();
	nodes_to_vtk((const MeshBase&)es.get_mesh(), _vtk_grid);
	cells_to_vtk((const MeshBase&)es.get_mesh(), _vtk_grid);

	// I'd like to write out meshdata, but this requires some coding, in
	// particular, non_const meshdata iterators
//   const MeshData& md = es.get_mesh_data();
//   if(es.has_mesh_data())
//      meshdata_to_vtk(md,_vtk_grid);
//   assert (soln.size() ==mesh.n_nodes()*names.size());
	solution_to_vtk(es,_vtk_grid);

#ifdef DEBUG
	if(true) // add some condition here, although maybe it is more sensible to give each vector a flag on whether it is to be written out or not
		system_vectors_to_vtk(es,_vtk_grid);
#endif
//   writer->SetNumberOfPieces(libMesh::n_processors());
//   writer->SetInput(libMesh::processor_id(),_vtk_grid);
	writer->SetInput(_vtk_grid);
	writer->SetFileName(fname.c_str());
	writer->Write();
#endif
}

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
  _vtk_grid = vtkUnstructuredGrid::New();
  vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
  nodes_to_vtk(mesh, _vtk_grid);
  cells_to_vtk(mesh, _vtk_grid);
  writer->SetInput(_vtk_grid);
  writer->SetFileName(name.c_str());
  writer->Write();
#endif // HAVE_VTK
}

//  vim: sw=3 ts=3  
