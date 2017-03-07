// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_config.h"
#include "libmesh/vtk_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"

#ifdef LIBMESH_HAVE_VTK

// I get a lot of "warning: extra ';' inside a class [-Wextra-semi]" from clang
// on VTK header files.
#include "libmesh/ignore_warnings.h"

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkConfigure.h"
#include "vtkDoubleArray.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"

#include "libmesh/restore_warnings.h"

// A convenient macro for comparing VTK versions.  Returns 1 if the
// current VTK version is < major.minor.subminor and zero otherwise.
//
// It relies on the VTK version numbers detected during configure.  Note that if
// LIBMESH_HAVE_VTK is not defined, none of the LIBMESH_DETECTED_VTK_VERSION_* variables will
// be defined either.
#define VTK_VERSION_LESS_THAN(major,minor,subminor)                     \
  ((LIBMESH_DETECTED_VTK_VERSION_MAJOR < (major) ||                     \
    (LIBMESH_DETECTED_VTK_VERSION_MAJOR == (major) && (LIBMESH_DETECTED_VTK_VERSION_MINOR < (minor) || \
                                                       (LIBMESH_DETECTED_VTK_VERSION_MINOR == (minor) && \
                                                        LIBMESH_DETECTED_VTK_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)

#endif // LIBMESH_HAVE_VTK



namespace libMesh
{

// Constructor for reading
VTKIO::VTKIO (MeshBase & mesh) :
  MeshInput<MeshBase> (mesh, /*is_parallel_format=*/true),
  MeshOutput<MeshBase>(mesh, /*is_parallel_format=*/true)
#ifdef LIBMESH_HAVE_VTK
  ,_compress(false)
#endif
{
}



// Constructor for writing
VTKIO::VTKIO (const MeshBase & mesh) :
  MeshOutput<MeshBase>(mesh, /*is_parallel_format=*/true)
#ifdef LIBMESH_HAVE_VTK
  ,_compress(false)
#endif
{
}



// Output the mesh without solutions to a .pvtu file
void VTKIO::write (const std::string & name)
{
  std::vector<Number> soln;
  std::vector<std::string> names;
  this->write_nodal_data(name, soln, names);
}



// The rest of the file is wrapped in ifdef LIBMESH_HAVE_VTK except for
// a couple of "stub" functions at the bottom.
#ifdef LIBMESH_HAVE_VTK

// Initialize the static _element_maps struct.
VTKIO::ElementMaps VTKIO::_element_maps = VTKIO::build_element_maps();

// Static function which constructs the ElementMaps object.
VTKIO::ElementMaps VTKIO::build_element_maps()
{
  // Object to be filled up
  ElementMaps em;

  em.associate(EDGE2, VTK_LINE);
  em.associate(EDGE3, VTK_QUADRATIC_EDGE);
  em.associate(TRI3, VTK_TRIANGLE);
  em.associate(TRI6, VTK_QUADRATIC_TRIANGLE);
  em.associate(QUAD4, VTK_QUAD);
  em.associate(QUAD8, VTK_QUADRATIC_QUAD);
  em.associate(TET4, VTK_TETRA);
  em.associate(TET10, VTK_QUADRATIC_TETRA);
  em.associate(HEX8, VTK_HEXAHEDRON);
  em.associate(HEX20, VTK_QUADRATIC_HEXAHEDRON);
  em.associate(HEX27, VTK_TRIQUADRATIC_HEXAHEDRON);
  em.associate(PRISM6, VTK_WEDGE);
  em.associate(PRISM15, VTK_QUADRATIC_WEDGE);
  em.associate(PRISM18, VTK_BIQUADRATIC_QUADRATIC_WEDGE);
  em.associate(PYRAMID5, VTK_PYRAMID);

  // VTK_BIQUADRATIC_QUAD has been around since VTK 5.0
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
  em.associate(QUAD9, VTK_BIQUADRATIC_QUAD);
#endif

  // TRI3SUBDIVISION is for writing only
  em.writing_map[TRI3SUBDIVISION] = VTK_TRIANGLE;

  return em;
}



void VTKIO::read (const std::string & name)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  // Keep track of what kinds of elements this file contains
  elems_of_dimension.clear();
  elems_of_dimension.resize(4, false);

  // Use a typedef, because these names are just crazy
  typedef vtkSmartPointer<vtkXMLUnstructuredGridReader> MyReader;
  MyReader reader = MyReader::New();

  // Pass the filename along to the reader
  reader->SetFileName(name.c_str());

  // Force reading
  reader->Update();

  // read in the grid
  _vtk_grid = reader->GetOutput();

  // Get a reference to the mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();

  // Clear out any pre-existing data from the Mesh
  mesh.clear();

  // Get the number of points from the _vtk_grid object
  const unsigned int vtk_num_points = static_cast<unsigned int>(_vtk_grid->GetNumberOfPoints());

  // always numbered nicely so we can loop like this
  for (unsigned int i=0; i<vtk_num_points; ++i)
    {
      // add to the id map
      // and add the actual point
      double pnt[3];
      _vtk_grid->GetPoint(static_cast<vtkIdType>(i), pnt);
      Point xyz(pnt[0], pnt[1], pnt[2]);
      mesh.add_point(xyz, i);
    }

  // Get the number of cells from the _vtk_grid object
  const unsigned int vtk_num_cells = static_cast<unsigned int>(_vtk_grid->GetNumberOfCells());

  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  for (unsigned int i=0; i<vtk_num_cells; ++i)
    {
      _vtk_grid->GetCell(i, cell);

      // Get the libMesh element type corresponding to this VTK element type.
      ElemType libmesh_elem_type = _element_maps.find(cell->GetCellType());
      Elem * elem = Elem::build(libmesh_elem_type).release();

      // get the straightforward numbering from the VTK cells
      for (unsigned int j=0; j<elem->n_nodes(); ++j)
        elem->set_node(j) =
          mesh.node_ptr(cast_int<dof_id_type>(cell->GetPointId(j)));

      // then get the connectivity
      std::vector<dof_id_type> conn;
      elem->connectivity(0, VTK, conn);

      // then reshuffle the nodes according to the connectivity, this
      // two-time-assign would evade the definition of the vtk_mapping
      for (std::size_t j=0; j<conn.size(); ++j)
        elem->set_node(j) = mesh.node_ptr(conn[j]);

      elem->set_id(i);

      elems_of_dimension[elem->dim()] = true;

      mesh.add_elem(elem);
    } // end loop over VTK cells

  // Set the mesh dimension to the largest encountered for an element
  for (unsigned char i=0; i!=4; ++i)
    if (elems_of_dimension[i])
      mesh.set_mesh_dimension(i);

#if LIBMESH_DIM < 3
  if (mesh.mesh_dimension() > LIBMESH_DIM)
    libmesh_error_msg("Cannot open dimension "  \
                      << mesh.mesh_dimension()              \
                      << " mesh file when configured without "  \
                      << mesh.mesh_dimension()                  \
                      << "D support.");
#endif // LIBMESH_DIM < 3
}



void VTKIO::write_nodal_data (const std::string & fname,
                              const std::vector<Number> & soln,
                              const std::vector<std::string> & names)
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // Warn that the .pvtu file extension should be used.  Paraview
  // recognizes this, and it works in both serial and parallel.  Only
  // warn about this once.
  if (fname.substr(fname.rfind("."), fname.size()) != ".pvtu")
    libmesh_do_once(libMesh::err << "The .pvtu extension should be used when writing VTK files in libMesh.");

  // If there are variable names being written, the solution vector
  // should not be empty, it should have been broadcast to all
  // processors by the MeshOutput base class, since VTK is a parallel
  // format.  Verify this before going further.
  if (!names.empty() && soln.empty())
    libmesh_error_msg("Empty soln vector in VTKIO::write_nodal_data().");

  // we only use Unstructured grids
  _vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

  // add nodes to the grid and update _local_node_map
  _local_node_map.clear();
  this->nodes_to_vtk();

  // add cells to the grid
  this->cells_to_vtk();

  // add nodal solutions to the grid, if solutions are given
  if (names.size() > 0)
    {
      std::size_t num_vars = names.size();
      dof_id_type num_nodes = mesh.n_nodes();

      for (std::size_t variable=0; variable<num_vars; ++variable)
        {
          vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
          data->SetName(names[variable].c_str());

          // number of local and ghost nodes
          data->SetNumberOfValues(_local_node_map.size());

          // loop over all nodes and get the solution for the current
          // variable, if the node is in the current partition
          for (dof_id_type k=0; k<num_nodes; ++k)
            {
              std::map<dof_id_type, dof_id_type>::iterator local_node_it = _local_node_map.find(k);
              if (local_node_it == _local_node_map.end())
                continue; // not a local node

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
              libmesh_do_once (libMesh::err << "Only writing the real part for complex numbers!\n"
                               << "if you need this support contact " << LIBMESH_PACKAGE_BUGREPORT
                               << std::endl);
              data->SetValue(local_node_it->second, soln[k*num_vars + variable].real());
#else
              data->SetValue(local_node_it->second, soln[k*num_vars + variable]);
#endif
            }
          _vtk_grid->GetPointData()->AddArray(data);
        }
    }

  // Tell the writer how many partitions exist and on which processor
  // we are currently
  writer->SetNumberOfPieces(MeshOutput<MeshBase>::mesh().n_processors());
  writer->SetStartPiece(MeshOutput<MeshBase>::mesh().processor_id());
  writer->SetEndPiece(MeshOutput<MeshBase>::mesh().processor_id());

  // partitions overlap by one node
  // FIXME: According to this document
  // http://paraview.org/Wiki/images/5/51/SC07_tut107_ParaView_Handouts.pdf
  // the ghosts are cells rather than nodes.
  writer->SetGhostLevel(1);

  // VTK 6 replaces SetInput() with SetInputData(). See
  // http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
  // for the full explanation.
#if VTK_VERSION_LESS_THAN(6,0,0)
  writer->SetInput(_vtk_grid);
#else
  writer->SetInputData(_vtk_grid);
#endif

  writer->SetFileName(fname.c_str());
  writer->SetDataModeToAscii();

  // compress the output, if desired (switches also to binary)
  if (this->_compress)
    {
#if !VTK_VERSION_LESS_THAN(5,6,0)
      writer->SetCompressorTypeToZLib();
#else
      libmesh_do_once(libMesh::err << "Compression not implemented with old VTK libs!" << std::endl;);
#endif
    }

  writer->Write();

}



vtkUnstructuredGrid * VTKIO::get_vtk_grid()
{
  return _vtk_grid;
}



void VTKIO::set_compression(bool b)
{
  this->_compress = b;
}



void VTKIO::nodes_to_vtk()
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // containers for points and coordinates of points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkDoubleArray> pcoords = vtkSmartPointer<vtkDoubleArray>::New();
  // if this grid is to be used in VTK then the dimesion of the points should be 3
  pcoords->SetNumberOfComponents(LIBMESH_DIM);
  pcoords->Allocate(3*mesh.n_local_nodes());
  points->SetNumberOfPoints(mesh.n_local_nodes()); // it seems that it needs this to prevent a segfault

  unsigned int local_node_counter = 0;

  MeshBase::const_node_iterator nd = mesh.local_nodes_begin();
  MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();
  for (; nd != nd_end; nd++, ++local_node_counter)
    {
      Node & node = **nd;

      double pnt[3] = {0, 0, 0};
      for (unsigned int i=0; i<LIBMESH_DIM; ++i)
        pnt[i] = node(i);

      // Fill mapping between global and local node numbers
      _local_node_map[node.id()] = local_node_counter;

      // add point
#if VTK_VERSION_LESS_THAN(7,1,0)
      pcoords->InsertNextTupleValue(pnt);
#else
      pcoords->InsertNextTuple(pnt);
#endif
    }

  // add coordinates to points
  points->SetData(pcoords);

  // add points to grid
  _vtk_grid->SetPoints(points);
}



void VTKIO::cells_to_vtk()
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();

  std::vector<int> types(mesh.n_active_local_elem());
  unsigned active_element_counter = 0;

  vtkSmartPointer<vtkIntArray> elem_id = vtkSmartPointer<vtkIntArray>::New();
  elem_id->SetName("libmesh_elem_id");
  elem_id->SetNumberOfComponents(1);

  vtkSmartPointer<vtkIntArray> subdomain_id = vtkSmartPointer<vtkIntArray>::New();
  subdomain_id->SetName("subdomain_id");
  subdomain_id->SetNumberOfComponents(1);

  vtkSmartPointer<vtkIntArray> elem_proc_id = vtkSmartPointer<vtkIntArray>::New();
  elem_proc_id->SetName("processor_id");
  elem_proc_id->SetNumberOfComponents(1);

  MeshBase::const_element_iterator it = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_local_elements_end();
  for (; it != end; ++it, ++active_element_counter)
    {
      Elem * elem = *it;

      pts->SetNumberOfIds(elem->n_nodes());

      // get the connectivity for this element
      std::vector<dof_id_type> conn;
      elem->connectivity(0, VTK, conn);

      for (std::size_t i=0; i<conn.size(); ++i)
        {
          // If the node ID is not found in the _local_node_map, we'll
          // add it to the _vtk_grid.  NOTE[JWP]: none of the examples
          // I have actually enters this section of code...
          if (_local_node_map.find(conn[i]) == _local_node_map.end())
            {
              dof_id_type global_node_id = elem->node_id(i);

              const Point & the_node = mesh.point(global_node_id);

              // InsertNextPoint accepts either a double or float array of length 3.
              double pt[3] = {0., 0., 0.};
              for (unsigned int d=0; d<LIBMESH_DIM; ++d)
                pt[d] = the_node(d);

              // Insert the point into the _vtk_grid
              vtkIdType local = _vtk_grid->GetPoints()->InsertNextPoint(pt);

              // Update the _local_node_map with the ID returned by VTK
              _local_node_map[global_node_id] =
                cast_int<dof_id_type>(local);
            }

          // Otherwise, the node ID was found in the _local_node_map, so
          // insert it into the vtkIdList.
          pts->InsertId(i, _local_node_map[conn[i]]);
        }

      vtkIdType vtkcellid = cells->InsertNextCell(pts);
      types[active_element_counter] = cast_int<int>(_element_maps.find(elem->type()));

      elem_id->InsertTuple1(vtkcellid, elem->id());
      subdomain_id->InsertTuple1(vtkcellid, elem->subdomain_id());
      elem_proc_id->InsertTuple1(vtkcellid, elem->processor_id());
    } // end loop over active elements

  _vtk_grid->SetCells(&types[0], cells);
  _vtk_grid->GetCellData()->AddArray(elem_id);
  _vtk_grid->GetCellData()->AddArray(subdomain_id);
  _vtk_grid->GetCellData()->AddArray(elem_proc_id);
}



/**
 * FIXME: This is known to write nonsense on AMR meshes
 * and it strips the imaginary parts of complex Numbers
 *
 * This function is not currently used by anything, so it is commented
 * out, and may eventually be removed entirely.
 */
// void VTKIO::system_vectors_to_vtk(const EquationSystems & es,
//                                   vtkUnstructuredGrid *& grid)
// {
//   if (MeshOutput<MeshBase>::mesh().processor_id() == 0)
//     {
//       std::map<std::string, std::vector<Number> > vecs;
//       for (unsigned int i=0; i<es.n_systems(); ++i)
//         {
//           const System & sys = es.get_system(i);
//           System::const_vectors_iterator v_end = sys.vectors_end();
//           System::const_vectors_iterator it = sys.vectors_begin();
//           for (; it!= v_end; ++it)
//             {
//               // for all vectors on this system
//               std::vector<Number> values;
//               // libMesh::out<<"it "<<it->first<<std::endl;
//
//               it->second->localize_to_one(values, 0);
//               // libMesh::out<<"finish localize"<<std::endl;
//               vecs[it->first] = values;
//             }
//         }
//
//       std::map<std::string, std::vector<Number> >::iterator it = vecs.begin();
//
//       for (; it!=vecs.end(); ++it)
//         {
//           vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
//           data->SetName(it->first.c_str());
//           libmesh_assert_equal_to (it->second.size(), es.get_mesh().n_nodes());
//           data->SetNumberOfValues(it->second.size());
//
//           for (std::size_t i=0; i<it->second.size(); ++i)
//             {
// #ifdef LIBMESH_USE_COMPLEX_NUMBERS
//               libmesh_do_once (libMesh::err << "Only writing the real part for complex numbers!\n"
//                                << "if you need this support contact " << LIBMESH_PACKAGE_BUGREPORT
//                                << std::endl);
//               data->SetValue(i, it->second[i].real());
// #else
//               data->SetValue(i, it->second[i]);
// #endif
//
//             }
//           grid->GetPointData()->AddArray(data);
//         }
//     }
// }



#else // !LIBMESH_HAVE_VTK

void VTKIO::read (const std::string & name)
{
  libmesh_error_msg("Cannot read VTK file: " << name \
                    << "\nYou must have VTK installed and correctly configured to read VTK meshes.");
}



void VTKIO::write_nodal_data (const std::string & fname,
                              const std::vector<Number> &,
                              const std::vector<std::string> &)
{
  libmesh_error_msg("Cannot write VTK file: " << fname                  \
                    << "\nYou must have VTK installed and correctly configured to read VTK meshes.");
}


#endif // LIBMESH_HAVE_VTK



} // namespace libMesh
