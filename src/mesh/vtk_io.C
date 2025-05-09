// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/vtk_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/enum_io_package.h"

#ifdef LIBMESH_HAVE_VTK

// I get a lot of "warning: extra ';' inside a class [-Wextra-semi]" from clang
// on VTK header files.
#include "libmesh/ignore_warnings.h"

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLPUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"

#ifdef LIBMESH_HAVE_MPI
#include "vtkMPI.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIController.h"
#endif

#include "libmesh/restore_warnings.h"

// C++ includes
#include <fstream>


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

// Initialize the static _element_maps map.
std::map<ElemMappingType, VTKIO::ElementMaps> VTKIO::_element_maps = VTKIO::build_element_maps();

// Static function which constructs the ElementMaps object.
std::map<ElemMappingType, VTKIO::ElementMaps> VTKIO::build_element_maps()
{
  // Object to be filled up
  std::map<ElemMappingType, VTKIO::ElementMaps> all_maps;
  ElementMaps em; // Lagrange element maps

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
  all_maps[ElemMappingType::LAGRANGE_MAP] = em;


  // VTK_BEZIER_* types were introduced in VTK 9.0
#if VTK_VERSION_LESS_THAN(9,0,0)
  // Revert back to previous behavior when using an older version of VTK
  all_maps[ElemMappingType::RATIONAL_BERNSTEIN_MAP] = em;
#else
  ElementMaps bem; // Rational Bernstein element maps
  bem.associate(EDGE2, VTK_LINE);
  bem.associate(EDGE3, VTK_BEZIER_CURVE);
  bem.associate(TRI3, VTK_TRIANGLE);
  bem.associate(TRI6, VTK_BEZIER_TRIANGLE);
  bem.associate(QUAD4, VTK_QUAD);
  bem.associate(QUAD8, VTK_QUADRATIC_QUAD);
  bem.associate(QUAD9, VTK_BEZIER_QUADRILATERAL);
  bem.associate(TET4, VTK_TETRA);
  bem.associate(TET10, VTK_QUADRATIC_TETRA);
  bem.associate(HEX8, VTK_HEXAHEDRON);
  bem.associate(HEX20, VTK_QUADRATIC_HEXAHEDRON);
  bem.associate(HEX27, VTK_BEZIER_HEXAHEDRON);
  bem.associate(PRISM6, VTK_WEDGE);
  bem.associate(PRISM15, VTK_QUADRATIC_WEDGE);
  bem.associate(PRISM18, VTK_BEZIER_WEDGE);
  bem.associate(PYRAMID5, VTK_PYRAMID);
  bem.writing_map[TRI3SUBDIVISION] = VTK_TRIANGLE;
  all_maps[ElemMappingType::RATIONAL_BERNSTEIN_MAP] = bem;
#endif

  return all_maps;
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
  typedef vtkSmartPointer<vtkXMLPUnstructuredGridReader> MyReader;
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

  // Try to preserve any libMesh ids and subdomain ids we find in the
  // file.  This will be null if there are none, e.g. if a non-libMesh
  // program wrote this file.

  vtkAbstractArray * abstract_elem_id =
    _vtk_grid->GetCellData()->GetAbstractArray("libmesh_elem_id");
  vtkAbstractArray * abstract_node_id =
    _vtk_grid->GetPointData()->GetAbstractArray("libmesh_node_id");
  vtkAbstractArray * abstract_subdomain_id =
    _vtk_grid->GetCellData()->GetAbstractArray("subdomain_id");

  // Get ids as integers.  This will be null if they are another data
  // type, e.g. if a non-libMesh program used the names we thought
  // were unique for different data.
  vtkIntArray * elem_id = vtkIntArray::SafeDownCast(abstract_elem_id);
  vtkIntArray * node_id = vtkIntArray::SafeDownCast(abstract_node_id);
  vtkIntArray * subdomain_id = vtkIntArray::SafeDownCast(abstract_subdomain_id);

  if (abstract_elem_id && !elem_id)
    libmesh_warning("Found non-integral libmesh_elem_id array; forced to ignore it.\n"
                    "This is technically valid but probably broken.");

  if (abstract_node_id && !node_id)
    libmesh_warning("Found non-integral libmesh_node_id array; forced to ignore it.\n"
                    "This is technically valid but probably broken.");

  if (abstract_subdomain_id && !subdomain_id)
    libmesh_warning("Found non-integral subdomain_id array; forced to ignore it.\n"
                    "This is technically valid but probably broken.");

  // Get the number of points from the _vtk_grid object
  const unsigned int vtk_num_points = static_cast<unsigned int>(_vtk_grid->GetNumberOfPoints());

  // Map from VTK indexing to libMesh id if necessary
  std::vector<dof_id_type> vtk_node_to_libmesh;
  if (node_id)
    vtk_node_to_libmesh.resize(vtk_num_points);

  // always numbered nicely so we can loop like this
  for (unsigned int i=0; i<vtk_num_points; ++i)
    {
      // add to the id map
      // and add the actual point
      double pnt[3];
      _vtk_grid->GetPoint(static_cast<vtkIdType>(i), pnt);
      Point xyz(pnt[0], pnt[1], pnt[2]);

      if (node_id)
        {
          auto id = node_id->GetValue(i);

          // It would nice to distinguish between "duplicate nodes
          // because one was ghosted in a parallel file segment" and
          // "duplicate nodes because there was a bug", but I'm not
          // sure how to do that with vtkXMLPUnstructuredGridReader
          if (!mesh.query_node_ptr(id))
            mesh.add_point(xyz, id);
          vtk_node_to_libmesh[i] = id;
        }
      else
        mesh.add_point(xyz, i);
    }

  // Get the number of cells from the _vtk_grid object
  const unsigned int vtk_num_cells = static_cast<unsigned int>(_vtk_grid->GetNumberOfCells());

  auto& element_map = libmesh_map_find(_element_maps, mesh.default_mapping_type());

  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  for (unsigned int i=0; i<vtk_num_cells; ++i)
    {
      _vtk_grid->GetCell(i, cell);

      // Get the libMesh element type corresponding to this VTK element type.
      ElemType libmesh_elem_type = element_map.find(cell->GetCellType());
      auto elem = Elem::build(libmesh_elem_type);

      // get the straightforward numbering from the VTK cells
      for (auto j : elem->node_index_range())
        {
          const auto vtk_point_id = cell->GetPointId(j);
          const dof_id_type libmesh_node_id = node_id ?
            vtk_node_to_libmesh[vtk_point_id] : vtk_point_id;

          elem->set_node(j, mesh.node_ptr(libmesh_node_id));
        }

      // then get the connectivity
      std::vector<dof_id_type> conn;
      elem->connectivity(0, VTK, conn);

      // then reshuffle the nodes according to the connectivity, this
      // two-time-assign would evade the definition of the vtk_mapping
      for (unsigned int j=0,
           n_conn = cast_int<unsigned int>(conn.size());
           j != n_conn; ++j)
        elem->set_node(j, mesh.node_ptr(conn[j]));

      if (elem_id)
        {
          auto id = elem_id->GetValue(i);
          libmesh_error_msg_if
            (mesh.query_elem_ptr(id), "Duplicate element id " << id <<
             " found in libmesh_elem_ids");
          elem->set_id(id);
        }
      else
        elem->set_id(i);

      if (subdomain_id)
        {
          auto sbdid = subdomain_id->GetValue(i);
          elem->subdomain_id() = sbdid;
        }

      elems_of_dimension[elem->dim()] = true;

      mesh.add_elem(std::move(elem));
    } // end loop over VTK cells

  // Set the mesh dimension to the largest encountered for an element
  for (unsigned char i=0; i!=4; ++i)
    if (elems_of_dimension[i])
      mesh.set_mesh_dimension(i);

#if LIBMESH_DIM < 3
  libmesh_error_msg_if(mesh.mesh_dimension() > LIBMESH_DIM,
                       "Cannot open dimension "
                       << mesh.mesh_dimension()
                       << " mesh file when configured without "
                       << mesh.mesh_dimension()
                       << "D support.");
#endif // LIBMESH_DIM < 3
}



void VTKIO::write_nodal_data (const std::string & fname,
                              const std::vector<Number> & soln,
                              const std::vector<std::string> & names)
{
  // Warn that the .pvtu file extension should be used.  Paraview
  // recognizes this, and it works in both serial and parallel.  Only
  // warn about this once.
  if (fname.substr(fname.rfind("."), fname.size()) != ".pvtu")
    libmesh_do_once(libMesh::err << "The .pvtu extension should be used when writing VTK files in libMesh.");

  // If there are variable names being written, the solution vector
  // should not be empty, it should have been broadcast to all
  // processors by the MeshOutput base class, since VTK is a parallel
  // format.  Verify this before going further.
  libmesh_error_msg_if(!names.empty() && soln.empty(),
                       "Empty soln vector in VTKIO::write_nodal_data().");

  // Get a reference to the mesh
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  // we only use Unstructured grids
  _vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
#ifdef LIBMESH_HAVE_MPI
  // Set VTK to the same communicator as libMesh
  vtkSmartPointer<vtkMPICommunicator> vtk_comm = vtkSmartPointer<vtkMPICommunicator>::New();
  MPI_Comm mpi_comm = mesh.comm().get();
  vtkMPICommunicatorOpaqueComm vtk_opaque_comm(&mpi_comm);
  vtk_comm->InitializeExternal(&vtk_opaque_comm);

  vtkSmartPointer<vtkMPIController> vtk_mpi_ctrl = vtkSmartPointer<vtkMPIController>::New();
  vtk_mpi_ctrl->SetCommunicator(vtk_comm);

  writer->SetController(vtk_mpi_ctrl);
#endif

  // add nodes to the grid and update _local_node_map
  _local_node_map.clear();
  this->nodes_to_vtk();

  // add cells to the grid
  this->cells_to_vtk();

  // add nodal solutions to the grid, if solutions are given
  if (names.size() > 0)
    {
      std::size_t num_vars = names.size();
      std::vector<Number> local_values;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      std::vector<Real> local_real_values;
#endif

      for (std::size_t variable=0; variable<num_vars; ++variable)
        {
          get_local_node_values(local_values, variable, soln, names);

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
          // write real part
          local_real_values.resize(local_values.size());
          std::transform(local_values.begin(), local_values.end(),
                         local_real_values.begin(),
                         [](Number x) { return x.real(); });
          node_values_to_vtk(names[variable] + "_real", local_real_values);

          // write imaginary part
          local_real_values.resize(local_values.size());
          std::transform(local_values.begin(), local_values.end(),
                         local_real_values.begin(),
                         [](Number x) { return x.imag(); });
          node_values_to_vtk(names[variable] + "_imag", local_real_values);
#else
          node_values_to_vtk(names[variable], local_values);
#endif
        }
    }

  // Tell the writer how many partitions exist and on which processor
  // we are currently
  writer->SetNumberOfPieces(mesh.n_processors());
  writer->SetStartPiece(mesh.processor_id());
  writer->SetEndPiece(mesh.processor_id());

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
      writer->SetDataModeToBinary();
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
  // if this grid is to be used in VTK then the dimension of the points should be 3
  pcoords->SetNumberOfComponents(LIBMESH_DIM);
  pcoords->Allocate(3*mesh.n_local_nodes());
  points->SetNumberOfPoints(mesh.n_local_nodes()); // it seems that it needs this to prevent a segfault

  // SetRationalWeights() was introduced in VTK 9.0
#if !VTK_VERSION_LESS_THAN(9,0,0)
  bool have_weights = false;
  int weight_index = 0;
  vtkSmartPointer<vtkDoubleArray> rational_weights;

  if (mesh.default_mapping_type() == ElemMappingType::RATIONAL_BERNSTEIN_MAP)
  {
    rational_weights = vtkSmartPointer<vtkDoubleArray>::New();
    rational_weights->SetName("RationalWeights");
    rational_weights->SetNumberOfComponents(1);
    weight_index = static_cast<int>(mesh.default_mapping_data());
    have_weights = true;
  }
#endif

  vtkSmartPointer<vtkIntArray> node_id = vtkSmartPointer<vtkIntArray>::New();
  node_id->SetName("libmesh_node_id");
  node_id->SetNumberOfComponents(1);

  unsigned int local_node_counter = 0;

  for (const auto & node_ptr : mesh.local_node_ptr_range())
    {
      const Node & node = *node_ptr;

      double pnt[3] = {0, 0, 0};
      for (unsigned int i=0; i<LIBMESH_DIM; ++i)
        pnt[i] = double(node(i));

      // Fill mapping between global and local node numbers
      _local_node_map[node.id()] = local_node_counter;

      // add point
#if VTK_VERSION_LESS_THAN(7,1,0)
      pcoords->InsertNextTupleValue(pnt);
#else
      pcoords->InsertNextTuple(pnt);
#endif
#if !VTK_VERSION_LESS_THAN(9,0,0)
      if (have_weights)
      {
        Real weight = node.get_extra_datum<Real>(weight_index);
        rational_weights->InsertTuple1(local_node_counter, double(weight));
      }
#endif

      node_id->InsertTuple1(local_node_counter, node.id());

      ++local_node_counter;
    }

  // add coordinates to points
  points->SetData(pcoords);

  // add points to grid
  _vtk_grid->SetPoints(points);
  _vtk_grid->GetPointData()->AddArray(node_id);

#if !VTK_VERSION_LESS_THAN(9,0,0)
  if (have_weights)
    _vtk_grid->GetPointData()->SetRationalWeights(rational_weights);

#endif
}



void VTKIO::cells_to_vtk()
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();

  auto& element_map = libmesh_map_find(_element_maps, mesh.default_mapping_type());
  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();

  std::vector<int> types(mesh.n_active_local_elem());

  // We already created this but we need to add more if we have any
  // ghost nodes
  vtkAbstractArray * abstract_node_id =
    _vtk_grid->GetPointData()->GetAbstractArray("libmesh_node_id");
  vtkIntArray * node_id = vtkIntArray::SafeDownCast(abstract_node_id);
  libmesh_assert(node_id);

  vtkSmartPointer<vtkIntArray> elem_id = vtkSmartPointer<vtkIntArray>::New();
  elem_id->SetName("libmesh_elem_id");
  elem_id->SetNumberOfComponents(1);

  vtkSmartPointer<vtkIntArray> subdomain_id = vtkSmartPointer<vtkIntArray>::New();
  subdomain_id->SetName("subdomain_id");
  subdomain_id->SetNumberOfComponents(1);

  vtkSmartPointer<vtkIntArray> elem_proc_id = vtkSmartPointer<vtkIntArray>::New();
  elem_proc_id->SetName("processor_id");
  elem_proc_id->SetNumberOfComponents(1);

  unsigned active_element_counter = 0;
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // When using rational bernstein these hold the weights
      if ( elem->type() == NODEELEM )
        continue;

      pts->SetNumberOfIds(elem->n_nodes());

      // get the connectivity for this element
      std::vector<dof_id_type> conn;
      elem->connectivity(0, VTK, conn);

      for (unsigned int i=0,
           n_conn = cast_int<unsigned int>(conn.size());
           i != n_conn; ++i)
        {
          // If the node ID is not found in the _local_node_map, we'll
          // add it to the _vtk_grid.  NOTE[JWP]: none of the examples
          // I have actually enters this section of code...
          if (!_local_node_map.count(conn[i]))
            {
              dof_id_type global_node_id = elem->node_id(i);

              const Point & the_node = mesh.point(global_node_id);

              // InsertNextPoint accepts either a double or float array of length 3.
              double pt[3] = {0., 0., 0.};
              for (unsigned int d=0; d<LIBMESH_DIM; ++d)
                pt[d] = double(the_node(d));

              // Insert the point into the _vtk_grid
              vtkIdType local = _vtk_grid->GetPoints()->InsertNextPoint(pt);

              // Update the _local_node_map with the ID returned by VTK
              _local_node_map[global_node_id] =
                cast_int<dof_id_type>(local);

              node_id->InsertTuple1(local, global_node_id);
            }

          // Otherwise, the node ID was found in the _local_node_map, so
          // insert it into the vtkIdList.
          pts->InsertId(i, _local_node_map[conn[i]]);
        }

      vtkIdType vtkcellid = cells->InsertNextCell(pts);
      types[active_element_counter] = cast_int<int>(element_map.find(elem->type()));

      elem_id->InsertTuple1(vtkcellid, elem->id());
      subdomain_id->InsertTuple1(vtkcellid, elem->subdomain_id());
      elem_proc_id->InsertTuple1(vtkcellid, elem->processor_id());
      ++active_element_counter;
    } // end loop over active elements

  _vtk_grid->SetCells(types.data(), cells);
  _vtk_grid->GetCellData()->AddArray(elem_id);
  _vtk_grid->GetCellData()->AddArray(subdomain_id);
  _vtk_grid->GetCellData()->AddArray(elem_proc_id);
}

void VTKIO::node_values_to_vtk(const std::string & name,
                               const std::vector<Real> & local_values)
{
  vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
  data->SetName(name.c_str());

  libmesh_assert_equal_to(_local_node_map.size(), local_values.size());

  // number of local and ghost nodes
  data->SetNumberOfValues(_local_node_map.size());

  // copy values into vtk
  for (auto i : index_range(local_values)) {
    data->SetValue(i, double(local_values[i]));
  }

  _vtk_grid->GetPointData()->AddArray(data);
}

void VTKIO::get_local_node_values(std::vector<Number> & local_values,
                                  std::size_t variable,
                                  const std::vector<Number> & soln,
                                  const std::vector<std::string> & names)
{
  const MeshBase & mesh = MeshOutput<MeshBase>::mesh();
  std::size_t num_vars = names.size();
  dof_id_type num_nodes = mesh.n_nodes();

  local_values.clear();
  local_values.resize(_local_node_map.size(), 0.0);

  // loop over all nodes and get the solution for the current
  // variable, if the node is in the current partition
  for (dof_id_type k=0; k<num_nodes; ++k)
    if (const auto local_node_it = _local_node_map.find(k);
        local_node_it != _local_node_map.end())
      local_values[local_node_it->second] = soln[k*num_vars + variable];
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
//       std::map<std::string, std::vector<Number>> vecs;
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
//       std::map<std::string, std::vector<Number>>::iterator it = vecs.begin();
//
//       for (; it!=vecs.end(); ++it)
//         {
//           vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
//           data->SetName(it->first.c_str());
//           libmesh_assert_equal_to (it->second.size(), es.get_mesh().n_nodes());
//           data->SetNumberOfValues(it->second.size());
//
//           for (auto i : index_range(it->second))
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
