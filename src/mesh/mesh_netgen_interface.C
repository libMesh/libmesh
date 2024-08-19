// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_NETGEN


// C++ includes
#include <sstream>

// Local includes
#include "libmesh/mesh_netgen_interface.h"

#include "libmesh/boundary_info.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/utility.h" // libmesh_map_find

namespace {

// RAII for exception safety
class WrappedNgMesh
{
public:
  WrappedNgMesh() {
    _ngmesh = nglib::Ng_NewMesh();
  }

  ~WrappedNgMesh() {
    nglib::Ng_DeleteMesh(_ngmesh);
  }

  void clear() {
    nglib::Ng_DeleteMesh(_ngmesh);
    _ngmesh = nglib::Ng_NewMesh();
  }

  operator nglib::Ng_Mesh* () {
    return _ngmesh;
  }

private:
  nglib::Ng_Mesh * _ngmesh;
};

}

namespace libMesh
{

//----------------------------------------------------------------------
// NetGenMeshInterface class members
NetGenMeshInterface::NetGenMeshInterface (UnstructuredMesh & mesh) :
  MeshTetInterface(mesh),
  _serializer(mesh)
{
}



void NetGenMeshInterface::triangulate ()
{
  using namespace nglib;

  // We're hoping to do volume_to_surface_mesh in parallel at least,
  // but then we'll need to serialize any hole meshes to rank 0 so it
  // can use them in serial.

  const BoundingBox mesh_bb =
    MeshTetInterface::volume_to_surface_mesh(this->_mesh);

  std::vector<MeshSerializer> hole_serializers;
  if (_holes)
    for (std::unique_ptr<UnstructuredMesh> & hole : *_holes)
      {
        const BoundingBox hole_bb =
          MeshTetInterface::volume_to_surface_mesh(*hole);

        libmesh_error_msg_if
          (!mesh_bb.contains(hole_bb),
           "Found hole with bounding box " << hole_bb <<
           "\nextending outside of mesh bounding box " << mesh_bb);

        hole_serializers.emplace_back
          (*hole, /* need_serial */ true,
           /* serial_only_needed_on_proc_0 */ true);
      }

  // If we're not rank 0, we're just going to wait for rank 0 to call
  // Netgen, then receive its data afterward, we're not going to hope
  // that Netgen does the exact same thing on every processor.
  if (this->_mesh.processor_id() != 0)
    {
      // We don't need our holes anymore.  Delete their serializers
      // first to avoid dereferencing dangling pointers.
      hole_serializers.clear();
      if (_holes)
        _holes->clear();

      // Receive the mesh data rank 0 will send later, then fix it up
      // together
      MeshCommunication().broadcast(this->_mesh);
      this->_mesh.prepare_for_use();
      return;
    }

  this->check_hull_integrity();

  Ng_Meshing_Parameters params;

  // Override any default parameters we might need to, to avoid
  // inserting nodes we don't want.
  params.uselocalh = false;
  params.minh = 0;
  params.elementsperedge = 1;
  params.elementspercurve = 1;
  params.closeedgeenable = false;
  params.closeedgefact = 0;
  params.minedgelenenable = false;
  params.minedgelen = 0;

  // Try to get a no-extra-nodes mesh if we're asked to, or try to
  // translate our desired volume into NetGen terms otherwise.
  //
  // Spoiler alert: all we can do is try; NetGen uses a marching front
  // algorithm that can insert extra nodes despite all my best
  // efforts.
  if (_desired_volume == 0) // shorthand for "no refinement"
    {
      params.maxh = std::numeric_limits<double>::max();
      params.fineness = 0; // "coarse" in the docs
      params.grading = 1;  // "aggressive local grading" to avoid smoothing??

      // Turning off optimization steps avoids another opportunity for
      // Netgen to try to add more nodes.
      params.optsteps_3d = 0;
    }
  else
    params.maxh = double(std::pow(_desired_volume, 1./3.));

  // Keep track of how NetGen copies of nodes map back to our original
  // nodes, so we can connect new elements to nodes correctly.
  std::unordered_map<int, dof_id_type> ng_to_libmesh_id;

  auto handle_ng_result = [](Ng_Result result) {
    static const std::vector<std::string> result_types =
      {"Netgen error", "Netgen success", "Netgen surface input error",
       "Netgen volume failure", "Netgen STL input error",
       "Netgen surface failure", "Netgen file not found"};

    if (result+1 >= 0 &&
        std::size_t(result+1) < result_types.size())
      libmesh_error_msg_if
        (result, "Ng_GenerateVolumeMesh failed: " <<
         result_types[result+1]);
    else
      libmesh_error_msg
        ("Ng_GenerateVolumeMesh failed with an unknown error code");
  };

  WrappedNgMesh ngmesh;

  // Create surface mesh in the WrappedNgMesh
  {
    // NetGen appears to use ONE-BASED numbering for its nodes, and
    // since it doesn't return an id when adding nodes we'll have to
    // track the numbering ourselves.
    int ng_id = 1;

    auto create_surface_component =
      [this, &ng_id, &ng_to_libmesh_id, &ngmesh]
      (UnstructuredMesh & srcmesh, bool hole_mesh)
    {
      // Keep track of what nodes we've already added to the Netgen
      // mesh vs what nodes we need to add.  We'll keep track by id,
      // not by point location.  I don't know if Netgen can handle
      // multiple nodes with the same point location, but if they can
      // it's not going to be *us* who breaks that feature.
      std::unordered_map<dof_id_type, int> libmesh_to_ng_id;

      // Keep track of what nodes we've already added to the main
      // mesh from a hole mesh.
      std::unordered_map<dof_id_type, dof_id_type> hole_to_main_mesh_id;

      // Use a separate array for passing points to NetGen, just in case
      // we're not using double-precision ourselves.
      std::array<double, 3> point_val;

      // And an array for element vertices
      std::array<int, 3> elem_nodes;

      for (const auto * elem : srcmesh.element_ptr_range())
        {
          // If someone has triangles we can't triangulate, we have a
          // problem
          if (elem->type() == TRI6 ||
              elem->type() == TRI7)
            libmesh_not_implemented_msg
              ("Netgen tetrahedralization currently only supports TRI3 boundaries");

          // If someone has non-triangles, let's just ignore them.
          if (elem->type() != TRI3)
            continue;

          for (int ni : make_range(3))
            {
              // Just using the "invert_trigs" option in NetGen params
              // doesn't work for me, so we'll have to have properly
              // oriented the tris earlier.
              auto & elem_node = hole_mesh ? elem_nodes[2-ni] : elem_nodes[ni];

              const Node & n = elem->node_ref(ni);
              auto n_id = n.id();
              if (hole_mesh)
                {
                  if (auto it = hole_to_main_mesh_id.find(n_id);
                      it != hole_to_main_mesh_id.end())
                    {
                      n_id = it->second;
                    }
                  else
                    {
                      Node * n_new = this->_mesh.add_point(n);
                      const dof_id_type n_new_id = n_new->id();
                      hole_to_main_mesh_id.emplace(n_id, n_new_id);
                      n_id = n_new_id;
                    }
                }

              if (auto it = libmesh_to_ng_id.find(n_id);
                  it != libmesh_to_ng_id.end())
                {
                  const int existing_ng_id = it->second;
                  elem_node = existing_ng_id;
                }
              else
                {
                  for (auto i : make_range(3))
                    point_val[i] = double(n(i));

                  Ng_AddPoint(ngmesh, point_val.data());

                  ng_to_libmesh_id[ng_id] = n_id;
                  libmesh_to_ng_id[n_id] = ng_id;
                  elem_node = ng_id;
                  ++ng_id;
                }
            }

          Ng_AddSurfaceElement(ngmesh, NG_TRIG, elem_nodes.data());
        }
    };

    create_surface_component(this->_mesh, false);

    if (_holes)
      for (const std::unique_ptr<UnstructuredMesh> & h : *_holes)
        create_surface_component(*h, true);
  }

  auto result = Ng_GenerateVolumeMesh(ngmesh, &params);
  handle_ng_result(result);

  const int n_elem = Ng_GetNE(ngmesh);

  libmesh_error_msg_if (n_elem <= 0,
                        "NetGen failed to generate any tetrahedra");

  const dof_id_type n_points = Ng_GetNP(ngmesh);
  const dof_id_type old_nodes = this->_mesh.n_nodes();

  // Netgen may have generated new interior nodes
  if (n_points != old_nodes)
    {
      std::array<double, 3> point_val;

      // We should only be getting new nodes if we asked for them
      if (!_desired_volume)
        {
          std::cout <<
            "NetGen output " << n_points <<
            " points when we gave it " <<
            old_nodes << " and disabled refinement\n" <<
            "If new interior points are acceptable in your mesh, please set\n" <<
            "a non-zero desired_volume to indicate that.  If new interior\n" <<
            "points are not acceptable in your mesh, you may need a different\n" <<
            "(non-advancing-front?) mesh generator." << std::endl;
          libmesh_error();
        }
      else
        for (auto i : make_range(old_nodes, n_points))
          {
            // i+1 since ng uses ONE-BASED numbering
            Ng_GetPoint (ngmesh, i+1, point_val.data());
            const Point p(point_val[0], point_val[1], point_val[2]);
            Node * n_new = this->_mesh.add_point(p);
            const dof_id_type n_new_id = n_new->id();
            ng_to_libmesh_id[i+1] = n_new_id;
          }
    }

  for (auto * elem : this->_mesh.element_ptr_range())
    this->_mesh.delete_elem(elem);

  BoundaryInfo * bi = & this->_mesh.get_boundary_info();

  for (auto i : make_range(n_elem))
  {
    // Enough data to return even a Tet10 without a segfault if nglib
    // went nuts
    int ngnodes[11];

    // i+1 since we must be 1-based with these ids too...
    Ng_Volume_Element_Type ngtype =
      Ng_GetVolumeElement(ngmesh, i+1, ngnodes);

    // But really nglib shouldn't go nuts
    libmesh_assert(ngtype == NG_TET);
    libmesh_ignore(ngtype);

    auto elem = this->_mesh.add_elem(Elem::build_with_id(TET4, i));
    for (auto n : make_range(4))
      {
        const dof_id_type node_id =
          libmesh_map_find(ng_to_libmesh_id, ngnodes[n]);
        elem->set_node(n) = this->_mesh.node_ptr(node_id);
      }

    // NetGen and we disagree about node numbering orientation
    elem->orient(bi);
  }

  // We don't need our holes anymore.  Delete their serializers
  // first to avoid dereferencing dangling pointers.
  hole_serializers.clear();
  if (_holes)
    _holes->clear();

  // Send our data to other ranks, then fix it up together
  MeshCommunication().broadcast(this->_mesh);
  this->_mesh.prepare_for_use();
}



} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_NETGEN
