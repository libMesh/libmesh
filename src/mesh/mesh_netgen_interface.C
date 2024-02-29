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

  this->check_hull_integrity();

  WrappedNgMesh ngmesh;

  Ng_Meshing_Parameters params;

  // Override any default parameters we might need to
  params.uselocalh = false;
  params.maxh = std::pow(_desired_volume, 1./3.);
  params.minh = 0;
  params.elementsperedge = 1;
  params.elementspercurve = 1;
  params.closeedgeenable = false;
  params.closeedgefact = 0;
  params.minedgelenenable = false;
  params.minedgelen = 0;

  // Keep track of how NetGen copies of nodes map back to our original
  // nodes, so we can connect new elements to nodes correctly.
  std::unordered_map<int, dof_id_type> ng_to_libmesh_id;

  auto create_surface_mesh = [this, &ng_to_libmesh_id](WrappedNgMesh & wngmesh,
                                                       bool invert_tris) {
    // Keep track of what nodes we've already added to the Netgen mesh
    // vs what nodes we need to add.  We'll keep track by id, not by
    // point location.  I don't know if Netgen can handle multiple
    // nodes with the same point location, but if they can it's not
    // going to be *us* who breaks that feature.
    std::unordered_map<dof_id_type, int> libmesh_to_ng_id;

    // NetGen appears to use ONE-BASED numbering for its nodes, and
    // since it doesn't return an id when adding nodes we'll have to
    // track the numbering ourselves.
    int ng_id = 1;

    // Use a separate array for passing points to NetGen, just in case
    // we're not using double-precision ourselves.
    std::array<double, 3> point_val;

    // And an array for element vertices
    std::array<int, 3> elem_nodes;

    for (const auto * elem : this->_mesh.element_ptr_range())
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
            // Just using the "invert_trigs" option in params doesn't
            // work for me, so we'll try inverting them manually if
            // need be.
            auto & elem_node = invert_tris ? elem_nodes[2-ni] : elem_nodes[ni];

            const Node & n = elem->node_ref(ni);
            auto it = libmesh_to_ng_id.find(n.id());
            if (it == libmesh_to_ng_id.end())
              {
                for (auto i : make_range(3))
                  point_val[i] = n(i);

                Ng_AddPoint(wngmesh, point_val.data());

                ng_to_libmesh_id[ng_id] = n.id();
                libmesh_to_ng_id[n.id()] = ng_id;
                elem_node = ng_id;
                ++ng_id;
              }
            else
              {
                const int existing_ng_id = it->second;
                elem_node = existing_ng_id;
              }
          }

        Ng_AddSurfaceElement(wngmesh, NG_TRIG, elem_nodes.data());
      }
  };

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

  // We want to support any orientation of input triangles, but NetGen
  // wants us to tell them the triangle orientation a priori or it
  // won't generate any volume elements.  Let's just cheat and try
  // twice.
  create_surface_mesh(ngmesh, /* invert_tris = */ false);
  auto result = Ng_GenerateVolumeMesh(ngmesh, &params);
  handle_ng_result(result);

  int n_elem = Ng_GetNE(ngmesh);

  if (!n_elem)
    {
      ngmesh.clear();
      create_surface_mesh(ngmesh, /* invert_tris = */ true);
      result = Ng_GenerateVolumeMesh(ngmesh, &params);
      handle_ng_result(result);
    }

  n_elem = Ng_GetNE(ngmesh);

  libmesh_error_msg_if (n_elem <= 0,
                        "NetGen failed to generate any tetrahedra");

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

  this->_mesh.prepare_for_use();
}



} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_NETGEN
