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
#include "libmesh/cell_tet10.h"
#include "libmesh/elem.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/threads.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/utility.h" // libmesh_map_find

namespace nglib {
#include "netgen/nglib/nglib.h"
}

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

  LOG_SCOPE("triangulate()", "NetGenMeshInterface");

  if (_elem_type != TET4 &&
      _elem_type != TET10 &&
      _elem_type != TET14)
    libmesh_not_implemented();

  // We're hoping to do volume_to_surface_mesh in parallel at least,
  // but then we'll need to serialize any hole meshes to rank 0 so it
  // can use them in serial.

  // If the user wants higher-order output, record midpoints from any
  // quadratic boundary elements before stripping them to TRI3 for
  // NetGen.  We key by sorted position pairs (not node IDs) so that
  // outer mesh and hole mesh node namespaces cannot conflict.
  std::map<std::pair<Point,Point>, Point> edge_midpoints;

  // TRI7 boundary faces additionally carry a face-centroid node (local index
  // 6) that TET14 output reproduces (its faces are TRI7).  Record those too,
  // keyed by the sorted triple of the face's three corner positions so the
  // key is independent of node id and rotation.
  std::map<std::array<Point,3>, Point> face_centroids;

  auto record_and_strip =
    [&edge_midpoints, &face_centroids](UnstructuredMesh & m)
    {
      bool has_quadratic = false;
      for (const auto & elem : m.element_ptr_range())
        {
          if (elem->type() != TRI6 && elem->type() != TRI7) continue;
          has_quadratic = true;
          // TRI6/TRI7: edge e runs node[e]→node[(e+1)%3], midpoint=node[e+3]
          for (auto e : make_range(3u))
            {
              const Point & pa = elem->point(e);
              const Point & pb = elem->point((e+1)%3);
              auto key = pa < pb
                ? std::make_pair(pa,pb) : std::make_pair(pb,pa);
              edge_midpoints[key] = elem->point(e+3);
            }
          // TRI7: node 6 is the face-centroid node.
          if (elem->type() == TRI7)
            {
              std::array<Point,3> corners =
                {elem->point(0), elem->point(1), elem->point(2)};
              std::sort(corners.begin(), corners.end());
              face_centroids[corners] = elem->point(6);
            }
        }
      if (has_quadratic)
        m.all_first_order();
    };

  const BoundingBox mesh_bb =
    MeshTetInterface::volume_to_surface_mesh(this->_mesh);

  if (_elem_type != TET4)
    record_and_strip(this->_mesh);

  std::vector<MeshSerializer> hole_serializers;
  if (_holes)
    for (std::unique_ptr<UnstructuredMesh> & hole : *_holes)
      {
        const BoundingBox hole_bb =
          MeshTetInterface::volume_to_surface_mesh(*hole);

        if (_elem_type != TET4)
          record_and_strip(*hole);

        libmesh_error_msg_if
          (!mesh_bb.contains(hole_bb),
           "Found hole with bounding box " << hole_bb <<
           "\nextending outside of mesh bounding box " << mesh_bb);

        hole_serializers.emplace_back
          (*hole, /* need_serial */ true,
           /* serial_only_needed_on_proc_0 */ true);
      }

  // Increasing the element order (all_second_order()/all_complete_order()
  // via increase_tet_order()) performs collective MPI communication, so it
  // must run on every rank in lockstep -- it cannot run on rank 0 alone
  // while the other ranks wait in broadcast().  We therefore broadcast
  // NetGen's TET4 result first, then increase the order and restore the
  // curved-boundary midpoints identically on all ranks.  This mirrors what
  // the 2D interfaces do (see TriangulatorInterface::increase_triangle_order,
  // which is likewise called on all ranks).  edge_midpoints was built while
  // the mesh was serialized on every rank, so it is identical everywhere and
  // this post-broadcast fixup is deterministic.
  auto increase_order_and_restore_midpoints =
    [this, &edge_midpoints, &face_centroids]()
    {
      if (_elem_type == TET4)
        return;

      // find_neighbors() is needed before all_second_order() can place
      // shared edge midpoints correctly.
      this->_mesh.find_neighbors();

      // Refresh the cached element dimensions.  We just replaced the 2D
      // TRI3 surface elements with 3D TET4 volume elements, but the mesh's
      // cached dimension is still that of the original surface (2).
      // all_second_order()/all_complete_order() derive the per-element
      // unique_id reservation width from mesh_dimension(): a stale value of
      // 2 reserves only 9-4=5 slots per element, too few for the 6 new edge
      // nodes of a TET10, so adjacent elements' unique_id ranges overlap and
      // collide.  cache_elem_data() recomputes the dimension to 3 (reserving
      // 27-8=19 slots) before the order increase runs.
      this->_mesh.cache_elem_data();

      this->increase_tet_order();

      // Move auto-placed geometric midpoints to the recorded positions,
      // preserving any curvature from the original quadratic boundary.
      if (!edge_midpoints.empty() || !face_centroids.empty())
        for (Elem * elem : _mesh.element_ptr_range())
          for (auto s : elem->side_index_range())
            {
              if (elem->neighbor_ptr(s)) continue;

              // build_side_ptr() returns a TRI6 (TET10) or TRI7 (TET14)
              // whose node pointers reference the actual mesh nodes; point()
              // assignments update mesh node coordinates in place.
              auto side = elem->build_side_ptr(s);
              for (auto e : make_range(3u))
                {
                  const Point & pa = side->point(e);
                  const Point & pb = side->point((e+1)%3);
                  auto key = pa < pb
                    ? std::make_pair(pa,pb) : std::make_pair(pb,pa);
                  if (auto it = edge_midpoints.find(key);
                      it != edge_midpoints.end())
                    side->point(e+3) = it->second;
                }

              // TRI7 faces (TET14 output) also carry a face-centroid node
              // at local index 6; restore its recorded curved position.
              if (side->type() == TRI7)
                {
                  std::array<Point,3> corners =
                    {side->point(0), side->point(1), side->point(2)};
                  std::sort(corners.begin(), corners.end());
                  if (auto it = face_centroids.find(corners);
                      it != face_centroids.end())
                    side->point(6) = it->second;
                }
            }
    };

  // This should probably only be done on rank 0, but the API is
  // designed with the hope that we'll parallelize it eventually
  auto integrity = this->improve_hull_integrity();
  this->process_hull_integrity_result(integrity);

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

      // Receive the TET4 mesh data rank 0 will send later.
      MeshCommunication().broadcast(this->_mesh);

      // If we got an empty mesh here then our tetrahedralization
      // failed.
      libmesh_error_msg_if (!this->_mesh.n_elem(),
                            "NetGen failed to generate any tetrahedra");

      // Increase the element order collectively, in lockstep with rank 0.
      increase_order_and_restore_midpoints();

      this->_mesh.prepare_for_use();
      return;
    }

  Ng_Meshing_Parameters params;

  Ng_SetNumThreads(cast_int<int>(libMesh::n_threads()));

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

  // Keep track of what boundary ids we want to assign to each new
  // triangle.  We'll give the outer boundary BC 0, and give holes ids
  // starting from 1.
  // We key on sorted tuples of node ids to identify a side.
  std::unordered_map<std::array<dof_id_type,3>,
                     boundary_id_type, libMesh::hash> side_boundary_id;

  auto insert_id = []
    (std::array<dof_id_type,3> & array,
     dof_id_type n_id)
  {
    libmesh_assert_less(n_id, DofObject::invalid_id);
    unsigned int i=0;
    while (array[i] < n_id)
      ++i;
    while (i < 3)
      std::swap(array[i++], n_id);
  };

  WrappedNgMesh ngmesh;

  // Create surface mesh in the WrappedNgMesh
  {
    // NetGen appears to use ONE-BASED numbering for its nodes, and
    // since it doesn't return an id when adding nodes we'll have to
    // track the numbering ourselves.
    int ng_id = 1;

    auto create_surface_component =
      [this, &ng_id, &ng_to_libmesh_id, &ngmesh, &side_boundary_id, &insert_id]
      (UnstructuredMesh & srcmesh, bool hole_mesh,
       boundary_id_type bcid)
    {
      LOG_SCOPE("create_surface_component()", "NetGenMeshInterface");

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

          std::array<dof_id_type,3> sorted_ids =
            {DofObject::invalid_id, DofObject::invalid_id,
             DofObject::invalid_id};

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

              insert_id(sorted_ids, n_id);
            }

          side_boundary_id[sorted_ids] = bcid;

          Ng_AddSurfaceElement(ngmesh, NG_TRIG, elem_nodes.data());
        }
    };

    // Number the outer boundary 0, and the holes starting from 1
    boundary_id_type bcid = 0;

    create_surface_component(this->_mesh, false, bcid);

    if (_holes)
      for (const std::unique_ptr<UnstructuredMesh> & h : *_holes)
        create_surface_component(*h, true, ++bcid);
  }

  {
    LOG_SCOPE("Ng_GenerateVolumeMesh()", "NetGenMeshInterface");

    auto result = Ng_GenerateVolumeMesh(ngmesh, &params);
    handle_ng_result(result);
  }

  const int n_elem = Ng_GetNE(ngmesh);

  // If Netgen fails us, we're likely to get n_elem <= 0.  This is a
  // common enough failure from bad setups that I want to make sure
  // it's thrown in parallel so as to not desynchronize any unit tests
  // that trigger it.  So we'll broadcast the empty mesh to indicate
  // the problem and enable throwing exceptions in parallel.
  if (n_elem <= 0)
    {
      this->_mesh.clear();
      MeshCommunication().broadcast(this->_mesh);
      libmesh_error_msg ("NetGen failed to generate any tetrahedra");
    }

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
        elem->set_node(n, this->_mesh.node_ptr(node_id));
      }

    // NetGen and we disagree about node numbering orientation
    elem->orient(bi);

    for (auto s : make_range(4))
      {
        std::array<dof_id_type,3> sorted_ids =
          {DofObject::invalid_id, DofObject::invalid_id,
           DofObject::invalid_id};

        std::vector<unsigned int> nos = elem->nodes_on_side(s);
        for (auto n : nos)
          insert_id(sorted_ids, elem->node_id(n));

        if (auto it = side_boundary_id.find(sorted_ids);
            it != side_boundary_id.end())
          bi->add_side(elem, s, it->second);
      }
  }

  // We don't need our holes anymore.  Delete their serializers
  // first to avoid dereferencing dangling pointers.
  hole_serializers.clear();
  if (_holes)
    _holes->clear();

  // Send NetGen's TET4 result to the other ranks, then increase the element
  // order collectively on all ranks together.  increase_tet_order() performs
  // collective communication, so it must not run on rank 0 alone -- doing so
  // would deadlock against the other ranks waiting in broadcast() above.
  MeshCommunication().broadcast(this->_mesh);
  increase_order_and_restore_midpoints();
  this->_mesh.prepare_for_use();
}



} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_NETGEN
