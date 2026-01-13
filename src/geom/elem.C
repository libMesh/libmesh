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
#include "libmesh/elem.h"

#include "libmesh/boundary_info.h"
#include "libmesh/fe_type.h"
#include "libmesh/fe_interface.h"
#include "libmesh/node_elem.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_edge4.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/face_c0polygon.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri3_subdivision.h"
#include "libmesh/face_tri3_shell.h"
#include "libmesh/face_tri6.h"
#include "libmesh/face_tri7.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_quad4_shell.h"
#include "libmesh/face_quad8.h"
#include "libmesh/face_quad8_shell.h"
#include "libmesh/face_quad9.h"
#include "libmesh/face_quad9_shell.h"
#include "libmesh/face_inf_quad4.h"
#include "libmesh/face_inf_quad6.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_tet10.h"
#include "libmesh/cell_tet14.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex20.h"
#include "libmesh/cell_hex27.h"
#include "libmesh/cell_inf_hex8.h"
#include "libmesh/cell_inf_hex16.h"
#include "libmesh/cell_inf_hex18.h"
#include "libmesh/cell_prism6.h"
#include "libmesh/cell_prism15.h"
#include "libmesh/cell_prism18.h"
#include "libmesh/cell_prism20.h"
#include "libmesh/cell_prism21.h"
#include "libmesh/cell_inf_prism6.h"
#include "libmesh/cell_inf_prism12.h"
#include "libmesh/cell_pyramid5.h"
#include "libmesh/cell_pyramid13.h"
#include "libmesh/cell_pyramid14.h"
#include "libmesh/cell_pyramid18.h"
#include "libmesh/fe_base.h"
#include "libmesh/mesh_base.h"
#include "libmesh/quadrature_nodal.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/remote_elem.h"
#include "libmesh/reference_elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/threads.h"
#include "libmesh/enum_elem_quality.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"
#include "libmesh/elem_internal.h"

#ifdef LIBMESH_ENABLE_PERIODIC
#include "libmesh/mesh.h"
#include "libmesh/periodic_boundaries.h"
#endif


// C++ includes
#include <algorithm> // for std::sort
#include <array>
#include <iterator>  // for std::ostream_iterator
#include <sstream>
#include <limits>    // for std::numeric_limits<>
#include <cmath>     // for std::sqrt()
#include <memory>
#include <regex>     // for exceptions in volume()


namespace libMesh
{

Threads::spin_mutex parent_indices_mutex;
Threads::spin_mutex parent_bracketing_nodes_mutex;

// Initialize static member variables
const unsigned int Elem::type_to_dim_map [] =
  {
    1,  // EDGE2
    1,  // EDGE3
    1,  // EDGE4

    2,  // TRI3
    2,  // TRI6

    2,  // QUAD4
    2,  // QUAD8
    2,  // QUAD9

    3,  // TET4
    3,  // TET10

    3,  // HEX8
    3,  // HEX20
    3,  // HEX27

    3,  // PRISM6
    3,  // PRISM15
    3,  // PRISM18

    3,  // PYRAMID5
    3,  // PYRAMID13
    3,  // PYRAMID14

    1,  // INFEDGE2

    2,  // INFQUAD4
    2,  // INFQUAD6

    3,  // INFHEX8
    3,  // INFHEX16
    3,  // INFHEX18

    3,  // INFPRISM6
    3,  // INFPRISM12

    0,  // NODEELEM

    0,  // REMOTEELEM

    2,  // TRI3SUBDIVISION
    2,  // TRISHELL3
    2,  // QUADSHELL4
    2,  // QUADSHELL8

    2,  // TRI7
    3,  // TET14
    3,  // PRISM20
    3,  // PRISM21
    3,  // PYRAMID18

    2,  // QUADSHELL9

    2,  // C0POLYGON
    3,  // C0POLYHEDRON
  };

const unsigned int Elem::max_n_nodes;

const unsigned int Elem::type_to_n_nodes_map [] =
  {
    2,  // EDGE2
    3,  // EDGE3
    4,  // EDGE4

    3,  // TRI3
    6,  // TRI6

    4,  // QUAD4
    8,  // QUAD8
    9,  // QUAD9

    4,  // TET4
    10, // TET10

    8,  // HEX8
    20, // HEX20
    27, // HEX27

    6,  // PRISM6
    15, // PRISM15
    18, // PRISM18

    5,  // PYRAMID5
    13, // PYRAMID13
    14, // PYRAMID14

    2,  // INFEDGE2

    4,  // INFQUAD4
    6,  // INFQUAD6

    8,  // INFHEX8
    16, // INFHEX16
    18, // INFHEX18

    6,  // INFPRISM6
    12, // INFPRISM12

    1,  // NODEELEM

    0,  // REMOTEELEM

    3,  // TRI3SUBDIVISION
    3,  // TRISHELL3
    4,  // QUADSHELL4
    8,  // QUADSHELL8

    7,  // TRI7
    14, // TET14
    20, // PRISM20
    21, // PRISM21
    18, // PYRAMID18

    9,  // QUADSHELL9

    invalid_uint,  // C0POLYGON
    invalid_uint,  // C0POLYHEDRON
  };

const unsigned int Elem::type_to_n_sides_map [] =
  {
    2,  // EDGE2
    2,  // EDGE3
    2,  // EDGE4

    3,  // TRI3
    3,  // TRI6

    4,  // QUAD4
    4,  // QUAD8
    4,  // QUAD9

    4,  // TET4
    4,  // TET10

    6,  // HEX8
    6,  // HEX20
    6,  // HEX27

    5,  // PRISM6
    5,  // PRISM15
    5,  // PRISM18

    5,  // PYRAMID5
    5,  // PYRAMID13
    5,  // PYRAMID14

    2,  // INFEDGE2

    3,  // INFQUAD4
    3,  // INFQUAD6

    5,  // INFHEX8
    5,  // INFHEX16
    5,  // INFHEX18

    4,  // INFPRISM6
    4,  // INFPRISM12

    0,  // NODEELEM

    0,  // REMOTEELEM

    3,  // TRI3SUBDIVISION
    3,  // TRISHELL3
    4,  // QUADSHELL4
    4,  // QUADSHELL8

    3,  // TRI7
    4,  // TET14
    5,  // PRISM20
    5,  // PRISM21
    5,  // PYRAMID18

    4,  // QUADSHELL9

    invalid_uint,  // C0POLYGON
    invalid_uint,  // C0POLYHEDRON
  };

const unsigned int Elem::type_to_n_edges_map [] =
  {
    0,  // EDGE2
    0,  // EDGE3
    0,  // EDGE4

    3,  // TRI3
    3,  // TRI6

    4,  // QUAD4
    4,  // QUAD8
    4,  // QUAD9

    6,  // TET4
    6,  // TET10

    12, // HEX8
    12, // HEX20
    12, // HEX27

    9,  // PRISM6
    9,  // PRISM15
    9,  // PRISM18

    8,  // PYRAMID5
    8,  // PYRAMID13
    8,  // PYRAMID14

    0,  // INFEDGE2

    3,  // INFQUAD4
    3,  // INFQUAD6

    8,  // INFHEX8
    8,  // INFHEX16
    8,  // INFHEX18

    6,  // INFPRISM6
    6,  // INFPRISM12

    0,  // NODEELEM

    0,  // REMOTEELEM

    3,  // TRI3SUBDIVISION
    3,  // TRISHELL3
    4,  // QUADSHELL4
    4,  // QUADSHELL8

    3,  // TRI7
    6,  // TET14
    9,  // PRISM20
    9,  // PRISM21
    8,  // PYRAMID18

    4,  // QUADSHELL9

    invalid_uint,  // C0POLYGON
    invalid_uint,  // C0POLYHEDRON
  };

const Order Elem::type_to_default_order_map [] =
  {
    FIRST,    // EDGE2
    SECOND,   // EDGE3
    THIRD,    // EDGE4

    FIRST,    // TRI3
    SECOND,   // TRI6

    FIRST,    // QUAD4
    SECOND,   // QUAD8
    SECOND,   // QUAD9

    FIRST,    // TET4
    SECOND,   // TET10

    FIRST,    // HEX8
    SECOND,   // HEX20
    SECOND,   // HEX27

    FIRST,    // PRISM6
    SECOND,   // PRISM15
    SECOND,   // PRISM18

    FIRST,    // PYRAMID5
    SECOND,   // PYRAMID13
    SECOND,   // PYRAMID14

    FIRST,    // INFEDGE2

    FIRST,    // INFQUAD4
    SECOND,   // INFQUAD6

    FIRST,    // INFHEX8
    SECOND,   // INFHEX16
    SECOND,   // INFHEX18

    FIRST,    // INFPRISM6
    SECOND,   // INFPRISM12

    CONSTANT, // NODEELEM

    INVALID_ORDER, // REMOTEELEM

    FIRST,    // TRI3SUBDIVISION
    FIRST,    // TRISHELL3
    FIRST,    // QUADSHELL4
    SECOND,   // QUADSHELL8

    THIRD,    // TRI7
    THIRD,    // TET14
    THIRD,    // PRISM20
    THIRD,    // PRISM21
    THIRD,    // PYRAMID18

    SECOND,   // QUADSHELL9

    FIRST,    // C0POLYGON
    FIRST,    // C0POLYHEDRON
  };

// ------------------------------------------------------------
// Elem class member functions
std::unique_ptr<Elem> Elem::disconnected_clone() const
{
  std::unique_ptr<Elem> returnval;

  switch (this->type())
    {
    case C0POLYGON:
      returnval = std::make_unique<C0Polygon>(this->n_sides());
      break;

    default:
      returnval = Elem::build(this->type());
    }

  returnval->set_id() = this->id();
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (this->valid_unique_id())
    returnval->set_unique_id(this->unique_id());
#endif

  const auto n_elem_ints = this->n_extra_integers();
  returnval->add_extra_integers(n_elem_ints);
  for (unsigned int i = 0; i != n_elem_ints; ++i)
    returnval->set_extra_integer(i, this->get_extra_integer(i));

  returnval->inherit_data_from(*this);

  return returnval;
}



std::unique_ptr<Elem> Elem::build(const ElemType type,
                                  Elem * p)
{
  switch (type)
    {
      // 0D elements
    case NODEELEM:
      return std::make_unique<NodeElem>(p);

      // 1D elements
    case EDGE2:
      return std::make_unique<Edge2>(p);
    case EDGE3:
      return std::make_unique<Edge3>(p);
    case EDGE4:
      return std::make_unique<Edge4>(p);

      // 2D elements
    case TRI3:
      return std::make_unique<Tri3>(p);
    case TRISHELL3:
      return std::make_unique<TriShell3>(p);
    case TRI3SUBDIVISION:
      return std::make_unique<Tri3Subdivision>(p);
    case TRI6:
      return std::make_unique<Tri6>(p);
    case TRI7:
      return std::make_unique<Tri7>(p);
    case QUAD4:
      return std::make_unique<Quad4>(p);
    case QUADSHELL4:
      return std::make_unique<QuadShell4>(p);
    case QUAD8:
      return std::make_unique<Quad8>(p);
    case QUADSHELL8:
      return std::make_unique<QuadShell8>(p);
    case QUAD9:
      return std::make_unique<Quad9>(p);
    case QUADSHELL9:
      return std::make_unique<QuadShell9>(p);

    // Well, a hexagon is *a* polygon...
    case C0POLYGON:
      return std::make_unique<C0Polygon>(6, p);

    // Building a polyhedron can't currently be done without creating
    // its nodes first
    case C0POLYHEDRON:
      libmesh_not_implemented_msg
        ("Polyhedra cannot be built via Elem::build()");

      // 3D elements
    case TET4:
      return std::make_unique<Tet4>(p);
    case TET10:
      return std::make_unique<Tet10>(p);
    case TET14:
      return std::make_unique<Tet14>(p);
    case HEX8:
      return std::make_unique<Hex8>(p);
    case HEX20:
      return std::make_unique<Hex20>(p);
    case HEX27:
      return std::make_unique<Hex27>(p);
    case PRISM6:
      return std::make_unique<Prism6>(p);
    case PRISM15:
      return std::make_unique<Prism15>(p);
    case PRISM18:
      return std::make_unique<Prism18>(p);
    case PRISM20:
      return std::make_unique<Prism20>(p);
    case PRISM21:
      return std::make_unique<Prism21>(p);
    case PYRAMID5:
      return std::make_unique<Pyramid5>(p);
    case PYRAMID13:
      return std::make_unique<Pyramid13>(p);
    case PYRAMID14:
      return std::make_unique<Pyramid14>(p);
    case PYRAMID18:
      return std::make_unique<Pyramid18>(p);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
      // 1D infinite elements
    case INFEDGE2:
      return std::make_unique<InfEdge2>(p);

      // 2D infinite elements
    case INFQUAD4:
      return std::make_unique<InfQuad4>(p);
    case INFQUAD6:
      return std::make_unique<InfQuad6>(p);

      // 3D infinite elements
    case INFHEX8:
      return std::make_unique<InfHex8>(p);
    case INFHEX16:
      return std::make_unique<InfHex16>(p);
    case INFHEX18:
      return std::make_unique<InfHex18>(p);
    case INFPRISM6:
      return std::make_unique<InfPrism6>(p);
    case INFPRISM12:
      return std::make_unique<InfPrism12>(p);
#endif

    default:
      libmesh_error_msg("ERROR: Undefined element type == " << Utility::enum_to_string(type));
    }
}



std::unique_ptr<Elem> Elem::build_with_id (const ElemType type,
                                           dof_id_type id)
{
  // Call the other build() method with nullptr parent, then set the
  // required id.
  auto temp = Elem::build(type, nullptr);
  temp->set_id(id);
  return temp;
}



const Elem * Elem::reference_elem () const
{
  return &(ReferenceElem::get(this->type()));
}



Point Elem::true_centroid() const
{
  // The base class implementation builds a finite element of the correct
  // order and computes the centroid, c=(cx, cy, cz), where:
  //
  // [cx]            [\int x dV]
  // [cy] := (1/V) * [\int y dV]
  // [cz]            [\int z dV]
  //
  // using quadrature. Note that we can expand "x" in the FE space as:
  //
  // x = \sum_i x_i \phi_i
  //
  // where x_i are the nodal positions of the element and \phi_i are the
  // associated Lagrange shape functions. This allows us to write the
  // integrals above as e.g.:
  //
  // \int x dV = \sum_i x_i \int \phi_i dV
  //
  // Defining:
  //
  // V_i := \int \phi_i dV
  //
  // we then have:
  //
  // [cx]           [\sum_i x_i V_i]
  // [cy] = (1/V) * [\sum_i y_i V_i]
  // [cz]           [\sum_i z_i V_i]
  //
  // where:
  // V = \sum_i V_i
  //
  // Derived element types can overload this method to compute
  // the centroid more efficiently when possible.

  // If this Elem has an elevated p_level, then we need to generate a
  // barebones copy of it with zero p_level and call true_centroid()
  // on that instead.  This workaround allows us to avoid issues with
  // calling FE::reinit() with a default_order() FEType, and then
  // having that order incorrectly boosted by p_level.
  if (this->p_level())
    {
      auto elem_copy = this->disconnected_clone();
#ifdef LIBMESH_ENABLE_AMR
      elem_copy->set_p_level(0);
#endif

      // Set node pointers
      for (auto n : this->node_index_range())
        elem_copy->set_node(n, _nodes[n]);

      return elem_copy->true_centroid();
    }

  const FEFamily mapping_family = FEMap::map_fe_type(*this);
  const FEType fe_type(this->default_order(), mapping_family);

  // Build FE and attach quadrature rule.  The default quadrature rule
  // integrates the mass matrix exactly, thus it is overkill to
  // integrate the basis functions, but this is convenient.
  std::unique_ptr<FEBase> fe = FEBase::build(this->dim(), fe_type);
  QGauss qrule (this->dim(), fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);

  // Pre-request required data
  const auto & JxW = fe->get_JxW();
  const auto & phi = fe->get_phi();

  // Re-compute element-specific values
  fe->reinit(this);

  // Number of basis functions
  auto N = phi.size();
  libmesh_assert_equal_to(N, this->n_nodes());

  // Compute V_i
  std::vector<Real> V(N);
  for (auto qp : index_range(JxW))
    for (auto i : make_range(N))
      V[i] += JxW[qp] * phi[i][qp];

  // Compute centroid
  Point cp;
  Real vol = 0.;

  for (auto i : make_range(N))
    {
      cp += this->point(i) * V[i];
      vol += V[i];
    }

  return cp / vol;
}

Point Elem::vertex_average() const
{
  Point cp;

  const auto n_vertices = this->n_vertices();

  for (unsigned int n=0; n<n_vertices; n++)
    cp.add (this->point(n));

  return (cp /= static_cast<Real>(n_vertices));
}



Real Elem::hmin() const
{
  Real h_min=std::numeric_limits<Real>::max();

  // Avoid calling a virtual a lot of times
  const auto n_vertices = this->n_vertices();

  for (unsigned int n_outer=0; n_outer<n_vertices; n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<n_vertices; n_inner++)
      {
        const auto diff = (this->point(n_outer) - this->point(n_inner));

        h_min = std::min(h_min, diff.norm_sq());
      }

  return std::sqrt(h_min);
}



Real Elem::hmax() const
{
  Real h_max=0;

  // Avoid calling a virtual a lot of times
  const auto n_vertices = this->n_vertices();

  for (unsigned int n_outer=0; n_outer<n_vertices; n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<n_vertices; n_inner++)
      {
        const auto diff = (this->point(n_outer) - this->point(n_inner));

        h_max = std::max(h_max, diff.norm_sq());
      }

  return std::sqrt(h_max);
}



Real Elem::length(const unsigned int n1,
                  const unsigned int n2) const
{
  libmesh_assert_less ( n1, this->n_vertices() );
  libmesh_assert_less ( n2, this->n_vertices() );

  return (this->point(n1) - this->point(n2)).norm();
}



dof_id_type Elem::key () const
{
  const unsigned short n_n = this->n_nodes();

  std::array<dof_id_type, Elem::max_n_nodes> node_ids;

  for (unsigned short n=0; n != n_n; ++n)
    node_ids[n] = this->node_id(n);

  // Always sort, so that different local node numberings hash to the
  // same value.
  std::sort (node_ids.begin(), node_ids.begin()+n_n);

  return Utility::hashword(node_ids.data(), n_n);
}



bool Elem::operator == (const Elem & rhs) const
{
  // If the elements aren't the same type, they aren't equal
  if (this->type() != rhs.type())
    return false;

  const unsigned short n_n = this->n_nodes();
  libmesh_assert_equal_to(n_n, rhs.n_nodes());

  // Make two sorted arrays of global node ids and compare them for
  // equality.
  std::array<dof_id_type, Elem::max_n_nodes> this_ids, rhs_ids;

  for (unsigned short n = 0; n != n_n; n++)
    {
      this_ids[n] = this->node_id(n);
      rhs_ids[n] = rhs.node_id(n);
    }

  // Sort the vectors to rule out different local node numberings.
  std::sort(this_ids.begin(), this_ids.begin()+n_n);
  std::sort(rhs_ids.begin(), rhs_ids.begin()+n_n);

  // If the node ids match, the elements are equal!
  for (unsigned short n = 0; n != n_n; ++n)
    if (this_ids[n] != rhs_ids[n])
      return false;
  return true;
}



bool Elem::topologically_equal (const Elem & rhs) const
{
  // If the elements aren't the same type, they aren't equal
  if (this->type() != rhs.type())
    return false;

  libmesh_assert_equal_to(this->n_nodes(), rhs.n_nodes());

  for (auto n : make_range(this->n_nodes()))
    if (this->node_id(n) != rhs.node_id(n))
      return false;

  for (auto neigh : make_range(this->n_neighbors()))
    {
      if (!this->neighbor_ptr(neigh))
        {
          if (rhs.neighbor_ptr(neigh))
            return false;
          continue;
        }
      if (!rhs.neighbor_ptr(neigh) ||
          this->neighbor_ptr(neigh)->id() !=
          rhs.neighbor_ptr(neigh)->id())
        return false;
    }

  if (this->parent())
    {
      if (!rhs.parent())
        return false;
      if (this->parent()->id() != rhs.parent()->id())
        return false;
    }
  else if (rhs.parent())
    return false;

  if (this->interior_parent())
    {
      if (!rhs.interior_parent())
        return false;
      if (this->interior_parent()->id() !=
          rhs.interior_parent()->id())
        return false;
    }
  else if (rhs.interior_parent())
    return false;

  return true;
}



bool Elem::is_semilocal(const processor_id_type my_pid) const
{
  std::set<const Elem *> point_neighbors;

  this->find_point_neighbors(point_neighbors);

  for (const auto & elem : point_neighbors)
    if (elem->processor_id() == my_pid)
      return true;

  return false;
}



unsigned int Elem::which_side_am_i (const Elem * e) const
{
  libmesh_assert(e);

  const unsigned int ns = this->n_sides();
  const unsigned int nn = this->n_nodes();

  const unsigned int en = e->n_nodes();

  // e might be on any side until proven otherwise
  std::vector<bool> might_be_side(ns, true);

  for (unsigned int i=0; i != en; ++i)
    {
      Point side_point = e->point(i);
      unsigned int local_node_id = libMesh::invalid_uint;

      // Look for a node of this that's contiguous with node i of
      // e. Note that the exact floating point comparison of Point
      // positions is intentional, see the class documentation for
      // this function.
      for (unsigned int j=0; j != nn; ++j)
        if (this->point(j) == side_point)
          local_node_id = j;

      // If a node of e isn't contiguous with some node of this, then
      // e isn't a side of this.
      if (local_node_id == libMesh::invalid_uint)
        return libMesh::invalid_uint;

      // If a node of e isn't contiguous with some node on side s of
      // this, then e isn't on side s.
      for (unsigned int s=0; s != ns; ++s)
        if (!this->is_node_on_side(local_node_id, s))
          might_be_side[s] = false;
    }

  for (unsigned int s=0; s != ns; ++s)
    if (might_be_side[s])
      {
#ifdef DEBUG
        for (unsigned int s2=s+1; s2 < ns; ++s2)
          libmesh_assert (!might_be_side[s2]);
#endif
        return s;
      }

  // Didn't find any matching side
  return libMesh::invalid_uint;
}



bool Elem::contains_vertex_of(const Elem * e, bool mesh_connection) const
{
  // Our vertices are the first numbered nodes
  const unsigned int nv = e->n_vertices();
  const unsigned int my_nv = this->n_vertices();

  // Check for vertex-to-vertex containment first; contains_point() is
  // expensive
  for (auto n : make_range(nv))
    {
      const Node * vertex = e->node_ptr(n);
      for (auto my_n : make_range(my_nv))
        if (&this->node_ref(my_n) == vertex)
          return true;
    }

  // If e is in our mesh, then we might be done testing
  if (mesh_connection)
    {
      const unsigned int l = this->level();
      const unsigned int el = e->level();

      if (l >= el)
        return false;

      // We could also return false for l==el-1 iff we knew we had no
      // triangular faces, but we don't have an API to check that.
    }

  // Our vertices are the first numbered nodes
  for (auto n : make_range(nv))
    if (this->contains_point(e->point(n)))
      return true;
  return false;
}



bool Elem::contains_edge_of(const Elem * e) const
{
  unsigned int num_contained_edges = 0;

  // Our vertices are the first numbered nodes
  for (auto n : make_range(e->n_vertices()))
    {
      if (this->contains_point(e->point(n)))
        {
          num_contained_edges++;
          if (num_contained_edges>=2)
            {
              return true;
            }
        }
    }
  return false;
}



void Elem::find_point_neighbors(const Point & p,
                                std::set<const Elem *> & neighbor_set) const
{
  libmesh_assert(this->contains_point(p));
  libmesh_assert(this->active());

  neighbor_set.clear();
  neighbor_set.insert(this);

  std::set<const Elem *> untested_set, next_untested_set;
  untested_set.insert(this);

#ifdef LIBMESH_ENABLE_AMR
  std::vector<const Elem *> active_neighbor_children;
#endif // #ifdef LIBMESH_ENABLE_AMR

  while (!untested_set.empty())
    {
      // Loop over all the elements in the patch that haven't already
      // been tested
      for (const auto & elem : untested_set)
          for (auto current_neighbor : elem->neighbor_ptr_range())
            {
              if (current_neighbor &&
                  current_neighbor != remote_elem)    // we have a real neighbor on this side
                {
                  if (current_neighbor->active())                // ... if it is active
                    {
                      auto it = neighbor_set.lower_bound(current_neighbor);
                      if ((it == neighbor_set.end() || *it != current_neighbor) &&
                          current_neighbor->contains_point(p))   // ... don't have and touches p
                        {
                          // Add it and test it
                          next_untested_set.insert(current_neighbor);
                          neighbor_set.emplace_hint(it, current_neighbor);
                        }
                    }
#ifdef LIBMESH_ENABLE_AMR
                  else                                 // ... the neighbor is *not* active,
                    {                                  // ... so add *all* neighboring
                                                       // active children that touch p
                      active_neighbor_children.clear();
                      current_neighbor->active_family_tree_by_neighbor
                        (active_neighbor_children, elem);

                      for (const auto & current_child : active_neighbor_children)
                        {
                          auto it = neighbor_set.lower_bound(current_child);
                          if ((it == neighbor_set.end() || *it != current_child) &&
                              current_child->contains_point(p))
                            {
                              // Add it and test it
                              next_untested_set.insert(current_child);
                              neighbor_set.emplace_hint(it, current_child);
                            }
                        }
                    }
#endif // #ifdef LIBMESH_ENABLE_AMR
                }
            }
      untested_set.swap(next_untested_set);
      next_untested_set.clear();
    }
}



void Elem::find_point_neighbors(std::set<const Elem *> & neighbor_set) const
{
  this->find_point_neighbors(neighbor_set, this);
}



void Elem::find_point_neighbors(std::set<const Elem *> & neighbor_set,
                                const Elem * start_elem) const
{
  ElemInternal::find_point_neighbors(this, neighbor_set, start_elem);
}



void Elem::find_point_neighbors(std::set<Elem *> & neighbor_set,
                                Elem * start_elem)
{
  ElemInternal::find_point_neighbors(this, neighbor_set, start_elem);
}



void Elem::find_edge_neighbors(const Point & p1,
                               const Point & p2,
                               std::set<const Elem *> & neighbor_set) const
{
  // Simple but perhaps suboptimal code: find elements containing the
  // first point, then winnow this set down by removing elements which
  // don't also contain the second point

  libmesh_assert(this->contains_point(p2));
  this->find_point_neighbors(p1, neighbor_set);

  std::set<const Elem *>::iterator        it = neighbor_set.begin();
  const std::set<const Elem *>::iterator end = neighbor_set.end();

  while (it != end)
    {
      // As of C++11, set::erase returns an iterator to the element
      // following the erased element, or end.
      if (!(*it)->contains_point(p2))
        it = neighbor_set.erase(it);
      else
        ++it;
    }
}



void Elem::find_edge_neighbors(std::set<const Elem *> & neighbor_set) const
{
  neighbor_set.clear();
  neighbor_set.insert(this);

  std::set<const Elem *> untested_set, next_untested_set;
  untested_set.insert(this);

  while (!untested_set.empty())
    {
      // Loop over all the elements in the patch that haven't already
      // been tested
      for (const auto & elem : untested_set)
        {
          for (auto current_neighbor : elem->neighbor_ptr_range())
            {
              if (current_neighbor &&
                  current_neighbor != remote_elem)    // we have a real neighbor on this side
                {
                  if (current_neighbor->active())                // ... if it is active
                    {
                      if (this->contains_edge_of(current_neighbor) // ... and touches us
                          || current_neighbor->contains_edge_of(this))
                        {
                          // Make sure we'll test it
                          if (!neighbor_set.count(current_neighbor))
                            next_untested_set.insert (current_neighbor);

                          // And add it
                          neighbor_set.insert (current_neighbor);
                        }
                    }
#ifdef LIBMESH_ENABLE_AMR
                  else                                 // ... the neighbor is *not* active,
                    {                                  // ... so add *all* neighboring
                                                       // active children
                      std::vector<const Elem *> active_neighbor_children;

                      current_neighbor->active_family_tree_by_neighbor
                        (active_neighbor_children, elem);

                      for (const auto & current_child : active_neighbor_children)
                        if (this->contains_edge_of(current_child) || current_child->contains_edge_of(this))
                          {
                            // Make sure we'll test it
                            if (!neighbor_set.count(current_child))
                              next_untested_set.insert (current_child);

                            neighbor_set.insert (current_child);
                          }
                    }
#endif // #ifdef LIBMESH_ENABLE_AMR
                }
            }
        }
      untested_set.swap(next_untested_set);
      next_untested_set.clear();
    }
}



void Elem::find_interior_neighbors(std::set<const Elem *> & neighbor_set) const
{
  ElemInternal::find_interior_neighbors(this, neighbor_set);
}



void Elem::find_interior_neighbors(std::set<Elem *> & neighbor_set)
{
  ElemInternal::find_interior_neighbors(this, neighbor_set);
}



const Elem * Elem::interior_parent () const
{
  // interior parents make no sense for full-dimensional elements.
  if (this->dim() >= LIBMESH_DIM)
      return nullptr;

  // they USED TO BE only good for level-0 elements, but we now
  // support keeping interior_parent() valid on refined boundary
  // elements.
  // if (this->level() != 0)
  // return this->parent()->interior_parent();

  // We store the interior_parent pointer after both the parent
  // neighbor and neighbor pointers
  Elem * interior_p = _elemlinks[1+this->n_sides()];

  // If we have an interior_parent, we USED TO assume it was a
  // one-higher-dimensional interior element, but we now allow e.g.
  // edge elements to have a 3D interior_parent with no
  // intermediate 2D element.
  // libmesh_assert (!interior_p ||
  //                interior_p->dim() == (this->dim()+1));
  libmesh_assert (!interior_p ||
                  (interior_p == remote_elem) ||
                  (interior_p->dim() > this->dim()));

  // If an element in a multi-dimensional mesh has an interior_parent
  // link, it should be at our level or coarser, just like a neighbor
  // link.  Our collect_families() code relies on this, but it might
  // be tempting for users to manually assign something that breaks
  // it.
  //
  // However, we *also* create temporary side elements, and we don't
  // bother with creating ancestors for those, so they can be at level
  // 0 even when they're sides of non-level-0 elements.
  libmesh_assert (!interior_p ||
                  (interior_p->level() <= this->level()) ||
                  (this->level() == 0 &&
                   this->id() == DofObject::invalid_id));

  return interior_p;
}



Elem * Elem::interior_parent ()
{
  // See the const version for comments
  if (this->dim() >= LIBMESH_DIM)
      return nullptr;

  Elem * interior_p = _elemlinks[1+this->n_sides()];

  libmesh_assert (!interior_p ||
                  (interior_p == remote_elem) ||
                  (interior_p->dim() > this->dim()));

  return interior_p;
}



void Elem::set_interior_parent (Elem * p)
{
  // interior parents make no sense for full-dimensional elements.
  libmesh_assert (!p ||
                  this->dim() < LIBMESH_DIM);

  // If we have an interior_parent, we USED TO assume it was a
  // one-higher-dimensional interior element, but we now allow e.g.
  // edge elements to have a 3D interior_parent with no
  // intermediate 2D element.
  // libmesh_assert (!p ||
  //                 p->dim() == (this->dim()+1));
  libmesh_assert (!p ||
                  (p == remote_elem) ||
                  (p->dim() > this->dim()));

  _elemlinks[1+this->n_sides()] = p;
}



#ifdef LIBMESH_ENABLE_PERIODIC

Elem * Elem::topological_neighbor (const unsigned int i,
                                   MeshBase & mesh,
                                   const PointLocatorBase & point_locator,
                                   const PeriodicBoundaries * pb)
{
  libmesh_assert_less (i, this->n_neighbors());

  Elem * neighbor_i = this->neighbor_ptr(i);
  if (neighbor_i != nullptr)
    return neighbor_i;

  if (pb)
    {
      // Since the neighbor is nullptr it must be on a boundary. We need
      // see if this is a periodic boundary in which case it will have a
      // topological neighbor
      std::vector<boundary_id_type> bc_ids;
      mesh.get_boundary_info().boundary_ids(this, cast_int<unsigned short>(i), bc_ids);
      for (const auto & id : bc_ids)
        if (pb->boundary(id))
          {
            // Since the point locator inside of periodic boundaries
            // returns a const pointer we will retrieve the proper
            // pointer directly from the mesh object.
            const Elem * const cn = pb->neighbor(id, point_locator, this, i);
            neighbor_i = const_cast<Elem *>(cn);

            // Since coarse elements do not have more refined
            // neighbors we need to make sure that we don't return one
            // of these types of neighbors.
            if (neighbor_i)
              while (level() < neighbor_i->level())
                neighbor_i = neighbor_i->parent();
            return neighbor_i;
          }
    }

  return nullptr;
}



const Elem * Elem::topological_neighbor (const unsigned int i,
                                         const MeshBase & mesh,
                                         const PointLocatorBase & point_locator,
                                         const PeriodicBoundaries * pb) const
{
  libmesh_assert_less (i, this->n_neighbors());

  const Elem * neighbor_i = this->neighbor_ptr(i);
  if (neighbor_i != nullptr)
    return neighbor_i;

  if (pb)
    {
      // Since the neighbor is nullptr it must be on a boundary. We need
      // see if this is a periodic boundary in which case it will have a
      // topological neighbor
      std::vector<boundary_id_type> bc_ids;
      mesh.get_boundary_info().boundary_ids(this, cast_int<unsigned short>(i), bc_ids);
      for (const auto & id : bc_ids)
        if (pb->boundary(id))
          {
            neighbor_i = pb->neighbor(id, point_locator, this, i);

            // Since coarse elements do not have more refined
            // neighbors we need to make sure that we don't return one
            // of these types of neighbors.
            if (neighbor_i)
              while (level() < neighbor_i->level())
                neighbor_i = neighbor_i->parent();
            return neighbor_i;
          }
    }

  return nullptr;
}


bool Elem::has_topological_neighbor (const Elem * elem,
                                     const MeshBase & mesh,
                                     const PointLocatorBase & point_locator,
                                     const PeriodicBoundaries * pb) const
{
  // First see if this is a normal "interior" neighbor
  if (has_neighbor(elem))
    return true;

  for (auto n : this->side_index_range())
    if (this->topological_neighbor(n, mesh, point_locator, pb))
      return true;

  return false;
}


#endif

#ifndef NDEBUG

void Elem::libmesh_assert_valid_node_pointers() const
{
  libmesh_assert(this->valid_id());
  for (auto n : this->node_index_range())
    {
      libmesh_assert(this->node_ptr(n));
      libmesh_assert(this->node_ptr(n)->valid_id());
    }
}



void Elem::libmesh_assert_valid_neighbors() const
{
  for (auto n : this->side_index_range())
    {
      const Elem * neigh = this->neighbor_ptr(n);

      // Any element might have a remote neighbor; checking
      // to make sure that's not inaccurate is tough.
      if (neigh == remote_elem)
        continue;

      if (neigh)
        {
          // Only subactive elements have subactive neighbors
          libmesh_assert (this->subactive() || !neigh->subactive());

          const Elem * elem = this;

          // If we're subactive but our neighbor isn't, its
          // return neighbor link will be to our first active
          // ancestor OR to our inactive ancestor of the same
          // level as neigh,
          if (this->subactive() && !neigh->subactive())
            {
              for (elem = this; !elem->active();
                   elem = elem->parent())
                libmesh_assert(elem);
            }
          else
            {
              unsigned int rev = neigh->which_neighbor_am_i(elem);
              libmesh_assert_less (rev, neigh->n_neighbors());

              if (this->subactive() && !neigh->subactive())
                {
                  while (neigh->neighbor_ptr(rev) != elem)
                    {
                      libmesh_assert(elem->parent());
                      elem = elem->parent();
                    }
                }
              else
                {
                  const Elem * nn = neigh->neighbor_ptr(rev);
                  libmesh_assert(nn);

                  for (; elem != nn; elem = elem->parent())
                    libmesh_assert(elem);
                }
            }
        }
      // If we don't have a neighbor and we're not subactive, our
      // ancestors shouldn't have any neighbors in this same
      // direction.
      else if (!this->subactive())
        {
          const Elem * my_parent = this->parent();
          if (my_parent &&
              // A parent with a different dimension isn't really one of
              // our ancestors, it means we're on a boundary mesh and this
              // is an interior mesh element for which we're on a side.
              // Nothing to test for in that case.
              (my_parent->dim() == this->dim()))
            libmesh_assert (!my_parent->neighbor_ptr(n));
        }
    }
}

#endif // !NDEBUG



void Elem::make_links_to_me_local(unsigned int n, unsigned int nn)
{
  Elem * neigh = this->neighbor_ptr(n);

  // Don't bother calling this function unless it's necessary
  libmesh_assert(neigh);
  libmesh_assert(!neigh->is_remote());

  // We never have neighbors more refined than us
  libmesh_assert_less_equal (neigh->level(), this->level());

  // We never have subactive neighbors of non subactive elements
  libmesh_assert(!neigh->subactive() || this->subactive());

  // If we have a neighbor less refined than us then it must not
  // have any more refined descendants we could have pointed to
  // instead.
  libmesh_assert((neigh->level() == this->level()) ||
                 (neigh->active() && !this->subactive()) ||
                 (!neigh->has_children() && this->subactive()));

  // If neigh is at our level, then its family might have
  // remote_elem neighbor links which need to point to us
  // instead, but if not, then we're done.
  if (neigh->level() != this->level())
    return;

  // What side of neigh are we on?  nn.
  //
  // We can't use the usual Elem method because we're in the middle of
  // restoring topology.  We can't compare side_ptr nodes because
  // users want to abuse neighbor_ptr to point to
  // not-technically-neighbors across mesh slits.  We can't compare
  // node locations because users want to move those
  // not-technically-neighbors until they're
  // not-even-geometrically-neighbors.

  // Find any elements that ought to point to elem
  std::vector<Elem *> neigh_family;
#ifdef LIBMESH_ENABLE_AMR
  if (this->active())
    neigh->family_tree_by_side(neigh_family, nn);
  else
#endif
    neigh_family.push_back(neigh);

  // And point them to elem
  for (auto & neigh_family_member : neigh_family)
    {
      // Only subactive elements point to other subactive elements
      if (this->subactive() && !neigh_family_member->subactive())
        continue;

      // Ideally, the neighbor link ought to either be correct
      // already or ought to be to remote_elem.
      //
      // However, if we're redistributing a newly created elem,
      // after an AMR step but before find_neighbors has fixed up
      // neighbor links, we might have an out of date neighbor
      // link to elem's parent instead.
#ifdef LIBMESH_ENABLE_AMR
      libmesh_assert((neigh_family_member->neighbor_ptr(nn) &&
                      (neigh_family_member->neighbor_ptr(nn)->active() ||
                       neigh_family_member->neighbor_ptr(nn)->is_ancestor_of(this))) ||
                     (neigh_family_member->neighbor_ptr(nn) == remote_elem) ||
                     ((this->refinement_flag() == JUST_REFINED) &&
                      (this->parent() != nullptr) &&
                      (neigh_family_member->neighbor_ptr(nn) == this->parent())));
#else
      libmesh_assert((neigh_family_member->neighbor_ptr(nn) == this) ||
                     (neigh_family_member->neighbor_ptr(nn) == remote_elem));
#endif

      neigh_family_member->set_neighbor(nn, this);
    }
}


void Elem::make_links_to_me_remote()
{
  libmesh_assert_not_equal_to (this, remote_elem);

  // We need to have handled any children first
#if defined(LIBMESH_ENABLE_AMR) && defined(DEBUG)
  if (this->has_children())
    for (auto & child : this->child_ref_range())
      libmesh_assert_equal_to (&child, remote_elem);
#endif

  // Remotify any neighbor links
  for (auto neigh : this->neighbor_ptr_range())
    {
      if (neigh && neigh != remote_elem)
        {
          // My neighbor should never be more refined than me; my real
          // neighbor would have been its parent in that case.
          libmesh_assert_greater_equal (this->level(), neigh->level());

          if (this->level() == neigh->level() &&
              neigh->has_neighbor(this))
            {
#ifdef LIBMESH_ENABLE_AMR
              // My neighbor may have descendants which also consider me a
              // neighbor
              std::vector<Elem *> family;
              neigh->total_family_tree_by_neighbor (family, this);

              // FIXME - There's a lot of ugly const_casts here; we
              // may want to make remote_elem non-const
              for (auto & n : family)
                {
                  libmesh_assert (n);
                  if (n == remote_elem)
                    continue;
                  unsigned int my_s = n->which_neighbor_am_i(this);
                  libmesh_assert_less (my_s, n->n_neighbors());
                  libmesh_assert_equal_to (n->neighbor_ptr(my_s), this);
                  n->set_neighbor(my_s, const_cast<RemoteElem *>(remote_elem));
                }
#else
              unsigned int my_s = neigh->which_neighbor_am_i(this);
              libmesh_assert_less (my_s, neigh->n_neighbors());
              libmesh_assert_equal_to (neigh->neighbor_ptr(my_s), this);
              neigh->set_neighbor(my_s, const_cast<RemoteElem *>(remote_elem));
#endif
            }
#ifdef LIBMESH_ENABLE_AMR
          // Even if my neighbor doesn't link back to me, it might
          // have subactive descendants which do
          else if (neigh->has_children())
            {
              // If my neighbor at the same level doesn't have me as a
              // neighbor, I must be subactive
              libmesh_assert(this->level() > neigh->level() ||
                             this->subactive());

              // My neighbor must have some ancestor of mine as a
              // neighbor
              Elem * my_ancestor = this->parent();
              libmesh_assert(my_ancestor);
              while (!neigh->has_neighbor(my_ancestor))
                {
                  my_ancestor = my_ancestor->parent();
                  libmesh_assert(my_ancestor);
                }

              // My neighbor may have descendants which consider me a
              // neighbor
              std::vector<Elem *> family;
              neigh->total_family_tree_by_subneighbor (family, my_ancestor, this);

              for (auto & n : family)
                {
                  libmesh_assert (n);
                  if (n->is_remote())
                    continue;
                  unsigned int my_s = n->which_neighbor_am_i(this);
                  libmesh_assert_less (my_s, n->n_neighbors());
                  libmesh_assert_equal_to (n->neighbor_ptr(my_s), this);
                  // TODO: we may want to make remote_elem non-const.
                  n->set_neighbor(my_s, const_cast<RemoteElem *>(remote_elem));
                }
            }
#endif
        }
    }

#ifdef LIBMESH_ENABLE_AMR
  // Remotify parent's child link
  Elem * my_parent = this->parent();
  if (my_parent &&
      // As long as it's not already remote
      my_parent != remote_elem &&
      // And it's a real parent, not an interior parent
      this->dim() == my_parent->dim())
    {
      unsigned int me = my_parent->which_child_am_i(this);
      libmesh_assert_equal_to (my_parent->child_ptr(me), this);
      my_parent->set_child(me, const_cast<RemoteElem *>(remote_elem));
    }
#endif
}


void Elem::remove_links_to_me()
{
  libmesh_assert_not_equal_to (this, remote_elem);

  // We need to have handled any children first
#ifdef LIBMESH_ENABLE_AMR
  libmesh_assert (!this->has_children());
#endif

  // Nullify any neighbor links
  for (auto neigh : this->neighbor_ptr_range())
    {
      if (neigh && neigh != remote_elem)
        {
          // My neighbor should never be more refined than me; my real
          // neighbor would have been its parent in that case.
          libmesh_assert_greater_equal (this->level(), neigh->level());

          if (this->level() == neigh->level() &&
              neigh->has_neighbor(this))
            {
#ifdef LIBMESH_ENABLE_AMR
              // My neighbor may have descendants which also consider me a
              // neighbor
              std::vector<Elem *> family;
              neigh->total_family_tree_by_neighbor (family, this);

              for (auto & n : family)
                {
                  libmesh_assert (n);
                  if (n->is_remote())
                    continue;
                  unsigned int my_s = n->which_neighbor_am_i(this);
                  libmesh_assert_less (my_s, n->n_neighbors());
                  libmesh_assert_equal_to (n->neighbor_ptr(my_s), this);
                  n->set_neighbor(my_s, nullptr);
                }
#else
              unsigned int my_s = neigh->which_neighbor_am_i(this);
              libmesh_assert_less (my_s, neigh->n_neighbors());
              libmesh_assert_equal_to (neigh->neighbor_ptr(my_s), this);
              neigh->set_neighbor(my_s, nullptr);
#endif
            }
#ifdef LIBMESH_ENABLE_AMR
          // Even if my neighbor doesn't link back to me, it might
          // have subactive descendants which do
          else if (neigh->has_children())
            {
              // If my neighbor at the same level doesn't have me as a
              // neighbor, I must be subactive
              libmesh_assert(this->level() > neigh->level() ||
                             this->subactive());

              // My neighbor must have some ancestor of mine as a
              // neighbor
              Elem * my_ancestor = this->parent();
              libmesh_assert(my_ancestor);
              while (!neigh->has_neighbor(my_ancestor))
                {
                  my_ancestor = my_ancestor->parent();
                  libmesh_assert(my_ancestor);
                }

              // My neighbor may have descendants which consider me a
              // neighbor
              std::vector<Elem *> family;
              neigh->total_family_tree_by_subneighbor (family, my_ancestor, this);

              for (auto & n : family)
                {
                  libmesh_assert (n);
                  if (n->is_remote())
                    continue;
                  unsigned int my_s = n->which_neighbor_am_i(this);
                  libmesh_assert_less (my_s, n->n_neighbors());
                  libmesh_assert_equal_to (n->neighbor_ptr(my_s), this);
                  n->set_neighbor(my_s, nullptr);
                }
            }
#endif
        }
    }

#ifdef LIBMESH_ENABLE_AMR
  // We can't currently delete a child with a parent!
  libmesh_assert (!this->parent());
#endif
}



void Elem::write_connectivity (std::ostream & out_stream,
                               const IOPackage iop) const
{
  libmesh_assert (out_stream.good());
  libmesh_assert(_nodes);
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        // This connectivity vector will be used repeatedly instead
        // of being reconstructed inside the loop.
        std::vector<dof_id_type> conn;
        for (auto sc : make_range(this->n_sub_elem()))
          {
            this->connectivity(sc, TECPLOT, conn);

            std::copy(conn.begin(),
                      conn.end(),
                      std::ostream_iterator<dof_id_type>(out_stream, " "));

            out_stream << '\n';
          }
        return;
      }

    case UCD:
      {
        for (auto i : this->node_index_range())
          out_stream << this->node_id(i)+1 << "\t";

        out_stream << '\n';
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



Real Elem::quality (const ElemQuality q) const
{
  switch (q)
    {
      // Return the maximum ratio of edge lengths, or zero if that
      // maximum would otherwise be infinity.
    case EDGE_LENGTH_RATIO:
      {
        if (this->dim() < 2)
          return 1;

        std::vector<Real> edge_lengths(this->n_edges());
        for (auto e : index_range(edge_lengths))
          edge_lengths[e] = (this->point(this->local_edge_node(e,1)) -
                             this->point(this->local_edge_node(e,0))).norm();

        auto [min, max] =
          std::minmax_element(edge_lengths.begin(), edge_lengths.end());

        if (*min == 0.)
          return 0.;
        else
          return *max / *min;
      }

    case MIN_ANGLE:
    case MAX_ANGLE:
      {
        // 1D elements don't have interior angles, so just return some
        // dummy value in that case.
        if (this->dim() < 2)
          return 0.;

        // Initialize return values
        Real min_angle = std::numeric_limits<Real>::max();
        Real max_angle = -std::numeric_limits<Real>::max();

        for (auto n : this->node_index_range())
          {
            // Get list of edge ids adjacent to this node.
            const auto adjacent_edge_ids = this->edges_adjacent_to_node(n);

            // Skip any nodes with fewer than 2 adjacent edges. You
            // need at least two adjacent edges to form an interior
            // element angle.
            auto N = adjacent_edge_ids.size();
            if (N < 2)
              continue;

            // Consider all possible pairs of edges adjacent to node n
            for (unsigned int first = 0; first < N-1; ++first)
              for (unsigned int second = first+1; second < N; ++second)
                {
                  // Get ids of first and second edges
                  auto first_edge = adjacent_edge_ids[first];
                  auto second_edge = adjacent_edge_ids[second];

                  // Get node ids of first and second edge
                  auto first_edge_node_0 = this->local_edge_node(first_edge, 0);
                  auto first_edge_node_1 = this->local_edge_node(first_edge, 1);
                  auto second_edge_node_0 = this->local_edge_node(second_edge, 0);
                  auto second_edge_node_1 = this->local_edge_node(second_edge, 1);

                  // Orient both edges so that edge node 0 == n
                  if (first_edge_node_0 != n)
                    std::swap(first_edge_node_0, first_edge_node_1);
                  if (second_edge_node_0 != n)
                    std::swap(second_edge_node_0, second_edge_node_1);

                  libmesh_assert_equal_to(first_edge_node_0, n);
                  libmesh_assert_equal_to(second_edge_node_0, n);

                  // Locally oriented edge vectors
                  Point
                    first_ev = this->point(first_edge_node_1) - this->point(first_edge_node_0),
                    second_ev = this->point(second_edge_node_1) - this->point(second_edge_node_0);

                  // Angle between them in the range [0, pi]
                  Real theta = std::acos(first_ev.unit() * second_ev.unit());

                  // Track min and max angles seen
                  min_angle = std::min(theta, min_angle);
                  max_angle = std::max(theta, max_angle);
                }
          }

        // Return requested extreme value (in degrees)
        return Real(180)/libMesh::pi * ((q == MIN_ANGLE) ? min_angle : max_angle);
      }

    case MIN_DIHEDRAL_ANGLE:
    case MAX_DIHEDRAL_ANGLE:
      {
        // For 1D and 2D elements, just call this function with MIN,MAX_ANGLE instead.
        if (this->dim() < 3)
          return this->quality((q == MIN_DIHEDRAL_ANGLE) ? MIN_ANGLE : MAX_ANGLE);

        // Initialize return values
        Real min_angle = std::numeric_limits<Real>::max();
        Real max_angle = -std::numeric_limits<Real>::max();

        // Algorithm for computing dihedral angles:
        // .) Loop over the edges, use the edge_sides_map to get the
        //    two sides adjacent to the edge.
        // .) For each adjacent side, use the side_nodes_map to get the
        //    list of three (or more) node ids on that side.
        // .) Use the node ids to compute the (approximate) inward
        //    normal for the side.  Note: this is approximate since 3D
        //    elements with four-sided faces and quadratic tetrahedra
        //    do not necessarily have planar faces.
        // .) Compute the dihedral angle for the current edge, compare
        //    it to the current min and max dihedral angles.
        for (auto e : this->edge_index_range())
          {
            // Get list of edge ids adjacent to this node.
            const auto adjacent_side_ids = this->sides_on_edge(e);

            // All 3D elements should have exactly two sides adjacent to each edge.
            libmesh_assert_equal_to(adjacent_side_ids.size(), 2);

            // Get lists of node ids on each side
            const auto side_0_node_ids = this->nodes_on_side(adjacent_side_ids[0]);
            const auto side_1_node_ids = this->nodes_on_side(adjacent_side_ids[1]);

            // All 3D elements have at least three nodes on each side
            libmesh_assert_greater_equal(side_0_node_ids.size(), 3);
            libmesh_assert_greater_equal(side_1_node_ids.size(), 3);

            // Construct (approximate) inward normal on each adjacent side
            const auto side_0_normal =
              (this->point(side_0_node_ids[2]) - this->point(side_0_node_ids[0])).cross
              (this->point(side_0_node_ids[1]) - this->point(side_0_node_ids[0])).unit();
            const auto side_1_normal =
              (this->point(side_1_node_ids[2]) - this->point(side_1_node_ids[0])).cross
              (this->point(side_1_node_ids[1]) - this->point(side_1_node_ids[0])).unit();

            // Compute dihedral angle between the planes.
            // Using the absolute value ensures that we get the same result
            // even if the orientation of one of the planes is flipped, i.e.
            // it always gives us a value in the range [0, pi/2].
            // https://en.wikipedia.org/wiki/Dihedral_angle
            Real theta = std::acos(std::abs(side_0_normal * side_1_normal));

            // Track min and max angles seen
            min_angle = std::min(theta, min_angle);
            max_angle = std::max(theta, max_angle);
          }

        // Return requested extreme value (in degrees)
        return Real(180)/libMesh::pi * ((q == MIN_DIHEDRAL_ANGLE) ? min_angle : max_angle);
      }

    case JACOBIAN:
    case SCALED_JACOBIAN:
      {
        // 1D elements don't have interior corners, so this metric
        // does not really apply to them.
        const auto N = this->dim();
        if (N < 2)
          return 1.;

        // Initialize return value
        Real min_node_area = std::numeric_limits<Real>::max();

        for (auto n : this->node_index_range())
          {
            // Get list of edge ids adjacent to this node.
            auto adjacent_edge_ids = this->edges_adjacent_to_node(n);

            // Skip any nodes that don't have dim() adjacent edges. In 2D,
            // you need at least two adjacent edges to compute the cross
            // product, and in 3D, you need at least three adjacent edges
            // to compute the scalar triple product. The only element
            // type which has an unusual topology in this regard is the
            // Pyramid element in 3D, where 4 edges meet at the apex node.
            // For now, we just skip this node when computing the JACOBIAN
            // metric for Pyramids.
            if (adjacent_edge_ids.size() != N)
              continue;

            // Construct oriented edges
            std::vector<Point> oriented_edges(N);
            for (auto i : make_range(N))
              {
                auto node_0 = this->local_edge_node(adjacent_edge_ids[i], 0);
                auto node_1 = this->local_edge_node(adjacent_edge_ids[i], 1);
                if (node_0 != n)
                  std::swap(node_0, node_1);
                oriented_edges[i] = this->point(node_1) - this->point(node_0);
              }

            // Compute unscaled area (2D) or volume (3D) using the
            // cross product (2D) or scalar triple product (3D) of the
            // oriented edges. We take the absolute value so that we
            // don't have to worry about the sign of the area.
            Real node_area = (N == 2) ?
              cross_norm(oriented_edges[0], oriented_edges[1]) :
              std::abs(triple_product(oriented_edges[0], oriented_edges[1], oriented_edges[2]));

            // Divide by (non-zero) edge lengths if computing scaled Jacobian.
            // If any length is zero, then set node_area to zero. This means that
            // e.g. degenerate quadrilaterals will have a zero SCALED_JACOBIAN metric.
            if (q == SCALED_JACOBIAN)
              for (auto i : make_range(N))
                {
                  Real len_i = oriented_edges[i].norm();
                  node_area = (len_i == 0.) ? 0. : (node_area / len_i);
                }

            // Update minimum
            min_node_area = std::min(node_area, min_node_area);
          }

        return min_node_area;
      }

      // Return 1 if we made it here
    default:
      {
        libmesh_do_once( libmesh_here();

                         libMesh::err << "ERROR: quality metric "
                         << Utility::enum_to_string(q)
                         << " not implemented on element type "
                         << Utility::enum_to_string(this->type())
                         << std::endl
                         << "Returning 1."
                         << std::endl; );

        return 1.;
      }
    }
}



bool Elem::ancestor() const
{
#ifdef LIBMESH_ENABLE_AMR

  // Use a fast, DistributedMesh-safe definition
  const bool is_ancestor =
    !this->active() && !this->subactive();

  // But check for inconsistencies if we have time
#ifdef DEBUG
  if (!is_ancestor && this->has_children())
    {
      for (auto & c : this->child_ref_range())
        {
          if (&c != remote_elem)
            {
              libmesh_assert(!c.active());
              libmesh_assert(!c.ancestor());
            }
        }
    }
#endif // DEBUG

  return is_ancestor;

#else
  return false;
#endif
}



#ifdef LIBMESH_ENABLE_AMR

void Elem::add_child (Elem * elem)
{
  const unsigned int nc = this->n_children();

  if (!_children)
    {
      _children = std::make_unique<Elem *[]>(nc);

      for (unsigned int c = 0; c != nc; c++)
        this->set_child(c, nullptr);
    }

  for (unsigned int c = 0; c != nc; c++)
    {
      if (this->_children[c] == nullptr || this->_children[c] == remote_elem)
        {
          libmesh_assert_equal_to (this, elem->parent());
          this->set_child(c, elem);
          return;
        }
    }

  libmesh_error_msg("Error: Tried to add a child to an element with full children array");
}



void Elem::add_child (Elem * elem, unsigned int c)
{
  if (!this->has_children())
    {
      const unsigned int nc = this->n_children();
      _children = std::make_unique<Elem *[]>(nc);

      for (unsigned int i = 0; i != nc; i++)
        this->set_child(i, nullptr);
    }

  libmesh_assert (this->_children[c] == nullptr || this->child_ptr(c) == remote_elem);
  libmesh_assert (elem == remote_elem || this == elem->parent());

  this->set_child(c, elem);
}



void Elem::replace_child (Elem * elem, unsigned int c)
{
  libmesh_assert(this->has_children());

  libmesh_assert(this->child_ptr(c));

  this->set_child(c, elem);
}



void Elem::family_tree (std::vector<const Elem *> & family,
                          bool reset) const
{
  ElemInternal::family_tree(this, family, reset);
}



void Elem::family_tree (std::vector<Elem *> & family,
                        bool reset)
{
  ElemInternal::family_tree(this, family, reset);
}



void Elem::total_family_tree (std::vector<const Elem *> & family,
                              bool reset) const
{
  ElemInternal::total_family_tree(this, family, reset);
}



void Elem::total_family_tree (std::vector<Elem *> & family,
                              bool reset)
{
  ElemInternal::total_family_tree(this, family, reset);
}



void Elem::active_family_tree (std::vector<const Elem *> & active_family,
                               bool reset) const
{
  ElemInternal::active_family_tree(this, active_family, reset);
}



void Elem::active_family_tree (std::vector<Elem *> & active_family,
                               bool reset)
{
  ElemInternal::active_family_tree(this, active_family, reset);
}



void Elem::family_tree_by_side (std::vector<const Elem *> & family,
                                unsigned int side,
                                bool reset) const
{
  ElemInternal::family_tree_by_side(this, family, side, reset);
}



void Elem:: family_tree_by_side (std::vector<Elem *> & family,
                                 unsigned int side,
                                 bool reset)
{
  ElemInternal::family_tree_by_side(this, family, side, reset);
}



void Elem::total_family_tree_by_side (std::vector<const Elem *> & family,
                                      unsigned int side,
                                      bool reset) const
{
  ElemInternal::total_family_tree_by_side(this, family, side, reset);
}



void Elem::total_family_tree_by_side (std::vector<Elem *> & family,
                                      unsigned int side,
                                      bool reset)
{
  ElemInternal::total_family_tree_by_side(this, family, side, reset);
}



void Elem::active_family_tree_by_side (std::vector<const Elem *> & family,
                                       unsigned int side,
                                       bool reset) const
{
  ElemInternal::active_family_tree_by_side(this, family, side, reset);
}



void Elem::active_family_tree_by_side (std::vector<Elem *> & family,
                                       unsigned int side,
                                       bool reset)
{
  ElemInternal::active_family_tree_by_side(this, family, side, reset);
}



void Elem::family_tree_by_neighbor (std::vector<const Elem *> & family,
                                    const Elem * neighbor,
                                    bool reset) const
{
  ElemInternal::family_tree_by_neighbor(this, family, neighbor, reset);
}



void Elem::family_tree_by_neighbor (std::vector<Elem *> & family,
                                    Elem * neighbor,
                                    bool reset)
{
  ElemInternal::family_tree_by_neighbor(this, family, neighbor, reset);
}



void Elem::total_family_tree_by_neighbor (std::vector<const Elem *> & family,
                                          const Elem * neighbor,
                                          bool reset) const
{
  ElemInternal::total_family_tree_by_neighbor(this, family, neighbor, reset);
}



void Elem::total_family_tree_by_neighbor (std::vector<Elem *> & family,
                                          Elem * neighbor,
                                          bool reset)
{
  ElemInternal::total_family_tree_by_neighbor(this, family, neighbor, reset);
}



void Elem::family_tree_by_subneighbor (std::vector<const Elem *> & family,
                                       const Elem * neighbor,
                                       const Elem * subneighbor,
                                       bool reset) const
{
  ElemInternal::family_tree_by_subneighbor(this, family, neighbor, subneighbor, reset);
}



void Elem::family_tree_by_subneighbor (std::vector<Elem *> & family,
                                       Elem * neighbor,
                                       Elem * subneighbor,
                                       bool reset)
{
  ElemInternal::family_tree_by_subneighbor(this, family, neighbor, subneighbor, reset);
}



void Elem::total_family_tree_by_subneighbor (std::vector<const Elem *> & family,
                                             const Elem * neighbor,
                                             const Elem * subneighbor,
                                             bool reset) const
{
  ElemInternal::total_family_tree_by_subneighbor(this, family, neighbor, subneighbor, reset);
}



void Elem::total_family_tree_by_subneighbor (std::vector<Elem *> & family,
                                             Elem * neighbor,
                                             Elem * subneighbor,
                                             bool reset)
{
  ElemInternal::total_family_tree_by_subneighbor(this, family, neighbor, subneighbor, reset);
}



void Elem::active_family_tree_by_neighbor (std::vector<const Elem *> & family,
                                           const Elem * neighbor,
                                           bool reset) const
{
  ElemInternal::active_family_tree_by_neighbor(this, family, neighbor, reset);
}



void Elem::active_family_tree_by_neighbor (std::vector<Elem *> & family,
                                           Elem * neighbor,
                                           bool reset)
{
  ElemInternal::active_family_tree_by_neighbor(this, family, neighbor, reset);
}



void Elem::active_family_tree_by_topological_neighbor (std::vector<const Elem *> & family,
                                                       const Elem * neighbor,
                                                       const MeshBase & mesh,
                                                       const PointLocatorBase & point_locator,
                                                       const PeriodicBoundaries * pb,
                                                       bool reset) const
{
  ElemInternal::active_family_tree_by_topological_neighbor(this, family, neighbor,
                                                           mesh, point_locator, pb,
                                                           reset);
}



void Elem::active_family_tree_by_topological_neighbor (std::vector<Elem *> & family,
                                                       Elem * neighbor,
                                                       const MeshBase & mesh,
                                                       const PointLocatorBase & point_locator,
                                                       const PeriodicBoundaries * pb,
                                                       bool reset)
{
  ElemInternal::active_family_tree_by_topological_neighbor(this, family, neighbor,
                                                           mesh, point_locator, pb,
                                                           reset);
}


bool Elem::is_child_on_edge(const unsigned int c,
                            const unsigned int e) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (e, this->n_edges());

  std::unique_ptr<const Elem> my_edge = this->build_edge_ptr(e);
  std::unique_ptr<const Elem> child_edge = this->child_ptr(c)->build_edge_ptr(e);

  // We're assuming that an overlapping child edge has the same
  // number and orientation as its parent
  return (child_edge->node_id(0) == my_edge->node_id(0) ||
          child_edge->node_id(1) == my_edge->node_id(1));
}



unsigned int Elem::min_p_level_by_neighbor(const Elem * neighbor_in,
                                           unsigned int current_min) const
{
  libmesh_assert(!this->subactive());
  libmesh_assert(neighbor_in->active());

  // If we're an active element this is simple
  if (this->active())
    return std::min(current_min, this->p_level());

  libmesh_assert(has_neighbor(neighbor_in));

  // The p_level() of an ancestor element is already the minimum
  // p_level() of its children - so if that's high enough, we don't
  // need to examine any children.
  if (current_min <= this->p_level())
    return current_min;

  unsigned int min_p_level = current_min;

  for (auto & c : this->child_ref_range())
    if (&c != remote_elem && c.has_neighbor(neighbor_in))
      min_p_level =
        c.min_p_level_by_neighbor(neighbor_in, min_p_level);

  return min_p_level;
}


unsigned int Elem::min_new_p_level_by_neighbor(const Elem * neighbor_in,
                                               unsigned int current_min) const
{
  libmesh_assert(!this->subactive());
  libmesh_assert(neighbor_in->active());

  // If we're an active element this is simple
  if (this->active())
    {
      unsigned int new_p_level = this->p_level();
      if (this->p_refinement_flag() == Elem::REFINE)
        new_p_level += 1;
      if (this->p_refinement_flag() == Elem::COARSEN)
        {
          libmesh_assert_greater (new_p_level, 0);
          new_p_level -= 1;
        }
      return std::min(current_min, new_p_level);
    }

  libmesh_assert(has_neighbor(neighbor_in));

  unsigned int min_p_level = current_min;

  for (auto & c : this->child_ref_range())
    if (&c != remote_elem && c.has_neighbor(neighbor_in))
      min_p_level =
        c.min_new_p_level_by_neighbor(neighbor_in, min_p_level);

  return min_p_level;
}



unsigned int Elem::as_parent_node (unsigned int child,
                                   unsigned int child_node) const
{
  const unsigned int nc = this->n_children();
  libmesh_assert_less(child, nc);

  // Cached return values, indexed first by embedding_matrix version,
  // then by child number, then by child node number.
  std::vector<std::vector<std::vector<signed char>>> &
    cached_parent_indices = this->_get_parent_indices_cache();

  unsigned int em_vers = this->embedding_matrix_version();

  // We may be updating the cache on one thread, and while that
  // happens we can't safely access the cache from other threads.
  Threads::spin_mutex::scoped_lock lock(parent_indices_mutex);

  if (em_vers >= cached_parent_indices.size())
    cached_parent_indices.resize(em_vers+1);

  if (child >= cached_parent_indices[em_vers].size())
    {
      const signed char nn = cast_int<signed char>(this->n_nodes());

      cached_parent_indices[em_vers].resize(nc);

      for (unsigned int c = 0; c != nc; ++c)
        {
          const unsigned int ncn = this->n_nodes_in_child(c);
          cached_parent_indices[em_vers][c].resize(ncn);
          for (unsigned int cn = 0; cn != ncn; ++cn)
            {
              for (signed char n = 0; n != nn; ++n)
                {
                  const Real em_val = this->embedding_matrix
                    (c, cn, n);
                  if (em_val == 1)
                    {
                      cached_parent_indices[em_vers][c][cn] = n;
                      break;
                    }

                  if (em_val != 0)
                    {
                      cached_parent_indices[em_vers][c][cn] =
                        -1;
                      break;
                    }

                  // We should never see an all-zero embedding matrix
                  // row
                  libmesh_assert_not_equal_to (n+1, nn);
                }
            }
        }
    }

  const signed char cache_val =
    cached_parent_indices[em_vers][child][child_node];
  if (cache_val == -1)
    return libMesh::invalid_uint;

  return cached_parent_indices[em_vers][child][child_node];
}



const std::vector<std::pair<unsigned char, unsigned char>> &
Elem::parent_bracketing_nodes(unsigned int child,
                              unsigned int child_node) const
{
  // Indexed first by embedding matrix type, then by child id, then by
  // child node, then by bracketing pair
  std::vector<std::vector<std::vector<std::vector<std::pair<unsigned char, unsigned char>>>>> &
    cached_bracketing_nodes = this->_get_bracketing_node_cache();

  const unsigned int em_vers = this->embedding_matrix_version();

  // We may be updating the cache on one thread, and while that
  // happens we can't safely access the cache from other threads.
  Threads::spin_mutex::scoped_lock lock(parent_bracketing_nodes_mutex);

  if (cached_bracketing_nodes.size() <= em_vers)
    cached_bracketing_nodes.resize(em_vers+1);

  const unsigned int nc = this->n_children();

  // If we haven't cached the bracketing nodes corresponding to this
  // embedding matrix yet, let's do so now.
  if (cached_bracketing_nodes[em_vers].size() < nc)
    {
      // If we're a second-order element but we're not a full-order
      // element, then some of our bracketing nodes may not exist
      // except on the equivalent full-order element.  Let's build an
      // equivalent full-order element and make a copy of its cache to
      // use.
      if (this->default_order() != FIRST &&
          second_order_equivalent_type(this->type(), /*full_ordered=*/ true) != this->type())
        {
          // Check that we really are the non-full-order type
          libmesh_assert_equal_to
            (second_order_equivalent_type (this->type(), false),
             this->type());

          // Build the full-order type
          ElemType full_type =
            second_order_equivalent_type(this->type(), /*full_ordered=*/ true);
          std::unique_ptr<Elem> full_elem = Elem::build(full_type);

          // This won't work for elements with multiple
          // embedding_matrix versions, but every such element is full
          // order anyways.
          libmesh_assert_equal_to(em_vers, 0);

          // Make sure its cache has been built.  We temporarily
          // release our mutex lock so that the inner call can
          // re-acquire it.
          lock.release();
          full_elem->parent_bracketing_nodes(0,0);

          // And then we need to lock again, so that if someone *else*
          // grabbed our lock before we did we don't risk accessing
          // cached_bracketing_nodes while they're working on it.
          // Threading is hard.
          lock.acquire(parent_bracketing_nodes_mutex);

          // Copy its cache
          cached_bracketing_nodes =
            full_elem->_get_bracketing_node_cache();

          // Now we don't need to build the cache ourselves.
          return cached_bracketing_nodes[em_vers][child][child_node];
        }

      cached_bracketing_nodes[em_vers].resize(nc);

      const unsigned int nn = this->n_nodes();

      // We have to examine each child
      for (unsigned int c = 0; c != nc; ++c)
        {
          const unsigned int ncn = this->n_nodes_in_child(c);

          cached_bracketing_nodes[em_vers][c].resize(ncn);

          // We have to examine each node in that child
          for (unsigned int n = 0; n != ncn; ++n)
            {
              // If this child node isn't a vertex or an infinite
              // child element's mid-infinite-edge node, then we need
              // to find bracketing nodes on the child.
              if (!this->is_vertex_on_child(c, n)
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                  && !this->is_mid_infinite_edge_node(n)
#endif
                  )
                {
                  // Use the embedding matrix to find the child node
                  // location in parent master element space
                  Point bracketed_pt;

                  for (unsigned int pn = 0; pn != nn; ++pn)
                    {
                      const Real em_val =
                        this->embedding_matrix(c,n,pn);

                      libmesh_assert_not_equal_to (em_val, 1);
                      if (em_val != 0.)
                        bracketed_pt.add_scaled(this->master_point(pn), em_val);
                    }

                  // Check each pair of nodes on the child which are
                  // also both parent nodes
                  for (unsigned int n1 = 0; n1 != ncn; ++n1)
                    {
                      if (n1 == n)
                        continue;

                      unsigned int parent_n1 =
                        this->as_parent_node(c,n1);

                      if (parent_n1 == libMesh::invalid_uint)
                        continue;

                      Point p1 = this->master_point(parent_n1);

                      for (unsigned int n2 = n1+1; n2 < nn; ++n2)
                        {
                          if (n2 == n)
                            continue;

                          unsigned int parent_n2 =
                            this->as_parent_node(c,n2);

                          if (parent_n2 == libMesh::invalid_uint)
                            continue;

                          Point p2 = this->master_point(parent_n2);

                          Point pmid = (p1 + p2)/2;

                          if (pmid == bracketed_pt)
                            {
                              cached_bracketing_nodes[em_vers][c][n].emplace_back(parent_n1, parent_n2);
                              break;
                            }
                          else
                            libmesh_assert(!pmid.absolute_fuzzy_equals(bracketed_pt));
                        }
                    }
                }
              // If this child node is a parent node, we need to
              // find bracketing nodes on the parent.
              else
                {
                  unsigned int parent_node = this->as_parent_node(c,n);

                  Point bracketed_pt;

                  // If we're not a parent node, use the embedding
                  // matrix to find the child node location in parent
                  // master element space
                  if (parent_node == libMesh::invalid_uint)
                    {
                      for (unsigned int pn = 0; pn != nn; ++pn)
                        {
                          const Real em_val =
                            this->embedding_matrix(c,n,pn);

                          libmesh_assert_not_equal_to (em_val, 1);
                          if (em_val != 0.)
                            bracketed_pt.add_scaled(this->master_point(pn), em_val);
                        }
                    }
                  // If we're a parent node then we need no arithmetic
                  else
                    bracketed_pt = this->master_point(parent_node);

                  for (unsigned int n1 = 0; n1 != nn; ++n1)
                    {
                      if (n1 == parent_node)
                        continue;

                      Point p1 = this->master_point(n1);

                      for (unsigned int n2 = n1+1; n2 < nn; ++n2)
                        {
                          if (n2 == parent_node)
                            continue;

                          Point pmid = (p1 + this->master_point(n2))/2;

                          if (pmid == bracketed_pt)
                            {
                              cached_bracketing_nodes[em_vers][c][n].emplace_back(n1, n2);
                              break;
                            }
                          else
                            libmesh_assert(!pmid.absolute_fuzzy_equals(bracketed_pt));
                        }
                    }
                }
            }
        }
    }

  return cached_bracketing_nodes[em_vers][child][child_node];
}


const std::vector<std::pair<dof_id_type, dof_id_type>>
Elem::bracketing_nodes(unsigned int child,
                       unsigned int child_node) const
{
  std::vector<std::pair<dof_id_type, dof_id_type>> returnval;

  const std::vector<std::pair<unsigned char, unsigned char>> & pbc =
    this->parent_bracketing_nodes(child,child_node);

  for (const auto & pb : pbc)
    {
      const unsigned short n_n = this->n_nodes();
      if (pb.first < n_n && pb.second < n_n)
        returnval.emplace_back(this->node_id(pb.first), this->node_id(pb.second));
      else
        {
          // We must be on a non-full-order higher order element...
          libmesh_assert_not_equal_to(this->default_order(), FIRST);
          libmesh_assert_not_equal_to
            (second_order_equivalent_type (this->type(), true),
             this->type());
          libmesh_assert_equal_to
            (second_order_equivalent_type (this->type(), false),
             this->type());

          // And that's a shame, because this is a nasty search:

          // Build the full-order type
          ElemType full_type =
            second_order_equivalent_type(this->type(), /*full_ordered=*/ true);
          std::unique_ptr<Elem> full_elem = Elem::build(full_type);

          dof_id_type pt1 = DofObject::invalid_id;
          dof_id_type pt2 = DofObject::invalid_id;

          // Find the bracketing nodes by figuring out what
          // already-created children will have them.

          // This only doesn't break horribly because we add children
          // and nodes in straightforward + hierarchical orders...
          for (unsigned int c=0; c <= child; ++c)
            {
              libmesh_assert(this->child_ptr(c));
              if (this->child_ptr(c) == remote_elem)
                continue;

              for (auto n : make_range(this->n_nodes_in_child(c)))
                {
                  if (c == child && n == child_node)
                    break;

                  if (pb.first == full_elem->as_parent_node(c,n))
                    {
                      // We should be consistent
                      if (pt1 != DofObject::invalid_id)
                        libmesh_assert_equal_to(pt1, this->child_ptr(c)->node_id(n));

                      pt1 = this->child_ptr(c)->node_id(n);
                    }

                  if (pb.second == full_elem->as_parent_node(c,n))
                    {
                      // We should be consistent
                      if (pt2 != DofObject::invalid_id)
                        libmesh_assert_equal_to(pt2, this->child_ptr(c)->node_id(n));

                      pt2 = this->child_ptr(c)->node_id(n);
                    }
                }
            }

          // We should *usually* find all bracketing nodes by the time
          // we query them (again, because of the child & node add
          // order)
          //
          // The exception is if we're a HEX20, in which case we will
          // find pairs of vertex nodes and edge nodes bracketing the
          // new central node but we *won't* find the pairs of face
          // nodes which we would have had on a HEX27.  In that case
          // we'll still have enough bracketing nodes for a
          // topological lookup, but we won't be able to make the
          // following assertions.
          if (this->type() != HEX20)
            {
              libmesh_assert_not_equal_to (pt1, DofObject::invalid_id);
              libmesh_assert_not_equal_to (pt2, DofObject::invalid_id);
            }

          if (pt1 != DofObject::invalid_id &&
              pt2 != DofObject::invalid_id)
            returnval.emplace_back(pt1, pt2);
        }
    }

  return returnval;
}
#endif // #ifdef LIBMESH_ENABLE_AMR




bool Elem::contains_point (const Point & p, Real tol) const
{
  // We currently allow the user to enlarge the bounding box by
  // providing a tol > TOLERANCE (so this routine is identical to
  // Elem::close_to_point()), but print a warning so that the
  // user can eventually switch his code over to calling close_to_point()
  // instead, which is intended to be used for this purpose.
  if (tol > TOLERANCE)
    {
      libmesh_do_once(libMesh::err
                      << "WARNING: Resizing bounding box to match user-specified tolerance!\n"
                      << "In the future, calls to Elem::contains_point() with tol > TOLERANCE\n"
                      << "will be more optimized, but should not be used\n"
                      << "to search for points 'close to' elements!\n"
                      << "Instead, use Elem::close_to_point() for this purpose.\n"
                      << std::endl;);
      return this->point_test(p, tol, tol);
    }
  else
    return this->point_test(p, TOLERANCE, tol);
}




bool Elem::close_to_point (const Point & p, Real tol) const
{
  // This test uses the user's passed-in tolerance for the
  // bounding box test as well, thereby allowing the routine to
  // find points which are not only "in" the element, but also
  // "nearby" to within some tolerance.
  return this->point_test(p, tol, tol);
}




bool Elem::point_test(const Point & p, Real box_tol, Real map_tol) const
{
  libmesh_assert_greater (box_tol, 0.);
  libmesh_assert_greater (map_tol, 0.);

  // This is a great optimization on first order elements, but it
  // could return false negatives on higher orders
  if (this->default_order() == FIRST)
    {
      // Check to make sure the element *could* contain this point, so we
      // can avoid an expensive inverse_map call if it doesn't.
      bool
#if LIBMESH_DIM > 2
        point_above_min_z = false,
        point_below_max_z = false,
#endif
#if LIBMESH_DIM > 1
        point_above_min_y = false,
        point_below_max_y = false,
#endif
        point_above_min_x = false,
        point_below_max_x = false;

      // For relative bounding box checks in physical space
      const Real my_hmax = this->hmax();

      for (auto & n : this->node_ref_range())
        {
          point_above_min_x = point_above_min_x || (n(0) - my_hmax*box_tol <= p(0));
          point_below_max_x = point_below_max_x || (n(0) + my_hmax*box_tol >= p(0));
#if LIBMESH_DIM > 1
          point_above_min_y = point_above_min_y || (n(1) - my_hmax*box_tol <= p(1));
          point_below_max_y = point_below_max_y || (n(1) + my_hmax*box_tol >= p(1));
#endif
#if LIBMESH_DIM > 2
          point_above_min_z = point_above_min_z || (n(2) - my_hmax*box_tol <= p(2));
          point_below_max_z = point_below_max_z || (n(2) + my_hmax*box_tol >= p(2));
#endif
        }

      if (
#if LIBMESH_DIM > 2
          !point_above_min_z ||
          !point_below_max_z ||
#endif
#if LIBMESH_DIM > 1
          !point_above_min_y ||
          !point_below_max_y ||
#endif
          !point_above_min_x ||
          !point_below_max_x)
        return false;
    }

  // To be on the safe side, we converge the inverse_map() iteration
  // to a slightly tighter tolerance than that requested by the
  // user...
  const Point mapped_point = FEMap::inverse_map(this->dim(),
                                                this,
                                                p,
                                                0.1*map_tol, // <- this is |dx| tolerance, the Newton residual should be ~ |dx|^2
                                                /*secure=*/ false,
                                                /*extra_checks=*/ false);

  // Check that the refspace point maps back to p!  This is only necessary
  // for 1D and 2D elements, 3D elements always live in 3D.
  //
  // TODO: The contains_point() function could most likely be implemented
  // more efficiently in the element sub-classes themselves, at least for
  // the linear element types.
  if (this->dim() < 3)
    {
      Point xyz = FEMap::map(this->dim(), this, mapped_point);

      // Compute the distance between the original point and the re-mapped point.
      // They should be in the same place.
      Real dist = (xyz - p).norm();


      // If dist is larger than some fraction of the tolerance, then return false.
      // This can happen when e.g. a 2D element is living in 3D, and
      // FEMap::inverse_map() maps p onto the projection of the element,
      // effectively "tricking" on_reference_element().
      if (dist > this->hmax() * map_tol)
        return false;
    }

  return this->on_reference_element(mapped_point, map_tol);
}




bool Elem::has_invertible_map(Real /*tol*/) const
{
  QNodal qnodal {this->dim()};
  FEMap fe_map;
  auto & jac = fe_map.get_jacobian();

  // We have a separate check for is_singular_node() below, so in this
  // case its "OK" to do nodal quadrature on pyramids.
  qnodal.allow_nodal_pyramid_quadrature = true;
  qnodal.init(*this);
  auto & qp = qnodal.get_points();
  libmesh_assert_equal_to(qp.size(), this->n_nodes());

  std::vector<Point> one_point(1);
  std::vector<Real> one_weight(1,1);
  for (auto i : index_range(qp))
    {
      if (this->is_singular_node(i))
        continue;

      one_point[0] = qp[i];

      fe_map.init_reference_to_physical_map(this->dim(), one_point, this);
      fe_map.compute_map(this->dim(), one_weight, this, false);

      if (jac[0] <= 0)
        return false;
    }

  return true;
}



void Elem::print_info (std::ostream & os) const
{
  os << this->get_info()
     << std::endl;
}



std::string Elem::get_info () const
{
  std::ostringstream oss;

  oss << "  Elem Information"                                      << '\n'
      << "   id()=";

  if (this->valid_id())
    oss << this->id();
  else
    oss << "invalid";

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  oss << ", unique_id()=";
  if (this->valid_unique_id())
    oss << this->unique_id();
  else
    oss << "invalid";
#endif

  oss << ", subdomain_id()=" << this->subdomain_id();
  oss << ", processor_id()=" << this->processor_id()               << '\n';

  oss << "   type()="    << Utility::enum_to_string(this->type())  << '\n'
      << "   dim()="     << this->dim()                            << '\n'
      << "   n_nodes()=" << this->n_nodes()                        << '\n';

  oss << "   mapping=" << Utility::enum_to_string(this->mapping_type()) << '\n';

  for (auto n : this->node_index_range())
    {
      oss << "    " << n << this->node_ref(n);
      if (this->mapping_type() == RATIONAL_BERNSTEIN_MAP)
        {
          const unsigned char datum_index = this->mapping_data();
          oss << "    weight=" <<
            this->node_ref(n).get_extra_datum<Real>(datum_index) << '\n';
        }
    }

  oss << "   n_sides()=" << this->n_sides()                        << '\n';

  for (auto s : this->side_index_range())
    {
      oss << "    neighbor(" << s << ")=";
      if (this->neighbor_ptr(s))
        oss << this->neighbor_ptr(s)->id() << '\n';
      else
        oss << "nullptr\n";
    }

  if (!this->infinite())
    {
    oss << "   hmin()=" << this->hmin()
        << ", hmax()=" << this->hmax()                             << '\n'
        << "   volume()=" << this->volume()                        << '\n';
    }
  oss << "   active()=" << this->active()
    << ", ancestor()=" << this->ancestor()
    << ", subactive()=" << this->subactive()
    << ", has_children()=" << this->has_children()               << '\n'
    << "   parent()=";
  if (this->parent())
    oss << this->parent()->id() << '\n';
  else
    oss << "nullptr\n";
  oss << "   level()=" << this->level()
      << ", p_level()=" << this->p_level()                         << '\n'
#ifdef LIBMESH_ENABLE_AMR
      << "   refinement_flag()=" << Utility::enum_to_string(this->refinement_flag())        << '\n'
      << "   p_refinement_flag()=" << Utility::enum_to_string(this->p_refinement_flag())    << '\n'
#endif
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
      << "   infinite()=" << this->infinite()    << '\n';
  if (this->infinite())
    oss << "   origin()=" << this->origin()    << '\n'
#endif
      ;

  oss << "   DoFs=";
  for (auto s : make_range(this->n_systems()))
    for (auto v : make_range(this->n_vars(s)))
      for (auto c : make_range(this->n_comp(s,v)))
        oss << '(' << s << '/' << v << '/' << this->dof_number(s,v,c) << ") ";


  return oss.str();
}



void Elem::nullify_neighbors ()
{
  // Tell any of my neighbors about my death...
  // Looks strange, huh?
  for (auto n : this->side_index_range())
    {
      Elem * current_neighbor = this->neighbor_ptr(n);
      if (current_neighbor && current_neighbor != remote_elem)
        {
          // Note:  it is possible that I see the neighbor
          // (which is coarser than me)
          // but they don't see me, so avoid that case.
          if (current_neighbor->level() == this->level())
            {
              const unsigned int w_n_a_i = current_neighbor->which_neighbor_am_i(this);
              libmesh_assert_less (w_n_a_i, current_neighbor->n_neighbors());
              current_neighbor->set_neighbor(w_n_a_i, nullptr);
              this->set_neighbor(n, nullptr);
            }
        }
    }
}




unsigned int Elem::n_second_order_adjacent_vertices (const unsigned int) const
{
  // for linear elements, always return 0
  return 0;
}



unsigned short int Elem::second_order_adjacent_vertex (const unsigned int,
                                                       const unsigned int) const
{
  // for linear elements, always return 0
  return 0;
}



std::pair<unsigned short int, unsigned short int>
Elem::second_order_child_vertex (const unsigned int) const
{
  // for linear elements, always return 0
  return std::pair<unsigned short int, unsigned short int>(0,0);
}



ElemType Elem::first_order_equivalent_type (const ElemType et)
{
  switch (et)
    {
    case NODEELEM:
      return NODEELEM;
    case EDGE2:
    case EDGE3:
    case EDGE4:
      return EDGE2;
    case TRI3:
    case TRI6:
    case TRI7:
      return TRI3;
    case TRISHELL3:
      return TRISHELL3;
    case QUAD4:
    case QUAD8:
    case QUAD9:
      return QUAD4;
    case QUADSHELL4:
    case QUADSHELL8:
    case QUADSHELL9:
      return QUADSHELL4;
    case TET4:
    case TET10:
    case TET14:
      return TET4;
    case HEX8:
    case HEX27:
    case HEX20:
      return HEX8;
    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      return PRISM6;
    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      return PYRAMID5;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

    case INFEDGE2:
      return INFEDGE2;
    case INFQUAD4:
    case INFQUAD6:
      return INFQUAD4;
    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
      return INFHEX8;
    case INFPRISM6:
    case INFPRISM12:
      return INFPRISM6;

#endif

    default:
      // unknown element
      return INVALID_ELEM;
    }
}



ElemType Elem::second_order_equivalent_type (const ElemType et,
                                             const bool full_ordered)
{
  switch (et)
    {
    case NODEELEM:
      return NODEELEM;
    case EDGE2:
    case EDGE3:
      {
        // full_ordered not relevant
        return EDGE3;
      }

    case EDGE4:
      {
        // full_ordered not relevant
        return EDGE4;
      }

    case TRI3:
    case TRI6:
      {
        // full_ordered not relevant
        return TRI6;
      }

    case TRI7:
      return TRI7;

      // Currently there is no TRISHELL6, so similarly to other types
      // where this is the case, we just return the input.
    case TRISHELL3:
      return TRISHELL3;

    case QUAD4:
    case QUAD8:
      {
        if (full_ordered)
          return QUAD9;
        else
          return QUAD8;
      }

    case QUADSHELL4:
    case QUADSHELL8:
      {
        if (full_ordered)
          return QUADSHELL9;
        else
          return QUADSHELL8;
      }

    case QUAD9:
      {
        // full_ordered not relevant
        return QUAD9;
      }

    case QUADSHELL9:
      {
        // full_ordered not relevant
        return QUADSHELL9;
      }

    case TET4:
    case TET10:
      {
        // full_ordered not relevant
        return TET10;
      }

    case TET14:
      return TET14;

    case HEX8:
    case HEX20:
      {
        // see below how this correlates with INFHEX8
        if (full_ordered)
          return HEX27;
        else
          return HEX20;
      }

    case HEX27:
      {
        // full_ordered not relevant
        return HEX27;
      }

    case PRISM6:
    case PRISM15:
      {
        if (full_ordered)
          return PRISM18;
        else
          return PRISM15;
      }

    // full_ordered not relevant, already fully second order
    case PRISM18:
      return PRISM18;
    case PRISM20:
      return PRISM20;
    case PRISM21:
      return PRISM21;
    case PYRAMID18:
      return PYRAMID18;

    case PYRAMID5:
    case PYRAMID13:
      {
        if (full_ordered)
          return PYRAMID14;
        else
          return PYRAMID13;
      }

    case PYRAMID14:
      {
        // full_ordered not relevant
        return PYRAMID14;
      }



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

      // infinite elements
    case INFEDGE2:
      {
        return INFEDGE2;
      }

    case INFQUAD4:
    case INFQUAD6:
      {
        // full_ordered not relevant
        return INFQUAD6;
      }

    case INFHEX8:
    case INFHEX16:
      {
        /*
         * Note that this matches with \p Hex8:
         * For full-ordered, \p InfHex18 and \p Hex27
         * belong together, and for not full-ordered,
         * \p InfHex16 and \p Hex20 belong together.
         */
        if (full_ordered)
          return INFHEX18;
        else
          return INFHEX16;
      }

    case INFHEX18:
      {
        // full_ordered not relevant
        return INFHEX18;
      }

    case INFPRISM6:
    case INFPRISM12:
      {
        // full_ordered not relevant
        return INFPRISM12;
      }

#endif


    default:
      {
        // what did we miss?
        libmesh_error_msg("No second order equivalent element type for et =  "
                          << Utility::enum_to_string(et));
      }
    }
}



ElemType Elem::complete_order_equivalent_type (const ElemType et)
{
  switch (et)
    {
    case NODEELEM:
      return NODEELEM;
    case EDGE2:
    case EDGE3:
      return EDGE3;

    case EDGE4:
      return EDGE4;

    case TRI3:
    case TRI6:
    case TRI7:
      return TRI7;

      // Currently there is no TRISHELL6, so similarly to other types
      // where this is the case, we just return the input.
    case TRISHELL3:
      return TRISHELL3;

    case QUAD4:
    case QUAD8:
    case QUAD9:
      return QUAD9;

    case QUADSHELL4:
    case QUADSHELL8:
    case QUADSHELL9:
      return QUADSHELL9;

    case TET4:
    case TET10:
    case TET14:
      return TET14;

    case HEX8:
    case HEX20:
    case HEX27:
      return HEX27;

    // We don't strictly need the 21st node for DoFs, but some code
    // depends on it for e.g. refinement patterns
    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      return PRISM21;

    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      return PYRAMID18;



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    case INFEDGE2:
      return INFEDGE2;

    case INFQUAD4:
    case INFQUAD6:
      return INFQUAD6;

    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
      return INFHEX18;

    // Probably ought to implement INFPRISM13 and/or 14; until then
    // we'll just return the not-complete-order ElemType, in
    // accordance which what seems like bad precedent...
    case INFPRISM6:
    case INFPRISM12:
      return INFPRISM12;
#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS

    default:
      {
        // what did we miss?
        libmesh_error_msg("No complete order equivalent element type for et =  "
                          << Utility::enum_to_string(et));
      }
    }
}



Elem::side_iterator Elem::boundary_sides_begin()
{
  Predicates::BoundarySide<SideIter> bsp;
  return side_iterator(this->_first_side(), this->_last_side(), bsp);
}




Elem::side_iterator Elem::boundary_sides_end()
{
  Predicates::BoundarySide<SideIter> bsp;
  return side_iterator(this->_last_side(), this->_last_side(), bsp);
}




Real Elem::volume () const
{
  // The default implementation builds a finite element of the correct
  // order and sums up the JxW contributions.  This can be expensive,
  // so the various element types can overload this method and compute
  // the volume more efficiently.
  const FEFamily mapping_family = FEMap::map_fe_type(*this);
  const FEType fe_type(this->default_order(), mapping_family);

  std::unique_ptr<FEBase> fe (FEBase::build(this->dim(),
                                            fe_type));

  // Use a map with a negative Jacobian tolerance, in case we're asked
  // for the net volume of a tangled element
  fe->get_fe_map().set_jacobian_tolerance(std::numeric_limits<Real>::lowest());

  const std::vector<Real> & JxW = fe->get_JxW();

  // The default quadrature rule should integrate the mass matrix,
  // thus it should be plenty to compute the area
  QGauss qrule (this->dim(), fe_type.default_quadrature_order());

  fe->attach_quadrature_rule(&qrule);

  fe->reinit(this);

  Real vol=0.;
  for (auto jxw : JxW)
    vol += jxw;

  return vol;

}



BoundingBox Elem::loose_bounding_box () const
{
  Point pmin = this->point(0);
  Point pmax = pmin;

  unsigned int n_points = this->n_nodes();
  for (unsigned int p=0; p != n_points; ++p)
    for (unsigned d=0; d<LIBMESH_DIM; ++d)
      {
        const Point & pt = this->point(p);
        if (pmin(d) > pt(d))
          pmin(d) = pt(d);

        if (pmax(d) < pt(d))
          pmax(d) = pt(d);
      }

  return BoundingBox(pmin, pmax);
}

Point
Elem::side_vertex_average_normal(const unsigned int s) const
{
  unsigned int dim = this->dim();
  const std::unique_ptr<const Elem> face = this->build_side_ptr(s);
  std::unique_ptr<libMesh::FEBase> fe(
      libMesh::FEBase::build(dim, libMesh::FEType(this->default_order())));
  const std::vector<Point> & normals = fe->get_normals();
  std::vector<Point> ref_side_vertex_average_v = {face->reference_elem()->vertex_average()};
  fe->reinit(this, s, TOLERANCE, &ref_side_vertex_average_v);
  return normals[0];
}

bool Elem::is_vertex_on_parent(unsigned int c,
                               unsigned int n) const
{
#ifdef LIBMESH_ENABLE_AMR

  unsigned int my_n_vertices = this->n_vertices();
  for (unsigned int n_parent = 0; n_parent != my_n_vertices;
       ++n_parent)
    if (this->node_ptr(n_parent) == this->child_ptr(c)->node_ptr(n))
      return true;
  return false;

#else

  // No AMR?
  libmesh_ignore(c,n);
  libmesh_error_msg("ERROR: AMR disabled, how did we get here?");
  return true;

#endif
}



unsigned int Elem::opposite_side(const unsigned int /*s*/) const
{
  // If the subclass didn't rederive this, using it is an error
  libmesh_not_implemented();
}



unsigned int Elem::opposite_node(const unsigned int /*n*/,
                                 const unsigned int /*s*/) const
{
  // If the subclass didn't rederive this, using it is an error
  libmesh_not_implemented();
}


unsigned int Elem::center_node_on_side(const unsigned short libmesh_dbg_var(side)) const
{
  libmesh_assert_less (side, this->n_sides());
  return invalid_uint;
}


void Elem::swap2boundarysides(unsigned short s1,
                              unsigned short s2,
                              BoundaryInfo * boundary_info) const
{
  std::vector<boundary_id_type> ids1, ids2;
  boundary_info->boundary_ids(this, s1, ids1);
  boundary_info->boundary_ids(this, s2, ids2);
  boundary_info->remove_side(this, s1);
  boundary_info->remove_side(this, s2);
  if (!ids1.empty())
    boundary_info->add_side(this, s2, ids1);
  if (!ids2.empty())
    boundary_info->add_side(this, s1, ids2);
}


void Elem::swap2boundaryedges(unsigned short e1,
                              unsigned short e2,
                              BoundaryInfo * boundary_info) const
{
  std::vector<boundary_id_type> ids1, ids2;
  boundary_info->edge_boundary_ids(this, e1, ids1);
  boundary_info->edge_boundary_ids(this, e2, ids2);
  boundary_info->remove_edge(this, e1);
  boundary_info->remove_edge(this, e2);
  if (!ids1.empty())
    boundary_info->add_edge(this, e2, ids1);
  if (!ids2.empty())
    boundary_info->add_edge(this, e1, ids2);
}

bool
Elem::is_internal(const unsigned int i) const
{
  switch (this->dim())
  {
    case 0:
      return false;

    case 1:
      return !this->is_vertex(i);

    case 2:
      return !this->is_vertex(i) && !this->is_edge(i);

    case 3:
      return !this->is_vertex(i) && !this->is_edge(i) && !this->is_face(i);

    default:
      libmesh_error_msg("impossible element dimension " << std::to_string(this->dim()));
      return 0;
  }
}

bool
Elem::positive_edge_orientation(const unsigned int i) const
{
  libmesh_assert_less (i, this->n_edges());

  return this->point(this->local_edge_node(i, 0)) >
         this->point(this->local_edge_node(i, 1));
}

bool
Elem::positive_face_orientation(const unsigned int i) const
{
  libmesh_assert_less (i, this->n_faces());

  // Get the number of vertices N of face i. Note that for 3d elements, i.e.
  // elements for which this->n_faces() > 0, the number of vertices on any of
  // its sides (or faces) is just the number of that face's sides (or edges).
  const unsigned int N = Elem::type_to_n_sides_map[this->side_type(i)];

  const std::vector<unsigned int> nodes = this->nodes_on_side(i);

  auto cmp = [&](const unsigned int & m, const unsigned int & n) -> bool
             { return this->point(m) < this->point(n); };

  const unsigned int v = std::distance(nodes.begin(),
                         std::min_element(nodes.begin(), nodes.begin() + N, cmp));

  return cmp(nodes[(v - 1 + N) % N], nodes[(v + 1) % N]);
}

bool
Elem::relative_edge_face_order(const unsigned int e, const unsigned int s) const
{
  libmesh_assert_less (e, this->n_edges());
  libmesh_assert_less (s, this->n_faces());
  libmesh_assert (is_edge_on_side(e, s));

  const unsigned int N = Elem::type_to_n_sides_map[this->side_type(s)];

  unsigned int v;
  for (v = 0; v < N; v++)
    if (this->local_side_node(s, v) == this->local_edge_node(e, 0))
      break;

  libmesh_assert (this->local_side_node(s, (v + 1) % N) ==
                  this->local_edge_node(e, 1) ||
                  this->local_side_node(s, (v - 1 + N) % N) ==
                  this->local_edge_node(e, 1));

  return this->local_side_node(s, (v + 1) % N) ==
         this->local_edge_node(e, 1);
}

} // namespace libMesh
