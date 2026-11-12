#ifndef MESH_ELEM_TEST_H
#define MESH_ELEM_TEST_H

#include "../geom/elem_test.h"

#include "libmesh/mesh_serializer.h"

using namespace libMesh;

template <ElemType elem_type>
class MeshPerElemTest : public PerElemTest<elem_type>
{
protected:

  bool meshes_equal_enough(Mesh & other_mesh, bool double_precision)
  {
    // We'll need to fix up processor_id() and unique_id() values
    // before we can operator== these meshes.  But worse: our gold
    // meshes might have been numbered differently to our generated
    // meshes.  Some of our generated mesh options practically
    // *require* renumbering (e.g. after interior HEX20 nodes are
    // deleted, ExodusII still wants to see a contiguous numbering),
    // but ReplicatedMesh and DistributedMesh renumber differently.
    //
    // So, let's renumber too.

    MeshSerializer serialthis(*this->_mesh);
    MeshSerializer serialother(other_mesh);

    const dof_id_type max_elem_id = this->_mesh->max_elem_id();
    const dof_id_type max_node_id = this->_mesh->max_node_id();

    CPPUNIT_ASSERT_EQUAL(max_elem_id, other_mesh.max_elem_id());
    CPPUNIT_ASSERT_EQUAL(max_node_id, other_mesh.max_node_id());

    auto locator = other_mesh.sub_point_locator();

    for (Elem * e1 : this->_mesh->element_ptr_range())
    {
      const Elem * e2c = (*locator)(e1->vertex_average());
      CPPUNIT_ASSERT(e2c);
      Elem & e2 = other_mesh.elem_ref(e2c->id());
      e1->processor_id() = 0;
      e2.processor_id() = 0;

      const dof_id_type e1_id = e1->id();
      const dof_id_type e2_id = e2.id();
      // Do a swap if necessary, using a free temporary id
      if (e1_id != e2_id)
        {
          other_mesh.renumber_elem(e1_id, max_elem_id);
          other_mesh.renumber_elem(e2_id, e1_id);
          other_mesh.renumber_elem(max_elem_id, e2_id);
        }

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      e2.set_unique_id(e1->unique_id());
#endif
    }

    for (Node * n1 : this->_mesh->node_ptr_range())
    {
      const Elem * e1c = (*locator)(*n1);
      Node * n2 = nullptr;
      for (const Node & n : e1c->node_ref_range())
        {
#if defined(LIBMESH_DEFAULT_QUADRUPLE_PRECISION) || defined(LIBMESH_DEFAULT_TRIPLE_PRECISION)
          if (double_precision)
            {
              const Point diff = Point(*n1)-Point(n);

              // We may be testing against ExodusII input, and if
              // we're in triple or quadruple precision that means our
              // lovely higher-precision node coordinates got
              // truncated to double to be written.  We need to adjust
              // ours or they won't satisfy operator== later.

              // We're *also* testing against gold files that were
              // calculated at double precision, so just casting a
              // higher precision calculation to double won't give the
              // exact same result, we have to account for error.
              if (diff.norm() < 1e-15)
                for (auto d : make_range(LIBMESH_DIM))
                  (*n1)(d) = double(n(d));
            }
#else
          libmesh_ignore(double_precision);
#endif
          if (Point(*n1) == Point(n))
            n2 = other_mesh.node_ptr(n.id());
        }
      CPPUNIT_ASSERT(n2);
      n1->processor_id() = 0;
      n2->processor_id() = 0;

      const dof_id_type n1_id = n1->id();
      const dof_id_type n2_id = n2->id();
      // Do a swap if necessary, using a free temporary id
      if (n1_id != n2_id)
        {
          other_mesh.renumber_node(n1_id,max_node_id);
          other_mesh.renumber_node(n2_id,n1_id);
          other_mesh.renumber_node(max_node_id, n2_id);
        }

#ifdef LIBMESH_ENABLE_UNIQUE_ID
      n2->set_unique_id(n1->unique_id());
#endif
    }

    // Setting all those processor ids to 0 changes our sets of local
    // subdomains too.
    other_mesh.cache_elem_data();
    this->_mesh->cache_elem_data();

#ifdef LIBMESH_ENABLE_UNIQUE_ID
    other_mesh.set_next_unique_id(this->_mesh->parallel_max_unique_id());
    this->_mesh->set_next_unique_id(this->_mesh->parallel_max_unique_id());
#endif

    return *this->_mesh == other_mesh;
  }
};

#endif // MESH_ELEM_TEST_H
