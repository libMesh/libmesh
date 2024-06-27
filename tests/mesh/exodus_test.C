#include "../geom/elem_test.h"

#ifdef LIBMESH_HAVE_EXODUS_API

#include "libmesh/enum_to_string.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_serializer.h"

using namespace libMesh;

template <ElemType elem_type>
class ExodusTest : public PerElemTest<elem_type> {

private:

  bool meshes_equal_enough(Mesh & other_mesh)
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
        if (Point(*n1) == Point(n))
          n2 = other_mesh.node_ptr(n.id());
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

#ifdef LIBMESH_ENABLE_UNIQUE_ID
    other_mesh.set_next_unique_id(this->_mesh->parallel_max_unique_id());
    this->_mesh->set_next_unique_id(this->_mesh->parallel_max_unique_id());
#endif

    return *this->_mesh == other_mesh;
  }

public:

  void test_read_gold()
  {
    LOG_UNIT_TEST;

    Mesh input_mesh(*TestCommWorld);

    ExodusII_IO exii(input_mesh);
    if (input_mesh.processor_id() == 0)
      exii.read("meshes/exodus_elements/read_exodus_" +
                Utility::enum_to_string(elem_type) + ".e");

    MeshCommunication().broadcast(input_mesh);
    input_mesh.prepare_for_use();

    CPPUNIT_ASSERT(this->meshes_equal_enough(input_mesh));
  }

  void test_write()
  {
    LOG_UNIT_TEST;

    // This is a *buffered* write; we use scope to make sure the
    // ExodusII_IO object gets destructed (and thus is guaranteed to
    // finish writing and close the file) before we try to read what
    // was written.
    {
      ExodusII_IO exii(*this->_mesh);
      exii.write("write_exodus_" +
                 Utility::enum_to_string(elem_type) + ".e");
    }

    Mesh input_mesh(*TestCommWorld);
    ExodusII_IO exii_input(input_mesh);
    if (input_mesh.processor_id() == 0)
      exii_input.read("write_exodus_" +
                      Utility::enum_to_string(elem_type) + ".e");

    MeshCommunication().broadcast(input_mesh);
    input_mesh.prepare_for_use();

    CPPUNIT_ASSERT(this->meshes_equal_enough(input_mesh));
  }
};

#define EXODUSTEST                              \
  CPPUNIT_TEST( test_read_gold );               \
  CPPUNIT_TEST( test_write );

#define INSTANTIATE_EXODUSTEST(elemtype)                        \
  class ExodusTest_##elemtype : public ExodusTest<elemtype> {   \
  public:                                                       \
  ExodusTest_##elemtype() :                                     \
    ExodusTest<elemtype>() {                                    \
    if (unitlog->summarized_logs_enabled())                     \
      this->libmesh_suite_name = "ExodusTest";                  \
    else                                                        \
      this->libmesh_suite_name = "ExodusTest_" #elemtype;       \
  }                                                             \
  CPPUNIT_TEST_SUITE( ExodusTest_##elemtype );                  \
  EXODUSTEST;                                                   \
  CPPUNIT_TEST_SUITE_END();                                     \
  };                                                            \
                                                                \
  CPPUNIT_TEST_SUITE_REGISTRATION( ExodusTest_##elemtype )

INSTANTIATE_EXODUSTEST(EDGE2);
INSTANTIATE_EXODUSTEST(EDGE3);
INSTANTIATE_EXODUSTEST(EDGE4);

#if LIBMESH_DIM > 1
INSTANTIATE_EXODUSTEST(TRI3);
INSTANTIATE_EXODUSTEST(TRISHELL3);
INSTANTIATE_EXODUSTEST(TRI6);
INSTANTIATE_EXODUSTEST(TRI7);

INSTANTIATE_EXODUSTEST(QUAD4);
INSTANTIATE_EXODUSTEST(QUADSHELL4);
INSTANTIATE_EXODUSTEST(QUAD8);
INSTANTIATE_EXODUSTEST(QUADSHELL8);
INSTANTIATE_EXODUSTEST(QUAD9);
INSTANTIATE_EXODUSTEST(QUADSHELL9);
#endif // LIBMESH_DIM > 1

#if LIBMESH_DIM > 2
INSTANTIATE_EXODUSTEST(TET4);
INSTANTIATE_EXODUSTEST(TET10);
INSTANTIATE_EXODUSTEST(TET14);

INSTANTIATE_EXODUSTEST(HEX8);
INSTANTIATE_EXODUSTEST(HEX20);
INSTANTIATE_EXODUSTEST(HEX27);

INSTANTIATE_EXODUSTEST(PRISM6);
INSTANTIATE_EXODUSTEST(PRISM15);
INSTANTIATE_EXODUSTEST(PRISM18);
INSTANTIATE_EXODUSTEST(PRISM20);
INSTANTIATE_EXODUSTEST(PRISM21);

// These tests use PointLocator, which uses contains_point(), which
// uses inverse_map(), which doesn't play nicely on Pyramids unless we
// have exceptions support
#ifdef LIBMESH_ENABLE_EXCEPTIONS
INSTANTIATE_EXODUSTEST(PYRAMID5);
INSTANTIATE_EXODUSTEST(PYRAMID13);
INSTANTIATE_EXODUSTEST(PYRAMID14);
INSTANTIATE_EXODUSTEST(PYRAMID18);
#endif
#endif // LIBMESH_DIM > 2

#endif // LIBMESH_HAVE_EXODUS_API
