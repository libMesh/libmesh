#include <libmesh/partitioner.h>
#include <libmesh/mesh.h>
#include <libmesh/auto_ptr.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/elem.h>
#include <libmesh/system.h>
#include <libmesh/dof_map.h>
#include <timpi/communicator.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class ContrivedPartitioner : public Partitioner
{
public:
  ContrivedPartitioner() = default;
  ContrivedPartitioner(const ContrivedPartitioner &) = default;
  ContrivedPartitioner(ContrivedPartitioner &&) = default;
  ContrivedPartitioner & operator=(const ContrivedPartitioner &) = default;
  ContrivedPartitioner & operator=(ContrivedPartitioner &&) = default;
  virtual ~ContrivedPartitioner() = default;

  std::unique_ptr<Partitioner> clone() const override
  {
    return libmesh_make_unique<ContrivedPartitioner>(*this);
  }

protected:
  void _do_partition(MeshBase & mesh, const unsigned int n) override
  {
    libmesh_assert(n > 1);

    unsigned int n_interior_elems = 0;
    for (Elem * const elem : as_range(mesh.elements_begin(), mesh.elements_end()))
    {
      bool internal_elem = true;
      for (const Elem * const neighbor : elem->neighbor_ptr_range())
        if (!neighbor)
        {
          internal_elem = false;
          break;
        }

      // The interior element will go on processor 1. Other elements will go on processor 0. For a
      // system/DofMap that has all nodal variables all elements on processor 0 will appear
      // evaluable. However, for a system/DofMap that has any elemental variables the interior
      // element will not appear evaluable on processor 0
      n_interior_elems += internal_elem;
      elem->processor_id() = internal_elem;
    }

    // This test is hand-crafted for one interior element
    libmesh_assert(n_interior_elems == 1);
  }
};

class MultiEvaluablePredTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(MultiEvaluablePredTest);

  CPPUNIT_TEST(test);

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void test()
  {
    const auto rank = TestCommWorld->rank();
    Mesh mesh(*TestCommWorld);
    mesh.partitioner() = libmesh_make_unique<ContrivedPartitioner>();
    MeshTools::Generation::build_square(mesh, 3, 3, 0, 3, 0, 3, QUAD4);

    EquationSystems es(mesh);
    auto & nodal_system = es.add_system<System>("nodal_system");
    nodal_system.add_variable("nodal", FIRST, LAGRANGE);
    auto & nodal_dof_map = nodal_system.get_dof_map();
    nodal_dof_map.remove_default_ghosting();
    auto & elem_system = es.add_system<System>("elem_system");
    elem_system.add_variable("elem", CONSTANT, MONOMIAL);
    auto & elem_dof_map = elem_system.get_dof_map();
    elem_dof_map.remove_default_ghosting();

    es.init();

    {
      auto n_evaluable = std::distance(mesh.evaluable_elements_begin(nodal_dof_map),
                                       mesh.evaluable_elements_end(nodal_dof_map));
      typedef decltype(n_evaluable) comp_type;
      switch (rank)
      {
        case 0:
          CPPUNIT_ASSERT_EQUAL(n_evaluable, comp_type(9));
          break;

        case 1:
          CPPUNIT_ASSERT_EQUAL(n_evaluable, comp_type(1));
          break;

        default:
          CPPUNIT_ASSERT_EQUAL(n_evaluable, comp_type(0));
          break;
      }
    }

    {
      auto n_evaluable = std::distance(mesh.evaluable_elements_begin(elem_dof_map),
                                       mesh.evaluable_elements_end(elem_dof_map));
      typedef decltype(n_evaluable) comp_type;
      switch (rank)
      {
        case 0:
          CPPUNIT_ASSERT_EQUAL(n_evaluable, comp_type(TestCommWorld->size() == 1 ? 9 : 8));
          break;

        case 1:
          CPPUNIT_ASSERT_EQUAL(n_evaluable, comp_type(1));
          break;

        default:
          CPPUNIT_ASSERT_EQUAL(n_evaluable, comp_type(0));
          break;
      }
    }

    {
      std::vector<const DofMap *> dof_maps = {&nodal_dof_map, &elem_dof_map};
      auto n_evaluable = std::distance(mesh.multi_evaluable_elements_begin(dof_maps),
                                       mesh.multi_evaluable_elements_end(dof_maps));
      typedef decltype(n_evaluable) comp_type;
      switch (rank)
      {
        case 0:
          CPPUNIT_ASSERT_EQUAL(n_evaluable, comp_type(TestCommWorld->size() == 1 ? 9 : 8));
          break;

        case 1:
          CPPUNIT_ASSERT_EQUAL(n_evaluable, comp_type(1));
          break;

        default:
          CPPUNIT_ASSERT_EQUAL(n_evaluable, comp_type(0));
          break;
      }
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(MultiEvaluablePredTest);
