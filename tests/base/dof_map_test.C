#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/elem.h>
#include <libmesh/dof_map.h>
#include <libmesh/partitioner.h>
#include <libmesh/replicated_mesh.h>

#include <timpi/parallel_implementation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <regex>
#include <string>

using namespace libMesh;

#ifdef LIBMESH_ENABLE_CONSTRAINTS
// This class is used by testConstraintLoopDetection
class MyConstraint : public System::Constraint
{
private:

  System & _sys;

public:

  MyConstraint( System & sys ) : Constraint(), _sys(sys) {}

  virtual ~MyConstraint() {}

  void constrain()
  {
    {
      const dof_id_type constrained_dof_index = 0;
      DofConstraintRow constraint_row;
      constraint_row[1] = 1.0;
      _sys.get_dof_map().add_constraint_row( constrained_dof_index, constraint_row, 0., true);
    }
    {
      const dof_id_type constrained_dof_index = 1;
      DofConstraintRow constraint_row;
      constraint_row[0] = 1.0;
      _sys.get_dof_map().add_constraint_row( constrained_dof_index, constraint_row, 0., true);
    }
  }
};

class GhostDofConstraint : public System::Constraint
{
private:

  System & _sys;

public:

  GhostDofConstraint(System & sys) : Constraint(), _sys(sys) {}

  void constrain() override
  {
    const unsigned int sys_num = _sys.number();

    // Constrain the interface DOF on node 3 by the distant DOF on node 0.
    // On the processor owning element 3, both DOFs are nonlocal, but only
    // node 3 is directly needed by the local element.
    const dof_id_type constrained_dof =
      _sys.get_mesh().node_ref(3).dof_number(sys_num, 0, 0);
    const dof_id_type dependency_dof =
      _sys.get_mesh().node_ref(0).dof_number(sys_num, 0, 0);

    // This represents constrained_dof = dependency_dof.
    DofConstraintRow constraint_row;
    constraint_row[dependency_dof] = 1.;
    _sys.get_dof_map().add_constraint_row(constrained_dof, constraint_row);
  }
};

class GhostDofConstraintPartitioner : public Partitioner
{
public:

  std::unique_ptr<Partitioner> clone() const override
  {
    return std::make_unique<GhostDofConstraintPartitioner>(*this);
  }

protected:

  void _do_partition(MeshBase & mesh, const unsigned int n) override
  {
    // Keep the first three elements on processor 0 and move the final
    // element to another processor, creating a processor interface at node 3.
    for (auto & elem : mesh.active_element_ptr_range())
      elem->processor_id() = 0;

    if (n > 1)
      mesh.elem_ref(3).processor_id() = n - 1;
  }
};
#endif


class DofMapTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( DofMapTest );

  CPPUNIT_TEST( testDofOwnerOnEdge3 );
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testDofOwnerOnQuad9 );
  CPPUNIT_TEST( testDofOwnerOnTri6 );
#endif
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testDofOwnerOnHex27 );
#endif

#if defined(LIBMESH_ENABLE_EXCEPTIONS)
  CPPUNIT_TEST( testBadElemFECombo );
#endif

#if defined(LIBMESH_ENABLE_CONSTRAINTS) && defined(LIBMESH_ENABLE_EXCEPTIONS) && LIBMESH_DIM > 1
  CPPUNIT_TEST( testConstraintLoopDetection );
#endif

#if defined(LIBMESH_ENABLE_CONSTRAINTS)
  CPPUNIT_TEST( testGhostConstraintSendList );
#endif

  CPPUNIT_TEST( testArrayDofIndices );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testDofOwner(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", THIRD, HIERARCHIC);

    const unsigned int n_elem_per_side = 3;
    const std::unique_ptr<Elem> test_elem = Elem::build(elem_type);
    const unsigned int ymax = test_elem->dim() > 1;
    const unsigned int zmax = test_elem->dim() > 2;
    const unsigned int ny = ymax * n_elem_per_side;
    const unsigned int nz = zmax * n_elem_per_side;

    MeshTools::Generation::build_cube (mesh,
                                       n_elem_per_side,
                                       ny,
                                       nz,
                                       0., 1.,
                                       0., ymax,
                                       0., zmax,
                                       elem_type);

    es.init();

    DofMap & dof_map = sys.get_dof_map();
    for (dof_id_type id = 0; id != dof_map.n_dofs(); ++id)
      {
        const processor_id_type pid = dof_map.dof_owner(id);
        CPPUNIT_ASSERT(dof_map.first_dof(pid) <= id);
        CPPUNIT_ASSERT(id < dof_map.end_dof(pid));
      }
  }



  void testDofOwnerOnEdge3() { LOG_UNIT_TEST; testDofOwner(EDGE3); }
  void testDofOwnerOnQuad9() { LOG_UNIT_TEST; testDofOwner(QUAD9); }
  void testDofOwnerOnTri6()  { LOG_UNIT_TEST; testDofOwner(TRI6); }
  void testDofOwnerOnHex27() { LOG_UNIT_TEST; testDofOwner(HEX27); }

#if defined(LIBMESH_ENABLE_EXCEPTIONS)
  void testBadElemFECombo()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", SECOND);

    MeshTools::Generation::build_square (mesh,4,4,-1., 1.,-1., 1., QUAD4);

    // We need at least one element per processor to make sure
    // everyone throws and we don't get out of sync before going on to
    // future tests.
    dof_id_type min_local_elem = mesh.n_local_elem();
    mesh.comm().min(min_local_elem);

    if (!min_local_elem)
      return;

    // We can't just CPPUNIT_ASSERT_THROW, because we want to make
    // sure we were thrown from the right place with the right error
    // message!
    bool threw_desired_exception = false;
    try {
      es.init();
    }
    catch (libMesh::LogicError & e) {
      std::regex msg_regex("only supports FEInterface::max_order");
      CPPUNIT_ASSERT(std::regex_search(e.what(), msg_regex));
      threw_desired_exception = true;
    }
    catch (...) {
      CPPUNIT_ASSERT_MESSAGE("Unexpected exception type thrown", false);
    }

    // If we have more than 4*4 processors, or a poor partitioner, we
    // might not get an exception on every processor
    mesh.comm().max(threw_desired_exception);

    CPPUNIT_ASSERT(threw_desired_exception);
  }
#endif

#if defined(LIBMESH_ENABLE_CONSTRAINTS) && defined(LIBMESH_ENABLE_EXCEPTIONS)
  void testConstraintLoopDetection()
  {
    LOG_UNIT_TEST;
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST);

    MyConstraint my_constraint(sys);
    sys.attach_constraint_object(my_constraint);

    MeshTools::Generation::build_square (mesh,4,4,-1., 1.,-1., 1., QUAD4);

    // Tell the dof_map to check for constraint loops
    DofMap & dof_map = sys.get_dof_map();
    dof_map.set_error_on_constraint_loop(true);

    CPPUNIT_ASSERT_THROW_MESSAGE("Constraint loop not detected", es.init(), libMesh::LogicError);
  }
#endif

#if defined(LIBMESH_ENABLE_CONSTRAINTS)
  void testGhostConstraintSendList()
  {
    LOG_UNIT_TEST;

    if (TestCommWorld->size() == 1)
      return;

    ReplicatedMesh mesh(*TestCommWorld);
    mesh.partitioner() = std::make_unique<GhostDofConstraintPartitioner>();

    // Build this four-element line:
    //
    //   node:  0 ----- 1 ----- 2 ----- 3 ----- 4
    //   elem:     0       1       2       3
    //   owner:    0       0       0     remote
    //
    // Node 3 remains owned by processor 0, so its DOF is ghosted on the
    // processor owning element 3. Node 0 is too distant to be ghosted by
    // the ordinary element-neighbor send-list construction.
    MeshTools::Generation::build_line(mesh, 4, 0., 4., EDGE2);
    const processor_id_type ghost_dof_processor = mesh.elem_ref(3).processor_id();
    mesh.node_ref(3).processor_id() = 0;

    EquationSystems es(mesh);
    // The default route should not ghost dependencies of constrained ghost DOFs.
    System & default_sys = es.add_system<System>("DefaultSystem");
    default_sys.add_variable("u", FIRST);
    // The enabled route should add those dependencies to the send list.
    System & ghost_sys = es.add_system<System>("GhostSystem");
    ghost_sys.add_variable("u", FIRST);

    GhostDofConstraint default_constraint(default_sys);
    default_sys.attach_constraint_object(default_constraint);
    GhostDofConstraint ghost_constraint(ghost_sys);
    ghost_sys.attach_constraint_object(ghost_constraint);

    // Enable the new behavior only for GhostSystem so DefaultSystem provides
    // a control route using the original send-list behavior.
    DofMap & ghost_dof_map = ghost_sys.get_dof_map();
    ghost_dof_map.ghost_constraints_needed();
    es.init();

    // Only the processor owning element 3 has the constrained DOF as a ghost
    // supported by one of its local elements.
    if (mesh.processor_id() == ghost_dof_processor)
      {
        const Elem & local_elem = mesh.elem_ref(3);
        CPPUNIT_ASSERT_EQUAL(mesh.processor_id(), local_elem.processor_id());

        auto verify_send_list = [&local_elem](System & sys, const bool dependency_is_ghosted)
        {
          DofMap & dof_map = sys.get_dof_map();
          const dof_id_type ghost_constrained_dof =
            sys.get_mesh().node_ref(3).dof_number(sys.number(), 0, 0);
          const dof_id_type dependency_dof =
            sys.get_mesh().node_ref(0).dof_number(sys.number(), 0, 0);

          CPPUNIT_ASSERT(!dof_map.local_index(ghost_constrained_dof));
          CPPUNIT_ASSERT(!dof_map.local_index(dependency_dof));

          // Verify that the nonlocal constrained DOF is supported by this
          // processor's local element, making it a constrained ghost DOF.
          std::vector<dof_id_type> local_dofs;
          dof_map.dof_indices(&local_elem, local_dofs);
          CPPUNIT_ASSERT(std::find(local_dofs.begin(), local_dofs.end(), ghost_constrained_dof) !=
                         local_dofs.end());

          const auto & send_list = dof_map.get_send_list();
          // The constrained DOF is always ghosted because the local element
          // uses it; its distant dependency is ghosted only when requested.
          CPPUNIT_ASSERT(std::binary_search(send_list.begin(), send_list.end(),
                                            ghost_constrained_dof));
          // This assertion is false for the control route and true for the
          // route with ghost constraint dependency expansion enabled.
          CPPUNIT_ASSERT_EQUAL(dependency_is_ghosted,
                               std::binary_search(send_list.begin(), send_list.end(),
                                                  dependency_dof));
        };

        verify_send_list(default_sys, false);
        verify_send_list(ghost_sys, true);
      }
  }
#endif

  void testArrayDofIndicesWithType(const FEType & fe_type)
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    auto & sys = es.add_system<System>("SimpleSystem");
    const auto last_var = sys.add_variable_array({"u0", "u1"}, fe_type);
    MeshTools::Generation::build_square (mesh,1,1,-1., 1.,-1., 1., QUAD9);
    es.init();
    auto & dof_map = sys.get_dof_map();
    const auto * const elem = mesh.query_elem_ptr(0);
    if (elem)
    {
      std::vector<dof_id_type> array_dofs, work_dofs, dofs;
      dof_map.array_dof_indices(elem, array_dofs, last_var);
      // Make sure there are no duplicates
      std::sort(array_dofs.begin(), array_dofs.end());
      auto it = std::unique(array_dofs.begin(), array_dofs.end());
      array_dofs.erase(it, array_dofs.end());
      CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(dof_map.n_dofs()), array_dofs.size());
      dof_map.dof_indices(elem, work_dofs, last_var);
      dofs = work_dofs;
      dof_map.dof_indices(elem, work_dofs, last_var - 1);
      dofs.insert(dofs.end(), work_dofs.begin(), work_dofs.end());
      std::sort(dofs.begin(), dofs.end());
      CPPUNIT_ASSERT(array_dofs == dofs);
    }
  }

  void testArrayDofIndices()
  {
    LOG_UNIT_TEST;
    for (int i = 1; i < 3; ++i)
    {
      testArrayDofIndicesWithType({i, LAGRANGE});
      testArrayDofIndicesWithType({i, HIERARCHIC});
      testArrayDofIndicesWithType({i, MONOMIAL});
      testArrayDofIndicesWithType({i, L2_LAGRANGE});
    }
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( DofMapTest );
