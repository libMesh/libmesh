#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/elem.h>
#include <libmesh/dof_map.h>

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
