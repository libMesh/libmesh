// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/elem.h>
#include <libmesh/dof_map.h>

#include "test_comm.h"

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

using namespace libMesh;

// This class is used by testCyclicConstraintDetection
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


class DofMapTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( DofMapTest );

  CPPUNIT_TEST( testDofOwnerOnEdge3 );
  CPPUNIT_TEST( testDofOwnerOnQuad9 );
  CPPUNIT_TEST( testDofOwnerOnTri6 );
  CPPUNIT_TEST( testDofOwnerOnHex27 );

#if defined(LIBMESH_ENABLE_CONSTRAINTS) && defined(LIBMESH_ENABLE_EXCEPTIONS)
  CPPUNIT_TEST( testCyclicConstraintDetection );
#endif

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

    const unsigned n_elem_per_side = 3;
    const std::unique_ptr<Elem> test_elem = Elem::build(elem_type);
    const Real ymax = test_elem->dim() > 1;
    const Real zmax = test_elem->dim() > 2;
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



  void testDofOwnerOnEdge3() { testDofOwner(EDGE3); }
  void testDofOwnerOnQuad9() { testDofOwner(QUAD9); }
  void testDofOwnerOnTri6()  { testDofOwner(TRI6); }
  void testDofOwnerOnHex27() { testDofOwner(HEX27); }

#if defined(LIBMESH_ENABLE_CONSTRAINTS) && defined(LIBMESH_ENABLE_EXCEPTIONS)
  void testCyclicConstraintDetection()
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST);

    MyConstraint my_constraint(sys);
    sys.attach_constraint_object(my_constraint);

    MeshTools::Generation::build_square (mesh,4,4,-1., 1.,-1., 1., QUAD4);

    // Tell the dof_map to check for cyclic constraints
    DofMap & dof_map = sys.get_dof_map();
    dof_map.set_error_on_cyclic_constraint(true);

    CPPUNIT_ASSERT_THROW_MESSAGE("Cyclic constraint not detected", es.init(), libMesh::LogicError);
  }
#endif

};

CPPUNIT_TEST_SUITE_REGISTRATION( DofMapTest );
