// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// libmesh includes
#include "libmesh/auto_ptr.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/system.h"

// test includes
#include "test_comm.h"

using namespace libMesh;

class ParsedFEMFunctionTest : public CppUnit::TestCase
{
public:
  void setUp() {
    mesh.reset(new Mesh(*TestCommWorld));
    MeshTools::Generation::build_cube(*mesh, 1, 1, 1);
    es.reset(new EquationSystems(*mesh));
    sys = &(es->add_system<System> ("SimpleSystem"));
    sys->add_variable("x2");
    sys->add_variable("x3");
    sys->add_variable("c05");
    sys->add_variable("y4");
    sys->add_variable("xy");
    sys->add_variable("yz");
    sys->add_variable("xyz");

    es->init();

    NumericVector<Number> & sol = *sys->solution;
    Elem *elem = mesh->elem(0);

    // Set x2 = 2*x
    sol.set(elem->get_node(1)->dof_number(0,0,0), 2);
    sol.set(elem->get_node(2)->dof_number(0,0,0), 2);
    sol.set(elem->get_node(4)->dof_number(0,0,0), 2);
    sol.set(elem->get_node(5)->dof_number(0,0,0), 2);

    // Set x3 = 3*x
    sol.set(elem->get_node(1)->dof_number(0,1,0), 3);
    sol.set(elem->get_node(2)->dof_number(0,1,0), 3);
    sol.set(elem->get_node(4)->dof_number(0,1,0), 3);
    sol.set(elem->get_node(5)->dof_number(0,1,0), 3);

    // Set c05 = 0.5
    sol.set(elem->get_node(0)->dof_number(0,2,0), 0.5);
    sol.set(elem->get_node(1)->dof_number(0,2,0), 0.5);
    sol.set(elem->get_node(2)->dof_number(0,2,0), 0.5);
    sol.set(elem->get_node(3)->dof_number(0,2,0), 0.5);
    sol.set(elem->get_node(4)->dof_number(0,2,0), 0.5);
    sol.set(elem->get_node(5)->dof_number(0,2,0), 0.5);
    sol.set(elem->get_node(6)->dof_number(0,2,0), 0.5);
    sol.set(elem->get_node(7)->dof_number(0,2,0), 0.5);

    // Set y4 = 4*y
    sol.set(elem->get_node(2)->dof_number(0,3,0), 4);
    sol.set(elem->get_node(3)->dof_number(0,3,0), 4);
    sol.set(elem->get_node(6)->dof_number(0,3,0), 4);
    sol.set(elem->get_node(7)->dof_number(0,3,0), 4);

    // Set xy = x*y
    sol.set(elem->get_node(2)->dof_number(0,4,0), 1);
    sol.set(elem->get_node(6)->dof_number(0,4,0), 1);

    // Set yz = y*z
    sol.set(elem->get_node(6)->dof_number(0,5,0), 1);
    sol.set(elem->get_node(7)->dof_number(0,5,0), 1);

    // Set xyz = x*y*z
    sol.set(elem->get_node(7)->dof_number(0,6,0), 1);

    sol.close();
    sys->update();
  }

  void tearDown() {
    es.reset();
    mesh.reset();
  }

  CPPUNIT_TEST_SUITE(ParsedFEMFunctionTest);

  CPPUNIT_TEST(testValues);

  CPPUNIT_TEST_SUITE_END();


private:
  AutoPtr<Mesh> mesh;
  AutoPtr<EquationSystems> es;
  System * sys;

  void testValues()
  {
    Elem *elem = mesh->elem(0);

    ParsedFEMFunction<Number> parsed_func(*sys, "x2*y4");

    FEMContext c(*sys);

    c.pre_fe_reinit(*sys, elem);
    c.elem_fe_reinit();
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(ParsedFEMFunctionTest);
