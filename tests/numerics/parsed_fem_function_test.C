// Ignore unused parameter warnings coming from cppunit headers
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
    Elem *elem = mesh->query_elem_ptr(0);

    if (elem && elem->processor_id() == TestCommWorld->rank())
      {
        // Set x2 = 2*x
        sol.set(elem->node_ref(1).dof_number(0,0,0), 2);
        sol.set(elem->node_ref(2).dof_number(0,0,0), 2);
        sol.set(elem->node_ref(5).dof_number(0,0,0), 2);
        sol.set(elem->node_ref(6).dof_number(0,0,0), 2);

        // Set x3 = 3*x
        sol.set(elem->node_ref(1).dof_number(0,1,0), 3);
        sol.set(elem->node_ref(2).dof_number(0,1,0), 3);
        sol.set(elem->node_ref(5).dof_number(0,1,0), 3);
        sol.set(elem->node_ref(6).dof_number(0,1,0), 3);

        // Set c05 = 0.5
        sol.set(elem->node_ref(0).dof_number(0,2,0), 0.5);
        sol.set(elem->node_ref(1).dof_number(0,2,0), 0.5);
        sol.set(elem->node_ref(2).dof_number(0,2,0), 0.5);
        sol.set(elem->node_ref(3).dof_number(0,2,0), 0.5);
        sol.set(elem->node_ref(4).dof_number(0,2,0), 0.5);
        sol.set(elem->node_ref(5).dof_number(0,2,0), 0.5);
        sol.set(elem->node_ref(6).dof_number(0,2,0), 0.5);
        sol.set(elem->node_ref(7).dof_number(0,2,0), 0.5);

        // Set y4 = 4*y
        sol.set(elem->node_ref(2).dof_number(0,3,0), 4);
        sol.set(elem->node_ref(3).dof_number(0,3,0), 4);
        sol.set(elem->node_ref(6).dof_number(0,3,0), 4);
        sol.set(elem->node_ref(7).dof_number(0,3,0), 4);

        // Set xy = x*y
        sol.set(elem->node_ref(2).dof_number(0,4,0), 1);
        sol.set(elem->node_ref(6).dof_number(0,4,0), 1);

        // Set yz = y*z
        sol.set(elem->node_ref(6).dof_number(0,5,0), 1);
        sol.set(elem->node_ref(7).dof_number(0,5,0), 1);

        // Set xyz = x*y*z
        sol.set(elem->node_ref(6).dof_number(0,6,0), 1);
      }

    sol.close();
    sys->update();

    c.reset(new FEMContext(*sys));
    s.reset(new FEMContext(*sys));
    if (elem && elem->processor_id() == TestCommWorld->rank())
      {
        c->pre_fe_reinit(*sys, elem);
        c->elem_fe_reinit();
        s->pre_fe_reinit(*sys, elem);
        s->side = 3;
        s->side_fe_reinit();
      }
  }

  void tearDown() {
    c.reset();
    s.reset();
    es.reset();
    mesh.reset();
  }

  CPPUNIT_TEST_SUITE(ParsedFEMFunctionTest);

  CPPUNIT_TEST(testValues);
  CPPUNIT_TEST(testGradients);
  CPPUNIT_TEST(testHessians);
  CPPUNIT_TEST(testInlineGetter);
  CPPUNIT_TEST(testInlineSetter);
  CPPUNIT_TEST(testNormals);

  CPPUNIT_TEST_SUITE_END();


private:
  UniquePtr<UnstructuredMesh> mesh;
  UniquePtr<EquationSystems> es;
  System * sys;
  UniquePtr<FEMContext> c, s;

  void testValues()
  {
    if (c->has_elem() &&
        c->get_elem().processor_id() == TestCommWorld->rank())
      {
        ParsedFEMFunction<Number> x2(*sys, "x2");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(x2(*c,Point(0.5,0.5,0.5))), 1.0, TOLERANCE*TOLERANCE);

        ParsedFEMFunction<Number> xy8(*sys, "x2*y4");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(xy8(*c,Point(0.5,0.5,0.5))), 2.0, TOLERANCE*TOLERANCE);
      }
  }


  void testGradients()
  {
    if (c->has_elem() &&
        c->get_elem().processor_id() == TestCommWorld->rank())
      {
        ParsedFEMFunction<Number> c2(*sys, "grad_x_x2");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(c2(*c,Point(0.35,0.45,0.55))), 2.0, TOLERANCE*TOLERANCE);

        ParsedFEMFunction<Number> xz(*sys, "grad_y_xyz");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(xz(*c,Point(0.25,0.35,0.75))), 0.1875, TOLERANCE*TOLERANCE);

        ParsedFEMFunction<Number> xyz(*sys, "grad_y_xyz*grad_x_xy");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(xyz(*c,Point(0.25,0.5,0.75))), 0.09375, TOLERANCE*TOLERANCE);
      }
  }


  void testHessians()
  {
    if (c->has_elem() &&
        c->get_elem().processor_id() == TestCommWorld->rank())
      {
        ParsedFEMFunction<Number> c1(*sys, "hess_xy_xy");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(c1(*c,Point(0.35,0.45,0.55))), 1.0, TOLERANCE*TOLERANCE);

        ParsedFEMFunction<Number> x(*sys, "hess_yz_xyz");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(x(*c,Point(0.25,0.35,0.55))), 0.25, TOLERANCE*TOLERANCE);

        ParsedFEMFunction<Number> xz(*sys, "hess_yz_xyz*grad_y_yz");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(xz(*c,Point(0.25,0.4,0.75))), 0.1875, TOLERANCE*TOLERANCE);
      }
  }

  void testInlineGetter()
  {
    if (c->has_elem() &&
        c->get_elem().processor_id() == TestCommWorld->rank())
      {
        ParsedFEMFunction<Number> ax2(*sys, "a:=4.5;a*x2");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(ax2(*c,Point(0.25,0.25,0.25))), 2.25, TOLERANCE*TOLERANCE);

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(ax2.get_inline_value("a")), 4.5, TOLERANCE*TOLERANCE);

        ParsedFEMFunction<Number> cxy8
          (*sys, "a := 4 ; b := a/2+1; c:=b-a+3.5; c*x2*y4");

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(cxy8(*c,Point(0.5,0.5,0.5))), 5.0, TOLERANCE*TOLERANCE);

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(cxy8.get_inline_value("b")), 3.0, TOLERANCE*TOLERANCE);

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(cxy8.get_inline_value("c")), 2.5, TOLERANCE*TOLERANCE);
      }
  }

  void testInlineSetter()
  {
    if (c->has_elem() &&
        c->get_elem().processor_id() == TestCommWorld->rank())
      {
        ParsedFEMFunction<Number> ax2(*sys, "a:=4.5;a*x2");
        ax2.set_inline_value("a", 2.5);

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(ax2(*c,Point(0.25,0.25,0.25))), 1.25, TOLERANCE*TOLERANCE);

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(ax2.get_inline_value("a")), 2.5, TOLERANCE*TOLERANCE);

        ParsedFEMFunction<Number> cxy8
          (*sys, "a := 4 ; b := a/2+1; c:=b-a+3.5; c*x2*y4");

        cxy8.set_inline_value("a", 2);

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(cxy8(*c,Point(0.5,0.5,0.5))), 7.0, TOLERANCE*TOLERANCE);

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(cxy8.get_inline_value("b")), 2.0, TOLERANCE*TOLERANCE);

        CPPUNIT_ASSERT_DOUBLES_EQUAL
          (libmesh_real(cxy8.get_inline_value("c")), 3.5, TOLERANCE*TOLERANCE);

      }
  }

  void testNormals()
  {
    if (s->has_elem() &&
        s->get_elem().processor_id() == TestCommWorld->rank())
      {
        ParsedFEMFunction<Number> nx(*sys, "n_x");

        ParsedFEMFunction<Number> ny(*sys, "n_y");

        ParsedFEMFunction<Number> nz(*sys, "n_z");

        const std::vector<Point> & xyz = s->get_side_fe(0)->get_xyz();

        // On side 3 of a hex the normal direction is +y
        for (std::size_t qp=0; qp != xyz.size(); ++qp)
          {
            CPPUNIT_ASSERT_DOUBLES_EQUAL
              (libmesh_real(nx(*s,xyz[qp])), 0.0, TOLERANCE*TOLERANCE);
            CPPUNIT_ASSERT_DOUBLES_EQUAL
              (libmesh_real(ny(*s,xyz[qp])), 1.0, TOLERANCE*TOLERANCE);
            CPPUNIT_ASSERT_DOUBLES_EQUAL
              (libmesh_real(nz(*s,xyz[qp])), 0.0, TOLERANCE*TOLERANCE);
          }
      }
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(ParsedFEMFunctionTest);
