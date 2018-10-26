// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/default_coupling.h>

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



Number cubic_default_coupling_test (const Point& p,
                                    const Parameters&,
                                    const std::string&,
                                    const std::string&)
{
  const Real & x = p(0);
  const Real & y = LIBMESH_DIM > 1 ? p(1) : 0;
  const Real & z = LIBMESH_DIM > 2 ? p(2) : 0;

  return x*(1-x)*(1-x) + x*x*(1-y) + x*(1-y)*(1-z) + y*(1-y)*z + z*(1-z)*(1-z);
}



class DefaultCouplingTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( DefaultCouplingTest );

#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testCouplingOnEdge3 );
#endif
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testCouplingOnQuad9 );
  CPPUNIT_TEST( testCouplingOnTri6 );
  CPPUNIT_TEST( testCouplingOnHex27 );
#endif

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testCoupling(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", THIRD, HIERARCHIC);
    sys.get_dof_map().default_algebraic_ghosting().set_n_levels(3);

    const unsigned int n_elem_per_side = 5;
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
    sys.project_solution(cubic_default_coupling_test, nullptr, es.parameters);

    std::set<dof_id_type> evaluable_elements;

    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        CPPUNIT_ASSERT(sys.get_dof_map().is_evaluable(*elem));
        evaluable_elements.insert(elem->id());

        for (unsigned int s1=0; s1 != elem->n_neighbors(); ++s1)
          {
            const Elem * n1 = elem->neighbor_ptr(s1);
            // Let's speed up this test by only checking the ghosted
            // elements which are most likely to break.
            if (!n1 ||
                n1->processor_id() == mesh.processor_id())
              continue;

            if (!evaluable_elements.count(n1->id()))
              {
                CPPUNIT_ASSERT(sys.get_dof_map().is_evaluable(*n1));
                evaluable_elements.insert(n1->id());
              }

            for (unsigned int s2=0; s2 != elem->n_neighbors(); ++s2)
              {
                const Elem * n2 = n1->neighbor_ptr(s2);
                if (!n2 ||
                    n2->processor_id() == mesh.processor_id())
                  continue;

                if (!evaluable_elements.count(n2->id()))
                  {
                    CPPUNIT_ASSERT(sys.get_dof_map().is_evaluable(*n2));
                    evaluable_elements.insert(n2->id());
                  }

                for (unsigned int s3=0; s3 != elem->n_neighbors(); ++s3)
                  {
                    const Elem * n3 = n2->neighbor_ptr(s3);
                    if (!n3 ||
                        n3->processor_id() == mesh.processor_id() ||
                        evaluable_elements.count(n3->id()))
                      continue;

                    CPPUNIT_ASSERT(sys.get_dof_map().is_evaluable(*n3));
                    evaluable_elements.insert(n3->id());

                    Point p = n3->centroid();

                    CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys.point_value(0,p,n3)),
                                                 libmesh_real(cubic_default_coupling_test(p,es.parameters,"","")),
                                                 TOLERANCE*TOLERANCE);
                  }
              }
          }
      }

    const std::size_t n_evaluable =
      std::distance(mesh.evaluable_elements_begin(sys.get_dof_map()),
                    mesh.evaluable_elements_end(sys.get_dof_map()));

    CPPUNIT_ASSERT_EQUAL(evaluable_elements.size(), n_evaluable);
  }



  void testCouplingOnEdge3() { testCoupling(EDGE3); }
  void testCouplingOnQuad9() { testCoupling(QUAD9); }
  void testCouplingOnTri6()  { testCoupling(TRI6); }
  void testCouplingOnHex27() { testCoupling(HEX27); }

};

CPPUNIT_TEST_SUITE_REGISTRATION( DefaultCouplingTest );
