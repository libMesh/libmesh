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
#include <libmesh/point_neighbor_coupling.h>

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



Number cubic_point_neighbor_coupling_test (const Point& p,
                                           const Parameters&,
                                           const std::string&,
                                           const std::string&)
{
  const Real & x = p(0);
  const Real & y = p(1);
  const Real & z = p(2);

  return x*(1-x)*(1-x) + x*x*(1-y) + x*(1-y)*(1-z) + y*(1-y)*z + z*(1-z)*(1-z);
}



class PointNeighborCouplingTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( PointNeighborCouplingTest );

  CPPUNIT_TEST( testCouplingOnEdge3 );
  CPPUNIT_TEST( testCouplingOnQuad9 );
  CPPUNIT_TEST( testCouplingOnTri6 );
  CPPUNIT_TEST( testCouplingOnHex27 );

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

    // Remove the default DoF ghosting functors
    sys.get_dof_map().remove_coupling_functor
      (sys.get_dof_map().default_coupling());
    sys.get_dof_map().remove_algebraic_ghosting_functor
      (sys.get_dof_map().default_algebraic_ghosting());

    // Create a replacement functor
    PointNeighborCoupling point_neighbor_coupling;

    // This just re-sets the default; real users may want a real
    // coupling matrix instead.
    point_neighbor_coupling.set_dof_coupling(NULL);

    point_neighbor_coupling.set_n_levels(3);

    sys.get_dof_map().add_algebraic_ghosting_functor
      (point_neighbor_coupling);

    const unsigned n_elem_per_side = 5;
    const UniquePtr<Elem> test_elem = Elem::build(elem_type);
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
    sys.project_solution(cubic_point_neighbor_coupling_test, NULL, es.parameters);

    for (MeshBase::const_element_iterator
           elem_it  = mesh.active_local_elements_begin(),
           elem_end = mesh.active_local_elements_end();
         elem_it != elem_end; ++elem_it)
      {
        const Elem * elem = *elem_it;
        for (unsigned int s1=0; s1 != elem->n_neighbors(); ++s1)
          {
            const Elem * n1 = elem->neighbor_ptr(s1);
            if (!n1)
              continue;

            libmesh_assert(sys.get_dof_map().is_evaluable(*n1, 0));

            // Let's speed up this test by only checking the ghosted
            // elements which are most likely to break.
            if (n1->processor_id() == mesh.processor_id())
              continue;

            for (unsigned int s2=0; s2 != elem->n_neighbors(); ++s2)
              {
                const Elem * n2 = elem->neighbor_ptr(s2);
                if (!n2 ||
                    n2->processor_id() == mesh.processor_id())
                  continue;

                libmesh_assert(sys.get_dof_map().is_evaluable(*n2, 0));

                for (unsigned int s3=0; s3 != elem->n_neighbors(); ++s3)
                  {
                    const Elem * n3 = elem->neighbor_ptr(s3);
                    if (!n3 ||
                        n3->processor_id() == mesh.processor_id())
                      continue;

                    libmesh_assert(sys.get_dof_map().is_evaluable(*n3, 0));

                    Point p = n3->centroid();

                    CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys.point_value(0,p,n3)),
                                                 libmesh_real(cubic_point_neighbor_coupling_test(p,es.parameters,"","")),
                                                 TOLERANCE*TOLERANCE);
                  }
              }
          }
      }
  }



  void testCouplingOnEdge3() { testCoupling(EDGE3); }
  void testCouplingOnQuad9() { testCoupling(QUAD9); }
  void testCouplingOnTri6()  { testCoupling(TRI6); }
  void testCouplingOnHex27() { testCoupling(HEX27); }

};

CPPUNIT_TEST_SUITE_REGISTRATION( PointNeighborCouplingTest );
