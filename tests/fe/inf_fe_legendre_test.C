// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// libmesh includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/inf_elem_builder.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"
#include "libmesh/inf_fe.h"

// unit test includes
#include "test_comm.h"

// C++ includes
#include <tuple>

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

class InfFELegendreTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( InfFELegendreTest );
  CPPUNIT_TEST( testDifferentOrders );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}

  struct TabulatedVal
  {
    TabulatedVal(unsigned int i_in, unsigned int qp_in, Real val_in)
      : i(i_in), qp(qp_in), val(val_in) {}
    unsigned int i;
    unsigned int qp;
    Real val;
  };

  struct TabulatedGrad
  {
    TabulatedGrad(unsigned int i_in, unsigned int qp_in, Point grad_in)
      : i(i_in), qp(qp_in), grad(grad_in) {}
    unsigned int i;
    unsigned int qp;
    Point grad;
  };

  typedef std::vector<TabulatedVal> TabulatedVals;
  typedef std::vector<TabulatedGrad> TabulatedGrads;

  void testDifferentOrders ()
  {
    // std::cout << "Called testDifferentOrders." << std::endl;

    // Arbitrarily selected (i, qp, phi) sample values for FIRST-order case.
    TabulatedVals first_vals =
      {
        {0,9,0.0550016}, {1,5,0.41674},   {2,8,0.0147376}, {3,7,0.111665},
        {4,6,0.0737011}, {5,1,0.0803773}, {6,11,0.275056}, {7,15,0.021537}
      };

    // Arbitrarily selected (i, qp, dphi) sample values for FIRST-order case.
    TabulatedGrads first_grads =
      {
        {0,9,Point(-0.0429458, -0.0115073, 0.)},
        {1,5,Point(0.177013, -0.177013, -0.483609)},
        {2,8,Point(0.0115073, 0.0115073, 0.00842393)},
        {3,7,Point(-0.177013, 0.0474305, 0.)},
        {4,6,Point(-0.031305, -0.116832,  0.10025)},
        {5,1,Point(0.0474191, -0.0474191, 0.872917)},
        {6,11,Point(0.0575466, 0.0575466, -0.11251)},
        {7,15,Point(-0.00353805, 0.000948017, 0.000111572)}
      };

    testSingleOrder(FIRST, first_vals, first_grads);

    // Arbitrarily selected (i, qp, phi) sample values for SECOND-order case.
    TabulatedVals second_vals =
      {
        {0,3,0.0425633}, {1,6,0.0343526}, {2,2,0.158848}, {3,7,0.128206},
        {4,12,0.220829}, {5,5,-0.136549}, {6,0,0.0149032}, {7,10,-0.0334936},
        {8,16,0.00399329}, {9,18,-0.00209733}, {10,2,0.0556194}, {11,17,-0.000561977}
      };

    // Arbitrarily selected (i, qp, dphi) sample values for SECOND-order case.
    TabulatedGrads second_grads =
      {
        {0,3,Point(-0.0959817, -0.0959817, 0.0702635)},
        {1,6,Point(0.0625228, -0.0625228, 0.0457699)},
        {2,2,Point(0.358209, 0.0959817, -2.77556e-17)},
        {3,7,Point(-0.233338, 0.0625228, -6.93889e-17)},
        {4,12,Point(-0.0323071, -0.0323071, -0.0729771)},
        {5,5,Point(0.248523, 0.0665915, -0.245097)},
        {6,0,Point(0.0336072, -0.00900502, 0.288589)},
        {7,10,Point(-0.0396234, 0.0396234, -0.0290064)},
        {8,16,Point(0.000443217, 0.000443217, 0.000333678)},
        {9,18,Point(-0.000232783, -6.23741e-05, 9.35433e-05)},
        {10,2,Point(-0.0336072, 0.0336072, 0.985214)},
        {11,17,Point(6.23741e-05, -6.23741e-05, -2.05961e-05)},
      };

    testSingleOrder(SECOND, second_vals, second_grads);

    // Arbitrarily selected (i, qp, phi) sample values for THIRD-order case.
    TabulatedVals third_vals =
      {
        {12,9,0.136657},
        {13,6,0.175034},
        {14,20,-0.0011016},
        {15,15,0.0428935}
      };

    // Arbitrarily selected (i, qp, dphi) sample values for THIRD-order case.
    TabulatedGrads third_grads =
      {
        {12,9,Point(0.0536552, 0.200244, -0.0849541)},
        {13,6,Point(-0.0921697, 0.0921697, 0.461056)},
        {14,20,Point(2.35811e-05, -8.80059e-05, 3.58959e-05)},
        {15,15,Point(-0.0386352, 0.0103523, -0.0197323)},
      };

    testSingleOrder(THIRD, third_vals, third_grads);

    // Arbitrarily selected (i, qp, phi) sample values for FOURTH-order case.
    TabulatedVals fourth_vals =
      {
        {16,3,0.00826618},
        {17,10,-0.547813},
        {18,14,0.311004},
        {19,27,-0.00192085}
      };

    // Arbitrarily selected (i, qp, dphi) sample values for FOURTH-order case.
    TabulatedGrads fourth_grads =
      {
        {16,3,Point(-0.0190603, 0.0051072, 0.308529)},
        {17,10,Point(0.244125, -0.244125, 0.140907)},
        {18,14,Point(-0.0985844, 0.0985844, -0.502591)},
        {19,27,Point(0.000115647, -3.09874e-05, 4.30775e-05)}
      };

    testSingleOrder(FOURTH, fourth_vals, fourth_grads);
  }

  void testSingleOrder (Order radial_order,
                        const TabulatedVals & val_table,
                        const TabulatedGrads & grad_table)
  {
    // Avoid warnings when not compiled with infinite element support.
    libmesh_ignore(radial_order, val_table, grad_table);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    // std::cout << "Called testSingleOrder with radial_order = "
    //           << radial_order
    //           << std::endl;

    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_cube
      (mesh,
       /*nx=*/1, /*ny=*/1, /*nz=*/1,
       /*xmin=*/-1., /*xmax=*/1.,
       /*ymin=*/-1., /*ymax=*/1.,
       /*zmin=*/-1., /*zmax=*/1.,
       HEX8);

    // Add infinite elements to the mesh. The optional verbose flag only
    // prints extra information in non-optimized mode.
    InfElemBuilder builder(mesh);
    builder.build_inf_elem(/*verbose=true*/);

    // Get pointer to the last infinite Elem. I originally intended to
    // get a pointer to the *first* infinite Elem but then I computed
    // all the reference values on the last Elem, so that's the one we
    // are using now...
    const Elem * infinite_elem = mesh.elem_ptr(mesh.n_elem() - 1);
    if (!infinite_elem || !infinite_elem->infinite())
      libmesh_error_msg("Error setting Elem pointer.");

    // We will construct FEs, etc. of the same dimension as the mesh elements.
    auto dim = mesh.mesh_dimension();

    FEType fe_type(/*Order*/FIRST,
                   /*FEFamily*/LAGRANGE,
                   radial_order,
                   /*radial_family*/LEGENDRE,
                   /*inf_map*/CARTESIAN);

    // Construct FE, quadrature rule, etc.
    std::unique_ptr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    inf_fe->attach_quadrature_rule(&qrule);
    const std::vector<std::vector<Real>> & phi = inf_fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = inf_fe->get_dphi();

    // Reinit on infinite elem
    inf_fe->reinit(infinite_elem);

    // Check phi values against reference values. The number of
    // quadrature points and shape functions seem to both increase by
    // 4 as the radial_order is increased.
    //
    // radial_order | n_qp | n_sf
    // --------------------------
    // FIRST        | 16   | 8
    // SECOND       | 20   | 12
    // THIRD        | 24   | 16
    // FOURTH       | 28   | 20
    // FIFTH        | 32   | 24

    // Debugging print values so we can tabulate them
    // auto n_qp = inf_fe->n_quadrature_points();
    // auto n_sf = inf_fe->n_shape_functions();
    // for (unsigned int i=0; i<n_sf; ++i)
    //   for (unsigned qp=0; qp<n_qp; ++qp)
    //     {
    //       libMesh::out << "phi[" << i << "][" << qp << "]=" << phi[i][qp] << std::endl;
    //       libMesh::out << "dphi[" << i << "][" << qp << "]=" << dphi[i][qp] << std::endl;
    //     }

    // Test whether the computed values match the tabulated values to
    // the specified accuracy.
    for (const auto & t : val_table)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(phi[t.i][t.qp], t.val, 1.e-4);
    for (const auto & t : grad_table)
      CPPUNIT_ASSERT_DOUBLES_EQUAL((dphi[t.i][t.qp] - t.grad).norm_sq(), 0., 1.e-4);

#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( InfFELegendreTest );
