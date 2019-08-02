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
#include <map>
#include <vector>

#include "libmesh_cppunit.h"

using namespace libMesh;

/**
 * This class is for unit testing the "radial" basis function
 * aspects of the InfFE. There are several possible choices
 * for these basis functions, including:
 * - JACOBI_20_00
 * - JACOBI_30_00
 * - LEGENDRE
 * - LAGRANGE
 * Only the first three are currently tested by this class.
 */
class InfFERadialTest : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( InfFERadialTest );
  CPPUNIT_TEST( testDifferentOrders );
  CPPUNIT_TEST_SUITE_END();

public:
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
  typedef std::pair<Order, FEFamily> KeyType;
  std::map<KeyType, TabulatedVals> tabulated_vals;
  std::map<KeyType, TabulatedGrads> tabulated_grads;

  void testDifferentOrders ()
  {
    testSingleOrder(FIRST, LEGENDRE);
    testSingleOrder(SECOND, LEGENDRE);
    testSingleOrder(THIRD, LEGENDRE);
    testSingleOrder(FOURTH, LEGENDRE);

    testSingleOrder(FIRST, JACOBI_20_00);
    testSingleOrder(SECOND, JACOBI_20_00);
    testSingleOrder(THIRD, JACOBI_20_00);
    testSingleOrder(FOURTH, JACOBI_20_00);

    testSingleOrder(FIRST, JACOBI_30_00);
    testSingleOrder(SECOND, JACOBI_30_00);
    testSingleOrder(THIRD, JACOBI_30_00);
    testSingleOrder(FOURTH, JACOBI_30_00);
  }

  void testSingleOrder (Order radial_order,
                        FEFamily radial_family)
  {
    // Avoid warnings when not compiled with infinite element support.
    libmesh_ignore(radial_order, radial_family);

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
                   radial_family,
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
    //       libMesh::out << "{" << i << "," << qp << "," << phi[i][qp] << "}," << std::endl;
    //     }
    // for (unsigned int i=0; i<n_sf; ++i)
    //   for (unsigned qp=0; qp<n_qp; ++qp)
    //     {
    //       libMesh::out << "{" << i << "," << qp << ",Point("
    //                    << dphi[i][qp](0)
    //                    << ","
    //                    << dphi[i][qp](1)
    //                    << ","
    //                    << dphi[i][qp](2)
    //                    << ")},"
    //                    << std::endl;
    //     }

    // Get reference vals for this Order/Family combination.
    const TabulatedVals & val_table =
      tabulated_vals[std::make_pair(radial_order, radial_family)];
    const TabulatedGrads & grad_table =
      tabulated_grads[std::make_pair(radial_order, radial_family)];

    // Test whether the computed values match the tabulated values to
    // the specified accuracy.
    for (const auto & t : val_table)
      LIBMESH_ASSERT_FP_EQUAL(t.val, phi[t.i][t.qp], 1.e-4);
    for (const auto & t : grad_table)
      LIBMESH_ASSERT_FP_EQUAL(0., (dphi[t.i][t.qp] - t.grad).norm_sq(), 1.e-4);

    // Make sure there are actually reference values
    if (val_table.empty())
      libmesh_error_msg("No tabulated values found!");
    if (grad_table.empty())
      libmesh_error_msg("No tabulated gradients found!");

#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
  }

  void setUp()
  {
    ////////////////////////////////////////////////////////////////////////////////
    // Arbitrarily selected LEGENDRE reference values.
    ////////////////////////////////////////////////////////////////////////////////
    tabulated_vals[std::make_pair(FIRST, LEGENDRE)] =
      {
        {0,9,0.0550016}, {1,5,0.41674},   {2,8,0.0147376}, {3,7,0.111665},
        {4,6,0.0737011}, {5,1,0.0803773}, {6,11,0.275056}, {7,15,0.021537}
      };

    tabulated_vals[std::make_pair(SECOND, LEGENDRE)] =
      {
        {0,3,0.0425633}, {1,6,0.0343526}, {2,2,0.158848}, {3,7,0.128206},
        {4,12,0.220829}, {5,5,-0.136549}, {6,0,0.0149032}, {7,10,-0.0334936},
        {8,16,0.00399329}, {9,18,-0.00209733}, {10,2,0.0556194}, {11,17,-0.000561977}
      };

    tabulated_vals[std::make_pair(THIRD, LEGENDRE)] =
      {
        {12,9,0.136657}, {13,6,0.175034}, {14,20,-0.0011016}, {15,15,0.0428935}
      };

    tabulated_vals[std::make_pair(FOURTH, LEGENDRE)] =
      {
        {16,3,0.00826618}, {17,10,-0.547813}, {18,14,0.311004}, {19,27,-0.00192085}
      };

    ////////////////////////////////////////////////////////////////////////////////
    // Arbitrarily selected LEGENDRE reference gradients.
    ////////////////////////////////////////////////////////////////////////////////
    tabulated_grads[std::make_pair(FIRST, LEGENDRE)] =
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

    tabulated_grads[std::make_pair(SECOND, LEGENDRE)] =
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

    tabulated_grads[std::make_pair(THIRD, LEGENDRE)] =
      {
        {12,9,Point(0.0536552, 0.200244, -0.0849541)},
        {13,6,Point(-0.0921697, 0.0921697, 0.461056)},
        {14,20,Point(2.35811e-05, -8.80059e-05, 3.58959e-05)},
        {15,15,Point(-0.0386352, 0.0103523, -0.0197323)},
      };

    tabulated_grads[std::make_pair(FOURTH, LEGENDRE)] =
      {
        {16,3,Point(-0.0190603, 0.0051072, 0.308529)},
        {17,10,Point(0.244125, -0.244125, 0.140907)},
        {18,14,Point(-0.0985844, 0.0985844, -0.502591)},
        {19,27,Point(0.000115647, -3.09874e-05, 4.30775e-05)}
      };

    ////////////////////////////////////////////////////////////////////////////////
    // Arbitrarily selected JACOBI_20_00 reference values.
    ////////////////////////////////////////////////////////////////////////////////
    tabulated_vals[std::make_pair(FIRST, JACOBI_20_00)] =
      {
        // These values are the same as Legendre
        {0,9,0.0550016}, {1,5,0.41674},   {2,8,0.0147376}, {3,7,0.111665},
        // These values are different from Legendre
        {4,6,0.147402}, {5,1,0.160755}, {6,11,0.550112}, {7,15,0.043074}
      };

    tabulated_vals[std::make_pair(SECOND, JACOBI_20_00)] =
      {
        // These values are the same as Legendre
        {0,3,0.0425633}, {1,6,0.0343526}, {2,2,0.158848}, {3,7,0.128206},
        // These values are different from Legendre
        {4,12,0.441658}, {5,5,-0.193445}, {6,0,0.0298063}, {7,10,-0.0279114},
        {8,16,0.00798659}, {9,18,0.0320146}, {10,2,0.111239}, {11,17,0.00857829}
      };

    tabulated_vals[std::make_pair(THIRD, JACOBI_20_00)] =
      {
        {12,9,0.0837876}, {13,6,0.350068}, {14,20,0.0244336}, {15,15,0.0181532}
      };

    tabulated_vals[std::make_pair(FOURTH, JACOBI_20_00)] =
      {
        {16,3,0.0165324}, {17,10,-0.720085}, {18,14,0.0777511}, {19,27,0.0453842}
      };

    ////////////////////////////////////////////////////////////////////////////////
    // Arbitrarily selected JACOBI_20_00 reference gradients.
    ////////////////////////////////////////////////////////////////////////////////
    tabulated_grads[std::make_pair(FIRST, JACOBI_20_00)] =
      {
        // These values are the same as Legendre
        {0,9,Point(-0.0429458, -0.0115073, 0.)},
        {1,5,Point(0.177013, -0.177013, -0.483609)},
        {2,8,Point(0.0115073, 0.0115073, 0.00842393)},
        {3,7,Point(-0.177013, 0.0474305, 0.)},
        // These values are different from Legendre
        {4,6,Point(-0.0626101, -0.233664,   0.2005)},
        {5,1,Point(0.0948382, -0.0948382,  1.74583)},
        {6,11,Point(0.115093, 0.115093, -0.22502)},
        {7,15,Point(-0.0070761, 0.00189603, 0.000223144)}
      };

    tabulated_grads[std::make_pair(SECOND, JACOBI_20_00)] =
      {
        // These values are the same as Legendre
        {0,3,Point(-0.0959817, -0.0959817, 0.0702635)},
        {1,6,Point(0.0625228, -0.0625228, 0.0457699)},
        {2,2,Point(0.358209, 0.0959817, -2.77556e-17)},
        {3,7,Point(-0.233338, 0.0625228, -6.93889e-17)},
        // These values are different from Legendre
        {4,12,Point(-0.0646142, -0.0646142, -0.145954)},
        {5,5,Point(0.352076, 0.0943384, -0.233431)},
        {6,0,Point(0.0672144, -0.01801, 0.577179)},
        {7,10,Point(-0.0330195, 0.0330195, 0.00373942)},
        {8,16,Point(0.000886435, 0.000886435, 0.000667355)},
        {9,18,Point(0.00355332, 0.000952108, 0.000319882)},
        {10,2,Point(-0.0672144, 0.0672144,  1.97043)},
        {11,17,Point(-0.000952108, 0.000952108, 0.000782704)}
      };

    tabulated_grads[std::make_pair(THIRD, JACOBI_20_00)] =
      {
        {12,9,Point(0.0328972,0.122774,-0.222472)},
        {13,6,Point(-0.184339,0.184339,0.922112)},
        {14,20,Point(-0.000523034,0.00195199,0.000121819)},
        {15,15,Point(-0.016351,0.00438123,0.0404817)}
      };

    tabulated_grads[std::make_pair(FOURTH, JACOBI_20_00)] =
      {
        {16,3,Point(-0.0381206,0.0102144,0.617059)},
        {17,10,Point(0.320895,-0.320895,0.641729)},
        {18,14,Point(-0.0246461,0.0246461,-0.300588)},
        {19,27,Point(-0.0027324,0.000732145,0.000328402)}
      };

    ////////////////////////////////////////////////////////////////////////////////
    // Arbitrarily selected JACOBI_30_00 reference values.
    ////////////////////////////////////////////////////////////////////////////////
    tabulated_vals[std::make_pair(FIRST, JACOBI_30_00)] =
      {
        // These values are the same as Legendre
        {0,9,0.0550016}, {1,5,0.41674},   {2,8,0.0147376}, {3,7,0.111665},
        // These values are different from Legendre
        {4,6,0.184253}, {5,1,0.200943}, {6,11,0.68764}, {7,15,0.0538426}
      };

    tabulated_vals[std::make_pair(SECOND, JACOBI_30_00)] =
      {
        // These values are the same as Legendre
        {0,3,0.0425633}, {1,6,0.0343526}, {2,2,0.158848}, {3,7,0.128206},
        // These values are different from Legendre
        {4,12,0.552072}, {5,5,-0.211652}, {6,0,0.0372579}, {7,10,-0.0167468},
        {8,16,0.00998323}, {9,18,0.0597236}, {10,2,0.139049}, {11,17,0.0160029}
      };

    tabulated_vals[std::make_pair(THIRD, JACOBI_30_00)] =
      {
        {12,9,0.0469849}, {13,6,0.437585}, {14,20,0.0450821}, {15,15,0.0469849}
      };

    tabulated_vals[std::make_pair(FOURTH, JACOBI_30_00)] =
      {
        {16,3,0.0206655}, {17,10,-0.748341}, {18,14,0}, {19,27,0.115995}
      };

    ////////////////////////////////////////////////////////////////////////////////
    // Arbitrarily selected JACOBI_30_00 reference gradients.
    ////////////////////////////////////////////////////////////////////////////////
    tabulated_grads[std::make_pair(FIRST, JACOBI_30_00)] =
      {
        // These values are the same as Legendre
        {0,9,Point(-0.0429458, -0.0115073, 0.)},
        {1,5,Point(0.177013, -0.177013, -0.483609)},
        {2,8,Point(0.0115073, 0.0115073, 0.00842393)},
        {3,7,Point(-0.177013, 0.0474305, 0.)},
        // These values are different from Legendre
        {4,6,Point(-0.0782626,-0.29208,0.250625)},
        {5,1,Point(0.118548,-0.118548,2.18229)},
        {6,11,Point(0.143866,0.143866,-0.281275)},
        {7,15,Point(-0.00884512,0.00237004,0.00027893)}
      };

    tabulated_grads[std::make_pair(SECOND, JACOBI_30_00)] =
      {
        // These values are the same as Legendre
        {0,3,Point(-0.0959817, -0.0959817, 0.0702635)},
        {1,6,Point(0.0625228, -0.0625228, 0.0457699)},
        {2,2,Point(0.358209, 0.0959817, -2.77556e-17)},
        {3,7,Point(-0.233338, 0.0625228, -6.93889e-17)},
        // These values are different from Legendre
        {4,12,Point(-0.0807678,-0.0807678,-0.182443)},
        {5,5,Point(0.385213,0.103218,-0.175079)},
        {6,0,Point(0.0840179,-0.0225125,0.721474)},
        {7,10,Point(-0.0198117,0.0198117,0.0357373)},
        {8,16,Point(0.00110804,0.00110804,0.000834194)},
        {9,18,Point(0.00662875,0.00177617,0.000482244)},
        {10,2,Point(-0.0840179,0.0840179,2.46303)},
        {11,17,Point(-0.00177617,0.00177617,0.00142946)},
      };

    tabulated_grads[std::make_pair(THIRD, JACOBI_30_00)] =
      {
        {12,9,Point(0.0184475,0.0688471,-0.254748)},
        {13,6,Point(-0.230424,0.230424,1.15264)},
        {14,20,Point(-0.000965042,0.00360159,0.000183379)},
        {15,15,Point(-0.0423204,0.0113397,0.12514)},
      };

    tabulated_grads[std::make_pair(FOURTH, JACOBI_30_00)] =
      {
        {16,3,Point(-0.0476508,0.012768,0.771323)},
        {17,10,Point(0.333487,-0.333487,1.01421)},
        {18,14,Point(0,0,0)},
        {19,27,Point(-0.00698362,0.00187126,0.000667664)},
      };
  }

  void tearDown() {}

};

CPPUNIT_TEST_SUITE_REGISTRATION( InfFERadialTest );
