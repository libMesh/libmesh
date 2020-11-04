// libmesh includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/inf_elem_builder.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"
#include "libmesh/inf_fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/jacobi_polynomials.h"

// unit test includes
#include "test_comm.h"

// C++ includes
#include <map>
#include <vector>

#include "libmesh_cppunit.h"
#include "libmesh/type_tensor.h"
#include "libmesh/fe.h"

#include "libmesh/face_tri6.h"

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
  CPPUNIT_TEST( testInfQuants );
  CPPUNIT_TEST( testSides );
  CPPUNIT_TEST( testInfQuants_numericDeriv );
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

  Point base_point(const Point physical_point, const Elem * inf_elem)
  {
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    libmesh_assert(inf_elem);

    // The strategy is:
    // compute the intersection of the line
    // physical_point - origin with the base element,
    // find its internal coordinatels using FEMap::inverse_map():
    // The radial part can then be computed directly later on.

    // 1.)
    // build a base element to do the map inversion in the base face
    std::unique_ptr<const Elem> base_elem = InfFEBase::build_elem (inf_elem);

    // the origin of the infinite element
    const Point o = inf_elem->origin();
    // 2.)
    // Now find the intersection of a plane represented by the base
    // element nodes and the line given by the origin of the infinite
    // element and the physical point.
    Point intersection;

    const Point xi ( base_elem->point(1) - base_elem->point(0));
    const Point eta( base_elem->point(2) - base_elem->point(0));
    const Point zeta( physical_point - o);

    // normal vector of the base elements plane
    Point n(xi.cross(eta));
    Real c_factor = (base_elem->point(0) - o)*n/(zeta*n) - 1.;

    // Check whether the above system is ill-posed.  It should
    // only happen when \p physical_point is not in \p
    // inf_elem. In this case, any point that is not in the
    // element is a valid answer.
    if (libmesh_isinf(c_factor))
      return Point(0., 0., -2.);

    // Compute the intersection with
    // {intersection} = {physical_point} + c*({physical_point}-{o}).
    intersection.add_scaled(physical_point,1.+c_factor);
    intersection.add_scaled(o,-c_factor);

    if (!base_elem->has_affine_map())
      {
        unsigned int iter_max = 20;

        // the number of shape functions needed for the base_elem
        unsigned int n_sf = FE<2,LAGRANGE>::n_shape_functions(base_elem->type(),base_elem->default_order());

        // guess base element coordinates: p=xi,eta,0
        // in first iteration, never run with 'secure' or
        // 'extra_checks' to avoid false warnings.
        Point ref_point= FEMap::inverse_map(2, base_elem.get(), intersection,
                                            TOLERANCE, false, false);

        // Newton iteration
        for(unsigned int it=0; it<iter_max; ++it)
          {
            // Get the shape function and derivative values at the reference coordinate
            // phi.size() == dphi.size()
            Point dxyz_dxi;
            Point dxyz_deta;

            Point intersection_guess;
            for(unsigned int i=0; i<n_sf; ++i)
              {

                intersection_guess += base_elem->node_ref(i) * FE<2,LAGRANGE>::shape(base_elem->type(),
                                                                                     base_elem->default_order(),
                                                                                     i,
                                                                                     ref_point);

                dxyz_dxi += base_elem->node_ref(i) * FE<2,LAGRANGE>::shape_deriv(base_elem->type(),
                                                                                 base_elem->default_order(),
                                                                                 i,
                                                                                 0, // d()/dxi
                                                                                 ref_point);

                dxyz_deta+= base_elem->node_ref(i) * FE<2,LAGRANGE>::shape_deriv(base_elem->type(),
                                                                                 base_elem->default_order(),
                                                                                 i,
                                                                                 1, // d()/deta
                                                                                 ref_point);
              } // for i

            TypeVector<Real> F(physical_point + c_factor*(physical_point-o) - intersection_guess);

            TypeTensor<Real> J;
            J(0,0) = (physical_point-o)(0);
            J(0,1) = -dxyz_dxi(0);
            J(0,2) = -dxyz_deta(0);
            J(1,0) = (physical_point-o)(1);
            J(1,1) = -dxyz_dxi(1);
            J(1,2) = -dxyz_deta(1);
            J(2,0) = (physical_point-o)(2);
            J(2,1) = -dxyz_dxi(2);
            J(2,2) = -dxyz_deta(2);

            // delta will be the newton step
            TypeVector<Real> delta;
            J.solve(F,delta);

            // check for convergence
            Real tol = std::min( TOLERANCE*0.01, TOLERANCE*base_elem->hmax() );
            if ( delta.norm() < tol )
              {
                // newton solver converged, now make sure it converged to a point on the base_elem
                if (base_elem->contains_point(intersection_guess,TOLERANCE*0.1))
                  {
                    intersection = intersection_guess;
                  }
                else
                  {
                    err<<" No! This inverse_map failed!"<<std::endl;
                  }
                break; // break out of 'for it'
              }
            else
              {
                c_factor     -= delta(0);
                ref_point(0) -= delta(1);
                ref_point(1) -= delta(2);
              }

          }

      }

    return intersection;
#else
    // lets make the compilers happy:
    return Point(0.,0.,-2.);
#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
  }

  typedef std::vector<TabulatedVal> TabulatedVals;
  typedef std::vector<TabulatedGrad> TabulatedGrads;
  typedef std::pair<Order, FEFamily> KeyType;
  std::map<KeyType, TabulatedVals> tabulated_vals;
  std::map<KeyType, TabulatedGrads> tabulated_grads;

  void testSides()
  {
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_cube
      (mesh,
       /*nx=*/2, /*ny=*/2, /*nz=*/2,
       /*xmin=*/-2., /*xmax=*/2.,
       /*ymin=*/-1., /*ymax=*/1.,
       /*zmin=*/-2., /*zmax=*/2.,
       PRISM18);

    const unsigned int n_fem =mesh.n_elem();

    InfElemBuilder builder(mesh);
    builder.build_inf_elem();

    // Get pointer to the first infinite Elem.
    Elem * infinite_elem = mesh.elem_ptr(n_fem+3);
    if (!infinite_elem || !infinite_elem->infinite())
      libmesh_error_msg("Error setting Elem pointer.");

    // We will construct FEs, etc. of the same dimension as the mesh elements.
    auto dim = mesh.mesh_dimension();

    //FEType fe_type(/*Order*/FIRST,
    FEType fe_type(/*Order*/SECOND,
                   /*FEFamily*/LAGRANGE,
                   /*radial order*/FIRST,
                   /*radial_family*/LAGRANGE,
                   /*inf_map*/CARTESIAN);


    // Construct FE, quadrature rule, etc.
    std::unique_ptr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
    std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
    std::unique_ptr<FEBase> fe2 (FEBase::build(dim-1, fe_type));
    QGauss qrule (dim-1, fe_type.default_quadrature_order());
    inf_fe->attach_quadrature_rule(&qrule);
    fe->attach_quadrature_rule(&qrule);
    fe2->attach_quadrature_rule(&qrule);

    //find the FE neighbour of \p infinite_elem:
    unsigned int side_num=7;
    MeshBase::const_element_iterator           el = mesh.elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.elements_end();
    Point side_pt = infinite_elem->side_ptr(0)->centroid();

    Elem * finite_elem = infinite_elem->neighbor_ptr(0);
    for(unsigned int i=0; i<finite_elem->n_sides(); ++i)
      {
        if (finite_elem->side_ptr(i)->contains_point(side_pt))
          {
            side_num = i;
            break;
          }
      }
    libmesh_assert(side_num < 7);
    // in particular, side_num = 0 (used later)

    // construct a 'side' element (specific to the current geometry):
    Tri6 side_elem;
    side_elem.set_node(0)=finite_elem->node_ptr(0);
    side_elem.set_node(1)=finite_elem->node_ptr(2);
    side_elem.set_node(2)=finite_elem->node_ptr(1);
    side_elem.set_node(3)=finite_elem->node_ptr(8);
    side_elem.set_node(4)=finite_elem->node_ptr(7);
    side_elem.set_node(5)=finite_elem->node_ptr(6);

    const std::vector<Point>& i_qpoint = inf_fe->get_xyz();
    const std::vector<Real> & i_weight = inf_fe->get_Sobolev_weight();
    const std::vector<Real> & i_JxW = inf_fe->get_JxWxdecay_sq();
    const std::vector<std::vector<Real> >& i_phi  = inf_fe->get_phi();
    const std::vector<Point> &                i_normal = inf_fe->get_normals();
    const std::vector<std::vector<Point> >& i_tangents = inf_fe->get_tangents();

    const std::vector<Point>& f_qpoint = fe->get_xyz();
    const std::vector<Real> & f_weight = fe->get_Sobolev_weight();
    const std::vector<Real> & f_JxW = fe->get_JxWxdecay_sq();
    const std::vector<std::vector<Real> >& f_phi  = fe->get_phi();
    const std::vector<Point> &                f_normal = fe->get_normals();
    const std::vector<std::vector<Point> >& f_tangents = fe->get_tangents();

    const std::vector<Point>& s_qpoint = fe2->get_xyz();
    const std::vector<Real> & s_JxW = fe2->get_JxWxdecay_sq();
    const std::vector<std::vector<Real> >& s_phi  = fe2->get_phi();

    // Reinit on infinite elem
    inf_fe->reinit(infinite_elem, (unsigned)0);
    fe->reinit(finite_elem, side_num);
    fe2->reinit(&side_elem);

    LIBMESH_ASSERT_FP_EQUAL(i_weight.size(), f_weight.size(), TOLERANCE);
    for(unsigned int qp =0 ; qp < i_weight.size() ; ++qp)
      {
        LIBMESH_ASSERT_FP_EQUAL(i_qpoint[qp](0), f_qpoint[qp](0), TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(i_qpoint[qp](1), f_qpoint[qp](1), TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(i_qpoint[qp](2), f_qpoint[qp](2), TOLERANCE);

        LIBMESH_ASSERT_FP_EQUAL(s_qpoint[qp](0), f_qpoint[qp](0), TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(s_qpoint[qp](1), f_qpoint[qp](1), TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(s_qpoint[qp](2), f_qpoint[qp](2), TOLERANCE);

        LIBMESH_ASSERT_FP_EQUAL(i_weight[qp], f_weight[qp], TOLERANCE); // this is both 1
        LIBMESH_ASSERT_FP_EQUAL(i_JxW[qp]   , f_JxW[qp]   , TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(s_JxW[qp]   , f_JxW[qp]   , TOLERANCE);

        // initialize index maps: (specific to current geometry:
        unsigned int inf_index[] ={0, 1, 2, 6, 7, 8};
        unsigned int fe_index[]  ={0, 2, 1, 8, 7, 6};
        unsigned int s_index[]   ={0, 1, 2, 3, 4, 5};
        for (unsigned int i=0; i<6; ++i)
          {
            LIBMESH_ASSERT_FP_EQUAL(i_phi[inf_index[i]][qp], f_phi[fe_index[i]][qp], TOLERANCE);
            LIBMESH_ASSERT_FP_EQUAL(i_phi[inf_index[i]][qp], s_phi[s_index[i]][qp], TOLERANCE);
          }

        // note that infinite element normals point OUTWARD of the element
        LIBMESH_ASSERT_FP_EQUAL(i_normal[qp](0), f_normal[qp](0), TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(i_normal[qp](1), f_normal[qp](1), TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(i_normal[qp](2), f_normal[qp](2), TOLERANCE);

        for (unsigned int i=0; i< i_tangents[0].size(); ++i)
          {
            LIBMESH_ASSERT_FP_EQUAL(i_tangents[qp][i](0), f_tangents[qp][i](0), TOLERANCE);
            LIBMESH_ASSERT_FP_EQUAL(i_tangents[qp][i](1), f_tangents[qp][i](1), TOLERANCE);
            LIBMESH_ASSERT_FP_EQUAL(i_tangents[qp][i](2), f_tangents[qp][i](2), TOLERANCE);
          }

      }

#endif
  }

  void testInfQuants ()
  {
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    // std::cout << "Called testSingleOrder with radial_order = "
    //           << radial_order
    //           << std::endl;

    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_cube
      (mesh,
       /*nx=*/3, /*ny=*/2, /*nz=*/5,
       /*xmin=*/0., /*xmax=*/2.,
       /*ymin=*/-1., /*ymax=*/1.,
       /*zmin=*/-4., /*zmax=*/-2.,
       TET10);

    InfElemBuilder::InfElemOriginValue com_x;
    InfElemBuilder::InfElemOriginValue com_y;
    InfElemBuilder::InfElemOriginValue com_z;
    com_x.first=true;
    com_y.first=true;
    com_z.first=true;
    com_x.second=0.0;
    com_y.second=-.5;
    com_z.second=0.0;

    const unsigned int n_fem =mesh.n_elem();

    InfElemBuilder builder(mesh);
    builder.build_inf_elem(com_x, com_y, com_z,
                           false,  false,  false,
                           true, libmesh_nullptr);

    // Get pointer to the first infinite Elem.
    const Elem * infinite_elem = mesh.elem_ptr(n_fem);
    if (!infinite_elem || !infinite_elem->infinite())
      libmesh_error_msg("Error setting Elem pointer.");

    // We will construct FEs, etc. of the same dimension as the mesh elements.
    auto dim = mesh.mesh_dimension();

    FEType fe_type(/*Order*/FIRST,
                   /*FEFamily*/LAGRANGE,
                   /*radial order*/FIRST,
                   /*radial_family*/LAGRANGE,
                   /*inf_map*/CARTESIAN);

    // Construct FE, quadrature rule, etc.
    std::unique_ptr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    inf_fe->attach_quadrature_rule(&qrule);

    const std::vector<Point>& _qpoint = inf_fe->get_xyz();
    const std::vector<Real> & _weight = inf_fe->get_Sobolev_weight();
    const std::vector<RealGradient>& dweight = inf_fe->get_Sobolev_dweight();
    const std::vector<Real>& dxidx = inf_fe->get_dxidx();
    const std::vector<Real>& dxidy = inf_fe->get_dxidy();
    const std::vector<Real>& dxidz = inf_fe->get_dxidz();
    const std::vector<Real>& detadx = inf_fe->get_detadx();
    const std::vector<Real>& detady = inf_fe->get_detady();
    const std::vector<Real>& detadz = inf_fe->get_detadz();
    const std::vector<Real>& dzetadx = inf_fe->get_dzetadx();
    const std::vector<Real>& dzetady = inf_fe->get_dzetady();
    const std::vector<Real>& dzetadz = inf_fe->get_dzetadz();
    const std::vector<RealGradient>& dphase  = inf_fe->get_dphase();
    // Reinit on infinite elem
    inf_fe->reinit(infinite_elem);

    for(unsigned int qp =0 ; qp < _weight.size() ; ++qp)
      {
        /**
         * dphase = r/|r| where r is the vector from the 'origin' to  the point of interest
         * Sob. weight = b*b/(r*r)  with b being the vector of the elements base to the origin.
         *                          Since Sobolev weight is 1 for finite elements, it gives a
         *                          continuous function.
         * its derivative, thus, -2*b*b*r/(r**4)
         */

        const Point b = base_point(_qpoint[qp], infinite_elem) - infinite_elem->origin();
        const Point p = _qpoint[qp] - infinite_elem ->origin();
        const Point rp = FEInterface::inverse_map(3, fe_type, infinite_elem, _qpoint[qp] );
        const Point again_qp = FEInterface::map(3, fe_type, infinite_elem, rp);
        //const Point rb = FEInterface::inverse_map(3, fe_type, infinite_elem, b+infinite_elem->origin());
        const Real v=rp(2);

        // check that map() does the opposite from inverse_map
        LIBMESH_ASSERT_FP_EQUAL(again_qp(0), _qpoint[qp](0), TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(again_qp(1), _qpoint[qp](1), TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(again_qp(2), _qpoint[qp](2), TOLERANCE);

        const Point normal(dzetadx[qp],
                           dzetady[qp],
                           dzetadz[qp]);

        // check that dweight is in direction of base elements normal
        // (holds when the base_elem has an affine map)
        LIBMESH_ASSERT_FP_EQUAL(normal*dweight[qp], -normal.norm()*dweight[qp].norm(), TOLERANCE);

        if (rp(2) < -.85)
          // close to the base, dphase is normal
          LIBMESH_ASSERT_FP_EQUAL(normal*dphase[qp], normal.norm()*dphase[qp].norm(), 3e-4);
        else if (rp(2) > .85)
          // 'very far away' dphase goes in radial direction
          LIBMESH_ASSERT_FP_EQUAL(p*dphase[qp]/p.norm(), dphase[qp].norm(), 1e-4);

        // check that the mapped radial coordinate corresponds to a/r
        //   - b: radius of the FEM-region (locally, i.e. distance of the projected base point)
        //   - p: distance of the point from origin
        LIBMESH_ASSERT_FP_EQUAL(b.norm()/p.norm(), (1.-v)/2., TOLERANCE);

        // this is how weight is defined;
        LIBMESH_ASSERT_FP_EQUAL(b.norm_sq()/p.norm_sq(), _weight[qp], TOLERANCE);

        const Point e_xi(dxidx[qp], dxidy[qp], dxidz[qp]);
        const Point e_eta(detadx[qp], detady[qp], detadz[qp]);

        LIBMESH_ASSERT_FP_EQUAL(0., b*e_xi,TOLERANCE);
        LIBMESH_ASSERT_FP_EQUAL(0., b*e_eta,TOLERANCE);

      }

#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
  }

  void testInfQuants_numericDeriv ()
  {
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
    ReplicatedMesh mesh(*TestCommWorld);
    MeshTools::Generation::build_sphere
      (mesh, /*rad*/ 1,
       /* nr */ 2, /*type*/ HEX27);

    InfElemBuilder::InfElemOriginValue com_x;
    InfElemBuilder::InfElemOriginValue com_y;
    InfElemBuilder::InfElemOriginValue com_z;
    com_x.first=true;
    com_y.first=true;
    com_z.first=true;
    com_x.second=0.7;
    com_y.second=-0.4;
    com_z.second=0.0;

    const unsigned int n_fem =mesh.n_elem();

    InfElemBuilder builder(mesh);
    builder.build_inf_elem(com_x, com_y, com_z,
                           false,  false,  false,
                           true, libmesh_nullptr);

    // Get pointer to the first infinite Elem.
    Elem * infinite_elem = mesh.elem_ptr(n_fem+1);
    if (!infinite_elem || !infinite_elem->infinite())
      libmesh_error_msg("Error setting Elem pointer.");
    // lets overemphasize that the base element has a non-affine map:
    for (unsigned int n=8; n<12; ++n)
      {
        Node* node=infinite_elem->node_ptr(n);
        // shift base-points on edges towards the center.
        // This leads to strong 'non-affine effects'
        *node -= 0.1*(*node-infinite_elem->origin()).unit();
      }

    // We will construct FEs, etc. of the same dimension as the mesh elements.
    auto dim = mesh.mesh_dimension();

    FEType fe_type(/*Order*/FIRST,
                   /*FEFamily*/LAGRANGE,
                   /*radial order*/FIRST,
                   /*radial_family*/LAGRANGE,
                   /*inf_map*/CARTESIAN);

    // Reinit on infinite elem
    //
    unsigned int num_pt=10;
    std::vector<Point> points(2*num_pt);
    points[0]=Point(-0.7, -0.5, -0.9);
    points[1]=Point(-0.1,  0.9, -0.9);
    points[2]=Point(-0.7, -0.5, -0.4);
    points[3]=Point(-0.1,  0.9, -0.4);
    points[4]=Point(-0.7, -0.5, -0.2);
    points[5]=Point(-0.1,  0.9, -0.2);
    points[6]=Point(-0.7, -0.5,  0.1);
    points[7]=Point(-0.1,  0.9,  0.1);
    points[8]=Point(-0.7, -0.5,  0.6);
    points[9]=Point(-0.1,  0.9,  0.6);

    //  Check derivatives along radial direction:
    Point delta(0.,0.,1e-3);
    for (unsigned int i=num_pt; i<2*num_pt; ++i)
      points[i]=points[i-num_pt]+delta;

    std::unique_ptr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
    const std::vector<Point> &                  q_point = inf_fe->get_xyz();
    const std::vector<Real>         &             sob_w = inf_fe->get_Sobolev_weightxR_sq();
    const std::vector<RealGradient> &            dsob_w = inf_fe->get_Sobolev_dweightxR_sq();
    const std::vector<Real> &                   sob_now = inf_fe->get_Sobolev_weight();
    const std::vector<RealGradient>&           dsob_now = inf_fe->get_Sobolev_dweight();
    const std::vector<RealGradient>&            dphase  = inf_fe->get_dphase();
    const std::vector<std::vector<RealGradient> >& dphi = inf_fe->get_dphi();
    const std::vector<std::vector<Real> >&         phi  = inf_fe->get_phi();
    const std::vector<std::vector<RealGradient> >& dphi_w = inf_fe->get_dphi_over_decayxR();
    const std::vector<std::vector<Real> >&         phi_w  = inf_fe->get_phi_over_decayxR();
    inf_fe->reinit(infinite_elem,&points);

    for(unsigned int qp =0 ; qp < num_pt ; ++qp)
      {
        const Point dxyz(q_point[qp+num_pt]-q_point[qp]);
        const Point b_i= base_point(q_point[qp], infinite_elem) - infinite_elem->origin();
        const Point b_o= base_point(q_point[qp+num_pt], infinite_elem) - infinite_elem->origin();

        LIBMESH_ASSERT_FP_EQUAL((b_i-b_o).norm_sq(), 0, TOLERANCE);

        Real weight_i = b_i.norm_sq()/(q_point[qp]-infinite_elem->origin()).norm_sq();
        Real weight_o = b_o.norm_sq()/(q_point[qp+num_pt]-infinite_elem->origin()).norm_sq();
        const Real phase_i = (q_point[qp]-infinite_elem->origin()).norm() - b_i.norm();
        const Real phase_o = (q_point[qp+num_pt]-infinite_elem->origin()).norm() - b_o.norm();

        Real tolerance = std::abs((dphase[qp+num_pt]-dphase[qp])*dxyz)+1e-10;
        Real deriv_mean = (dphase[qp]*dxyz + dphase[qp+num_pt]*dxyz)*0.5;
        LIBMESH_ASSERT_FP_EQUAL(phase_o - phase_i, deriv_mean, tolerance*.5);

        tolerance = std::abs((dsob_now[qp+num_pt]-dsob_now[qp])*dxyz)+1e-10;
        deriv_mean = (dsob_now[qp]*dxyz +  dsob_now[qp+num_pt]*dxyz)* 0.5;
        LIBMESH_ASSERT_FP_EQUAL(sob_now[qp+num_pt] - sob_now[qp], deriv_mean, tolerance*.5);

        LIBMESH_ASSERT_FP_EQUAL(sob_w[qp+num_pt]*weight_o - sob_w[qp]*weight_i, dsob_w[qp]*dxyz*weight_i, tolerance);

        for (unsigned int i=0; i< phi.size(); ++i)
          {
            tolerance = std::abs((dphi[i][qp+num_pt]-dphi[i][qp])*dxyz)+1e-10;
            deriv_mean = (dphi[i][qp]*dxyz + dphi[i][qp+num_pt]*dxyz)*0.5;
            LIBMESH_ASSERT_FP_EQUAL(phi[i][qp+num_pt] - phi[i][qp], deriv_mean, tolerance*.5);

            deriv_mean = 0.5*(dphi_w[i][qp]*dxyz*sqrt(weight_i) +dphi_w[i][qp+num_pt]*dxyz*sqrt(weight_o));
            LIBMESH_ASSERT_FP_EQUAL(phi_w[i][qp+num_pt]*sqrt(weight_o) - phi_w[i][qp]*sqrt(weight_i),
                                    deriv_mean, tolerance*.5);
          }

      }

    //  Check derivatives along angular direction:
    points[0 ]=Point(-0.7, -0.5, -0.9);
    points[2 ]=Point(-0.1,  0.9, -0.9);
    points[4 ]=Point(-0.7, -0.5, -0.4);
    points[6 ]=Point(-0.1,  0.9, -0.4);
    points[8 ]=Point(-0.7, -0.5, -0.2);
    points[10]=Point(-0.1,  0.9, -0.2);
    points[12]=Point(-0.7, -0.5,  0.1);
    points[14]=Point(-0.1,  0.9,  0.1);
    points[16]=Point(-0.7, -0.5,  0.6);
    points[18]=Point(-0.1,  0.9,  0.6);

    delta = Point(1.2e-4,-2.7e-4,0.);
    for (unsigned int i=0; i<2*num_pt; i+=2)
      points[i+1]=points[i]+delta;

    const std::vector<Real>& dzetadx = inf_fe->get_dzetadx();
    const std::vector<Real>& dzetady = inf_fe->get_dzetady();
    const std::vector<Real>& dzetadz = inf_fe->get_dzetadz();

    inf_fe->reinit(infinite_elem,&points);

    for(unsigned int qp =0 ; qp < 2*num_pt ; qp+=2)
      {
        const Point dxyz(q_point[qp+1]-q_point[qp]);
        const Point b_i= base_point(q_point[qp], infinite_elem) - infinite_elem->origin();
        const Point b_o= base_point(q_point[qp+1], infinite_elem) - infinite_elem->origin();
        Real weight_i = b_i.norm_sq()/(q_point[qp]-infinite_elem->origin()).norm_sq();
        Real weight_o = b_o.norm_sq()/(q_point[qp+1]-infinite_elem->origin()).norm_sq();
        const Real phase_i = (q_point[qp]-infinite_elem->origin()).norm() - b_i.norm();
        const Real phase_o = (q_point[qp+1]-infinite_elem->origin()).norm() - b_o.norm();

        LIBMESH_ASSERT_FP_EQUAL(weight_o ,weight_i, TOLERANCE);

        Point normal(dzetadx[qp],
                     dzetady[qp],
                     dzetadz[qp]);

        Real a_i= (q_point[qp]-infinite_elem->origin()).norm()*0.5*(1.-points[qp](2));

        LIBMESH_ASSERT_FP_EQUAL(a_i, b_i.norm(), TOLERANCE*TOLERANCE);

        // the elements normal should be orthogonal to dxyz, but in practice deviates
        // approx. 1 degree. So we must account for this error in direction.
        Real err_direct = 1.5*std::abs(dxyz*normal)/(dxyz.norm()*normal.norm());

        Real tolerance = std::abs((dphase[qp+1]-dphase[qp])*dxyz)*0.5 + err_direct*dxyz.norm()*dphase[qp].norm();
        Real deriv_mean = (dphase[qp] + dphase[qp+1])*dxyz*0.5;
        LIBMESH_ASSERT_FP_EQUAL(phase_o - phase_i, deriv_mean, tolerance );

        LIBMESH_ASSERT_FP_EQUAL(normal.cross(dsob_now[qp]).norm_sq(), 0.0, TOLERANCE*TOLERANCE);

        tolerance = std::abs((dsob_now[qp+1]-dsob_now[qp])*dxyz)*.5+1e-10 + err_direct*dxyz.norm()*dsob_now[qp].norm();
        deriv_mean = (dsob_now[qp]*dxyz + dsob_now[qp+1]*dxyz)*0.5;
        LIBMESH_ASSERT_FP_EQUAL(sob_now[qp+1] - sob_now[qp], deriv_mean, tolerance);

        deriv_mean = (dsob_w[qp]*dxyz*weight_i +dsob_w[qp+1]*dxyz*weight_o)*0.5;
        LIBMESH_ASSERT_FP_EQUAL(sob_w[qp+1]*weight_o - sob_w[qp]*weight_i, deriv_mean, tolerance);

        for (unsigned int i=0; i< phi.size(); ++i)
          {
            tolerance = std::abs((dphi[i][qp+1]-dphi[i][qp])*dxyz)*0.5+1e-10 + err_direct*dxyz.norm()*dphi[i][qp].norm();
            deriv_mean = (dphi[i][qp]*dxyz + dphi[i][qp+1]*dxyz)*.5;
            LIBMESH_ASSERT_FP_EQUAL(phi[i][qp+1] - phi[i][qp], deriv_mean, tolerance);
          }
      }

#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
  }

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
    //           << " radial_family: "<<radial_family
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
    libmesh_error_msg_if(!infinite_elem || !infinite_elem->infinite(),
                         "Error setting Elem pointer.");

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
      LIBMESH_ASSERT_FP_EQUAL(t.val, phi[t.i][t.qp], TOLERANCE);
    for (const auto & t : grad_table)
      LIBMESH_ASSERT_FP_EQUAL(0., (dphi[t.i][t.qp] - t.grad).norm_sq(), 1e-10);

    // Make sure there are actually reference values
    libmesh_error_msg_if(val_table.empty(), "No tabulated values found!");
    libmesh_error_msg_if(grad_table.empty(), "No tabulated gradients found!");

#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
  }

  void setUp()
  {
    ////////////////////////////////////////////////////////////////////////////////
    // Arbitrarily selected LEGENDRE reference values.
    ////////////////////////////////////////////////////////////////////////////////
    tabulated_vals[std::make_pair(FIRST, LEGENDRE)] = /* LEGENDRE = 14 */
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
        {0,9,Point(-0.0429458,-0.0115073,0.)},
        {1,5,Point(0.177013,-0.177013,-0.483609)},
        {2,8,Point(0.0115073,0.0115073,0.00842393)},
        {3,7,Point(-0.177013,0.0474305,0.)},
        {4,6,Point(-0.031305,-0.116832,0.10025)},
        {5,1,Point(0.0474191,-0.0474191,0.872917)},
        {6,11,Point(0.0575466,0.0575466,-0.11251)},
        {7,15,Point(-0.00353805,0.000948017,0.000111572)},
      };

    tabulated_grads[std::make_pair(SECOND, LEGENDRE)] =
      {
        {0,3,Point(-0.0959817, -0.0959817, 0.0702635)},
        {1,6,Point(0.0625228, -0.0625228, 0.0457699)},
        {2,2,Point(0.358209, 0.0959817, 0.)},
        {3,7,Point(-0.233338, 0.0625228, 0.)},
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
    tabulated_vals[std::make_pair(FIRST, JACOBI_20_00)] = /* JACOBI_20_00 = 12 */
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
        {2,2,Point(0.358209, 0.0959817, 0.)},
        {3,7,Point(-0.233338, 0.0625228, 0.)},
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
    tabulated_vals[std::make_pair(FIRST, JACOBI_30_00)] = /* JACOBI_30_00 = 13 */
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
        {2,2,Point(0.358209, 0.0959817, 0.)},
        {3,7,Point(-0.233338, 0.0625228, 0.)},
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
        {15,15,Point(-0.0423204,0.0113397,0.12514)}
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
