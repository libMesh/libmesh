#include "fe_test.h" // For linear_test(), headers, FETestBase, etc

// Helper functions for restricting normal components
namespace {

template <typename MyTensor, typename MyVector>
MyTensor tangential_hessian (const MyTensor hess,
                             const MyVector normal)
{
  auto full_hess_n = hess*normal;
  MyTensor tang_hess = hess;

  const MyTensor n_n = outer_product(normal, normal);
  for (unsigned int k=0; k != 3; ++k)
    {
      MyVector ek {0};
      ek(k) = 1;
      const auto Hkn = ek * full_hess_n;
      const auto Hkn_ek = Hkn*ek;
      tang_hess -= outer_product(Hkn_ek, normal);
      tang_hess -= outer_product(normal, Hkn_ek);
      tang_hess += 2.*(Hkn_ek*normal) * n_n;
    }
  tang_hess -= normal*(full_hess_n) * n_n;

  return tang_hess;
}

} // anonymous namespace

#define SIDEFETEST                              \
  CPPUNIT_TEST( testU );                        \
  CPPUNIT_TEST( testGradU );                    \
  CPPUNIT_TEST( testGradUComp );                \
  CPPUNIT_TEST( testHessU );                    \
  CPPUNIT_TEST( testHessUComp );

const unsigned int N_x=2;

template <Order order, FEFamily family, ElemType elem_type>
class FESideTest : public FETestBase<order, family, elem_type, N_x> {

private:
  std::unique_ptr<FEBase> _fe_side;
  std::unique_ptr<QGauss> _qrule_side;

public:
  void setUp()
  {
    FETestBase<order,family,elem_type,N_x>::setUp();

    FEType fe_type = this->_sys->variable_type(0);
    _fe_side = FEBase::build(this->_dim, fe_type);

    _qrule_side = std::make_unique<QGauss>(this->_dim-1, fe_type.default_quadrature_order());
    _fe_side->attach_quadrature_rule(this->_qrule_side.get());

    _fe_side->get_xyz();
    _fe_side->get_normals();
    _fe_side->get_phi();
    _fe_side->get_dphi();
    _fe_side->get_dphidx();
#if LIBMESH_DIM > 1
    _fe_side->get_dphidy();
#endif
#if LIBMESH_DIM > 2
    _fe_side->get_dphidz();
#endif

#if LIBMESH_ENABLE_SECOND_DERIVATIVES

    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

    _fe_side->get_d2phi();
    _fe_side->get_d2phidx2();
#if LIBMESH_DIM > 1
    _fe_side->get_d2phidxdy();
    _fe_side->get_d2phidy2();
#endif
#if LIBMESH_DIM > 2
    _fe_side->get_d2phidxdz();
    _fe_side->get_d2phidydz();
    _fe_side->get_d2phidz2();
#endif

#endif
  }

  void tearDown()
  {
    FETestBase<order,family,elem_type,N_x>::tearDown();
  }

  template <typename Functor>
  void testSideLoop(Functor f)
  {
    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    for (auto e : this->_mesh->active_local_element_ptr_range())
      {
        this->_elem = e;

        this->_sys->get_dof_map().dof_indices(this->_elem, this->_dof_indices);

        for (unsigned int s=0; s != this->_elem->n_sides(); ++s)
          {
            _fe_side->reinit(this->_elem, s);

            f();
          }
      }
  }

  void testU()
  {
    LOG_UNIT_TEST;

    auto f = [this]() {
      Parameters dummy;

      auto phi = this->_fe_side->get_phi();
      CPPUNIT_ASSERT(phi.size());
      CPPUNIT_ASSERT_EQUAL(phi.size(), this->_dof_indices.size());

      std::size_t n_qp = phi[0].size();
      CPPUNIT_ASSERT(n_qp > 0);

      auto xyz = this->_fe_side->get_xyz();
      CPPUNIT_ASSERT_EQUAL(n_qp, xyz.size());

      for (auto qp : index_range(phi[0]))
        {
          Number u = 0;
          for (auto d : index_range(phi))
            u += phi[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);

          const Point p = xyz[qp];
          Real x = p(0),
               y = LIBMESH_DIM > 1 ? p(1) : 0,
               z = LIBMESH_DIM > 2 ? p(2) : 0;

          if (family == RATIONAL_BERNSTEIN && order > 1)
            LIBMESH_ASSERT_FP_EQUAL
              (libmesh_real(rational_test(p, dummy, "", "")),
               libmesh_real(u),
               this->_value_tol);
          else if (order > 3)
            LIBMESH_ASSERT_FP_EQUAL
              (libmesh_real(fe_quartic_test(p, dummy, "", "")),
               libmesh_real(u),
               this->_value_tol);
          else if (order > 2)
            LIBMESH_ASSERT_FP_EQUAL
              (libmesh_real(fe_cubic_test(p, dummy, "", "")),
               libmesh_real(u),
               this->_value_tol);
          else if (order > 1)
            LIBMESH_ASSERT_FP_EQUAL
              (libmesh_real(x*x + 0.5*y*y + 0.25*z*z + 0.125*x*y +
                            0.0625*x*z + 0.03125*y*z),
               libmesh_real(u),
               this->_value_tol);
          else
            LIBMESH_ASSERT_FP_EQUAL
              (libmesh_real(x + 0.25*y + 0.0625*z),
               libmesh_real(u),
               this->_value_tol);
        }
    };

    testSideLoop(f);
  }

  void testGradU()
  {
    LOG_UNIT_TEST;

    auto f = [this]() {
      Parameters dummy;

      auto dphi = this->_fe_side->get_dphi();
      CPPUNIT_ASSERT(dphi.size() > 0);
      CPPUNIT_ASSERT_EQUAL(dphi.size(), this->_dof_indices.size());

      std::size_t n_qp = dphi[0].size();
      CPPUNIT_ASSERT(n_qp > 0);

      auto xyz = this->_fe_side->get_xyz();
      CPPUNIT_ASSERT_EQUAL(n_qp, xyz.size());

      auto normals = this->_fe_side->get_normals();
      CPPUNIT_ASSERT_EQUAL(n_qp, normals.size());

      for (auto qp : index_range(dphi[0]))
        {
          Gradient grad_u = 0;
          for (auto d : index_range(dphi))
            grad_u += dphi[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);

          const Point p = xyz[qp];

          // We can only trust the tangential component of gradient to
          // match!  The gradient in the normal direction should be 0
          // for a shape-preserving master->physical element
          // transformation, and could be just about anything for a
          // skew transformation (as occurs when we build meshes with
          // triangles, tets, prisms, pyramids...)
          Point n = normals[qp];

          const Gradient tangential_grad_u = grad_u - ((grad_u * n) * n);

          if (family == RATIONAL_BERNSTEIN && order > 1)
            {
              // TODO: Yeah we'll test the ugly expressions later.
              return;
            }

          RealGradient true_grad = this->true_gradient(p);
          if (order == 0)
            true_grad.zero();

          const Gradient true_tangential_grad = true_grad - ((true_grad * n) * n);

          LIBMESH_ASSERT_FP_EQUAL(libmesh_real(tangential_grad_u(0)),
                                  libmesh_real(true_tangential_grad(0)),
                                  this->_grad_tol);
          if (this->_dim > 1)
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(tangential_grad_u(1)),
                                    libmesh_real(true_tangential_grad(1)),
                                    this->_grad_tol);
          if (this->_dim > 2)
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(tangential_grad_u(2)),
                                    libmesh_real(true_tangential_grad(2)),
                                    this->_grad_tol);
        }
    };

    testSideLoop(f);
  }

  void testGradUComp()
  {
    LOG_UNIT_TEST;

    auto f = [this]() {
      Parameters dummy;

      auto dphidx = this->_fe_side->get_dphidx();
      CPPUNIT_ASSERT(dphidx.size() > 0);
      CPPUNIT_ASSERT_EQUAL(dphidx.size(), this->_dof_indices.size());
#if LIBMESH_DIM > 1
      auto dphidy = this->_fe_side->get_dphidy();
      CPPUNIT_ASSERT_EQUAL(dphidy.size(), this->_dof_indices.size());
#endif
#if LIBMESH_DIM > 2
      auto dphidz = this->_fe_side->get_dphidz();
      CPPUNIT_ASSERT_EQUAL(dphidz.size(), this->_dof_indices.size());
#endif

      std::size_t n_qp = dphidx[0].size();
      CPPUNIT_ASSERT(n_qp > 0);

      auto xyz = this->_fe_side->get_xyz();
      CPPUNIT_ASSERT_EQUAL(n_qp, xyz.size());

      auto normals = this->_fe_side->get_normals();
      CPPUNIT_ASSERT_EQUAL(n_qp, normals.size());

      for (auto qp : index_range(dphidx[0]))
        {
          Number grad_u_x = 0, grad_u_y = 0, grad_u_z = 0;
          for (auto d : index_range(dphidx))
            {
              grad_u_x += dphidx[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#if LIBMESH_DIM > 1
              grad_u_y += dphidy[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#endif
#if LIBMESH_DIM > 2
              grad_u_z += dphidz[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#endif
            }

          const Point p = xyz[qp];

          // We can only trust the tangential component of gradient to
          // match!  The gradient in the normal direction should be 0
          // for a shape-preserving master->physical element
          // transformation, and could be just about anything for a
          // skew transformation (as occurs when we build meshes with
          // triangles, tets, prisms, pyramids...)
          Point n = normals[qp];

          const Gradient grad_u {grad_u_x, grad_u_y, grad_u_z};
          const Gradient tangential_grad_u = grad_u - ((grad_u * n) * n);

          if (family == RATIONAL_BERNSTEIN && order > 1)
            {
              // TODO: Yeah we'll test the ugly expressions later.
              return;
            }

          RealGradient true_grad = this->true_gradient(p);
          if (order == 0)
            true_grad.zero();

          const Gradient true_tangential_grad = true_grad - ((true_grad * n) * n);

          LIBMESH_ASSERT_FP_EQUAL(libmesh_real(tangential_grad_u(0)),
                                  libmesh_real(true_tangential_grad(0)),
                                  this->_grad_tol);
          if (this->_dim > 1)
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(tangential_grad_u(1)),
                                    libmesh_real(true_tangential_grad(1)),
                                    this->_grad_tol);
          if (this->_dim > 2)
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(tangential_grad_u(2)),
                                    libmesh_real(true_tangential_grad(2)),
                                    this->_grad_tol);
        }
    };

    testSideLoop(f);
  }


  void testHessU()
  {
    LOG_UNIT_TEST;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

    auto f = [this]() {
      Parameters dummy;

      auto d2phi = this->_fe_side->get_d2phi();
      CPPUNIT_ASSERT(d2phi.size() > 0);
      CPPUNIT_ASSERT_EQUAL(d2phi.size(), this->_dof_indices.size());

      std::size_t n_qp = d2phi[0].size();
      CPPUNIT_ASSERT(n_qp > 0);

      auto xyz = this->_fe_side->get_xyz();
      CPPUNIT_ASSERT_EQUAL(n_qp, xyz.size());

      auto normals = this->_fe_side->get_normals();
      CPPUNIT_ASSERT_EQUAL(n_qp, normals.size());

      for (auto qp : index_range(d2phi[0]))
        {
          Tensor hess_u;
          for (auto d : index_range(d2phi))
            hess_u += d2phi[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);

          // We can only trust the tangential-tangential components of
          // hessians, t1'*H*t2 for tangential vectors t1 and t2, to
          // match!  The hessian components in the normal direction
          // should be 0 for a shape-preserving master->physical
          // element transformation, but could be just about anything
          // for a skew transformation (as occurs when we build meshes
          // with triangles, tets, prisms, pyramids...)
          Point n = normals[qp];
          hess_u = tangential_hessian(hess_u, n);

          if (family == RATIONAL_BERNSTEIN && order > 1)
            {
              // TODO: Yeah we'll test the ugly expressions later.
              return;
            }

          const Point p = xyz[qp];
          RealTensor true_hess = this->true_hessian(p);

          // Subtract off values only relevant in normal directions
          RealTensor tangential_hess = tangential_hessian(true_hess, n);

          LIBMESH_ASSERT_FP_EQUAL(tangential_hess(0,0), libmesh_real(hess_u(0,0)),
                                  this->_hess_tol);
          if (this->_dim > 1)
            {
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,1)), libmesh_real(hess_u(1,0)),
                                      this->_hess_tol);
              LIBMESH_ASSERT_FP_EQUAL(tangential_hess(0,1), libmesh_real(hess_u(0,1)),
                                      this->_hess_tol);
              LIBMESH_ASSERT_FP_EQUAL(tangential_hess(1,1), libmesh_real(hess_u(1,1)),
                                      this->_hess_tol);
            }
          if (this->_dim > 2)
            {
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,2)), libmesh_real(hess_u(2,0)),
                                      this->_hess_tol);
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(1,2)), libmesh_real(hess_u(2,1)),
                                      this->_hess_tol);
              LIBMESH_ASSERT_FP_EQUAL(tangential_hess(0,2), libmesh_real(hess_u(0,2)),
                                      this->_hess_tol);
              LIBMESH_ASSERT_FP_EQUAL(tangential_hess(1,2), libmesh_real(hess_u(1,2)),
                                      this->_hess_tol);
              LIBMESH_ASSERT_FP_EQUAL(tangential_hess(2,2), libmesh_real(hess_u(2,2)),
                                      this->_hess_tol);
            }
        }
    };

    testSideLoop(f);
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
  }

  void testHessUComp()
  {
    LOG_UNIT_TEST;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

    auto f = [this]() {
      Parameters dummy;

      auto d2phidx2 = this->_fe_side->get_d2phidx2();
      CPPUNIT_ASSERT(d2phidx2.size() > 0);
      CPPUNIT_ASSERT_EQUAL(d2phidx2.size(), this->_dof_indices.size());
#if LIBMESH_DIM > 1
      auto d2phidxdy = this->_fe_side->get_d2phidxdy();
      CPPUNIT_ASSERT_EQUAL(d2phidxdy.size(), this->_dof_indices.size());
      auto d2phidy2 = this->_fe_side->get_d2phidy2();
      CPPUNIT_ASSERT_EQUAL(d2phidy2.size(), this->_dof_indices.size());
#endif
#if LIBMESH_DIM > 2
      auto d2phidxdz = this->_fe_side->get_d2phidxdz();
      CPPUNIT_ASSERT_EQUAL(d2phidxdz.size(), this->_dof_indices.size());
      auto d2phidydz = this->_fe_side->get_d2phidydz();
      CPPUNIT_ASSERT_EQUAL(d2phidydz.size(), this->_dof_indices.size());
      auto d2phidz2 = this->_fe_side->get_d2phidz2();
      CPPUNIT_ASSERT_EQUAL(d2phidz2.size(), this->_dof_indices.size());
#endif

      std::size_t n_qp = d2phidx2[0].size();
      CPPUNIT_ASSERT(n_qp > 0);

      auto xyz = this->_fe_side->get_xyz();
      CPPUNIT_ASSERT_EQUAL(n_qp, xyz.size());

      auto normals = this->_fe_side->get_normals();
      CPPUNIT_ASSERT_EQUAL(n_qp, normals.size());

      for (auto qp : index_range(d2phidx2[0]))
        {
          Tensor hess_u;
          for (auto d : index_range(d2phidx2))
            {
              hess_u(0,0) += d2phidx2[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#if LIBMESH_DIM > 1
              hess_u(0,1) += d2phidxdy[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
              hess_u(1,1) += d2phidy2[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#endif
#if LIBMESH_DIM > 2
              hess_u(0,2) += d2phidxdz[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
              hess_u(1,2) += d2phidydz[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
              hess_u(2,2) += d2phidz2[d][qp] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#endif
            }

#if LIBMESH_DIM > 1
          hess_u(1,0) = hess_u(0,1);
#endif
#if LIBMESH_DIM > 2
          hess_u(2,0) = hess_u(0,2);
          hess_u(2,1) = hess_u(1,2);
#endif

          // We can only trust the tangential-tangential components of
          // hessians, t1'*H*t2 for tangential vectors t1 and t2, to
          // match!  The hessian components in the normal direction
          // should be 0 for a shape-preserving master->physical
          // element transformation, but could be just about anything
          // for a skew transformation (as occurs when we build meshes
          // with triangles, tets, prisms, pyramids...)
          Point n = normals[qp];
          hess_u = tangential_hessian(hess_u, n);

          if (family == RATIONAL_BERNSTEIN && order > 1)
            return;

          const Point p = xyz[qp];
          RealTensor true_hess = this->true_hessian(p);

          // Subtract off values only relevant in normal directions
          RealTensor tangential_hess = tangential_hessian(true_hess, n);

          LIBMESH_ASSERT_FP_EQUAL(tangential_hess(0,0), libmesh_real(hess_u(0,0)),
                                  this->_hess_tol);
          if (this->_dim > 1)
            {
              LIBMESH_ASSERT_FP_EQUAL(tangential_hess(0,1), libmesh_real(hess_u(0,1)),
                                      this->_hess_tol);
              LIBMESH_ASSERT_FP_EQUAL(tangential_hess(1,1), libmesh_real(hess_u(1,1)),
                                      this->_hess_tol);
            }
          if (this->_dim > 2)
            {
              LIBMESH_ASSERT_FP_EQUAL( tangential_hess(0,2), libmesh_real(hess_u(0,2)),
                                       this->_hess_tol);
              LIBMESH_ASSERT_FP_EQUAL( tangential_hess(1,2), libmesh_real(hess_u(1,2)),
                                       this->_hess_tol);
              LIBMESH_ASSERT_FP_EQUAL( tangential_hess(2,2), libmesh_real(hess_u(2,2)),
                                       this->_hess_tol);
            }
        }
    };

    testSideLoop(f);
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
  }

};


#define INSTANTIATE_FESIDETEST(order, family, elemtype)                 \
  class FESideTest_##order##_##family##_##elemtype : public FESideTest<order, family, elemtype> { \
  public:                                                               \
  FESideTest_##order##_##family##_##elemtype() :                        \
    FESideTest<order,family,elemtype>() {                               \
    if (unitlog->summarized_logs_enabled())                             \
      this->libmesh_suite_name = "FESideTest";                          \
    else                                                                \
      this->libmesh_suite_name = "FESideTest_" #order "_" #family "_" #elemtype; \
  }                                                                     \
  CPPUNIT_TEST_SUITE( FESideTest_##order##_##family##_##elemtype );     \
  SIDEFETEST                                                            \
  CPPUNIT_TEST_SUITE_END();                                             \
  };                                                                    \
                                                                        \
  CPPUNIT_TEST_SUITE_REGISTRATION( FESideTest_##order##_##family##_##elemtype );


INSTANTIATE_FESIDETEST(CONSTANT, SIDE_HIERARCHIC, EDGE3);
INSTANTIATE_FESIDETEST(FIRST, SIDE_HIERARCHIC, EDGE3);
INSTANTIATE_FESIDETEST(SECOND, SIDE_HIERARCHIC, EDGE3);
INSTANTIATE_FESIDETEST(THIRD, SIDE_HIERARCHIC, EDGE3);
INSTANTIATE_FESIDETEST(FOURTH, SIDE_HIERARCHIC, EDGE3);

#if LIBMESH_DIM > 1
INSTANTIATE_FESIDETEST(CONSTANT, SIDE_HIERARCHIC, QUAD9);
INSTANTIATE_FESIDETEST(FIRST, SIDE_HIERARCHIC, QUAD9);
INSTANTIATE_FESIDETEST(SECOND, SIDE_HIERARCHIC, QUAD9);
INSTANTIATE_FESIDETEST(THIRD, SIDE_HIERARCHIC, QUAD9);
INSTANTIATE_FESIDETEST(FOURTH, SIDE_HIERARCHIC, QUAD9);

INSTANTIATE_FESIDETEST(CONSTANT, SIDE_HIERARCHIC, QUAD8);
INSTANTIATE_FESIDETEST(FIRST, SIDE_HIERARCHIC, QUAD8);
INSTANTIATE_FESIDETEST(SECOND, SIDE_HIERARCHIC, QUAD8);
INSTANTIATE_FESIDETEST(THIRD, SIDE_HIERARCHIC, QUAD8);
INSTANTIATE_FESIDETEST(FOURTH, SIDE_HIERARCHIC, QUAD8);

INSTANTIATE_FESIDETEST(CONSTANT, SIDE_HIERARCHIC, TRI6);
INSTANTIATE_FESIDETEST(FIRST, SIDE_HIERARCHIC, TRI6);
INSTANTIATE_FESIDETEST(SECOND, SIDE_HIERARCHIC, TRI6);
INSTANTIATE_FESIDETEST(THIRD, SIDE_HIERARCHIC, TRI6);
INSTANTIATE_FESIDETEST(FOURTH, SIDE_HIERARCHIC, TRI6);

INSTANTIATE_FESIDETEST(CONSTANT, SIDE_HIERARCHIC, TRI7);
INSTANTIATE_FESIDETEST(FIRST, SIDE_HIERARCHIC, TRI7);
INSTANTIATE_FESIDETEST(SECOND, SIDE_HIERARCHIC, TRI7);
INSTANTIATE_FESIDETEST(THIRD, SIDE_HIERARCHIC, TRI7);
INSTANTIATE_FESIDETEST(FOURTH, SIDE_HIERARCHIC, TRI7);
#endif

#if LIBMESH_DIM > 2
INSTANTIATE_FESIDETEST(CONSTANT, SIDE_HIERARCHIC, HEX27);
INSTANTIATE_FESIDETEST(FIRST, SIDE_HIERARCHIC, HEX27);
INSTANTIATE_FESIDETEST(SECOND, SIDE_HIERARCHIC, HEX27);
INSTANTIATE_FESIDETEST(THIRD, SIDE_HIERARCHIC, HEX27);
INSTANTIATE_FESIDETEST(FOURTH, SIDE_HIERARCHIC, HEX27);
//
INSTANTIATE_FESIDETEST(CONSTANT, SIDE_HIERARCHIC, TET14);
INSTANTIATE_FESIDETEST(FIRST, SIDE_HIERARCHIC, TET14);
INSTANTIATE_FESIDETEST(SECOND, SIDE_HIERARCHIC, TET14);
INSTANTIATE_FESIDETEST(THIRD, SIDE_HIERARCHIC, TET14);
INSTANTIATE_FESIDETEST(FOURTH, SIDE_HIERARCHIC, TET14);
//
INSTANTIATE_FESIDETEST(CONSTANT, SIDE_HIERARCHIC, PRISM20);
INSTANTIATE_FESIDETEST(FIRST, SIDE_HIERARCHIC, PRISM20);
INSTANTIATE_FESIDETEST(SECOND, SIDE_HIERARCHIC, PRISM20);
INSTANTIATE_FESIDETEST(THIRD, SIDE_HIERARCHIC, PRISM20);
INSTANTIATE_FESIDETEST(FOURTH, SIDE_HIERARCHIC, PRISM20);
//
INSTANTIATE_FESIDETEST(CONSTANT, SIDE_HIERARCHIC, PRISM21);
INSTANTIATE_FESIDETEST(FIRST, SIDE_HIERARCHIC, PRISM21);
INSTANTIATE_FESIDETEST(SECOND, SIDE_HIERARCHIC, PRISM21);
INSTANTIATE_FESIDETEST(THIRD, SIDE_HIERARCHIC, PRISM21);
INSTANTIATE_FESIDETEST(FOURTH, SIDE_HIERARCHIC, PRISM21);
#endif
