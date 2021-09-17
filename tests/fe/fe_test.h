#ifndef __fe_test_h__
#define __fe_test_h__

#include "test_comm.h"

#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe_base.h>
#include <libmesh/fe_interface.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/system.h>
#include <libmesh/quadrature_gauss.h>

#include <vector>

#include "libmesh_cppunit.h"

#define FETEST                                  \
  CPPUNIT_TEST( testU );                        \
  CPPUNIT_TEST( testGradU );                    \
  CPPUNIT_TEST( testGradUComp );                \
  CPPUNIT_TEST( testHessU );                    \
  CPPUNIT_TEST( testHessUComp );                \
  CPPUNIT_TEST( testDualDoesntScreamAndDie );   \
  CPPUNIT_TEST( testCustomReinit );

using namespace libMesh;

inline
Number linear_test (const Point& p,
                    const Parameters&,
                    const std::string&,
                    const std::string&)
{
  const Real & x = p(0);
  const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
  const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

  return x + 0.25*y + 0.0625*z;
}

inline
Gradient linear_test_grad (const Point&,
                           const Parameters&,
                           const std::string&,
                           const std::string&)
{
  Gradient grad = 1;
  if (LIBMESH_DIM > 1)
    grad(1) = 0.25;
  if (LIBMESH_DIM > 2)
    grad(2) = 0.0625;

  return grad;
}


inline
Number quadratic_test (const Point& p,
                       const Parameters&,
                       const std::string&,
                       const std::string&)
{
  const Real & x = p(0);
  const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
  const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

  return x*x + 0.5*y*y + 0.25*z*z + 0.125*x*y + 0.0625*x*z + 0.03125*y*z;
}

inline
Gradient quadratic_test_grad (const Point & p,
                              const Parameters&,
                              const std::string&,
                              const std::string&)
{
  const Real & x = p(0);
  const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
  const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

  Gradient grad = 2*x + 0.125*y + 0.0625*z;
  if (LIBMESH_DIM > 1)
    grad(1) = y + 0.125*x + 0.03125*z;
  if (LIBMESH_DIM > 2)
    grad(2) = 0.5*z + 0.0625*x + 0.03125*y;

  return grad;
}


inline
Number fe_cubic_test (const Point& p,
                      const Parameters&,
                      const std::string&,
                      const std::string&)
{
  const Real & x = p(0);
  const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
  const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

  return x*(1-x)*(1-x) + x*x*(1-y) + x*(1-y)*(1-z) + y*(1-y)*z + z*(1-z)*(1-z);
}

inline
Gradient fe_cubic_test_grad (const Point & p,
                             const Parameters&,
                             const std::string&,
                             const std::string&)
{
  const Real & x = p(0);
  const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
  const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

  Gradient grad = 3*x*x-4*x+1 + 2*x*(1-y) + (1-y)*(1-z);
  if (LIBMESH_DIM > 1)
    grad(1) = -x*x - x*(1-z) + (1-2*y)*z;
  if (LIBMESH_DIM > 2)
    grad(2) = -x*(1-y) + y*(1-y) + 3*z*z-4*z+1;

  return grad;
}


// Higher order rational bases need uniform weights to exactly
// represent linears; we can easily try out other functions on
// tensor product elements.
static const Real rational_w = 0.75;

inline
Number rational_test (const Point& p,
                      const Parameters&,
                      const std::string&,
                      const std::string&)
{
  const Real & x = p(0);
  const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
  const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

  const Real denom = ((1-x)*(1-x)+x*x+2*rational_w*x*(1-x))*
                     ((1-y)*(1-y)+y*y+2*rational_w*y*(1-y))*
                     ((1-z)*(1-z)+z*z+2*rational_w*z*(1-z));

  return (x + 0.25*y + 0.0625*z)/denom;
}

inline
Gradient rational_test_grad (const Point& p,
                             const Parameters&,
                             const std::string&,
                             const std::string&)
{
  const Real & x = p(0);
  const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
  const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

  const Real xpoly = (1-x)*(1-x)+x*x+2*rational_w*x*(1-x);
  const Real xderiv = -2*(1-x)+2*x+2*rational_w*(1-2*x);
  const Real ypoly = (1-y)*(1-y)+y*y+2*rational_w*y*(1-y);
  const Real yderiv = -2*(1-y)+2*y+2*rational_w*(1-2*y);
  const Real zpoly = (1-z)*(1-z)+z*z+2*rational_w*z*(1-z);
  const Real zderiv = -2*(1-z)+2*z+2*rational_w*(1-2*z);

  const Real denom = xpoly * ypoly * zpoly;

  const Real numer = (x + 0.25*y + 0.0625*z);

  Gradient grad_n = 1, grad_d = xderiv * ypoly * zpoly;
  if (LIBMESH_DIM > 1)
    {
      grad_n(1) = 0.25;
      grad_d(1) = xpoly * yderiv * zpoly;
    }
  if (LIBMESH_DIM > 2)
    {
      grad_n(2) = 0.0625;
      grad_d(2) = xpoly * ypoly * zderiv;
    }

  Gradient grad = (grad_n - numer * grad_d / denom) / denom;

  return grad;
}


#define FE_CAN_TEST_CUBIC \
  (((family != LAGRANGE && family != L2_LAGRANGE) || elem_type != TRI7) && order > 2)



template <Order order, FEFamily family, ElemType elem_type, unsigned int build_nx>
class FETestBase : public CppUnit::TestCase {

protected:
  unsigned int _dim, _nx, _ny, _nz;
  Elem *_elem;
  std::vector<dof_id_type> _dof_indices;
  FEBase * _fe;
  Mesh * _mesh;
  System * _sys;
  EquationSystems * _es;

  Real _value_tol, _grad_tol, _hess_tol;

  QGauss * _qrule;

public:
  void setUp()
  {
    _mesh = new Mesh(*TestCommWorld);
    const std::unique_ptr<Elem> test_elem = Elem::build(elem_type);
    _dim = test_elem->dim();
    const unsigned int build_ny = (_dim > 1) * build_nx;
    const unsigned int build_nz = (_dim > 2) * build_nx;

    unsigned char weight_index = 0;

    if (family == RATIONAL_BERNSTEIN)
      {
        // Make sure we can handle non-zero weight indices
        _mesh->add_node_integer("buffer integer");

        // Default weight to 1.0 so we don't get NaN from contains_point
        // checks with a default GhostPointNeighbors ghosting
        const Real default_weight = 1.0;
        weight_index = cast_int<unsigned char>
          (_mesh->add_node_datum<Real>("rational_weight", true,
                                       &default_weight));
        libmesh_assert_not_equal_to(weight_index, 0);

        // Set mapping data but not type, since here we're testing
        // rational basis functions in the FE space but testing then
        // with Lagrange bases for the mapping space.
        _mesh->set_default_mapping_data(weight_index);
      }

    MeshTools::Generation::build_cube (*_mesh,
                                       build_nx, build_ny, build_nz,
                                       0., 1., 0., 1., 0., 1.,
                                       elem_type);

    // Permute our elements randomly and rotate our mesh so we test
    // all sorts of orientations ... except with Hermite elements,
    // which are only designed to support meshes with a single
    // orientation shared by all elements.  We're also not rotating
    // the rational elements, since our test solution was designed for
    // a specific weighted mesh.
    if (family != HERMITE &&
        family != RATIONAL_BERNSTEIN)
      {
        MeshTools::Modification::permute_elements(*_mesh);

        // Not yet testing manifolds embedded in higher-D space
        if (_dim > 1)
          MeshTools::Modification::rotate(*_mesh, 4,
                                          8*(_dim>2), 16*(_dim>2));
      }

    // Set rational weights so we can exactly match our test solution
    if (family == RATIONAL_BERNSTEIN)
      {
        for (auto elem : _mesh->active_element_ptr_range())
          {
            const unsigned int nv = elem->n_vertices();
            const unsigned int nn = elem->n_nodes();
            // We want interiors in lower dimensional elements treated
            // like edges or faces as appropriate.
            const unsigned int n_edges =
              (elem->type() == EDGE3) ? 1 : elem->n_edges();
            const unsigned int n_faces =
              (elem->type() == QUAD9) ? 1 : elem->n_faces();
            const unsigned int nve = std::min(nv + n_edges, nn);
            const unsigned int nvef = std::min(nve + n_faces, nn);

            for (unsigned int i = 0; i != nv; ++i)
              elem->node_ref(i).set_extra_datum<Real>(weight_index, 1.);
            for (unsigned int i = nv; i != nve; ++i)
              elem->node_ref(i).set_extra_datum<Real>(weight_index, rational_w);
            const Real w2 = rational_w * rational_w;
            for (unsigned int i = nve; i != nvef; ++i)
              elem->node_ref(i).set_extra_datum<Real>(weight_index, w2);
            const Real w3 = rational_w * w2;
            for (unsigned int i = nvef; i != nn; ++i)
              elem->node_ref(i).set_extra_datum<Real>(weight_index, w3);
          }
      }

    _es = new EquationSystems(*_mesh);
    _sys = &(_es->add_system<System> ("SimpleSystem"));
    _sys->add_variable("u", order, family);
    _es->init();

    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    if (family == RATIONAL_BERNSTEIN && order > 1)
      {
        _sys->project_solution(rational_test, rational_test_grad, _es->parameters);
      }
    // Lagrange "cubic" on Tet7 only supports a bubble function, not
    // all of P^3
    else if (FE_CAN_TEST_CUBIC)
      {
        _sys->project_solution(fe_cubic_test, fe_cubic_test_grad, _es->parameters);
      }
    else if (order > 1)
      {
        _sys->project_solution(quadratic_test, quadratic_test_grad, _es->parameters);
      }
    else
      {
        _sys->project_solution(linear_test, linear_test_grad, _es->parameters);
      }

    FEType fe_type = _sys->variable_type(0);
    _fe = FEBase::build(_dim, fe_type).release();

    // Create quadrature rule for use in computing dual shape coefficients
    _qrule = new QGauss(_dim, fe_type.default_quadrature_order());
    _fe->attach_quadrature_rule(_qrule);

    auto rng = _mesh->active_local_element_ptr_range();
    this->_elem = rng.begin() == rng.end() ? nullptr : *(rng.begin());

    _sys->get_dof_map().dof_indices(this->_elem, _dof_indices);

    _nx = 10;
    _ny = (_dim > 1) ? _nx : 1;
    _nz = (_dim > 2) ? _nx : 1;

    this->_value_tol = TOLERANCE * sqrt(TOLERANCE);

    // We see 6.5*tol*sqrt(tol) errors on cubic Hermites with the fe_cubic
    // hermite test function
    // On Tri7 we see 10*tol*sqrt(tol) errors, even!
    this->_grad_tol = 12 * TOLERANCE * sqrt(TOLERANCE);

    this->_hess_tol = sqrt(TOLERANCE); // FIXME: we see some ~1e-5 errors?!?

    // Prerequest everything we'll want to calculate later.
    _fe->get_phi();
    _fe->get_dphi();
    _fe->get_dphidx();
#if LIBMESH_DIM > 1
    _fe->get_dphidy();
#endif
#if LIBMESH_DIM > 2
    _fe->get_dphidz();
#endif

#if LIBMESH_ENABLE_SECOND_DERIVATIVES

    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

    _fe->get_d2phi();
    _fe->get_d2phidx2();
#if LIBMESH_DIM > 1
    _fe->get_d2phidxdy();
    _fe->get_d2phidy2();
#endif
#if LIBMESH_DIM > 2
    _fe->get_d2phidxdz();
    _fe->get_d2phidydz();
    _fe->get_d2phidz2();
#endif

#endif
  }

  void tearDown()
  {
    delete _fe;
    delete _es;
    delete _mesh;
    delete _qrule;
  }
};



template <Order order, FEFamily family, ElemType elem_type>
class FETest : public FETestBase<order, family, elem_type, 1> {

public:

  template <typename Functor>
  void testLoop(Functor f)
  {
    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Handle the "more processors than elements" case
    if (!this->_elem)
      return;

    // These tests require exceptions to be enabled because a
    // TypeTensor::solve() call down in Elem::contains_point()
    // actually throws a non-fatal exception for a certain Point which
    // is not in the Elem. When exceptions are not enabled, this test
    // simply aborts.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
    for (unsigned int i=0; i != this->_nx; ++i)
      for (unsigned int j=0; j != this->_ny; ++j)
        for (unsigned int k=0; k != this->_nz; ++k)
          {
            Real x = (Real(i)/this->_nx), y = 0, z = 0;
            Point p = x;
            if (j > 0)
              p(1) = y = (Real(j)/this->_ny);
            if (k > 0)
              p(2) = z = (Real(k)/this->_nz);
            if (!this->_elem->contains_point(p))
              continue;

            // If at a singular node, cannot use FEMap::map
            if (this->_elem->local_singular_node(p) != invalid_uint)
              continue;

            std::vector<Point> master_points
              (1, FEMap::inverse_map(this->_dim, this->_elem, p));

            // Reinit at point to test against analytic solution
            this->_fe->reinit(this->_elem, &master_points);

            f(p, x, y, z);
          }
#endif // LIBMESH_ENABLE_EXCEPTIONS
  }


  void testU()
  {
    auto f = [this](Point p, Real x, Real y, Real z)
      {
        Parameters dummy;

        Number u = 0;
        for (std::size_t d = 0; d != this->_dof_indices.size(); ++d)
          u += this->_fe->get_phi()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);

        if (family == RATIONAL_BERNSTEIN && order > 1)
          LIBMESH_ASSERT_FP_EQUAL
            (libmesh_real(rational_test(p, dummy, "", "")),
             libmesh_real(u),
             this->_value_tol);
        else if (FE_CAN_TEST_CUBIC)
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
      };

    testLoop(f);
  }

  void testDualDoesntScreamAndDie()
  {
    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Handle the "more processors than elements" case
    if (!this->_elem)
      return;

    // Request dual calculations
    this->_fe->get_dual_phi();

    // reinit using the default quadrature rule in order to calculate the dual coefficients
    this->_fe->reinit(this->_elem);
  }


  void testGradU()
  {
    auto f = [this](Point p, Real x, Real y, Real z)
      {
        Parameters dummy;

        Gradient grad_u = 0;
        for (std::size_t d = 0; d != this->_dof_indices.size(); ++d)
          grad_u += this->_fe->get_dphi()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);

        if (family == RATIONAL_BERNSTEIN && order > 1)
          {
            const Gradient rat_grad =
              rational_test_grad(p, dummy, "", "");

            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(0)),
                                    libmesh_real(rat_grad(0)),
                                    this->_grad_tol);
            if (this->_dim > 1)
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(1)),
                                      libmesh_real(rat_grad(1)),
                                      this->_grad_tol);
            if (this->_dim > 2)
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(2)),
                                      libmesh_real(rat_grad(2)),
                                      this->_grad_tol);
          }
        else if (FE_CAN_TEST_CUBIC)
          {
            const Gradient cub_grad =
              fe_cubic_test_grad(p, dummy, "", "");

            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(0)),
                                    libmesh_real(cub_grad(0)),
                                    this->_grad_tol);
            if (this->_dim > 1)
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(1)),
                                      libmesh_real(cub_grad(1)),
                                      this->_grad_tol);
            if (this->_dim > 2)
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(2)),
                                      libmesh_real(cub_grad(2)),
                                      this->_grad_tol);
          }
        else if (order > 1)
          {
            const Real & x = p(0);
            const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
            const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

            LIBMESH_ASSERT_FP_EQUAL(2*x+0.125*y+0.0625*z, libmesh_real(grad_u(0)),
                                    this->_grad_tol);
            if (this->_dim > 1)
              LIBMESH_ASSERT_FP_EQUAL(y+0.125*x+0.03125*z, libmesh_real(grad_u(1)),
                                      this->_grad_tol);
            if (this->_dim > 2)
              LIBMESH_ASSERT_FP_EQUAL(0.5*z+0.0625*x+0.03125*y, libmesh_real(grad_u(2)),
                                      this->_grad_tol);
          }
        else
          {
            LIBMESH_ASSERT_FP_EQUAL(1.0, libmesh_real(grad_u(0)),
                                    this->_grad_tol);
            if (this->_dim > 1)
              LIBMESH_ASSERT_FP_EQUAL(0.25, libmesh_real(grad_u(1)),
                                      this->_grad_tol);
            if (this->_dim > 2)
              LIBMESH_ASSERT_FP_EQUAL(0.0625, libmesh_real(grad_u(2)),
                                      this->_grad_tol);
          }
      };

    testLoop(f);
  }

  void testGradUComp()
  {
    auto f = [this](Point p, Real x, Real y, Real z)
      {
        Parameters dummy;

        Number grad_u_x = 0, grad_u_y = 0, grad_u_z = 0;
        for (std::size_t d = 0; d != this->_dof_indices.size(); ++d)
          {
            grad_u_x += this->_fe->get_dphidx()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#if LIBMESH_DIM > 1
            grad_u_y += this->_fe->get_dphidy()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#endif
#if LIBMESH_DIM > 2
            grad_u_z += this->_fe->get_dphidz()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#endif
          }

        if (family == RATIONAL_BERNSTEIN && order > 1)
          {
            const Gradient rat_grad =
              rational_test_grad(p, dummy, "", "");

            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_x),
                                    libmesh_real(rat_grad(0)),
                                    this->_grad_tol);
            if (this->_dim > 1)
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_y),
                                      libmesh_real(rat_grad(1)),
                                      this->_grad_tol);
            if (this->_dim > 2)
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_z),
                                      libmesh_real(rat_grad(2)),
                                      this->_grad_tol);
          }
        else if (FE_CAN_TEST_CUBIC)
          {
            const Gradient cub_grad =
              fe_cubic_test_grad(p, dummy, "", "");

            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_x),
                                    libmesh_real(cub_grad(0)),
                                    this->_grad_tol);
            if (this->_dim > 1)
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_y),
                                      libmesh_real(cub_grad(1)),
                                      this->_grad_tol);
            if (this->_dim > 2)
              LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_z),
                                      libmesh_real(cub_grad(2)),
                                      this->_grad_tol);
          }
        else if (order > 1)
          {
            const Real & x = p(0);
            const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
            const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

            LIBMESH_ASSERT_FP_EQUAL(2*x+0.125*y+0.0625*z, libmesh_real(grad_u_x),
                                    this->_grad_tol);
            if (this->_dim > 1)
              LIBMESH_ASSERT_FP_EQUAL(y+0.125*x+0.03125*z, libmesh_real(grad_u_y),
                                      this->_grad_tol);
            if (this->_dim > 2)
              LIBMESH_ASSERT_FP_EQUAL(0.5*z+0.0625*x+0.03125*y, libmesh_real(grad_u_z),
                                      this->_grad_tol);
          }
        else
          {
            LIBMESH_ASSERT_FP_EQUAL(1.0, libmesh_real(grad_u_x),
                                    this->_grad_tol);
            if (this->_dim > 1)
              LIBMESH_ASSERT_FP_EQUAL(0.25, libmesh_real(grad_u_y),
                                      this->_grad_tol);
            if (this->_dim > 2)
              LIBMESH_ASSERT_FP_EQUAL(0.0625, libmesh_real(grad_u_z),
                                      this->_grad_tol);
          }
      };

    testLoop(f);
  }


  void testHessU()
  {
    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    auto f = [this](Point p, Real x, Real y, Real z)
      {
        Parameters dummy;

        Tensor hess_u;
        for (std::size_t d = 0; d != this->_dof_indices.size(); ++d)
          hess_u += this->_fe->get_d2phi()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);

        if (family == RATIONAL_BERNSTEIN && order > 1)
          {
            // TODO: Yeah we'll test the ugly expressions later.
          }
        else if (FE_CAN_TEST_CUBIC)
          {
            const Real & x = p(0);
            const Real & y = LIBMESH_DIM > 1 ? p(1) : 0;
            const Real & z = LIBMESH_DIM > 2 ? p(2) : 0;
            const RealTensor full_hess { 6*x-4+2*(1-y), -2*x+z-1,     y-1,
                                              -2*x+z-1,     -2*z, x+1-2*y,
                                                   y-1,  x+1-2*y, 6*z-4 };

            LIBMESH_ASSERT_FP_EQUAL(full_hess(0,0), libmesh_real(hess_u(0,0)),
                                    this->_hess_tol);
            if (this->_dim > 1)
              {
                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,1)), libmesh_real(hess_u(1,0)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(full_hess(0,1), libmesh_real(hess_u(0,1)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(full_hess(1,1), libmesh_real(hess_u(1,1)),
                                        this->_hess_tol);
              }
            if (this->_dim > 2)
              {
                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,2)), libmesh_real(hess_u(2,0)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(1,2)), libmesh_real(hess_u(2,1)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(full_hess(0,2), libmesh_real(hess_u(0,2)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(full_hess(1,2), libmesh_real(hess_u(1,2)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(full_hess(2,2), libmesh_real(hess_u(2,2)),
                                        this->_hess_tol);
              }

          }
        else if (order > 1)
          {
            LIBMESH_ASSERT_FP_EQUAL(2, libmesh_real(hess_u(0,0)),
                                    this->_hess_tol);
            if (this->_dim > 1)
              {
                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,1)), libmesh_real(hess_u(1,0)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0.125, libmesh_real(hess_u(0,1)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(1, libmesh_real(hess_u(1,1)),
                                        this->_hess_tol);
              }
            if (this->_dim > 2)
              {
                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,2)), libmesh_real(hess_u(2,0)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(1,2)), libmesh_real(hess_u(2,1)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0.0625, libmesh_real(hess_u(0,2)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0.03125, libmesh_real(hess_u(1,2)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0.5, libmesh_real(hess_u(2,2)),
                                        this->_hess_tol);
              }
          }
        else
          {
            LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(0,0)),
                                    this->_hess_tol);
            if (this->_dim > 1)
              {
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(0,1)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(1,0)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(1,1)),
                                        this->_hess_tol);
              }
            if (this->_dim > 2)
              {
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(0,2)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(1,2)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(2,0)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(2,1)),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(2,2)),
                                        this->_hess_tol);
              }
          }
      };

    testLoop(f);
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
  }

  void testHessUComp()
  {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

    auto f = [this](Point p, Real x, Real y, Real z)
      {
        Parameters dummy;

        Number hess_u_xx = 0, hess_u_xy = 0, hess_u_yy = 0,
               hess_u_xz = 0, hess_u_yz = 0, hess_u_zz = 0;
        for (std::size_t d = 0; d != this->_dof_indices.size(); ++d)
          {
            hess_u_xx += this->_fe->get_d2phidx2()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#if LIBMESH_DIM > 1
            hess_u_xy += this->_fe->get_d2phidxdy()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
            hess_u_yy += this->_fe->get_d2phidy2()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#endif
#if LIBMESH_DIM > 2
            hess_u_xz += this->_fe->get_d2phidxdz()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
            hess_u_yz += this->_fe->get_d2phidydz()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
            hess_u_zz += this->_fe->get_d2phidz2()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);
#endif
          }

        if (family == RATIONAL_BERNSTEIN && order > 1)
          {
            // TODO: tedious calculus
          }
        else if (FE_CAN_TEST_CUBIC)
          {
            const Real & x = p(0);
            const Real & y = LIBMESH_DIM > 1 ? p(1) : 0;
            const Real & z = LIBMESH_DIM > 2 ? p(2) : 0;
            const RealTensor full_hess { 6*x-4+2*(1-y), -2*x+z-1,     y-1,
                                              -2*x+z-1,     -2*z, x+1-2*y,
                                                   y-1,  x+1-2*y, 6*z-4 };

            LIBMESH_ASSERT_FP_EQUAL(full_hess(0,0), libmesh_real(hess_u_xx),
                                    this->_hess_tol);
            if (this->_dim > 1)
              {
                LIBMESH_ASSERT_FP_EQUAL(full_hess(0,1), libmesh_real(hess_u_xy),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(full_hess(1,1), libmesh_real(hess_u_yy),
                                        this->_hess_tol);
              }
            if (this->_dim > 2)
              {
                LIBMESH_ASSERT_FP_EQUAL(full_hess(0,2), libmesh_real(hess_u_xz),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(full_hess(1,2), libmesh_real(hess_u_yz),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(full_hess(2,2), libmesh_real(hess_u_zz),
                                        this->_hess_tol);
              }

          }
        else if (order > 1)
          {
            LIBMESH_ASSERT_FP_EQUAL(2, libmesh_real(hess_u_xx),
                                    this->_hess_tol);
            if (this->_dim > 1)
              {
                LIBMESH_ASSERT_FP_EQUAL(0.125, libmesh_real(hess_u_xy),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(1, libmesh_real(hess_u_yy),
                                        this->_hess_tol);
              }
            if (this->_dim > 2)
              {
                LIBMESH_ASSERT_FP_EQUAL( 0.0625, libmesh_real(hess_u_xz),
                                         this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL( 0.03125, libmesh_real(hess_u_yz),
                                         this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL( 0.5, libmesh_real(hess_u_zz),
                                         this->_hess_tol);
              }
          }
        else
          {
            LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_xx),
                                    this->_hess_tol);
            if (this->_dim > 1)
              {
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_xy),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_yy),
                                        this->_hess_tol);
              }
            if (this->_dim > 2)
              {
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_xz),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_yz),
                                        this->_hess_tol);
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_zz),
                                        this->_hess_tol);
              }
          }
      };

    testLoop(f);
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
  }

  void testCustomReinit()
  {
    std::vector<Point> q_points;
    std::vector<Real> weights;
    q_points.resize(3); weights.resize(3);
    q_points[0](0) = 0.0; q_points[0](1) = 0.0; weights[0] = Real(1)/6;
    q_points[1](0) = 1.0; q_points[1](1) = 0.0; weights[1] = weights[0];
    q_points[2](0) = 0.0; q_points[2](1) = 1.0; weights[2] = weights[0];

    FEType fe_type = this->_sys->variable_type(0);
    std::unique_ptr<FEBase> fe (FEBase::build(this->_dim, fe_type));
    const int extraorder = 3;
    std::unique_ptr<QBase> qrule (fe_type.default_quadrature_rule (this->_dim, extraorder));
    fe->attach_quadrature_rule (qrule.get());

    const std::vector<Point> & q_pos = fe->get_xyz();

    for (const auto & elem : this->_mesh->active_local_element_ptr_range()) {
      fe->reinit (elem, &q_points, &weights);
      CPPUNIT_ASSERT_EQUAL(q_points.size(), std::size_t(3));
      CPPUNIT_ASSERT_EQUAL(q_pos.size(), std::size_t(3));  // 6? bug?
    }
  }

};


#define INSTANTIATE_FETEST(order, family, elemtype)                     \
  class FETest_##order##_##family##_##elemtype : public FETest<order, family, elemtype> { \
  public:                                                               \
  CPPUNIT_TEST_SUITE( FETest_##order##_##family##_##elemtype );         \
  FETEST                                                                \
  CPPUNIT_TEST_SUITE_END();                                             \
  };                                                                    \
                                                                        \
  CPPUNIT_TEST_SUITE_REGISTRATION( FETest_##order##_##family##_##elemtype );

#endif // #ifdef __fe_test_h__
