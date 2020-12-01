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
  CPPUNIT_TEST( testDualDoesntScreamAndDie );

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




template <Order order, FEFamily family, ElemType elem_type>
class FETest : public CppUnit::TestCase {

private:
  unsigned int _dim, _nx, _ny, _nz;
  Elem *_elem;
  std::vector<dof_id_type> _dof_indices;
  FEBase * _fe;
  Mesh * _mesh;
  System * _sys;
  EquationSystems * _es;

  Real value_tol, grad_tol, hess_tol;

  QGauss * _qrule;

public:
  void setUp()
  {
    _mesh = new Mesh(*TestCommWorld);
    const std::unique_ptr<Elem> test_elem = Elem::build(elem_type);
    _dim = test_elem->dim();
    const unsigned int ny = _dim > 1;
    const unsigned int nz = _dim > 2;

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
                                       1, ny, nz,
                                       0., 1., 0., ny, 0., nz,
                                       elem_type);

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
    _elem = rng.begin() == rng.end() ? nullptr : *(rng.begin());

    _sys->get_dof_map().dof_indices(_elem, _dof_indices);

    _nx = 10;
    _ny = (_dim > 1) ? _nx : 1;
    _nz = (_dim > 2) ? _nx : 1;

    // TOLERANCE * TOLERANCE doesn't work for 3D cubic Hermite?
    value_tol = TOLERANCE * sqrt(TOLERANCE);

    // TOLERANCE * sqrt(TOLERANCE) too low for 3D quartic Hierarchic
    grad_tol = 2 * TOLERANCE * sqrt(TOLERANCE);

    hess_tol = sqrt(TOLERANCE); // FIXME: we see some ~1e-5 errors?!?

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

  void testU()
  {
    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Handle the "more processors than elements" case
    if (!_elem)
      return;

    Parameters dummy;

    // These tests require exceptions to be enabled because a
    // TypeTensor::solve() call down in Elem::contains_point()
    // actually throws a non-fatal exception for a certain Point which
    // is not in the Elem. When exceptions are not enabled, this test
    // simply aborts.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
    for (unsigned int i=0; i != _nx; ++i)
      for (unsigned int j=0; j != _ny; ++j)
        for (unsigned int k=0; k != _nz; ++k)
          {
            Real x = (Real(i)/_nx), y = 0, z = 0;
            Point p = x;
            if (j > 0)
              p(1) = y = (Real(j)/_ny);
            if (k > 0)
              p(2) = z = (Real(k)/_nz);
            if (!_elem->contains_point(p))
              continue;

            // If at a singular node, cannot use FEMap::map
            if (_elem->local_singular_node(p) != invalid_uint)
              continue;

            std::vector<Point> master_points
              (1, FEMap::inverse_map(_dim, _elem, p));

            // Reinit at point to test against analytic solution
            _fe->reinit(_elem, &master_points);

            Number u = 0;
            for (std::size_t d = 0; d != _dof_indices.size(); ++d)
              u += _fe->get_phi()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);

            if (family == RATIONAL_BERNSTEIN && order > 1)
              LIBMESH_ASSERT_FP_EQUAL
                (libmesh_real(rational_test(p, dummy, "", "")),
                 libmesh_real(u),
                 value_tol);
            else if (order > 1)
              LIBMESH_ASSERT_FP_EQUAL
                (libmesh_real(x*x + 0.5*y*y + 0.25*z*z + 0.125*x*y +
                              0.0625*x*z + 0.03125*y*z),
                 libmesh_real(u),
                 value_tol);
            else
              LIBMESH_ASSERT_FP_EQUAL
                (libmesh_real(x + 0.25*y + 0.0625*z),
                 libmesh_real(u),
                 value_tol);
          }
#endif
  }

  void testDualDoesntScreamAndDie()
  {
    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Handle the "more processors than elements" case
    if (!_elem)
      return;

    // Request dual calculations
    _fe->get_dual_phi();

    // reinit using the default quadrature rule in order to calculate the dual coefficients
    _fe->reinit(_elem);
  }


  void testGradU()
  {
    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Handle the "more processors than elements" case
    if (!_elem)
      return;

    Parameters dummy;

    // These tests require exceptions to be enabled because a
    // TypeTensor::solve() call down in Elem::contains_point()
    // actually throws a non-fatal exception for a certain Point which
    // is not in the Elem. When exceptions are not enabled, this test
    // simply aborts.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
    for (unsigned int i=0; i != _nx; ++i)
      for (unsigned int j=0; j != _ny; ++j)
        for (unsigned int k=0; k != _nz; ++k)
          {
            Point p(Real(i)/_nx);
            if (j > 0)
              p(1) = Real(j)/_ny;
            if (k > 0)
              p(2) = Real(k)/_ny;
            if (!_elem->contains_point(p))
              continue;

            // If at a singular node, cannot use FEMap::map
            if (_elem->local_singular_node(p) != invalid_uint)
              continue;

            std::vector<Point> master_points
              (1, FEMap::inverse_map(_dim, _elem, p));

            _fe->reinit(_elem, &master_points);

            Gradient grad_u = 0;
            for (std::size_t d = 0; d != _dof_indices.size(); ++d)
              grad_u += _fe->get_dphi()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);

            if (family == RATIONAL_BERNSTEIN && order > 1)
              {
                const Gradient rat_grad =
                  rational_test_grad(p, dummy, "", "");

                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(0)),
                                        libmesh_real(rat_grad(0)),
                                        grad_tol);
                if (_dim > 1)
                  LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(1)),
                                          libmesh_real(rat_grad(1)),
                                          grad_tol);
                if (_dim > 2)
                  LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(2)),
                                          libmesh_real(rat_grad(2)),
                                          grad_tol);
              }
            else if (order > 1)
              {
                const Real & x = p(0);
                const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
                const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

                LIBMESH_ASSERT_FP_EQUAL(2*x+0.125*y+0.0625*z, libmesh_real(grad_u(0)),
                                        grad_tol);
                if (_dim > 1)
                  LIBMESH_ASSERT_FP_EQUAL(y+0.125*x+0.03125*z, libmesh_real(grad_u(1)),
                                          grad_tol);
                if (_dim > 2)
                  LIBMESH_ASSERT_FP_EQUAL(0.5*z+0.0625*x+0.03125*y, libmesh_real(grad_u(2)),
                                          grad_tol);
              }
            else
              {
                LIBMESH_ASSERT_FP_EQUAL(1.0, libmesh_real(grad_u(0)),
                                        grad_tol);
                if (_dim > 1)
                  LIBMESH_ASSERT_FP_EQUAL(0.25, libmesh_real(grad_u(1)),
                                          grad_tol);
                if (_dim > 2)
                  LIBMESH_ASSERT_FP_EQUAL(0.0625, libmesh_real(grad_u(2)),
                                          grad_tol);
              }
          }
#endif
  }

  void testGradUComp()
  {
    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Handle the "more processors than elements" case
    if (!_elem)
      return;

    Parameters dummy;

    // These tests require exceptions to be enabled because a
    // TypeTensor::solve() call down in Elem::contains_point()
    // actually throws a non-fatal exception for a certain Point which
    // is not in the Elem. When exceptions are not enabled, this test
    // simply aborts.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
    for (unsigned int i=0; i != _nx; ++i)
      for (unsigned int j=0; j != _ny; ++j)
        for (unsigned int k=0; k != _nz; ++k)
          {
            Point p(Real(i)/_nx);
            if (j > 0)
              p(1) = Real(j)/_ny;
            if (k > 0)
              p(2) = Real(k)/_ny;
            if (!_elem->contains_point(p))
              continue;

            // If at a singular node, cannot use FEMap::map
            if (_elem->local_singular_node(p) != invalid_uint)
              continue;

            std::vector<Point> master_points
              (1, FEMap::inverse_map(_dim, _elem, p));

            _fe->reinit(_elem, &master_points);

            Number grad_u_x = 0, grad_u_y = 0, grad_u_z = 0;
            for (std::size_t d = 0; d != _dof_indices.size(); ++d)
              {
                grad_u_x += _fe->get_dphidx()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);
#if LIBMESH_DIM > 1
                grad_u_y += _fe->get_dphidy()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);
#endif
#if LIBMESH_DIM > 2
                grad_u_z += _fe->get_dphidz()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);
#endif
              }

            if (family == RATIONAL_BERNSTEIN && order > 1)
              {
                const Gradient rat_grad =
                  rational_test_grad(p, dummy, "", "");

                LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_x),
                                        libmesh_real(rat_grad(0)),
                                        grad_tol);
                if (_dim > 1)
                  LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_y),
                                          libmesh_real(rat_grad(1)),
                                          grad_tol);
                if (_dim > 2)
                  LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_z),
                                          libmesh_real(rat_grad(2)),
                                          grad_tol);
              }
            else if (order > 1)
              {
                const Real & x = p(0);
                const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
                const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

                LIBMESH_ASSERT_FP_EQUAL(2*x+0.125*y+0.0625*z, libmesh_real(grad_u_x),
                                        grad_tol);
                if (_dim > 1)
                  LIBMESH_ASSERT_FP_EQUAL(y+0.125*x+0.03125*z, libmesh_real(grad_u_y),
                                          grad_tol);
                if (_dim > 2)
                  LIBMESH_ASSERT_FP_EQUAL(0.5*z+0.0625*x+0.03125*y, libmesh_real(grad_u_z),
                                          grad_tol);
              }
            else
              {
                LIBMESH_ASSERT_FP_EQUAL(1.0, libmesh_real(grad_u_x),
                                        grad_tol);
                if (_dim > 1)
                  LIBMESH_ASSERT_FP_EQUAL(0.25, libmesh_real(grad_u_y),
                                          grad_tol);
                if (_dim > 2)
                  LIBMESH_ASSERT_FP_EQUAL(0.0625, libmesh_real(grad_u_z),
                                          grad_tol);
              }
          }
#endif
  }


  void testHessU()
  {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

    // Handle the "more processors than elements" case
    if (!_elem)
      return;

    Parameters dummy;

    // These tests require exceptions to be enabled because a
    // TypeTensor::solve() call down in Elem::contains_point()
    // actually throws a non-fatal exception for a certain Point which
    // is not in the Elem. When exceptions are not enabled, this test
    // simply aborts.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
    for (unsigned int i=0; i != _nx; ++i)
      for (unsigned int j=0; j != _ny; ++j)
        for (unsigned int k=0; k != _nz; ++k)
          {
            Point p(Real(i)/_nx);
            if (j > 0)
              p(1) = Real(j)/_ny;
            if (k > 0)
              p(2) = Real(k)/_ny;
            if (!_elem->contains_point(p))
              continue;

            // If at a singular node, cannot use FEMap::map
            if (_elem->local_singular_node(p) != invalid_uint)
              continue;

            std::vector<Point> master_points
              (1, FEMap::inverse_map(_dim, _elem, p));

            _fe->reinit(_elem, &master_points);

            Tensor hess_u;
            for (std::size_t d = 0; d != _dof_indices.size(); ++d)
              hess_u += _fe->get_d2phi()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);

            if (family == RATIONAL_BERNSTEIN && order > 1)
              {
                // TODO: Yeah we'll test the ugly expressions later.
              }
            else if (order > 1)
              {
                LIBMESH_ASSERT_FP_EQUAL(2, libmesh_real(hess_u(0,0)),
                                        hess_tol);
                if (_dim > 1)
                  {
                    LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,1)), libmesh_real(hess_u(1,0)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0.125, libmesh_real(hess_u(0,1)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(1, libmesh_real(hess_u(1,1)),
                                            hess_tol);
                  }
                if (_dim > 2)
                  {
                    LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,2)), libmesh_real(hess_u(2,0)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(1,2)), libmesh_real(hess_u(2,1)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0.0625, libmesh_real(hess_u(0,2)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0.03125, libmesh_real(hess_u(1,2)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0.5, libmesh_real(hess_u(2,2)),
                                            hess_tol);
                  }
              }
            else
              {
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(0,0)),
                                        hess_tol);
                if (_dim > 1)
                  {
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(0,1)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(1,0)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(1,1)),
                                            hess_tol);
                  }
                if (_dim > 2)
                  {
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(0,2)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(1,2)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(2,0)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(2,1)),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u(2,2)),
                                            hess_tol);
                  }
              }
          }
#endif // LIBMESH_ENABLE_EXCEPTIONS
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
  }

  void testHessUComp()
  {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    // Clough-Tocher elements still don't work multithreaded
    if (family == CLOUGH && libMesh::n_threads() > 1)
      return;

    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

    // Handle the "more processors than elements" case
    if (!_elem)
      return;

    Parameters dummy;

    // Why are we seeing O(1e-5) errors??
    const Real hess_tol = sqrt(TOLERANCE);

    // These tests require exceptions to be enabled because a
    // TypeTensor::solve() call down in Elem::contains_point()
    // actually throws a non-fatal exception for a certain Point which
    // is not in the Elem. When exceptions are not enabled, this test
    // simply aborts.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
    for (unsigned int i=0; i != _nx; ++i)
      for (unsigned int j=0; j != _ny; ++j)
        for (unsigned int k=0; k != _nz; ++k)
          {
            Point p(Real(i)/_nx);
            if (j > 0)
              p(1) = Real(j)/_ny;
            if (k > 0)
              p(2) = Real(k)/_ny;
            if (!_elem->contains_point(p))
              continue;

            // If at a singular node, cannot use FEMap::map
            if (_elem->local_singular_node(p) != invalid_uint)
              continue;

            std::vector<Point> master_points
              (1, FEMap::inverse_map(_dim, _elem, p));

            _fe->reinit(_elem, &master_points);

            Number hess_u_xx = 0, hess_u_xy = 0, hess_u_yy = 0,
                   hess_u_xz = 0, hess_u_yz = 0, hess_u_zz = 0;
            for (std::size_t d = 0; d != _dof_indices.size(); ++d)
              {
                hess_u_xx += _fe->get_d2phidx2()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);
#if LIBMESH_DIM > 1
                hess_u_xy += _fe->get_d2phidxdy()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);
                hess_u_yy += _fe->get_d2phidy2()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);
#endif
#if LIBMESH_DIM > 2
                hess_u_xz += _fe->get_d2phidxdz()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);
                hess_u_yz += _fe->get_d2phidydz()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);
                hess_u_zz += _fe->get_d2phidz2()[d][0] * (*_sys->current_local_solution)(_dof_indices[d]);
#endif
              }

            if (family == RATIONAL_BERNSTEIN && order > 1)
              {
                // TODO: tedious calculus
              }
            else if (order > 1)
              {
                LIBMESH_ASSERT_FP_EQUAL(2, libmesh_real(hess_u_xx),
                                        hess_tol);
                if (_dim > 1)
                  {
                    LIBMESH_ASSERT_FP_EQUAL(0.125, libmesh_real(hess_u_xy),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(1, libmesh_real(hess_u_yy),
                                            hess_tol);
                  }
                if (_dim > 2)
                  {
                    LIBMESH_ASSERT_FP_EQUAL( 0.0625, libmesh_real(hess_u_xz),
                                             hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL( 0.03125, libmesh_real(hess_u_yz),
                                             hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL( 0.5, libmesh_real(hess_u_zz),
                                             hess_tol);
                  }
              }
            else
              {
                LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_xx),
                                        hess_tol);
                if (_dim > 1)
                  {
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_xy),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_yy),
                                            hess_tol);
                  }
                if (_dim > 2)
                  {
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_xz),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_yz),
                                            hess_tol);
                    LIBMESH_ASSERT_FP_EQUAL(0, libmesh_real(hess_u_zz),
                                            hess_tol);
                  }
              }
          }
#endif // LIBMESH_ENABLE_EXCEPTIONS
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
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
