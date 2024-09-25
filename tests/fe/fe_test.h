#ifndef FE_TEST_H
#define FE_TEST_H

#include "test_comm.h"

#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe_base.h>
#include <libmesh/fe_interface.h>
#include <libmesh/function_base.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/system.h>
#include <libmesh/quadrature_gauss.h>

#include <vector>

#include "libmesh_cppunit.h"

#define FETEST                                  \
  CPPUNIT_TEST( testFEInterface );              \
  CPPUNIT_TEST( testU );                        \
  CPPUNIT_TEST( testPartitionOfUnity );         \
  CPPUNIT_TEST( testGradU );                    \
  CPPUNIT_TEST( testGradUComp );                \
  CPPUNIT_TEST( testHessU );                    \
  CPPUNIT_TEST( testHessUComp );                \
  CPPUNIT_TEST( testDualDoesntScreamAndDie );   \
  CPPUNIT_TEST( testCustomReinit );

using namespace libMesh;


class SkewFunc : public FunctionBase<Real>
{
  std::unique_ptr<FunctionBase<Real>> clone () const override
  { return std::make_unique<SkewFunc>(); }

  Real operator() (const Point &,
                   const Real = 0.) override
  { libmesh_not_implemented(); } // scalar-only API

  // Skew in x based on y, y based on z
  void operator() (const Point & p,
                   const Real,
                   DenseVector<Real> & output)
  {
    output.resize(3);
    output(0) = p(0);
#if LIBMESH_DIM > 0
    output(0) += .1*p(1);
    output(1) = p(1);
#endif
#if LIBMESH_DIM > 1
    output(1) += .2*p(2);
    output(2) = p(2);
#endif
  }
};


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


inline
Number fe_quartic_test (const Point& p,
                        const Parameters&,
                        const std::string&,
                        const std::string&)
{
  const Real & x = p(0);
  const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
  const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

  return x*x*(1-x)*(1-x) + x*x*z*(1-y) + x*(1-x)*(1-y)*(1-z) + (1-x)*y*(1-y)*z + z*z*(1-z)*(1-z);
}

inline
Gradient fe_quartic_test_grad (const Point & p,
                               const Parameters&,
                               const std::string&,
                               const std::string&)
{
  const Real & x = p(0);
  const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
  const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

  Gradient grad = 4*x*x*x-6*x*x+2*x + 2*x*z*(1-y) + (1-2*x)*(1-y)*(1-z) - y*(1-y)*z;
  if (LIBMESH_DIM > 1)
    grad(1) = -x*x*z - x*(1-x)*(1-z) + (1-x)*(1-2*y)*z;
  if (LIBMESH_DIM > 2)
    grad(2) = x*x*(1-y) - x*(1-x)*(1-y) + (1-x)*y*(1-y) + 4*z*z*z-6*z*z+2*z;

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
  (((family != LAGRANGE && family != L2_LAGRANGE) || \
    (elem_type != TRI7 && elem_type != TET14 && \
     elem_type != PRISM20 && elem_type != PRISM21 && \
     elem_type != PYRAMID18)) && order > 2)


template <Order order, FEFamily family, ElemType elem_type, unsigned int build_nx>
class FETestBase : public CppUnit::TestCase {

protected:
  std::string libmesh_suite_name;

  unsigned int _dim, _nx, _ny, _nz;
  Elem *_elem;
  std::vector<dof_id_type> _dof_indices;
  System * _sys;
  std::unique_ptr<Mesh> _mesh;
  std::unique_ptr<EquationSystems> _es;
  std::unique_ptr<FEBase> _fe;
  std::unique_ptr<QGauss> _qrule;

  Real _value_tol, _grad_tol, _hess_tol;


  static RealGradient true_gradient(Point p)
  {
    Parameters dummy;

    Gradient true_grad;
    RealGradient returnval;

    if (family == RATIONAL_BERNSTEIN && order > 1)
      true_grad = rational_test_grad(p, dummy, "", "");
    else if (order > 3)
      true_grad = fe_quartic_test_grad(p, dummy, "", "");
    else if (FE_CAN_TEST_CUBIC)
      true_grad = fe_cubic_test_grad(p, dummy, "", "");
    else if (order > 1)
      {
        const Real & x = p(0);
        const Real & y = (LIBMESH_DIM > 1) ? p(1) : 0;
        const Real & z = (LIBMESH_DIM > 2) ? p(2) : 0;

        true_grad = Gradient(2*x+0.125*y+0.0625*z,
                             y+0.125*x+0.03125*z,
                             0.5*z+0.0625*x+0.03125*y);
      }
    else
      true_grad = Gradient(1.0, 0.25, 0.0625);

    for (unsigned int d=0; d != LIBMESH_DIM; ++d)
      {
        CPPUNIT_ASSERT(true_grad(d) ==
                       Number(libmesh_real(true_grad(d))));

        returnval(d) = libmesh_real(true_grad(d));
      }

    return returnval;
  }


  static RealTensor true_hessian(Point p)
  {
    const Real & x = p(0);
    const Real & y = LIBMESH_DIM > 1 ? p(1) : 0;
    const Real & z = LIBMESH_DIM > 2 ? p(2) : 0;

    if (order > 3)
      return RealTensor
        { 12*x*x-12*x+2+2*z*(1-y)-2*(1-y)*(1-z), -2*x*z-(1-2*x)*(1-z)-(1-2*y)*z, 2*x*(1-y)-(1-2*x)*(1-y)-y*(1-y),
                 -2*x*z-(1-2*x)*(1-z)-(1-2*y)*z,                     -2*(1-x)*z,      -x*x+x*(1-x)+(1-x)*(1-2*y),
                2*x*(1-y)-(1-2*x)*(1-y)-y*(1-y),     -x*x+x*(1-x)+(1-x)*(1-2*y),                   12*z*z-12*z+2 };
    else if (FE_CAN_TEST_CUBIC)
      return RealTensor
        { 6*x-4+2*(1-y), -2*x+z-1,     y-1,
               -2*x+z-1,     -2*z, x+1-2*y,
                    y-1,  x+1-2*y, 6*z-4 };
    else if (order > 1)
      return RealTensor
        { 2, 0.125, 0.0625,
          0.125, 1, 0.03125,
          0.0625, 0.03125, 0.5 };

    return RealTensor
      { 0, 0, 0,
        0, 0, 0,
        0, 0, 0 };
  }

public:
  void setUp()
  {
    _mesh = std::make_unique<Mesh>(*TestCommWorld);
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

    // For debugging purposes it can be helpful to only consider one
    // element even when we're using an element type that requires
    // more than one element to fill out a square or cube.
#if 0
    for (dof_id_type i = 0; i != _mesh->max_elem_id(); ++i)
      {
        Elem * elem = _mesh->query_elem_ptr(i);
        if (elem && elem->id())
          _mesh->delete_elem(elem);
      }
    _mesh->prepare_for_use();
    CPPUNIT_ASSERT_EQUAL(_mesh->n_elem(), dof_id_type(1));
#endif

    // Permute our elements randomly and rotate and skew our mesh so
    // we test all sorts of orientations ... except with Hermite
    // elements, which are only designed to support meshes with a
    // single orientation shared by all elements.  We're also not
    // rotating and/or skewing the rational elements, since our test
    // solution was designed for a specific weighted mesh.
    if (family != HERMITE &&
        family != RATIONAL_BERNSTEIN)
      {
        MeshTools::Modification::permute_elements(*_mesh);

        // Not yet testing manifolds embedded in higher-D space
        if (_dim > 1)
          MeshTools::Modification::rotate(*_mesh, 4,
                                          8*(_dim>2), 16*(_dim>2));

        SkewFunc skew_func;
        MeshTools::Modification::redistribute(*_mesh, skew_func);
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

    _es = std::make_unique<EquationSystems>(*_mesh);
    _sys = &(_es->add_system<System> ("SimpleSystem"));
    _sys->add_variable("u", order, family);
    _es->init();

    if (family == RATIONAL_BERNSTEIN && order > 1)
      {
        _sys->project_solution(rational_test, rational_test_grad, _es->parameters);
      }
    else if (order > 3)
      {
        _sys->project_solution(fe_quartic_test, fe_quartic_test_grad, _es->parameters);
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
    _fe = FEBase::build(_dim, fe_type);

    // Create quadrature rule for use in computing dual shape coefficients
    _qrule = std::make_unique<QGauss>(_dim, fe_type.default_quadrature_order());
    _fe->attach_quadrature_rule(_qrule.get());

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
    // On Tet14 Monomial gives us at least 13*tol*sqrt(tol) in some
    // cases
    this->_grad_tol = 15 * TOLERANCE * sqrt(TOLERANCE);

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

  void tearDown() {}
};



template <Order order, FEFamily family, ElemType elem_type>
class FETest : public FETestBase<order, family, elem_type, 1> {

public:

  template <typename Functor>
  void testLoop(Functor f)
  {
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
            Point p = Real(i)/this->_nx;
            if (j > 0)
              p(1) = Real(j)/this->_ny;
            if (k > 0)
              p(2) = Real(k)/this->_nz;
            if (!this->_elem->contains_point(p))
              continue;

            // If at a singular node, cannot use FEMap::map
            if (this->_elem->local_singular_node(p) != invalid_uint)
              continue;

            std::vector<Point> master_points
              (1, FEMap::inverse_map(this->_dim, this->_elem, p));

            // Reinit at point to test against analytic solution
            this->_fe->reinit(this->_elem, &master_points);

            f(p);
          }
#endif // LIBMESH_ENABLE_EXCEPTIONS
  }

  void testPartitionOfUnity()
    {
      if (!this->_elem)
        return;

      this->_fe->reinit(this->_elem);

      bool satisfies_partition_of_unity = true;
      for (const auto qp : make_range(this->_qrule->n_points()))
      {
        Real phi_sum = 0;
        for (std::size_t d = 0; d != this->_dof_indices.size(); ++d)
          phi_sum += this->_fe->get_phi()[d][qp];
        if (phi_sum < (1 - TOLERANCE) || phi_sum > (1 + TOLERANCE))
        {
          satisfies_partition_of_unity = false;
          break;
        }
      }

      switch (this->_fe->get_family())
      {
        case MONOMIAL:
        {
          switch (this->_fe->get_order())
          {
            case CONSTANT:
              CPPUNIT_ASSERT(satisfies_partition_of_unity);
              break;

            default:
              CPPUNIT_ASSERT(!satisfies_partition_of_unity);
              break;
          }
          break;
        }
        case SZABAB:
        case HIERARCHIC:
        case L2_HIERARCHIC:
        {
          switch (this->_fe->get_order())
          {
            case FIRST:
              CPPUNIT_ASSERT(satisfies_partition_of_unity);
              break;

            default:
              CPPUNIT_ASSERT(!satisfies_partition_of_unity);
              break;
          }
          break;
        }

        case XYZ:
        case CLOUGH:
        case HERMITE:
        {
          CPPUNIT_ASSERT(!satisfies_partition_of_unity);
          break;
        }

        case LAGRANGE:
        case L2_LAGRANGE:
        case BERNSTEIN:
        case RATIONAL_BERNSTEIN:
        {
          CPPUNIT_ASSERT(satisfies_partition_of_unity);
          break;
        }

        default:
          CPPUNIT_FAIL("Uncovered FEFamily");
      }
    }


  void testFEInterface()
  {
    LOG_UNIT_TEST;

    // Handle the "more processors than elements" case
    if (!this->_elem)
      return;

    this->_fe->reinit(this->_elem);

    const FEType fe_type = this->_sys->variable_type(0);

    CPPUNIT_ASSERT_EQUAL(
      FEInterface::n_shape_functions(fe_type, this->_elem),
      this->_fe->n_shape_functions());

    CPPUNIT_ASSERT_EQUAL(
      FEInterface::get_continuity(fe_type),
      this->_fe->get_continuity());

    CPPUNIT_ASSERT_EQUAL(
      FEInterface::is_hierarchic(fe_type),
      this->_fe->is_hierarchic());
  }

  void testU()
  {
    LOG_UNIT_TEST;

    auto f = [this](Point p)
      {
        Parameters dummy;

        Number u = 0;
        for (std::size_t d = 0; d != this->_dof_indices.size(); ++d)
          u += this->_fe->get_phi()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);

        Number true_u;

        if (family == RATIONAL_BERNSTEIN && order > 1)
          true_u = rational_test(p, dummy, "", "");
        else if (order > 3)
          true_u = fe_quartic_test(p, dummy, "", "");
        else if (FE_CAN_TEST_CUBIC)
          true_u = fe_cubic_test(p, dummy, "", "");
        else if (order > 1)
          true_u = p(0)*p(0) + 0.5*p(1)*p(1) + 0.25*p(2)*p(2) +
            0.125*p(0)*p(1) + 0.0625*p(0)*p(2) + 0.03125*p(1)*p(2);
        else
          true_u = p(0) + 0.25*p(1) + 0.0625*p(2);

        LIBMESH_ASSERT_FP_EQUAL
          (libmesh_real(true_u), libmesh_real(u), this->_value_tol);
      };

    testLoop(f);
  }

  void testDualDoesntScreamAndDie()
  {
    LOG_UNIT_TEST;

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
    LOG_UNIT_TEST;

    auto f = [this](Point p)
      {
        Parameters dummy;

        Gradient grad_u = 0;
        for (std::size_t d = 0; d != this->_dof_indices.size(); ++d)
          grad_u += this->_fe->get_dphi()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);

        RealGradient true_grad = this->true_gradient(p);

        LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(0)),
                                true_grad(0), this->_grad_tol);
        if (this->_dim > 1)
          LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(1)),
                                  true_grad(1), this->_grad_tol);
        if (this->_dim > 2)
          LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u(2)),
                                  true_grad(2), this->_grad_tol);
      };

    testLoop(f);
  }

  void testGradUComp()
  {
    LOG_UNIT_TEST;

    auto f = [this](Point p)
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

        RealGradient true_grad = this->true_gradient(p);

        LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_x),
                                true_grad(0), this->_grad_tol);
        if (this->_dim > 1)
          LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_y),
                                  true_grad(1), this->_grad_tol);
        if (this->_dim > 2)
          LIBMESH_ASSERT_FP_EQUAL(libmesh_real(grad_u_z),
                                  true_grad(2), this->_grad_tol);
      };

    testLoop(f);
  }


  void testHessU()
  {
    LOG_UNIT_TEST;

    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    auto f = [this](Point p)
      {
        Tensor hess_u;
        for (std::size_t d = 0; d != this->_dof_indices.size(); ++d)
          hess_u += this->_fe->get_d2phi()[d][0] * (*this->_sys->current_local_solution)(this->_dof_indices[d]);

        // TODO: Yeah we'll test the ugly expressions later.
        if (family == RATIONAL_BERNSTEIN && order > 1)
          return;

        RealTensor true_hess = this->true_hessian(p);

        LIBMESH_ASSERT_FP_EQUAL(true_hess(0,0), libmesh_real(hess_u(0,0)),
                                this->_hess_tol);
        if (this->_dim > 1)
          {
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,1)), libmesh_real(hess_u(1,0)),
                                    this->_hess_tol);
            LIBMESH_ASSERT_FP_EQUAL(true_hess(0,1), libmesh_real(hess_u(0,1)),
                                    this->_hess_tol);
            LIBMESH_ASSERT_FP_EQUAL(true_hess(1,1), libmesh_real(hess_u(1,1)),
                                    this->_hess_tol);
          }
        if (this->_dim > 2)
          {
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(0,2)), libmesh_real(hess_u(2,0)),
                                    this->_hess_tol);
            LIBMESH_ASSERT_FP_EQUAL(libmesh_real(hess_u(1,2)), libmesh_real(hess_u(2,1)),
                                    this->_hess_tol);
            LIBMESH_ASSERT_FP_EQUAL(true_hess(0,2), libmesh_real(hess_u(0,2)),
                                    this->_hess_tol);
            LIBMESH_ASSERT_FP_EQUAL(true_hess(1,2), libmesh_real(hess_u(1,2)),
                                    this->_hess_tol);
            LIBMESH_ASSERT_FP_EQUAL(true_hess(2,2), libmesh_real(hess_u(2,2)),
                                    this->_hess_tol);
          }
      };

    testLoop(f);
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
  }

  void testHessUComp()
  {
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    LOG_UNIT_TEST;

    // Szabab elements don't have second derivatives yet
    if (family == SZABAB)
      return;

    auto f = [this](Point p)
      {
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

        // TODO: Yeah we'll test the ugly expressions later.
        if (family == RATIONAL_BERNSTEIN && order > 1)
          return;

        RealTensor true_hess = this->true_hessian(p);

        LIBMESH_ASSERT_FP_EQUAL(true_hess(0,0), libmesh_real(hess_u_xx),
                                this->_hess_tol);
        if (this->_dim > 1)
          {
            LIBMESH_ASSERT_FP_EQUAL(true_hess(0,1), libmesh_real(hess_u_xy),
                                    this->_hess_tol);
            LIBMESH_ASSERT_FP_EQUAL(true_hess(1,1), libmesh_real(hess_u_yy),
                                    this->_hess_tol);
          }
        if (this->_dim > 2)
          {
            LIBMESH_ASSERT_FP_EQUAL(true_hess(0,2), libmesh_real(hess_u_xz),
                                    this->_hess_tol);
            LIBMESH_ASSERT_FP_EQUAL(true_hess(1,2), libmesh_real(hess_u_yz),
                                    this->_hess_tol);
            LIBMESH_ASSERT_FP_EQUAL(true_hess(2,2), libmesh_real(hess_u_zz),
                                    this->_hess_tol);
          }
      };

    testLoop(f);
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
  }

  void testCustomReinit()
  {
    LOG_UNIT_TEST;

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
  FETest_##order##_##family##_##elemtype() :                            \
    FETest<order,family,elemtype>() {                                   \
    if (unitlog->summarized_logs_enabled())                             \
      this->libmesh_suite_name = "FETest";                              \
    else                                                                \
      this->libmesh_suite_name = "FETest_" #order "_" #family "_" #elemtype; \
  }                                                                     \
  CPPUNIT_TEST_SUITE( FETest_##order##_##family##_##elemtype );         \
  FETEST                                                                \
  CPPUNIT_TEST_SUITE_END();                                             \
  };                                                                    \
                                                                        \
  CPPUNIT_TEST_SUITE_REGISTRATION( FETest_##order##_##family##_##elemtype );

#endif
