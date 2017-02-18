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

// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <vector>

#define FETEST                   \
  CPPUNIT_TEST( testU );         \
  CPPUNIT_TEST( testGradU );     \
  CPPUNIT_TEST( testGradUComp );

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

public:
  void setUp()
  {
    _mesh = new Mesh(*TestCommWorld);
    const UniquePtr<Elem> test_elem = Elem::build(elem_type);
    _dim = test_elem->dim();
    const unsigned int ny = _dim > 1;
    const unsigned int nz = _dim > 2;

    MeshTools::Generation::build_cube (*_mesh,
                                      1, ny, nz,
                                      0., 1., 0., ny, 0., nz,
                                      elem_type);

    _es = new EquationSystems(*_mesh);
    _sys = &(_es->add_system<System> ("SimpleSystem"));
    _sys->add_variable("u", order, family);
    _es->init();
    _sys->project_solution(linear_test, linear_test_grad, _es->parameters);

    FEType fe_type = _sys->variable_type(0);
    _fe = FEBase::build(_dim, fe_type).release();
    _fe->get_phi();
    _fe->get_dphi();
    _fe->get_dphidx();
#if LIBMESH_DIM > 1
    _fe->get_dphidy();
#endif
#if LIBMESH_DIM > 2
    _fe->get_dphidz();
#endif

    MeshBase::const_element_iterator
      elem_it  = _mesh->active_local_elements_begin(),
      elem_end = _mesh->active_local_elements_end();

    if (elem_it == elem_end)
      _elem = libmesh_nullptr;
    else
      _elem = *elem_it;

    _sys->get_dof_map().dof_indices(_elem, _dof_indices);

    _nx = 10;
    _ny = _nx * (_dim > 1);
    _nz = _nx * (_dim > 2);
  }

  void tearDown()
  {
    delete _fe;
    delete _es;
    delete _mesh;
  }

  void testU()
  {
    // Handle the "more processors than elements" case
    if (!_elem)
      return;

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

            std::vector<Point> master_points
              (1, FEInterface::inverse_map(_dim, _fe->get_fe_type(), _elem, p));

            _fe->reinit(_elem, &master_points);

            Number u = 0;
            for (std::size_t d = 0; d != _dof_indices.size(); ++d)
              u += _fe->get_phi()[d][0] *
                   (*_sys->current_local_solution)(_dof_indices[d]);

            CPPUNIT_ASSERT_DOUBLES_EQUAL
              (libmesh_real(u),
               libmesh_real(x + 0.25*y + 0.0625*z),
               TOLERANCE*TOLERANCE);
          }
  }

  void testGradU()
  {
    // Handle the "more processors than elements" case
    if (!_elem)
      return;

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

            std::vector<Point> master_points
              (1, FEInterface::inverse_map(_dim, _fe->get_fe_type(), _elem, p));

            _fe->reinit(_elem, &master_points);

            Gradient grad_u = 0;
            for (std::size_t d = 0; d != _dof_indices.size(); ++d)
              grad_u += _fe->get_dphi()[d][0] *
                        (*_sys->current_local_solution)(_dof_indices[d]);

            CPPUNIT_ASSERT_DOUBLES_EQUAL
              (libmesh_real(grad_u(0)), 1.0,
               TOLERANCE*sqrt(TOLERANCE));
            if (_dim > 1)
              CPPUNIT_ASSERT_DOUBLES_EQUAL
                (libmesh_real(grad_u(1)), 0.25,
                 TOLERANCE*sqrt(TOLERANCE));
            if (_dim > 2)
              CPPUNIT_ASSERT_DOUBLES_EQUAL
                (libmesh_real(grad_u(2)), 0.0625,
                 TOLERANCE*sqrt(TOLERANCE));
          }
  }

  void testGradUComp()
  {
    // Handle the "more processors than elements" case
    if (!_elem)
      return;

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

            std::vector<Point> master_points
              (1, FEInterface::inverse_map(_dim, _fe->get_fe_type(), _elem, p));

            _fe->reinit(_elem, &master_points);

            Number grad_u_x = 0, grad_u_y = 0, grad_u_z = 0;
            for (std::size_t d = 0; d != _dof_indices.size(); ++d)
              {
                grad_u_x += _fe->get_dphidx()[d][0] *
                            (*_sys->current_local_solution)(_dof_indices[d]);
#if LIBMESH_DIM > 1
                grad_u_y += _fe->get_dphidy()[d][0] *
                            (*_sys->current_local_solution)(_dof_indices[d]);
#endif
#if LIBMESH_DIM > 2
                grad_u_z += _fe->get_dphidz()[d][0] *
                            (*_sys->current_local_solution)(_dof_indices[d]);
#endif
              }

            CPPUNIT_ASSERT_DOUBLES_EQUAL
              (libmesh_real(grad_u_x), 1.0,
               TOLERANCE*sqrt(TOLERANCE));
            if (_dim > 1)
              CPPUNIT_ASSERT_DOUBLES_EQUAL
                (libmesh_real(grad_u_y), 0.25,
                 TOLERANCE*sqrt(TOLERANCE));
            if (_dim > 2)
              CPPUNIT_ASSERT_DOUBLES_EQUAL
                (libmesh_real(grad_u_z), 0.0625,
                 TOLERANCE*sqrt(TOLERANCE));
          }
  }

};

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We'll put an
// ignore_warnings at the end of this file so it's the last warnings
// related header that our including code sees.
#include <libmesh/ignore_warnings.h>

#define INSTANTIATE_FETEST(order, family, elemtype) \
class FETest_##order##_##family##_##elemtype : public FETest<order, family, elemtype> { \
public: \
  CPPUNIT_TEST_SUITE( FETest_##order##_##family##_##elemtype ); \
  FETEST \
  CPPUNIT_TEST_SUITE_END(); \
}; \
 \
CPPUNIT_TEST_SUITE_REGISTRATION( FETest_##order##_##family##_##elemtype );

#endif // #ifdef __fe_test_h__
