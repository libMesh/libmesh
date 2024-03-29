#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/fem_system.h>
#include <libmesh/quadrature.h>
#include <libmesh/diff_solver.h>
#include <libmesh/newmark_solver.h>
#include <libmesh/euler_solver.h>
#include <libmesh/euler2_solver.h>

#include "solvers/time_solver_test_common.h"


//! Implements ODE: 3.14\ddot{u} = 2.71, u(0) = 0, \dot{u}(0) = 0
template<typename SystemBase>
class ConstantSecondOrderODE : public SystemBase
{
public:
  ConstantSecondOrderODE(EquationSystems & es,
                         const std::string & name_in,
                         const unsigned int number_in)
    : SystemBase(es, name_in, number_in)
  {}

  virtual Number F( FEMContext & /*context*/, unsigned int /*qp*/ ) override
  { return -Real(271)/100; }

  virtual Number C( FEMContext & /*context*/, unsigned int /*qp*/ ) override
  { return 0.0; }

  virtual Number M( FEMContext & /*context*/, unsigned int /*qp*/ ) override
  { return Real(314)/100; }

  virtual Number u( Real t ) override
  { return Real(271)/Real(314)*0.5*t*t; }
};

//! Implements ODE: 1.0\ddot{u} = 6.0*t+2.0, u(0) = 0, \dot{u}(0) = 0
template<typename SystemBase>
class LinearTimeSecondOrderODE : public SystemBase
{
public:
  LinearTimeSecondOrderODE(EquationSystems & es,
                           const std::string & name_in,
                           const unsigned int number_in)
    : SystemBase(es, name_in, number_in)
  {}

  virtual Number F( FEMContext & context, unsigned int /*qp*/ ) override
  { return -6.0*context.get_time()-2.0; }

  virtual Number C( FEMContext & /*context*/, unsigned int /*qp*/ ) override
  { return 0.0; }

  virtual Number M( FEMContext & /*context*/, unsigned int /*qp*/ ) override
  { return 1.0; }

  virtual Number u( Real t ) override
  { return t*t*t+t*t; }
};

class NewmarkSolverTestBase : public TimeSolverTestImplementation<NewmarkSolver>
{
public:
  NewmarkSolverTestBase()
    : TimeSolverTestImplementation<NewmarkSolver>(),
    _beta(0.25)
  {}

protected:

  virtual void aux_time_solver_init( NewmarkSolver & time_solver ) override
  { time_solver.set_beta(_beta);
    time_solver.compute_initial_accel(); }

  void set_beta( Real beta )
  { _beta = beta; }

  Real _beta;
};

class NewmarkSolverTest : public CppUnit::TestCase,
                          public NewmarkSolverTestBase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( NewmarkSolverTest );

#ifdef LIBMESH_HAVE_SOLVER
  CPPUNIT_TEST( testNewmarkSolverConstantSecondOrderODESecondOrderStyle );
  CPPUNIT_TEST( testNewmarkSolverLinearTimeSecondOrderODESecondOrderStyle );
  CPPUNIT_TEST( testNewmarkSolverConstantSecondOrderODEFirstOrderStyle );
  CPPUNIT_TEST( testNewmarkSolverLinearTimeSecondOrderODEFirstOrderStyle );
#endif

  CPPUNIT_TEST_SUITE_END();

public:

  void testNewmarkSolverConstantSecondOrderODESecondOrderStyle()
  {
    LOG_UNIT_TEST;

    this->run_test_with_exact_soln<ConstantSecondOrderODE<SecondOrderScalarSystemSecondOrderTimeSolverBase>>(0.5,10);
  }

  void testNewmarkSolverLinearTimeSecondOrderODESecondOrderStyle()
  {
    LOG_UNIT_TEST;

    // For \beta = 1/6, we have the "linear acceleration method" for which
    // we should be able to exactly integrate linear (in time) acceleration
    // functions.
    this->set_beta(Real(1)/Real(6));
    this->run_test_with_exact_soln<LinearTimeSecondOrderODE<SecondOrderScalarSystemSecondOrderTimeSolverBase>>(0.5,10);
  }

  void testNewmarkSolverConstantSecondOrderODEFirstOrderStyle()
  {
    LOG_UNIT_TEST;

    this->run_test_with_exact_soln<ConstantSecondOrderODE<SecondOrderScalarSystemFirstOrderTimeSolverBase>>(0.5,10);
  }

  void testNewmarkSolverLinearTimeSecondOrderODEFirstOrderStyle()
  {
    LOG_UNIT_TEST;

    // For \beta = 1/6, we have the "linear acceleration method" for which
    // we should be able to exactly integrate linear (in time) acceleration
    // functions.
    this->set_beta(Real(1)/Real(6));
    this->run_test_with_exact_soln<LinearTimeSecondOrderODE<SecondOrderScalarSystemFirstOrderTimeSolverBase>>(0.5,10);
  }

};

template<typename TimeSolverType>
class ThetaSolverTestBase : public TimeSolverTestImplementation<TimeSolverType>
{
public:
  ThetaSolverTestBase()
    : TimeSolverTestImplementation<TimeSolverType>(),
    _theta(1.0)
  {}

protected:

  virtual void aux_time_solver_init( TimeSolverType & time_solver )
  { time_solver.theta = _theta; }

  void set_theta( Real theta )
  { _theta = theta; }

  Real _theta;
};

class EulerSolverSecondOrderTest : public CppUnit::TestCase,
                                   public ThetaSolverTestBase<EulerSolver>
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( EulerSolverSecondOrderTest );

#ifdef LIBMESH_HAVE_SOLVER
  CPPUNIT_TEST( testEulerSolverConstantSecondOrderODE );
#endif

  CPPUNIT_TEST_SUITE_END();

public:

  void testEulerSolverConstantSecondOrderODE()
  {
    LOG_UNIT_TEST;

    this->set_theta(0.5);
    this->run_test_with_exact_soln<ConstantSecondOrderODE<SecondOrderScalarSystemFirstOrderTimeSolverBase>>(0.5,10);
  }

};

class Euler2SolverSecondOrderTest : public CppUnit::TestCase,
                                    public ThetaSolverTestBase<Euler2Solver>
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( Euler2SolverSecondOrderTest );

#ifdef LIBMESH_HAVE_SOLVER
  CPPUNIT_TEST( testEuler2SolverConstantSecondOrderODE );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void testEuler2SolverConstantSecondOrderODE()
  {
    LOG_UNIT_TEST;

    this->set_theta(0.5);
    this->run_test_with_exact_soln<ConstantSecondOrderODE<SecondOrderScalarSystemFirstOrderTimeSolverBase>>(0.5,10);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( NewmarkSolverTest );
CPPUNIT_TEST_SUITE_REGISTRATION( EulerSolverSecondOrderTest );
CPPUNIT_TEST_SUITE_REGISTRATION( Euler2SolverSecondOrderTest );
