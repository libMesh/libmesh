// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/fem_system.h>
#include <libmesh/quadrature.h>
#include <libmesh/diff_solver.h>
#include <libmesh/newmark_solver.h>
#include <libmesh/euler_solver.h>
#include <libmesh/euler2_solver.h>

#include "test_comm.h"

#include "solvers/time_solver_test_common.h"

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

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

  virtual Number F( FEMContext & /*context*/, unsigned int /*qp*/ ) libmesh_override
  { return -2.71; }

  virtual Number C( FEMContext & /*context*/, unsigned int /*qp*/ ) libmesh_override
  { return 0.0; }

  virtual Number M( FEMContext & /*context*/, unsigned int /*qp*/ ) libmesh_override
  { return 3.14; }

  virtual Number u( Real t ) libmesh_override
  { return 2.71/3.14*0.5*t*t; }
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

  virtual Number F( FEMContext & context, unsigned int /*qp*/ ) libmesh_override
  { return -6.0*context.get_time()-2.0; }

  virtual Number C( FEMContext & /*context*/, unsigned int /*qp*/ ) libmesh_override
  { return 0.0; }

  virtual Number M( FEMContext & /*context*/, unsigned int /*qp*/ ) libmesh_override
  { return 1.0; }

  virtual Number u( Real t ) libmesh_override
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

  virtual void aux_time_solver_init( NewmarkSolver & time_solver ) libmesh_override
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
  CPPUNIT_TEST_SUITE( NewmarkSolverTest );

  CPPUNIT_TEST( testNewmarkSolverConstantSecondOrderODESecondOrderStyle );
  CPPUNIT_TEST( testNewmarkSolverLinearTimeSecondOrderODESecondOrderStyle );
  CPPUNIT_TEST( testNewmarkSolverConstantSecondOrderODEFirstOrderStyle );
  CPPUNIT_TEST( testNewmarkSolverLinearTimeSecondOrderODEFirstOrderStyle );

  CPPUNIT_TEST_SUITE_END();

public:

  void testNewmarkSolverConstantSecondOrderODESecondOrderStyle()
  {
    this->run_test_with_exact_soln<ConstantSecondOrderODE<SecondOrderScalarSystemSecondOrderTimeSolverBase> >(0.5,10);
  }

  void testNewmarkSolverLinearTimeSecondOrderODESecondOrderStyle()
  {
    // For \beta = 1/6, we have the "linear acceleration method" for which
    // we should be able to exactly integrate linear (in time) acceleration
    // functions.
    this->set_beta(1.0/6.0);
    this->run_test_with_exact_soln<LinearTimeSecondOrderODE<SecondOrderScalarSystemSecondOrderTimeSolverBase> >(0.5,10);
  }

  void testNewmarkSolverConstantSecondOrderODEFirstOrderStyle()
  {
    this->run_test_with_exact_soln<ConstantSecondOrderODE<SecondOrderScalarSystemFirstOrderTimeSolverBase> >(0.5,10);
  }

  void testNewmarkSolverLinearTimeSecondOrderODEFirstOrderStyle()
  {
    // For \beta = 1/6, we have the "linear acceleration method" for which
    // we should be able to exactly integrate linear (in time) acceleration
    // functions.
    this->set_beta(1.0/6.0);
    this->run_test_with_exact_soln<LinearTimeSecondOrderODE<SecondOrderScalarSystemFirstOrderTimeSolverBase> >(0.5,10);
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
  CPPUNIT_TEST_SUITE( EulerSolverSecondOrderTest );

  CPPUNIT_TEST( testEulerSolverConstantSecondOrderODE );

  CPPUNIT_TEST_SUITE_END();

public:

  void testEulerSolverConstantSecondOrderODE()
  {
    this->set_theta(0.5);
    this->run_test_with_exact_soln<ConstantSecondOrderODE<SecondOrderScalarSystemFirstOrderTimeSolverBase> >(0.5,10);
  }

};

class Euler2SolverSecondOrderTest : public CppUnit::TestCase,
                                    public ThetaSolverTestBase<Euler2Solver>
{
public:
  CPPUNIT_TEST_SUITE( Euler2SolverSecondOrderTest );

  CPPUNIT_TEST( testEuler2SolverConstantSecondOrderODE );

  CPPUNIT_TEST_SUITE_END();

public:
  void testEuler2SolverConstantSecondOrderODE()
  {
    this->set_theta(0.5);
    this->run_test_with_exact_soln<ConstantSecondOrderODE<SecondOrderScalarSystemFirstOrderTimeSolverBase> >(0.5,10);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( NewmarkSolverTest );
CPPUNIT_TEST_SUITE_REGISTRATION( EulerSolverSecondOrderTest );
CPPUNIT_TEST_SUITE_REGISTRATION( Euler2SolverSecondOrderTest );
