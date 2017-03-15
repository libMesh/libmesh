// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/quadrature.h>
#include <libmesh/diff_solver.h>
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


//! Implements ODE: 2.1\dot{u} = 5, u(0) = 0;
class ConstantFirstOrderODE : public FirstOrderScalarSystemBase
{
public:
  ConstantFirstOrderODE(EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in)
    : FirstOrderScalarSystemBase(es, name_in, number_in)
  {}

  virtual Number F( FEMContext & /*context*/, unsigned int /*qp*/ )
  { return 5.0; }

  virtual Number M( FEMContext & /*context*/, unsigned int /*qp*/ )
  { return 2.1; }

  virtual Number u( Real t )
  { return 5.0/2.1*t; }
};

//! Implements ODE: 5.0\dot{u} = 2.0t, u(0) = 0;
class LinearTimeFirstOrderODE : public FirstOrderScalarSystemBase
{
public:
  LinearTimeFirstOrderODE(EquationSystems & es,
                          const std::string & name_in,
                          const unsigned int number_in)
    : FirstOrderScalarSystemBase(es, name_in, number_in)
  {}

  virtual Number F( FEMContext & context, unsigned int /*qp*/ )
  { return 2.0*context.get_time(); }

  virtual Number M( FEMContext & /*context*/, unsigned int /*qp*/ )
  { return 5.0; }

  virtual Number u( Real t )
  { return 1.0/5.0*t*t; }
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

class EulerSolverTest : public CppUnit::TestCase,
                        public ThetaSolverTestBase<EulerSolver>
{
public:
  CPPUNIT_TEST_SUITE( EulerSolverTest );

  CPPUNIT_TEST( testEulerSolverConstantFirstOrderODE );
  CPPUNIT_TEST( testEulerSolverLinearTimeFirstOrderODE );

  CPPUNIT_TEST_SUITE_END();

public:

  void testEulerSolverConstantFirstOrderODE()
  {
    this->set_theta(1.0);
    this->run_test_with_exact_soln<ConstantFirstOrderODE>(0.5,10);

    this->set_theta(0.5);
    this->run_test_with_exact_soln<ConstantFirstOrderODE>(0.5,10);
  }

  void testEulerSolverLinearTimeFirstOrderODE()
  {
    // Need \theta = 0.5 since this has t in F.
    this->set_theta(0.5);
    this->run_test_with_exact_soln<LinearTimeFirstOrderODE>(0.5,10);
  }

};

class Euler2SolverTest : public CppUnit::TestCase,
                         public ThetaSolverTestBase<Euler2Solver>
{
public:
  CPPUNIT_TEST_SUITE( Euler2SolverTest );

  CPPUNIT_TEST( testEuler2SolverConstantFirstOrderODE );
  CPPUNIT_TEST( testEuler2SolverLinearTimeFirstOrderODE );

  CPPUNIT_TEST_SUITE_END();

public:
  void testEuler2SolverConstantFirstOrderODE()
  {
    this->set_theta(1.0);
    this->run_test_with_exact_soln<ConstantFirstOrderODE>(0.5,10);

    this->set_theta(0.5);
    this->run_test_with_exact_soln<ConstantFirstOrderODE>(0.5,10);
  }

  void testEuler2SolverLinearTimeFirstOrderODE()
  {
    // Need \theta = 0.5 since this has t in F.
    this->set_theta(0.5);
    this->run_test_with_exact_soln<LinearTimeFirstOrderODE>(0.5,10);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( EulerSolverTest );
CPPUNIT_TEST_SUITE_REGISTRATION( Euler2SolverTest );
