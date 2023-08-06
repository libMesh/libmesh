#include "libmesh/shell_matrix.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/sparse_shell_matrix.h"
#include "libmesh/solver_configuration.h"
#ifdef LIBMESH_HAVE_PETSC
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_shell_matrix.h"

#include "libmesh_cppunit.h"
#include "test_comm.h"

using namespace libMesh;

class ShellMatrixTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(ShellMatrixTest);

  CPPUNIT_TEST(testPetscShell);
  CPPUNIT_TEST(testSparseShell_petsc);

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
    _comm = TestCommWorld;
  }

  void tearDown() {}

  void testPetscShell()
  {
    LOG_UNIT_TEST;

    // Depending on the config: Use different vector-types here!
    PetscLinearSolver<Number> linear_solver(*_comm);
    linear_solver.init();
    NoPcSolverConfiguration solver_conf(linear_solver);
    linear_solver.set_solver_configuration(solver_conf);
    KSPSetInitialGuessNonzero(linear_solver.ksp(), PETSC_TRUE);

    unsigned int maxits=20;
    Real tol = 1e-3;

    UnityShellMat mat = (*_comm);

    // Depending on the config: Use different vector-types here!
    PetscVector<Number> rhs(*_comm, 20);
    rhs.set(5,2.);
    rhs.close();
    auto solution = rhs.clone();
    solution->set(7,2.);
    solution->close();
    linear_solver.solve(mat, *solution, rhs, tol, maxits);

    for (numeric_index_type m=0; m < rhs.size(); ++m)
      CPPUNIT_ASSERT_EQUAL(solution->el(m), rhs.el(m));
  }

  void testSparseShell_petsc()
  {
    LOG_UNIT_TEST;

    // Depending on the config: Use different vector-types here!
    PetscLinearSolver<Number> linear_solver(*_comm);
    linear_solver.init();
    NoPcSolverConfiguration solver_conf(linear_solver);
    linear_solver.set_solver_configuration(solver_conf);
    KSPSetInitialGuessNonzero(linear_solver.ksp(), PETSC_TRUE);

    PetscMatrix<Number> petsc_mat(*_comm);
    petsc_mat.init(20,20,20,20,1,0,1); // hard-coded dimension.
    set_unity(petsc_mat);
    petsc_mat.close();
    SparseShellMatrix<Number> mat(petsc_mat);

    unsigned int maxits=20;
    Real tol = 1e-3;

    // Depending on the config: Use different vector-types here!
    PetscVector<Number> rhs(*_comm, 20);
    rhs.set(5,2.);
    rhs.close();
    auto solution = rhs.clone();
    solution->set(7,2.);
    solution->close();
    linear_solver.solve(mat, *solution, rhs, tol, maxits);

    for (numeric_index_type m=0; m < rhs.size(); ++m)
      CPPUNIT_ASSERT_EQUAL(solution->el(m), rhs.el(m));

    // Make the matrix singular:
    // a) set an irrelevant element to '0':
    solution->set(7,2.);
    petsc_mat.set(2,2,0.);
    petsc_mat.close();
    auto rval = linear_solver.solve(mat, *solution, rhs, tol, maxits);
    // In this case, we can solve the equation.
    CPPUNIT_ASSERT_EQUAL(rval.second, 0.00);
    // b) set the relevant element to '0':
    petsc_mat.set(5,5,0.);
    petsc_mat.close();
    rval = linear_solver.solve(mat, *solution, rhs, tol, maxits);
    // In this case, we can not solve the equation.
    //   The best solution is the 0-vector and we have an error of 2.0
    CPPUNIT_ASSERT_EQUAL(rval.second, 2.00);

  }

private:
  Parallel::Communicator * _comm;

  void set_unity(SparseMatrix<Number> & M)
  {
    for (libMesh::numeric_index_type n=0; n < M.m(); ++n)
      M.set(n,n,1.);
  }


  class NoPcSolverConfiguration : public libMesh::SolverConfiguration
  {
  public:

    NoPcSolverConfiguration(libMesh::PetscLinearSolver<libMesh::Number> & petsc_linear_solver) :
      _petsc_linear_solver(petsc_linear_solver)
    {
    }

    // Shell-matrices are implemented only for few solvers/preconditioners.
    // The SolverConfiguration overrides input-arguments, so here we force GMRES without precond.
    virtual void configure_solver()
    {
      PetscErrorCode ierr = 0;
      ierr = KSPSetType (_petsc_linear_solver.ksp(), const_cast<KSPType>(KSPGMRES));
      CHKERRABORT(_petsc_linear_solver.comm().get(), ierr);
      ierr = PCSetType (_petsc_linear_solver.pc(), const_cast<PCType>(PCNONE));
      CHKERRABORT(_petsc_linear_solver.comm().get(), ierr);
    }

    // The linear solver object that we are configuring
    libMesh::PetscLinearSolver<libMesh::Number> & _petsc_linear_solver;

  };


  class UnityShellMat : public libMesh::ShellMatrix<libMesh::Number>
  {
  public:
    UnityShellMat(Parallel::Communicator & comm):
      ShellMatrix(comm){}

    virtual ~UnityShellMat() = default;

    virtual libMesh::numeric_index_type m () const override
    {
      return 20;
    }
    virtual libMesh::numeric_index_type n () const override
    {
      return 20;
    }

    virtual void vector_mult(libMesh::NumericVector<libMesh::Number> &dest,
                             const libMesh::NumericVector< libMesh::Number > &arg) const override
    {
      dest = arg;
    }

    virtual void vector_mult_add (libMesh::NumericVector<libMesh::Number> &dest,
                                  const libMesh::NumericVector< libMesh::Number > &arg) const override
    {
      dest += arg;
    }

    virtual void get_diagonal (libMesh::NumericVector<libMesh::Number> &dest) const override
    {
      for (numeric_index_type m=0; m < dest.local_size(); ++m)
        dest.set(m, 1.);
    }

  };

};

CPPUNIT_TEST_SUITE_REGISTRATION(ShellMatrixTest);

#endif // #ifdef LIBMESH_HAVE_PETSC
