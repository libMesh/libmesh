#include <libmesh/dof_map.h>
#include <libmesh/fem_system.h>

using namespace libMesh;

template<typename TimeSolverType>
class TimeSolverTestImplementation
{
protected:

  // Any specialized initialization that's needed for the test
  virtual void aux_time_solver_init( TimeSolverType & /*time_solver*/ ){}

  // Implementation for solving ODE of SystemType
  // Note this test assumes that the time integrator gets the *exact* solution
  // to within floating point tolerance.
  template<typename SystemType>
  void run_test_with_exact_soln(Real deltat, unsigned int n_timesteps)
  {
    Mesh mesh(*TestCommWorld);
    MeshTools::Generation::build_point(mesh);
    EquationSystems es(mesh);
    SystemType & system = es.add_system<SystemType>("ScalarSystem");

    system.time_solver.reset(new TimeSolverType(system));

    es.init();

    DiffSolver & solver = *(system.time_solver->diff_solver().get());
    solver.relative_step_tolerance = std::numeric_limits<Real>::epsilon()*10;
    solver.relative_residual_tolerance = std::numeric_limits<Real>::epsilon()*10;
    solver.absolute_residual_tolerance = std::numeric_limits<Real>::epsilon()*10;

    system.deltat = deltat;

    TimeSolverType * time_solver = cast_ptr<TimeSolverType *>(system.time_solver.get());
    this->aux_time_solver_init(*time_solver);

    // We're going to want to check our solution, and when we run
    // "make check" with LIBMESH_RUN='mpirun -np N" for N>1 then we'll
    // need to avoid checking on the processors that are just
    // twiddling their thumbs, not owning our mesh point.
    std::vector<dof_id_type> solution_index;
    solution_index.push_back(0);
    const bool has_solution = system.get_dof_map().all_semilocal_indices(solution_index);

    for (unsigned int t_step=0; t_step != n_timesteps; ++t_step)
      {
        system.solve();
        system.time_solver->advance_timestep();

        if (has_solution)
          {
            // Use relative error for comparison, so "exact" is 0
            Number exact_soln = system.u(system.time);
            Real rel_error =  std::abs((exact_soln - (*system.solution)(0))/exact_soln);

            CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0,
                                          rel_error,
                                          std::numeric_limits<Real>::epsilon()*10 );
          }
      }
  }

};

//! FEMSystem-based class for testing of TimeSolvers using first order SCALARs
/**
 *  We're assuming the ODEs are only dependent on time, so no Jacobian
 *  functions are needed, just F and M.
 */
class FirstOrderScalarSystemBase : public FEMSystem
{
public:
  FirstOrderScalarSystemBase(EquationSystems & es,
                             const std::string & name_in,
                             const unsigned int number_in)
    : FEMSystem(es, name_in, number_in)
  {}


  //! Value of F(u)
  virtual Number F( FEMContext & context, unsigned int qp ) =0;

  //! Value of M(u).
  virtual Number M( FEMContext & context, unsigned int qp ) =0;

  //! Exact solution as a function of time t.
  virtual Number u( Real t ) =0;

  virtual void init_data ()
  {
    _u_var = this->add_variable ("u", FIRST, LAGRANGE);
    this->time_evolving(_u_var,1);
    FEMSystem::init_data();
  }

  //! Note the nonlinear residual is F(u)-M(u)\dot{u}
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context) libmesh_override
  {
    FEMContext & c = cast_ref<FEMContext &>(context);
    DenseSubVector<Number> & Fu = c.get_elem_residual(_u_var);

    const unsigned int n_u_dofs = c.get_dof_indices(_u_var).size();
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        Number Fval = this->F(c,qp);

        for (unsigned int i=0; i != n_u_dofs; i++)
          Fu(i) += Fval;
      }

    return request_jacobian;
  }

  //! Note the nonlinear residual is F(u)-M(u)\dot{u}
  virtual bool mass_residual (bool request_jacobian,
                              DiffContext & context) libmesh_override
  {
    FEMContext & c = cast_ref<FEMContext &>(context);
    DenseSubVector<Number> & Fu = c.get_elem_residual(_u_var);
    DenseSubMatrix<Number> & Kuu = c.get_elem_jacobian(_u_var, _u_var);

    const unsigned int n_u_dofs = c.get_dof_indices(_u_var).size();
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number udot;
        c.interior_rate(_u_var,qp,udot);

        Number Mval = this->M(c,qp);

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) -= Mval*udot;

            if (request_jacobian)
              for (unsigned int j=0; j != n_u_dofs; j++)
                Kuu(i,j) -= Mval*context.get_elem_solution_rate_derivative();
          }
      }

    return request_jacobian;
  }

protected:
  unsigned int _u_var;
};

//! FEMSystem-based class for testing of TimeSolvers using second order SCALARs
/**
 *  This is for solving second order systems using second order time solvers.
 *  We're assuming the ODEs are only dependent on time, so no Jacobian
 *  functions are needed, just F, C, and M.
 */
class SecondOrderScalarSystemSecondOrderTimeSolverBase : public FirstOrderScalarSystemBase
{
public:
  SecondOrderScalarSystemSecondOrderTimeSolverBase(EquationSystems & es,
                                                   const std::string & name_in,
                                                   const unsigned int number_in)
    : FirstOrderScalarSystemBase(es, name_in, number_in)
  {}

  virtual void init_data ()
  {
    _u_var = this->add_variable ("u", FIRST, LAGRANGE);
    this->time_evolving(_u_var,2);
    FEMSystem::init_data();
  }

  //! Value of C(u).
  virtual Number C( FEMContext & context, unsigned int qp ) =0;

  //! Note the nonlinear residual is M(u)\ddot{u} + C(u)\dot{u} + F(u)
  virtual bool damping_residual (bool request_jacobian,
                                 DiffContext & context) libmesh_override
  {
    FEMContext & c = cast_ref<FEMContext &>(context);
    DenseSubVector<Number> & Fu = c.get_elem_residual(_u_var);
    DenseSubMatrix<Number> & Kuu = c.get_elem_jacobian(_u_var, _u_var);

    const unsigned int n_u_dofs = c.get_dof_indices(_u_var).size();
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number udot;
        c.interior_rate(_u_var,qp,udot);

        Number Cval = this->C(c,qp);

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += Cval*udot;

            if (request_jacobian)
              for (unsigned int j=0; j != n_u_dofs; j++)
                Kuu(i,j) += Cval*context.get_elem_solution_rate_derivative();
          }
      }

    return request_jacobian;
  }

  virtual bool mass_residual (bool request_jacobian,
                              DiffContext & context) libmesh_override
  {
    FEMContext & c = cast_ref<FEMContext &>(context);
    DenseSubVector<Number> & Fu = c.get_elem_residual(_u_var);
    DenseSubMatrix<Number> & Kuu = c.get_elem_jacobian(_u_var, _u_var);

    const unsigned int n_u_dofs = c.get_dof_indices(_u_var).size();
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number uddot;
        c.interior_accel(_u_var,qp,uddot);

        Number Mval = this->M(c,qp);

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += Mval*uddot;

            if (request_jacobian)
              for (unsigned int j=0; j != n_u_dofs; j++)
                Kuu(i,j) += Mval*context.get_elem_solution_accel_derivative();
          }
      }

    return request_jacobian;
  }
};


//! FEMSystem-based class for testing of TimeSolvers using second order SCALARs
/**
 *  This is for solving second order systems using *first* order *or* second order
 *  time solvers. We're assuming the ODEs are only dependent on time, so no Jacobian
 *  functions are needed, just F, C, and M.
 */
class SecondOrderScalarSystemFirstOrderTimeSolverBase : public SecondOrderScalarSystemSecondOrderTimeSolverBase
{
public:
  SecondOrderScalarSystemFirstOrderTimeSolverBase(EquationSystems & es,
                                                  const std::string & name_in,
                                                  const unsigned int number_in)
    : SecondOrderScalarSystemSecondOrderTimeSolverBase(es, name_in, number_in)
  {}

  //! Note the nonlinear residual is M(u)\dot{v} + C(u)v + F(u)
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context) libmesh_override
  {
    FEMContext & c = cast_ref<FEMContext &>(context);

    unsigned int v_var = this->get_second_order_dot_var(_u_var);

    DenseSubVector<Number> & Fv = c.get_elem_residual(v_var);

    const unsigned int n_u_dofs = c.get_dof_indices(_u_var).size();
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        Number Fval = this->F(c,qp);

        for (unsigned int i=0; i != n_u_dofs; i++)
          Fv(i) += Fval;
      }

    return request_jacobian;
  }

  //! Note the nonlinear residual is M(u)\dot{v} + C(u)v + F(u)
  virtual bool damping_residual (bool request_jacobian,
                                 DiffContext & context) libmesh_override
  {
    FEMContext & c = cast_ref<FEMContext &>(context);

    unsigned int v_var = this->get_second_order_dot_var(_u_var);

    DenseSubVector<Number> & Fv = c.get_elem_residual(v_var);

    DenseSubMatrix<Number> & Kvv = c.get_elem_jacobian(v_var, v_var);

    const unsigned int n_u_dofs = c.get_dof_indices(_u_var).size();
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number udot;
        c.interior_rate(v_var,qp,udot);

        Number Cval = this->C(c,qp);

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fv(i) += Cval*udot;

            if (request_jacobian)
              for (unsigned int j=0; j != n_u_dofs; j++)
                Kvv(i,j) += Cval*context.get_elem_solution_rate_derivative();
          }
      }

    return request_jacobian;
  }

  virtual bool mass_residual (bool request_jacobian,
                              DiffContext & context) libmesh_override
  {
    FEMContext & c = cast_ref<FEMContext &>(context);

    unsigned int v_var = this->get_second_order_dot_var(_u_var);

    DenseSubVector<Number> & Fv = c.get_elem_residual(v_var);
    DenseSubMatrix<Number> & Kvv = c.get_elem_jacobian(v_var, v_var);

    const unsigned int n_u_dofs = c.get_dof_indices(_u_var).size();
    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number uddot;
        c.interior_accel(v_var,qp,uddot);

        Number Mval = this->M(c,qp);

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fv(i) += Mval*uddot;

            if (request_jacobian)
              for (unsigned int j=0; j != n_u_dofs; j++)
                Kvv(i,j) += Mval*context.get_elem_solution_accel_derivative();
          }
      }

    return request_jacobian;
  }
};
