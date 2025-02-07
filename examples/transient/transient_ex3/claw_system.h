// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef CLAW_SYSTEM_H
#define CLAW_SYSTEM_H

// libMesh includes
#include "libmesh/linear_implicit_system.h" // base class

namespace libMesh
{

/**
 * This class encapsulates functionality that allows us to
 * solve conservation laws. Here we focus on 2D, and our approach
 * is to "pre-assemble" all the discretization matrices rather
 * than doing explicit elements loops in each time step.
 *
 * @author David J. Knezevic, 2012
 * @author John W. Peterson, 2025 (modernization and libmesh example)
 */
class ClawSystem : public LinearImplicitSystem
{
public:

  /**
   * Define an enumeration for temporal discretization types
   * ForwardEuler = 1st-order Explicit Euler
   * RK4 = 4th-order explicit Runge-Kutta
   */
  enum TemporalDiscretizationType {ForwardEuler = 0, RK4};

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  ClawSystem (EquationSystems & es,
              const std::string & name,
              const unsigned int number);

  /**
   * Destructor. Defaulted out-of-line.
   */
  virtual ~ClawSystem ();

  /**
   * The type of the parent.
   */
  typedef LinearImplicitSystem Parent;

  /**
   * @returns a string indicating the type of the system.
   */
  virtual std::string system_type () const override;

  /**
   * Convert between string and enum for TemporalDiscretizationType.
   */
  static TemporalDiscretizationType string_to_enum(std::string string_type);
  static std::string enum_to_string(TemporalDiscretizationType enum_type);

  /**
   * Set parameters for this system (e.g. LxF_constant, delta_t) by reading
   * in from file.
   */
  virtual void process_parameters_file (const std::string& parameters_filename);

  /**
   * Solve the conservation law.
   */
  void solve_conservation_law();

  /**
   * Assemble the right-hand side vector. This should
   * be implemented in subclasses
   */
  virtual void assemble_claw_rhs(NumericVector<Number> &) = 0;

  /**
   * Assemble the matrices we need to solve a conservation law.
   */
  void assemble_all_matrices();

  /**
   * Get a reference to one of the discretization matrices.
   * For APIs that take a "dim" argument, this function returns
   * the matrix that corresponds to that direction, i.e. 0==x, 1==y.
   */
  SparseMatrix<Number> & get_mass_matrix();
  SparseMatrix<Number> & get_advection_matrix(unsigned int dim);
  SparseMatrix<Number> & get_avg_matrix(unsigned int dim);
  SparseMatrix<Number> & get_jump_matrix();
  SparseMatrix<Number> & get_boundary_condition_matrix(unsigned int dim);

  /**
   * Set/get the time-step size.
   */
  void set_delta_t(Real delta_t_in);
  Real get_delta_t();

  /**
   * Set/get the number of time-steps.
   */
  void set_n_time_steps(unsigned int n_time_steps_in);
  unsigned int get_n_time_steps();

  /**
   * Set/get the Lax-Friedrichs constant.
   * If not set by the user, the default value is 1.0.
   */
  void set_LxF_constant(Real LxF_constant_in);
  Real get_LxF_constant();

  /**
   * Set/get the temporal discretization type.
   */
  void set_temporal_discretization_type(TemporalDiscretizationType td_in);
  TemporalDiscretizationType get_temporal_discretization_type();

  /**
   * Set/get write_interval.
   */
  void set_write_interval(unsigned int write_interval_in);
  unsigned int get_write_interval();

  /**
   * Set/get the current time in the conservation law solve.
   */
  void set_time(Real time_in);
  Real get_time();

  /**
   * Print out some info about the system's configuration.
   */
  virtual void print_info();

  /**
   * Print discretization matrices to file in MATLAB-readable format.
   */
  void write_out_discretization_matrices();

  /**
   * Initialize the system (e.g. add variables)
   */
  virtual void init_data() override;

private:

  /**
   * Helper functions called by assemble_all_matrices().
   */
  void assemble_mass_matrix();
  void assemble_advection_matrices();
  void assemble_avg_coupling_matrices();
  void assemble_jump_coupling_matrix();
  void assemble_boundary_condition_matrices();

  /**
   * The mass matrix.
   */
  std::unique_ptr<SparseMatrix<Number>> _mass_matrix;

  /**
   * The "advection" matrices.
   */
  std::vector<std::unique_ptr<SparseMatrix<Number>>> _advection_matrices;

  /**
   * The "flux average" element boundary matrices.
   */
  std::vector<std::unique_ptr<SparseMatrix<Number>>> _avg_matrices;

  /**
   * The "jump" element boundary matrix.
   */
  std::unique_ptr<SparseMatrix<Number>> _jump_matrix;

  /**
   * The "boundary condition" matrices. There is one matrix per space
   * dimension in the problem.
   */
  std::vector<std::unique_ptr<SparseMatrix<Number>>> _boundary_condition_matrices;

  /**
   * String that defines the type of temporal discretization
   * that we use.
   */
  TemporalDiscretizationType _temporal_discretization_type;

  /**
   * The time step size
   */
  Real _delta_t;

  /**
   * The number of time steps
   */
  unsigned int _n_time_steps;

  /**
   * The constant C in the Lax-Friedrichs flux.
   */
  Real _LxF_constant;

  /**
   * The time step interval between writing out solutions
   */
  unsigned int _write_interval;

  /**
   * The current time in the conservation law solve
   */
  Real _time;
};

} // namespace libMesh

#endif
