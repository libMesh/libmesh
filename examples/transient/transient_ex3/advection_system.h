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

#ifndef ADVECTION_SYSTEM_H
#define ADVECTION_SYSTEM_H

// Application includes
#include "claw_system.h" // base class

namespace libMesh
{

/**
 * This class extends ClawSystem to implement pure advection
 * in 2D.
 *
 * @author David J. Knezevic, 2012
 * @author John W. Peterson, 2025 (modernization and libmesh example)
 */
class AdvectionSystem : public ClawSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  AdvectionSystem (EquationSystems& es,
                   const std::string& name,
                   const unsigned int number);

  /**
   * Destructor. Defaulted out-of-line.
   */
  virtual ~AdvectionSystem ();

  /**
   * The type of system.
   */
  typedef AdvectionSystem sys_type;

  /**
   * @return a reference to ourself... I'm not sure what the point of
   * this is. The LinearImplicitSystem base class has a similar API,
   * but it's not virtual, so we aren't overriding it by doing
   * this. This could probably be removed.
   */
  sys_type & system () { return *this; }

  /**
   * The type of the parent.
   */
  typedef ClawSystem Parent;

  /**
   * Right-hand side assembly. This function is called from the
   * base class function, ClawSystem::solve_conservation_law().
   * The input vector is either the initial condition (at the
   * first timestep) or the solution vector from the previous
   * timestep.
   */
  virtual void assemble_claw_rhs(NumericVector<Number> & q) override;

  /**
   * Initialize the system (e.g. add variables)
   */
  virtual void init_data() override;

  /**
   * Read in data from input file.
   * Calls the base class version of this function and then reads in
   * some additional parameters.
   */
  virtual void process_parameters_file (const std::string & parameters_filename) override;

  /**
   * Print out some info about the system's configuration.
   */
  virtual void print_info() override;

protected:

  /**
   * Computes the vectors "uq" and "vq" where (u,v) are the (given,
   * constant) advective velocity coefficients, and q is either the
   * initial condition or the previous solution vector.
   */
  void update_Fh(NumericVector<Number> & q);

  /**
   * The advective flux vectors.  The outer std::vector length
   * corresponds to the number of space dimensions in the problem.
   */
  std::vector<std::unique_ptr<NumericVector<Number>>> _Fh;

  /**
   * Variable number for q1
   */
  unsigned int _q1_var;

  /**
   * The (assumed constant) advection velocity
   */
  Point _u;

  /**
   * Store the FE Order and family specified by the input file.
   * Defaults to (CONSTANT, MONOMIAL).
   */
  Order _fe_order;
  FEFamily _fe_family;
};

} // namespace libMesh

#endif
