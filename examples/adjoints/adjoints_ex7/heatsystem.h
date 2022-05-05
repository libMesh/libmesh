// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/enum_fe_family.h"
#include "libmesh/fem_system.h"
#include "libmesh/parameter_pointer.h"
#include "libmesh/parameter_vector.h"

using namespace libMesh;

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class HeatSystem : public FEMSystem
{
public:
  // Constructor
  HeatSystem(EquationSystems & es,
             const std::string & name_in,
             const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
      _k(1.0),
      _fe_family("LAGRANGE"),
      _fe_order(1),
      _analytic_jacobians(true)
  { this->init_qois(2); QoI_time_instant.resize(2);}

  Real & k() { return _k; }
  bool & analytic_jacobians() { return _analytic_jacobians; }

  // A vector to specify which QoIs are instantaneous and what time instant they are evaluated at
  std::vector<Real> QoI_time_instant;

protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context);

  // Library evaluation of QoIs
  virtual void element_qoi (DiffContext & context, const QoISet & qois);

  // RHS for adjoint problem
  virtual void element_qoi_derivative (DiffContext & context,
                                       const QoISet & /* qois */);

  // Parameters associated with the system
  std::vector<Number> parameters;

  // The ParameterVector object that will contain pointers to
  // the system parameters
  ParameterVector parameter_vector;

  // The parameters to solve for
  Real _k;

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;

  // Use averaged model
  bool _averaged_model;

};
