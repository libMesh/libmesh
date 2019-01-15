// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// DiffSystem framework files
#include "libmesh/fem_system.h"

using namespace libMesh;

// The Navier-Stokes system class.
// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class NavierSystem : public FEMSystem
{
public:
  // Constructor
  NavierSystem(EquationSystems & es,
               const std::string & name_in,
               const unsigned int number_in)
    : FEMSystem(es, name_in, number_in), Reynolds(1.), application(0) {}

  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context(DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context);

  // Constraint parts
  virtual bool element_constraint (bool request_jacobian,
                                   DiffContext & context);
  virtual bool side_constraint (bool request_jacobian,
                                DiffContext & context);

  // Mass matrix part
  virtual bool mass_residual (bool request_jacobian,
                              DiffContext & context);

  // Postprocessed output
  virtual void postprocess ();

  // Indices for each variable;
  unsigned int p_var, u_var, v_var, w_var;

  // The Reynolds number to solve for
  Real Reynolds;

  // The application number controls what boundary conditions and/or
  // forcing functions are applied.  Current options are:
  // 0 - discontinuous lid velocity driven cavity
  // 1 - homogeneous Dirichlet BC with smooth forcing
  unsigned int application;

  // Returns the value of a forcing function at point p.  This value
  // depends on which application is being used.
  Point forcing(const Point & p);
};
