/* $Id: naviersystem.h,v 1.2 2006-12-11 23:17:54 roystgnr Exp $ */

/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

// DiffSystem framework files
#include "fem_system.h"

// The Navier-Stokes system class.
// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class NavierSystem : public FEMSystem
{
public:
  // Constructor
  NavierSystem(EquationSystems& es,
               const std::string& name,
               const unsigned int number)
  : FEMSystem(es, name, number), Reynolds(1.) {}

  // System initialization
  virtual void init_data ();

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian);

  // Constraint parts
  virtual bool element_constraint (bool request_jacobian);
  virtual bool side_constraint (bool request_jacobian);

  // Indices for each variable;
  unsigned int p_var, u_var, v_var, w_var;

  // Finite elements for the velocity and pressure on element interiors
  FEBase *fe_velocity, *fe_pressure;

  // Finite element for the velocity on element sides
  FEBase *fe_side_vel;

  // The Reynolds number to solve for
  Real Reynolds;
};
