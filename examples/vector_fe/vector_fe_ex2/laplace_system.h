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
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/dirichlet_boundaries.h"

#include "solution_function.h"

using namespace libMesh;

#ifndef LAPLACE_SYSTEM_H
#define LAPLACE_SYSTEM_H

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class LaplaceSystem : public FEMSystem
{
public:
  // Constructor
  LaplaceSystem(EquationSystems & es,
                const std::string & name,
                const unsigned int number);

  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context(DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context);

  // Constraint part
  // virtual bool side_constraint (bool request_jacobian,
  //                               DiffContext & context);

protected:
  // Indices for each variable;
  unsigned int u_var;

  void init_dirichlet_bcs();

  // Returns the value of a forcing function at point p.
  RealGradient forcing(const Point & p);

  LaplaceExactSolution exact_solution;
};

#endif // LAPLACE_SYSTEM_H
