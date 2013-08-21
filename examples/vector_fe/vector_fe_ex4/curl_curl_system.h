/* The libMesh Finite Element Library. */
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
#include "libmesh/fem_system.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/dirichlet_boundaries.h"

#include "curl_curl_exact_solution.h"

using namespace libMesh;

#ifndef __curl_curl_system_h__
#define __curl_curl_system_h__

// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class CurlCurlSystem : public FEMSystem
{

public:

  // Constructor
  CurlCurlSystem( EquationSystems& es,
		 const std::string& name,
		 const unsigned int number );

  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context(DiffContext &context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext& context);

  virtual bool side_time_derivative(bool request_jacobian,
				    DiffContext& context);

protected:
  // Indices for each variable;
  unsigned int u_var;

  void init_dirichlet_bcs();

  // Returns the value of a forcing function at point p.
  RealGradient forcing(const Point& p);

  // Returns the value of the exact solution for this model problem.
  RealGradient exact_solution( const Point& p );

  CurlCurlExactSolution soln;

};

#endif //__curl_curl_system_h__
