/* The libMesh Finite Element Library. */
/* Copyright (C) 2002-2013 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner */

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




#ifndef __fem_parameters_h__
#define __fem_parameters_h__

#include "libmesh/libmesh_common.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/function_base.h"
#include "libmesh/getpot.h"
#include "libmesh/id_types.h"
#include "libmesh/periodic_boundaries.h"

#include <limits>
#include <map>
#include <string>
#include <vector>

// Bring in everything from the libMesh namespace
using namespace libMesh;

class FEMParameters
{
public:
  FEMParameters();

  ~FEMParameters();

  void read(GetPot &input);

  // Parameters applicable to entire EquationSystems:

  unsigned int initial_timestep, n_timesteps;
  bool transient;
  unsigned int deltat_reductions;
  std::string timesolver_core;
  Real end_time, deltat, timesolver_theta,
    timesolver_maxgrowth, timesolver_tolerance,
    timesolver_upper_tolerance, steadystate_tolerance;
  std::vector<FEMNormType> timesolver_norm;

  //   Mesh generation

  unsigned int dimension;
  std::string domaintype, domainfile, elementtype;
  Real elementorder;
  Real domain_xmin, domain_ymin, domain_zmin;
  Real domain_edge_width, domain_edge_length, domain_edge_height;
  unsigned int coarsegridx, coarsegridy, coarsegridz;
  unsigned int coarserefinements, extrarefinements;

  //   Mesh refinement

  unsigned int nelem_target;
  Real global_tolerance;
  Real refine_fraction, coarsen_fraction, coarsen_threshold;
  unsigned int max_adaptivesteps;
  unsigned int initial_adaptivesteps;

  //   Output

  unsigned int write_interval;
  bool write_gmv_error, write_tecplot_error, output_xda, output_xdr,
    output_bz2, output_gz, output_gmv, output_tecplot;

  // Types of Systems to create

  std::vector<std::string> system_types;

  // Parameters applicable to each system:

  //   Boundary and initial conditions

  //std::vector<PeriodicBoundary> periodic_boundaries;

  std::map<subdomain_id_type, FunctionBase<Number> *> initial_conditions;
  std::map<boundary_id_type, FunctionBase<Number> *>  dirichlet_conditions,
    neumann_conditions;
  std::map<boundary_id_type, std::vector<unsigned int> > dirichlet_condition_variables,
    neumann_condition_variables;
  std::map<int, std::map<subdomain_id_type, FunctionBase<Number> *> >
  other_interior_functions;
  std::map<int, std::map<boundary_id_type, FunctionBase<Number> *> >
  other_boundary_functions;

  //   Execution type

  bool run_simulation, run_postprocess;

  //   Approximation type

  std::vector<std::string> fe_family;
  std::vector<unsigned int> fe_order;
  int extra_quadrature_order;

  //   Assembly options

  bool analytic_jacobians;
  Real verify_analytic_jacobians;
  Real numerical_jacobian_h;

  bool print_solution_norms, print_solutions,
    print_residual_norms, print_residuals,
    print_jacobian_norms, print_jacobians,
    print_element_jacobians, print_element_residuals;

  //   Solver options

  bool use_petsc_snes;
  bool time_solver_quiet, solver_quiet, solver_verbose, require_residual_reduction;
  Real min_step_length;
  unsigned int max_linear_iterations, max_nonlinear_iterations;
  Real relative_step_tolerance, relative_residual_tolerance,
    initial_linear_tolerance, minimum_linear_tolerance,
    linear_tolerance_multiplier;

  // Initialization

  unsigned int initial_sobolev_order;
  unsigned int initial_extra_quadrature;

  //   Error indicators

  std::string indicator_type;
  unsigned int sobolev_order;

  //   System-specific parameters file

  std::string system_config_file;
};

#endif // __fem_parameters_h__
