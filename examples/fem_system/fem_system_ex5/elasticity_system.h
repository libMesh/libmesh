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

// boundary IDs
extern boundary_id_type
  boundary_id_min_x, boundary_id_max_x,
  boundary_id_min_y, boundary_id_max_y,
  boundary_id_min_z, boundary_id_max_z;
const boundary_id_type node_boundary_id = 10;
const boundary_id_type edge_boundary_id = 20;
static const boundary_id_type & traction_boundary_id = boundary_id_max_x;
const boundary_id_type pressure_boundary_id = 30;
const boundary_id_type fixed_u_boundary_id = 40;
const boundary_id_type fixed_v_boundary_id = 50;

// The Navier-Stokes system class.
// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class ElasticitySystem : public FEMSystem
{
public:
  // Constructor
  ElasticitySystem(EquationSystems & es,
                   const std::string & name_in,
                   const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
      _dim(3),
      _fe_type(),
      _rho(1.0)
  {}

  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context(DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context);

  virtual bool side_time_derivative (bool request_jacobian,
                                     DiffContext & context);

  // Mass matrix part
  virtual bool mass_residual (bool request_jacobian,
                              DiffContext & context);

  void set_dim(unsigned int dim) { _dim = dim; }

  void set_fe_type(const FEType & fe_type) { _fe_type = fe_type; }

private:

  unsigned int _dim;

  FEType _fe_type;

  // Indices for each variable;
  unsigned int _u_var, _v_var, _w_var;

  Real _rho;

  Real kronecker_delta(unsigned int i, unsigned int j)
  {
    return i == j ? 1. : 0.;
  }

  Real elasticity_tensor(unsigned int i, unsigned int j, unsigned int k, unsigned int l);
};
