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

// libMesh includes
#include "libmesh/enum_fe_family.h"
#include "libmesh/fdm_gradient.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/fem_system.h"
#include "libmesh/libmesh_common.h"

// C++ includes
#include <map>
#include <memory>


// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class HilbertSystem : public libMesh::FEMSystem
{
public:
  // Constructor
  HilbertSystem(libMesh::EquationSystems & es,
                const std::string & name,
                const unsigned int number)
  : libMesh::FEMSystem(es, name, number),
    input_system(nullptr),
    _fe_family("LAGRANGE"),
    _fe_order(1),
    _hilbert_order(0),
    _subdomains_list() {}

  // Default destructor
  ~HilbertSystem();

  std::string & fe_family() { return _fe_family; }
  unsigned int & fe_order() { return _fe_order; }
  std::set<libMesh::subdomain_id_type> & subdomains_list() { return _subdomains_list; }

  unsigned int & hilbert_order() { return _hilbert_order; }
  void set_fdm_eps(libMesh::Real eps) {
    _fdm_eps = eps;
    if (_goal_func.get())
      _goal_grad = std::make_unique<libMesh::FDMGradient<libMesh::Gradient>>(*_goal_func, _fdm_eps);
  }

  void set_goal_func(libMesh::FEMFunctionBase<libMesh::Number> & goal)
  {
    _goal_func = goal.clone();
    _goal_grad = std::make_unique<libMesh::FDMGradient<libMesh::Gradient>>(*_goal_func, _fdm_eps);
  }

  // We want to be able to project functions based on *other* systems'
  // values.  For that we need not only a FEMFunction but also a
  // reference to the system where it applies and a separate context
  // object (or multiple separate context objects, in the threaded
  // case) for that system.
  libMesh::System * input_system;

protected:
  std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number> > _goal_func;

  std::map<libMesh::FEMContext *, std::unique_ptr<libMesh::FEMContext>>
    input_contexts;

  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (libMesh::DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        libMesh::DiffContext & context);

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // The Hilbert order our subclass will project with
  unsigned int _hilbert_order;

  // The function we will call to finite difference our goal
  // function
  std::unique_ptr<libMesh::FDMGradient<libMesh::Gradient>> _goal_grad;

  // The perturbation we will use when finite differencing our goal
  // function
  libMesh::Real _fdm_eps;

  // Which subdomains to integrate on (all subdomains, if empty())
  std::set<libMesh::subdomain_id_type> _subdomains_list;
};
