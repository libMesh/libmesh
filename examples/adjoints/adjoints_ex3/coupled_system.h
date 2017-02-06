// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef COUPLED_SYSTEM_H
#define COUPLED_SYSTEM_H

// DiffSystem framework files
#include "libmesh/fem_function_base.h"
#include "libmesh/fem_system.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/parameter_vector.h"

using namespace libMesh;

// The Navier-Stokes system class.
// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class CoupledSystem : public FEMSystem
{
public:
  // Constructor
  CoupledSystem(EquationSystems & es,
                const std::string & name_in,
                const unsigned int number_in)
    : FEMSystem(es, name_in, number_in), Peclet(1.) {qoi.resize(1);}

  // Function to get computed QoI values

  Number & get_QoI_value()
  {
    return computed_QoI;
  }

  Number & get_parameter_value(unsigned int parameter_index)
  {
    return parameters[parameter_index];
  }

  ParameterVector & get_parameter_vector()
  {
    parameter_vector.resize(parameters.size());
    for (std::size_t i = 0; i != parameters.size(); ++i)
      parameter_vector[i] = &parameters[i];

    return parameter_vector;
  }

  Real & get_Pe()
  {
    return Peclet;
  }

protected:

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

  // Postprocessed output
  virtual void postprocess ();

  // Parameters associated with the system
  std::vector<Number> parameters;

  // Indices for each variable;
  unsigned int p_var, u_var, v_var, C_var;

  // The ParameterVector object that will contain pointers to
  // the system parameters
  ParameterVector parameter_vector;

  // The Peclet number for the species transport
  Real Peclet;

  // The functionals to be computed as QoIs
  Number computed_QoI;
};


class CoupledFEMFunctionsx : public FEMFunctionBase<Number>
{
public:
  // Constructor
  CoupledFEMFunctionsx(System & /* sys */,
                       unsigned int var_number)
  {var = var_number;}

  // Destructor
  virtual ~CoupledFEMFunctionsx () {}

  virtual UniquePtr<FEMFunctionBase<Number> > clone () const
  {
    return UniquePtr<FEMFunctionBase<Number> >(new CoupledFEMFunctionsx(*this));
  }

  virtual void operator() (const FEMContext &,
                           const Point &,
                           const Real,
                           DenseVector<Number> &)
  { libmesh_not_implemented(); }

  virtual Number operator() (const FEMContext &,
                             const Point & p,
                             const Real time = 0.);

private:
  unsigned int var;
};


class CoupledFEMFunctionsy : public FEMFunctionBase<Number>
{
public:
  // Constructor
  CoupledFEMFunctionsy(System & /* sys */,
                       unsigned int var_number)
  { var = var_number; }

  // Destructor
  virtual ~CoupledFEMFunctionsy () {}

  virtual UniquePtr<FEMFunctionBase<Number> > clone () const
  {
    return UniquePtr<FEMFunctionBase<Number> >(new CoupledFEMFunctionsy(*this));
  }

  virtual void operator() (const FEMContext &,
                           const Point &,
                           const Real,
                           DenseVector<Number> &)
  { libmesh_not_implemented(); }

  virtual Number operator() (const FEMContext &,
                             const Point & p,
                             const Real time = 0.);

private:
  unsigned int var;
};

#endif // COUPLED_SYSTEM_H
