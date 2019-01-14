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



// \author Robert Weidlich
// \date Copyright 2012

#ifndef SOLID_SYSTEM_H_
#define SOLID_SYSTEM_H_

#include "libmesh/fem_system.h"

using namespace libMesh;

class SolidSystem: public FEMSystem
{
public:
  // Constructor
  SolidSystem(EquationSystems & es,
              const std::string & name,
              const unsigned int number);

  // System initialization
  virtual void init_data();

  // Context initialization
  virtual void init_context(DiffContext & context);

  // Element residual and jacobian calculations
  virtual bool element_time_derivative(bool request_jacobian,
                                       DiffContext & context);

  // Contributions for adding boundary conditions
  virtual bool side_time_derivative(bool request_jacobian,
                                    DiffContext & context);

  virtual bool eulerian_residual(bool, DiffContext &)
  {
    return false;
  }

  // Simulation parameters
  GetPot args;

  // Custom Identifier
  virtual std::string system_type() const
  {
    return "Solid";
  }

  // override method to update mesh also
  virtual void update();

  // save the undeformed mesh to an auxiliary system
  void save_initial_mesh();

  // variable numbers of primary variables in the current system
  unsigned int var[3];

  // variable numbers of primary variables in auxiliary system (for accessing the
  // undeformed configuration)
  unsigned int undefo_var[3];
};

#endif /* SOLID_SYSTEM_H_ */
