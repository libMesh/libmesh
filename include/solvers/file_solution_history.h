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



#ifndef LIBMESH_FILE_SOLUTION_HISTORY_H
#define LIBMESH_FILE_SOLUTION_HISTORY_H

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/solution_history.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/libmesh.h"
#include "libmesh/diff_system.h"

// C++ includes
#include <list>
#include <memory>

namespace libMesh
{

/**
 * Subclass of Solution History that stores the solutions
 * and other important vectors onto disk.
 *
 * \author Vikram Garg
 * \date 2020
 * \brief Stores past solutions onto disk.
 */
class FileSolutionHistory : public SolutionHistory
{
public:

  /**
   * Constructor, reference to system to be passed by user, set the
   * stored_sols iterator to some initial value
   */
  FileSolutionHistory(DifferentiableSystem & system_);
  /**
   * Destructor
   */
  ~FileSolutionHistory();

  /**
   * Virtual function store which we will be overriding to store timesteps
   */
  virtual void store(bool is_adjoint_solve, Real time) override;

  /**
   * Virtual function retrieve which we will be overriding to retrieve timesteps
   */
  virtual void retrieve(bool is_adjoint_solve, Real time) override;

  /**
   * Definition of the clone function needed for the setter function
   */
  virtual std::unique_ptr<SolutionHistory > clone() const override
  {
    return std::make_unique<FileSolutionHistory>(_system);
  }

private:

  // A system reference
  DifferentiableSystem & _system ;

  /**
   * A vector of pointers to adjoint and old adjoint solutions at the last time step.
   * These are used to prevent the zeroing of the adjoint and old adjoint by es::read.
   */
  std::vector< std::unique_ptr<NumericVector<Number>> > dual_solution_copies;
  std::vector< std::unique_ptr<NumericVector<Number>> > old_dual_solution_copies;
};

} // end namespace libMesh

#endif // LIBMESH_FILE_SOLUTION_HISTORY_H
