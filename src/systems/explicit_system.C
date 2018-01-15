// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// C++ includes

// Local includes
#include "libmesh/explicit_system.h"
#include "libmesh/numeric_vector.h"

namespace libMesh
{


// ------------------------------------------------------------
// ExplicitSystem implementation
ExplicitSystem::ExplicitSystem (EquationSystems & es,
                                const std::string & name_in,
                                const unsigned int number_in) :
  Parent (es, name_in, number_in),
  rhs(libmesh_nullptr)

{
  // Add the system RHS.
  this->add_system_rhs ();
}



void ExplicitSystem::clear ()
{
  // Clear the parent data
  Parent::clear();

  // Restore us to a "basic" state
  this->add_system_rhs ();
}



void ExplicitSystem::assemble_qoi (const QoISet & qoi_indices)
{
  // The user quantity of interest assembly gets to expect to
  // accumulate on initially zero values
  for (std::size_t i=0; i != qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      qoi[i] = 0;

  Parent::assemble_qoi (qoi_indices);
}



void ExplicitSystem::assemble_qoi_derivative (const QoISet & qoi_indices,
                                              bool include_liftfunc,
                                              bool apply_constraints)
{
  // The user quantity of interest derivative assembly gets to expect
  // to accumulate on initially zero vectors
  for (std::size_t i=0; i != qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      this->add_adjoint_rhs(i).zero();

  Parent::assemble_qoi_derivative (qoi_indices, include_liftfunc,
                                   apply_constraints);
}



void ExplicitSystem::solve ()
{
  // Assemble the linear system
  this->assemble ();

  // Update the system after the solve
  this->update();
}



void ExplicitSystem::add_system_rhs ()
{
  // Possible that we cleared the _vectors but
  // forgot to NULL-out the rhs?
  if (this->n_vectors() == 0)
    rhs = libmesh_nullptr;


  // Only need to add the rhs if it isn't there
  // already!
  if (rhs == libmesh_nullptr)
    rhs = &(this->add_vector ("RHS Vector", false));

  libmesh_assert(rhs);
}

} // namespace libMesh
