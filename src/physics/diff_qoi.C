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


#include "libmesh/diff_qoi.h"

namespace libMesh
{

DifferentiableQoI::DifferentiableQoI () :
  assemble_qoi_sides(false),
  assemble_qoi_internal_sides(false),
  assemble_qoi_elements(true)
{
}

void DifferentiableQoI::thread_join( std::vector<Number> & qoi,
                                     const std::vector<Number> & other_qoi,
                                     const QoISet &)
{
  for (std::size_t i=0; i != qoi.size(); ++i)
    qoi[i] += other_qoi[i];
}

void DifferentiableQoI::parallel_op(const Parallel::Communicator & communicator,
                                    std::vector<Number> & sys_qoi,
                                    std::vector<Number> & local_qoi,
                                    const QoISet &)
{
  // Sum everything into local_qoi
  communicator.sum(local_qoi);

  // Now put into system qoi
  sys_qoi = local_qoi;
}

} // namespace libMesh
