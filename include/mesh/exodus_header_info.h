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

#ifndef LIBMESH_EXODUS_HEADER_INFO_H
#define LIBMESH_EXODUS_HEADER_INFO_H

// TIMPI includes
#include "timpi/communicator.h"
#include "timpi/parallel_implementation.h"

// C++ includes
#include <vector>

namespace libMesh
{

/**
 * This class is used as both an external data structure for passing
 * around Exodus file header information, and for storing information
 * internally in ExodusII_IO_Helper. It can't really be defined in
 * either of those existing headers, so it has its own header...
 */
class ExodusHeaderInfo
{
public:
  // All special functions can be defaulted for this simple class.
  ExodusHeaderInfo () = default;
  ExodusHeaderInfo (ExodusHeaderInfo &&) = default;
  ExodusHeaderInfo (const ExodusHeaderInfo &) = default;
  ExodusHeaderInfo & operator= (const ExodusHeaderInfo &) = default;
  ExodusHeaderInfo & operator= (ExodusHeaderInfo &&) = default;
  ~ExodusHeaderInfo() = default;

  /**
   * Broadcasts data from processor 0 to other procs using the
   * provided Communicator.
   */
  void broadcast(const Parallel::Communicator & comm)
  {
    // broadcast vector<char> separately
    comm.broadcast(title);

    // Pack individual integers into vector
    std::vector<int> buffer =
      {num_dim, num_elem, num_elem_blk, num_node_sets,
       num_side_sets, num_edge_blk, num_edge};

    // broadcast integers
    comm.broadcast(buffer);

    // unpack
    num_dim       = buffer[0];
    num_elem      = buffer[1];
    num_elem_blk  = buffer[2];
    num_node_sets = buffer[3];
    num_side_sets = buffer[4];
    num_edge_blk  = buffer[5];
    num_edge      = buffer[6];
  }

  std::vector<char> title;
  int num_dim;
  int num_nodes;
  int num_elem;
  int num_elem_blk;
  int num_node_sets;
  int num_side_sets;
  int num_edge_blk;
  int num_edge;
};

} // namespace libMesh

#endif
