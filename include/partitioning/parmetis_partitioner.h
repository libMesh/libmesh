// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PARMETIS_PARTITIONER_H
#define LIBMESH_PARMETIS_PARTITIONER_H

// Local Includes -----------------------------------
#include "libmesh/id_types.h"
#include "libmesh/partitioner.h"
#include "libmesh/vectormap.h"

// C++ Includes   -----------------------------------
#include <cstddef>
#include <vector>

namespace libMesh
{



/**
 * The \p ParmetisPartitioner uses the Parmetis graph partitioner
 * to partition the elements.
 */

// ------------------------------------------------------------
// ParmetisPartitioner class definition
class ParmetisPartitioner : public Partitioner
{
 public:

  /**
   * Constructor.
   */
  ParmetisPartitioner () {}

  /**
   * Creates a new partitioner of this type and returns it in
   * an \p AutoPtr.
   */
  virtual AutoPtr<Partitioner> clone () const {
    AutoPtr<Partitioner> cloned_partitioner
      (new ParmetisPartitioner());
    return cloned_partitioner;
  }


protected:

  /**
   * Parmetis can handle dynamically repartitioning a mesh such
   * that the redistribution costs are minimized.  This method
   * takes a previously partitioned domain (which may have
   * then been adaptively refined) and repartitions it.
   */
  virtual void _do_repartition (MeshBase& mesh,
				const unsigned int n);

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase& mesh,
			      const unsigned int n);

private:

// These methods & data only need to be available if the
// ParMETIS library is available.
#ifdef LIBMESH_HAVE_PARMETIS

  /**
   * Initialize data structures.
   */
  void initialize (const MeshBase& mesh, const unsigned int n_sbdmns);

  /**
   * Build the graph.
   */
  void build_graph (const MeshBase& mesh);

  /**
   * Assign the computed partitioning to the mesh.
   */
  void assign_partitioning (MeshBase& mesh);

  /**
   * The number of active elements on each processor.  Note that
   * ParMETIS requires that each processor have some active elements,
   * it will abort if any processor passes a NULL _part array.
   */
  std::vector<dof_id_type> _n_active_elem_on_proc;

  /**
   * Maps active element ids into a contiguous range, as needed by ParMETIS.
   */
  vectormap<dof_id_type, dof_id_type> _global_index_by_pid_map;

  /**
   * Data structures used by ParMETIS to describe the connectivity graph
   * of the mesh.  Consult the ParMETIS documentation.
   */
  std::vector<int>    _vtxdist;
  std::vector<int>    _xadj;
  std::vector<int>    _adjncy;
  std::vector<int>    _part;
  std::vector<float>  _tpwgts;
  std::vector<float>  _ubvec;
  std::vector<int>    _options;
  std::vector<int>    _vwgt;

  int _wgtflag;
  int _ncon;
  int _numflag;
  int _nparts;
  int _edgecut;

#endif
};


} // namespace libMesh



#endif // LIBMESH_PARMETIS_PARTITIONER_H
