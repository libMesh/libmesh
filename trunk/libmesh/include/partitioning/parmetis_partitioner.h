// $Id: parmetis_partitioner.h,v 1.2 2004-01-03 15:37:42 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __parmetis_partitioner_h__
#define __parmetis_partitioner_h__

// C++ Includes   -----------------------------------
#include <string>
#include <vector>

// Local Includes -----------------------------------
#include "partitioner.h"



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
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void partition (MeshBase& mesh,
			  const unsigned int n = libMesh::n_processors());

  /**
   * Parmetis can handle dynamically repartitioning a mesh such
   * that the redistribution costs are minimized.  This method
   * takes a previously partitioned domain (which may have
   * then been adaptively refined) and repartitions it.
   */
  virtual void repartition (MeshBase& mesh,
			    const unsigned int n = libMesh::n_processors());

private:

// These methods & data only need to be available if the
// ParMETIS library is available.
#ifdef HAVE_PARMETIS

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
   * Maps active element ids into a contiguous range, as needed by ParMETIS.
   */
  std::vector<unsigned int> _forward_map;
  unsigned int _first_local_elem;
  
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




#endif // #define __parmetis_partitioner_h__
