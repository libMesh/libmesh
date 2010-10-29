// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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


#ifndef __sparsity_pattern_h__
#define __sparsity_pattern_h__

// C++ includes
#include <vector>

// Local Includes
#include "threads_allocators.h"
#include "elem_range.h"

namespace libMesh
{

  // Forward declaractions
  class MeshBase;
  class DofMap;
  class CouplingMatrix;
  
// ------------------------------------------------------------
// Sparsity Pattern

/**
 * This defines the sparsity pattern, or graph, of a sparse matrix.
 * The format is quite simple -- the global indices of the nonzero entries
 * in each row are packed into a vector.  The global indices (i,j) of the 
 * nth nonzero entry of row i are given by j = sparsity_pattern[i][n];
 */
namespace SparsityPattern // use a namespace so member classes can be forward-declared.
{
  typedef std::vector<unsigned int, Threads::scalable_allocator<unsigned int> > Row;
  class Graph : public std::vector<Row> {};
  
  /**
   * Splices the two sorted ranges [begin,middle) and [middle,end)
   * into one sorted range [begin,end).  This method is much like
   * std::inplace_merge except it assumes the intersection
   * of the two sorted ranges is empty and that any element in
   * each range occurs only once in that range.  Additionally,
   * this sort occurs in-place, while std::inplace_merge may
   * use a temporary buffer.
   */ 
  template<typename BidirectionalIterator>
  static void sort_row (const BidirectionalIterator begin,
			BidirectionalIterator       middle,
			const BidirectionalIterator end);

  /**
   * This helper class can be called on multiple threads to compute 
   * the sparsity pattern (or graph) of the sparse matrix resulting
   * from the discretization.  This pattern may be used directly by
   * a particular sparse matrix format (e.g. \p LaspackMatrix)
   * or indirectly (e.g. \p PetscMatrix).  In the latter case the
   * number of nonzeros per row of the matrix is needed for efficient 
   * preallocation.  In this case it suffices to provide estimate
   * (but bounding) values, and in this case the threaded method can
   * take some short-cuts for efficiency.
   */
  class Build
  {
  private:
    const MeshBase &mesh;
    const DofMap &dof_map;
    const CouplingMatrix *dof_coupling;
    const bool implicit_neighbor_dofs;
    const bool need_full_sparsity_pattern;

  public:

    SparsityPattern::Graph sparsity_pattern;
    std::vector<unsigned int> n_nz;
    std::vector<unsigned int> n_oz;
    
    Build (const MeshBase &mesh_in,
	   const DofMap &dof_map_in,
	   const CouplingMatrix *dof_coupling_in,
	   const bool implicit_neighbor_dofs_in,
	   const bool need_full_sparsity_pattern_in) :
      mesh(mesh_in),
      dof_map(dof_map_in),
      dof_coupling(dof_coupling_in),
      implicit_neighbor_dofs(implicit_neighbor_dofs_in),
      need_full_sparsity_pattern(need_full_sparsity_pattern_in)
    {}

    Build (Build &other, Threads::split) :
      mesh(other.mesh),
      dof_map(other.dof_map),
      dof_coupling(other.dof_coupling),
      implicit_neighbor_dofs(other.implicit_neighbor_dofs),
      need_full_sparsity_pattern(other.need_full_sparsity_pattern)
    {}

    void operator()(const ConstElemRange &range);

    void join (const Build &other);
  };

#if defined(__GNUC__) && (__GNUC__ < 4) && !defined(__INTEL_COMPILER)
  /**
   * Dummy function that does nothing but can be used to prohibit
   * compiler optimization in some situations where some compilers
   * have optimization bugs.
   */
  void _dummy_function(void);
#endif
  
}



// ------------------------------------------------------------
// SparsityPattern inline member functions
template<typename BidirectionalIterator>
inline
void SparsityPattern::sort_row (const BidirectionalIterator begin,
				BidirectionalIterator       middle,
				const BidirectionalIterator end)
{
  if ((begin == middle) || (middle == end)) return;

  libmesh_assert (std::distance (begin,  middle) > 0);
  libmesh_assert (std::distance (middle, end)    > 0);
  libmesh_assert (std::unique (begin,  middle) == middle);
  libmesh_assert (std::unique (middle, end)    == end);
  
  while (middle != end)
    {
      BidirectionalIterator
	b = middle,
	a = b-1;

      // Bubble-sort the middle value downward
      while (!(*a < *b)) // *a & *b are less-than comparable, so use <
	{
	  std::swap (*a, *b);

#if defined(__GNUC__) && (__GNUC__ < 4) && !defined(__INTEL_COMPILER)
	  /* Prohibit optimization at this point since gcc 3.3.5 seems
	     to have a bug.  */
	  SparsityPattern::_dummy_function();
#endif

	  if (a == begin) break;
	  
	  b=a;
	  --a;
	}
      
      ++middle;
    }

  // Assure the algorithm worked if we are in DEBUG mode
#ifdef DEBUG
  {
    // SGI STL extension!
    // libmesh_assert (std::is_sorted(begin,end));
    
    BidirectionalIterator
      prev  = begin,
      first = begin;

    for (++first; first != end; prev=first, ++first)
      if (*first < *prev)
	libmesh_assert(false);
  }  
#endif

  // Make sure the two ranges did not contain any common elements
  libmesh_assert (std::unique (begin, end) == end);
}

} // namespace libMesh

#endif // __sparsity_pattern_h__
