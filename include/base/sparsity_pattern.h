// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_SPARSITY_PATTERN_H
#define LIBMESH_SPARSITY_PATTERN_H

// Local Includes
#include "libmesh/elem_range.h"
#include "libmesh/threads_allocators.h"
#include "libmesh/parallel_object.h"

// C++ includes
#include <vector>
#include <unordered_set>

namespace libMesh
{

// Forward declarations
class DofMap;
class CouplingMatrix;

/**
 * This defines the sparsity pattern, or graph, of a sparse matrix.
 * The format is quite simple -- the global indices of the nonzero entries
 * in each row are packed into a vector.  The global indices (i,j) of the
 * nth nonzero entry of row i are given by j = sparsity_pattern[i][n];
 *
 * \author Roy Stogner
 * \date 2010
 * \brief Defines the sparsity pattern of a sparse matrix.
 */
namespace SparsityPattern // use a namespace so member classes can be forward-declared.
{
typedef std::vector<dof_id_type, Threads::scalable_allocator<dof_id_type>> Row;
class Graph : public std::vector<Row> {};

class NonlocalGraph : public std::map<dof_id_type, Row> {};

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
   * Abstract base class to be used to add user-defined implicit
   * degree of freedom couplings.
   */
  class AugmentSparsityPattern
  {
  public:
    virtual ~AugmentSparsityPattern () {}

    /**
     * User-defined function to augment the sparsity pattern.
     */
    virtual void augment_sparsity_pattern (SparsityPattern::Graph & sparsity,
                                           std::vector<dof_id_type> & n_nz,
                                           std::vector<dof_id_type> & n_oz) = 0;
  };



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
class Build : public ParallelObject
{
public:
  Build (const DofMap & dof_map_in,
         const CouplingMatrix * dof_coupling_in,
         const std::set<GhostingFunctor *> & coupling_functors_in,
         const bool implicit_neighbor_dofs_in,
         const bool need_full_sparsity_pattern_in);

  /**
   * Special functions.
   * - The Build object holds references to DofMap and coupling
   *   functors, therefore it can't be assigned.
   */
  Build (const Build &) = default;
  Build & operator= (const Build &) = delete;
  Build (Build &&) = default;
  Build & operator= (Build &&) = delete;
  ~Build () = default;

  /**
   * Splitting constructor, for use in multithreaded loops.
   */
  Build (Build & other, Threads::split);

  /**
   * Add entries from a range of elements to this object's sparsity
   * pattern.
   */
  void operator()(const ConstElemRange & range);

  /**
   * Combine the sparsity pattern in \p other with this object's
   * sparsity pattern.  Useful in multithreaded loops.
   */
  void join (const Build & other);

  /**
   * Send sparsity pattern data relevant to other processors to those
   * processors, and receive and incorporate data relevant to us.
   */
  void parallel_sync ();

  /**
   * Rows of sparse matrix indices, indexed by the offset from the
   * first DoF on this processor.
   */
  const SparsityPattern::Graph & get_sparsity_pattern() const
  { return sparsity_pattern; }

  /**
   * Rows of sparse matrix indices, mapped from global DoF number,
   * which belong on other processors.  Stored here only temporarily
   * until a parallel_sync() sends them where they belong.
   */
  const SparsityPattern::NonlocalGraph & get_nonlocal_pattern() const
  { return nonlocal_pattern; }

  /**
   * The number of on-processor nonzeros in my portion of the
   * global matrix.
   */
  const std::vector<dof_id_type> & get_n_nz() const
  { return n_nz; }

  /**
   * The number of off-processor nonzeros in my portion of the
   * global matrix.
   */
  const std::vector<dof_id_type> & get_n_oz() const
  { return n_oz; }

  /**
   * Let a user-provided AugmentSparsityPattern subclass modify our
   * sparsity structure.
   */
  void apply_extra_sparsity_object(SparsityPattern::AugmentSparsityPattern & asp);

  /**
   * Let a user-provided function modify our sparsity structure.
   */
  void apply_extra_sparsity_function(void (*func)(SparsityPattern::Graph & sparsity,
                                                   std::vector<dof_id_type> & n_nz,
                                                   std::vector<dof_id_type> & n_oz,
                                                   void * context),
                                      void * context)
  { func(sparsity_pattern, n_nz, n_oz, context); }

  /**
   * Clear the "full" details of our sparsity structure, leaving only
   * the counts of non-zero entries.
   */
  void clear_full_sparsity()
  {
    sparsity_pattern.clear();
    nonlocal_pattern.clear();
  }

private:
  const DofMap & dof_map;
  const CouplingMatrix * dof_coupling;
  const std::set<GhostingFunctor *> & coupling_functors;
  const bool implicit_neighbor_dofs;
  const bool need_full_sparsity_pattern;

  // If there are "spider" nodes in the mesh (i.e. a single node which
  // is connected to many 1D elements) and Constraints, we can end up
  // sorting the same set of DOFs multiple times in handle_vi_vj(),
  // only to find that it has no net effect on the final sparsity. In
  // such cases it is much faster to keep track of (element_dofs_i,
  // element_dofs_j) pairs which have already been handled and not
  // repeat the computation. We use this data structure to keep track
  // of hashes of sets of dofs we have already seen, thus avoiding
  // unnecessary caluclations.
  std::unordered_set<dof_id_type> hashed_dof_sets;

  void handle_vi_vj(const std::vector<dof_id_type> & element_dofs_i,
                    const std::vector<dof_id_type> & element_dofs_j);

  void sorted_connected_dofs(const Elem * elem,
                             std::vector<dof_id_type> & dofs_vi,
                             unsigned int vi);

#ifndef LIBMESH_ENABLE_DEPRECATED
private:
#endif

  SparsityPattern::Graph sparsity_pattern;

  SparsityPattern::NonlocalGraph nonlocal_pattern;

  std::vector<dof_id_type> n_nz;

  std::vector<dof_id_type> n_oz;
};

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

  libmesh_assert_greater (std::distance (begin,  middle), 0);
  libmesh_assert_greater (std::distance (middle, end), 0);
  libmesh_assert (std::unique (begin,  middle) == middle);
  libmesh_assert (std::unique (middle, end) == end);

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

#endif // LIBMESH_SPARSITY_PATTERN_H
