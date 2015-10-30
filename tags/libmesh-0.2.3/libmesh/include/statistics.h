// $Id: statistics.h,v 1.6 2003-02-03 03:51:49 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __statistics_h__
#define __statistics_h__

// C++ includes
#include <vector>


// Local includes
#include "mesh_common.h"

/**
 * The StatisticsVector class is derived from
 * the std::vector<> and therefore has all of
 * its useful features.  It was designed 
 * to not have any internal state, i.e. no
 * public or private data members.  Also,
 * it was only designed for classes and
 * types for which the operators +,*,/ have
 * meaining, specifically floats, doubles,
 * ints, etc.  The main reason for this design
 * decision was to allow a std::vector<> to
 * be successfully cast to a StatisticsVector,
 * thereby enabling its additional functionality.
 * We do not anticipate any problems with 
 * deriving from an stl container which lacks
 * a virtual destructor in this case.
 *
 * Where manipulation of the data set was necessary
 * (for example sorting) two versions of member
 * functions have been implemented.  The non-const
 * versions perform sorting directly in the
 * data set, invalidating pointers and changing
 * the entries.  const versions of the same
 * functions are generally available, and will be
 * automatically invoked on const StatisticsVector objects.
 * A draw-back to the const versions is that they
 * simply make a copy of the original object and
 * therefore double the original memory requirement
 * for the data set.
 *
 * Most of the actual code was copied or adapted from
 * the GNU Scientific Library (GSL). More precisely,
 * the recursion relations for computing the mean were
 * implemented in order to avoid possible problems with
 * buffer overruns.
 *
 * @author John W. Peterson, 2002
 */

// ------------------------------------------------------------
// StatisticsVector class definition
template <typename T>
class StatisticsVector : public std::vector<T>
{
 public:

  /**
   * Call the std::vector constructor.
   */
  StatisticsVector(unsigned int i=0) { this->resize(i); };
  
  /**
   * Returns the mean value of the
   * data set using a recurrence relation.
   * Source: GNU Scientific Library
   */
  Real mean() const;

  /**
   * Returns the median (e.g. the middle)
   * value of the data set.
   * This function modifies the
   * original data by sorting, so it
   * can't be called on const objects.
   * Source: GNU Scientific Library
   */
  Real median();

  /**
   * A const version of the median funtion.
   * Requires twice the memory of original
   * data set but does not change the original.
   */
  Real median() const;

  /**
   * Computes the variance of the data set.
   * Uses a recurrence relation to prevent
   * data overflow for large sums.
   * Note: The variance is equal to the
   * standard deviation squared.  The variance
   * is normalized by N in this case.
   * Source: GNU Scientific Library
   */
  Real variance() const;

  /**
   * Computes and returns a histogram with
   * n_bins bins for the data set.
   * For simplicity, the bins are assumed to
   * be of uniform size.  The return values
   * are unsigned integers which give the number
   * of members in each bin.
   * Source: GNU Scientific Library
   */
  std::vector<unsigned int> histogram(unsigned int n_bins=10);

  /**
   * A const version of the median function.
   * Requires twice the memory of the data set
   * but does not change the original.
   */
  std::vector<unsigned int> histogram(unsigned int n_bins=10) const;

  /**
   * Returns a vector of unsigned ints which correspond
   * to the indices of every member of the data set
   * below the cutoff value cut.
   */
  std::vector<unsigned int> cut_below(Real cut) const;

  /**
   * Returns a vector of unsigned ints which correspond
   * to the indices of every member of the data set
   * above the cutoff value cut.  I chose not to combine
   * these two functions since the interface is cleaner
   * with one passed parameter instead of two.
   */
  std::vector<unsigned int> cut_above(Real cut) const;
  
  
 private:
  
};

#endif
