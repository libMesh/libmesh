// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_STATISTICS_H
#define LIBMESH_STATISTICS_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/id_types.h"

// C++ includes
#include <vector>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>

namespace libMesh
{

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
 * \author John W. Peterson
 * \date 2002
 */
template <typename T>
class StatisticsVector : public std::vector<T>
{
public:

  /**
   * Call the std::vector constructor.
   */
  explicit
  StatisticsVector(dof_id_type i=0) : std::vector<T> (i) {}

  /**
   * Call the std::vector constructor, fill each entry with \p val
   */
  StatisticsVector(dof_id_type i, T val) : std::vector<T> (i,val) {}

  /**
   * Destructor.  Virtual so we can derive from the \p StatisticsVector
   */
  virtual ~StatisticsVector () {}


  /**
   * Returns the l2 norm of the data set.
   */
  virtual Real l2_norm() const;

  /**
   * Returns the minimum value in the data set.
   */
  virtual T minimum() const;

  /**
   * Returns the maximum value in the data set.
   */
  virtual T maximum() const;

  /**
   * Returns the mean value of the
   * data set using a recurrence relation.
   * Source: GNU Scientific Library
   */
  virtual Real mean() const;

  /**
   * Returns the median (e.g. the middle)
   * value of the data set.
   * This function modifies the
   * original data by sorting, so it
   * can't be called on const objects.
   * Source: GNU Scientific Library
   */
  virtual Real median();

  /**
   * A const version of the median funtion.
   * Requires twice the memory of original
   * data set but does not change the original.
   */
  virtual Real median() const;

  /**
   * Computes the variance of the data set.
   * Uses a recurrence relation to prevent
   * data overflow for large sums.
   * Note: The variance is equal to the
   * standard deviation squared.
   * Source: GNU Scientific Library
   */
  virtual Real variance() const
  { return this->variance(this->mean()); }

  /**
   * Computes the variance of the data set
   * where the \p mean is provided. This is useful
   * for efficiency when you have already calculated
   * the mean. Uses a recurrence relation to prevent
   * data overflow for large sums.
   * Note: The variance is equal to the
   * standard deviation squared.
   * Source: GNU Scientific Library
   */
  virtual Real variance(const Real known_mean) const;

  /**
   * Computes the standard deviation of
   * the data set, which is simply the square-root
   * of the variance.
   */
  virtual Real stddev() const
  { return std::sqrt(this->variance()); }

  /**
   * Computes the standard deviation of
   * the data set, which is simply the square-root
   * of the variance.  This method can be used for
   * efficiency when the \p mean has already been computed.
   */
  virtual Real stddev(const Real known_mean) const
  { return std::sqrt(this->variance(known_mean)); }

  /**
   * Divides all entries by the largest entry and
   * stores the result
   */
  void normalize();

  /**
   * Computes and returns a histogram with n_bins bins for the data
   * set.  For simplicity, the bins are assumed to be of uniform size.
   * Upon return, the bin_members vector will contain unsigned
   * integers which give the number of members in each bin.
   * WARNING: This non-const function sorts the vector, changing its
   * order.
   * Source: GNU Scientific Library
   */
  virtual void histogram (std::vector<dof_id_type> & bin_members,
                          unsigned int n_bins=10);

  /**
   * Generates a Matlab/Octave style file which can be used to
   * make a plot of the histogram having the desired number of bins.
   * Uses the histogram(...) function in this class
   * WARNING: The histogram(...) function is non-const, and changes
   * the order of the vector.
   */
  void plot_histogram(const processor_id_type my_procid,
                      const std::string & filename,
                      unsigned int n_bins);

  /**
   * A const version of the histogram function.
   */
  virtual void histogram (std::vector<dof_id_type> & bin_members,
                          unsigned int n_bins=10) const;

  /**
   * Returns a vector of dof_id_types which correspond
   * to the indices of every member of the data set
   * below the cutoff value "cut".
   */
  virtual std::vector<dof_id_type> cut_below(Real cut) const;

  /**
   * Returns a vector of dof_id_types which correspond
   * to the indices of every member of the data set
   * above the cutoff value cut.  I chose not to combine
   * these two functions since the interface is cleaner
   * with one passed parameter instead of two.
   */
  virtual std::vector<dof_id_type> cut_above(Real cut) const;
};

} // namespace libMesh

#endif // LIBMESH_STATISTICS_H
