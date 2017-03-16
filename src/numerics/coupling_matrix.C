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


#include "libmesh/coupling_matrix.h"


namespace libMesh {

CouplingMatrix & CouplingMatrix::operator&= (const CouplingMatrix & other)
{
  const std::size_t max_size = std::numeric_limits<std::size_t>::max();

  rc_type::iterator start_range = this->_ranges.begin();

  rc_type::const_iterator     other_range = other._ranges.begin();
  const rc_type::const_iterator other_end = other._ranges.end();

  for (; other_range != other_end; ++other_range)
    {
      std::size_t other_range_start = other_range->first;
      std::size_t other_range_end = other_range->second;

      // Find our range that might contain the start of the other
      // range.
      // lower_bound isn't *quite* what we want.
      // Because the other._ranges is sorted, we can contract this
      // search as we proceed, beginning with lb rather than at
      // begin() every time.
      rc_type::iterator lb =
        std::upper_bound (start_range, this->_ranges.end(),
                          std::make_pair(other_range_start, max_size));
      if (lb!=start_range)
        --lb;
      else
        lb=this->_ranges.end();

      start_range = lb;

      // If no range might contain the start of the new range then
      // we can just break out of here and start appending any
      // remaining ranges.
      if (lb == this->_ranges.end())
        break;

      // We did find a range which might contain the start of the new
      // range.
      const std::size_t lastloc  = lb->second;
      libmesh_assert_less_equal(lb->first, lastloc);
      libmesh_assert_less_equal(lb->first, other_range_start);

#ifdef DEBUG
      {
        CouplingMatrix::rc_type::const_iterator next = lb;
        next++;
        if (next != this->_ranges.end())
          {
            // Ranges should be sorted and should not touch
            libmesh_assert_greater(next->first, lastloc+1);
          }
      }
#endif

      CouplingMatrix::rc_type::iterator next = lb;
      next++;

      // We might need to extend this range or add a new range.
      // Merge contiguous ranges
      if (other_range_start <= lastloc)
        lb->second = other_range_end;
      // Or insert a new range.  This invalidates existing iterators,
      // but that's okay; the only iterator we need is for the newly
      // inserted range.
      else
        start_range = lb = this->_ranges.insert
          (next, std::make_pair(other_range_start, other_range_end));

      // At this point we have a range lb that may potentially overlap
      // subsequent existing ranges, in which case we need to merge
      // some.

      // First expand our range as necessary while finding ranges
      // which will be redundant later
      for (const std::size_t nextloc =
             (next == this->_ranges.end()) ?
             std::numeric_limits<std::size_t>::max() : next->first;
           nextloc <= lb->second; ++next)
        {
          // Ranges should be sorted and should not have been touching
          // initially
          libmesh_assert_greater(nextloc, lastloc+1);

          lb->second = std::max(lb->second, next->second);
        }

      CouplingMatrix::rc_type::iterator oldnext = lb;
      oldnext++;

      // Finally remove the redundant ranges
      this->_ranges.erase(oldnext, next);
    }

  // If we broke out early but we still have more other_ranges then
  // we can safely just append them to our ranges.
  for (; other_range != other_end; ++other_range)
    this->_ranges.push_back(*other_range);

  // Behave like a standard modification operator
  return *this;
}

}
