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



// Local Includes
#include "libmesh/qoi_set.h"
#include "libmesh/system.h"

// C++ Includes
#include <vector>


namespace libMesh
{

QoISet::QoISet(const System & sys) : _indices(sys.qoi.size(), true) {}



unsigned int QoISet::size (const System & sys) const
{
  unsigned int qoi_count = 0;
  for (std::size_t i=0; i != sys.qoi.size(); ++i)
    if (this->has_index(i))
      qoi_count++;
  return qoi_count;
}



void QoISet::add_indices(const std::vector<unsigned int> & indices)
{
  unsigned int max_size = 0;
  for (std::vector<unsigned int>::const_iterator i = indices.begin();
       i != indices.end(); ++i)
    max_size = std::max(max_size, *i + 1);

  _indices.resize(max_size);

  for (std::vector<unsigned int>::const_iterator i = indices.begin();
       i != indices.end(); ++i)
    _indices[*i] = true;
}



inline
void QoISet::remove_indices(const std::vector<unsigned int> & indices)
{
  for (std::vector<unsigned int>::const_iterator i = indices.begin();
       i != indices.end(); ++i)
    _indices[*i] = false;
}

} // namespace libMesh
