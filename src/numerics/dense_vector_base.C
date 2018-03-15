// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ includes
#include <iostream>
#include <iomanip> // for std::setw


// Local Includes
#include "libmesh/dense_vector_base.h"

namespace libMesh
{

template<typename T>
void DenseVectorBase<T>::print_scientific (std::ostream & os, unsigned precision) const
{
  // save the initial format flags
  std::ios_base::fmtflags os_flags = os.flags();

  // Print the vector entries.
  for (unsigned int i=0; i<this->size(); i++)
    os << std::setw(10)
       << std::scientific
       << std::setprecision(precision)
       << this->el(i)
       << std::endl;

  // reset the original format flags
  os.flags(os_flags);
}



template<typename T>
void DenseVectorBase<T>::print (std::ostream & os) const
{
  for (unsigned int i=0; i<this->size(); i++)
    os << std::setw(8)
       << this->el(i)
       << std::endl;
}


//--------------------------------------------------------------
// Explicit instantiations
template class DenseVectorBase<Real>;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template class DenseVectorBase<Complex>;
#endif

} // namespace libMesh
