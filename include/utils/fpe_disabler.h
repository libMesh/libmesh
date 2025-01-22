// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FPE_DISABLER_H
#define LIBMESH_FPE_DISABLER_H


// Local includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <cfenv>

namespace libMesh
{

/**
 * The \p FPEDisabler class puts Floating-Point Exception (FPE)
 * trapping on hold during its lifetime, restoring the FE environment
 * upon destruction.  This allows code to ignore known/expected FPE
 * sources, e.g. to work around a bug in HDF5 1.14.3 or to test for
 * errors in a more sophisticated way and throw a C++ exception
 * instead.
 *
 * FE environment changes should not be made during the lifetime of
 * this object, as they will not be preserved when the destructor
 * restores the original environment.
 *
 * \author Roy Stogner
 * \date 2024
 * \brief Temporarily disables Floating-Point Exception (FPE) trapping
 */
struct FPEDisabler
{
  FPEDisabler()
  {
    std::feholdexcept(&old_env);
  }

  ~FPEDisabler()
  {
    std::fesetenv(&old_env);
  }

  std::fenv_t old_env;
};

} // namespace libMesh

#endif // LIBMESH_PERFLOG_H
