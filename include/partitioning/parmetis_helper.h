// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PARMETIS_HELPER_H
#define LIBMESH_PARMETIS_HELPER_H

// Local Includes
#include "libmesh/libmesh_config.h"

// C++ Includes
#include <vector>

// Include the ParMETIS header files.  We need this so we can use
// ParMetis' idx_t and real_t types directly.
#ifdef LIBMESH_HAVE_PARMETIS
namespace Parmetis {
extern "C" {
#     include "libmesh/ignore_warnings.h"
#     include "parmetis.h"
#     include "libmesh/restore_warnings.h"
}
}
#endif // LIBMESH_HAVE_PARMETIS


namespace libMesh
{

/**
 * The \p ParmetisHelper class allows us to use a 'pimpl' strategy in
 * the ParmetisPartitioner class.  Since we don't include the
 * parmetis.h header file here, we don't have to install it, either.
 * This class is empty when Parmetis is not available, otherwise it is
 * simply a data container.
 *
 * \author John W. Peterson
 * \date 2015
 * \brief Pointer-to-implementation class used by ParmetisPartitioner.
 */
class ParmetisHelper
{
public:
  /**
   * Defaulted constructors, assignment operators, and destructor.
   */
  ParmetisHelper () = default;
  ParmetisHelper (const ParmetisHelper &) = default;
  ParmetisHelper (ParmetisHelper &&) = default;
  ParmetisHelper & operator= (const ParmetisHelper &) = default;
  ParmetisHelper & operator= (ParmetisHelper &&) = default;
  ~ParmetisHelper () = default;

#ifdef LIBMESH_HAVE_PARMETIS

  /**
   * Data structures used by ParMETIS to describe the connectivity graph
   * of the mesh. Consult the ParMETIS documentation.
   */
  std::vector<Parmetis::idx_t>  vtxdist;
  std::vector<Parmetis::idx_t>  xadj;
  std::vector<Parmetis::idx_t>  adjncy;

  // We use dof_id_type for part so we can pass it directly to
  // Partitioner:: methods expecting that type.
  std::vector<dof_id_type>  part;

  // But we plan to pass a pointer to part as a buffer to ParMETIS, so
  // it had better be using a simply reinterpretable type!
  static_assert(sizeof(Parmetis::idx_t) == sizeof(dof_id_type),
                "ParMETIS and libMesh ID sizes must match!");

  std::vector<Parmetis::real_t> tpwgts;
  std::vector<Parmetis::real_t> ubvec;
  std::vector<Parmetis::idx_t>  options;
  std::vector<Parmetis::idx_t>  vwgt;

  Parmetis::idx_t wgtflag;
  Parmetis::idx_t ncon;
  Parmetis::idx_t numflag;
  Parmetis::idx_t nparts;
  Parmetis::idx_t edgecut;

#endif // LIBMESH_HAVE_PARMETIS
};

} // namespace libMesh

#endif // LIBMESH_PARMETIS_HELPER_H
