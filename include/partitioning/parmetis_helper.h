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

#ifndef LIBMESH_PARMETIS_HELPER_H
#define LIBMESH_PARMETIS_HELPER_H

// Local Includes -----------------------------------
#include "libmesh/libmesh_config.h"

// C++ Includes   -----------------------------------
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
 * the ParmetisPartitioner class.  Since we don't install this header
 * file, apps do not include parmetis.h, and consequently we don't
 * have to install it, either.  This class is empty when Parmetis
 * is not available, when it is it is simply a data container.
 */
class ParmetisHelper
{
public:
  /**
   * Constructor.
   */
  ParmetisHelper () {}

#ifdef LIBMESH_HAVE_PARMETIS

  /**
   * Data structures used by ParMETIS to describe the connectivity graph
   * of the mesh.  Consult the ParMETIS documentation.
   */
  std::vector<Parmetis::idx_t>  vtxdist;
  std::vector<Parmetis::idx_t>  xadj;
  std::vector<Parmetis::idx_t>  adjncy;
  std::vector<Parmetis::idx_t>  part;
  std::vector<Parmetis::real_t> tpwgts;
  std::vector<Parmetis::real_t> ubvec;
  std::vector<Parmetis::idx_t>  options;
  std::vector<Parmetis::idx_t>  vwgt;

  Parmetis::idx_t wgtflag;
  Parmetis::idx_t ncon;
  Parmetis::idx_t numflag;
  Parmetis::idx_t nparts;
  Parmetis::idx_t edgecut;

#endif
};

}

#endif
