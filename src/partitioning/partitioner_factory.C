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



// Local Includes
#include "libmesh/libmesh_config.h"
#include "libmesh/centroid_partitioner.h"
#include "libmesh/metis_partitioner.h"
#include "libmesh/parmetis_partitioner.h"
#include "libmesh/linear_partitioner.h"
#include "libmesh/hilbert_sfc_partitioner.h"
#include "libmesh/morton_sfc_partitioner.h"
#include "libmesh/factory.h"

namespace libMesh
{


//-------------------------------------------------
// Full specialization for the Factory<Partitioner>
template<>
std::map<std::string, Factory<Partitioner> *> &
Factory<Partitioner>::factory_map()
{
  static std::map<std::string, Factory<Partitioner> *> _map;
  return _map;
}



// ------------------------------------------------------------
// Register Partitioning classes with the factory.  These will never
// be called from user code, they just need to get instantiated.  Hide
// them in an anonymous namespace to prevent name clashes
namespace {

#ifdef LIBMESH_HAVE_METIS
FactoryImp<MetisPartitioner,      Partitioner> metis    ("Metis");
#endif

#ifdef LIBMESH_HAVE_PARMETIS
FactoryImp<ParmetisPartitioner,   Partitioner> parmetis ("Parmetis");
#endif

#ifdef LIBMESH_HAVE_SFCURVES
FactoryImp<SFCPartitioner,        Partitioner> sfc      ("SFCurves");
FactoryImp<HilbertSFCPartitioner, Partitioner> hilbert  ("Hilbert");
FactoryImp<MortonSFCPartitioner,  Partitioner> morton   ("Morton");
#endif

FactoryImp<LinearPartitioner,     Partitioner> linear   ("Linear");
FactoryImp<CentroidPartitioner,   Partitioner> centroid ("Centroid");

}

} // namespace libMesh
