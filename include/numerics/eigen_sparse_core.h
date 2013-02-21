// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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




#ifndef LIBMESH_EIGEN_SPARSE_CORE_H
#define LIBMESH_EIGEN_SPARSE_CORE_H



#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_EIGEN

// Local includes
#include "libmesh/id_types.h"

// C++ includes

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/SparseCore>



namespace libMesh
{

  // must be a signed type!!
#if LIBMESH_DOF_ID_BYTES == 1
  typedef int8_t eigen_idx_type;
#elif LIBMESH_DOF_ID_BYTES == 2
  typedef int16_t eigen_idx_type;
#elif LIBMESH_DOF_ID_BYTES == 8
  typedef int64_t eigen_idx_type;
#else // LIBMESH_DOF_ID_BYTES = 4 (default)
  typedef int32_t eigen_idx_type;
#endif

  typedef Eigen::SparseMatrix<Number, Eigen::ColMajor, eigen_idx_type> EigenSM;
  //typedef Eigen::SparseVector<Number, Eigen::ColMajor, eigen_idx_type> EigenSV;
  typedef Eigen::Matrix<Number, Eigen::Dynamic, 1> EigenSV;
} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_EIGEN
#endif // LIBMESH_EIGEN_SPARSE_CORE_H
