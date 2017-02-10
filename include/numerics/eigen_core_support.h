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




#ifndef LIBMESH_EIGEN_CORE_SUPPORT_H
#define LIBMESH_EIGEN_CORE_SUPPORT_H



#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_EIGEN

// Local includes
#include "libmesh/id_types.h"

// C++ includes

// // hack to avoid MatType collision...
// #undef libMeshSaveMatType
// #ifdef MatType
// #  define MatType libMeshSaveMatType
// #  undef  MatType
// #endif

// Eigen uses deprecated std::binder1st/2nd classes, which GCC warns about.
// GCC 6.1 also warns about misleading indentation from if blocks in
// Eigen.
#include "libmesh/ignore_warnings.h"

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "libmesh/restore_warnings.h"

// #ifdef libMeshSaveMatType
// #  define libMeshSaveMatType MatType
// #  undef  libMeshSaveMatType
// #endif


namespace libMesh
{

// must be a signed type!!
#if LIBMESH_DOF_ID_BYTES == 1
// Workaround for Eigen bug
// typedef int8_t eigen_idx_type;
typedef int32_t eigen_idx_type;
#elif LIBMESH_DOF_ID_BYTES == 2
// Workaround for Eigen bug
// typedef int16_t eigen_idx_type;
typedef int32_t eigen_idx_type;
#elif LIBMESH_DOF_ID_BYTES == 8
typedef int64_t eigen_idx_type;
#else // LIBMESH_DOF_ID_BYTES = 4 (default)
typedef int32_t eigen_idx_type;
#endif

// We have to use RowMajor SparseMatrix storage for our preallocation to work
typedef Eigen::SparseMatrix<Number, Eigen::RowMajor, eigen_idx_type> EigenSM;
typedef Eigen::Matrix<Number, Eigen::Dynamic, 1> EigenSV;
} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_EIGEN
#endif // LIBMESH_EIGEN_CORE_SUPPORT_H
