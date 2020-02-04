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



// Local includes
#include "libmesh/fe_xyz_impl.h"
#include "libmesh/fe_xyz_shape_0D_impl.h"
#include "libmesh/fe_xyz_shape_1D_impl.h"
#include "libmesh/fe_xyz_shape_2D_impl.h"
#include "libmesh/fe_xyz_shape_3D_impl.h"

namespace libMesh
{
// Explicit instantiations for Real
template struct FEShim<0,XYZ,Real>;
template struct FEShim<1,XYZ,Real>;
template struct FEShim<2,XYZ,Real>;
template struct FEShim<3,XYZ,Real>;

// Explicit instantiations for non-static FEXYZ member functions.
// These non-static member functions map more naturally to explicit
// instantiations than the functions above:
//
// 1.)  Since they are member functions, they rely on
// private/protected member data, and therefore do not work well
// with the "anonymous function call" model we've used above for
// the specializations.
//
// 2.) There is (IMHO) less chance of the linker calling the
// wrong version of one of these member functions, since there is
// only one FEXYZ.
template void  FEXYZ<0,Real>::init_shape_functions(const std::vector<Point> &, const Elem *);
template void  FEXYZ<1,Real>::init_shape_functions(const std::vector<Point> &, const Elem *);
template void  FEXYZ<2,Real>::init_shape_functions(const std::vector<Point> &, const Elem *);
template void  FEXYZ<3,Real>::init_shape_functions(const std::vector<Point> &, const Elem *);

template void  FEXYZ<0,Real>::compute_shape_functions(const Elem *,const std::vector<Point> &);
template void  FEXYZ<1,Real>::compute_shape_functions(const Elem *,const std::vector<Point> &);
template void  FEXYZ<2,Real>::compute_shape_functions(const Elem *,const std::vector<Point> &);
template void  FEXYZ<3,Real>::compute_shape_functions(const Elem *,const std::vector<Point> &);

} // namespace libMesh
