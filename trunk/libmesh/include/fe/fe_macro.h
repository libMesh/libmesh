// $Id: fe_macro.h,v 1.9 2005-05-07 00:22:29 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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


#ifndef __fe_macro_h__
#define __fe_macro_h__


// Local includes


/**
 * This macro helps in instantiating specific versions
 * of the \p FE class.  Simply include this file, and
 * instantiate at the end for the desired dimension.
 */

#ifndef ENABLE_HIGHER_ORDER_SHAPES

#define INSTANTIATE_FE(_dim)   template class FE< (_dim), CLOUGH>; \
                               template class FE< (_dim), HIERARCHIC>;   \
                               template class FE< (_dim), LAGRANGE>;   \
                               template class FE< (_dim), MONOMIAL>;   \
                               template class FE< (_dim), XYZ>

#define INSTANTIATE_IMAP(_dim) \
  template void FE<_dim,CLOUGH>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&); \
  template void FE<_dim,HIERARCHIC>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&);\
  template void FE<_dim,LAGRANGE>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&);\
  template void FE<_dim,MONOMIAL>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&);\
  template void FE<_dim,XYZ>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&)
 
#else //ENABLE_HIGHER_ORDER_SHAPES

#define INSTANTIATE_FE(_dim)   template class FE< (_dim), CLOUGH>; \
                               template class FE< (_dim), HIERARCHIC>;   \
                               template class FE< (_dim), LAGRANGE>;   \
                               template class FE< (_dim), MONOMIAL>;   \
                               template class FE< (_dim), SZABAB>;     \
                               template class FE< (_dim), XYZ>

#define INSTANTIATE_IMAP(_dim) \
  template void FE<_dim,CLOUGH>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&); \
  template void FE<_dim,HIERARCHIC>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&);\
  template void FE<_dim,LAGRANGE>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&);\
  template void FE<_dim,MONOMIAL>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&);\
  template void FE<_dim,SZABAB>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&);\
  template void FE<_dim,XYZ>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&)

#endif //ENABLE_HIGHER_ORDER_SHAPES

#define INSTANTIATE_MBRF(_dim,_t) \
  template unsigned int FE<_dim,_t>::n_dofs_at_node(ElemType,Order,unsigned int); \
  template unsigned int FE<_dim,_t>::n_dofs(ElemType,Order);\
  template bool         FE<_dim,_t>::shapes_need_reinit() const;\
  template FEContinuity FE<_dim,_t>::get_continuity() const;\
  template unsigned int FE<_dim,_t>::n_dofs_per_elem(ElemType,Order);\
  template void         FE<_dim,_t>::nodal_soln(const Elem*,const Order,const std::vector<Number>&,std::vector<Number>&)


#endif
