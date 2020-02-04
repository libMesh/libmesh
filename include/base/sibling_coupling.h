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



#ifndef LIBMESH_SIBLING_COUPLING_H
#define LIBMESH_SIBLING_COUPLING_H

// Local Includes
#include "libmesh/ghosting_functor.h"

namespace libMesh
{

/**
 * This class adds coupling (for use in send_list construction)
 * between active elements and all descendants of their parent.
 *
 * This may be useful when constructing projections onto the parent
 * element for error indicator evaluation.
 *
 * \author Roy H. Stogner
 * \date 2016
 */
template <typename RealType = Real>
class SiblingCouplingTempl : public GhostingFunctorTempl<RealType>
{
public:
  typedef SiblingCouplingTempl<RealType> SiblingCoupling;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef ElemTempl<RealType> Elem;
  typedef NodeTempl<RealType> Node;
  using map_type = std::unordered_map<const ElemTempl<RealType> *, const CouplingMatrix*>;

  /**
   * Constructor.
   */
  SiblingCouplingTempl() :
    _dof_coupling(nullptr) {}

  // Change coupling matrix after construction
  void set_dof_coupling(const CouplingMatrix * dof_coupling)
  { _dof_coupling = dof_coupling; }

  /**
   * For the specified range of active elements, find any sibling
   * elements which should be evaluable too.
   */
  virtual void operator() (const typename MeshBaseTempl<RealType>::const_element_iterator & range_begin,
                           const typename MeshBaseTempl<RealType>::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements);

private:

  const CouplingMatrix * _dof_coupling;
};

typedef SiblingCouplingTempl<Real> SiblingCoupling;

} // namespace libMesh

#endif // LIBMESH_SIBLING_COUPLING_H
