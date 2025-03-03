// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public  License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local Includes
#include "libmesh/overlap_coupling.h"

#include "libmesh/elem.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/remote_elem.h"
#include "libmesh/libmesh_logging.h"

// C++ Includes
#include <unordered_set>

namespace libMesh
{

OverlapCoupling::OverlapCoupling() :
  _dof_coupling(nullptr)
{
  _q_rules[0] = QBase::build(QSIMPSON, 1, SECOND);
  _q_rules[1] = QBase::build(QSIMPSON, 2, SECOND);
  _q_rules[2] = QBase::build(QSIMPSON, 3, SECOND);
}



OverlapCoupling::OverlapCoupling(const OverlapCoupling & other) :
  GhostingFunctor(other),
  _dof_coupling(other._dof_coupling)
{
  for (auto i : make_range(2))
    if (other._q_rules[i].get())
      _q_rules[i] = other._q_rules[i]->clone();
}



void OverlapCoupling::set_quadrature_rule
  (std::unique_ptr<QBase> new_q_rule)
{
  libmesh_assert(new_q_rule.get());
  const unsigned int dim = new_q_rule->get_dim();
  _q_rules[dim-1] = std::move(new_q_rule);
}



void OverlapCoupling::mesh_reinit()
{
  // We'll need a master point locator, so we'd better have a mesh
  // to build it on.
  libmesh_assert(_mesh);

  // Make sure an up-to-date master point locator has been
  // constructed; we'll need to grab sub-locators soon.
  _mesh->sub_point_locator();
}



void OverlapCoupling::operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   map_type & coupled_elements)
{
  LOG_SCOPE("operator()", "OverlapCoupling");

  std::unique_ptr<PointLocatorBase> point_locator
    = _mesh->sub_point_locator();

  const std::vector<Point> & xyz = this->_fe_map.get_xyz();

  for (const auto & elem : as_range(range_begin, range_end))
  {
    const unsigned int dim = elem->dim();

    QBase & qrule = *_q_rules[dim-1];
    qrule.init(*elem);

    _fe_map.init_reference_to_physical_map(dim, qrule.get_points(), elem);
    _fe_map.compute_map(dim, qrule.get_weights(), elem, /*d2phi=*/ false);

    std::set<const Elem *> overlapping_elements;
    for (const Point & qp : xyz)
      (*point_locator)(qp, overlapping_elements);

    // We should at least find ourselves
    libmesh_assert(overlapping_elements.count(elem));

    for (const Elem * e : overlapping_elements)
      if (e->processor_id() != p)
        coupled_elements.emplace(e, _dof_coupling);
  }
}


} // namespace libMesh
