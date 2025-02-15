// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_OVERLAP_COUPLING_H
#define LIBMESH_OVERLAP_COUPLING_H

// Local Includes
#include "libmesh/ghosting_functor.h"

#include "libmesh/fe_map.h"
#include "libmesh/quadrature.h"

// C++ Includes
#include <memory>

namespace libMesh
{

// Forward declarations
class PeriodicBoundaries;


/**
 * This class implements ghosting of elements that overlap or touch at
 * at least one sampled point, even if no topological connection
 * between the elements exists.  This may be useful for immersed
 * methods, or for integration along both sides of mesh "slits" with
 * no neighbor information.
 *
 * \author Roy H. Stogner
 * \date 2025
 */
class OverlapCoupling : public GhostingFunctor
{
public:

  /**
   * Constructor.  Defaults to using a Simpson's Rule quadrature to
   * choose sampling points.  Users may wish to attach custom
   * quadrature rules matching those used for their system
   * integration.
   */
  OverlapCoupling();

  /**
   * Copy constructor.
   */
  OverlapCoupling(const OverlapCoupling & other) :
    GhostingFunctor(other),
    _dof_coupling(other._dof_coupling)
  {
    for (auto i : make_range(2))
      if (other._q_rules[i].get())
        _q_rules[i] = other._q_rules[i]->clone();
  }

  /**
   * A clone() is needed because GhostingFunctor can not be shared between
   * different meshes. The operations in GhostingFunctor are mesh dependent.
   */
  virtual std::unique_ptr<GhostingFunctor> clone () const override
  { return std::make_unique<OverlapCoupling>(*this); }

  // Change coupling matrix after construction
  void set_dof_coupling(const CouplingMatrix * dof_coupling)
  { _dof_coupling = dof_coupling; }

  // Change quadrature rule after construction.  This function takes
  // ownership of the rule object and uses it for elements of the
  // appropriate dimension.
  void set_quadrature_rule(std::unique_ptr<QBase> new_q_rule)
  {
    libmesh_assert(new_q_rule.get());
    const unsigned int dim = new_q_rule->get_dim();
    _q_rules[dim-1] = std::move(new_q_rule);
  }

  /**
   * We need an updated point locator to see what other elements might
   * share each quadrature point.
   */
  virtual void mesh_reinit () override;

  virtual void redistribute () override
  { this->mesh_reinit(); }

  virtual void delete_remote_elements() override
  { this->mesh_reinit(); }

  /**
   * For the specified range of active elements, find the elements
   * which will be coupled to them in the sparsity pattern.
   *
   * This will include whatever the point locator finds at each
   * quadrature point in the range.
   */
  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements) override;

private:

  const CouplingMatrix * _dof_coupling;

  // Quadrature rules for different element dimensions
  std::array<std::unique_ptr<QBase>, 3> _q_rules;

  // FE Map from quadrature space to physical space
  FEMap _fe_map;
};

} // namespace libMesh

#endif // LIBMESH_OVERLAP_COUPLING_H
