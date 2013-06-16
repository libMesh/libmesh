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



#ifndef LIBMESH_DG_FEM_CONTEXT_H
#define LIBMESH_DG_FEM_CONTEXT_H

// Local Includes
#include "libmesh/fem_context.h"


namespace libMesh
{

/**
 * This class extends FEMContext in order to provide extra data
 * required to perform local element residual and Jacobian assembly
 * in the case of a discontinuous Galerkin (DG) discretization.
 */

// ------------------------------------------------------------
// DGFEMContext class definition

class DGFEMContext : public FEMContext
{
public:

  /**
   * Constructor.  Allocates some but fills no data structures.
   */
  explicit
  DGFEMContext (const System &sys);

  /**
   * Destructor.
   */
  virtual ~DGFEMContext ();

  /**
   * Override pre_fe_reinit to set a boolean flag so that by
   * default DG terms are assumed to be inactive.
   */
  virtual void pre_fe_reinit(const System& sys, const Elem *e);

  /**
   * Const accessor for neighbor residual.
   */
  const DenseVector<Number>& get_neighbor_residual() const
  { return _neighbor_residual; }

  /**
   * Non-const accessor for neighbor residual.
   */
  DenseVector<Number>& get_neighbor_residual()
  { return _neighbor_residual; }

  /**
   * Const accessor for neighbor residual of a particular variable corresponding
   * to the variable index argument.
   */
  const DenseSubVector<Number>& get_neighbor_residual( unsigned int var ) const
  { return *(_neighbor_residuals[var]); }

  /**
   * Non-const accessor for neighbor residual of a particular variable corresponding
   * to the variable index argument.
   */
  DenseSubVector<Number>& get_neighborresidual( unsigned int var )
  { return *(_neighbor_residuals[var]); }

  /**
   * Const accessor for element-neighbor Jacobian.
   */
  const DenseVector<Number>& get_elem_neighbor_jacobian() const
  { return _elem_neighbor_jacobian; }

  /**
   * Non-const accessor for element Jacobian.
   */
  DenseVector<Number>& get_elem_neighbor_jacobian()
  { return _elem_neighbor_jacobian; }

  /**
   * Const accessor for element-neighbor Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  const DenseSubMatrix<Number>& get_elem_neighbor_jacobian( unsigned int var1, unsigned int var2 ) const
  { return *(_elem_neighbor_subjacobians[var1][var2]); }

  /**
   * Non-const accessor for element-neighbor Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  DenseSubMatrix<Number>& get_elem_neighbor_jacobian( unsigned int var1, unsigned int var2 )
  { return *(_elem_neighbor_subjacobians[var1][var2]); }

  /**
   * Const accessor for element-neighbor Jacobian.
   */
  const DenseVector<Number>& get_neighbor_elem_jacobian() const
  { return _neighbor_elem_jacobian; }

  /**
   * Non-const accessor for element Jacobian.
   */
  DenseVector<Number>& get_neighbor_elem_jacobian()
  { return _neighbor_elem_jacobian; }

  /**
   * Const accessor for neighbor-element Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  const DenseSubMatrix<Number>& get_neighbor_elem_jacobian( unsigned int var1, unsigned int var2 ) const
  { return *(neighbor_elem_subjacobians[var1][var2]); }

  /**
   * Non-const accessor for neighbor-element Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  DenseSubMatrix<Number>& get_neighbor_elem_jacobian( unsigned int var1, unsigned int var2 )
  { return *(_neighbor_elem_subjacobians[var1][var2]); }

  /**
   * Const accessor for element-neighbor Jacobian.
   */
  const DenseVector<Number>& get_neighbor_neighbor_jacobian() const
  { return _neighbor_neighbor_jacobian; }

  /**
   * Non-const accessor for element Jacobian.
   */
  DenseVector<Number>& get_neighbor_neighbor_jacobian()
  { return _neighbor_neighbor_jacobian; }

  /**
   * Const accessor for neighbor-neighbor Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  const DenseSubMatrix<Number>& get_neighbor_neighbor_jacobian( unsigned int var1, unsigned int var2 ) const
  { return *(_neighbor_neighbor_subjacobians[var1][var2]); }

  /**
   * Non-const accessor for neighbor-neighbor Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  DenseSubMatrix<Number>& get_neighbor_neighbor_jacobian( unsigned int var1, unsigned int var2 )
  { return *(_neighbor_neighbor_subjacobians[var1][var2]); }

  /**
   * Set the neighbor element which we will use to assemble DG terms.
   * Note that we do not assume that this element is get_elem().neighbor(side)
   * because we also need to be able to handle the special case of DG terms on
   * "cracks" in a mesh to model certain types of interface conditions. In this
   * case, we need to be able to specify the neighbor element manually.
   */
  void set_neighbor(Elem& neighbor)
  { _neighbor = &neighbor; }

  /**
   * Accessor for current neighbor Elem object for assembling DG terms.
   */
  const Elem& get_neighbor() const
  { return *_neighbor; }

  /**
   * Are the DG terms active, i.e. have they been assembled?
   */
  bool are_dg_terms_active() const
  { return _dg_terms_active; }

private:

  /**
   * Current neighbor element for assembling DG terms.
   */
  const Elem *_neighbor;

  /**
   * Residual vector of the neighbor component.
   */
  DenseVector<Number> _neighbor_residual;

  /**
   * The element-neighbor Jacobian terms.
   * Test functions are from "element", trial functions are from "neighbor".
   */
  DenseMatrix<Number> _elem_neighbor_jacobian;

  /**
   * The neighbor-element Jacobian terms.
   * Test functions are from "neighbor", trial functions are from "element".
   */
  DenseMatrix<Number> _neighbor_elem_jacobian;
  
  /**
   * The neighbor-neighbor Jacobian terms.
   * Both trial and test functions are from "neighbor".
   */
  DenseMatrix<Number> _neighbor_neighbor_jacobian;

  /**
   * Element residual subvectors and Jacobian submatrices
   */
  std::vector<DenseSubVector<Number> *> _neighbor_subresiduals;
  std::vector<std::vector<DenseSubMatrix<Number> *> > _elem_neighbor_subjacobians;
  std::vector<std::vector<DenseSubMatrix<Number> *> > _neighbor_elem_subjacobians;
  std::vector<std::vector<DenseSubMatrix<Number> *> > _neighbor_neighbor_subjacobians;

  /**
   * Global Degree of freedom index lists for the neighbor element
   */
  std::vector<dof_id_type> _neighbor_dof_indices;
  std::vector<std::vector<dof_id_type> > _neighbor_dof_indices_var;

  /**
   * Finite element objects for each variable's
   * sides on the neighbor element.
   * We do not need FE objects for neighbor element
   * interior since we just need to handle DG interface
   * terms here.
   */
  std::map<FEType, FEBase *> _neighbor_side_fe;

  /**
   * Pointers to the same finite element objects on the neighbor element,
   * but indexed by variable number
   */
  std::vector<FEBase *> _neighbor_side_fe_var;

  /**
   * Quadrature rules for neighbor sides
   * Analogous to side_qrule in FEMContext.
   */
  QBase *_neighbor_side_qrule;

  /**
   * Boolean flag to indicate whether or not the DG terms have been
   * assembled and should be used in the global matrix assembly.
   */
  bool _dg_terms_active;
};


} // namespace libMesh

#endif // LIBMESH_FEM_CONTEXT_H
