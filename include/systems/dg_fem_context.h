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
 *
 * \author David Knezevic
 * \date 2015
 * \brief Extends FEMContext to work for DG problems.
 */
class DGFEMContext : public FEMContext
{
public:

  /**
   * Constructor.  Allocates some but fills no data structures.
   */
  explicit
  DGFEMContext (const System & sys);

  /**
   * Destructor.
   */
  virtual ~DGFEMContext ();

  /**
   * Override side_fe_reinit to set a boolean flag so that by
   * default DG terms are assumed to be inactive. DG terms are
   * only active if neighbor_side_fe_reinit is called.
   */
  virtual void side_fe_reinit () override;

  /**
   * Initialize neighbor side data needed to assemble DG terms.
   * The neighbor element is determined by the current value of
   * get_neighbor().
   */
  void neighbor_side_fe_reinit ();

  /**
   * Accessor for neighbor dof indices
   */
  const std::vector<dof_id_type> & get_neighbor_dof_indices() const
  { return _neighbor_dof_indices; }

  /**
   * Accessor for element dof indices of a particular variable corresponding
   * to the index argument.
   */
  const std::vector<dof_id_type> & get_neighbor_dof_indices( unsigned int var ) const
  { return _neighbor_dof_indices_var[var]; }

  /**
   * Const accessor for neighbor residual.
   */
  const DenseVector<Number> & get_neighbor_residual() const
  { return _neighbor_residual; }

  /**
   * Non-const accessor for neighbor residual.
   */
  DenseVector<Number> & get_neighbor_residual()
  { return _neighbor_residual; }

  /**
   * Const accessor for neighbor residual of a particular variable corresponding
   * to the variable index argument.
   */
  const DenseSubVector<Number> & get_neighbor_residual( unsigned int var ) const
  { return *(_neighbor_subresiduals[var]); }

  /**
   * Non-const accessor for neighbor residual of a particular variable corresponding
   * to the variable index argument.
   */
  DenseSubVector<Number> & get_neighbor_residual( unsigned int var )
  { return *(_neighbor_subresiduals[var]); }

  /**
   * Const accessor for element-element Jacobian.
   */
  const DenseMatrix<Number> & get_elem_elem_jacobian() const
  { return _elem_elem_jacobian; }

  /**
   * Non-const accessor for element-element Jacobian.
   */
  DenseMatrix<Number> & get_elem_elem_jacobian()
  { return _elem_elem_jacobian; }

  /**
   * Const accessor for element-element Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  const DenseSubMatrix<Number> & get_elem_elem_jacobian( unsigned int var1, unsigned int var2 ) const
  { return *(_elem_elem_subjacobians[var1][var2]); }

  /**
   * Non-const accessor for element-element Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  DenseSubMatrix<Number> & get_elem_elem_jacobian( unsigned int var1, unsigned int var2 )
  { return *(_elem_elem_subjacobians[var1][var2]); }

  /**
   * Const accessor for element-neighbor Jacobian.
   */
  const DenseMatrix<Number> & get_elem_neighbor_jacobian() const
  { return _elem_neighbor_jacobian; }

  /**
   * Non-const accessor for element -neighborJacobian.
   */
  DenseMatrix<Number> & get_elem_neighbor_jacobian()
  { return _elem_neighbor_jacobian; }

  /**
   * Const accessor for element-neighbor Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  const DenseSubMatrix<Number> & get_elem_neighbor_jacobian( unsigned int var1, unsigned int var2 ) const
  { return *(_elem_neighbor_subjacobians[var1][var2]); }

  /**
   * Non-const accessor for element-neighbor Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  DenseSubMatrix<Number> & get_elem_neighbor_jacobian( unsigned int var1, unsigned int var2 )
  { return *(_elem_neighbor_subjacobians[var1][var2]); }

  /**
   * Const accessor for element-neighbor Jacobian.
   */
  const DenseMatrix<Number> & get_neighbor_elem_jacobian() const
  { return _neighbor_elem_jacobian; }

  /**
   * Non-const accessor for element Jacobian.
   */
  DenseMatrix<Number> & get_neighbor_elem_jacobian()
  { return _neighbor_elem_jacobian; }

  /**
   * Const accessor for neighbor-element Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  const DenseSubMatrix<Number> & get_neighbor_elem_jacobian( unsigned int var1, unsigned int var2 ) const
  { return *(_neighbor_elem_subjacobians[var1][var2]); }

  /**
   * Non-const accessor for neighbor-element Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  DenseSubMatrix<Number> & get_neighbor_elem_jacobian( unsigned int var1, unsigned int var2 )
  { return *(_neighbor_elem_subjacobians[var1][var2]); }

  /**
   * Const accessor for element-neighbor Jacobian.
   */
  const DenseMatrix<Number> & get_neighbor_neighbor_jacobian() const
  { return _neighbor_neighbor_jacobian; }

  /**
   * Non-const accessor for element Jacobian.
   */
  DenseMatrix<Number> & get_neighbor_neighbor_jacobian()
  { return _neighbor_neighbor_jacobian; }

  /**
   * Const accessor for neighbor-neighbor Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  const DenseSubMatrix<Number> & get_neighbor_neighbor_jacobian( unsigned int var1, unsigned int var2 ) const
  { return *(_neighbor_neighbor_subjacobians[var1][var2]); }

  /**
   * Non-const accessor for neighbor-neighbor Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  DenseSubMatrix<Number> & get_neighbor_neighbor_jacobian( unsigned int var1, unsigned int var2 )
  { return *(_neighbor_neighbor_subjacobians[var1][var2]); }

  /**
   * Set the neighbor element which we will use to assemble DG terms.
   *
   * \note We do not assume that this element is get_elem().neighbor(side)
   * because we also need to be able to handle the special case of DG terms on
   * "cracks" in a mesh to model certain types of interface conditions. In this
   * case, we need to be able to specify the neighbor element manually.
   * Also, this should give us more flexibility to handle non-conforming meshes.
   */
  void set_neighbor(const Elem & neighbor)
  { _neighbor = &neighbor; }

  /**
   * Accessor for current neighbor Elem object for assembling DG terms.
   */
  const Elem & get_neighbor() const
  { return *_neighbor; }

  /**
   * Are the DG terms active, i.e. have they been assembled?
   */
  bool dg_terms_are_active() const
  { return _dg_terms_active; }

  /**
   * Accessor for neighbor edge/face (2D/3D) finite element object for variable var.
   */
  template<typename OutputShape>
  void get_neighbor_side_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const;

private:

  /**
   * Current neighbor element for assembling DG terms.
   */
  const Elem * _neighbor;

  /**
   * Residual vector of the neighbor component.
   */
  DenseVector<Number> _neighbor_residual;

  /**
   * The DG Jacobian terms.
   * Trial and test functions come from either element or neighbor.
   */
  DenseMatrix<Number> _elem_elem_jacobian;
  DenseMatrix<Number> _elem_neighbor_jacobian;
  DenseMatrix<Number> _neighbor_elem_jacobian;
  DenseMatrix<Number> _neighbor_neighbor_jacobian;

  /**
   * Element residual subvectors and Jacobian submatrices
   */
  std::vector<std::unique_ptr<DenseSubVector<Number>>> _neighbor_subresiduals;
  std::vector<std::vector<std::unique_ptr<DenseSubMatrix<Number>>>> _elem_elem_subjacobians;
  std::vector<std::vector<std::unique_ptr<DenseSubMatrix<Number>>>> _elem_neighbor_subjacobians;
  std::vector<std::vector<std::unique_ptr<DenseSubMatrix<Number>>>> _neighbor_elem_subjacobians;
  std::vector<std::vector<std::unique_ptr<DenseSubMatrix<Number>>>> _neighbor_neighbor_subjacobians;

  /**
   * Global Degree of freedom index lists for the neighbor element
   */
  std::vector<dof_id_type> _neighbor_dof_indices;
  std::vector<std::vector<dof_id_type>> _neighbor_dof_indices_var;

  /**
   * Finite element objects for each variable's
   * sides on the neighbor element.
   * We do not need FE objects for neighbor element
   * interior since we just need to handle DG interface
   * terms here.
   */
  std::map<FEType, std::unique_ptr<FEAbstract<>>> _neighbor_side_fe;

  /**
   * Pointers to the same finite element objects on the neighbor element,
   * but indexed by variable number
   */
  std::vector<FEAbstract<> *> _neighbor_side_fe_var;

  /**
   * Boolean flag to indicate whether or not the DG terms have been
   * assembled and should be used in the global matrix assembly.
   */
  bool _dg_terms_active;
};

template<typename OutputShape>
inline
void DGFEMContext::get_neighbor_side_fe( unsigned int var, FEGenericBase<OutputShape> *& fe ) const
{
  libmesh_assert_less ( var, _neighbor_side_fe_var.size() );
  fe = cast_ptr<FEGenericBase<OutputShape> *>( _neighbor_side_fe_var[var] );
}

} // namespace libMesh

#endif // LIBMESH_FEM_CONTEXT_H
