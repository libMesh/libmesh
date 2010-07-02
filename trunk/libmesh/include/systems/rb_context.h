// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef __rb_context_h__
#define __rb_context_h__

// C++ includes

// Local Includes
#include "dense_matrix.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "dense_vector.h"
#include "elem.h"
#include "vector_value.h"
#include "fe_type.h"

namespace libMesh
{

// Forward declaration
class RBSystem;
class FEBase;
class QBase;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBContext borrows code liberally from libMesh's DiffContext
 * and FEMContext in order to implement an analogous context
 * class that can be used to faciliate assembly of reduced
 * basis affine operators.
 *
 * @author David J. Knezevic and John W. Peterson, 2009
 */

// ------------------------------------------------------------
// RBContext class definition

class RBContext
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBContext (const RBSystem &);

  /**
   * Destructor.
   */
  virtual ~RBContext ();

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element interior
   */
  Number interior_value(unsigned int var, unsigned int qp);

  /**
   * Returns the value of the solution variable \p var at the quadrature
   * point \p qp on the current element side
   */
  Number side_value(unsigned int var, unsigned int qp);

  /**
   * Returns the value of the solution variable \p var at the physical
   * point \p p on the current element
   */
  Number point_value(unsigned int var, Point &p);

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element interior
   */
  Gradient interior_gradient(unsigned int var, unsigned int qp);

  /**
   * Returns the gradient of the solution variable \p var at the quadrature
   * point \p qp on the current element side
   */
  Gradient side_gradient(unsigned int var, unsigned int qp);

  /**
   * Reinitialize all the context data on a given
   * element for the given system
   */
  virtual void reinit(RBSystem&, Elem*);

  /**
   * Reinitializes interior FE objects on the current geometric element
   */
  void elem_fe_reinit();

  /**
   * Reinitializes side FE objects on the current geometric element
   */
  void elem_side_fe_reinit();

  /**
   * Element by element components of nonlinear_solution
   * as adjusted by a time_solver
   */
  DenseVector<Number> elem_solution;
  std::vector<DenseSubVector<Number> *> elem_subsolutions;

  /**
   * Element residual vector
   */
  DenseVector<Number> elem_vector;

  /**
   * Element jacobian: derivatives of elem_residual with respect to
   * elem_solution
   */
  DenseMatrix<Number> elem_matrix;

  /**
   * Element residual subvectors and Jacobian submatrices
   */
  std::vector<DenseSubVector<Number> *> elem_subvectors;
  std::vector<std::vector<DenseSubMatrix<Number> *> > elem_submatrices;

  /**
   * Global Degree of freedom index lists
   */
  std::vector<unsigned int> dof_indices;
  std::vector<std::vector<unsigned int> > dof_indices_var;

  /**
   * Finite element objects for each variable's interior and sides.
   */
  std::map<FEType, FEBase *> element_fe;
  std::map<FEType, FEBase *> side_fe;

  /**
   * Pointers to the same finite element objects, but indexed
   * by variable number
   */
  std::vector<FEBase *> element_fe_var;
  std::vector<FEBase *> side_fe_var;

  /**
   * Quadrature rules for element interior and sides.
   * The FEM system will try to find a quadrature rule that
   * correctly integrates all variables
   */
  QBase *element_qrule;
  QBase *side_qrule;

  /**
   * Current element for element_* to examine
   */
  Elem *elem;

  /**
   * Current side for side_* to examine
   */
  unsigned char side;

  /**
   * Cached dimension of elements in this mesh
   */
  unsigned char dim;
};

} // namespace libMesh


#endif
