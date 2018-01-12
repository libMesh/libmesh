// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// \author Robert Weidlich
// \date Copyright 2012

#ifndef NONLINEAR_NEOHOOKE_CC_H_
#define NONLINEAR_NEOHOOKE_CC_H_

#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/getpot.h"

using namespace libMesh;

/**
 * This class implements a constitutive formulation for an Neo-Hookean elastic solid
 * in terms of the current configuration. This implementation is suitable for computing
 * large deformations. See e.g. "Nonlinear finite element methods" (P. Wriggers, 2008,
 * Springer) for details.
 */
class NonlinearNeoHookeCurrentConfig
{
public:
  NonlinearNeoHookeCurrentConfig(const std::vector<std::vector<RealGradient>> & dphi_in,
                                 GetPot & args,
                                 bool calculate_linearized_stiffness_in) :
    calculate_linearized_stiffness(calculate_linearized_stiffness_in),
    dphi(dphi_in)
  {
    E = args("material/neohooke/e_modulus", 10000.0);
    nu = args("material/neohooke/nu", 0.3);
  }

  /**
   * Initialize the class for the given displacement gradient at the
   * specified quadrature point.
   */
  void init_for_qp(VectorValue<Gradient> & grad_u,
                   unsigned int qp);

  /**
   * Return the residual vector for the current state.
   */
  void get_residual(DenseVector<Real> & residuum,
                    unsigned int & i);

  /**
   * Return the stiffness matrix for the current state.
   */
  void get_linearized_stiffness(DenseMatrix<Real> & stiffness,
                                unsigned int & i,
                                unsigned int & j);

  /**
   * Flag to indicate if it is necessary to calculate values for stiffness
   * matrix during initialization.
   */
  bool calculate_linearized_stiffness;
private:
  void build_b_0_mat(int i, DenseMatrix<Real> & b_l_mat);
  void calculate_stress();
  void calculate_tangent();
  static void tensor_to_voigt(const RealTensor & tensor,
                              DenseVector<Real> & vec);

  unsigned int current_qp;
  const std::vector<std::vector<RealGradient>> & dphi;

  DenseMatrix<Real> C_mat;
  Real E;
  Real nu;
  RealTensor F, S, tau, sigma;
  DenseMatrix<Real> B_L;
  DenseMatrix<Real> B_K;
};

#endif // NONLINEAR_NEOHOOKE_CC_H_
