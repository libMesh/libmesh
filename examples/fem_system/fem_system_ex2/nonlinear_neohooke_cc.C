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

#include "nonlinear_neohooke_cc.h"

void NonlinearNeoHookeCurrentConfig::init_for_qp(VectorValue<Gradient> & grad_u,
                                                 unsigned int qp)
{
  this->current_qp = qp;
  F.zero();
  S.zero();

  {
    RealTensor invF;
    invF.zero();
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        invF(i, j) += libmesh_real(grad_u(i)(j));

    F.add(invF.inverse());

    libmesh_assert_greater (F.det(), -TOLERANCE);
  }

  if (this->calculate_linearized_stiffness)
    this->calculate_tangent();

  this->calculate_stress();
}



void NonlinearNeoHookeCurrentConfig::calculate_tangent()
{
  Real mu = E / (2 * (1 + nu));
  Real lambda = E * nu / ((1 + nu) * (1 - 2 * nu));

  Real detF = F.det();

  C_mat.resize(6, 6);
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      if (i == j) {
        C_mat(i, j) = 2 * mu + lambda;
        C_mat(i + 3, j + 3) = mu - 0.5 * lambda * (detF * detF - 1);
      } else {
        C_mat(i, j) = lambda * detF * detF;
      }
    }
  }
}


void NonlinearNeoHookeCurrentConfig::calculate_stress()
{
  double mu = E / (2.0 * (1.0 + nu));
  double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));

  Real detF = F.det();
  RealTensor Ft = F.transpose();

  RealTensor C = Ft * F;
  RealTensor b = F * Ft;
  RealTensor identity;
  identity(0, 0) = 1.0; identity(1, 1) = 1.0; identity(2, 2) = 1.0;
  RealTensor invC = C.inverse();

  S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC);

  tau = (F * S) * Ft;
  sigma = 1/detF * tau;
}

void NonlinearNeoHookeCurrentConfig::get_residual(DenseVector<Real> & residuum,
                                                  unsigned int & i)
{
  B_L.resize(3, 6);
  DenseVector<Real> sigma_voigt(6);

  this->build_b_0_mat(i, B_L);

  tensor_to_voigt(sigma, sigma_voigt);

  B_L.vector_mult(residuum, sigma_voigt);
}

void NonlinearNeoHookeCurrentConfig::tensor_to_voigt(const RealTensor & tensor,
                                                     DenseVector<Real> & vec)
{
  vec(0) = tensor(0, 0);
  vec(1) = tensor(1, 1);
  vec(2) = tensor(2, 2);
  vec(3) = tensor(0, 1);
  vec(4) = tensor(1, 2);
  vec(5) = tensor(0, 2);

}

void NonlinearNeoHookeCurrentConfig::get_linearized_stiffness(DenseMatrix<Real> & stiffness,
                                                              unsigned int & i,
                                                              unsigned int & j)
{
  stiffness.resize(3, 3);

  double G_IK = (sigma * dphi[i][current_qp]) * dphi[j][current_qp];
  stiffness(0, 0) += G_IK;
  stiffness(1, 1) += G_IK;
  stiffness(2, 2) += G_IK;

  B_L.resize(3, 6);
  this->build_b_0_mat(i, B_L);
  B_K.resize(3, 6);
  this->build_b_0_mat(j, B_K);

  B_L.right_multiply(C_mat);
  B_L.right_multiply_transpose(B_K);
  B_L *= 1/F.det();

  stiffness += B_L;
}

void NonlinearNeoHookeCurrentConfig::build_b_0_mat(int i, DenseMatrix<Real> & b_0_mat)
{
  for (unsigned int ii = 0; ii < 3; ++ii)
    b_0_mat(ii, ii) = dphi[i][current_qp](ii);

  b_0_mat(0, 3) = dphi[i][current_qp](1);
  b_0_mat(1, 3) = dphi[i][current_qp](0);
  b_0_mat(1, 4) = dphi[i][current_qp](2);
  b_0_mat(2, 4) = dphi[i][current_qp](1);
  b_0_mat(0, 5) = dphi[i][current_qp](2);
  b_0_mat(2, 5) = dphi[i][current_qp](0);
}
