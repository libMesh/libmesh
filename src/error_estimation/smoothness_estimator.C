// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// libmesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/smoothness_estimator.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/error_vector.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/enum_to_string.h"

// C++ includes
#include <algorithm> // for std::fill
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>     // for std::sqrt std::pow std::abs

namespace libMesh
{


//-----------------------------------------------------------------
// SmoothnessEstimator implementations
SmoothnessEstimator::SmoothnessEstimator() :
    _extra_order(1)
{
}



std::vector<Real> SmoothnessEstimator::legepoly(const unsigned int dim,
                                                const Order order,
                                                const Point p,
                                                const unsigned int matsize)
{
  std::vector<Real> psi;
  psi.reserve(matsize);

  // Evaluate 1D Legendre polynomials at x, y, z up to order
  const int npols = static_cast<int>(order) + 1;

  std::vector<Real> Lx(npols, 0.), Ly(npols, 1.), Lz(npols, 1.);
  const Real x = p(0);
  Lx[0] = 1.;
  if (npols > 1)
    Lx[1] = x;
  for (int n = 2; n < npols; ++n)
    Lx[n] = ((2. * n - 1) * x * Lx[n - 1] - (n - 1) * Lx[n - 2]) / n;

  if (dim > 1)
  {
    const Real y = p(1);
    Ly[0] = 1.;
    if (npols > 1)
      Ly[1] = y;
    for (int n = 2; n < npols; ++n)
      Ly[n] = ((2. * n - 1) * y * Ly[n - 1] - (n - 1) * Ly[n - 2]) / n;
  }

  if (dim > 2)
  {
    const Real z = p(2);
    Lz[0] = 1.;
    if (npols > 1)
      Lz[1] = z;
    for (int n = 2; n < npols; ++n)
      Lz[n] = ((2. * n - 1) * z * Lz[n - 1] - (n - 1) * Lz[n - 2]) / n;
  }

  // Loop over total degree and build tensor-product Legendre basis
  for (unsigned int poly_deg = 0; poly_deg <= static_cast<unsigned int>(order); poly_deg++)
  {
    switch (dim)
    {
      case 3:
        for (int i = poly_deg; i >= 0; --i)
          for (int j = poly_deg - i; j >= 0; --j)
          {
            int k = poly_deg - i - j;
            psi.push_back(Lx[i] * Ly[j] * Lz[k]);
          }
        break;

      case 2:
        for (int i = poly_deg; i >= 0; --i)
        {
          int j = poly_deg - i;
          psi.push_back(Lx[i] * Ly[j]);
        }
        break;

      case 1:
        psi.push_back(Lx[poly_deg]);
        break;

      default:
        libmesh_error_msg("Invalid dimension dim " << dim);
    }
  }

  return psi;
}

Real SmoothnessEstimator::compute_slope(int N, Real Sx, Real Sy, Real Sxx, Real Sxy)
{
    const Real denom = (N * Sxx - Sx * Sx);
    // Triggers for first order polynomials when there only 1 point
    // available at log |c_k| and log |k| space to fit.
    if (std::abs(denom) < std::numeric_limits<Real>::epsilon()) {
        return std::numeric_limits<Real>::max();
    }
    return (N * Sxy - Sx * Sy) / denom;
}

void SmoothnessEstimator::reduce_smoothness (std::vector<ErrorVectorReal> & error_per_cell,
                                   const Parallel::Communicator & comm)
{
  // Aggregates element-wise contributions computed
  // in parallel across all processors
  comm.sum(error_per_cell);
}

void SmoothnessEstimator::estimate_smoothness (const System & system,
                                                  ErrorVector & smoothness_per_cell,
                                                  const NumericVector<Number> * solution_vector)
{
  LOG_SCOPE("estimate_error()", "SmoothnessEstimator");

  // The current mesh
  const MeshBase & mesh = system.get_mesh();

  // Resize the smoothness_per_cell vector to be
  // the number of elements, initialize it to 0.
  smoothness_per_cell.resize (mesh.max_elem_id());
  std::fill (smoothness_per_cell.begin(), smoothness_per_cell.end(), 0.);

  // Prepare current_local_solution to localize a non-standard
  // solution vector if necessary
  if (solution_vector && solution_vector != system.solution.get())
    {
      NumericVector<Number> * newsol =
        const_cast<NumericVector<Number> *>(solution_vector);
      System & sys = const_cast<System &>(system);
      newsol->swap(*sys.solution);
      sys.update();
    }

  //------------------------------------------------------------
  // Iterate over all the active elements in the mesh
  // that live on this processor.
  Threads::parallel_for (ConstElemRange(mesh.active_local_elements_begin(),
                                        mesh.active_local_elements_end(),
                                        200),
                         EstimateSmoothness(system,
                                       *this,
                                       smoothness_per_cell)
                         );

  // Each processor has now computed the smoothness for its local
  // elements, and smoothness_per_cell contains 0 for all the
  // non-local elements.  Summing the vector will provide the true
  // value for each element, local or remote
  this->reduce_smoothness (smoothness_per_cell, system.comm());

  // If we used a non-standard solution before, now is the time to fix
  // the current_local_solution
  if (solution_vector && solution_vector != system.solution.get())
    {
      NumericVector<Number> * newsol = const_cast<NumericVector<Number> *>(solution_vector);
      System & sys = const_cast<System &>(system);
      newsol->swap(*sys.solution);
      sys.update();
    }
}


void SmoothnessEstimator::EstimateSmoothness::operator()(const ConstElemRange & range) const
{
  // The current mesh
  const MeshBase & mesh = system.get_mesh();

  // The dimensionality of the mesh
  const unsigned int dim = mesh.mesh_dimension();

  // The number of variables in the system
  const unsigned int n_vars = system.n_vars();

  // The DofMap for this system
  const DofMap & dof_map = system.get_dof_map();

  //------------------------------------------------------------
  // Iterate over all the elements in the range.
for (const auto & elem : range)
{
  const dof_id_type e_id = elem->id();
  // Loop over each variable
  for (unsigned int var=0; var<n_vars; var++)
  {
    const FEType & fe_type = dof_map.variable_type(var);
    const Order element_order = fe_type.order + elem->p_level();

    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
    std::unique_ptr<QBase> qrule(fe_type.default_quadrature_rule(dim, smoothness_estimator._extra_order));
    fe->attach_quadrature_rule(qrule.get());

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<Point> & q_point = fe->get_xyz();
    const std::vector<std::vector<Real>> * phi = &(fe->get_phi());

    std::vector<dof_id_type> dof_indices;

    unsigned int matsize = element_order + 1;
    if (dim > 1)
    {
      matsize *= (element_order + 2);
      matsize /= 2;
    }
    if (dim > 2)
    {
      matsize *= (element_order + 3);
      matsize /= 3;
    }

    DenseMatrix<Number> Kp(matsize, matsize);
    DenseVector<Number> F;
    DenseVector<Number> Pu_h;

    F.resize(matsize);
    Pu_h.resize(matsize);


    fe->reinit(elem);

    dof_map.dof_indices(elem, dof_indices, var);
    libmesh_assert_equal_to(dof_indices.size(), phi->size());

    const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());
    const unsigned int n_qp   = qrule->n_points();

    for (unsigned int qp = 0; qp < n_qp; qp++)
    {
      std::vector<Real> psi = legepoly(dim, element_order, q_point[qp], matsize);
      const unsigned int psi_size = cast_int<unsigned int>(psi.size());

      for (unsigned int i = 0; i < matsize; i++)
        for (unsigned int j = 0; j < matsize; j++)
          Kp(i, j) += JxW[qp] * psi[i] * psi[j];

      Number u_h = libMesh::zero;
      for (unsigned int i = 0; i < n_dofs; i++)
        u_h += (*phi)[i][qp] * system.current_solution(dof_indices[i]);

      for (unsigned int i = 0; i < psi_size; i++)
        F(i) += JxW[qp] * u_h * psi[i];
    }

    Kp.lu_solve(F, Pu_h);

    // Generate index to total degree map. Total_degree_per_index[i] gives the degree ith basis function
    std::vector<unsigned int> total_degree_per_index;

    for (unsigned int poly_deg = 0; poly_deg <= static_cast<unsigned int>(element_order); ++poly_deg)
    {
      switch (dim)
      {
        case 3:
          for (int i = poly_deg; i >= 0; --i)
            for (int j = poly_deg - i; j >= 0; --j)
            {
              int k = poly_deg - i - j;
              total_degree_per_index.push_back(i + j + k);
            }
          break;

        case 2:
          for (int i = poly_deg; i >= 0; --i)
          {
            int j = poly_deg - i;
            total_degree_per_index.push_back(i + j);
          }
          break;

        case 1:
          total_degree_per_index.push_back(poly_deg);
          break;

        default:
          libmesh_error_msg("Invalid dimension dim " << dim);
      }
    }

    // Group coefficients |c_k| by degree k
    std::map<unsigned int, std::vector<Number>> coeff_by_degree;

    for (unsigned int i = 0; i < Pu_h.size(); ++i)
    {
      unsigned int degree = total_degree_per_index[i];
      // Excluding the constant mode i.e, zeroth order coefficient as we use ln|order| later
      // Constant order term often dominates the scale and doesn't inform about regularity.
      if (degree > 0)
        coeff_by_degree[degree].push_back(Pu_h(i));
    }

    // Compute L2 norm of each group |c_k|
    std::vector<unsigned int> degrees;
    std::vector<double> norms;

    for (const auto & pair : coeff_by_degree)
    {
      unsigned int deg = pair.first;
      const std::vector<Number> & coeffs = pair.second;

      double norm = 0.;
      for (const auto & c : coeffs)
        norm += std::norm(c);

      norm = std::sqrt(norm);

      degrees.push_back(deg);
      norms.push_back(norm);
    }


    // Collect log |c_k| and log |k|
    std::vector<double> log_norms, log_degrees;

    for (unsigned int i = 0; i < degrees.size(); ++i)
    {
      // filter tiny terms
      if (norms[i] > 1e-12)
      {
        log_degrees.push_back(std::log(static_cast<double>(degrees[i])));
        log_norms.push_back(std::log(norms[i]));
      }
    }


    // Fit log-log line - we use least-squares fit
    double Sx = 0., Sy = 0., Sxx = 0., Sxy = 0.;
    const size_t N = log_degrees.size();

    for (size_t i = 0; i < N; ++i)
    {
      Sx += log_degrees[i];
      Sy += log_norms[i];
      Sxx += log_degrees[i] * log_degrees[i];
      Sxy += log_degrees[i] * log_norms[i];
    }

    const double regularity = -compute_slope(N, Sx, Sy, Sxx, Sxy);
    // const double intercept = (Sy - slope * Sx) / N;

    Threads::spin_mutex::scoped_lock acquire(Threads::spin_mtx);
    smoothness_per_cell[e_id] = regularity;
  } // end variables loop

} // end element loop

} // End () operator definition

} // namespace libMesh