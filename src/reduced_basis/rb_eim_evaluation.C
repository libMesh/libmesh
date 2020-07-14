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

// C++ includes
#include <sstream>
#include <fstream>

// rbOOmit includes
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_eim_theta.h"
#include "libmesh/rb_parametrized_function.h"

// libMesh includes
#include "libmesh/xdr_cxx.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/elem.h"
#include "timpi/parallel_implementation.h"

namespace libMesh
{

RBEIMEvaluation::RBEIMEvaluation(const Parallel::Communicator & comm)
:
ParallelObject(comm),
evaluate_eim_error_bound(true)
{
}

RBEIMEvaluation::~RBEIMEvaluation()
{
}

void RBEIMEvaluation::clear()
{
  _interpolation_points_xyz.clear();
  _interpolation_points_comp.clear();
  _interpolation_points_subdomain_id.clear();
  _interpolation_points_elem_id.clear();
  _interpolation_points_qp.clear();

  _interpolation_matrix.resize(0,0);

  // Delete any RBTheta objects that were created
  _rb_eim_theta_objects.clear();
}

void RBEIMEvaluation::resize_data_structures(const unsigned int Nmax)
{
  // Resize the data structures relevant to the EIM system
  _interpolation_points_xyz.clear();
  _interpolation_points_comp.clear();
  _interpolation_points_subdomain_id.clear();
  _interpolation_points_elem_id.clear();
  _interpolation_points_qp.clear();

  _interpolation_matrix.resize(Nmax,Nmax);
}

void RBEIMEvaluation::set_parametrized_function(std::unique_ptr<RBParametrizedFunction> pf)
{
  _parametrized_function = std::move(pf);
}

RBParametrizedFunction & RBEIMEvaluation::get_parametrized_function()
{
  if(!_parametrized_function)
    libmesh_error_msg("Parametrized function not initialized yet");

  return *_parametrized_function;
}

Real RBEIMEvaluation::rb_eim_solve(unsigned int N)
{
  LOG_SCOPE("rb_eim_solve()", "RBEIMEvaluation");

  if (N > get_n_basis_functions())
    libmesh_error_msg("Error: N cannot be larger than the number of basis functions in rb_solve");

  if (N==0)
    libmesh_error_msg("Error: N must be greater than 0 in rb_solve");

  // Get the rhs by sampling parametrized_function
  // at the first N interpolation_points
  DenseVector<Number> EIM_rhs(N);
  for (unsigned int i=0; i<N; i++)
    {
      EIM_rhs(i) = get_parametrized_function().evaluate(get_parameters(),
                                                        _interpolation_points_comp[i],
                                                        _interpolation_points_xyz[i],
                                                        _interpolation_points_subdomain_id[i]);
    }

  DenseMatrix<Number> interpolation_matrix_N;
  _interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  interpolation_matrix_N.lu_solve(EIM_rhs, _rb_eim_solution);

  // Optionally evaluate an a posteriori error bound. The EIM error estimate
  // recommended in the literature is based on using "next" EIM point, so
  // we skip this if N == get_n_basis_functions()
  if (evaluate_eim_error_bound && (N != get_n_basis_functions()))
    {
      // Compute the a posteriori error bound
      // First, sample the parametrized function at x_{N+1}
      Number g_at_next_x = get_parametrized_function().evaluate(get_parameters(),
                                                        _interpolation_points_comp[N],
                                                        _interpolation_points_xyz[N],
                                                        _interpolation_points_subdomain_id[N]);

      // Next, evaluate the EIM approximation at x_{N+1}
      Number EIM_approx_at_next_x = 0.;
      for (unsigned int j=0; j<N; j++)
        {
          EIM_approx_at_next_x += _rb_eim_solution(j) * _interpolation_matrix(N,j);
        }

      Real error_estimate = std::abs(g_at_next_x - EIM_approx_at_next_x);
      return error_estimate;
    }
  else // Don't evaluate an error bound
    {
      return -1.;
    }
}

void RBEIMEvaluation::rb_eim_solve(DenseVector<Number> & EIM_rhs)
{
  LOG_SCOPE("rb_eim_solve()", "RBEIMEvaluation");

  if (EIM_rhs.size() > get_n_basis_functions())
    libmesh_error_msg("Error: N cannot be larger than the number of basis functions in rb_solve");

  if (EIM_rhs.size()==0)
    libmesh_error_msg("Error: N must be greater than 0 in rb_solve");

  const unsigned int N = EIM_rhs.size();
  DenseMatrix<Number> interpolation_matrix_N;
  _interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  interpolation_matrix_N.lu_solve(EIM_rhs, _rb_eim_solution);
}

unsigned int RBEIMEvaluation::get_n_basis_functions() const
{
  return _local_eim_basis_functions.size();
}

void RBEIMEvaluation::set_n_basis_functions(unsigned int n_bfs)
{
  _local_eim_basis_functions.resize(n_bfs);
}

void RBEIMEvaluation::decrement_vector(std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & v,
                                       const DenseVector<Number> & coeffs)
{
  if(get_n_basis_functions() != coeffs.size())
    libmesh_error_msg("Error: Number of coefficients should match number of basis functions");

  for (const auto v_it : v)
    {
      dof_id_type elem_id = v_it.first;
      const auto & v_comp_and_qp = v_it.second;
      
      for (const auto & comp : index_range(v_comp_and_qp))
        for (unsigned int qp : index_range(v_comp_and_qp[comp]))
          for (unsigned int i : index_range(_local_eim_basis_functions))
            {
              // Check that entry (elem_id,comp,qp) exists in _local_eim_basis_functions so that
              // we get a clear error message if there is any missing data
              auto basis_it = _local_eim_basis_functions[i].find(elem_id);
              if (basis_it == _local_eim_basis_functions[i].end())
                libmesh_error_msg("Error: Missing elem_id");

              const auto & basis_comp_and_qp = basis_it->second;
              if(comp >= basis_comp_and_qp.size())
                libmesh_error_msg("Error: Invalid comp");
              if(qp >= basis_comp_and_qp[comp].size())
                libmesh_error_msg("Error: Invalid qp");

              v[elem_id][comp][qp] -= get_rb_eim_solution()(i) * basis_comp_and_qp[comp][qp];
            }
    }

}

void RBEIMEvaluation::initialize_eim_theta_objects()
{
  // Initialize the rb_theta objects that access the solution from this rb_eim_evaluation
  _rb_eim_theta_objects.clear();
  for (unsigned int i=0; i<get_n_basis_functions(); i++)
    _rb_eim_theta_objects.emplace_back(build_eim_theta(i));
}

std::vector<std::unique_ptr<RBTheta>> & RBEIMEvaluation::get_eim_theta_objects()
{
  return _rb_eim_theta_objects;
}

std::unique_ptr<RBTheta> RBEIMEvaluation::build_eim_theta(unsigned int index)
{
  return libmesh_make_unique<RBEIMTheta>(*this, index);
}

void RBEIMEvaluation::get_parametrized_function_values_at_qps(
  const std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & pf,
  dof_id_type elem_id,
  unsigned int comp,
  std::vector<Number> & values)
{
  LOG_SCOPE("get_parametrized_function_values_at_qps()", "RBEIMConstruction");

  values.clear();

  const auto it = pf.find(elem_id);
  if(it != pf.end())
  {
    const auto & comps_and_qps_on_elem = it->second;
    if(comp >= comps_and_qps_on_elem.size())
    {
      libmesh_error_msg("Invalid comp index: " + std::to_string(comp));
    }

    values = comps_and_qps_on_elem[comp];
  }
}

Number RBEIMEvaluation::get_parametrized_function_value(
  const Parallel::Communicator & comm,
  const std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & pf,
  dof_id_type elem_id,
  unsigned int comp,
  unsigned int qp)
{
  std::vector<Number> values;
  get_parametrized_function_values_at_qps(pf, elem_id, comp, values);

  // In parallel, values should only be non-empty on one processor
  Number value = 0.;
  if(!values.empty())
  {
    if(qp >= values.size())
      libmesh_error_msg("Error: Invalid qp index");

    value = values[qp];
  }
  comm.sum(value);

  return value;
}

void RBEIMEvaluation::get_eim_basis_function_values_at_qps(unsigned int basis_function_index,
                                                           dof_id_type elem_id,
                                                           unsigned int comp,
                                                           std::vector<Number> & values) const
{
  if(basis_function_index >= _local_eim_basis_functions.size())
  {
    libmesh_error_msg("Invalid basis function index: " + std::to_string(basis_function_index));
  }

  get_parametrized_function_values_at_qps(
    _local_eim_basis_functions[basis_function_index],
    elem_id,
    comp,
    values);
}

Number RBEIMEvaluation::get_eim_basis_function_value(unsigned int basis_function_index,
                                                     dof_id_type elem_id,
                                                     unsigned int comp,
                                                     unsigned int qp) const
{
  if(basis_function_index >= _local_eim_basis_functions.size())
  {
    libmesh_error_msg("Invalid basis function index: " + std::to_string(basis_function_index));
  }

  return get_parametrized_function_value(
    comm(),
    _local_eim_basis_functions[basis_function_index],
    elem_id,
    comp,
    qp);
}

const std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> &
  RBEIMEvaluation::get_basis_function(unsigned int i) const
{
  return _local_eim_basis_functions[i];
}

const DenseVector<Number> & RBEIMEvaluation::get_rb_eim_solution() const
{
  return _rb_eim_solution;
}

void RBEIMEvaluation::add_interpolation_points_xyz(Point p)
{
  _interpolation_points_xyz.emplace_back(p);
}

void RBEIMEvaluation::add_interpolation_points_comp(unsigned int comp)
{
  _interpolation_points_comp.emplace_back(comp);
}

void RBEIMEvaluation::add_interpolation_points_subdomain_id(subdomain_id_type sbd_id)
{
  _interpolation_points_subdomain_id.emplace_back(sbd_id);
}

void RBEIMEvaluation::add_interpolation_points_elem_id(dof_id_type elem_id)
{
  _interpolation_points_elem_id.emplace_back(elem_id);
}

void RBEIMEvaluation::add_interpolation_points_qp(unsigned int qp)
{
  _interpolation_points_qp.emplace_back(qp);
}

Point RBEIMEvaluation::get_interpolation_points_xyz(unsigned int index) const
{
  if(index >= _interpolation_points_xyz.size())
    libmesh_error_msg("Error: Invalid index");

  return _interpolation_points_xyz[index];
}

unsigned int RBEIMEvaluation::get_interpolation_points_comp(unsigned int index) const
{
  if(index >= _interpolation_points_comp.size())
    libmesh_error_msg("Error: Invalid index");

  return _interpolation_points_comp[index];
}

subdomain_id_type RBEIMEvaluation::get_interpolation_points_subdomain_id(unsigned int index) const
{
  if(index >= _interpolation_points_subdomain_id.size())
    libmesh_error_msg("Error: Invalid index");

  return _interpolation_points_subdomain_id[index];
}

dof_id_type RBEIMEvaluation::get_interpolation_points_elem_id(unsigned int index) const
{
  if(index >= _interpolation_points_elem_id.size())
    libmesh_error_msg("Error: Invalid index");

  return _interpolation_points_elem_id[index];
}

unsigned int RBEIMEvaluation::get_interpolation_points_qp(unsigned int index) const
{
  if(index >= _interpolation_points_qp.size())
    libmesh_error_msg("Error: Invalid index");

  return _interpolation_points_qp[index];
}

void RBEIMEvaluation::set_interpolation_matrix_entry(unsigned int i, unsigned int j, Number value)
{
  if( (i >= _interpolation_matrix.m()) || (j >= _interpolation_matrix.n()) )
    libmesh_error_msg("Error: Invalid matrix indices");

  _interpolation_matrix(i,j) = value;
}

const DenseMatrix<Number> & RBEIMEvaluation::get_interpolation_matrix() const
{
  return _interpolation_matrix;
}

void RBEIMEvaluation::add_basis_function_and_interpolation_data(
  const std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & bf,
  Point p,
  unsigned int comp,
  dof_id_type elem_id,
  subdomain_id_type subdomain_id,
  unsigned int qp)
{
  _local_eim_basis_functions.emplace_back(bf);

  _interpolation_points_xyz.emplace_back(p);
  _interpolation_points_comp.emplace_back(comp);
  _interpolation_points_elem_id.emplace_back(elem_id);
  _interpolation_points_subdomain_id.emplace_back(subdomain_id);
  _interpolation_points_qp.emplace_back(qp);
}

}
