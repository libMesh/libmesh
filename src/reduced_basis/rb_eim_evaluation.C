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

// rbOOmit includes
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_eim_theta.h"
#include "libmesh/rb_parametrized_function.h"
#include "libmesh/utility.h" // Utility::mkdir

// libMesh includes
#include "libmesh/xdr_cxx.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/system.h"
#include "timpi/parallel_implementation.h"

// C++ includes
#include <sstream>
#include <fstream>
#include <numeric> // std::accumulate
#include <iterator> // std::advance

namespace libMesh
{

RBEIMEvaluation::RBEIMEvaluation(const Parallel::Communicator & comm)
:
ParallelObject(comm),
evaluate_eim_error_bound(true),
_rb_eim_solves_N(0)
{
}

RBEIMEvaluation::~RBEIMEvaluation() = default;

void RBEIMEvaluation::clear()
{
  _interpolation_points_xyz.clear();
  _interpolation_points_comp.clear();
  _interpolation_points_subdomain_id.clear();
  _interpolation_points_xyz_perturbations.clear();
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
  _interpolation_points_xyz_perturbations.clear();
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
  libmesh_error_msg_if(!_parametrized_function, "Parametrized function not initialized yet");

  return *_parametrized_function;
}

Real RBEIMEvaluation::rb_eim_solve(unsigned int N)
{
  LOG_SCOPE("rb_eim_solve()", "RBEIMEvaluation");

  libmesh_error_msg_if(N > get_n_basis_functions(), "Error: N cannot be larger than the number of basis functions in rb_solve");
  libmesh_error_msg_if(N==0, "Error: N must be greater than 0 in rb_solve");

  // Get the rhs by sampling parametrized_function
  // at the first N interpolation_points
  DenseVector<Number> EIM_rhs(N);
  for (unsigned int i=0; i<N; i++)
    {
      EIM_rhs(i) =
        get_parametrized_function().evaluate_comp(get_parameters(),
                                                  _interpolation_points_comp[i],
                                                  _interpolation_points_xyz[i],
                                                  _interpolation_points_subdomain_id[i],
                                                  _interpolation_points_xyz_perturbations[i]);
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
      Number g_at_next_x =
        get_parametrized_function().evaluate_comp(get_parameters(),
                                                  _interpolation_points_comp[N],
                                                  _interpolation_points_xyz[N],
                                                  _interpolation_points_subdomain_id[N],
                                                  _interpolation_points_xyz_perturbations[N]);

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

  libmesh_error_msg_if(EIM_rhs.size() > get_n_basis_functions(),
                       "Error: N cannot be larger than the number of basis functions in rb_solve");

  libmesh_error_msg_if(EIM_rhs.size()==0, "Error: N must be greater than 0 in rb_solve");

  const unsigned int N = EIM_rhs.size();
  DenseMatrix<Number> interpolation_matrix_N;
  _interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  interpolation_matrix_N.lu_solve(EIM_rhs, _rb_eim_solution);
}

void RBEIMEvaluation::rb_eim_solves(const std::vector<RBParameters> & mus,
                                    unsigned int N)
{
  libmesh_error_msg_if(N > get_n_basis_functions(),
    "Error: N cannot be larger than the number of basis functions in rb_eim_solves");
  libmesh_error_msg_if(N==0, "Error: N must be greater than 0 in rb_eim_solves");

  // If mus and N are the same as before, then we return early
  if ((_rb_eim_solves_mus == mus) && (_rb_eim_solves_N == N))
    return;

  LOG_SCOPE("rb_eim_solves()", "RBEIMEvaluation");

  _rb_eim_solves_mus = mus;
  _rb_eim_solves_N = N;

  if (get_parametrized_function().is_lookup_table)
    {
      _rb_eim_solutions.resize(mus.size());
      for (auto mu_index : index_range(mus))
        {
          Real lookup_table_param =
            mus[mu_index].get_value(get_parametrized_function().lookup_table_param_name);

          // Cast lookup_table_param to an unsigned integer so that we can use
          // it as an index into the EIM rhs values obtained from the lookup table.
          unsigned int lookup_table_index =
            cast_int<unsigned int>(std::round(lookup_table_param));

          DenseVector<Number> values;
          eim_solutions[lookup_table_index].get_principal_subvector(N, values);
          _rb_eim_solutions[mu_index] = values;
        }

      return;
    }

  // output all comps indexing is as follows:
  //   mu index --> interpolation point index --> component index --> value.
  std::vector<std::vector<std::vector<Number>>> output_all_comps;
  get_parametrized_function().vectorized_evaluate(mus,
                                                  _interpolation_points_xyz,
                                                  _interpolation_points_subdomain_id,
                                                  _interpolation_points_xyz_perturbations,
                                                  output_all_comps);

  std::vector<std::vector<Number>> evaluated_values_at_interp_points(output_all_comps.size());

  for(unsigned int mu_index : index_range(evaluated_values_at_interp_points))
    {
      evaluated_values_at_interp_points[mu_index].resize(N);
      for(unsigned int interp_pt_index=0; interp_pt_index<N; interp_pt_index++)
        {
          unsigned int comp = _interpolation_points_comp[interp_pt_index];

          evaluated_values_at_interp_points[mu_index][interp_pt_index] =
            output_all_comps[mu_index][interp_pt_index][comp];
        }
    }

  DenseMatrix<Number> interpolation_matrix_N;
  _interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  _rb_eim_solutions.resize(mus.size());
  for(unsigned int mu_index : index_range(mus))
    {
      DenseVector<Number> EIM_rhs(N);
      for (unsigned int i=0; i<N; i++)
        {
          EIM_rhs(i) = evaluated_values_at_interp_points[mu_index][i];
        }

      interpolation_matrix_N.lu_solve(EIM_rhs, _rb_eim_solutions[mu_index]);
    }
}

unsigned int RBEIMEvaluation::get_n_basis_functions() const
{
  return _local_eim_basis_functions.size();
}

void RBEIMEvaluation::set_n_basis_functions(unsigned int n_bfs)
{
  _local_eim_basis_functions.resize(n_bfs);
}

void RBEIMEvaluation::decrement_vector(QpDataMap & v,
                                       const DenseVector<Number> & coeffs)
{
  LOG_SCOPE("decrement_vector()", "RBEIMEvaluation");

  libmesh_error_msg_if(get_n_basis_functions() != coeffs.size(),
                       "Error: Number of coefficients should match number of basis functions");

  for (auto & pr : v)
    {
      dof_id_type elem_id = pr.first;
      auto & v_comp_and_qp = pr.second;

      for (const auto & comp : index_range(v_comp_and_qp))
        for (unsigned int qp : index_range(v_comp_and_qp[comp]))
          for (unsigned int i : index_range(_local_eim_basis_functions))
            {
              // Check that entry (elem_id,comp,qp) exists in _local_eim_basis_functions so that
              // we get a clear error message if there is any missing data
              const auto & basis_comp_and_qp = libmesh_map_find(_local_eim_basis_functions[i], elem_id);

              libmesh_error_msg_if(comp >= basis_comp_and_qp.size(), "Error: Invalid comp");
              libmesh_error_msg_if(qp >= basis_comp_and_qp[comp].size(), "Error: Invalid qp");

              v_comp_and_qp[comp][qp] -= coeffs(i) * basis_comp_and_qp[comp][qp];
            }
    }

}

void RBEIMEvaluation::initialize_eim_theta_objects()
{
  // Initialize the rb_theta objects that access the solution from this rb_eim_evaluation
  _rb_eim_theta_objects.clear();
  for (auto i : make_range(get_n_basis_functions()))
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
  const QpDataMap & pf,
  dof_id_type elem_id,
  unsigned int comp,
  std::vector<Number> & values)
{
  LOG_SCOPE("get_parametrized_function_values_at_qps()", "RBEIMConstruction");

  values.clear();

  const auto it = pf.find(elem_id);
  if (it != pf.end())
  {
    const auto & comps_and_qps_on_elem = it->second;
    libmesh_error_msg_if(comp >= comps_and_qps_on_elem.size(),
                         "Invalid comp index: " << comp);

    values = comps_and_qps_on_elem[comp];
  }
}

Number RBEIMEvaluation::get_parametrized_function_value(
  const Parallel::Communicator & comm,
  const QpDataMap & pf,
  dof_id_type elem_id,
  unsigned int comp,
  unsigned int qp)
{
  std::vector<Number> values;
  get_parametrized_function_values_at_qps(pf, elem_id, comp, values);

  // In parallel, values should only be non-empty on one processor
  Number value = 0.;
  if (!values.empty())
  {
    libmesh_error_msg_if(qp >= values.size(), "Error: Invalid qp index");

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
  libmesh_error_msg_if(basis_function_index >= _local_eim_basis_functions.size(),
                       "Invalid basis function index: " << basis_function_index);

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
  libmesh_error_msg_if(basis_function_index >= _local_eim_basis_functions.size(),
                       "Invalid basis function index: " + basis_function_index);

  return get_parametrized_function_value(
    comm(),
    _local_eim_basis_functions[basis_function_index],
    elem_id,
    comp,
    qp);
}

const RBEIMEvaluation::QpDataMap &
RBEIMEvaluation::get_basis_function(unsigned int i) const
{
  return _local_eim_basis_functions[i];
}

const DenseVector<Number> & RBEIMEvaluation::get_rb_eim_solution() const
{
  return _rb_eim_solution;
}

std::vector<Number> RBEIMEvaluation::get_rb_eim_solutions_entries(unsigned int index) const
{
  LOG_SCOPE("get_rb_eim_solutions_entries()", "RBEIMEvaluation");

  std::vector<Number> rb_eim_solutions_entries(_rb_eim_solutions.size());
  for (unsigned int mu_index : index_range(_rb_eim_solutions))
    {
      libmesh_error_msg_if(index >= _rb_eim_solutions[mu_index].size(), "Error: Invalid index");
      rb_eim_solutions_entries[mu_index] = _rb_eim_solutions[mu_index](index);
    }

  return rb_eim_solutions_entries;
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

void RBEIMEvaluation::add_interpolation_points_xyz_perturbations(const std::vector<Point> & perturbs)
{
  _interpolation_points_xyz_perturbations.emplace_back(perturbs);
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
  libmesh_error_msg_if(index >= _interpolation_points_xyz.size(), "Error: Invalid index");

  return _interpolation_points_xyz[index];
}

unsigned int RBEIMEvaluation::get_interpolation_points_comp(unsigned int index) const
{
  libmesh_error_msg_if(index >= _interpolation_points_comp.size(), "Error: Invalid index");

  return _interpolation_points_comp[index];
}

subdomain_id_type RBEIMEvaluation::get_interpolation_points_subdomain_id(unsigned int index) const
{
  libmesh_error_msg_if(index >= _interpolation_points_subdomain_id.size(), "Error: Invalid index");

  return _interpolation_points_subdomain_id[index];
}

const std::vector<Point> & RBEIMEvaluation::get_interpolation_points_xyz_perturbations(unsigned int index) const
{
  libmesh_error_msg_if(index >= _interpolation_points_xyz_perturbations.size(), "Error: Invalid index");

  return _interpolation_points_xyz_perturbations[index];
}

dof_id_type RBEIMEvaluation::get_interpolation_points_elem_id(unsigned int index) const
{
  libmesh_error_msg_if(index >= _interpolation_points_elem_id.size(), "Error: Invalid index");

  return _interpolation_points_elem_id[index];
}

unsigned int RBEIMEvaluation::get_interpolation_points_qp(unsigned int index) const
{
  libmesh_error_msg_if(index >= _interpolation_points_qp.size(), "Error: Invalid index");

  return _interpolation_points_qp[index];
}

void RBEIMEvaluation::set_interpolation_matrix_entry(unsigned int i, unsigned int j, Number value)
{
  libmesh_error_msg_if((i >= _interpolation_matrix.m()) || (j >= _interpolation_matrix.n()),
                       "Error: Invalid matrix indices");

  _interpolation_matrix(i,j) = value;
}

const DenseMatrix<Number> & RBEIMEvaluation::get_interpolation_matrix() const
{
  return _interpolation_matrix;
}

void RBEIMEvaluation::add_basis_function_and_interpolation_data(
  const QpDataMap & bf,
  Point p,
  unsigned int comp,
  dof_id_type elem_id,
  subdomain_id_type subdomain_id,
  unsigned int qp,
  const std::vector<Point> & perturbs)
{
  _local_eim_basis_functions.emplace_back(bf);

  _interpolation_points_xyz.emplace_back(p);
  _interpolation_points_comp.emplace_back(comp);
  _interpolation_points_elem_id.emplace_back(elem_id);
  _interpolation_points_subdomain_id.emplace_back(subdomain_id);
  _interpolation_points_qp.emplace_back(qp);
  _interpolation_points_xyz_perturbations.emplace_back(perturbs);
}

void RBEIMEvaluation::
write_out_basis_functions(const std::string & directory_name,
                          bool write_binary_basis_functions)
{
  LOG_SCOPE("write_out_basis_functions()", "RBEIMEvaluation");

  // Quick return if there is no work to do. Note: make sure all procs
  // agree there is no work to do.
  bool is_empty = _local_eim_basis_functions.empty();
  this->comm().verify(is_empty);

  if (is_empty)
    return;

  // Gather basis function data from other procs, storing it in
  // _local_eim_basis_functions, so that we can then print everything
  // from processor 0.
  this->gather_bfs();

  // Write values from processor 0 only.
  if (this->processor_id() == 0)
    {
      // Make a directory to store all the data files
      Utility::mkdir(directory_name.c_str());

      // Create filename
      std::ostringstream file_name;
      const std::string basis_function_suffix = (write_binary_basis_functions ? ".xdr" : ".dat");
      file_name << directory_name << "/" << "bf_data" << basis_function_suffix;

      // Create XDR writer object
      Xdr xdr(file_name.str(), write_binary_basis_functions ? ENCODE : WRITE);

      // Write number of basis functions to file. Note: the
      // Xdr::data() function takes non-const references, so you can't
      // pass e.g. vec.size() to that interface.
      auto n_bf = _local_eim_basis_functions.size();
      xdr.data(n_bf, "# Number of basis functions");

      // We assume that each basis function has data for the same
      // number of elements as basis function 0, which is equal to the
      // size of the map.
      auto n_elem = _local_eim_basis_functions[0].size();
      xdr.data(n_elem, "# Number of elements");

      // We assume that each element has the same number of variables,
      // and we get the number of vars from the first element of the
      // first basis function.
      auto n_vars = _local_eim_basis_functions[0].begin()->second.size();
      xdr.data(n_vars, "# Number of variables");

      // We assume that the list of elements for each basis function
      // is the same as basis function 0. We also assume that all vars
      // have the same number of qps.
      std::vector<unsigned int> n_qp_per_elem;
      n_qp_per_elem.reserve(n_elem);
      dof_id_type expected_elem_id = 0;
      for (const auto & pr : _local_eim_basis_functions[0])
        {
          // Note: Currently we require that the Elems are numbered
          // contiguously from [0..n_elem).  This allows us to avoid
          // writing the Elem ids to the Xdr file, but if we need to
          // generalize this assumption later, we can.
          const auto & actual_elem_id = pr.first;

          libmesh_error_msg_if(actual_elem_id != expected_elem_id++,
                               "RBEIMEvaluation currently assumes a contiguous Elem numbering starting from 0.");

          // array[n_vars][n_qp] per Elem. We get the number of QPs
          // for variable 0, assuming they are all the same.
          const auto & array = pr.second;
          n_qp_per_elem.push_back(array[0].size());
        }
      xdr.data(n_qp_per_elem, "# Number of QPs per Elem");

      // The total amount of qp data for each var is the sum of the
      // entries in the "n_qp_per_elem" array.
      auto n_qp_data =
        std::accumulate(n_qp_per_elem.begin(),
                        n_qp_per_elem.end(),
                        0u);

      // Reserve space to store continguous vectors of qp data for each var
      std::vector<std::vector<Number>> qp_data(n_vars);
      for (auto var : index_range(qp_data))
        qp_data[var].reserve(n_qp_data);

      // Now we construct a vector for each basis function, for each
      // variable which is ordered according to:
      // [ [qp vals for Elem 0], [qp vals for Elem 1], ... [qp vals for Elem N] ]
      // and write it to file.
      for (auto bf : index_range(_local_eim_basis_functions))
        {
          // Clear any data from previous bf
          for (auto var : index_range(qp_data))
            qp_data[var].clear();

          for (const auto & pr : _local_eim_basis_functions[bf])
            {
              // array[n_vars][n_qp] per Elem
              const auto & array = pr.second;
              for (auto var : index_range(array))
                {
                  // Insert all qp values for this var
                  qp_data[var].insert(/*insert at*/qp_data[var].end(),
                                      /*data start*/array[var].begin(),
                                      /*data end*/array[var].end());
                }
            }

          // Write all the var values for this bf
          for (auto var : index_range(qp_data))
            xdr.data_stream(qp_data[var].data(), qp_data[var].size(), /*line_break=*/qp_data[var].size());
        }
    }
}

void RBEIMEvaluation::
read_in_basis_functions(const System & sys,
                        const std::string & directory_name,
                        bool read_binary_basis_functions)
{
  LOG_SCOPE("read_in_basis_functions()", "RBEIMEvaluation");

  // Read values on processor 0 only.
  if (sys.comm().rank() == 0)
    {
      // Create filename
      std::ostringstream file_name;
      const std::string basis_function_suffix = (read_binary_basis_functions ? ".xdr" : ".dat");
      file_name << directory_name << "/" << "bf_data" << basis_function_suffix;

      // Create XDR reader object
      Xdr xdr(file_name.str(), read_binary_basis_functions ? DECODE : READ);

      // Read in the number of basis functions. The comment parameter
      // is ignored when reading.
      std::size_t n_bf;
      xdr.data(n_bf);

      // Read in the number of elements
      std::size_t n_elem;
      xdr.data(n_elem);

      // Read in the number of variables.
      std::size_t n_vars;
      xdr.data(n_vars);

      // Read in vector containing the number of QPs per elem. We can
      // create this vector with the required size or let it be read
      // from the file and sized for us.
      std::vector<unsigned int> n_qp_per_elem(n_elem);
      xdr.data(n_qp_per_elem);

      // The total amount of qp data for each var is the sum of the
      // entries in the "n_qp_per_elem" array.
      auto n_qp_data =
        std::accumulate(n_qp_per_elem.begin(),
                        n_qp_per_elem.end(),
                        0u);

      // Allocate space to store all required basis functions,
      // clearing any data that may have been there previously.
      //
      // TODO: Do we need to also write out/read in Elem ids?
      // Or can we assume they will always be contiguously
      // numbered (at least on proc 0)?
      _local_eim_basis_functions.clear();
      _local_eim_basis_functions.resize(n_bf);
      for (auto i : index_range(_local_eim_basis_functions))
        for (std::size_t elem_id=0; elem_id<n_elem; ++elem_id)
          {
            auto & array = _local_eim_basis_functions[i][elem_id];
            array.resize(n_vars);
          }

      // Allocate temporary storage for one var's worth of qp data.
      std::vector<Number> qp_data;

      // Read in data for each basis function
      for (auto i : index_range(_local_eim_basis_functions))
        {
          // Reference to the data map for the current basis function.
          auto & bf_map = _local_eim_basis_functions[i];

          for (std::size_t var=0; var<n_vars; ++var)
            {
              qp_data.clear();
              qp_data.resize(n_qp_data);

              // Read data using data_stream() since that is
              // (currently) how we write it out. The "line_break"
              // parameter of data_stream() is ignored while reading.
              xdr.data_stream(qp_data.data(), qp_data.size());

              // Iterate over the qp_data vector, filling in the
              // "small" vectors for each Elem.
              auto cursor = qp_data.begin();
              for (std::size_t elem_id=0; elem_id<n_elem; ++elem_id)
                {
                  // Get reference to the [n_vars][n_qp] array for
                  // this Elem. We assign() into the vector of
                  // quadrature point values, which allocates space if
                  // it doesn't already exist.
                  auto & array = bf_map[elem_id];
                  array[var].assign(cursor, cursor + n_qp_per_elem[elem_id]);
                  std::advance(cursor, n_qp_per_elem[elem_id]);
                }
            } // end for (var)
        } // end for (i)
    } // end if processor 0

  // Distribute the basis function information to the processors that require it
  this->distribute_bfs(sys);
}

void RBEIMEvaluation::print_local_eim_basis_functions() const
{
  for (auto bf : index_range(_local_eim_basis_functions))
    {
      libMesh::out << "Basis function " << bf << std::endl;
      for (const auto & pr : _local_eim_basis_functions[bf])
        {
          libMesh::out << "Elem " << pr.first << std::endl;
          const auto & array = pr.second;
          for (auto var : index_range(array))
            {
              libMesh::out << "Variable " << var << std::endl;
              for (auto qp : index_range(array[var]))
                libMesh::out << array[var][qp] << " ";
              libMesh::out << std::endl;
            }
        }
    }
}

void RBEIMEvaluation::gather_bfs()
{
  // We need to gather _local_eim_basis_functions data from other
  // procs for printing.
  //
  // Ideally, this could be accomplished by simply calling:
  // this->comm().gather(/*root_id=*/0, _local_eim_basis_functions);
  //
  // but the data structure seems to be too complicated for this to
  // work automatically. (I get some error about the function called
  // being "private within this context".) Therefore, we have to
  // gather the information manually.

  // So we can avoid calling this many times below
  auto n_procs = this->n_processors();

  // In serial there's nothing to gather
  if (n_procs == 1)
    return;

  // Current assumption is that the number of basis functions stored on
  // each processor is the same, the only thing that differs is the number
  // of elements, so make sure that is the case now.
  auto n_bf = _local_eim_basis_functions.size();
  this->comm().verify(n_bf);

  // This function should never be called if there are no basis
  // functions, so if it was, something went wrong.
  libmesh_error_msg_if(!n_bf, "RBEIMEvaluation::gather_bfs() should not be called with 0 basis functions.");

  // The number of variables should be the same on all processors
  // and we can get this from _local_eim_basis_functions. However,
  // it may be that some processors have no local elements, so on
  // those processors we cannot look up the size from
  // _local_eim_basis_functions. As a result we use comm().max(n_vars)
  // to make sure all processors agree on the final value.
  std::size_t n_vars =
    _local_eim_basis_functions[0].empty() ? 0 : _local_eim_basis_functions[0].begin()->second.size();
  this->comm().max(n_vars);

  // Gather list of Elem ids stored on each processor to proc 0.  We
  // use basis function 0 as an example and assume all the basis
  // functions are distributed similarly.
  std::vector<dof_id_type> elem_ids;
  elem_ids.reserve(_local_eim_basis_functions[0].size());
  for (const auto & pr : _local_eim_basis_functions[0])
    elem_ids.push_back(pr.first);
  this->comm().gather(/*root_id=*/0, elem_ids);

  // Store the number of qps per Elem on this processor. Again, use
  // basis function 0 (and variable 0) to get this information, then
  // apply it to all basis functions.
  std::vector<unsigned int> n_qp_per_elem;
  n_qp_per_elem.reserve(_local_eim_basis_functions[0].size());
  for (const auto & pr : _local_eim_basis_functions[0])
    {
      // array[n_vars][n_qp] per Elem. We get the number of QPs
      // for variable 0, assuming they are all the same.
      const auto & array = pr.second;
      n_qp_per_elem.push_back(array[0].size());
    }

  // Before gathering, compute the total amount of local qp data for
  // each var, which is the sum of the entries in the "n_qp_per_elem" array.
  // This will be used to reserve space in a vector below.
  auto n_local_qp_data =
    std::accumulate(n_qp_per_elem.begin(),
                    n_qp_per_elem.end(),
                    0u);

  // Gather the number of qps per Elem for each processor onto processor 0.
  this->comm().gather(/*root_id=*/0, n_qp_per_elem);

  // Sanity check: On processor 0, this checks that we have gathered the same number
  // of elem ids and qp counts.
  libmesh_error_msg_if(elem_ids.size() != n_qp_per_elem.size(),
                       "Must gather same number of Elem ids as qps per Elem.");

  // Reserve space to store contiguous vectors of qp data for each var
  std::vector<std::vector<Number>> gathered_qp_data(n_vars);
  for (auto var : index_range(gathered_qp_data))
    gathered_qp_data[var].reserve(n_local_qp_data);

  // Now we construct a vector for each basis function, for each
  // variable, which is ordered according to:
  // [ [qp vals for Elem 0], [qp vals for Elem 1], ... [qp vals for Elem N] ]
  // and gather it to processor 0.
  for (auto bf : index_range(_local_eim_basis_functions))
    {
      // Clear any data from previous bf
      for (auto var : index_range(gathered_qp_data))
        gathered_qp_data[var].clear();

      for (const auto & pr : _local_eim_basis_functions[bf])
        {
          // array[n_vars][n_qp] per Elem
          const auto & array = pr.second;
          for (auto var : index_range(array))
            {
              // Insert all qp values for this var
              gathered_qp_data[var].insert(/*insert at*/gathered_qp_data[var].end(),
                                           /*data start*/array[var].begin(),
                                           /*data end*/array[var].end());
            }
        }

      // Reference to the data map for the current basis function.
      auto & bf_map = _local_eim_basis_functions[bf];

      for (auto var : index_range(gathered_qp_data))
        {
          // For each var, gather gathered_qp_data[var] onto processor
          // 0. There apparently is not a gather overload for
          // vector-of-vectors...
          this->comm().gather(/*root_id=*/0, gathered_qp_data[var]);

          // On processor 0, iterate over the gathered_qp_data[var]
          // vector we just gathered, filling in the "small" vectors
          // for each Elem. Note: here we ignore the fact that we
          // already have the data on processor 0 and just overwrite
          // it, this makes the indexing logic a bit simpler.
          if (this->processor_id() == 0)
            {
              auto cursor = gathered_qp_data[var].begin();
              for (auto i : index_range(elem_ids))
                {
                  auto elem_id = elem_ids[i];
                  auto n_qp_this_elem = n_qp_per_elem[i];

                  // Get reference to the [n_vars][n_qp] array for
                  // this Elem. We assign() into the vector of
                  // quadrature point values, which allocates space if
                  // it doesn't already exist.
                  auto & array = bf_map[elem_id];

                  // Possibly allocate space if this is data for a new
                  // element we haven't seen before.
                  if (array.empty())
                    array.resize(n_vars);

                  array[var].assign(cursor, cursor + n_qp_this_elem);
                  std::advance(cursor, n_qp_this_elem);
                }
            }
        }
    } // end loop over basis functions
}



void RBEIMEvaluation::distribute_bfs(const System & sys)
{
  // So we can avoid calling these many times below
  auto n_procs = sys.comm().size();
  auto rank = sys.comm().rank();

  // In serial there's nothing to distribute
  if (n_procs == 1)
    return;

  // Broadcast the number of basis functions from proc 0. After
  // distributing, all procs should have the same number of basis
  // functions.
  auto n_bf = _local_eim_basis_functions.size();
  sys.comm().broadcast(n_bf);

  // Allocate enough space to store n_bf basis functions on non-zero ranks
  if (rank != 0)
    _local_eim_basis_functions.resize(n_bf);

  // Broadcast the number of variables from proc 0. After
  // distributing, all procs should have the same number of variables.
  auto n_vars = _local_eim_basis_functions[0].begin()->second.size();
  sys.comm().broadcast(n_vars);

  // Construct lists of elem ids owned by different processors
  const MeshBase & mesh = sys.get_mesh();

  std::vector<dof_id_type> gathered_local_elem_ids;
  gathered_local_elem_ids.reserve(mesh.n_elem());
  for (const auto & elem : mesh.active_local_element_ptr_range())
    gathered_local_elem_ids.push_back(elem->id());

  // I _think_ the local elem ids are likely to already be sorted in
  // ascending order, since that is how they are stored on the Mesh,
  // but we can always just guarantee this to be on the safe side as
  // well.
  std::sort(gathered_local_elem_ids.begin(), gathered_local_elem_ids.end());

  // Gather the number of local elems from all procs to proc 0
  auto n_local_elems = gathered_local_elem_ids.size();
  std::vector<std::size_t> gathered_n_local_elems = {n_local_elems};
  sys.comm().gather(/*root_id=*/0, gathered_n_local_elems);

  // Gather the elem ids owned by each processor onto processor 0.
  sys.comm().gather(/*root_id=*/0, gathered_local_elem_ids);

  // Construct vectors of "start" and "one-past-the-end" indices into
  // the gathered_local_elem_ids vector for each proc. Only valid on
  // processor 0.
  std::vector<std::size_t> start_elem_ids_index, end_elem_ids_index;

  if (rank == 0)
    {
      start_elem_ids_index.resize(n_procs);
      start_elem_ids_index[0] = 0;
      for (processor_id_type p=1; p<n_procs; ++p)
        start_elem_ids_index[p] = start_elem_ids_index[p-1] + gathered_n_local_elems[p-1];

      end_elem_ids_index.resize(n_procs);
      end_elem_ids_index[n_procs - 1] = gathered_local_elem_ids.size();
      for (processor_id_type p=0; p<n_procs - 1; ++p)
        end_elem_ids_index[p] = start_elem_ids_index[p+1];
    }

  // On processor 0, using basis function 0 and variable 0, prepare a
  // vector with the number of qps per Elem.  Then scatter this vector
  // out to the processors that require it. The order of this vector
  // matches the gathered_local_elem_ids ordering. The counts will be
  // gathered_n_local_elems, since there will be one qp count per Elem.
  std::vector<unsigned int> n_qp_per_elem_data;

  // On rank 0, the "counts" vector holds the number of floating point values that
  // are to be scattered to each proc. It is only required on proc 0.
  std::vector<int> counts;

  if (rank == 0)
    {
      n_qp_per_elem_data.reserve(gathered_local_elem_ids.size());
      counts.resize(n_procs);

      auto & bf_map = _local_eim_basis_functions[0];

      for (processor_id_type p=0; p<n_procs; ++p)
        {
          for (auto e : make_range(start_elem_ids_index[p], end_elem_ids_index[p]))
            {
              auto elem_id = gathered_local_elem_ids[e];

              // Get reference to array[n_vars][n_qp] for current Elem.
              // Throws an error if the required elem_id is not found.
              const auto & array = libmesh_map_find(bf_map, elem_id);

              auto n_qps = array[0].size();

              // We use var==0 to set the number of qps for all vars
              n_qp_per_elem_data.push_back(n_qps);

              // Accumulate the count for this proc
              counts[p] += n_qps;
            } // end for (e)
        } // end for proc_id
    } // if (rank == 0)

  // Now scatter the n_qp_per_elem_data to all procs (must call the
  // scatter on all procs, it is a collective).
  {
    std::vector<unsigned int> recv;
    std::vector<int> tmp(gathered_n_local_elems.begin(), gathered_n_local_elems.end());
    sys.comm().scatter(n_qp_per_elem_data, tmp, recv, /*root_id=*/0);

    // Now swap n_qp_per_elem_data and recv. All processors now have a
    // vector of length n_local_elems containing the number of
    // quadarature points per Elem.
    n_qp_per_elem_data.swap(recv);
  }

  // For each basis function and each variable, build a vector
  // of qp data in the Elem ordering given by the
  // gathered_local_elem_ids, then call
  //
  // sys.comm().scatter(data, counts, recv, /*root_id=*/0);
  std::vector<std::vector<Number>> qp_data(n_vars);
  if (rank == 0)
    {
      // The total amount of qp data is given by summing the entries
      // of the "counts" vector.
      auto n_qp_data =
        std::accumulate(counts.begin(), counts.end(), 0u);

      // On processor 0, reserve enough space to hold all the qp
      // data for a single basis function for each var.
      for (auto var : index_range(qp_data))
        qp_data[var].reserve(n_qp_data);
    }

  // The recv_qp_data vector will be used on the receiving end of all
  // the scatters below.
  std::vector<Number> recv_qp_data;

  // Loop from 0..n_bf on _all_ procs, since the scatters inside this
  // loop are collective.
  for (auto bf : make_range(n_bf))
    {
      // Prepare data for scattering (only on proc 0)
      if (rank == 0)
        {
          // Reference to the data map for the current basis function.
          auto & bf_map = _local_eim_basis_functions[bf];

          // Clear any data from previous bf
          for (auto var : index_range(qp_data))
            qp_data[var].clear();

          for (processor_id_type p=0; p<n_procs; ++p)
            {
              for (auto e : make_range(start_elem_ids_index[p], end_elem_ids_index[p]))
                {
                  auto elem_id = gathered_local_elem_ids[e];

                  // Get reference to array[n_vars][n_qp] for current Elem.
                  // Throws an error if the required elem_id is not found.
                  const auto & array = libmesh_map_find(bf_map, elem_id);

                  for (auto var : index_range(array))
                    {
                      // Insert all qp values for this var
                      qp_data[var].insert(/*insert at*/qp_data[var].end(),
                                          /*data start*/array[var].begin(),
                                          /*data end*/array[var].end());
                    } // end for (var)
                } // end for (e)
            } // end for proc_id
        } // end if rank==0

      // Perform the scatters (all procs)
      for (auto var : make_range(n_vars))
        {
          // Do the scatter for the current var
          sys.comm().scatter(qp_data[var], counts, recv_qp_data, /*root_id=*/0);

          if (rank != 0)
            {
              // Store the scattered data we received in _local_eim_basis_functions[bf]
              auto & bf_map = _local_eim_basis_functions[bf];
              auto cursor = recv_qp_data.begin();

              for (auto i : index_range(gathered_local_elem_ids))
                {
                  auto elem_id = gathered_local_elem_ids[i];
                  auto n_qp_this_elem = n_qp_per_elem_data[i];
                  auto & array = bf_map[elem_id];

                  // Create space to store the data if it doesn't already exist.
                  if (array.empty())
                    array.resize(n_vars);

                  array[var].assign(cursor, cursor + n_qp_this_elem);
                  std::advance(cursor, n_qp_this_elem);
                }
            } // if (rank != 0)
        } // end for (var)
    } // end for (bf)

  // Now that the scattering is done, delete non-local Elem
  // information from processor 0's _local_eim_basis_functions data
  // structure.
  if (rank == 0)
    {
      for (processor_id_type p=1; p<n_procs; ++p)
        {
          for (auto e : make_range(start_elem_ids_index[p], end_elem_ids_index[p]))
            {
              auto elem_id = gathered_local_elem_ids[e];

              // Delete this Elem's information from every basis function.
              for (auto & bf_map : _local_eim_basis_functions)
                bf_map.erase(elem_id);
            } // end for (e)
        } // end for proc_id
    } // if (rank == 0)
}

} // namespace libMesh
