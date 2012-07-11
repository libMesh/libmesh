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
#include "rb_eim_evaluation.h"
#include "rb_eim_theta.h"
#include "rb_parametrized_function.h"

// libMesh includes
#include "o_string_stream.h"
#include "xdr_cxx.h"
#include "libmesh_logging.h"

namespace libMesh
{

RBEIMEvaluation::RBEIMEvaluation()
  :
  _previous_N(0),
  _previous_error_bound(-1)
{
  // Indicate that we need to compute the RB
  // inner product matrix in this case
  compute_RB_inner_product = true;

  // initialize to the empty RBThetaExpansion object
  set_rb_theta_expansion(_empty_rb_theta_expansion);
}

RBEIMEvaluation::~RBEIMEvaluation()
{
  this->clear();
}

void RBEIMEvaluation::clear()
{
  Parent::clear();

  interpolation_points.clear();
  interpolation_points_var.clear();

  // Delete any RBTheta objects that were created
  for(unsigned int i=0; i<_rb_eim_theta_objects.size(); i++)
  {
    delete _rb_eim_theta_objects[i];
  }
  _rb_eim_theta_objects.clear();
}

void RBEIMEvaluation::resize_data_structures(const unsigned int Nmax,
                                             bool resize_error_bound_data)
{
  Parent::resize_data_structures(Nmax, resize_error_bound_data);

  // Resize the data structures relevant to the EIM system
  interpolation_points.clear();
  interpolation_points_var.clear();
  interpolation_matrix.resize(Nmax,Nmax);

  // Resize the "extra" row due to the "extra Greedy step"
  extra_interpolation_matrix_row.resize(Nmax);
}

void RBEIMEvaluation::attach_paramerized_function(RBParametrizedFunction* pf)
{
  _parametrized_functions.push_back(pf);
}

unsigned int RBEIMEvaluation::get_n_parametrized_functions() const
{
  return _parametrized_functions.size();
}

Number RBEIMEvaluation::evaluate_parametrized_function(unsigned int var_index, const Point& p)
{
  if(var_index >= get_n_parametrized_functions())
  {
    libMesh::err << "Error: We must have var_index < get_n_parametrized_functions() in evaluate_parametrized_function."
                 << std::endl;
    libmesh_error();
  }

  return _parametrized_functions[var_index]->evaluate(get_parameters(), p);
}

Real RBEIMEvaluation::rb_solve(unsigned int N)
{
  // Short-circuit if we are using the same parameters and value of N
  if( (_previous_parameters == get_parameters()) && 
      (_previous_N == N) )
  {
    return _previous_error_bound;
  }
  
  // Otherwise, update _previous parameters, _previous_N
  _previous_parameters = get_parameters();
  _previous_N = N;
  
  START_LOG("rb_solve()", "RBEIMEvaluation");

  if(N > get_n_basis_functions())
  {
    libMesh::err << "ERROR: N cannot be larger than the number "
                 << "of basis functions in rb_solve" << std::endl;
    libmesh_error();
  }
  if(N==0)
  {
    libMesh::err << "ERROR: N must be greater than 0 in rb_solve" << std::endl;
    libmesh_error();
  }

  // Get the rhs by sampling parametrized_function
  // at the first N interpolation_points
  DenseVector<Number> EIM_rhs(N);
  for(unsigned int i=0; i<N; i++)
  {
    EIM_rhs(i) = evaluate_parametrized_function(interpolation_points_var[i], interpolation_points[i]);
  }



  DenseMatrix<Number> interpolation_matrix_N;
  interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  interpolation_matrix_N.lu_solve(EIM_rhs, RB_solution);

  // Evaluate an a posteriori error bound
  if(evaluate_RB_error_bound)
  {
    // Compute the a posteriori error bound
    // First, sample the parametrized function at x_{N+1}
    Number g_at_next_x;
    if(N == get_n_basis_functions())
      g_at_next_x = evaluate_parametrized_function(extra_interpolation_point_var, extra_interpolation_point);
    else
      g_at_next_x = evaluate_parametrized_function(interpolation_points_var[N], interpolation_points[N]);

    // Next, evaluate the EIM approximation at x_{N+1}
    Number EIM_approx_at_next_x = 0.;
    for(unsigned int j=0; j<N; j++)
    {
      if(N == get_n_basis_functions())
      {
        EIM_approx_at_next_x += RB_solution(j) * extra_interpolation_matrix_row(j);
      }
      else
      {
        EIM_approx_at_next_x += RB_solution(j) * interpolation_matrix(N,j);
      }
    }

    Real error_estimate = std::abs(g_at_next_x - EIM_approx_at_next_x);

    STOP_LOG("rb_solve()", "RBEIMEvaluation");

    _previous_error_bound = error_estimate;
    return error_estimate;
  }
  else // Don't evaluate an error bound
  {
    STOP_LOG("rb_solve()", "RBEIMEvaluation");
    _previous_error_bound = -1.;
    return -1.;
  }

}

void RBEIMEvaluation::rb_solve(DenseVector<Number>& EIM_rhs)
{
  START_LOG("rb_solve()", "RBEIMEvaluation");

  if(EIM_rhs.size() > get_n_basis_functions())
  {
    libMesh::err << "ERROR: N cannot be larger than the number "
                 << "of basis functions in rb_solve" << std::endl;
    libmesh_error();
  }
  if(EIM_rhs.size()==0)
  {
    libMesh::err << "ERROR: N must be greater than 0 in rb_solve" << std::endl;
    libmesh_error();
  }

  const unsigned int N = EIM_rhs.size();
  DenseMatrix<Number> interpolation_matrix_N;
  interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  interpolation_matrix_N.lu_solve(EIM_rhs, RB_solution);

  STOP_LOG("rb_solve()", "RBEIMEvaluation");
}

void RBEIMEvaluation::initialize_eim_theta_objects()
{
  // Initialize the rb_theta objects that access the solution from this rb_eim_evaluation
  _rb_eim_theta_objects.clear();
  for(unsigned int i=0; i<get_n_basis_functions(); i++)
  {
    _rb_eim_theta_objects.push_back(new RBEIMTheta(*this, i));
  }
}

std::vector<RBTheta*> RBEIMEvaluation::get_eim_theta_objects()
{
  return _rb_eim_theta_objects;
}

void RBEIMEvaluation::write_offline_data_to_files(const std::string& directory_name,
                                                  const bool read_binary_data)
{
  START_LOG("write_offline_data_to_files()", "RBEIMEvaluation");

  Parent::write_offline_data_to_files(directory_name);

  // Get the number of basis functions
  unsigned int n_bfs = get_n_basis_functions();

  // The writing mode: ENCODE for binary, WRITE for ASCII
  XdrMODE mode = read_binary_data ? ENCODE : WRITE;

  // The suffix to use for all the files that are written out
  const std::string suffix = read_binary_data ? ".xdr" : ".dat";

  if(libMesh::processor_id() == 0)
  {
    OStringStream file_name;
    
    // Next write out the interpolation_matrix
    file_name.str("");
    file_name << directory_name << "/interpolation_matrix" << suffix;
    Xdr interpolation_matrix_out(file_name.str(), mode);
    
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<=i; j++)
      {
        interpolation_matrix_out << interpolation_matrix(i,j);
      }
    }

    // Also, write out the "extra" row
    file_name.str("");
    file_name << directory_name << "/extra_interpolation_matrix_row" << suffix;
    Xdr extra_interpolation_matrix_row_out(file_name.str(), mode);
    
    for(unsigned int j=0; j<n_bfs; j++)
    {
      extra_interpolation_matrix_row_out << extra_interpolation_matrix_row(j);
    }
    extra_interpolation_matrix_row_out.close();

    // Next write out interpolation_points
    file_name.str("");
    file_name << directory_name << "/interpolation_points" << suffix;
    Xdr interpolation_points_out(file_name.str(), mode);
    
    for(unsigned int i=0; i<n_bfs; i++)
    {
      interpolation_points_out << interpolation_points[i](0);

      if(LIBMESH_DIM >= 2)
        interpolation_points_out << interpolation_points[i](1);

      if(LIBMESH_DIM >= 3)
        interpolation_points_out << interpolation_points[i](2);
    }
    interpolation_points_out.close();

    // Also, write out the "extra" interpolation point
    file_name.str("");
    file_name << directory_name << "/extra_interpolation_point" << suffix;
    Xdr extra_interpolation_point_out(file_name.str(), mode);
    
    extra_interpolation_point_out << extra_interpolation_point(0);

    if(LIBMESH_DIM >= 2)
      extra_interpolation_point_out << extra_interpolation_point(1);

    if(LIBMESH_DIM >= 3)
      extra_interpolation_point_out << extra_interpolation_point(2);

    extra_interpolation_point_out.close();

    // Next write out interpolation_points_var
    file_name.str("");
    file_name << directory_name << "/interpolation_points_var" << suffix;
    Xdr interpolation_points_var_out(file_name.str(), mode);

    for(unsigned int i=0; i<n_bfs; i++)
    {
      interpolation_points_var_out << interpolation_points_var[i];
    }
    interpolation_points_var_out.close();

    // Also, write out the "extra" interpolation variable
    file_name.str("");
    file_name << directory_name << "/extra_interpolation_point_var" << suffix;
    Xdr extra_interpolation_point_var_out(file_name.str(), mode);
    
    extra_interpolation_point_var_out << extra_interpolation_point_var;
    extra_interpolation_point_var_out.close();
  }

  STOP_LOG("write_offline_data_to_files()", "RBEIMEvaluation");
}

void RBEIMEvaluation::read_offline_data_from_files(const std::string& directory_name,
                                                   bool read_error_bound_data,
                                                   const bool read_binary_data)
{
  START_LOG("read_offline_data_from_files()", "RBEIMEvaluation");

  Parent::read_offline_data_from_files(directory_name, read_error_bound_data);

  // First, find out how many basis functions we had when Greedy terminated
  // This was set in RBSystem::read_offline_data_from_files
  unsigned int n_bfs = this->get_n_basis_functions();

  // The writing mode: DECODE for binary, READ for ASCII
  XdrMODE mode = read_binary_data ? DECODE : READ;

  // The suffix to use for all the files that are written out
  const std::string suffix = read_binary_data ? ".xdr" : ".dat";

  // Stream for creating file names
  OStringStream file_name;

  // Read in the interpolation matrix
  file_name.str("");
  file_name << directory_name << "/interpolation_matrix" << suffix;
  Xdr interpolation_matrix_in(file_name.str(), mode);

  for(unsigned int i=0; i<n_bfs; i++)
  {
    for(unsigned int j=0; j<=i; j++)
    {
      Number value;
      interpolation_matrix_in >> value;
      interpolation_matrix(i,j) = value;
    }
  }
  interpolation_matrix_in.close();

  // Also, read in the "extra" row
  file_name.str("");
  file_name << directory_name << "/extra_interpolation_matrix_row" << suffix;
  Xdr extra_interpolation_matrix_row_in(file_name.str(), mode);
  
  for(unsigned int j=0; j<n_bfs; j++)
  {
    Number value;
    extra_interpolation_matrix_row_in >> value;
    extra_interpolation_matrix_row(j) = value;
  }
  extra_interpolation_matrix_row_in.close();

  // Next read in interpolation_points
  file_name.str("");
  file_name << directory_name << "/interpolation_points" << suffix;
  Xdr interpolation_points_in(file_name.str(), mode);
  
  for(unsigned int i=0; i<n_bfs; i++)
  {
    Real x_val, y_val, z_val = 0.;
    interpolation_points_in >> x_val;
    
    if(LIBMESH_DIM >= 2)
      interpolation_points_in >> y_val;

    if(LIBMESH_DIM >= 3)
      interpolation_points_in >> z_val;

    Point p(x_val, y_val, z_val);
    interpolation_points.push_back(p);
  }
  interpolation_points_in.close();

  // Also, read in the extra interpolation point
  file_name.str("");
  file_name << directory_name << "/extra_interpolation_point" << suffix;
  Xdr extra_interpolation_point_in(file_name.str(), mode);
  
  for(unsigned int i=0; i<n_bfs; i++)
  {
    Real x_val, y_val, z_val = 0.;
    extra_interpolation_point_in >> x_val;

    if(LIBMESH_DIM >= 2)
      extra_interpolation_point_in >> y_val;

    if(LIBMESH_DIM >= 3)
      extra_interpolation_point_in >> z_val;

    Point p(x_val, y_val, z_val);
    extra_interpolation_point = p;
  }
  extra_interpolation_point_in.close();


  // Next read in interpolation_points_var
  file_name.str("");
  file_name << directory_name << "/interpolation_points_var" << suffix;
  Xdr interpolation_points_var_in(file_name.str(), mode);
  
  for(unsigned int i=0; i<=n_bfs; i++)
  {
    unsigned int var;
    interpolation_points_var_in >> var;
    interpolation_points_var.push_back(var);
  }
  interpolation_points_var_in.close();

  // Also, read in extra_interpolation_point_var
  file_name.str("");
  file_name << directory_name << "/extra_interpolation_point_var" << suffix;
  Xdr extra_interpolation_point_var_in(file_name.str(), mode);
  
  for(unsigned int i=0; i<=n_bfs; i++)
  {
    unsigned int var;
    extra_interpolation_point_var_in >> var;
    extra_interpolation_point_var = var;
  }
  extra_interpolation_point_var_in.close();

  STOP_LOG("read_offline_data_from_files()", "RBEIMEvaluation");
}

}
