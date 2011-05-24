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

#include "rb_eim_evaluation.h"
#include "rb_eim_system.h"
#include "libmesh_logging.h"
#include "rb_eim_theta.h"

namespace libMesh
{

RBEIMEvaluation::RBEIMEvaluation(RBEIMSystem& rb_eim_sys_in)
  :
  rb_eim_sys(rb_eim_sys_in)
{}


void RBEIMEvaluation::clear()
{
  Parent::clear();
  
  interpolation_points.clear();
  interpolation_points_var.clear();
}

void RBEIMEvaluation::initialize(const unsigned int Nmax)
{
  Parent::initialize(Nmax);

  // Resize the data structures relevant to the EIM system
  interpolation_points.clear();
  interpolation_points_var.clear();
  interpolation_matrix.resize(Nmax,Nmax);
  
  // Resize the "extra" row due to the "extra Greedy step"
  extra_interpolation_matrix_row.resize(Nmax);
}

Real RBEIMEvaluation::RB_solve(unsigned int N)
{
  START_LOG("RB_solve()", "RBEIMEvaluation");

  if(N > get_n_basis_functions())
  {
    libMesh::err << "ERROR: N cannot be larger than the number "
                 << "of basis functions in RB_solve" << std::endl;
    libmesh_error();
  }
  if(N==0)
  {
    libMesh::err << "ERROR: N must be greater than 0 in RB_solve" << std::endl;
    libmesh_error();
  }
  
  // This function uses rb_eim_sys, so set rb_eim_sys's parameter
  rb_eim_sys.set_current_parameters( get_current_parameters() );

  // Get the rhs by sampling parametrized_function
  // at the first N interpolation_points
  DenseVector<Number> EIM_rhs(N);
  for(unsigned int i=0; i<N; i++)
  {
    EIM_rhs(i) = rb_eim_sys.evaluate_parametrized_function(interpolation_points_var[i], interpolation_points[i]);
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
      g_at_next_x = rb_eim_sys.evaluate_parametrized_function(extra_interpolation_point_var, extra_interpolation_point);
    else
      g_at_next_x = rb_eim_sys.evaluate_parametrized_function(interpolation_points_var[N], interpolation_points[N]);

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

    STOP_LOG("RB_solve()", "RBEIMEvaluation");
  
    return error_estimate;
  }
  else // Don't evaluate an error bound
  {
    STOP_LOG("RB_solve()", "RBEIMEvaluation");
    return -1.;
  }

}

void RBEIMEvaluation::RB_solve(DenseVector<Number>& EIM_rhs)
{
  START_LOG("RB_solve()", "RBEIMEvaluation");
  
  if(EIM_rhs.size() > get_n_basis_functions())
  {
    libMesh::err << "ERROR: N cannot be larger than the number "
                 << "of basis functions in RB_solve" << std::endl;
    libmesh_error();
  }
  if(EIM_rhs.size()==0)
  {
    libMesh::err << "ERROR: N must be greater than 0 in RB_solve" << std::endl;
    libmesh_error();
  }
  
  const unsigned int N = EIM_rhs.size();
  DenseMatrix<Number> interpolation_matrix_N;
  interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);
  
  interpolation_matrix_N.lu_solve(EIM_rhs, RB_solution);
  
  STOP_LOG("RB_solve()", "RBEIMEvaluation");
}

AutoPtr<RBTheta> RBEIMEvaluation::build_rb_eim_theta(unsigned int index)
{
  return AutoPtr<RBTheta>(new RBEIMTheta(*this, index));
}

void RBEIMEvaluation::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "RBEIMEvaluation");

  Parent::write_offline_data_to_files(directory_name);

  const unsigned int n_bfs = get_n_basis_functions();

  const unsigned int precision_level = 14;

  if(libMesh::processor_id() == 0)
  {
    // Next write out the interpolation_matrix
    std::ofstream interpolation_matrix_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/interpolation_matrix.dat";
      interpolation_matrix_out.open(file_name.str().c_str());
    }
    if ( !interpolation_matrix_out.good() )
    {
      libMesh::err << "Error opening interpolation_matrix.dat" << std::endl;
      libmesh_error();
    }
    interpolation_matrix_out.precision(precision_level);
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<=i; j++)
      {
        interpolation_matrix_out << std::scientific
          << interpolation_matrix(i,j) << " ";
      }
    }
    
    // Also, write out the "extra" row
    std::ofstream extra_interpolation_matrix_row_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/extra_interpolation_matrix_row.dat";
      extra_interpolation_matrix_row_out.open(file_name.str().c_str());
    }
    if ( !extra_interpolation_matrix_row_out.good() )
    {
      libMesh::err << "Error opening extra_interpolation_matrix_row.dat" << std::endl;
      libmesh_error();
    }
    extra_interpolation_matrix_row_out.precision(precision_level);
    for(unsigned int j=0; j<n_bfs; j++)
      extra_interpolation_matrix_row_out << std::scientific
          << extra_interpolation_matrix_row(j) << " ";
    extra_interpolation_matrix_row_out.close();
    
    // Next write out interpolation_points
    std::ofstream interpolation_points_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/interpolation_points.dat";
      interpolation_points_out.open(file_name.str().c_str());
    }
    if ( !interpolation_points_out.good() )
    {
      libMesh::err << "Error opening interpolation_points.dat" << std::endl;
      libmesh_error();
    }
    interpolation_points_out.precision(precision_level);
    for(unsigned int i=0; i<n_bfs; i++)
      interpolation_points_out << std::scientific
          << interpolation_points[i](0) << " "
          << interpolation_points[i](1) << " "
          << interpolation_points[i](2) << " ";

    // Also, write out the "extra" interpolation point
    std::ofstream extra_interpolation_point_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/extra_interpolation_point.dat";
      extra_interpolation_point_out.open(file_name.str().c_str());
    }
    if ( !extra_interpolation_point_out.good() )
    {
      libMesh::err << "Error opening extra_interpolation_point.dat" << std::endl;
      libmesh_error();
    }
    extra_interpolation_point_out.precision(precision_level);
    extra_interpolation_point_out << std::scientific
          << extra_interpolation_point(0) << " "
          << extra_interpolation_point(1) << " "
          << extra_interpolation_point(2) << " ";
    extra_interpolation_point_out.close();
    
    // Next write out interpolation_points_var
    std::ofstream interpolation_points_var_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/interpolation_points_var.dat";
      interpolation_points_var_out.open(file_name.str().c_str());
    }
    if ( !interpolation_points_var_out.good() )
    {
      libMesh::err << "Error opening interpolation_points_var.dat" << std::endl;
      libmesh_error();
    }
    interpolation_points_var_out.precision(precision_level);
    for(unsigned int i=0; i<n_bfs; i++)
      interpolation_points_var_out << std::scientific
          << interpolation_points_var[i] << " ";

    // Also, write out the "extra" interpolation variable
    std::ofstream extra_interpolation_point_var_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/extra_interpolation_point_var.dat";
      extra_interpolation_point_var_out.open(file_name.str().c_str());
    }
    if ( !extra_interpolation_point_var_out.good() )
    {
      libMesh::err << "Error opening extra_interpolation_point_var.dat" << std::endl;
      libmesh_error();
    }
    extra_interpolation_point_var_out.precision(precision_level);
    extra_interpolation_point_var_out << std::scientific
          << extra_interpolation_point_var << " ";
    extra_interpolation_point_var_out.close();
  }

  STOP_LOG("write_offline_data_to_files()", "RBEIMEvaluation");
}

void RBEIMEvaluation::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "RBEIMEvaluation");

  Parent::read_offline_data_from_files(directory_name);
  
  // First, find out how many basis functions we had when Greedy terminated
  // This was set in RBSystem::read_offline_data_from_files
  unsigned int n_bfs = this->get_n_basis_functions();
  
  // Read in the interpolation matrix
  std::ifstream interpolation_matrix_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/interpolation_matrix.dat";
    interpolation_matrix_in.open(file_name.str().c_str());
  }
  if ( !interpolation_matrix_in.good() )
  {
    libMesh::err << "Error opening interpolation_matrix.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<n_bfs; i++)
  {
    for(unsigned int j=0; j<=i; j++)
    {
      Number value;
      interpolation_matrix_in >> value;
      interpolation_matrix(i,j) = value;
    }
  }

  // Also, read in the "extra" row
  std::ifstream extra_interpolation_matrix_row_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/extra_interpolation_matrix_row.dat";
    extra_interpolation_matrix_row_in.open(file_name.str().c_str());
  }
  if ( !extra_interpolation_matrix_row_in.good() )
  {
    libMesh::err << "Error opening extra_interpolation_matrix_row.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int j=0; j<n_bfs; j++)
  {
    Number value;
    extra_interpolation_matrix_row_in >> value;
    extra_interpolation_matrix_row(j) = value;
  }
  extra_interpolation_matrix_row_in.close();

  // Next read in interpolation_points
  std::ifstream interpolation_points_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/interpolation_points.dat";
    interpolation_points_in.open(file_name.str().c_str());
  }
  if ( !interpolation_points_in.good() )
  {
    libMesh::err << "Error opening interpolation_points.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<n_bfs; i++)
  {
    Real x_val, y_val, z_val;
    interpolation_points_in >> x_val;
    interpolation_points_in >> y_val;
    interpolation_points_in >> z_val;
    Point p(x_val, y_val, z_val);
    interpolation_points.push_back(p);
  }
  interpolation_points_in.close();
  
  // Also, read in the extra interpolation point
  std::ifstream extra_interpolation_point_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/extra_interpolation_point.dat";
    extra_interpolation_point_in.open(file_name.str().c_str());
  }
  if ( !extra_interpolation_point_in.good() )
  {
    libMesh::err << "Error opening extra_interpolation_point.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<n_bfs; i++)
  {
    Real x_val, y_val, z_val;
    extra_interpolation_point_in >> x_val;
    extra_interpolation_point_in >> y_val;
    extra_interpolation_point_in >> z_val;
    Point p(x_val, y_val, z_val);
    extra_interpolation_point = p;
  }
  extra_interpolation_point_in.close();
  

  // Next read in interpolation_points_var
  std::ifstream interpolation_points_var_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/interpolation_points_var.dat";
    interpolation_points_var_in.open(file_name.str().c_str());
  }
  if ( !interpolation_points_var_in.good() )
  {
    libMesh::err << "Error opening interpolation_points_var.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<=n_bfs; i++)
  {
    unsigned int var;
    interpolation_points_var_in >> var;
    interpolation_points_var.push_back(var);
  }
  interpolation_points_var_in.close();
  
  // Also, read in extra_interpolation_point_var
  std::ifstream extra_interpolation_point_var_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/extra_interpolation_point_var.dat";
    extra_interpolation_point_var_in.open(file_name.str().c_str());
  }
  if ( !extra_interpolation_point_var_in.good() )
  {
    libMesh::err << "Error opening extra_interpolation_point_var.dat" << std::endl;
    libmesh_error();
  }
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
