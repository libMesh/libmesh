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

namespace libMesh
{

RBEIMEvaluation::RBEIMEvaluation(const libMesh::Parallel::Communicator & comm_in)
  :
  RBEvaluation(comm_in),
  _previous_N(0),
  _previous_error_bound(-1),
  _interpolation_points_mesh(comm_in)
{
  // Indicate that we need to compute the RB
  // inner product matrix in this case
  compute_RB_inner_product = true;

  // initialize to the empty RBThetaExpansion object
  set_rb_theta_expansion(_empty_rb_theta_expansion);

  // Let's not renumber the _interpolation_points_mesh
  _interpolation_points_mesh.allow_renumbering(false);
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
  interpolation_points_elem.clear();
  _interpolation_points_mesh.clear();

  // Delete any RBTheta objects that were created
  for (std::size_t i=0; i<_rb_eim_theta_objects.size(); i++)
    delete _rb_eim_theta_objects[i];
  _rb_eim_theta_objects.clear();
}

void RBEIMEvaluation::resize_data_structures(const unsigned int Nmax,
                                             bool resize_error_bound_data)
{
  Parent::resize_data_structures(Nmax, resize_error_bound_data);

  // Resize the data structures relevant to the EIM system
  interpolation_points.clear();
  interpolation_points_var.clear();
  interpolation_points_elem.clear();
  interpolation_matrix.resize(Nmax,Nmax);
}

void RBEIMEvaluation::attach_parametrized_function(RBParametrizedFunction * pf)
{
  _parametrized_functions.push_back(pf);
}

unsigned int RBEIMEvaluation::get_n_parametrized_functions() const
{
  return cast_int<unsigned int>
    (_parametrized_functions.size());
}

ReplicatedMesh & RBEIMEvaluation::get_interpolation_points_mesh()
{
  return _interpolation_points_mesh;
}

Number RBEIMEvaluation::evaluate_parametrized_function(unsigned int var_index,
                                                       const Point & p,
                                                       const Elem & elem)
{
  if(var_index >= get_n_parametrized_functions())
    libmesh_error_msg("Error: We must have var_index < get_n_parametrized_functions() in evaluate_parametrized_function.");

  return _parametrized_functions[var_index]->evaluate(get_parameters(), p, elem);
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

  LOG_SCOPE("rb_solve()", "RBEIMEvaluation");

  if(N > get_n_basis_functions())
    libmesh_error_msg("ERROR: N cannot be larger than the number of basis functions in rb_solve");

  if(N==0)
    libmesh_error_msg("ERROR: N must be greater than 0 in rb_solve");

  // Get the rhs by sampling parametrized_function
  // at the first N interpolation_points
  DenseVector<Number> EIM_rhs(N);
  for(unsigned int i=0; i<N; i++)
    {
      EIM_rhs(i) = evaluate_parametrized_function(interpolation_points_var[i],
                                                  interpolation_points[i],
                                                  *interpolation_points_elem[i]);
    }



  DenseMatrix<Number> interpolation_matrix_N;
  interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  interpolation_matrix_N.lu_solve(EIM_rhs, RB_solution);

  // Optionally evaluate an a posteriori error bound. The EIM error estimate
  // recommended in the literature is based on using "next" EIM point, so
  // we skip this if N == get_n_basis_functions()
  if(evaluate_RB_error_bound && (N != get_n_basis_functions()))
    {
      // Compute the a posteriori error bound
      // First, sample the parametrized function at x_{N+1}
      Number g_at_next_x = evaluate_parametrized_function(interpolation_points_var[N],
                                                          interpolation_points[N],
                                                          *interpolation_points_elem[N]);

      // Next, evaluate the EIM approximation at x_{N+1}
      Number EIM_approx_at_next_x = 0.;
      for(unsigned int j=0; j<N; j++)
        {
          EIM_approx_at_next_x += RB_solution(j) * interpolation_matrix(N,j);
        }

      Real error_estimate = std::abs(g_at_next_x - EIM_approx_at_next_x);

      _previous_error_bound = error_estimate;
      return error_estimate;
    }
  else // Don't evaluate an error bound
    {
      _previous_error_bound = -1.;
      return -1.;
    }

}

void RBEIMEvaluation::rb_solve(DenseVector<Number> & EIM_rhs)
{
  LOG_SCOPE("rb_solve()", "RBEIMEvaluation");

  if(EIM_rhs.size() > get_n_basis_functions())
    libmesh_error_msg("ERROR: N cannot be larger than the number of basis functions in rb_solve");

  if(EIM_rhs.size()==0)
    libmesh_error_msg("ERROR: N must be greater than 0 in rb_solve");

  const unsigned int N = EIM_rhs.size();
  DenseMatrix<Number> interpolation_matrix_N;
  interpolation_matrix.get_principal_submatrix(N, interpolation_matrix_N);

  interpolation_matrix_N.lu_solve(EIM_rhs, RB_solution);
}

Real RBEIMEvaluation::get_error_bound_normalization()
{
  // Just set the normalization factor to 1 in this case.
  // Users can override this method if specific behavior
  // is required.

  return 1.;
}

void RBEIMEvaluation::initialize_eim_theta_objects()
{
  // Initialize the rb_theta objects that access the solution from this rb_eim_evaluation
  _rb_eim_theta_objects.clear();
  for(unsigned int i=0; i<get_n_basis_functions(); i++)
    {
      _rb_eim_theta_objects.push_back( build_eim_theta(i).release() );
    }
}

std::vector<RBTheta *> RBEIMEvaluation::get_eim_theta_objects()
{
  return _rb_eim_theta_objects;
}

UniquePtr<RBTheta> RBEIMEvaluation::build_eim_theta(unsigned int index)
{
  return UniquePtr<RBTheta>( new RBEIMTheta(*this, index) );
}

void RBEIMEvaluation::legacy_write_offline_data_to_files(const std::string & directory_name,
                                                         const bool read_binary_data)
{
  LOG_SCOPE("legacy_write_offline_data_to_files()", "RBEIMEvaluation");

  Parent::legacy_write_offline_data_to_files(directory_name);

  // Get the number of basis functions
  unsigned int n_bfs = get_n_basis_functions();

  // The writing mode: ENCODE for binary, WRITE for ASCII
  XdrMODE mode = read_binary_data ? ENCODE : WRITE;

  // The suffix to use for all the files that are written out
  const std::string suffix = read_binary_data ? ".xdr" : ".dat";

  if(this->processor_id() == 0)
    {
      std::ostringstream file_name;

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

      // Next write out interpolation_points_var
      file_name.str("");
      file_name << directory_name << "/interpolation_points_var" << suffix;
      Xdr interpolation_points_var_out(file_name.str(), mode);

      for(unsigned int i=0; i<n_bfs; i++)
        {
          interpolation_points_var_out << interpolation_points_var[i];
        }
      interpolation_points_var_out.close();
    }

  // Write out the elements associated with the interpolation points.
  // This uses mesh I/O, hence we have to do it on all processors.
  legacy_write_out_interpolation_points_elem(directory_name);
}

void RBEIMEvaluation::legacy_write_out_interpolation_points_elem(const std::string & directory_name)
{
  _interpolation_points_mesh.clear();

  // Maintain a set of node IDs to make sure we don't insert
  // the same node into _interpolation_points_mesh more than once
  std::set<dof_id_type> node_ids;
  std::map<dof_id_type, dof_id_type> node_id_map;

  unsigned int new_node_id = 0;
  for (std::size_t i=0; i<interpolation_points_elem.size(); i++)
    {
      Elem * old_elem = interpolation_points_elem[i];

      for(unsigned int n=0; n<old_elem->n_nodes(); n++)
        {
          Node & node_ref = old_elem->node_ref(n);
          dof_id_type old_node_id = node_ref.id();

          // Check if this node has already been added. This
          // could happen if some of the elements are neighbors.
          if( node_ids.find(old_node_id) == node_ids.end() )
            {
              node_ids.insert(old_node_id);
              _interpolation_points_mesh.add_point(node_ref, new_node_id, /* proc_id */ 0);

              node_id_map[old_node_id] = new_node_id;

              new_node_id++;
            }
        }
    }

  // Maintain a map of elem IDs to make sure we don't insert
  // the same elem into _interpolation_points_mesh more than once
  std::map<dof_id_type,dof_id_type> elem_id_map;
  std::vector<dof_id_type> interpolation_elem_ids(interpolation_points_elem.size());
  dof_id_type new_elem_id = 0;
  for (std::size_t i=0; i<interpolation_elem_ids.size(); i++)
    {
      Elem * old_elem = interpolation_points_elem[i];

      dof_id_type old_elem_id = old_elem->id();

      // Only insert the element into the mesh if it hasn't already been inserted
      std::map<dof_id_type,dof_id_type>::iterator id_it = elem_id_map.find(old_elem_id);
      if(id_it == elem_id_map.end())
        {
          Elem * new_elem = Elem::build(old_elem->type(), /*parent*/ libmesh_nullptr).release();
          new_elem->subdomain_id() = old_elem->subdomain_id();

          // Assign all the nodes
          for(unsigned int n=0; n<new_elem->n_nodes(); n++)
            {
              dof_id_type old_node_id = old_elem->node_id(n);
              new_elem->set_node(n) =
                _interpolation_points_mesh.node_ptr( node_id_map[old_node_id] );
            }

          // Just set all proc_ids to 0
          new_elem->processor_id() = 0;

          // Add the element to the mesh
          _interpolation_points_mesh.add_elem(new_elem);

          // Set the id of new_elem appropriately
          new_elem->set_id(new_elem_id);
          interpolation_elem_ids[i] = new_elem->id();
          elem_id_map[old_elem_id] = new_elem->id();

          new_elem_id++;
        }
      else
        {
          interpolation_elem_ids[i] = id_it->second;
        }

    }

  libmesh_assert(new_elem_id == _interpolation_points_mesh.n_elem());

  _interpolation_points_mesh.write(directory_name + "/interpolation_points_mesh.xda");

  // Also, write out the vector that tells us which element each entry
  // of interpolation_points_elem corresponds to. This allows us to handle
  // the case in which elements are repeated in interpolation_points_elem.
  if(processor_id() == 0)
    {
      // These are just integers, so no need for a binary format here
      std::ofstream interpolation_elem_ids_out
        ((directory_name + "/interpolation_elem_ids.dat").c_str(), std::ofstream::out);

      for (std::size_t i=0; i<interpolation_elem_ids.size(); i++)
        interpolation_elem_ids_out << interpolation_elem_ids[i] << std::endl;

      interpolation_elem_ids_out.close();
    }
}

void RBEIMEvaluation::legacy_read_offline_data_from_files(const std::string & directory_name,
                                                          bool read_error_bound_data,
                                                          const bool read_binary_data)
{
  LOG_SCOPE("legacy_read_offline_data_from_files()", "RBEIMEvaluation");

  Parent::legacy_read_offline_data_from_files(directory_name, read_error_bound_data);

  // First, find out how many basis functions we had when Greedy terminated
  // This was set in RBSystem::read_offline_data_from_files
  unsigned int n_bfs = this->get_n_basis_functions();

  // The writing mode: DECODE for binary, READ for ASCII
  XdrMODE mode = read_binary_data ? DECODE : READ;

  // The suffix to use for all the files that are written out
  const std::string suffix = read_binary_data ? ".xdr" : ".dat";

  // Stream for creating file names
  std::ostringstream file_name;

  // Read in the interpolation matrix
  file_name.str("");
  file_name << directory_name << "/interpolation_matrix" << suffix;
  assert_file_exists(file_name.str());

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

  // Next read in interpolation_points
  file_name.str("");
  file_name << directory_name << "/interpolation_points" << suffix;
  assert_file_exists(file_name.str());

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

  // Next read in interpolation_points_var
  file_name.str("");
  file_name << directory_name << "/interpolation_points_var" << suffix;
  assert_file_exists(file_name.str());

  Xdr interpolation_points_var_in(file_name.str(), mode);

  for(unsigned int i=0; i<n_bfs; i++)
    {
      unsigned int var;
      interpolation_points_var_in >> var;
      interpolation_points_var.push_back(var);
    }
  interpolation_points_var_in.close();

  // Read in the elements corresponding to the interpolation points
  legacy_read_in_interpolation_points_elem(directory_name);
}

void RBEIMEvaluation::legacy_read_in_interpolation_points_elem(const std::string & directory_name)
{
  _interpolation_points_mesh.read(directory_name + "/interpolation_points_mesh.xda");

  // We have an element for each EIM basis function
  unsigned int n_bfs = this->get_n_basis_functions();

  std::vector<dof_id_type> interpolation_elem_ids;
  {
    // These are just integers, so no need for a binary format here
    std::ifstream interpolation_elem_ids_in
      ((directory_name + "/interpolation_elem_ids.dat").c_str(), std::ifstream::in);

    if (!interpolation_elem_ids_in)
      libmesh_error_msg("RB data missing: " + directory_name + "/interpolation_elem_ids.dat");

    for(unsigned int i=0; i<n_bfs; i++)
      {
        dof_id_type elem_id;
        interpolation_elem_ids_in >> elem_id;
        interpolation_elem_ids.push_back(elem_id);
      }
    interpolation_elem_ids_in.close();
  }

  interpolation_points_elem.resize(n_bfs);
  for(unsigned int i=0; i<n_bfs; i++)
    {
      interpolation_points_elem[i] =
        _interpolation_points_mesh.elem_ptr(interpolation_elem_ids[i]);
    }
}

}
