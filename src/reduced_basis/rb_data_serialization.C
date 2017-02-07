// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010, 2015 David J. Knezevic

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

#include "libmesh/libmesh_config.h"
#if defined(LIBMESH_HAVE_CAPNPROTO)

//libMesh includes
#include "libmesh/rb_data_serialization.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/transient_rb_theta_expansion.h"
#include "libmesh/rb_evaluation.h"
#include "libmesh/transient_rb_evaluation.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_scm_evaluation.h"
#include "libmesh/elem.h"

// Cap'n'Proto includes
#include <capnp/serialize.h>

// C++ includes
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>

namespace libMesh
{

namespace
{

/**
 * Helper function that sets either real or complex numbers, based on
 * the libMesh config options.
 */
template <typename T, typename U>
void set_scalar_in_list(T list, unsigned int i, U value)
{
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  list[i].setReal(std::real(value));
  list[i].setImag(std::imag(value));
#else
  list.set(i, value);
#endif
}

}

namespace RBDataSerialization
{

// ---- RBEvaluationSerialization (BEGIN) ----

RBEvaluationSerialization::RBEvaluationSerialization(RBEvaluation & rb_eval)
  :
  _rb_eval(rb_eval)
{
}

RBEvaluationSerialization::~RBEvaluationSerialization()
{
}

void RBEvaluationSerialization::write_to_file(const std::string & path)
{
  LOG_SCOPE("write_to_file()", "RBEvaluationSerialization");

  if(_rb_eval.comm().rank() == 0)
    {
      capnp::MallocMessageBuilder message;

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
      RBData::RBEvaluationReal::Builder rb_eval_builder =
        message.initRoot<RBData::RBEvaluationReal>();
#else
      RBData::RBEvaluationComplex::Builder rb_eval_builder =
        message.initRoot<RBData::RBEvaluationComplex>();
#endif

      add_rb_evaluation_data_to_builder(_rb_eval, rb_eval_builder);

      int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
      if (!fd)
        libmesh_error_msg("Error opening a write-only file descriptor to " + path);

      capnp::writeMessageToFd(fd, message);

      int error = close(fd);
      if (error)
        libmesh_error_msg("Error closing a write-only file descriptor to " + path);
    }
}

// ---- RBEvaluationSerialization (END) ----


// ---- TransientRBEvaluationSerialization (BEGIN) ----

TransientRBEvaluationSerialization::
TransientRBEvaluationSerialization(TransientRBEvaluation & trans_rb_eval) :
  _trans_rb_eval(trans_rb_eval)
{
}

TransientRBEvaluationSerialization::~TransientRBEvaluationSerialization()
{
}

void TransientRBEvaluationSerialization::write_to_file(const std::string & path)
{
  LOG_SCOPE("write_to_file()", "TransientRBEvaluationSerialization");

  if(_trans_rb_eval.comm().rank() == 0)
    {
      capnp::MallocMessageBuilder message;

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
      RBData::TransientRBEvaluationReal::Builder trans_rb_eval_builder =
        message.initRoot<RBData::TransientRBEvaluationReal>();
      RBData::RBEvaluationReal::Builder rb_eval_builder =
        trans_rb_eval_builder.initRbEvaluation();
#else
      RBData::TransientRBEvaluationComplex::Builder trans_rb_eval_builder =
        message.initRoot<RBData::TransientRBEvaluationComplex>();
      RBData::RBEvaluationComplex::Builder rb_eval_builder =
        trans_rb_eval_builder.initRbEvaluation();
#endif

      add_transient_rb_evaluation_data_to_builder(_trans_rb_eval,
                                                  rb_eval_builder,
                                                  trans_rb_eval_builder);

      int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
      if (!fd)
        libmesh_error_msg("Error opening a write-only file descriptor to " + path);

      capnp::writeMessageToFd(fd, message);

      int error = close(fd);
      if (error)
        libmesh_error_msg("Error closing a write-only file descriptor to " + path);
    }
}

// ---- TransientRBEvaluationSerialization (END) ----


// ---- RBEIMEvaluationSerialization (BEGIN) ----

RBEIMEvaluationSerialization::RBEIMEvaluationSerialization(RBEIMEvaluation & rb_eim_eval)
  :
  _rb_eim_eval(rb_eim_eval)
{
}

RBEIMEvaluationSerialization::~RBEIMEvaluationSerialization()
{
}

void RBEIMEvaluationSerialization::write_to_file(const std::string & path)
{
  LOG_SCOPE("write_to_file()", "RBEIMEvaluationSerialization");

  if(_rb_eim_eval.comm().rank() == 0)
    {
      capnp::MallocMessageBuilder message;

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
      RBData::RBEIMEvaluationReal::Builder rb_eim_eval_builder =
        message.initRoot<RBData::RBEIMEvaluationReal>();
      RBData::RBEvaluationReal::Builder rb_eval_builder =
        rb_eim_eval_builder.initRbEvaluation();
#else
      RBData::RBEIMEvaluationComplex::Builder rb_eim_eval_builder =
        message.initRoot<RBData::RBEIMEvaluationComplex>();
      RBData::RBEvaluationComplex::Builder rb_eval_builder =
        rb_eim_eval_builder.initRbEvaluation();
#endif

      add_rb_eim_evaluation_data_to_builder(_rb_eim_eval,
                                            rb_eval_builder,
                                            rb_eim_eval_builder);

      int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
      if (!fd)
        libmesh_error_msg("Error opening a write-only file descriptor to " + path);

      capnp::writeMessageToFd(fd, message);

      int error = close(fd);
      if (error)
        libmesh_error_msg("Error closing a write-only file descriptor to " + path);
    }
}

// ---- RBEIMEvaluationSerialization (END) ----


// ---- RBSCMEvaluationSerialization (BEGIN) ----

#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

RBSCMEvaluationSerialization::RBSCMEvaluationSerialization(RBSCMEvaluation & rb_scm_eval)
  :
  _rb_scm_eval(rb_scm_eval)
{
}

RBSCMEvaluationSerialization::~RBSCMEvaluationSerialization()
{
}

void RBSCMEvaluationSerialization::write_to_file(const std::string & path)
{
  LOG_SCOPE("write_to_file()", "RBSCMEvaluationSerialization");

  if(_rb_scm_eval.comm().rank() == 0)
    {
      capnp::MallocMessageBuilder message;

      RBData::RBSCMEvaluation::Builder rb_scm_eval_builder =
        message.initRoot<RBData::RBSCMEvaluation>();

      add_rb_scm_evaluation_data_to_builder(_rb_scm_eval, rb_scm_eval_builder);

      int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
      if (!fd)
        libmesh_error_msg("Error opening a write-only file descriptor to " + path);

      capnp::writeMessageToFd(fd, message);

      int error = close(fd);
      if (error)
        libmesh_error_msg("Error closing a write-only file descriptor to " + path);
    }
}

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

// ---- RBSCMEvaluationSerialization (END) ----


// ---- Helper functions for adding data to capnp Builders (BEGIN) ----

void add_parameter_ranges_to_builder(const RBParametrized & rb_evaluation,
                                     RBData::ParameterRanges::Builder & parameter_ranges_list,
                                     RBData::DiscreteParameterList::Builder & discrete_parameters_list)
{
  // Continuous parameters
  {
    unsigned int n_continuous_parameters = rb_evaluation.get_n_continuous_params();
    auto names = parameter_ranges_list.initNames(n_continuous_parameters);
    auto mins = parameter_ranges_list.initMinValues(n_continuous_parameters);
    auto maxs = parameter_ranges_list.initMaxValues(n_continuous_parameters);

    std::set<std::string> parameter_names = rb_evaluation.get_parameter_names();
    const RBParameters & parameters_min = rb_evaluation.get_parameters_min();
    const RBParameters & parameters_max = rb_evaluation.get_parameters_max();

    unsigned int count = 0;
    for(const auto & parameter_name : parameter_names)
      {
        if(!rb_evaluation.is_discrete_parameter(parameter_name))
          {
            names.set(count, parameter_name);
            mins.set(count, parameters_min.get_value(parameter_name));
            maxs.set(count, parameters_max.get_value(parameter_name));

            ++count;
          }
      }

    if (count != n_continuous_parameters)
      libmesh_error_msg("Mismatch in number of continuous parameters");
  }

  // Discrete parameters
  {
    unsigned int n_discrete_parameters = rb_evaluation.get_n_discrete_params();
    auto names = discrete_parameters_list.initNames(n_discrete_parameters);
    auto values_outer = discrete_parameters_list.initValues(n_discrete_parameters);

    const std::map< std::string, std::vector<Real> > & discrete_parameters =
      rb_evaluation.get_discrete_parameter_values();

    unsigned int count = 0;
    for(const auto & discrete_parameter : discrete_parameters)
      {
        names.set(count, discrete_parameter.first);

        const std::vector<Real> & values = discrete_parameter.second;
        unsigned int n_values = values.size();

        values_outer.init(count, n_values);
        auto values_inner = values_outer[count];
        for(unsigned int i=0; i<n_values; ++i)
          {
            values_inner.set(i, values[i]);
          }

        ++count;
      }

    if (count != n_discrete_parameters)
      libmesh_error_msg("Mismatch in number of discrete parameters");
  }
}

template <typename RBEvaluationBuilderNumber>
void add_rb_evaluation_data_to_builder(RBEvaluation & rb_eval,
                                       RBEvaluationBuilderNumber & rb_evaluation_builder)
{
  const RBThetaExpansion & rb_theta_expansion = rb_eval.get_rb_theta_expansion();

  unsigned int n_F_terms = rb_theta_expansion.get_n_F_terms();
  unsigned int n_A_terms = rb_theta_expansion.get_n_A_terms();

  // Number of basis functions
  unsigned int n_bfs = rb_eval.get_n_basis_functions();
  rb_evaluation_builder.setNBfs(n_bfs);

  // Fq representor inner-product data
  {
    unsigned int Q_f_hat = n_F_terms*(n_F_terms+1)/2;

    auto fq_innerprods_list = rb_evaluation_builder.initFqInnerprods(Q_f_hat);

    for(unsigned int i=0; i<Q_f_hat; i++)
      set_scalar_in_list(fq_innerprods_list,
                         i,
                         rb_eval.Fq_representor_innerprods[i]);
  }

  // FqAq representor inner-product data
  {
    auto fq_aq_innerprods_list =
      rb_evaluation_builder.initFqAqInnerprods(n_F_terms*n_A_terms*n_bfs);

    for(unsigned int q_f=0; q_f < n_F_terms; ++q_f)
      for(unsigned int q_a=0; q_a < n_A_terms; ++q_a)
        for(unsigned int i=0; i < n_bfs; ++i)
          {
            unsigned int offset = q_f*n_A_terms*n_bfs + q_a*n_bfs + i;
            set_scalar_in_list(
                               fq_aq_innerprods_list, offset,
                               rb_eval.Fq_Aq_representor_innerprods[q_f][q_a][i]);
          }
  }

  // AqAq representor inner-product data
  {
    unsigned int Q_a_hat = n_A_terms*(n_A_terms+1)/2;
    auto aq_aq_innerprods_list =
      rb_evaluation_builder.initAqAqInnerprods(Q_a_hat*n_bfs*n_bfs);

    for(unsigned int i=0; i < Q_a_hat; ++i)
      for(unsigned int j=0; j < n_bfs; ++j)
        for(unsigned int l=0; l < n_bfs; ++l)
          {
            unsigned int offset = i*n_bfs*n_bfs + j*n_bfs + l;
            set_scalar_in_list(
                               aq_aq_innerprods_list,
                               offset,
                               rb_eval.Aq_Aq_representor_innerprods[i][j][l]);
          }
  }

  // Output dual inner-product data, and output vectors
  {
    unsigned int n_outputs = rb_theta_expansion.get_n_outputs();
    auto output_innerprod_outer = rb_evaluation_builder.initOutputDualInnerprods(n_outputs);
    auto output_vector_outer = rb_evaluation_builder.initOutputVectors(n_outputs);

    for(unsigned int output_id=0; output_id < n_outputs; ++output_id)
      {
        unsigned int n_output_terms = rb_theta_expansion.get_n_output_terms(output_id);

        {
          unsigned int Q_l_hat = n_output_terms*(n_output_terms+1)/2;
          auto output_innerprod_inner = output_innerprod_outer.init(output_id, Q_l_hat);
          for(unsigned int q=0; q < Q_l_hat; ++q)
            {
              set_scalar_in_list(
                                 output_innerprod_inner, q, rb_eval.output_dual_innerprods[output_id][q]);
            }
        }

        {
          auto output_vector_middle = output_vector_outer.init(output_id, n_output_terms);
          for(unsigned int q_l=0; q_l<n_output_terms; ++q_l)
            {
              auto output_vector_inner = output_vector_middle.init(q_l, n_bfs);
              for(unsigned int j=0; j<n_bfs; ++j)
                {
                  set_scalar_in_list(
                                     output_vector_inner, j, rb_eval.RB_output_vectors[output_id][q_l](j));
                }
            }
        }
      }
  }

  // Fq vectors and Aq matrices
  {
    unsigned int n_F_terms = rb_theta_expansion.get_n_F_terms();
    unsigned int n_A_terms = rb_theta_expansion.get_n_A_terms();

    auto rb_fq_vectors_outer_list = rb_evaluation_builder.initRbFqVectors(n_F_terms);
    for(unsigned int q_f=0; q_f < n_F_terms; ++q_f)
      {
        auto rb_fq_vectors_inner_list = rb_fq_vectors_outer_list.init(q_f, n_bfs);
        for(unsigned int i=0; i<n_bfs; i++)
          set_scalar_in_list(rb_fq_vectors_inner_list, i, rb_eval.RB_Fq_vector[q_f](i));
      }

    auto rb_Aq_matrices_outer_list = rb_evaluation_builder.initRbAqMatrices(n_A_terms);
    for(unsigned int q_a=0; q_a < n_A_terms; ++q_a)
      {
        auto rb_Aq_matrices_inner_list = rb_Aq_matrices_outer_list.init(q_a, n_bfs*n_bfs);
        for(unsigned int i=0; i < n_bfs; ++i)
          for(unsigned int j=0; j < n_bfs; ++j)
            {
              unsigned int offset = i*n_bfs+j;
              set_scalar_in_list(rb_Aq_matrices_inner_list, offset, rb_eval.RB_Aq_vector[q_a](i,j));
            }
      }
  }

  // Inner-product matrix
  if(rb_eval.compute_RB_inner_product)
    {
      auto rb_inner_product_matrix_list =
        rb_evaluation_builder.initRbInnerProductMatrix(n_bfs*n_bfs);

      for(unsigned int i=0; i < n_bfs; ++i)
        for(unsigned int j=0; j < n_bfs; ++j)
          {
            unsigned int offset = i*n_bfs + j;
            set_scalar_in_list(
                               rb_inner_product_matrix_list,
                               offset,
                               rb_eval.RB_inner_product_matrix(i,j) );
          }
    }

  auto parameter_ranges_list =
    rb_evaluation_builder.initParameterRanges();
  auto discrete_parameters_list =
    rb_evaluation_builder.initDiscreteParameters();
  add_parameter_ranges_to_builder(rb_eval,
                                  parameter_ranges_list,
                                  discrete_parameters_list);
}

template <typename RBEvaluationBuilderNumber, typename TransRBEvaluationBuilderNumber>
void add_transient_rb_evaluation_data_to_builder(TransientRBEvaluation & trans_rb_eval,
                                                 RBEvaluationBuilderNumber & rb_eval_builder,
                                                 TransRBEvaluationBuilderNumber & trans_rb_eval_builder)
{
  add_rb_evaluation_data_to_builder(trans_rb_eval, rb_eval_builder);

  trans_rb_eval_builder.setDeltaT(trans_rb_eval.get_delta_t());
  trans_rb_eval_builder.setEulerTheta(trans_rb_eval.get_euler_theta());
  trans_rb_eval_builder.setNTimeSteps(trans_rb_eval.get_n_time_steps());
  trans_rb_eval_builder.setTimeStep(trans_rb_eval.get_time_step());

  unsigned int n_bfs = trans_rb_eval.get_n_basis_functions();

  // L2-inner-product matrix
  {
    auto rb_L2_matrix_list =
      trans_rb_eval_builder.initRbL2Matrix(n_bfs*n_bfs);

    for(unsigned int i=0; i<n_bfs; ++i)
      for(unsigned int j=0; j<n_bfs; ++j)
        {
          unsigned int offset = i*n_bfs + j;
          set_scalar_in_list(rb_L2_matrix_list,
                             offset,
                             trans_rb_eval.RB_L2_matrix(i,j));
        }
  }

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(trans_rb_eval.get_rb_theta_expansion());
  unsigned int n_M_terms = trans_theta_expansion.get_n_M_terms();
  // Mq matrices
  {
    auto rb_Mq_matrices_outer_list = trans_rb_eval_builder.initRbMqMatrices(n_M_terms);
    for(unsigned int q_m=0; q_m < n_M_terms; ++q_m)
      {
        auto rb_Mq_matrices_inner_list = rb_Mq_matrices_outer_list.init(q_m, n_bfs*n_bfs);
        for(unsigned int i=0; i < n_bfs; ++i)
          for(unsigned int j=0; j < n_bfs; ++j)
            {
              unsigned int offset = i*n_bfs+j;
              set_scalar_in_list(rb_Mq_matrices_inner_list,
                                 offset,
                                 trans_rb_eval.RB_M_q_vector[q_m](i,j));
            }
      }
  }

  // The initial condition and L2 error at t=0.
  // We store the values for each RB space of dimension (0,...,n_basis_functions).
  {
    auto initial_l2_errors_builder =
      trans_rb_eval_builder.initInitialL2Errors(n_bfs);
    auto initial_conditions_outer_list =
      trans_rb_eval_builder.initInitialConditions(n_bfs);

    for(unsigned int i=0; i<n_bfs; i++)
      {
        initial_l2_errors_builder.set(i, trans_rb_eval.initial_L2_error_all_N[i]);

        auto initial_conditions_inner_list =
          initial_conditions_outer_list.init(i, i+1);
        for(unsigned int j=0; j<=i; j++)
          {
            set_scalar_in_list(initial_conditions_inner_list,
                               j,
                               trans_rb_eval.RB_initial_condition_all_N[i](j));
          }
      }
  }

  // FqMq representor inner-product data
  {
    unsigned int n_F_terms = trans_theta_expansion.get_n_F_terms();
    auto fq_mq_innerprods_list =
      trans_rb_eval_builder.initFqMqInnerprods(n_F_terms*n_M_terms*n_bfs);

    for(unsigned int q_f=0; q_f<n_F_terms; ++q_f)
      for(unsigned int q_m=0; q_m<n_M_terms; ++q_m)
        for(unsigned int i=0; i<n_bfs; ++i)
          {
            unsigned int offset = q_f*n_M_terms*n_bfs + q_m*n_bfs + i;
            set_scalar_in_list(fq_mq_innerprods_list,
                               offset,
                               trans_rb_eval.Fq_Mq_representor_innerprods[q_f][q_m][i]);
          }
  }

  // MqMq representor inner-product data
  {
    unsigned int Q_m_hat = n_M_terms*(n_M_terms+1)/2;
    auto mq_mq_innerprods_list =
      trans_rb_eval_builder.initMqMqInnerprods(Q_m_hat*n_bfs*n_bfs);

    for(unsigned int i=0; i < Q_m_hat; ++i)
      for(unsigned int j=0; j < n_bfs; ++j)
        for(unsigned int l=0; l < n_bfs; ++l)
          {
            unsigned int offset = i*n_bfs*n_bfs + j*n_bfs + l;
            set_scalar_in_list(mq_mq_innerprods_list,
                               offset,
                               trans_rb_eval.Mq_Mq_representor_innerprods[i][j][l]);
          }
  }

  // AqMq representor inner-product data
  {
    unsigned int n_A_terms = trans_theta_expansion.get_n_A_terms();

    auto aq_mq_innerprods_list =
      trans_rb_eval_builder.initAqMqInnerprods(n_A_terms*n_M_terms*n_bfs*n_bfs);

    for(unsigned int q_a=0; q_a<n_A_terms; q_a++)
      for(unsigned int q_m=0; q_m<n_M_terms; q_m++)
        for(unsigned int i=0; i<n_bfs; i++)
          for(unsigned int j=0; j<n_bfs; j++)
            {
              unsigned int offset =
                q_a*(n_M_terms*n_bfs*n_bfs) + q_m*(n_bfs*n_bfs) + i*n_bfs + j;
              set_scalar_in_list(aq_mq_innerprods_list,
                                 offset,
                                 trans_rb_eval.Aq_Mq_representor_innerprods[q_a][q_m][i][j]);
            }
  }

}

template <typename RBEvaluationBuilderNumber, typename RBEIMEvaluationBuilderNumber>
void add_rb_eim_evaluation_data_to_builder(RBEIMEvaluation & rb_eim_evaluation,
                                           RBEvaluationBuilderNumber & rb_evaluation_builder,
                                           RBEIMEvaluationBuilderNumber & rb_eim_evaluation_builder)
{
  add_rb_evaluation_data_to_builder(rb_eim_evaluation, rb_evaluation_builder);

  unsigned int n_bfs = rb_eim_evaluation.get_n_basis_functions();

  // EIM interpolation matrix
  {
    // We store the lower triangular part of an NxN matrix, the size of which is given by
    // (N(N + 1))/2
    unsigned int half_matrix_size = n_bfs*(n_bfs+1)/2;

    auto interpolation_matrix_list =
      rb_eim_evaluation_builder.initInterpolationMatrix(half_matrix_size);
    for(unsigned int i=0; i < n_bfs; ++i)
      for(unsigned int j=0; j <= i; ++j)
        {
          unsigned int offset = i*(i+1)/2 + j;
          set_scalar_in_list(interpolation_matrix_list,
                             offset,
                             rb_eim_evaluation.interpolation_matrix(i,j));
        }
  }

  // Interpolation points
  {
    auto interpolation_points_list =
      rb_eim_evaluation_builder.initInterpolationPoints(n_bfs);
    for(unsigned int i=0; i < n_bfs; ++i)
      add_point_to_builder(rb_eim_evaluation.interpolation_points[i],
                           interpolation_points_list[i]);
  }

  // Interpolation points variables
  {
    auto interpolation_points_var_list =
      rb_eim_evaluation_builder.initInterpolationPointsVar(n_bfs);
    for(unsigned int i=0; i<n_bfs; ++i)
      interpolation_points_var_list.set(i,
                                        rb_eim_evaluation.interpolation_points_var[i]);
  }

  // Interpolation elements
  {
    unsigned int n_interpolation_elems =
      rb_eim_evaluation.interpolation_points_elem.size();
    auto interpolation_points_elem_list =
      rb_eim_evaluation_builder.initInterpolationPointsElems(n_interpolation_elems);

    if (n_interpolation_elems != n_bfs)
      libmesh_error_msg("The number of elements should match the number of basis functions");

    for(unsigned int i=0; i<n_interpolation_elems; ++i)
      {
        const libMesh::Elem & elem = *rb_eim_evaluation.interpolation_points_elem[i];
        auto mesh_elem_builder = interpolation_points_elem_list[i];
        add_elem_to_builder(elem, mesh_elem_builder);
      }
  }
}

#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)
void add_rb_scm_evaluation_data_to_builder(RBSCMEvaluation & rb_scm_eval,
                                           RBData::RBSCMEvaluation::Builder & rb_scm_eval_builder)
{
  auto parameter_ranges_list =
    rb_scm_eval_builder.initParameterRanges();
  auto discrete_parameters_list =
    rb_scm_eval_builder.initDiscreteParameters();
  add_parameter_ranges_to_builder(rb_scm_eval,
                                  parameter_ranges_list,
                                  discrete_parameters_list);

  {
    if (rb_scm_eval.B_min.size() != rb_scm_eval.get_rb_theta_expansion().get_n_A_terms())
      libmesh_error_msg("Size error while writing B_min");
    auto b_min_list = rb_scm_eval_builder.initBMin( rb_scm_eval.B_min.size() );
    for (std::size_t i=0; i<rb_scm_eval.B_min.size(); i++)
      b_min_list.set(i, rb_scm_eval.get_B_min(i));
  }

  {
    if (rb_scm_eval.B_max.size() != rb_scm_eval.get_rb_theta_expansion().get_n_A_terms())
      libmesh_error_msg("Size error while writing B_max");

    auto b_max_list = rb_scm_eval_builder.initBMax( rb_scm_eval.B_max.size() );
    for (std::size_t i=0; i<rb_scm_eval.B_max.size(); i++)
      b_max_list.set(i, rb_scm_eval.get_B_max(i));
  }

  {
    auto cj_stability_vector =
      rb_scm_eval_builder.initCJStabilityVector( rb_scm_eval.C_J_stability_vector.size() );
    for (std::size_t i=0; i<rb_scm_eval.C_J_stability_vector.size(); i++)
      cj_stability_vector.set(i, rb_scm_eval.get_C_J_stability_constraint(i));
  }

  {
    auto cj_parameters_outer =
      rb_scm_eval_builder.initCJ( rb_scm_eval.C_J.size() );

    for (std::size_t i=0; i<rb_scm_eval.C_J.size(); i++)
      {
        auto cj_parameters_inner =
          cj_parameters_outer.init(i, rb_scm_eval.C_J[i].n_parameters());

        RBParameters::const_iterator it     = rb_scm_eval.C_J[i].begin();
        RBParameters::const_iterator it_end = rb_scm_eval.C_J[i].end();

        unsigned int count = 0;
        for( ; it != it_end; ++it)
          {
            cj_parameters_inner[count].setName( it->first );
            cj_parameters_inner[count].setValue( it->second );
            count++;
          }

      }
  }

  {
    unsigned int n_C_J_values = rb_scm_eval.C_J.size();
    unsigned int n_A_terms = rb_scm_eval.get_rb_theta_expansion().get_n_A_terms();
    unsigned int n_values = n_C_J_values*n_A_terms;
    auto scm_ub_vectors =
      rb_scm_eval_builder.initScmUbVectors( n_values );

    for(unsigned int i=0; i<n_C_J_values; i++)
      for(unsigned int j=0; j<n_A_terms; j++)
        {
          unsigned int offset = i*n_A_terms + j;
          scm_ub_vectors.set(offset, rb_scm_eval.get_SCM_UB_vector(i,j));
        }
  }
}
#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

void add_point_to_builder(const Point & point, RBData::Point3D::Builder point_builder)
{
  point_builder.setX(point(0));

  if (LIBMESH_DIM >= 2)
    point_builder.setY(point(1));

  if (LIBMESH_DIM >= 3)
    point_builder.setZ(point(2));
}

void add_elem_to_builder(const libMesh::Elem & elem, RBData::MeshElem::Builder mesh_elem_builder)
{
  std::string elem_type_string = libMesh::Utility::enum_to_string(elem.type());

  mesh_elem_builder.setType(elem_type_string.c_str());
  mesh_elem_builder.setSubdomainId(elem.subdomain_id());

  unsigned int n_points = elem.n_nodes();
  auto mesh_elem_point_list = mesh_elem_builder.initPoints(n_points);

  for(unsigned int j=0; j < n_points; ++j)
    {
      add_point_to_builder(elem.node_ref(j), mesh_elem_point_list[j]);
    }
}

// ---- Helper functions for adding data to capnp Builders (END) ----

} // namespace RBDataSerialization

} // namespace libMesh

#endif // #if defined(LIBMESH_HAVE_CAPNPROTO)
