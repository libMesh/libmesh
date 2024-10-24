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

// libMesh includes
#include "libmesh/rb_data_serialization.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/transient_rb_theta_expansion.h"
#include "libmesh/rb_evaluation.h"
#include "libmesh/transient_rb_evaluation.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_scm_evaluation.h"
#include "libmesh/elem.h"
#include "libmesh/int_range.h"
#include "libmesh/rb_parametrized_function.h"
#include "libmesh/libmesh_logging.h"

// Cap'n'Proto includes
#include <capnp/serialize.h>

// C++ includes
#include <iostream>
#include <fstream>
#ifdef LIBMESH_HAVE_UNISTD_H
#include <unistd.h> // for close()
#endif
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

RBEvaluationSerialization::~RBEvaluationSerialization() = default;

void RBEvaluationSerialization::write_to_file(const std::string & path)
{
  LOG_SCOPE("write_to_file()", "RBEvaluationSerialization");

  if (_rb_eval.comm().rank() == 0)
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
      libmesh_error_msg_if(!fd, "Error opening a write-only file descriptor to " + path);

      capnp::writeMessageToFd(fd, message);

      int error = close(fd);
      libmesh_error_msg_if(error, "Error closing a write-only file descriptor to " + path);
    }
}

// ---- RBEvaluationSerialization (END) ----


// ---- TransientRBEvaluationSerialization (BEGIN) ----

TransientRBEvaluationSerialization::
TransientRBEvaluationSerialization(TransientRBEvaluation & trans_rb_eval) :
  _trans_rb_eval(trans_rb_eval)
{
}

TransientRBEvaluationSerialization::~TransientRBEvaluationSerialization() = default;

void TransientRBEvaluationSerialization::write_to_file(const std::string & path)
{
  LOG_SCOPE("write_to_file()", "TransientRBEvaluationSerialization");

  if (_trans_rb_eval.comm().rank() == 0)
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
      libmesh_error_msg_if(!fd, "Error opening a write-only file descriptor to " + path);

      capnp::writeMessageToFd(fd, message);

      int error = close(fd);
      libmesh_error_msg_if(error, "Error closing a write-only file descriptor to " + path);
    }
}

// ---- TransientRBEvaluationSerialization (END) ----


// ---- RBEIMEvaluationSerialization (BEGIN) ----

RBEIMEvaluationSerialization::RBEIMEvaluationSerialization(RBEIMEvaluation & rb_eim_eval)
  :
  _rb_eim_eval(rb_eim_eval)
{
}

RBEIMEvaluationSerialization::~RBEIMEvaluationSerialization() = default;

void RBEIMEvaluationSerialization::write_to_file(const std::string & path)
{
  LOG_SCOPE("write_to_file()", "RBEIMEvaluationSerialization");

  if (_rb_eim_eval.comm().rank() == 0)
    {
      capnp::MallocMessageBuilder message;

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
      RBData::RBEIMEvaluationReal::Builder rb_eim_eval_builder =
        message.initRoot<RBData::RBEIMEvaluationReal>();
#else
      RBData::RBEIMEvaluationComplex::Builder rb_eim_eval_builder =
        message.initRoot<RBData::RBEIMEvaluationComplex>();
#endif

      add_rb_eim_evaluation_data_to_builder(_rb_eim_eval,
                                            rb_eim_eval_builder);

      int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
      libmesh_error_msg_if(!fd, "Error opening a write-only file descriptor to " + path);

      capnp::writeMessageToFd(fd, message);

      int error = close(fd);
      libmesh_error_msg_if(error, "Error closing a write-only file descriptor to " + path);
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

RBSCMEvaluationSerialization::~RBSCMEvaluationSerialization() = default;

void RBSCMEvaluationSerialization::write_to_file(const std::string & path)
{
  LOG_SCOPE("write_to_file()", "RBSCMEvaluationSerialization");

  if (_rb_scm_eval.comm().rank() == 0)
    {
      capnp::MallocMessageBuilder message;

      RBData::RBSCMEvaluation::Builder rb_scm_eval_builder =
        message.initRoot<RBData::RBSCMEvaluation>();

      add_rb_scm_evaluation_data_to_builder(_rb_scm_eval, rb_scm_eval_builder);

      int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
      libmesh_error_msg_if(!fd, "Error opening a write-only file descriptor to " + path);

      capnp::writeMessageToFd(fd, message);

      int error = close(fd);
      libmesh_error_msg_if(error, "Error closing a write-only file descriptor to " + path);
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

    const RBParameters & parameters_min = rb_evaluation.get_parameters_min();
    const RBParameters & parameters_max = rb_evaluation.get_parameters_max();

    // We could loop over either parameters_min or parameters_max, they should have the same keys.
    unsigned int count = 0;
    for (const auto & [key, val] : parameters_min)
      if (!rb_evaluation.is_discrete_parameter(key))
        {
          names.set(count, key);
          mins.set(count, parameters_min.get_value(key));
          maxs.set(count, parameters_max.get_value(key));
          ++count;
        }

    libmesh_error_msg_if(count != n_continuous_parameters, "Mismatch in number of continuous parameters");
  }

  // Discrete parameters
  {
    unsigned int n_discrete_parameters = rb_evaluation.get_n_discrete_params();
    auto names = discrete_parameters_list.initNames(n_discrete_parameters);
    auto values_outer = discrete_parameters_list.initValues(n_discrete_parameters);

    const std::map<std::string, std::vector<Real>> & discrete_parameters =
      rb_evaluation.get_discrete_parameter_values();

    unsigned int count = 0;
    for (const auto & discrete_parameter : discrete_parameters)
      {
        names.set(count, discrete_parameter.first);

        const std::vector<Real> & values = discrete_parameter.second;
        unsigned int n_values = values.size();

        values_outer.init(count, n_values);
        auto values_inner = values_outer[count];
        for (unsigned int i=0; i<n_values; ++i)
          {
            values_inner.set(i, values[i]);
          }

        ++count;
      }

    libmesh_error_msg_if(count != n_discrete_parameters, "Mismatch in number of discrete parameters");
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

    for (unsigned int i=0; i<Q_f_hat; i++)
      set_scalar_in_list(fq_innerprods_list,
                         i,
                         rb_eval.Fq_representor_innerprods[i]);
  }

  // FqAq representor inner-product data
  {
    auto fq_aq_innerprods_list =
      rb_evaluation_builder.initFqAqInnerprods(n_F_terms*n_A_terms*n_bfs);

    for (unsigned int q_f=0; q_f < n_F_terms; ++q_f)
      for (unsigned int q_a=0; q_a < n_A_terms; ++q_a)
        for (unsigned int i=0; i < n_bfs; ++i)
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

    for (unsigned int i=0; i < Q_a_hat; ++i)
      for (unsigned int j=0; j < n_bfs; ++j)
        for (unsigned int l=0; l < n_bfs; ++l)
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

    for (unsigned int output_id=0; output_id < n_outputs; ++output_id)
      {
        unsigned int n_output_terms = rb_theta_expansion.get_n_output_terms(output_id);

        {
          unsigned int Q_l_hat = n_output_terms*(n_output_terms+1)/2;
          auto output_innerprod_inner = output_innerprod_outer.init(output_id, Q_l_hat);
          for (unsigned int q=0; q < Q_l_hat; ++q)
            {
              set_scalar_in_list(
                                 output_innerprod_inner, q, rb_eval.output_dual_innerprods[output_id][q]);
            }
        }

        {
          auto output_vector_middle = output_vector_outer.init(output_id, n_output_terms);
          for (unsigned int q_l=0; q_l<n_output_terms; ++q_l)
            {
              auto output_vector_inner = output_vector_middle.init(q_l, n_bfs);
              for (unsigned int j=0; j<n_bfs; ++j)
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
    for (unsigned int q_f=0; q_f < n_F_terms; ++q_f)
      {
        auto rb_fq_vectors_inner_list = rb_fq_vectors_outer_list.init(q_f, n_bfs);
        for (unsigned int i=0; i<n_bfs; i++)
          set_scalar_in_list(rb_fq_vectors_inner_list, i, rb_eval.RB_Fq_vector[q_f](i));
      }

    auto rb_Aq_matrices_outer_list = rb_evaluation_builder.initRbAqMatrices(n_A_terms);
    for (unsigned int q_a=0; q_a < n_A_terms; ++q_a)
      {
        auto rb_Aq_matrices_inner_list = rb_Aq_matrices_outer_list.init(q_a, n_bfs*n_bfs);
        for (unsigned int i=0; i < n_bfs; ++i)
          for (unsigned int j=0; j < n_bfs; ++j)
            {
              unsigned int offset = i*n_bfs+j;
              set_scalar_in_list(rb_Aq_matrices_inner_list, offset, rb_eval.RB_Aq_vector[q_a](i,j));
            }
      }
  }

  // Inner-product matrix
  if (rb_eval.compute_RB_inner_product)
    {
      auto rb_inner_product_matrix_list =
        rb_evaluation_builder.initRbInnerProductMatrix(n_bfs*n_bfs);

      for (unsigned int i=0; i < n_bfs; ++i)
        for (unsigned int j=0; j < n_bfs; ++j)
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

    for (unsigned int i=0; i<n_bfs; ++i)
      for (unsigned int j=0; j<n_bfs; ++j)
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
    for (unsigned int q_m=0; q_m < n_M_terms; ++q_m)
      {
        auto rb_Mq_matrices_inner_list = rb_Mq_matrices_outer_list.init(q_m, n_bfs*n_bfs);
        for (unsigned int i=0; i < n_bfs; ++i)
          for (unsigned int j=0; j < n_bfs; ++j)
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

    for (unsigned int i=0; i<n_bfs; i++)
      {
        initial_l2_errors_builder.set(i, trans_rb_eval.initial_L2_error_all_N[i]);

        auto initial_conditions_inner_list =
          initial_conditions_outer_list.init(i, i+1);
        for (unsigned int j=0; j<=i; j++)
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

    for (unsigned int q_f=0; q_f<n_F_terms; ++q_f)
      for (unsigned int q_m=0; q_m<n_M_terms; ++q_m)
        for (unsigned int i=0; i<n_bfs; ++i)
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

    for (unsigned int i=0; i < Q_m_hat; ++i)
      for (unsigned int j=0; j < n_bfs; ++j)
        for (unsigned int l=0; l < n_bfs; ++l)
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

    for (unsigned int q_a=0; q_a<n_A_terms; q_a++)
      for (unsigned int q_m=0; q_m<n_M_terms; q_m++)
        for (unsigned int i=0; i<n_bfs; i++)
          for (unsigned int j=0; j<n_bfs; j++)
            {
              unsigned int offset =
                q_a*(n_M_terms*n_bfs*n_bfs) + q_m*(n_bfs*n_bfs) + i*n_bfs + j;
              set_scalar_in_list(aq_mq_innerprods_list,
                                 offset,
                                 trans_rb_eval.Aq_Mq_representor_innerprods[q_a][q_m][i][j]);
            }
  }

}

template <typename RBEIMEvaluationBuilderNumber>
void add_rb_eim_evaluation_data_to_builder(RBEIMEvaluation & rb_eim_evaluation,
                                           RBEIMEvaluationBuilderNumber & rb_eim_evaluation_builder)
{
  // Number of basis functions
  unsigned int n_bfs = rb_eim_evaluation.get_n_basis_functions();
  rb_eim_evaluation_builder.setNBfs(n_bfs);

  // EIM interpolation matrix
  {
    // We store the lower triangular part of an NxN matrix, the size of which is given by
    // (N(N + 1))/2
    unsigned int half_matrix_size = n_bfs*(n_bfs+1)/2;

    auto interpolation_matrix_list =
      rb_eim_evaluation_builder.initInterpolationMatrix(half_matrix_size);
    for (unsigned int i=0; i < n_bfs; ++i)
      for (unsigned int j=0; j <= i; ++j)
        {
          unsigned int offset = i*(i+1)/2 + j;
          set_scalar_in_list(interpolation_matrix_list,
                             offset,
                             rb_eim_evaluation.get_interpolation_matrix()(i,j));
        }
  }

  // We check use_eim_error_indicator() and the number of interpolation
  // points in order to set use_error_indicator. This is because even if
  // use_error_indicator() is true, we may not be using an error indicator,
  // e.g. if there are no parameters in the model then we will not use an
  // error indicator, and we can detect this by comparing the number of
  // interpolation points with n_bfs.
  bool use_error_indicator =
    (rb_eim_evaluation.use_eim_error_indicator() &&
    (rb_eim_evaluation.get_n_interpolation_points() > n_bfs));

  if (use_error_indicator)
  {
    auto error_indicator_data_list =
      rb_eim_evaluation_builder.initEimErrorIndicatorInterpData(n_bfs);
    const auto & error_indicator_row =
      rb_eim_evaluation.get_error_indicator_interpolation_row();
    for (unsigned int i=0; i < n_bfs; ++i)
      {
        set_scalar_in_list(error_indicator_data_list,
                           i,
                           error_indicator_row(i));
      }
  }

  // If we're using the EIM error indicator then we store one extra
  // interpolation point and associated data, hence we increment n_bfs
  // here so that we write out the extra data below.
  //
  // However the interpolation matrix and _error_indicator_interpolation_row
  // does not include this extra point, which is why we write it out above
  // before n_bfs is incremented.
  if (use_error_indicator)
    n_bfs++;

  libmesh_error_msg_if(n_bfs != rb_eim_evaluation.get_n_interpolation_points(),
    "Number of basis functions should match number of interpolation points");

  auto parameter_ranges_list =
    rb_eim_evaluation_builder.initParameterRanges();
  auto discrete_parameters_list =
    rb_eim_evaluation_builder.initDiscreteParameters();
  add_parameter_ranges_to_builder(rb_eim_evaluation,
                                  parameter_ranges_list,
                                  discrete_parameters_list);

  // Interpolation points
  {
    auto interpolation_points_list =
      rb_eim_evaluation_builder.initInterpolationXyz(n_bfs);
    for (unsigned int i=0; i < n_bfs; ++i)
      add_point_to_builder(rb_eim_evaluation.get_interpolation_points_xyz(i),
                           interpolation_points_list[i]);
  }

  // Interpolation points comps
  {
    auto interpolation_points_comp_list =
      rb_eim_evaluation_builder.initInterpolationComp(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      interpolation_points_comp_list.set(i,
                                         rb_eim_evaluation.get_interpolation_points_comp(i));
  }

  // Interpolation points subdomain IDs
  {
    auto interpolation_points_subdomain_id_list =
      rb_eim_evaluation_builder.initInterpolationSubdomainId(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      interpolation_points_subdomain_id_list.set(i,
                                                 rb_eim_evaluation.get_interpolation_points_subdomain_id(i));
  }

  // Interpolation points boundary IDs, relevant if the parametrized function is defined
  // on mesh sides or nodesets
  if (rb_eim_evaluation.get_parametrized_function().on_mesh_sides() ||
      rb_eim_evaluation.get_parametrized_function().on_mesh_nodes())
  {
    auto interpolation_points_boundary_id_list =
      rb_eim_evaluation_builder.initInterpolationBoundaryId(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      interpolation_points_boundary_id_list.set(i,
                                                 rb_eim_evaluation.get_interpolation_points_boundary_id(i));
  }

  // Interpolation points element IDs
  {
    auto interpolation_points_elem_id_list =
      rb_eim_evaluation_builder.initInterpolationElemId(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      interpolation_points_elem_id_list.set(i,
                                            rb_eim_evaluation.get_interpolation_points_elem_id(i));
  }

  // Interpolation points node IDs, relevant if the parametrized function is defined on mesh sides
  if (rb_eim_evaluation.get_parametrized_function().on_mesh_nodes())
  {
    auto interpolation_points_node_id_list =
      rb_eim_evaluation_builder.initInterpolationNodeId(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      interpolation_points_node_id_list.set(i,
                                            rb_eim_evaluation.get_interpolation_points_node_id(i));
  }

  // Interpolation points side indices, relevant if the parametrized function is defined on mesh sides
  if (rb_eim_evaluation.get_parametrized_function().on_mesh_sides())
  {
    auto interpolation_points_side_index_list =
      rb_eim_evaluation_builder.initInterpolationSideIndex(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      interpolation_points_side_index_list.set(i,
                                               rb_eim_evaluation.get_interpolation_points_side_index(i));
  }

  // Interpolation points quadrature point indices
  {
    auto interpolation_points_qp_list =
      rb_eim_evaluation_builder.initInterpolationQp(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      interpolation_points_qp_list.set(i,
                                       rb_eim_evaluation.get_interpolation_points_qp(i));
  }

  // Interpolation points perturbations
  {
    auto interpolation_points_list_outer =
      rb_eim_evaluation_builder.initInterpolationXyzPerturb(n_bfs);
    for (unsigned int i=0; i < n_bfs; ++i)
      {
        const std::vector<Point> & perturbs = rb_eim_evaluation.get_interpolation_points_xyz_perturbations(i);
        auto interpolation_points_list_inner = interpolation_points_list_outer.init(i, perturbs.size());

        for (unsigned int j : index_range(perturbs))
          {
            add_point_to_builder(perturbs[j], interpolation_points_list_inner[j]);
          }
      }
  }

  unsigned int n_elems = rb_eim_evaluation.get_n_elems();
  rb_eim_evaluation_builder.setNElems(n_elems);
  // Elem id to local index map
  {
    auto elem_id_to_local_index_list =
      rb_eim_evaluation_builder.initElemIdToLocalIndex(n_elems);
    unsigned int counter = 0;
    for (auto const & [elem_id, local_index] : rb_eim_evaluation.get_elem_id_to_local_index_map())
      {
        elem_id_to_local_index_list[counter].setFirst(elem_id);
        elem_id_to_local_index_list[counter].setSecond(local_index);
        counter++;
      }
  }

  // Interpolation points JxW values at each qp
  {
    auto interpolation_points_list_outer =
      rb_eim_evaluation_builder.initInterpolationJxWAllQp(n_elems);
    for (unsigned int i=0; i < n_elems; ++i)
      {
        const std::vector<Real> & JxW = rb_eim_evaluation.get_interpolation_points_JxW_all_qp(i);
        auto interpolation_points_list_inner = interpolation_points_list_outer.init(i, JxW.size());

        for (unsigned int j : index_range(JxW))
          {
            // Here we can use set() instead of set_scalar_in_list() because
            // phi stores real-valued data only.
            interpolation_points_list_inner.set(j, JxW[j]);
          }
      }
  }

  // Interpolation points phi values at each qp
  {
    auto interpolation_points_list_outer =
      rb_eim_evaluation_builder.initInterpolationPhiValuesAllQp(n_elems);

    for (unsigned int i=0; i < n_elems; ++i)
      {
        const auto & phi_i_all_qp = rb_eim_evaluation.get_interpolation_points_phi_i_all_qp(i);
        auto interpolation_points_list_middle = interpolation_points_list_outer.init(i, phi_i_all_qp.size());

        for (unsigned int j : index_range(phi_i_all_qp))
          {
            auto interpolation_points_list_inner = interpolation_points_list_middle.init(j, phi_i_all_qp[j].size());

            for (auto k : index_range(phi_i_all_qp[j]))
                interpolation_points_list_inner.set(k, phi_i_all_qp[j][k]);
          }
      }
  }

  // Dxyzdxi at the element center for the element that contains each interpolation point.
  {
    auto dxyzdxi_elem_center =
      rb_eim_evaluation_builder.initInterpolationDxyzDxiElem(n_elems);
    for (auto i : make_range(n_elems))
      {
        add_point_to_builder(rb_eim_evaluation.get_elem_center_dxyzdxi(i), dxyzdxi_elem_center[i]);
      }
  }

  // Dxyzdeta at the element center for the element that contains each interpolation point.
  {
    auto dxyzdeta_elem_center =
      rb_eim_evaluation_builder.initInterpolationDxyzDetaElem(n_elems);
    for (auto i : make_range(n_elems))
      {
        add_point_to_builder(rb_eim_evaluation.get_elem_center_dxyzdeta(i), dxyzdeta_elem_center[i]);
      }
  }

  // Quadrature rule order associated to the element that contains each interpolation point.
  {
    auto interpolation_points_qrule_order_list =
      rb_eim_evaluation_builder.initInterpolationQruleOrder(n_elems);
    for (unsigned int i=0; i<n_elems; ++i)
      interpolation_points_qrule_order_list.set(i,
                                                rb_eim_evaluation.get_interpolation_points_qrule_order(i));
  }

  // Element type for the element that contains each interpolation point
  {
    auto interpolation_points_elem_type_list =
      rb_eim_evaluation_builder.initInterpolationElemType(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      interpolation_points_elem_type_list.set(i,
                                              rb_eim_evaluation.get_interpolation_points_elem_type(i));
  }

  // Optionally store EIM solutions for the training set
  if (rb_eim_evaluation.get_parametrized_function().is_lookup_table)
    {
      const std::vector<DenseVector<Number>> & eim_solutions = rb_eim_evaluation.get_eim_solutions_for_training_set();

      auto eim_rhs_list_outer =
        rb_eim_evaluation_builder.initEimSolutionsForTrainingSet(eim_solutions.size());
      for (auto i : make_range(eim_solutions.size()))
        {
          const DenseVector<Number> & values = eim_solutions[i];
          auto eim_rhs_list_inner = eim_rhs_list_outer.init(i, values.size());

          for (auto j : index_range(values))
            {
              set_scalar_in_list(eim_rhs_list_inner,
                                 j,
                                 values(j));
            }
        }
    }

  // The shape function values at the interpolation points. This can be used to evaluate nodal data
  // at EIM interpolation points, which are at quadrature points.
  {
    auto interpolation_points_list_outer =
      rb_eim_evaluation_builder.initInterpolationPhiValues(n_bfs);
    for (unsigned int i=0; i < n_bfs; ++i)
      {
        const std::vector<Real> & phi_i_qp_vec = rb_eim_evaluation.get_interpolation_points_phi_i_qp(i);
        auto interpolation_points_list_inner = interpolation_points_list_outer.init(i, phi_i_qp_vec.size());

        for (unsigned int j : index_range(phi_i_qp_vec))
          {
            // Here we can use set() instead of set_scalar_in_list() because
            // phi stores real-valued data only.
            interpolation_points_list_inner.set(j, phi_i_qp_vec[j]);
          }
      }
  }

  // Initialize EIM spatial indices for the interpolation points, and if there are
  // any spatial indices then we store them in the buffer.
  rb_eim_evaluation.initialize_interpolation_points_spatial_indices();
  if (rb_eim_evaluation.get_n_interpolation_points_spatial_indices() > 0)
  {
    libmesh_error_msg_if(n_bfs != rb_eim_evaluation.get_n_interpolation_points_spatial_indices(),
                         "Error: Number of spatial indices should match number of EIM basis functions");

    auto interpolation_points_spatial_indices_builder =
      rb_eim_evaluation_builder.initInterpolationSpatialIndices(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      {
        const std::vector<unsigned int> & spatial_indices =
          rb_eim_evaluation.get_interpolation_points_spatial_indices(i);
        unsigned int n_spatial_indices = spatial_indices.size();

        auto spatial_indices_builder =
          interpolation_points_spatial_indices_builder.init(i, n_spatial_indices);

        for (auto j : make_range(n_spatial_indices))
          {
            spatial_indices_builder.set(j, spatial_indices[j]);
          }
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
    libmesh_error_msg_if(rb_scm_eval.B_min.size() != rb_scm_eval.get_rb_theta_expansion().get_n_A_terms(),
                         "Size error while writing B_min");
    auto b_min_list = rb_scm_eval_builder.initBMin( rb_scm_eval.B_min.size() );
    for (auto i : index_range(rb_scm_eval.B_min))
      b_min_list.set(i, rb_scm_eval.get_B_min(i));
  }

  {
    libmesh_error_msg_if(rb_scm_eval.B_max.size() != rb_scm_eval.get_rb_theta_expansion().get_n_A_terms(),
                         "Size error while writing B_max");

    auto b_max_list = rb_scm_eval_builder.initBMax( rb_scm_eval.B_max.size() );
    for (auto i : index_range(rb_scm_eval.B_max))
      b_max_list.set(i, rb_scm_eval.get_B_max(i));
  }

  {
    auto cj_stability_vector =
      rb_scm_eval_builder.initCJStabilityVector( rb_scm_eval.C_J_stability_vector.size() );
    for (auto i : index_range(rb_scm_eval.C_J_stability_vector))
      cj_stability_vector.set(i, rb_scm_eval.get_C_J_stability_constraint(i));
  }

  {
    auto cj_parameters_outer =
      rb_scm_eval_builder.initCJ( rb_scm_eval.C_J.size() );

    for (auto i : index_range(rb_scm_eval.C_J))
      {
        auto cj_parameters_inner =
          cj_parameters_outer.init(i, rb_scm_eval.C_J[i].n_parameters());

        unsigned int count = 0;
        for (const auto & [key, val] : rb_scm_eval.C_J[i])
          {
            cj_parameters_inner[count].setName(key);
            cj_parameters_inner[count].setValue(val);
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

    for (unsigned int i=0; i<n_C_J_values; i++)
      for (unsigned int j=0; j<n_A_terms; j++)
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

// ---- Helper functions for adding data to capnp Builders (END) ----

} // namespace RBDataSerialization

} // namespace libMesh

#endif // #if defined(LIBMESH_HAVE_CAPNPROTO)
