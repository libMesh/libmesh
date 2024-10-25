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
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/elem.h"
#include "libmesh/mesh.h"
#include "libmesh/rb_data_deserialization.h"
#include "libmesh/transient_rb_theta_expansion.h"
#include "libmesh/transient_rb_evaluation.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_scm_evaluation.h"
#include "libmesh/rb_parametrized_function.h"
#include "libmesh/libmesh_logging.h"

// Cap'n'Proto includes
#include "capnp/serialize.h"

// C++ includes
#ifdef LIBMESH_HAVE_UNISTD_H
#include <unistd.h> // for close()
#endif
#include <iostream>
#include <fstream>
#include <fcntl.h>

namespace libMesh
{

namespace
{

/**
 * Helper function that reads either real or complex numbers, based on
 * the libMesh config options.
 */
template <typename T>
Number load_scalar_value(const T & value)
{
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
  return Number(value.getReal(), value.getImag());
#else
  return value;
#endif
}

}

namespace RBDataDeserialization
{

// ---- RBEvaluationDeserialization (BEGIN) ----

RBEvaluationDeserialization::RBEvaluationDeserialization(RBEvaluation & rb_eval)
  :
  _rb_eval(rb_eval)
{}

RBEvaluationDeserialization::~RBEvaluationDeserialization() = default;

void RBEvaluationDeserialization::read_from_file(const std::string & path,
                                                 bool read_error_bound_data)
{
  LOG_SCOPE("read_from_file()", "RBEvaluationDeserialization");

  int fd = open(path.c_str(), O_RDONLY);
  libmesh_error_msg_if(!fd, "Couldn't open the buffer file: " + path);

  // Turn off the limit to the amount of data we can read in
  capnp::ReaderOptions reader_options;
  reader_options.traversalLimitInWords = std::numeric_limits<uint64_t>::max();

  std::unique_ptr<capnp::StreamFdMessageReader> message;
  libmesh_try
    {
      message = std::make_unique<capnp::StreamFdMessageReader>(fd, reader_options);
    }
  libmesh_catch(...)
    {
      libmesh_error_msg("Failed to open capnp buffer");
    }

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  RBData::RBEvaluationReal::Reader rb_eval_reader =
    message->getRoot<RBData::RBEvaluationReal>();
#else
  RBData::RBEvaluationComplex::Reader rb_eval_reader =
    message->getRoot<RBData::RBEvaluationComplex>();
#endif

  load_rb_evaluation_data(_rb_eval, rb_eval_reader, read_error_bound_data);

  int error = close(fd);
  libmesh_error_msg_if(error, "Error closing a read-only file descriptor: " + path);
}

// ---- RBEvaluationDeserialization (END) ----


// ---- TransientRBEvaluationDeserialization (BEGIN) ----

TransientRBEvaluationDeserialization::
TransientRBEvaluationDeserialization(TransientRBEvaluation & trans_rb_eval) :
  _trans_rb_eval(trans_rb_eval)
{}

TransientRBEvaluationDeserialization::~TransientRBEvaluationDeserialization() = default;

void TransientRBEvaluationDeserialization::read_from_file(const std::string & path,
                                                          bool read_error_bound_data)
{
  LOG_SCOPE("read_from_file()", "TransientRBEvaluationDeserialization");

  int fd = open(path.c_str(), O_RDONLY);
  libmesh_error_msg_if(!fd, "Couldn't open the buffer file: " + path);

  // Turn off the limit to the amount of data we can read in
  capnp::ReaderOptions reader_options;
  reader_options.traversalLimitInWords = std::numeric_limits<uint64_t>::max();

  std::unique_ptr<capnp::StreamFdMessageReader> message;
  libmesh_try
    {
      message = std::make_unique<capnp::StreamFdMessageReader>(fd, reader_options);
    }
  libmesh_catch(...)
    {
      libmesh_error_msg("Failed to open capnp buffer");
    }

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  RBData::TransientRBEvaluationReal::Reader trans_rb_eval_reader =
    message->getRoot<RBData::TransientRBEvaluationReal>();
  RBData::RBEvaluationReal::Reader rb_eval_reader =
    trans_rb_eval_reader.getRbEvaluation();
#else
  RBData::TransientRBEvaluationComplex::Reader trans_rb_eval_reader =
    message->getRoot<RBData::TransientRBEvaluationComplex>();
  RBData::RBEvaluationComplex::Reader rb_eval_reader =
    trans_rb_eval_reader.getRbEvaluation();
#endif

  load_transient_rb_evaluation_data(_trans_rb_eval,
                                    rb_eval_reader,
                                    trans_rb_eval_reader,
                                    read_error_bound_data);

  int error = close(fd);
  libmesh_error_msg_if(error, "Error closing a read-only file descriptor: " + path);
}

// ---- TransientRBEvaluationDeserialization (END) ----


// ---- RBEIMEvaluationDeserialization (BEGIN) ----

RBEIMEvaluationDeserialization::
RBEIMEvaluationDeserialization(RBEIMEvaluation & rb_eim_eval) :
  _rb_eim_eval(rb_eim_eval)
{}

RBEIMEvaluationDeserialization::~RBEIMEvaluationDeserialization() = default;

void RBEIMEvaluationDeserialization::read_from_file(const std::string & path)
{
  LOG_SCOPE("read_from_file()", "RBEIMEvaluationDeserialization");

  int fd = open(path.c_str(), O_RDONLY);
  libmesh_error_msg_if(!fd, "Couldn't open the buffer file: " + path);

  // Turn off the limit to the amount of data we can read in
  capnp::ReaderOptions reader_options;
  reader_options.traversalLimitInWords = std::numeric_limits<uint64_t>::max();

  std::unique_ptr<capnp::StreamFdMessageReader> message;
  libmesh_try
    {
      message = std::make_unique<capnp::StreamFdMessageReader>(fd, reader_options);
    }
  libmesh_catch(...)
    {
      libmesh_error_msg("Failed to open capnp buffer");
    }

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  RBData::RBEIMEvaluationReal::Reader rb_eim_eval_reader =
    message->getRoot<RBData::RBEIMEvaluationReal>();
#else
  RBData::RBEIMEvaluationComplex::Reader rb_eim_eval_reader =
    message->getRoot<RBData::RBEIMEvaluationComplex>();
#endif

  load_rb_eim_evaluation_data(_rb_eim_eval,
                              rb_eim_eval_reader);

  int error = close(fd);
  libmesh_error_msg_if(error, "Error closing a read-only file descriptor: " + path);
}

// ---- RBEIMEvaluationDeserialization (END) ----


// ---- RBSCMEvaluationDeserialization (BEGIN) ----
// RBSCMEvaluationDeserialization is only available if both SLEPC and
// GLPK are available.

#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

RBSCMEvaluationDeserialization::
RBSCMEvaluationDeserialization(RBSCMEvaluation & rb_scm_eval) :
  _rb_scm_eval(rb_scm_eval)
{}

RBSCMEvaluationDeserialization::~RBSCMEvaluationDeserialization() = default;

void RBSCMEvaluationDeserialization::read_from_file(const std::string & path)
{
  LOG_SCOPE("read_from_file()", "RBSCMEvaluationDeserialization");

  int fd = open(path.c_str(), O_RDONLY);
  libmesh_error_msg_if(!fd, "Couldn't open the buffer file: " + path);

  // Turn off the limit to the amount of data we can read in
  capnp::ReaderOptions reader_options;
  reader_options.traversalLimitInWords = std::numeric_limits<uint64_t>::max();

  std::unique_ptr<capnp::StreamFdMessageReader> message;
  libmesh_try
    {
      message = std::make_unique<capnp::StreamFdMessageReader>(fd, reader_options);
    }
  libmesh_catch(...)
    {
      libmesh_error_msg("Failed to open capnp buffer");
    }

  RBData::RBSCMEvaluation::Reader rb_scm_eval_reader =
    message->getRoot<RBData::RBSCMEvaluation>();

  load_rb_scm_evaluation_data(_rb_scm_eval,
                              rb_scm_eval_reader);

  int error = close(fd);
  libmesh_error_msg_if(error, "Error closing a read-only file descriptor: " + path);
}

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

// ---- RBSCMEvaluationDeserialization (END) ----


// ---- Helper functions for loading data from buffers (BEGIN) ----

void load_parameter_ranges(RBParametrized & rb_evaluation,
                           RBData::ParameterRanges::Reader & parameter_ranges,
                           RBData::DiscreteParameterList::Reader & discrete_parameters_list)
{
  // Continuous parameters
  RBParameters parameters_min;
  RBParameters parameters_max;

  // Discrete parameters
  std::map<std::string, std::vector<Real>> discrete_parameter_values;

  // The RBData::ParameterRanges::Reader (CapnProto class) will throw
  // an exception if there is a problem with the data file (e.g. if it
  // doesn't exist). We don't use the libmesh_try/catch() macros here
  // since they don't really allow you to do any processing of the
  // exception object.
#ifdef LIBMESH_ENABLE_EXCEPTIONS
  try
#endif
  {
    const auto & parameter_names = parameter_ranges.getNames();
    for (auto i : make_range(parameter_names.size()))
      {
        parameters_min.set_value(parameter_names[i], parameter_ranges.getMinValues()[i]);
        parameters_max.set_value(parameter_names[i], parameter_ranges.getMaxValues()[i]);
      }


    const auto & discrete_names = discrete_parameters_list.getNames();
    for (auto i : make_range(discrete_names.size()))
      {
        const auto & value_list = discrete_parameters_list.getValues()[i];
        unsigned int n_values = value_list.size();
        std::vector<Real> values(n_values);
        for (unsigned int j=0; j<n_values; ++j)
          values[j] = value_list[j];

        discrete_parameter_values[discrete_names[i]] = values;
      }
  }
#ifdef LIBMESH_ENABLE_EXCEPTIONS
  catch (std::exception & e)
  {
    libmesh_error_msg("Error loading parameter ranges from capnp reader.\n"
                      "This usually means that the training data either doesn't exist or is out of date.\n"
                      "Detailed information about the error is below:\n\n" << e.what() << "\n");
  }
#endif

  rb_evaluation.initialize_parameters(parameters_min,
                                      parameters_max,
                                      discrete_parameter_values);
}

template <typename RBEvaluationReaderNumber>
void load_rb_evaluation_data(RBEvaluation & rb_evaluation,
                             RBEvaluationReaderNumber & rb_evaluation_reader,
                             bool read_error_bound_data)
{
  // Set number of basis functions
  unsigned int n_bfs = rb_evaluation_reader.getNBfs();
  rb_evaluation.set_n_basis_functions(n_bfs);

  rb_evaluation.resize_data_structures(n_bfs, read_error_bound_data);

  auto parameter_ranges =
    rb_evaluation_reader.getParameterRanges();
  auto discrete_parameters_list =
    rb_evaluation_reader.getDiscreteParameters();

  load_parameter_ranges(rb_evaluation,
                        parameter_ranges,
                        discrete_parameters_list);

  const RBThetaExpansion & rb_theta_expansion = rb_evaluation.get_rb_theta_expansion();

  unsigned int n_F_terms = rb_theta_expansion.get_n_F_terms();
  unsigned int n_A_terms = rb_theta_expansion.get_n_A_terms();

  if (read_error_bound_data)
    {

      // Fq representor inner-product data
      {
        unsigned int Q_f_hat = n_F_terms*(n_F_terms+1)/2;

        auto fq_innerprods_list = rb_evaluation_reader.getFqInnerprods();
        libmesh_error_msg_if(fq_innerprods_list.size() != Q_f_hat,
                             "Size error while reading Fq representor norm data from buffer.");

        for (unsigned int i=0; i < Q_f_hat; ++i)
          rb_evaluation.Fq_representor_innerprods[i] = load_scalar_value(fq_innerprods_list[i]);
      }

      // Fq_Aq representor inner-product data
      {
        auto fq_aq_innerprods_list = rb_evaluation_reader.getFqAqInnerprods();
        libmesh_error_msg_if(fq_aq_innerprods_list.size() != n_F_terms*n_A_terms*n_bfs,
                             "Size error while reading Fq-Aq representor norm data from buffer.");

        for (unsigned int q_f=0; q_f<n_F_terms; ++q_f)
          for (unsigned int q_a=0; q_a<n_A_terms; ++q_a)
            for (unsigned int i=0; i<n_bfs; ++i)
              {
                unsigned int offset = q_f*n_A_terms*n_bfs + q_a*n_bfs + i;
                rb_evaluation.Fq_Aq_representor_innerprods[q_f][q_a][i] =
                  load_scalar_value(fq_aq_innerprods_list[offset]);
              }
      }

      // Aq_Aq representor inner-product data
      {
        unsigned int Q_a_hat = n_A_terms*(n_A_terms+1)/2;
        auto aq_aq_innerprods_list = rb_evaluation_reader.getAqAqInnerprods();
        libmesh_error_msg_if(aq_aq_innerprods_list.size() != Q_a_hat*n_bfs*n_bfs,
                             "Size error while reading Aq-Aq representor norm data from buffer.");

        for (unsigned int i=0; i<Q_a_hat; ++i)
          for (unsigned int j=0; j<n_bfs; ++j)
            for (unsigned int l=0; l<n_bfs; ++l)
              {
                unsigned int offset = i*n_bfs*n_bfs + j*n_bfs + l;
                rb_evaluation.Aq_Aq_representor_innerprods[i][j][l] =
                  load_scalar_value(aq_aq_innerprods_list[offset]);
              }
      }

      // Output dual inner-product data
      {
        unsigned int n_outputs = rb_theta_expansion.get_n_outputs();
        auto output_innerprod_outer = rb_evaluation_reader.getOutputDualInnerprods();

        libmesh_error_msg_if(output_innerprod_outer.size() != n_outputs,
                             "Incorrect number of outputs detected in the buffer");

        for (unsigned int output_id=0; output_id<n_outputs; ++output_id)
          {
            unsigned int n_output_terms = rb_theta_expansion.get_n_output_terms(output_id);

            unsigned int Q_l_hat = n_output_terms*(n_output_terms+1)/2;
            auto output_innerprod_inner = output_innerprod_outer[output_id];

            libmesh_error_msg_if(output_innerprod_inner.size() != Q_l_hat,
                                 "Incorrect number of output terms detected in the buffer");

            for (unsigned int q=0; q<Q_l_hat; ++q)
              {
                rb_evaluation.output_dual_innerprods[output_id][q] =
                  load_scalar_value(output_innerprod_inner[q]);
              }
          }
      }
    }

  // Output vectors
  {
    unsigned int n_outputs = rb_theta_expansion.get_n_outputs();
    auto output_vector_outer = rb_evaluation_reader.getOutputVectors();

    libmesh_error_msg_if(output_vector_outer.size() != n_outputs,
                         "Incorrect number of outputs detected in the buffer");

    for (unsigned int output_id=0; output_id<n_outputs; ++output_id)
      {
        unsigned int n_output_terms = rb_theta_expansion.get_n_output_terms(output_id);

        auto output_vector_middle = output_vector_outer[output_id];
        libmesh_error_msg_if(output_vector_middle.size() != n_output_terms,
                             "Incorrect number of output terms detected in the buffer");

        for (unsigned int q_l=0; q_l<n_output_terms; ++q_l)
          {
            auto output_vectors_inner_list = output_vector_middle[q_l];

            libmesh_error_msg_if(output_vectors_inner_list.size() != n_bfs,
                                 "Incorrect number of output terms detected in the buffer");

            for (unsigned int j=0; j<n_bfs; ++j)
              {
                rb_evaluation.RB_output_vectors[output_id][q_l](j) =
                  load_scalar_value(output_vectors_inner_list[j]);
              }
          }
      }
  }

  // Fq vectors and Aq matrices
  {
    auto rb_fq_vectors_outer_list = rb_evaluation_reader.getRbFqVectors();
    libmesh_error_msg_if(rb_fq_vectors_outer_list.size() != n_F_terms,
                         "Incorrect number of Fq vectors detected in the buffer");

    for (unsigned int q_f=0; q_f<n_F_terms; ++q_f)
      {
        auto rb_fq_vectors_inner_list = rb_fq_vectors_outer_list[q_f];
        libmesh_error_msg_if(rb_fq_vectors_inner_list.size() != n_bfs,
                             "Incorrect Fq vector size detected in the buffer");

        for (unsigned int i=0; i < n_bfs; ++i)
          {
            rb_evaluation.RB_Fq_vector[q_f](i) =
              load_scalar_value(rb_fq_vectors_inner_list[i]);
          }
      }

    auto rb_Aq_matrices_outer_list = rb_evaluation_reader.getRbAqMatrices();
    libmesh_error_msg_if(rb_Aq_matrices_outer_list.size() != n_A_terms,
                         "Incorrect number of Aq matrices detected in the buffer");

    for (unsigned int q_a=0; q_a<n_A_terms; ++q_a)
      {
        auto rb_Aq_matrices_inner_list = rb_Aq_matrices_outer_list[q_a];
        libmesh_error_msg_if(rb_Aq_matrices_inner_list.size() != n_bfs*n_bfs,
                             "Incorrect Aq matrix size detected in the buffer");

        for (unsigned int i=0; i<n_bfs; ++i)
          for (unsigned int j=0; j<n_bfs; ++j)
            {
              unsigned int offset = i*n_bfs+j;
              rb_evaluation.RB_Aq_vector[q_a](i,j) =
                load_scalar_value(rb_Aq_matrices_inner_list[offset]);
            }
      }
  }

  // Inner-product matrix
  if (rb_evaluation.compute_RB_inner_product)
    {
      auto rb_inner_product_matrix_list =
        rb_evaluation_reader.getRbInnerProductMatrix();

      libmesh_error_msg_if(rb_inner_product_matrix_list.size() != n_bfs*n_bfs,
                           "Size error while reading the inner product matrix.");

      for (unsigned int i=0; i<n_bfs; ++i)
        for (unsigned int j=0; j<n_bfs; ++j)
          {
            unsigned int offset = i*n_bfs + j;
            rb_evaluation.RB_inner_product_matrix(i,j) =
              load_scalar_value(rb_inner_product_matrix_list[offset]);
          }
    }
}

template <typename RBEvaluationReaderNumber, typename TransRBEvaluationReaderNumber>
void load_transient_rb_evaluation_data(TransientRBEvaluation & trans_rb_eval,
                                       RBEvaluationReaderNumber & rb_eval_reader,
                                       TransRBEvaluationReaderNumber & trans_rb_eval_reader,
                                       bool read_error_bound_data)
{
  load_rb_evaluation_data(trans_rb_eval,
                          rb_eval_reader,
                          read_error_bound_data);

  trans_rb_eval.set_delta_t( trans_rb_eval_reader.getDeltaT() );
  trans_rb_eval.set_euler_theta( trans_rb_eval_reader.getEulerTheta() );
  trans_rb_eval.set_n_time_steps( trans_rb_eval_reader.getNTimeSteps() );
  trans_rb_eval.set_time_step( trans_rb_eval_reader.getTimeStep() );

  unsigned int n_bfs = rb_eval_reader.getNBfs();
  unsigned int n_F_terms = trans_rb_eval.get_rb_theta_expansion().get_n_F_terms();
  unsigned int n_A_terms = trans_rb_eval.get_rb_theta_expansion().get_n_A_terms();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(trans_rb_eval.get_rb_theta_expansion());
  unsigned int n_M_terms = trans_theta_expansion.get_n_M_terms();

  // L2 matrix
  {
    auto rb_L2_matrix_list =
      trans_rb_eval_reader.getRbL2Matrix();

    libmesh_error_msg_if(rb_L2_matrix_list.size() != n_bfs*n_bfs,
                         "Size error while reading the L2 matrix.");

    for (unsigned int i=0; i<n_bfs; ++i)
      for (unsigned int j=0; j<n_bfs; ++j)
        {
          unsigned int offset = i*n_bfs + j;
          trans_rb_eval.RB_L2_matrix(i,j) =
            load_scalar_value(rb_L2_matrix_list[offset]);
        }
  }

  // Mq matrices
  {
    auto rb_Mq_matrices_outer_list = trans_rb_eval_reader.getRbMqMatrices();

    libmesh_error_msg_if(rb_Mq_matrices_outer_list.size() != n_M_terms,
                         "Incorrect number of Mq matrices detected in the buffer");

    for (unsigned int q_m=0; q_m < n_M_terms; ++q_m)
      {
        auto rb_Mq_matrices_inner_list = rb_Mq_matrices_outer_list[q_m];
        libmesh_error_msg_if(rb_Mq_matrices_inner_list.size() != n_bfs*n_bfs,
                             "Incorrect Mq matrix size detected in the buffer");

        for (unsigned int i=0; i<n_bfs; ++i)
          for (unsigned int j=0; j<n_bfs; ++j)
            {
              unsigned int offset = i*n_bfs+j;
              trans_rb_eval.RB_M_q_vector[q_m](i,j) =
                load_scalar_value(rb_Mq_matrices_inner_list[offset]);
            }
      }
  }

  // The initial condition and L2 error at t=0.
  {
    auto initial_l2_errors_reader =
      trans_rb_eval_reader.getInitialL2Errors();
    libmesh_error_msg_if(initial_l2_errors_reader.size() != n_bfs,
                         "Incorrect number of initial L2 error terms detected in the buffer");

    auto initial_conditions_outer_list =
      trans_rb_eval_reader.getInitialConditions();
    libmesh_error_msg_if(initial_conditions_outer_list.size() != n_bfs,
                         "Incorrect number of outer initial conditions detected in the buffer");

    for (unsigned int i=0; i<n_bfs; i++)
      {
        trans_rb_eval.initial_L2_error_all_N[i] =
          initial_l2_errors_reader[i];

        auto initial_conditions_inner_list = initial_conditions_outer_list[i];
        libmesh_error_msg_if(initial_conditions_inner_list.size() != (i+1),
                             "Incorrect number of inner initial conditions detected in the buffer");

        for (unsigned int j=0; j<=i; j++)
          {
            trans_rb_eval.RB_initial_condition_all_N[i](j) =
              load_scalar_value(initial_conditions_inner_list[j]);
          }
      }
  }


  if (read_error_bound_data)
    {
      // Fq_Mq data
      {
        auto fq_mq_innerprods_list = trans_rb_eval_reader.getFqMqInnerprods();
        libmesh_error_msg_if(fq_mq_innerprods_list.size() != n_F_terms*n_M_terms*n_bfs,
                             "Size error while reading Fq-Mq representor data from buffer.");

        for (unsigned int q_f=0; q_f<n_F_terms; ++q_f)
          for (unsigned int q_m=0; q_m<n_M_terms; ++q_m)
            for (unsigned int i=0; i<n_bfs; ++i)
              {
                unsigned int offset = q_f*n_M_terms*n_bfs + q_m*n_bfs + i;
                trans_rb_eval.Fq_Mq_representor_innerprods[q_f][q_m][i] =
                  load_scalar_value(fq_mq_innerprods_list[offset]);
              }
      }


      // Mq_Mq representor inner-product data
      {
        unsigned int Q_m_hat = n_M_terms*(n_M_terms+1)/2;
        auto mq_mq_innerprods_list = trans_rb_eval_reader.getMqMqInnerprods();
        libmesh_error_msg_if(mq_mq_innerprods_list.size() != Q_m_hat*n_bfs*n_bfs,
                             "Size error while reading Mq-Mq representor data from buffer.");

        for (unsigned int i=0; i<Q_m_hat; ++i)
          for (unsigned int j=0; j<n_bfs; ++j)
            for (unsigned int l=0; l<n_bfs; ++l)
              {
                unsigned int offset = i*n_bfs*n_bfs + j*n_bfs + l;
                trans_rb_eval.Mq_Mq_representor_innerprods[i][j][l] =
                  load_scalar_value(mq_mq_innerprods_list[offset]);
              }
      }

      // Aq_Mq representor inner-product data
      {
        auto aq_mq_innerprods_list =
          trans_rb_eval_reader.getAqMqInnerprods();
        libmesh_error_msg_if(aq_mq_innerprods_list.size() != n_A_terms*n_M_terms*n_bfs*n_bfs,
                             "Size error while reading Aq-Mq representor data from buffer.");

        for (unsigned int q_a=0; q_a<n_A_terms; q_a++)
          for (unsigned int q_m=0; q_m<n_M_terms; q_m++)
            for (unsigned int i=0; i<n_bfs; i++)
              for (unsigned int j=0; j<n_bfs; j++)
                {
                  unsigned int offset =
                    q_a*(n_M_terms*n_bfs*n_bfs) + q_m*(n_bfs*n_bfs) + i*n_bfs + j;

                  trans_rb_eval.Aq_Mq_representor_innerprods[q_a][q_m][i][j] =
                    load_scalar_value(aq_mq_innerprods_list[offset]);
                }
      }

    }
}

template <typename RBEIMEvaluationReaderNumber>
void load_rb_eim_evaluation_data(RBEIMEvaluation & rb_eim_evaluation,
                                 RBEIMEvaluationReaderNumber & rb_eim_evaluation_reader)
{
  // Set number of basis functions
  unsigned int n_bfs = rb_eim_evaluation_reader.getNBfs();
  rb_eim_evaluation.set_n_basis_functions(n_bfs);

  rb_eim_evaluation.resize_data_structures(n_bfs);

  // EIM interpolation matrix
  {
    auto interpolation_matrix_list = rb_eim_evaluation_reader.getInterpolationMatrix();

    libmesh_error_msg_if(interpolation_matrix_list.size() != n_bfs*(n_bfs+1)/2,
                         "Size error while reading the eim inner product matrix.");

    for (unsigned int i=0; i<n_bfs; ++i)
      for (unsigned int j=0; j<=i; ++j)
        {
          unsigned int offset = i*(i+1)/2 + j;
          rb_eim_evaluation.set_interpolation_matrix_entry(
            i,j, load_scalar_value(interpolation_matrix_list[offset]));
        }
  }

  // Note that we do not check rb_eim_evaluation.use_eim_error_indicator()
  // below because it can be false in some cases when we're reading in
  // data (e.g. if we're using a base class RBEIMEvaluation to do the
  // reading). The check below will still give the right behavior
  // regardless of the value of use_eim_error_indicator().
  if (rb_eim_evaluation_reader.getInterpolationXyz().size() > n_bfs)
  {
    auto error_indicator_data_list = rb_eim_evaluation_reader.getEimErrorIndicatorInterpData();

    libmesh_error_msg_if(error_indicator_data_list.size() != n_bfs,
                         "Size error while reading the eim error indicator data.");

    DenseVector<Number> eim_error_indicator_data(n_bfs);
    for (unsigned int i=0; i<n_bfs; ++i)
      eim_error_indicator_data(i) = load_scalar_value(error_indicator_data_list[i]);

    rb_eim_evaluation.set_error_indicator_interpolation_row(eim_error_indicator_data);
  }

  // If we're using the EIM error indicator then we store one extra
  // interpolation point and associated data, hence we increment n_bfs
  // here so that we write out the extra data below.
  //
  // However the interpolation matrix and _error_indicator_interpolation_row
  // does not include this extra point, which is why we read it in above
  // before n_bfs is incremented.
  if (rb_eim_evaluation_reader.getInterpolationXyz().size() > n_bfs)
    n_bfs++;

  auto parameter_ranges =
    rb_eim_evaluation_reader.getParameterRanges();
  auto discrete_parameters_list =
    rb_eim_evaluation_reader.getDiscreteParameters();

  load_parameter_ranges(rb_eim_evaluation,
                        parameter_ranges,
                        discrete_parameters_list);

  // Interpolation points
  {
    auto interpolation_points_list =
      rb_eim_evaluation_reader.getInterpolationXyz();

    libmesh_error_msg_if(interpolation_points_list.size() != n_bfs,
                         "Size error while reading the eim interpolation points.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        Point p;
        load_point(interpolation_points_list[i], p);
        rb_eim_evaluation.add_interpolation_points_xyz(p);
      }
  }

  // Interpolation points components
  {
    auto interpolation_points_comp_list =
      rb_eim_evaluation_reader.getInterpolationComp();

    libmesh_error_msg_if(interpolation_points_comp_list.size() != n_bfs,
                         "Size error while reading the eim interpolation components.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        rb_eim_evaluation.add_interpolation_points_comp(interpolation_points_comp_list[i]);
      }
  }

  // Interpolation points subdomain IDs
  {
    auto interpolation_points_subdomain_id_list =
      rb_eim_evaluation_reader.getInterpolationSubdomainId();

    libmesh_error_msg_if(interpolation_points_subdomain_id_list.size() != n_bfs,
                         "Size error while reading the eim interpolation subdomain IDs.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        rb_eim_evaluation.add_interpolation_points_subdomain_id(interpolation_points_subdomain_id_list[i]);
      }
  }

  // Interpolation points side indices, relevant if the parametrized function is defined
  // on mesh sides or nodesets
  if (rb_eim_evaluation.get_parametrized_function().on_mesh_sides() ||
      rb_eim_evaluation.get_parametrized_function().on_mesh_nodes())
  {
    auto interpolation_points_boundary_id_list =
      rb_eim_evaluation_reader.getInterpolationBoundaryId();

    libmesh_error_msg_if(interpolation_points_boundary_id_list.size() != n_bfs,
                         "Size error while reading the eim interpolation boundary IDs.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        rb_eim_evaluation.add_interpolation_points_boundary_id(interpolation_points_boundary_id_list[i]);
      }
  }

  // Interpolation points element IDs
  {
    auto interpolation_points_elem_id_list =
      rb_eim_evaluation_reader.getInterpolationElemId();

    libmesh_error_msg_if(interpolation_points_elem_id_list.size() != n_bfs,
                         "Size error while reading the eim interpolation element IDs.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        rb_eim_evaluation.add_interpolation_points_elem_id(interpolation_points_elem_id_list[i]);
      }
  }

  // Interpolation points node IDs, relevant if the parametrized function is defined on mesh sides
  if (rb_eim_evaluation.get_parametrized_function().on_mesh_nodes())
  {
    auto interpolation_points_node_id_list =
      rb_eim_evaluation_reader.getInterpolationNodeId();

    libmesh_error_msg_if(interpolation_points_node_id_list.size() != n_bfs,
                         "Size error while reading the eim interpolation node IDs.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        rb_eim_evaluation.add_interpolation_points_node_id(interpolation_points_node_id_list[i]);
      }
  }

  // Interpolation points side indices, relevant if the parametrized function is defined on mesh sides
  if (rb_eim_evaluation.get_parametrized_function().on_mesh_sides())
  {
    auto interpolation_points_side_index_list =
      rb_eim_evaluation_reader.getInterpolationSideIndex();

    libmesh_error_msg_if(interpolation_points_side_index_list.size() != n_bfs,
                         "Size error while reading the eim interpolation side indices.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        rb_eim_evaluation.add_interpolation_points_side_index(interpolation_points_side_index_list[i]);
      }
  }

  // Interpolation points quad point indices
  {
    auto interpolation_points_qp_list =
      rb_eim_evaluation_reader.getInterpolationQp();

    libmesh_error_msg_if(interpolation_points_qp_list.size() != n_bfs,
                         "Size error while reading the eim interpolation qps.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        rb_eim_evaluation.add_interpolation_points_qp(interpolation_points_qp_list[i]);
      }
  }

  unsigned int n_elems = rb_eim_evaluation_reader.getNElems();
  // Elem id to local index map
  {
    auto elem_id_to_local_index_list =
      rb_eim_evaluation_reader.getElemIdToLocalIndex();

    if (elem_id_to_local_index_list.size() > 0)
      {
        libmesh_error_msg_if(elem_id_to_local_index_list.size() != n_elems,
                             "Size error while reading the eim elem id to local index map.");

        for (unsigned int i=0; i<n_elems; ++i)
          {
            rb_eim_evaluation.add_elem_id_local_index_map_entry(elem_id_to_local_index_list[i].getFirst(), elem_id_to_local_index_list[i].getSecond());
          }
      }
  }

  // Interpolation points JxW all qp values
  {
    auto interpolation_points_list_outer =
      rb_eim_evaluation_reader.getInterpolationJxWAllQp();

    if (interpolation_points_list_outer.size() > 0)
      {
        libmesh_error_msg_if(interpolation_points_list_outer.size() != n_elems,
                             "Size error while reading the eim JxW values.");

        for (unsigned int i=0; i<n_elems; ++i)
          {
            auto interpolation_points_list_inner = interpolation_points_list_outer[i];

            std::vector<Real> JxW(interpolation_points_list_inner.size());
            for (unsigned int j=0; j<JxW.size(); j++)
              {
                JxW[j] = interpolation_points_list_inner[j];
              }
            rb_eim_evaluation.add_interpolation_points_JxW_all_qp(JxW);
          }
      }
  }

  // Interpolation points phi_i all qp values
  {
    auto interpolation_points_list_outer =
      rb_eim_evaluation_reader.getInterpolationPhiValuesAllQp();

    if (interpolation_points_list_outer.size() > 0)
      {
        libmesh_error_msg_if(interpolation_points_list_outer.size() != n_elems,
                             "Size error while reading the eim phi_i all qp values.");

        for (unsigned int i=0; i<n_elems; ++i)
          {
            auto interpolation_points_list_middle = interpolation_points_list_outer[i];

            std::vector<std::vector<Real>> phi_i_all_qp(interpolation_points_list_middle.size());
            for (auto j : index_range(phi_i_all_qp))
              {
                auto interpolation_points_list_inner = interpolation_points_list_middle[j];

                phi_i_all_qp[j].resize(interpolation_points_list_inner.size());
                for (auto k : index_range(phi_i_all_qp[j]))
                  phi_i_all_qp[j][k] = interpolation_points_list_inner[k];
              }
            rb_eim_evaluation.add_interpolation_points_phi_i_all_qp(phi_i_all_qp);
          }
      }
  }

  // Dxyzdxi at the element center for the element that contains each interpolation point.
  {
    auto elem_center_dxyzdxi =
      rb_eim_evaluation_reader.getInterpolationDxyzDxiElem();

    if (elem_center_dxyzdxi.size() > 0)
      {
        libmesh_error_msg_if(elem_center_dxyzdxi.size() != n_elems,
                             "Size error while reading the eim elem center tangent derivative dxyzdxi.");

        Point dxyzdxi_buffer;
        for (const auto i : make_range(n_elems))
          {
            load_point(elem_center_dxyzdxi[i], dxyzdxi_buffer);
            rb_eim_evaluation.add_elem_center_dxyzdxi(dxyzdxi_buffer);
          }
      }
  }

  // Dxyzdxi at the element center for the element that contains each interpolation point.
  {
    auto elem_center_dxyzdeta =
      rb_eim_evaluation_reader.getInterpolationDxyzDetaElem();

    if (elem_center_dxyzdeta.size() > 0)
      {
        libmesh_error_msg_if(elem_center_dxyzdeta.size() != n_elems,
                             "Size error while reading the eim elem center tangent derivative dxyzdeta.");

        Point dxyzdeta_buffer;
        for (const auto i : make_range(n_elems))
          {
            load_point(elem_center_dxyzdeta[i], dxyzdeta_buffer);
            rb_eim_evaluation.add_elem_center_dxyzdeta(dxyzdeta_buffer);
          }
      }
  }

  // Quadrature rule order associated to the element that contains each interpolation point.
  {
    auto interpolation_points_qrule_order_list =
      rb_eim_evaluation_reader.getInterpolationQruleOrder();

    if (interpolation_points_qrule_order_list.size() > 0)
      {
        libmesh_error_msg_if(interpolation_points_qrule_order_list.size() != n_elems,
                             "Size error while reading the eim elem qrule order.");

        for (unsigned int i=0; i<n_elems; ++i)
          {
            rb_eim_evaluation.add_interpolation_points_qrule_order(static_cast<Order>(interpolation_points_qrule_order_list[i]));
          }
      }
  }

  // Interpolation points element types
  {
    auto interpolation_points_elem_type_list =
      rb_eim_evaluation_reader.getInterpolationElemType();

    if (interpolation_points_elem_type_list.size() > 0)
      {
        libmesh_error_msg_if(interpolation_points_elem_type_list.size() != n_bfs,
                             "Size error while reading the eim interpolation element types.");

        for (unsigned int i=0; i<n_bfs; ++i)
          {
            rb_eim_evaluation.add_interpolation_points_elem_type(static_cast<ElemType>(interpolation_points_elem_type_list[i]));
          }
      }
  }

  // Interpolation points perturbations
  {
    auto interpolation_points_list_outer =
      rb_eim_evaluation_reader.getInterpolationXyzPerturb();

    libmesh_error_msg_if(interpolation_points_list_outer.size() != n_bfs,
                         "Size error while reading the eim interpolation points.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        auto interpolation_points_list_inner = interpolation_points_list_outer[i];

        std::vector<Point> perturbs(interpolation_points_list_inner.size());
        for (unsigned int j=0; j<perturbs.size(); j++)
          {
            load_point(interpolation_points_list_inner[j], perturbs[j]);
          }
        rb_eim_evaluation.add_interpolation_points_xyz_perturbations(perturbs);
      }
  }

  // Optionally load EIM solutions for the training set
  if (rb_eim_evaluation.get_parametrized_function().is_lookup_table)
    {
      auto eim_rhs_list_outer =
        rb_eim_evaluation_reader.getEimSolutionsForTrainingSet();

      std::vector<DenseVector<Number>> & eim_solutions = rb_eim_evaluation.get_eim_solutions_for_training_set();
      eim_solutions.clear();
      eim_solutions.resize(eim_rhs_list_outer.size());

      for (auto i : make_range(eim_rhs_list_outer.size()))
        {
          auto eim_rhs_list_inner = eim_rhs_list_outer[i];

          DenseVector<Number> values(eim_rhs_list_inner.size());
          for (auto j : index_range(values))
            {
              values(j) = load_scalar_value(eim_rhs_list_inner[j]);
            }
          eim_solutions[i] = values;
        }
    }

  // The shape function values at the interpolation points. This can be used to evaluate nodal data
  // at EIM interpolation points, which are at quadrature points.
  {
    auto interpolation_points_list_outer =
      rb_eim_evaluation_reader.getInterpolationPhiValues();

    libmesh_error_msg_if(interpolation_points_list_outer.size() != n_bfs,
                         "Size error while reading the eim interpolation points.");

    for (unsigned int i=0; i<n_bfs; ++i)
      {
        auto interpolation_points_list_inner = interpolation_points_list_outer[i];

        std::vector<Real> phi_i_qp(interpolation_points_list_inner.size());
        for (unsigned int j=0; j<phi_i_qp.size(); j++)
          {
            phi_i_qp[j] = interpolation_points_list_inner[j];
          }
        rb_eim_evaluation.add_interpolation_points_phi_i_qp(phi_i_qp);
      }
  }

  // Read in the spatial indices at interpolation points. This data is optional
  // so this may be empty.
  {
    auto interpolation_points_list_outer =
      rb_eim_evaluation_reader.getInterpolationSpatialIndices();

    for (unsigned int i=0; i<interpolation_points_list_outer.size(); ++i)
      {
        auto interpolation_points_list_inner = interpolation_points_list_outer[i];

        std::vector<unsigned int> spatial_indices(interpolation_points_list_inner.size());
        for (unsigned int j=0; j<spatial_indices.size(); j++)
          {
            spatial_indices[j] = interpolation_points_list_inner[j];
          }
        rb_eim_evaluation.add_interpolation_points_spatial_indices(spatial_indices);
      }

    rb_eim_evaluation.initialize_param_fn_spatial_indices();
  }
}

#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)
void load_rb_scm_evaluation_data(RBSCMEvaluation & rb_scm_eval,
                                 RBData::RBSCMEvaluation::Reader & rb_scm_evaluation_reader)
{
  auto parameter_ranges =
    rb_scm_evaluation_reader.getParameterRanges();
  auto discrete_parameters_list =
    rb_scm_evaluation_reader.getDiscreteParameters();
  load_parameter_ranges(rb_scm_eval,
                        parameter_ranges,
                        discrete_parameters_list);

  unsigned int n_A_terms = rb_scm_eval.get_rb_theta_expansion().get_n_A_terms();

  {
    auto b_min_list = rb_scm_evaluation_reader.getBMin();

    libmesh_error_msg_if(b_min_list.size() != n_A_terms,
                         "Size error while reading B_min");

    rb_scm_eval.B_min.clear();
    for (unsigned int i=0; i<n_A_terms; i++)
      rb_scm_eval.B_min.push_back(b_min_list[i]);
  }

  {
    auto b_max_list = rb_scm_evaluation_reader.getBMax();

    libmesh_error_msg_if(b_max_list.size() != n_A_terms,
                         "Size error while reading B_max");

    rb_scm_eval.B_max.clear();
    for (unsigned int i=0; i<n_A_terms; i++)
      rb_scm_eval.B_max.push_back(b_max_list[i]);
  }

  {
    auto cJ_stability_vector =
      rb_scm_evaluation_reader.getCJStabilityVector();

    rb_scm_eval.C_J_stability_vector.clear();
    for (const auto & val : cJ_stability_vector)
      rb_scm_eval.C_J_stability_vector.push_back(val);
  }

  {
    auto cJ_parameters_outer =
      rb_scm_evaluation_reader.getCJ();

    rb_scm_eval.C_J.resize(cJ_parameters_outer.size());
    for (auto i : make_range(cJ_parameters_outer.size()))
      {
        auto cJ_parameters_inner =
          cJ_parameters_outer[i];

        for (auto j : make_range(cJ_parameters_inner.size()))
          {
            std::string param_name = cJ_parameters_inner[j].getName();
            Real param_value = cJ_parameters_inner[j].getValue();
            rb_scm_eval.C_J[i].set_value(param_name, param_value);
          }
      }
  }

  {
    auto scm_ub_vectors =
      rb_scm_evaluation_reader.getScmUbVectors();

    // The number of UB vectors is the same as the number of C_J values
    unsigned int n_C_J_values = rb_scm_evaluation_reader.getCJ().size();

    libmesh_error_msg_if(scm_ub_vectors.size() != n_C_J_values*n_A_terms,
                         "Size mismatch in SCB UB vectors");

    rb_scm_eval.SCM_UB_vectors.resize( n_C_J_values );
    for (unsigned int i=0; i<n_C_J_values; i++)
      {
        rb_scm_eval.SCM_UB_vectors[i].resize(n_A_terms);
        for (unsigned int j=0; j<n_A_terms; j++)
          {
            unsigned int offset = i*n_A_terms + j;
            rb_scm_eval.SCM_UB_vectors[i][j] = scm_ub_vectors[offset];

          }
      }
  }

}
#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK



void load_point(RBData::Point3D::Reader point_reader, Point & point)
{
  point(0) = point_reader.getX();

  if (LIBMESH_DIM >= 2)
    point(1) = point_reader.getY();

  if (LIBMESH_DIM >= 3)
    point(2) = point_reader.getZ();
}

// ---- Helper functions for adding data to capnp Builders (END) ----

} // namespace RBDataSerialization

} // namespace libMesh

#endif // #if defined(LIBMESH_HAVE_CAPNPROTO)
