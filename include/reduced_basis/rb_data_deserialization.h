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

#ifndef RB_DATA_DESERIALIZATION_H
#define RB_DATA_DESERIALIZATION_H

// This class is only available if we have Cap'n Proto
#include "libmesh/libmesh_config.h"
#if defined(LIBMESH_HAVE_CAPNPROTO)

// libMesh/reduced_basis includes
#include "libmesh/rb_data.capnp.h"

// Cap'n'Proto includes
#include "capnp/message.h"

// C++ includes
#include <string>

namespace libMesh
{

// Forward declarations
class RBEvaluation;
class TransientRBEvaluation;
class RBEIMEvaluation;
class RBSCMEvaluation;
class RBParametrized;
class ReplicatedMesh;
class Elem;
class Point;

namespace RBDataDeserialization
{

/**
 * This class de-serializes an RBEvaluation object
 * using the Cap'n Proto library.
 *
 * \author David Knezevic
 * \date 2015
 * \brief Deserializes RBEvaluation objects using Cap'n Proto.
 */
class RBEvaluationDeserialization
{
public:

  /**
   * Initialize a new buffer using the structure from the Cap'n'Proto schema
   * described in rb_data.capnp.
   */
  RBEvaluationDeserialization(RBEvaluation & rb_eval);

  /**
   * Destructor.
   */
  virtual ~RBEvaluationDeserialization();

  /**
   * Write the Cap'n'Proto buffer to disk.
   */
  void read_from_file(const std::string & path, bool read_error_bound_data);

private:

  /**
   * The RBEvaluation object that we will read into.
   */
  RBEvaluation & _rb_eval;

};

/**
 * This class de-serializes a TransientRBEvaluation object
 * using the Cap'n Proto library.
 */
class TransientRBEvaluationDeserialization
{
public:

  /**
   * Initialize a new buffer using the structure from the Cap'n'Proto schema
   * described in rb_data.capnp.
   */
  TransientRBEvaluationDeserialization(TransientRBEvaluation & trans_rb_eval);

  /**
   * Destructor.
   */
  virtual ~TransientRBEvaluationDeserialization();

  /**
   * Write the Cap'n'Proto buffer to disk.
   */
  void read_from_file(const std::string & path, bool read_error_bound_data);

private:

  /**
   * The TransientRBEvaluation object that we will read into.
   */
  TransientRBEvaluation & _trans_rb_eval;

};

/**
 * This class de-serializes a RBEIMEvaluation object
 * using the Cap'n Proto library.
 */
class RBEIMEvaluationDeserialization
{
public:

  /**
   * Initialize a new buffer using the structure from the Cap'n'Proto schema
   * described in rb_data.capnp.
   */
  RBEIMEvaluationDeserialization(RBEIMEvaluation & trans_rb_eval);

  /**
   * Destructor.
   */
  virtual ~RBEIMEvaluationDeserialization();

  /**
   * Write the Cap'n'Proto buffer to disk.
   */
  void read_from_file(const std::string & path);

private:

  /**
   * The RBEIMEvaluation object we will read into.
   */
  RBEIMEvaluation & _rb_eim_eval;

};

// RBSCMEvaluation should only be available
// if SLEPc and GLPK support is enabled.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

/**
 * This class de-serializes a RBSCMEvaluation object
 * using the Cap'n Proto library.
 */
class RBSCMEvaluationDeserialization
{
public:

  /**
   * Initialize a new buffer using the structure from the Cap'n'Proto schema
   * described in rb_data.capnp.
   */
  RBSCMEvaluationDeserialization(RBSCMEvaluation & trans_rb_eval);

  /**
   * Destructor.
   */
  virtual ~RBSCMEvaluationDeserialization();

  /**
   * Write the Cap'n'Proto buffer to disk.
   */
  void read_from_file(const std::string & path);

private:

  /**
   * The RBSCMEvaluation object we will read into.
   */
  RBSCMEvaluation & _rb_scm_eval;

};
#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

/**
 * Load parameter ranges and discrete parameter values into an RBEvaluation
 * from the corresponding structure in the buffer.
 */
void load_parameter_ranges(RBParametrized & rb_evaluation,
                           RBData::ParameterRanges::Reader & parameter_ranges,
                           RBData::DiscreteParameterList::Reader & discrete_parameters_list);

/**
 * Load an RB evaluation from a corresponding reader structure in the buffer.
 */
template <typename RBEvaluationReaderNumber>
void load_rb_evaluation_data(RBEvaluation & rb_evaluation,
                             RBEvaluationReaderNumber & rb_evaluation_reader,
                             bool read_error_bound_data);

/**
 * Load an RB evaluation from a corresponding reader structure in the buffer.
 * Templated to deal with both Real and Complex numbers.
 */
template <typename RBEvaluationReaderNumber, typename TransRBEvaluationReaderNumber>
void load_transient_rb_evaluation_data(TransientRBEvaluation & trans_rb_eval,
                                       RBEvaluationReaderNumber & rb_evaluation_reader,
                                       TransRBEvaluationReaderNumber & trans_rb_eval_reader,
                                       bool read_error_bound_data);

/**
 * Load an EIM RB evaluation from a corresponding reader structure in the buffer.
 * Templated to deal with both Real and Complex numbers.
 */
template <typename RBEIMEvaluationReaderNumber>
void load_rb_eim_evaluation_data(RBEIMEvaluation & rb_eim_eval,
                                 RBEIMEvaluationReaderNumber & rb_eim_eval_reader);

#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)
/**
 * Load an SCM RB evaluation from a corresponding reader structure in the buffer.
 * Unlike the other functions above, this does not need
 * to be templated because an RBSCMEvaluation only stores
 * Real values, and hence doesn't depend on whether we're
 * using complex numbers or not.
 */
void load_rb_scm_evaluation_data(RBSCMEvaluation & rb_scm_eval,
                                 RBData::RBSCMEvaluation::Reader & rb_scm_eval_reader);
#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

/**
 * Helper function that loads point data.
 */
void load_point(RBData::Point3D::Reader point_reader, Point & point);

} // namespace RBDataDeserialization

} // namespace libMesh

#endif // #if defined(LIBMESH_HAVE_CAPNPROTO)

#endif // RB_COMPONENT_DATA_DESERIALIZATION_H
