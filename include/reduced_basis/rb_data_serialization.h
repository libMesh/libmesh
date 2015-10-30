#ifndef RB_DATA_SERIALIZATION_H
#define RB_DATA_SERIALIZATION_H

#include "libmesh/libmesh_config.h"
#if defined(LIBMESH_HAVE_CAPNPROTO)

// C++ includes
#include <string>

// libMesh/reduced_basis includes
#include "libmesh/rb_evaluation.h"
#include "libmesh/transient_rb_evaluation.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/rb_scm_evaluation.h"
#include "libmesh/rb_data.capnp.h"

// Cap'n'Proto includes
#include "capnp/message.h"

namespace libMesh
{

namespace RBDataSerialization
{

/**
 * This class serializes an RBEvaluation object
 * using the Cap'n Proto library.
 */
class RBEvaluationSerialization
{
public:

  /**
   * Initialize a new buffer using the structure from the Cap'n'Proto schema
   * described in rb_data.capnp.
   */
  RBEvaluationSerialization(RBEvaluation& rb_eval);

  /**
   * Destructor.
   */
  virtual ~RBEvaluationSerialization();

  /**
   * Write the Cap'n'Proto buffer to disk.
   */
  void write_to_file(const std::string& path);

private:

  /**
   * The RBEvaluation object that will be written to disk.
   */
  RBEvaluation& _rb_eval;

};

/**
 * This class serializes a TransientRBEvaluation object
 * using the Cap'n Proto library.
 */
class TransientRBEvaluationSerialization
{
public:

  /**
   * Initialize a new buffer using the structure from the Cap'n'Proto schema
   * described in rb_data.capnp.
   */
  TransientRBEvaluationSerialization(TransientRBEvaluation& rb_eval);

  /**
   * Destructor.
   */
  virtual ~TransientRBEvaluationSerialization();

  /**
   * Write the Cap'n'Proto buffer to disk.
   */
  void write_to_file(const std::string& path);

private:

  /**
   * The RBEvaluation object that will be written to disk.
   */
  TransientRBEvaluation& _trans_rb_eval;

};

/**
 * This class serializes an RBEIMEvaluation object
 * using the Cap'n Proto library.
 */
class RBEIMEvaluationSerialization
{
public:

  /**
   * Initialize a new buffer using the structure from the Cap'n'Proto schema
   * described in rb_data.capnp.
   */
  RBEIMEvaluationSerialization(RBEIMEvaluation& rb_eval);

  /**
   * Destructor.
   */
  virtual ~RBEIMEvaluationSerialization();

  /**
   * Write the Cap'n'Proto buffer to disk.
   */
  void write_to_file(const std::string& path);

private:

  /**
   * The RBEvaluation object that will be written to disk.
   */
  RBEIMEvaluation& _rb_eim_eval;

};

// RBSCMEvaluation should only be available
// if SLEPc and GLPK support is enabled.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

/**
 * This class serializes an RBSCMEvaluation object
 * using the Cap'n Proto library.
 */
class RBSCMEvaluationSerialization
{
public:

  /**
   * Initialize a new buffer using the structure from the Cap'n'Proto schema
   * described in rb_data.capnp.
   */
  RBSCMEvaluationSerialization(RBSCMEvaluation& rb_eval);

  /**
   * Destructor.
   */
  virtual ~RBSCMEvaluationSerialization();

  /**
   * Write the Cap'n'Proto buffer to disk.
   */
  void write_to_file(const std::string& path);

private:

  /**
   * The RBEvaluation object that will be written to disk.
   */
  RBSCMEvaluation& _rb_scm_eval;

};
#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

/**
 * Add parameter ranges for continuous and discrete parameters.
 */
void add_parameter_ranges_to_builder(const RBParametrized& rb_evaluation,
                                     RBData::ParameterRanges::Builder& parameter_ranges,
                                     RBData::DiscreteParameterList::Builder& discrete_parameters_list);

/**
 * Add data for an RBEvaluation to the builder.
 */
template <typename RBEvaluationBuilderNumber>
void add_rb_evaluation_data_to_builder(RBEvaluation& rb_eval,
                                       RBEvaluationBuilderNumber& rb_eval_builder);

/**
 * Add data for a TransientRBEvaluation to the builder.
 * Templated to deal with both Real and Complex numbers.
 */
template <typename RBEvaluationBuilderNumber, typename TransRBEvaluationBuilderNumber>
void add_transient_rb_evaluation_data_to_builder(TransientRBEvaluation& trans_rb_eval,
                                                 RBEvaluationBuilderNumber& rb_eval_builder,
                                                 TransRBEvaluationBuilderNumber& trans_rb_eval_builder);

/**
 * Add data for an RBEIMEvaluation to the builder.
 * Templated to deal with both Real and Complex numbers.
 */
template <typename RBEvaluationBuilderNumber, typename RBEIMEvaluationBuilderNumber>
void add_rb_eim_evaluation_data_to_builder(RBEIMEvaluation& rb_eim_eval,
                                           RBEvaluationBuilderNumber& rb_eval_builder,
                                           RBEIMEvaluationBuilderNumber& rb_eim_eval_builder);

#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)
/**
 * Add data for an RBSCMEvaluation to the builder.
 * Unlike the other functions above, this does not need
 * to be templated because an RBSCMEvaluation only stores
 * Real values, and hence doesn't depend on whether we're
 * using complex numbers or not.
 */
void add_rb_scm_evaluation_data_to_builder(RBSCMEvaluation& rb_scm_eval,
                                           RBData::RBSCMEvaluation::Builder& rb_scm_eval_builder);
#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

/**
 * Helper function that adds point data.
 */
void add_point_to_builder(const Point& point,
                          RBData::Point3D::Builder point_builder);

/**
 * Helper function that adds element data.
 */
void add_elem_to_builder(const libMesh::Elem& elem,
                         RBData::MeshElem::Builder mesh_elem_builder);

} // namespace RBDataSerialization

} // namespace libMesh

#endif // #if defined(LIBMESH_HAVE_CAPNPROTO)

#endif // RB_DATA_SERIALIZATION_H
