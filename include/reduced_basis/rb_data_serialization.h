#ifndef RB_DATA_SERIALIZATION_H
#define RB_DATA_SERIALIZATION_H

// C++ includes
#include <string>

// libMesh/reduced_basis includes
#include "libmesh/rb_evaluation.h"
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

///**
// * This class serializes a TransientRBEvaluation object
// * using the Cap'n Proto library.
// */ 
//class TransientRBEvaluationSerialization
//{
//public:
//
//  /**
//   * Initialize a new buffer using the structure from the Cap'n'Proto schema
//   * described in rb_data.capnp.
//   */
//  TransientRBEvaluationSerialization(TransientRBEvaluation& rb_eval);
//
//  /**
//   * Destructor.
//   */
//  virtual ~TransientRBEvaluationSerialization();
//  
//  /**
//   * Write the Cap'n'Proto buffer to disk.
//   */
//  void write_to_file(const std::string& path);
//
//private:
//
//  /**
//   * The RBEvaluation object that will be written to disk.
//   */
//  TransientRBEvaluation& _trans_rb_eval;
//
//};
//
///**
// * This class serializes an RBEIMEvaluation object
// * using the Cap'n Proto library.
// */ 
//class RBEIMEvaluationSerialization
//{
//public:
//
//  /**
//   * Initialize a new buffer using the structure from the Cap'n'Proto schema
//   * described in rb_data.capnp.
//   */
//  RBEIMEvaluationSerialization(RBEIMEvaluation& rb_eval);
//
//  /**
//   * Destructor.
//   */
//  virtual ~RBEIMEvaluationSerialization();
//  
//  /**
//   * Write the Cap'n'Proto buffer to disk.
//   */
//  void write_to_file(const std::string& path);
//
//private:
//
//  /**
//   * The RBEvaluation object that will be written to disk.
//   */
//  RBEIMEvaluation& _rb_eim_eval;
//
//};
//
///**
// * This class serializes an RBSCMEvaluation object
// * using the Cap'n Proto library.
// */ 
//class RBSCMEvaluationSerialization
//{
//public:
//
//  /**
//   * Initialize a new buffer using the structure from the Cap'n'Proto schema
//   * described in rb_data.capnp.
//   */
//  RBSCMEvaluationSerialization(RBSCMEvaluation& rb_eval);
//
//  /**
//   * Destructor.
//   */
//  virtual ~RBSCMEvaluationSerialization();
//  
//  /**
//   * Write the Cap'n'Proto buffer to disk.
//   */
//  void write_to_file(const std::string& path);
//
//private:
//
//  /**
//   * The RBEvaluation object that will be written to disk.
//   */
//  RBEIMEvaluation& _rb_scm_eval;
//
//};

/**
 * Add parameter ranges for continuous and discrete parameters.
 */
void add_parameter_ranges_to_builder(
  const RBParametrized& rb_evaluation,
  RBData::ParameterRanges::Builder& parameter_ranges,
  RBData::DiscreteParameterList::Builder& discrete_parameters_list);

/**
 * Add data for an RBEvaluation to the builder.
 */
template <typename RBEvaluationBuilderNumber>
void add_rb_evaluation_data_to_builder(
  RBEvaluation& rb_eval,
  RBEvaluationBuilderNumber& rb_eval_builder);

///**
// * Add data for a TransientRBEvaluation to the builder.
// */
//void add_transient_rb_evaluation_data_to_builder(
//  TransientRBEvaluation& trans_rb_eval,
//  RBData::RBEvaluation::Builder& rb_eval_builder,
//  RBData::TransientRBEvaluation::Builder& trans_rb_eval_builder);
//
///**
// * Add data for an RBEIMEvaluation to the builder.
// */
//void add_rb_eim_evaluation_data_to_builder(
//  RBEIMEvaluation& rb_eim_eval,
//  RBData::RBEvaluation::Builder& rb_eval_builder,
//  RBData::RBEIMEvaluation::Builder& rb_eim_eval_builder);
//
///**
// * Add data for an RBSCMEvaluation to the builder.
// */
//void add_rb_scm_evaluation_data_to_builder(
//  SCMRBEvaluation& rb_scm_eval,
//  RBData::RBSCMEvaluation::Builder& rb_scm_eval_builder);

} // namespace RBDataSerialization

} // namespace libMesh

#endif // RB_DATA_SERIALIZATION_H
