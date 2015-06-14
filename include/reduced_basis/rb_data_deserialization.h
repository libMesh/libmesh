//#ifndef RB_DATA_DESERIALIZATION_H
//#define RB_DATA_DESERIALIZATION_H
//
//// libMesh/reduced_basis includes
//#include "libmesh/transient_rb_evaluation.h"
//#include "libmesh/rb_eim_evaluation.h"
//#include "libmesh/rb_scm_evaluation.h"
//#include "libmesh/rb_data_real.capnp.h"
//
//// Cap'n'Proto includes
//#include "capnp/common.h"
//
//// C++ includes
//
//namespace libMesh
//{
//
//namespace RBDataDeserialization
//{
//
///**
// * This class de-serializes an RBEvaluation object
// * using the Cap'n Proto library.
// */
//class RBEvaluationDeserialization
//{
//public:
//
//  /**
//   * Initialize a new buffer using the structure from the Cap'n'Proto schema
//   * described in rb_data.capnp.
//   */
//  RBEvaluationDeserialization(RBEvaluation& rb_eval);
//
//  /**
//   * Destructor.
//   */
//  virtual ~RBEvaluationDeserialization();
//  
//  /**
//   * Write the Cap'n'Proto buffer to disk.
//   */
//  void read_from_file(const std::string& path, bool read_error_bound_data);
//
//private:
//
//  /**
//   * The RBEvaluation object that we will read into.
//   */
//  RBEvaluation& _rb_eval;
//
//};
//
///**
// * This class de-serializes a TransientRBEvaluation object
// * using the Cap'n Proto library.
// */
//class TransientRBEvaluationDeserialization
//{
//public:
//
//  /**
//   * Initialize a new buffer using the structure from the Cap'n'Proto schema
//   * described in rb_data.capnp.
//   */
//  TransientRBEvaluationDeserialization(TransientRBEvaluation& trans_rb_eval);
//
//  /**
//   * Destructor.
//   */
//  virtual ~TransientRBEvaluationDeserialization();
//  
//  /**
//   * Write the Cap'n'Proto buffer to disk.
//   */
//  void read_from_file(const std::string& path, bool read_error_bound_data);
//
//private:
//
//  /**
//   * The TransientRBEvaluation object that we will read into.
//   */
//  TransientRBEvaluation& _trans_rb_eval;
//
//};
//
///**
// * This class de-serializes a RBEIMEvaluation object
// * using the Cap'n Proto library.
// */
//class RBEIMEvaluationDeserialization
//{
//public:
//
//  /**
//   * Initialize a new buffer using the structure from the Cap'n'Proto schema
//   * described in rb_data.capnp.
//   */
//  RBEIMEvaluationDeserialization(RBEIMEvaluation& trans_rb_eval);
//
//  /**
//   * Destructor.
//   */
//  virtual ~RBEIMEvaluationDeserialization();
//  
//  /**
//   * Write the Cap'n'Proto buffer to disk.
//   */
//  void read_from_file(const std::string& path, bool read_error_bound_data);
//
//private:
//
//  /**
//   * The RBEIMEvaluation object we will read into.
//   */
//  RBEIMEvaluation& _rb_eim_eval;
//
//};
//
///**
// * This class de-serializes a RBSCMEvaluation object
// * using the Cap'n Proto library.
// */
//class RBSCMEvaluationDeserialization
//{
//public:
//
//  /**
//   * Initialize a new buffer using the structure from the Cap'n'Proto schema
//   * described in rb_data.capnp.
//   */
//  RBSCMEvaluationDeserialization(RBSCMEvaluation& trans_rb_eval);
//
//  /**
//   * Destructor.
//   */
//  virtual ~RBSCMEvaluationDeserialization();
//  
//  /**
//   * Write the Cap'n'Proto buffer to disk.
//   */
//  void read_from_file(const std::string& path, bool read_error_bound_data);
//
//private:
//
//  /**
//   * The RBSCMEvaluation object we will read into.
//   */
//  RBSCMEvaluation& _rb_scm_eval;
//
//};
//
///**
// * Load an RB evaluation from a corresponding reader structure in the buffer.
// */
//void load_rb_evaluation_data(
//  RBEvaluation& rb_evaluation,
//  RBData::RBEvaluation::Reader& rb_evaluation_reader,
//  bool read_error_bound_data);
//
///**
// * Load an RB evaluation from a corresponding reader structure in the buffer.
// */
//void load_transient_rb_evaluation_data(
//  TransientRBEvaluation& trans_rb_eval,
//  RBData::RBEvaluation::Reader& rb_evaluation_reader,
//  RBData::TransientRBEvaluation::Reader& trans_rb_eval_reader,
//  bool read_error_bound_data);
//
///**
// * Load an EIM RB evaluation from a corresponding reader structure in the buffer.
// */
//void load_rb_eim_evaluation_data(
//  RBEIMEvaluation& rb_eim_eval,
//  RBData::RBEvaluation::Reader& rb_evaluation_reader,
//  RBData::RBEIMEvaluation::Reader& rb_eim_eval_reader,
//  bool read_error_bound_data);
//
///**
// * Load an SCM RB evaluation from a corresponding reader structure in the buffer.
// */
//void load_rb_scm_evaluation_data(
//  RBSCMEvaluation& rb_scm_eval,
//  RBData::RBSCMEvaluation::Reader& rb_scm_eval_reader,
//  bool read_error_bound_data);
//
///**
// * Load parameter ranges and discrete parameter values into an RBEvaluation
// * from the corresponding structure in the buffer.
// */
//void load_parameter_ranges(
//  RBParametrized& rb_evaluation,
//  RBData::ParameterRanges::Reader& parameter_ranges,
//  RBData::DiscreteParameterList::Reader& discrete_parameters_list);
//
//} // namespace RBDataDeserialization
//
//} // namespace libMesh
//
//#endif // RB_COMPONENT_DATA_DESERIALIZATION_H
