#ifndef ONCE_FPARSERAD_H_
#define ONCE_FPARSERAD_H_

#include "fparser.hh"
#include <exception>

template<typename Value_t>
class ADImplementation;

template<typename Value_t>
class FunctionParserADBase : public FunctionParserBase<Value_t>
{
public:
  FunctionParserADBase();
  FunctionParserADBase(const FunctionParserADBase& cpy);
  virtual ~FunctionParserADBase();

  /**
   * auto-differentiate for var
   */
  int AutoDiff(const std::string & var_name);

  /**
   * add another variable
   */
  bool AddVariable(const std::string & var_name);

  /**
   * check if the function is equal to 0
   * This is a common case for vanishing derivatives. This relies on the
   * function to be optimized.
   */
  bool isZero();

  /**
   * set the bytecode of this function to return constant zero.
   * this provides a well defined state in case AutoDiff fails
   */
  void setZero();

  // feature flags for this parser
  enum ADFlags {
    /**
     * In certain applications derivatives are built proactively and may never be used.
     * We silence all AutoDiff exceptions in that case to avoid confusing the user.
     */
    ADSilenceErrors = 1,
    /**
     * Immediately apply the optimizer grammars to the derivative tree structure.
     * This saves a round trip to byte code.
     */
    ADAutoOptimize = 2,
    /**
     * Use cached JIT compiled function files, This bypasses the compilation stage
     * for all further runs, after the JIT compilation ran successfully at least once.
     */
    ADJITCache = 4,
    /**
     * Use cached bytecode for the derivatives. This bypasses the automatic differentiation,
     * (optimization), and byte code synthesis.
     */
    ADCacheDerivatives = 8
  };

  /**
   * Set the feature flags for this parser (this way gives us better control over default values)
   */
  void SetADFlags(int flags, bool turnon = true) {
    if (turnon)
      mADFlags |= flags;
    else
      mADFlags &= ~flags;
  }
  void UnsetADFlags(int flags) { mADFlags &= ~flags; }
  void ClearADFlags() { mADFlags = 0; }

  /**
   * compile the current function, or load a previously compiled copy.
   * Warning: When re-using an FParser function object by parsing a new expression
   *          the previously JIT compiled function will continue to be Evaled until the
   *          JITCompile method is called again.
   */
  bool JITCompile();

  /**
   * wrap Optimize of the parent class to check for a JIT compiled version and redo
   * the compilation after Optimization
   */
  void Optimize();

  /**
   * write the full state of the current FParser object to a stream
   */
  void Serialize(std::ostream &);

  /**
   * restore the full state of the current FParser object from a stream
   */
  void Unserialize(std::istream &);

#if LIBMESH_HAVE_FPARSER_JIT
  /**
   * Overwrite the Exec function with one that tests for a JIT compiled version
   * and uses that if it exists
   */
  Value_t Eval(const Value_t* Vars);
#endif

  /**
   * look up the opcode number for a given variable name
   * throws UnknownVariableException if the variable is not found
   */
  unsigned int LookUpVarOpcode(const std::string & var_name);

  /**
   * register a dependency between variables so that da/db = c
   */
  void RegisterDerivative(const std::string & a, const std::string & b, const std::string & c);

private:
  /// helper function to perform the JIT compilation (needs the Value_t typename as a string)
  bool JITCompileHelper(const std::string &);

  /// JIT function pointer
  Value_t (*compiledFunction)(const Value_t *, const Value_t *, const Value_t);

  /// pointer to the mImmed values (or NULL if the mImmed vector is empty)
  Value_t * pImmed;

  // user function plog
  static Value_t fp_plog(const Value_t * params);

  // user function erf
  static Value_t fp_erf(const Value_t * params);

  // function ID for the plog function
  unsigned int mFPlog;

  // function ID for the erf function
  unsigned int mFErf;

  // flags that control cache bahavior, optimization, and error reporting
  int mADFlags;

  // registered derivative table, and entry structure
  struct VariableDerivative {
    unsigned int var, dependence, derivative;
  };
  std::vector<VariableDerivative> mRegisteredDerivatives;

  // private implementaion of the automatic differentiation algorithm
  ADImplementation<Value_t> * ad;

  // the firewalled implementation class of the AD algorithm has full access to the FParser object
  friend class ADImplementation<Value_t>;

  // Exceptions
  class UnknownVariable : public std::exception {
    virtual const char* what() const throw() { return "Unknown variable"; }
  } UnknownVariableException;
  class UnknownSerializationVersion : public std::exception {
    virtual const char* what() const throw() { return "Unknown serialization file version"; }
  } UnknownSerializationVersionException;
};


class FunctionParserAD: public FunctionParserADBase<double> {};
class FunctionParserAD_f: public FunctionParserADBase<float> {};

#endif //ONCE_FPARSERAD_H_
