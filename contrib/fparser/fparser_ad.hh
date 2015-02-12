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

  /**
   * Do not report any error messages to the console
   */
  void silenceAutoDiffErrors(bool _silence = true) { mSilenceErrors = _silence; }

  /**
   * compile the current function, or load a previously compiled copy.
   * Warning: When re-using an FParser function object by parsing a new expression
   *          the previously JIT compiled function will continue to be Evaled until the
   *          JITCompile method is called again.
   */
  bool JITCompile(bool cacheFunction = true);

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
  typename FunctionParserBase<Value_t>::Data * mData;

  /// helper function to perform the JIT compilation (needs the Value_t typename as a string)
  bool JITCompileHelper(const std::string &, bool);

  /// JIT function pointer
  Value_t (*compiledFunction)(const Value_t *, const Value_t *, const Value_t);

  /**
   * In certain applications derivatives are built proactively and may never be used.
   * We silence all AutoDiff exceptions in that case to avoid confusing the user.
   */
  bool mSilenceErrors;

  // user function plog
  static Value_t fp_plog(const Value_t * params);

  // function ID for the plog function
  unsigned int mFPlog;

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
};


class FunctionParserAD: public FunctionParserADBase<double> {};
class FunctionParserAD_f: public FunctionParserADBase<float> {};

#endif //ONCE_FPARSERAD_H_
