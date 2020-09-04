#ifndef ONCE_FPARSERAD_H_
#define ONCE_FPARSERAD_H_

#include "fparser.hh"
#include <exception>
#include <iostream>
#include <fstream>

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
   * This class manages its own memory, so the compiler-generated copy
   * assignment, move assignment, and move constructor implementations
   * are not safe to use.  We therefore explicitly delete them so they
   * can't be called accidentally.
   */
  FunctionParserADBase (FunctionParserADBase &&) = delete;
  FunctionParserADBase & operator= (const FunctionParserADBase &) = delete;
  FunctionParserADBase & operator= (FunctionParserADBase &&) = delete;

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
   * check if the function's byte code is empty.
   */
  bool isEmpty() { return this->mData->mByteCode.empty(); }

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

protected:

#if LIBMESH_HAVE_FPARSER_JIT
  /// return a SHA1 hash for the current bytecode and value type name
  std::string JITCodeHash(const std::string & value_type_name);

  /// write generated C++ code to stream
  bool JITCodeGen(std::ostream & ccout, const std::string & fname, const std::string & Value_t_name);

  /// helper function to perform the JIT compilation (needs the Value_t typename as a string)
  bool JITCompileHelper(const std::string & Value_t_name,
                        const std::string & extra_options = "",
                        const std::string & extra_headers = "");
#endif // LIBMESH_HAVE_FPARSER_JIT

  /// function pointer type alias. This permits a Real Value_t function to be compiled
  /// to support dual numbers
  template <typename ActualValue_t>
  using CompiledFunctionPtr = void (*)(ActualValue_t *, const ActualValue_t *,
                                       const Value_t *, const Value_t);

  /// update pointer to immediate data
  void updatePImmed();

  /// clear the runtime evaluation error flag
  void clearEvalError() { this->mData->mEvalErrorType = 0; }

  /// JIT function pointer
  void *compiledFunction;

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

#ifdef LIBMESH_HAVE_FPARSER_JIT

/// Forward declare SHA1 hash object
class SHA1;

/// Namespacing the utility classes (rather than nesting them in a templated class)
namespace FParserJIT
{
/// Simplified C++ interface to lib SHA1
class Hash
{
public:
  Hash();
  ~Hash();

  template <typename T>
  void addData(const T & v);

  std::string get();

protected:
  /// the actual lib SHA1 call is in the helper so that we don't need to make lib/sha1.h available
  void addDataHelper(const char * start, std::size_t size);

  SHA1 * _sha1;
};

template <typename T>
void Hash::addData(const T & v)
{
  auto start = v.data();
  std::size_t size = v.size() * sizeof(*start);
  if (size > 0)
    addDataHelper(reinterpret_cast<const char *>(start), size);
}

/// Handle compilation, caching, and temporary files
class Compiler
{
public:
  Compiler(const std::string & master_hash = "");
  ~Compiler();
  std::ostream & source();

  bool probeCache();
  bool run(const std::string & compiler_options = "");
  void * getFunction(const std::string & fname);

protected:
  std::ofstream _ccout;
  void * _lib;
  const std::string _jitdir;
  std::string _ccname;
  std::string _objectname;
  std::string _object_so;
  bool _success;

  const std::string _master_hash;
  const bool _use_cache;
};
} // namespace FParserJIT

#endif // LIBMESH_HAVE_FPARSER_JIT

class FunctionParserAD: public FunctionParserADBase<double> {};
class FunctionParserAD_f: public FunctionParserADBase<float> {};

#endif //ONCE_FPARSERAD_H_
