#ifndef ONCE_FPARSERAD_H_
#define ONCE_FPARSERAD_H_

#include "fparser.hh"
#include <exception>

template<typename Value_t>
class FunctionParserADBase : public FunctionParserBase<Value_t>
{
public:
  FunctionParserADBase();
  FunctionParserADBase(const FunctionParserADBase& cpy);

  /**
   * auto-differentiate for var
   */
  int AutoDiff(const std::string& var);

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

protected:
  /**
   * A list of opcodes and immediate values
   */
  struct OpcodePacket {
    unsigned first, index;
    Value_t second;
    OpcodePacket() : first(0), second(Value_t()), index(0) {}
    OpcodePacket(unsigned _first, Value_t _second, unsigned _index) : first(_first), second(_second), index(_index) {}
  };
  struct OpcodePlain : OpcodePacket {
    OpcodePlain(unsigned _first) : OpcodePacket(_first, Value_t(), 0) {}
  };
  struct OpcodeImmediate : OpcodePacket {
    OpcodeImmediate(Value_t _second);
  };
  struct OpcodeFCall : OpcodePacket {
    OpcodeFCall(unsigned _index);
  };

  typedef std::vector<OpcodePacket> DiffProgramFragment;
  typedef std::pair<typename DiffProgramFragment::const_iterator,
                    typename DiffProgramFragment::const_iterator> Interval;

  /**
   * Recursively differentiate functions from the outside (end of program)
   * to the inside.
   */
  DiffProgramFragment DiffFunction(const DiffProgramFragment & orig);

  /**
   * how much does the current opcode move the stack pointer
   */
  int OpcodeSize(const OpcodePacket & p);

private:
  typename FunctionParserBase<Value_t>::Data * mData;

  /// get the preceeding argument
  Interval GetArgument(const DiffProgramFragment & code);

  /// get argument n before
  Interval GetArgument(const DiffProgramFragment & code, unsigned int index);

  /// remove the code fragments that generate the n previous arguments
  //std::vector<DiffProgramFragment> PopArguments(DiffProgramFragment & code, unsigned int index);

  /// variable to diff for
  unsigned mVarOp;

  /// write the DiffProgramFragement into the internal bytecode storage
  void Commit(const DiffProgramFragment & code);

  /// the function indices of the plog dplog and step functions
  unsigned mFStep, mFPlog;

  /**
   * In certain applications derivatives are built proactively and may never be used.
   * We silence all AutoDiff exceptions in that case to avoid confusing the user.
   */
  bool mSilenceErrors;

  // user functions plog and step
  static Value_t fp_step(const Value_t * params);
  static Value_t fp_plog(const Value_t * params);

  // Exceptions
  class UnsupportedOpcode : public std::exception {
    virtual const char* what() const throw() { return "Unsupported opcode"; }
  } UnsupportedOpcodeException;
  class StackExhausted : public std::exception {
    virtual const char* what() const throw() { return "Stack exhausted."; }
  } StackExhaustedException;
  class EmptyProgram : public std::exception {
    virtual const char* what() const throw() { return "Empty programm passed in."; }
  } EmptyProgramException;
  class UnsupportedArgumentCount : public std::exception {
    virtual const char* what() const throw() { return "Unsupported argument count."; }
  } UnsupportedArgumentCountException;
};


class FunctionParserAD: public FunctionParserADBase<double> {};
class FunctionParserAD_f: public FunctionParserADBase<float> {};

#endif //ONCE_FPARSERAD_H_
