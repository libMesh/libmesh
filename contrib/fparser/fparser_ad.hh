#ifndef ONCE_FPARSERAD_H_
#define ONCE_FPARSERAD_H_

#include "fparser.hh"
//#include "libmesh_common.h"
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

protected:
  /**
   * A list of opcodes and immediate values
   */
  typedef std::pair<unsigned, Value_t> OpcodeDataPair;
  typedef std::vector<OpcodeDataPair> DiffProgramFragment;
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
  int OpcodeSize(unsigned op);

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

  // write the DiffProgramFragement into the internal bytecode storage
  void Commit(const DiffProgramFragment & code);

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
