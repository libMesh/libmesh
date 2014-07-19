#include "fparser_ad.hh"
#include "extrasrc/fpaux.hh"
#include "extrasrc/fptypes.hh"
#include <stdlib.h>

using namespace FUNCTIONPARSERTYPES;

#include <iostream>

template<typename Value_t>
FunctionParserADBase<Value_t>::FunctionParserADBase() :
    FunctionParserBase<Value_t>(),
    mData(this->getParserData()),
    mSilenceErrors(false)
{
  mFPlog  = mData->mFuncPtrs.size();
  this->AddFunction("plog", fp_plog, 2);
}

template<typename Value_t>
FunctionParserADBase<Value_t>::FunctionParserADBase(const FunctionParserADBase& cpy) :
    FunctionParserBase<Value_t>(cpy),
    mData(this->getParserData()),
    mSilenceErrors(cpy.mSilenceErrors)
{
  mFPlog  = cpy.mFPlog;
}

template<typename Value_t>
Value_t FunctionParserADBase<Value_t>::fp_step(const Value_t * params)
{
  return params[0] > Value_t(0.0) ? Value_t(1.0) : Value_t(0.0);
}

template<typename Value_t>
Value_t FunctionParserADBase<Value_t>::fp_plog(const Value_t * params)
{
  return params[0] < params[1] ? fp_log(params[1]) + (params[0] - params[1]) / params[1] : fp_log(params[0]);
}

template<typename Value_t>
int FunctionParserADBase<Value_t>::OpcodeSize(const OpcodePacket & p)
{
  unsigned op = p.first, index = p.index;

  if (int(op) >= VarBegin || op == cImmed) {
    return -1;
  } else if (op == cFCall && index == mFPlog) {
    return mData->mFuncPtrs[index].mParams;
  } else {
    switch(op)
    {
      // these opcode takes one argument off the stack
      case cInv: case cNeg:
      case cSqr: case cAbs:
      case cSqrt: case cRSqrt: case cCbrt:
      case cExp: case cExp2:
      case cLog: case cLog2: case cLog10:
      case cSin: case cCos: case cTan:
      case cAsin: case cAcos: case cAtan:
      case cInt: case cFloor: case cCeil: case cTrunc:
        return 1;

      // these opcode takes two arguments off the stack
      case cAdd: case cSub: case cRSub:
      case cMul: case cDiv: case cRDiv:
      case cPow: case cHypot:
      case cLess: case cLessOrEq: case cGreater: case cGreaterOrEq:
        return 2;

      default:
        throw UnsupportedOpcodeException;
    }
  }
}

template<typename Value_t>
typename FunctionParserADBase<Value_t>::Interval
FunctionParserADBase<Value_t>::GetArgument(const DiffProgramFragment & orig)
{
  // count the number of elements on the stack (needs to reach 1)
  int stack_size = 0;

  // Extract the opcode sequence(s) that is responsible for the
  // top stack entry
  typename DiffProgramFragment::const_iterator ip = orig.end();
  do {
    // take one step back
    if (ip == orig.begin())
      throw StackExhaustedException;
    ip--;

    // a size two opcode needs two elements from the stack, but also puts one back on!
    int this_size = OpcodeSize(*ip);
    stack_size -= this_size > 0 ? this_size-1 : this_size;
  } while (stack_size < 1);

  return Interval(ip, orig.end());
}

template<typename Value_t>
typename FunctionParserADBase<Value_t>::Interval
FunctionParserADBase<Value_t>::GetArgument(const DiffProgramFragment & orig, unsigned int index)
{
  // argument loction
  Interval arg;
  DiffProgramFragment head;
  std::copy(orig.begin(), orig.end()-1, std::back_inserter(head));

  // iterate over the past arguments
  int i = index;
  do
  {
    arg = GetArgument(head);
    head.resize(std::distance<typename DiffProgramFragment::const_iterator>(head.begin(), arg.first));
  } while (--i >= 0);

  return arg;
}

template<typename Value_t>
typename FunctionParserADBase<Value_t>::DiffProgramFragment
FunctionParserADBase<Value_t>::DiffFunction(const DiffProgramFragment & orig)
{
  // check for empty DiffProgramFragments
  if (orig.empty())
    throw EmptyProgramException;

  // this stores the opcode of the differentiated function
  DiffProgramFragment outer, prog_a, prog_b, prog_da, prog_db;

  // size of the current end opcode
  int op_size = OpcodeSize(orig.back());

  // current opcode
  unsigned op = orig.back().first;
  unsigned findex = orig.back().index;

  // variable or immediate value
  if (op_size == -1)
  {
    if (op == cImmed)
      outer.push_back(OpcodeImmediate(Value_t(0)));
    else
    {
      if (op == mVarOp)
        outer.push_back(OpcodeImmediate(Value_t(1)));
      else
        outer.push_back(OpcodeImmediate(Value_t(0)));
    }

    return outer;
  }

  // Extract the opcode sequence(s) that generates
  // the argument(s) for the current opcode

  // get the opcode sequence preceeding the current opcode
  DiffProgramFragment head, head2;
  std::copy(orig.begin(), orig.end()-1, std::back_inserter(head));

  // get the opcode interval of the argument immediately preceeding the current opcode
  Interval arg = GetArgument(head);

  if (op_size == 1)
    std::copy(arg.first, arg.second, std::back_inserter(prog_a));
  else if (op_size == 2)
  {
    std::copy(arg.first, arg.second, std::back_inserter(prog_b));
    // get opcode squence before that argument interval
    std::copy<typename DiffProgramFragment::const_iterator>(head.begin(), arg.first, std::back_inserter(head2));
    arg = GetArgument(head2);
    std::copy(arg.first, arg.second, std::back_inserter(prog_a));
  }
  else
    throw UnsupportedArgumentCountException;

  // create derivatives
  switch (op)
  {
    case cRSub:
      // db - da
    case cAdd:
      // da + db
    case cSub:
      // da - db
      prog_da = DiffFunction(prog_a);
      prog_db = DiffFunction(prog_b);
      outer = prog_da;
      outer.insert(outer.end(), prog_db.begin(), prog_db.end());
      outer.push_back(OpcodePlain(op));
      return outer;

    case cMul:
      // a*db + da*b
      prog_da = DiffFunction(prog_a);
      prog_db = DiffFunction(prog_b);
      outer = prog_a;
      outer.insert(outer.end(), prog_db.begin(), prog_db.end());
      outer.push_back(OpcodePlain(cMul));
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.push_back(OpcodePlain(cMul));
      outer.push_back(OpcodePlain(cAdd));
      return outer;

    case cSqr:
      // a*da*2
      prog_da = DiffFunction(prog_a);
      outer = prog_a;
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.push_back(OpcodePlain(cMul));
      outer.push_back(OpcodeImmediate(2));
      outer.push_back(OpcodePlain(cMul));
      return outer;

    // also capture cRdiv here but switch a and b first!
    case cRDiv:
      std::swap(prog_a, prog_b);
    case cDiv:
      // db/a - a*db/b^2
      prog_da = DiffFunction(prog_a);
      prog_db = DiffFunction(prog_b);
      outer = prog_db;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cDiv));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.insert(outer.end(), prog_db.begin(), prog_db.end());
      outer.push_back(OpcodePlain(cMul));
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.push_back(OpcodePlain(cSqr));
      outer.push_back(OpcodePlain(cDiv));
      outer.push_back(OpcodePlain(cSub));
      return outer;

    case cLog2:
    case cLog10:
      // *1/ln(base)
    case cLog:
      // da/a
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cDiv));
      if (op == cLog2) {
        outer.push_back(OpcodeImmediate(std::log(Value_t(2))));
        outer.push_back(OpcodePlain(cDiv));
      }
      else if (op == cLog10) {
        outer.push_back(OpcodeImmediate(std::log(Value_t(10))));
        outer.push_back(OpcodePlain(cDiv));
      }
      return outer;

    case cExp:
    case cExp2:
      // da*exp(a)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(op));
      outer.push_back(OpcodePlain(cMul));
      if (op == cExp2) {
        outer.push_back(OpcodeImmediate(std::log(Value_t(2))));
        outer.push_back(OpcodePlain(cMul));
      }
      return outer;

    case cNeg:
      // -da
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodePlain(cNeg));
      return outer;

    case cInv:
      // -da/a^2
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cSqr));
      outer.push_back(OpcodePlain(cDiv));
      outer.push_back(OpcodePlain(cNeg));
      return outer;

    case cPow:
      // a**b * ( db*log(a) + b*da/a)
      prog_da = DiffFunction(prog_a);
      prog_db = DiffFunction(prog_b);
      outer = prog_a;
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.push_back(OpcodePlain(cPow));
      outer.insert(outer.end(), prog_db.begin(), prog_db.end());
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cLog));
      outer.push_back(OpcodePlain(cMul));
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.push_back(OpcodePlain(cMul));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cDiv));
      outer.push_back(OpcodePlain(cAdd));
      outer.push_back(OpcodePlain(cMul));
      return outer;

    case cSin :
      // da*cos(a)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cCos));
      outer.push_back(OpcodePlain(cMul));
      return outer;

    case cCos :
      // -da*sin(a)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cSin));
      outer.push_back(OpcodePlain(cMul));
      outer.push_back(OpcodePlain(cNeg));
      return outer;

    case cTan :
      // (tan(a)^2+1)*da
      prog_da = DiffFunction(prog_a);
      outer = prog_a;
      outer.push_back(OpcodePlain(cTan));
      outer.push_back(OpcodePlain(cSqr));
      outer.push_back(OpcodeImmediate(Value_t(1)));
      outer.push_back(OpcodePlain(cAdd));
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.push_back(OpcodePlain(cMul));
      return outer;

    case cAtan :
      // da/(a^2+1)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cSqr));
      outer.push_back(OpcodeImmediate(Value_t(1)));
      outer.push_back(OpcodePlain(cAdd));
      outer.push_back(OpcodePlain(cDiv));
      return outer;

    case cAsin:
      // da/sqrt(1-a^2)
    case cAcos:
      // -da/sqrt(1-a^2)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodeImmediate(Value_t(1)));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cSqr));
      outer.push_back(OpcodePlain(cSub));
      outer.push_back(OpcodePlain(cSqrt));
      outer.push_back(OpcodePlain(cDiv));
      if (op == cAcos)
        outer.push_back(OpcodePlain(cNeg));
      return outer;

    case cSqrt :
      // da/(2*sqrt(a))
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodeImmediate(Value_t(2)));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cSqrt));
      outer.push_back(OpcodePlain(cMul));
      outer.push_back(OpcodePlain(cDiv));
      return outer;

    case cRSqrt :
      // -da/(2*a^(3/2)))
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodeImmediate(Value_t(2)));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeImmediate(Value_t(3)));
      outer.push_back(OpcodePlain(cPow));
      outer.push_back(OpcodePlain(cSqrt));
      outer.push_back(OpcodePlain(cMul));
      outer.push_back(OpcodePlain(cDiv));
      outer.push_back(OpcodePlain(cNeg));
      return outer;

    case cCbrt :
      // da/(3*a^(2/3))
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodeImmediate(Value_t(3)));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cCbrt));
      outer.push_back(OpcodePlain(cSqr));
      outer.push_back(OpcodePlain(cMul));
      outer.push_back(OpcodePlain(cDiv));
      return outer;

    case cAbs:
      // da*a/|a|
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodePlain(cDup));
      outer.push_back(OpcodePlain(cAbs));
      outer.push_back(OpcodePlain(cDiv));
      outer.push_back(OpcodePlain(cMul));
      return outer;

    case cHypot:
      // (a*da+b*db)/hypot(a,b)
      prog_da = DiffFunction(prog_a);
      prog_db = DiffFunction(prog_b);
      outer = prog_a;
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.push_back(OpcodePlain(cMul));
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.insert(outer.end(), prog_db.begin(), prog_db.end());
      outer.push_back(OpcodePlain(cMul));
      outer.push_back(OpcodePlain(cAdd));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.push_back(OpcodePlain(cHypot));
      outer.push_back(OpcodePlain(cDiv));
      return outer;

    // no idea why anyone would like to diff those
    // (I'll pretend the discontinuities don't exist -
    // FP could not deal with them in a useful way in any case)
    case cInt:
    case cFloor:
    case cCeil:
    case cTrunc:
    case cLess:
    case cLessOrEq:
    case cGreater:
    case cGreaterOrEq:
      outer.push_back(OpcodeImmediate(Value_t(0)));
      return outer;
    case cFCall:
      if (findex == mFPlog)
      {
        // we assume that the second argument to plog is constant for now (TODO?).
        // da/a
        prog_da = DiffFunction(prog_a);
        outer = prog_da;
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.push_back(OpcodePlain(cDiv));
        // a>b
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cGreater));
        // *
        outer.push_back(OpcodePlain(cMul));
        // 1/b
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cInv));
        // a<=
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cLessOrEq));
        // *
        outer.push_back(OpcodePlain(cMul));
        // *
        outer.push_back(OpcodePlain(cAdd));
        return outer;
      }
      break;
  }

  // we encountered an unsupported opcode (this should not happen here)
  throw UnsupportedOpcodeException;
}

template<typename Value_t>
void FunctionParserADBase<Value_t>::Commit(const DiffProgramFragment & diff)
{
  // loop over diff and fill in mByteCode and mImmed
  mData->mByteCode.clear();
  mData->mImmed.clear();
  mData->mStackSize = 0;
  int stack_size = 0;

  // compressed immediate data representation
  for (unsigned int i = 0; i < diff.size(); ++i)
  {
    int op_size = OpcodeSize(diff[i]);
    stack_size -= op_size > 0 ? op_size-1 : op_size;
    if (stack_size > mData->mStackSize) mData->mStackSize = stack_size;

    mData->mByteCode.push_back(diff[i].first);

    // handle immediate value
    if (diff[i].first == cImmed)
      mData->mImmed.push_back(diff[i].second);

    // expand function call opcode
    if (diff[i].first == cFCall)
      mData->mByteCode.push_back(diff[i].index);
  }

#ifndef FP_USE_THREAD_SAFE_EVAL
  mData->mStack.resize(mData->mStackSize);
#endif
}

template<typename Value_t>
bool FunctionParserADBase<Value_t>::isZero()
{
  // determine if the program is a single cImmed 0
  return (mData->mByteCode.size() == 1  && mData->mImmed.size() == 1 &&
          mData->mByteCode[0] == cImmed && mData->mImmed[0] == Value_t(0));
}

template<typename Value_t>
void FunctionParserADBase<Value_t>::setZero()
{
  this->ForceDeepCopy();
  mData = this->getParserData();

  // set program to a single cImmed 0
  mData->mByteCode.resize(1);
  mData->mImmed.resize(1);
  mData->mByteCode[0] = cImmed;
  mData->mImmed[0] = Value_t(0);
}

template<typename Value_t>
int FunctionParserADBase<Value_t>::AutoDiff(const std::string& var)
{
  this->ForceDeepCopy();
  mData = this->getParserData();

  // get c string and length of var argument
  const unsigned len = unsigned(var.size());
  const char* name = var.c_str();

  // reset opcode of the variable we diff for
  mVarOp = 0;

  // figure out the opcode number that corresponds to 'var', the variable we diff for
  typename FUNCTIONPARSERTYPES::NamePtrsMap<Value_t> & NamePtrs = mData->mNamePtrs;
  typename FUNCTIONPARSERTYPES::NamePtrsMap<Value_t>::iterator vi;
  for (vi = NamePtrs.begin(); vi != NamePtrs.end(); ++vi)
  {
    if (len == vi->first.nameLength &&
        std::memcmp(name, vi->first.name, len) == 0) {
      mVarOp = vi->second.index;
      break;
    }
  }

  // invalid var argument, variable not found
  if (mVarOp == 0) return 1;

  // get a reference to the stored bytecode
  const std::vector<unsigned>& ByteCode = mData->mByteCode;

  // get a reference to the immediate values
  const std::vector<Value_t>& Immed = mData->mImmed;

  // uncompressed immediate data representation
  // we also expand a few multiopcodes into elementary opcodes
  // to keep the diff rules above simple (cDup must be expanded in this step)
  DiffProgramFragment orig, diff;
  unsigned op;
  unsigned int nImmed = 0;

  try
  {
    for (unsigned int i = 0; i < ByteCode.size(); ++i)
    {
      op = ByteCode[i];
      if (op == cImmed)
        orig.push_back(OpcodeImmediate(Immed[nImmed++]));
      else if (op == cDup)
      {
        // substitute full code for cDup opcodes
        Interval arg = GetArgument(orig);
        orig.insert(orig.end(), arg.first, arg.second);
      }
      else if (op == cFCall)
        orig.push_back(OpcodeFCall(ByteCode[++i]));
      else if (op == cJump)
        throw UnsupportedOpcodeException;
      else if (op == cFetch)
      {
        // get index and advance instruction counter
        throw UnsupportedOpcodeException; // for now

        unsigned index = ByteCode[++i];
        Interval arg = GetArgument(orig, index); // maybe off by one
        orig.insert(orig.end(), arg.first, arg.second);
      }
      else if (op == cSinCos)
      {
        // this instruction puts two values on the stack!
        Interval arg = GetArgument(orig);
        DiffProgramFragment sub(arg.first, arg.second);
        orig.push_back(OpcodePlain(cSin));
        orig.insert(orig.end(), sub.begin(), sub.end());
        orig.push_back(OpcodePlain(cCos));
      }
      else if (op == cCsc)
      {
        orig.push_back(OpcodePlain(cSin));
        orig.push_back(OpcodePlain(cInv));
      }
      else if (op == cSec)
      {
        orig.push_back(OpcodePlain(cSin));
        orig.push_back(OpcodePlain(cInv));
      }
      else if (op == cCot)
      {
        orig.push_back(OpcodePlain(cTan));
        orig.push_back(OpcodePlain(cInv));
      }
  #ifdef FP_SUPPORT_OPTIMIZER
      else if (op == cNop)
        continue;
  #endif
      else
        orig.push_back(OpcodePlain(op));
    }

    diff  = DiffFunction(orig);

    // create compressed program representation
    Commit(diff);
  }
  catch(std::exception &e)
  {
    static bool printed_error = false;
    if (!printed_error && !mSilenceErrors)
    {
      std::cerr << "AutoDiff exception: " << e.what() << " (this message will only be shown once per process)"<< std::endl;
      printed_error = true;
    }
    setZero();
    return 0;
  }

  return -1;
}

template<typename Value_t>
FunctionParserADBase<Value_t>::OpcodeImmediate::OpcodeImmediate(Value_t _second) :
    OpcodePacket(FUNCTIONPARSERTYPES::cImmed, _second, 0)
{
}

template<typename Value_t>
FunctionParserADBase<Value_t>::OpcodeFCall::OpcodeFCall(unsigned _index) :
    OpcodePacket(FUNCTIONPARSERTYPES::cFCall, Value_t(), _index)
{
}

#define FUNCTIONPARSERAD_INSTANTIATE_CLASS(type) \
    template class FunctionParserADBase< type >;

#ifndef FP_DISABLE_DOUBLE_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(double)
#endif

#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(std::complex<double>)
#endif

#ifdef FP_SUPPORT_FLOAT_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(float)
#endif

#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(std::complex<float>)
#endif

#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(long double)
#endif

#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(std::complex<long double>)
#endif
