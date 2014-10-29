#include "fparser_ad.hh"
#include "extrasrc/fpaux.hh"
#include "extrasrc/fptypes.hh"
#include "extrasrc/fpaux.hh"
#include <stdlib.h>

using namespace FUNCTIONPARSERTYPES;

#include <iostream>

#if LIBMESH_HAVE_FPARSER_JIT
#  include <fstream>
#  include <cstdio>
#  include <unistd.h>
#  include <dlfcn.h>
#  include "lib/sha1.h"
#  include <errno.h>
#  include <sys/stat.h>
#endif

template<typename Value_t>
FunctionParserADBase<Value_t>::FunctionParserADBase() :
    FunctionParserBase<Value_t>(),
    mData(this->getParserData()),
    compiledFunction(NULL),
    mSilenceErrors(false),
    mFPlog(mData->mFuncPtrs.size())
{
  this->AddFunction("plog", fp_plog, 2);
}

template<typename Value_t>
FunctionParserADBase<Value_t>::FunctionParserADBase(const FunctionParserADBase& cpy) :
    FunctionParserBase<Value_t>(cpy),
    mData(this->getParserData()),
    compiledFunction(cpy.compiledFunction),
    mSilenceErrors(cpy.mSilenceErrors),
    mFPlog(cpy.mFPlog)
{
}

template<typename Value_t>
Value_t FunctionParserADBase<Value_t>::fp_plog(const Value_t * params)
{
  const Value_t x = params[0];
  const Value_t a = params[1];
  // return x < a ? fp_log(a) + (params[0] - a) / a : fp_log(params[0]);
  // return x < a ? fp_log(a) - Value_t(1.5) + Value_t(2.0)/a * x - Value_t(0.5)/(a*a) * x*x : fp_log(x);
  return x < a ? fp_log(a) + (x-a)/a - (x-a)*(x-a)/(Value_t(2)*a*a) + (x-a)*(x-a)*(x-a)/(Value_t(3)*a*a*a) : fp_log(x);
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
      case cEqual: case cNEqual:
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
      outer.insert(outer.end(), prog_a.begin(), prog_a.end()); // insert twice, as cDup is not supported
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
      outer.push_back(OpcodeImmediate(Value_t(0)));
      return outer;

    // we return those undiffed to keep conditionals intact (for piecewise functions like:  (x<0) * 1 + (x>=0) * (1+x^2) )
    case cEqual:
    case cNEqual:
    case cLess:
    case cLessOrEq:
    case cGreater:
    case cGreaterOrEq:
      outer = prog_a;
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.push_back(OpcodePlain(op));
      return outer;

    case cFCall:
      if (findex == mFPlog)
      {
        // we assume that the second argument to plog is constant for now (TODO?).
        /*
        prog_da = DiffFunction(prog_a);
        // da (inner derivative)
        outer = prog_da;
        // 1/a
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.push_back(OpcodePlain(cInv));
        // a>b
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cGreater));
        // *
        outer.push_back(OpcodePlain(cMul));
        // 1/b
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cInv));
        // a<=b
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cLessOrEq));
        // *
        outer.push_back(OpcodePlain(cMul));
        // +
        outer.push_back(OpcodePlain(cAdd));
        // multiply by inner derivative da
        outer.push_back(OpcodePlain(cMul));
        */

        /*
        // x<e ? (2/e - 1/(e*e)*x) : 1/x
        prog_da = DiffFunction(prog_a);
        // da (inner derivative)
        outer = prog_da;
        // 1/a
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.push_back(OpcodePlain(cInv));
        // a>b
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cGreater));
        // *
        outer.push_back(OpcodePlain(cMul));
        // 2/b
        outer.push_back(OpcodeImmediate(Value_t(2)));
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cDiv));
        /// 1/(e*e) = e^-2
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodeImmediate(Value_t(-2)));
        outer.push_back(OpcodePlain(cPow));
        // * x
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.push_back(OpcodePlain(cMul));
        // -
        outer.push_back(OpcodePlain(cSub));
        // a<=b
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cLessOrEq));
        // *
        outer.push_back(OpcodePlain(cMul));
        // +
        outer.push_back(OpcodePlain(cAdd));
        // multiply by inner derivative da
        outer.push_back(OpcodePlain(cMul));
        */

        // x<e ? (1/e - 1/(e*e)*(x-e) + 1/e^3*(x-e)^2) : 1/x
        prog_da = DiffFunction(prog_a);
        // da (inner derivative)
        outer = prog_da;
        // 1/a
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        // actually do 1/(a + (a==0)). This avoids a 1/0 and a=0 should always make this the false-branch!
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.push_back(OpcodeImmediate(Value_t(0)));
        outer.push_back(OpcodePlain(cEqual));
        outer.push_back(OpcodePlain(cAdd));

        outer.push_back(OpcodePlain(cInv));
        // a>b
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cGreater));
        // *
        outer.push_back(OpcodePlain(cMul));
        // 1/b
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cInv));
        /// 1/(b*b) = b^-2
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodeImmediate(Value_t(-2)));
        outer.push_back(OpcodePlain(cPow));
        // * (a-b)
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cSub));
        outer.push_back(OpcodePlain(cMul));
        // -
        outer.push_back(OpcodePlain(cSub));
        // 1/(e*e*e) = e^-3
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodeImmediate(Value_t(-3)));
        outer.push_back(OpcodePlain(cPow));
        // * (x-e)^2
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cSub));
        outer.push_back(OpcodeImmediate(Value_t(2)));
        outer.push_back(OpcodePlain(cPow));
        outer.push_back(OpcodePlain(cMul));
        // +
        outer.push_back(OpcodePlain(cAdd));
        // a<=b
        outer.insert(outer.end(), prog_a.begin(), prog_a.end());
        outer.insert(outer.end(), prog_b.begin(), prog_b.end());
        outer.push_back(OpcodePlain(cLessOrEq));
        // *
        outer.push_back(OpcodePlain(cMul));
        // +
        outer.push_back(OpcodePlain(cAdd));
        // multiply by inner derivative da
        outer.push_back(OpcodePlain(cMul));

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
    // mData->mStackSize is unsigned, stack_size should be an int
    // since we subtract from it, so cast mStackSize to int for
    // comparison to avoid compiler warnings.
    if (stack_size > static_cast<int>(mData->mStackSize))
      mData->mStackSize = stack_size;

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
typename FunctionParserADBase<Value_t>::DiffProgramFragment
FunctionParserADBase<Value_t>::Expand()
{
  // get a reference to the stored bytecode
  const std::vector<unsigned>& ByteCode = mData->mByteCode;

  // get a reference to the immediate values
  const std::vector<Value_t>& Immed = mData->mImmed;

  DiffProgramFragment orig;
  unsigned op;
  unsigned int nImmed = 0;

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
      throw UnsupportedOpcodeException;
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
      orig.push_back(OpcodePlain(cCos));
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

  return orig;
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

  // uncompressed immediate data representation
  // we also expand a few multiopcodes into elementary opcodes
  // to keep the diff rules above simple (cDup must be expanded in this step)
  DiffProgramFragment orig, diff;

  try
  {
    // expand the internal byte code into a more convenient form
    orig = Expand();

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
bool FunctionParserADBase<Value_t>::JITCompile(bool)
{
  // JIT compile attempted for an unsupported value type
  return false;
}

#if LIBMESH_HAVE_FPARSER_JIT
// JIT compile for supported types
template<>
bool FunctionParserADBase<double>::JITCompile(bool cacheFunction) { return JITCompileHelper("double", cacheFunction); }
template<>
bool FunctionParserADBase<float>::JITCompile(bool cacheFunction) { return JITCompileHelper("float", cacheFunction); }
template<>
bool FunctionParserADBase<long double>::JITCompile(bool cacheFunction) { return JITCompileHelper("long double", cacheFunction); }

template<typename Value_t>
Value_t FunctionParserADBase<Value_t>::Eval(const Value_t* Vars)
{
  if (compiledFunction == NULL)
    return FunctionParserBase<Value_t>::Eval(Vars);
  else
    return (*compiledFunction)(Vars, &(mData->mImmed[0]), Epsilon<Value_t>::value);
}

template<typename Value_t>
bool FunctionParserADBase<Value_t>::JITCompileHelper(const std::string & Value_t_name, bool cacheFunction)
{
  // get a reference to the stored bytecode
  const std::vector<unsigned>& ByteCode = mData->mByteCode;

  // generate a sha1 hash of the current program and the Value type name
  SHA1 *sha1 = new SHA1();
  char result[41];
  sha1->addBytes((char*) &(ByteCode[0]), ByteCode.size() * sizeof(unsigned));
  sha1->addBytes(Value_t_name.c_str(), Value_t_name.size());
  unsigned char* digest = sha1->getDigest();
  for (unsigned int i = 0; i<20; ++i)
    sprintf(&(result[i*2]), "%02x", digest[i]);
  delete sha1;

  // function name
  std::string fnname = "f_";
  fnname += result;

  // cache file name
  std::string jitdir = ".jitcache";
  std::string libname = jitdir + "/" + result + ".so";

  // attempt to open the cache file
  void * lib;
  if (cacheFunction)
  {
    lib = dlopen(libname.c_str(), RTLD_NOW);
    if (lib != NULL) {
      // fetch function pointer
      *(void **) (&compiledFunction) = dlsym(lib, fnname.c_str());
      if (dlerror() == NULL)  {
        // success
        return true;
      }
    }
  }
  // opening the cached file did nbot work. (re)build it.

  // tmp filenames
  char ccname[] = "./tmp_jit_XXXXXX";
  if (mkstemp(ccname) == -1)
  {
    std::cerr << "Error creating JIT tmp file " << ccname << ".\n";
    return false;
  }
  char object[] = "./tmp_jit_XXXXXX";
  if (mkstemp(object) == -1)
  {
    std::cerr << "Error creating JIT tmp file " << object << ".\n";
    std::remove(ccname);
    return false;
  }

  int status;
  std::vector<bool> jumpTarget(ByteCode.size());
  for (unsigned int i = 0; i < ByteCode.size(); ++i) jumpTarget[i] = false;

  // determine all jump targets in the current program
  unsigned long ip;
  for (unsigned int i = 0; i < ByteCode.size(); ++i)
    switch(ByteCode[i])
    {
      case cJump: case cIf: case cAbsIf:
        ip = ByteCode[i+1] + 1;
        if (ip < ByteCode.size())
          jumpTarget[ip] = true;
#ifdef FP_SUPPORT_OPTIMIZER
      case cPopNMov:
#endif
        i += 2;
        break;

      case cFetch: case cFCall:  case cPCall:
        i++;
        break;
    }

  // stack pointer recovery at jump targets
  std::vector<int> stackAtTarget(ByteCode.size());
  for (unsigned int i = 0; i < ByteCode.size(); ++i) stackAtTarget[i] = -2;

  // save source
  std::ofstream ccfile;
  ccfile.open(ccname);
  ccfile << "#define _USE_MATH_DEFINES\n";
  ccfile << "#include <cmath>\n";
  ccfile << "extern \"C\" " << Value_t_name << ' '
         << fnname << "(const " << Value_t_name << " *params, const " << Value_t_name << " *immed, const " << Value_t_name << " eps) {\n";
  ccfile << Value_t_name << " r, s[" << mData->mStackSize << "];\n";
  int nImmed = 0, sp = -1, op;
  for (unsigned int i = 0; i < ByteCode.size(); ++i)
  {
    // place a label?
    if (jumpTarget[i])
    {
      ccfile << "l" << i << ":\n";
      if (stackAtTarget[i] > -2)
        sp = stackAtTarget[i];
    }

    // translate opcode into C code
    switch (op = ByteCode[i])
    {
      case cImmed:
        ++sp; ccfile << "s[" << sp << "] = immed[" << nImmed++ << "];\n"; break;
      case cAdd:
        --sp; ccfile << "s[" << sp << "] += s[" << (sp+1) << "];\n"; break;
      case cSub:
        --sp; ccfile << "s[" << sp << "] -= s[" << (sp+1) << "];\n"; break;
      case cRSub:
        --sp; ccfile << "s[" << sp << "] = s[" << (sp+1) << "] - s[" << sp << "];\n"; break;
      case cMul:
        --sp; ccfile << "s[" << sp << "] *= s[" << (sp+1) << "];\n"; break;
      case cDiv:
        --sp; ccfile << "s[" << sp << "] /= s[" << (sp+1) << "];\n"; break;
      case cMod:
        --sp; ccfile << "s[" << sp << "] = std::fmod(s[" << sp << "], s[" << (sp+1) << "]);\n"; break;
      case cRDiv:
        --sp; ccfile << "s[" << sp << "] = s[" << (sp+1) << "] / s[" << sp << "];\n"; break;

      case cSin:
        ccfile << "s[" << sp << "] = std::sin(s[" << sp << "]);\n"; break;
      case cCos:
        ccfile << "s[" << sp << "] = std::cos(s[" << sp << "]);\n"; break;
      case cTan:
        ccfile << "s[" << sp << "] = std::tan(s[" << sp << "]);\n"; break;
      case cSinh:
        ccfile << "s[" << sp << "] = std::sinh(s[" << sp << "]);\n"; break;
      case cCosh:
        ccfile << "s[" << sp << "] = std::cosh(s[" << sp << "]);\n"; break;
      case cTanh:
        ccfile << "s[" << sp << "] = std::tanh(s[" << sp << "]);\n"; break;
      // TODO: div by zero -> mData->mEvalErrorType=1; return Value_t(0);
      case cCsc:
        ccfile << "s[" << sp << "] = 1.0/std::sin(s[" << sp << "]);\n"; break;
      case cSec:
        ccfile << "s[" << sp << "] = 1.0/std::cos(s[" << sp << "]);\n"; break;
      case cCot:
        ccfile << "s[" << sp << "] = 1.0/std::tan(s[" << sp << "]);\n"; break;
      case cSinCos:
        ccfile << "s[" << (sp+1) << "] = std::cos(s[" << sp << "]);\n";
        ccfile << "s[" << sp << "] = std::sin(s[" << sp << "]);\n";
        ++sp;
        break;
      case cSinhCosh:
        ccfile << "s[" << (sp+1) << "] = std::cosh(s[" << sp << "]);\n";
        ccfile << "s[" << sp << "] = std::sinh(s[" << sp << "]);\n";
        ++sp;
        break;
      case cAsin:
        ccfile << "s[" << sp << "] = std::asin(s[" << sp << "]);\n"; break;
      case cAcos:
        ccfile << "s[" << sp << "] = std::acos(s[" << sp << "]);\n"; break;
      case cAsinh:
        ccfile << "s[" << sp << "] = std::asinh(s[" << sp << "]);\n"; break;
      case cAcosh:
        ccfile << "s[" << sp << "] = std::acosh(s[" << sp << "]);\n"; break;
      case cAtan:
        ccfile << "s[" << sp << "] = std::atan(s[" << sp << "]);\n"; break;
      case cAtanh:
        ccfile << "s[" << sp << "] = std::atanh(s[" << sp << "]);\n"; break;
      case cAtan2:
        --sp; ccfile << "s[" << sp << "] = std::atan2(s[" << sp << "], s[" << (sp+1) << "]);\n"; break;
      case cHypot:
        --sp; ccfile << "s[" << sp << "] = std::sqrt(s[" << sp << "]*s[" << sp << "] + s[" << (sp+1) << "]*s[" << (sp+1) << "]);\n"; break;

      case cAbs:
        ccfile << "s[" << sp << "] = std::abs(s[" << sp << "]);\n"; break;
      case cMax:
        --sp; ccfile << "s[" << sp << "] = s[" << sp << "] > s[" << (sp+1) << "] ? s[" << sp << "] : s[" << (sp+1) << "];\n"; break;
      case cMin:
        --sp; ccfile << "s[" << sp << "] = s[" << sp << "] < s[" << (sp+1) << "] ? s[" << sp << "] : s[" << (sp+1) << "];\n"; break;
      case cTrunc:
        ccfile << "s[" << sp << "] = s[" << sp << "] < 0 ? std::ceil(s[" << sp << "]) : std::floor(s[" << sp << "]);\n"; break;
      case cCeil:
        ccfile << "s[" << sp << "] = std::ceil(s[" << sp << "]);\n"; break;
      case cFloor:
        ccfile << "s[" << sp << "] = std::floor(s[" << sp << "]);\n"; break;
      case cInt:
        ccfile << "s[" << sp << "] = s[" << sp << "] < 0 ? std::ceil(s[" << sp << "] - 0.5) : std::floor(s[" << sp << "] + 0.5);\n"; break;

      case cEqual:
        //--sp; ccfile << "s[" << sp << "] = s[" << sp << "] == s[" << (sp+1) << "];\n"; break;
        --sp; ccfile << "s[" << sp << "] = std::abs(s[" << sp << "] - s[" << (sp+1) << "]) <= eps;\n"; break;
      case cNEqual:
        //--sp; ccfile << "s[" << sp << "] = s[" << sp << "] != s[" << (sp+1) << "];\n"; break;
        --sp; ccfile << "s[" << sp << "] = std::abs(s[" << sp << "] - s[" << (sp+1) << "]) > eps;\n"; break;
      case cLess:
        --sp; ccfile << "s[" << sp << "] = s[" << sp << "] < (s[" << (sp+1) << "] - eps);\n"; break;
      case cLessOrEq:
        --sp; ccfile << "s[" << sp << "] = s[" << sp << "] <= (s[" << (sp+1) << "] + eps);\n"; break;
      case cGreater:
        --sp; ccfile << "s[" << sp << "] = (s[" << sp << "] - eps) > s[" << (sp+1) << "];\n"; break;
      case cGreaterOrEq:
        --sp; ccfile << "s[" << sp << "] = (s[" << sp << "] + eps) >= s[" << (sp+1) << "];\n"; break;
      case cNot:
        ccfile << "s[" << sp << "] = std::abs(s[" << sp << "]) < 0.5;\n"; break;
      case cNotNot:
        ccfile << "s[" << sp << "] = std::abs(s[" << sp << "]) >= 0.5;\n"; break;
      case cAbsNot:
        ccfile << "s[" << sp << "] = s[" << sp << "] < 0.5;\n"; break;
      case cAbsNotNot:
        ccfile << "s[" << sp << "] = s[" << sp << "] >= 0.5;\n"; break;
      case cOr:
        --sp; ccfile << "s[" << sp << "] = (std::abs(s[" << sp << "]) >= 0.5) || (std::abs(s[" << (sp+1) << "]) >= 0.5);\n"; break;
      case cAbsOr:
        --sp; ccfile << "s[" << sp << "] = (s[" << sp << "] >= 0.5) || (s[" << (sp+1) << "] >= 0.5);\n"; break;
      case cAnd:
        --sp; ccfile << "s[" << sp << "] = (std::abs(s[" << sp << "]) >= 0.5) && (std::abs(s[" << (sp+1) << "]) >= 0.5);\n"; break;
      case cAbsAnd:
        --sp; ccfile << "s[" << sp << "] = (s[" << sp << "] >= 0.5) && (s[" << (sp+1) << "] >= 0.5);\n"; break;

      case cLog:
        ccfile << "s[" << sp << "] = std::log(s[" << sp << "]);\n"; break;
      case cLog2:
#ifdef FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS
        ccfile << "s[" << (sp-1) << "] = std::log2(s[" << (sp-1) << "]);\n";
#else
        ccfile << "s[" << sp << "] = std::log(s[" << sp << "])/log(2.0);\n";
#endif
        break;
      case cLog10:
        ccfile << "s[" << sp << "] = std::log10(s[" << sp << "]);\n"; break;

      case cNeg:
        ccfile << "s[" << sp << "] = -s[" << sp << "];\n"; break;
      case cInv:
        ccfile << "s[" << sp << "] = 1.0/s[" << sp << "];\n"; break;
      case cDeg:
        ccfile << "s[" << sp << "] *= 180.0/M_PI;\n"; break;
      case cRad:
        ccfile << "s[" << sp << "] /= 180.0/M_PI;\n"; break;

      case cFetch:
        ++sp; ccfile << "s[" << sp << "] = s[" << ByteCode[++i] << "];\n"; break;
      case cDup:
        ++sp; ccfile << "s[" << sp << "] = s[" << (sp-1) << "];\n"; break;

      case cFCall:
      {
        unsigned function = ByteCode[++i];
        if (function == mFPlog)
        {
          // --sp; ccfile << "s[" << sp << "] = s[" << sp << "] < s[" << (sp+1) << "] ? std::log(s[" << (sp+1) << "]) + (s[" << sp << "] - s[" << (sp+1) << "]) / s[" << (sp+1) << "] : std::log(s[" << sp << "]);\n";
          // --sp; ccfile << "s[" << sp << "] = s[" << sp << "] < s[" << (sp+1) << "] ? std::log(s[" << (sp+1) << "]) - 1.5 + 2.0/s[" << (sp+1) << "] * s[" << sp << "] - 0.5/(s[" << (sp+1) << "]*s[" << (sp+1) << "]) * s[" << sp << "]*s[" << sp << "] : std::log(s[" << sp << "]);\n";
          --sp; ccfile << "s[" << sp << "] = s[" << sp << "] < s[" << (sp+1) << "] ? std::log(s[" << (sp+1) << "])  +  (s[" << sp << "]-s[" << (sp+1) << "])/s[" << (sp+1) << "]  -  std::pow((s[" << sp << "]-s[" << (sp+1) << "])/s[" << (sp+1) << "],2.0)/2.0  +  std::pow((s[" << sp << "]-s[" << (sp+1) << "])/s[" << (sp+1) << "],3.0)/3.0 : std::log(s[" << sp << "]);\n";
        }
        else
        {
          std::cerr << "Function call not supported by JIT.\n";
          ccfile.close();
          std::remove(ccname);
          return false;
        }
        break;
      }

#ifdef FP_SUPPORT_OPTIMIZER
      case cPopNMov:
      {
        int dst = ByteCode[++i],
            src = ByteCode[++i];
        ccfile << "s[" << dst << "] = s[" << src << "];\n";
        sp = dst;
        break;
      }
      case cLog2by:
        --sp; ccfile << "s[" << sp << "] = std::log(s[" << sp << "])/log(2.0) * s[" << (sp+1) << "];\n"; break;
      case cNop:
        break;
#endif

      case cSqr:
        ccfile << "s[" << sp << "] *= s[" << sp << "];\n"; break;
      case cSqrt:
        ccfile << "s[" << sp << "] = std::sqrt(s[" << sp << "]);\n"; break;
      case cRSqrt:
        ccfile << "s[" << sp << "] = std::pow(s[" << sp << "], " << Value_t_name << "(-0.5));\n"; break;
      case cPow:
        --sp; ccfile << "s[" << sp << "] = std::pow(s[" << sp << "], s[" << (sp+1) << "]);\n"; break;
      case cExp:
        ccfile << "s[" << sp << "] = std::exp(s[" << sp << "]);\n"; break;
      case cExp2:
        ccfile << "s[" << sp << "] = std::pow(" << Value_t_name << "(2.0), s[" << sp << "]);\n"; break;
      case cCbrt:
#ifdef FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS
        ccfile << "s[" << sp << "] = std::cbrt(s[" << sp << "]);\n"; break;
#else
        ccfile << "s[" << sp << "] = s[" << sp << "] == 0 ? 0 : (s[" << sp << "] > 0 ? std::exp(std::log(s[" << sp << "])/3.0) : -std::exp(std::log(-s[" << sp << "])/3.0));\n"; break;
#endif

      case cJump:
      case cIf:
      case cAbsIf:
      {
        unsigned long ip = ByteCode[++i] + 1;

        if (op == cIf)
          ccfile << "if (std::abs(s[" << sp-- << "]) < 0.5) ";
        if (op == cAbsIf)
          ccfile << "if (s[" << sp-- << "] < 0.5) ";

        if (ip >= ByteCode.size())
          ccfile << "return s[" << sp << "];\n";
        else
        {
          ccfile << "goto l" << ip << ";\n";
          stackAtTarget[ip] = sp;
        }

        ++i;
        break;
      }

      default:
        if ( op>= VarBegin)
        {
          // load variable
          ++sp; ccfile << "s[" << sp << "] = params[" << (op-VarBegin) << "];\n";
        }
        else
        {
          std::cerr << "Opcode not supported by JIT.\n";
          ccfile.close();
          std::remove(ccname);
          return false;
        }
    }
  }
  ccfile << "return s[" << sp << "]; }\n";
  ccfile.close();

  // add a .cc extension to the source (needed by the compiler)
  std::string ccname_cc = ccname;
  ccname_cc += ".cc";
  status = std::rename(ccname, ccname_cc.c_str());
  if (status != 0)
  {
    std::cerr << "Unable to rename JIT source code file\n";
    std::remove(ccname);
    return false;
  }

  // run compiler
  std::string command = FPARSER_JIT_COMPILER" -O2 -shared -rdynamic -fPIC ";
  command += ccname_cc + " -o " + object;
  status = system(command.c_str());
  std::remove(ccname_cc.c_str());
  if (status != 0) {
    std::cerr << "JIT compile failed.\n";
    return false;
  }

  // add a .so extension to the object (needed by dlopen on mac)
  std::string object_so = object;
  object_so += ".so";
  status = std::rename(object, object_so.c_str());
  if (status != 0)
  {
    std::cerr << "Unable to rename JIT compiled function object\n";
    std::remove(object);
    return false;
  }

  // load compiled object
  lib = dlopen(object_so.c_str(), RTLD_NOW);
  if (lib == NULL) {
    std::cerr << "JIT object load failed.\n";
    std::remove(object_so.c_str());
    return false;
  }

  // fetch function pointer
  *(void **) (&compiledFunction) = dlsym(lib, fnname.c_str());
  char * error;
  if ((error = dlerror()) != NULL)  {
    std::cerr << "Error binding JIT compiled function\n" << error << '\n';
    compiledFunction = NULL;
    std::remove(object_so.c_str());
    return false;
  }

  // clear evalerror (this will not get set again by the JIT code)
  mData->mEvalErrorType = 0;

  // rename successfully compiled obj to cache file
  if (cacheFunction && (mkdir(jitdir.c_str(), 0700) == 0 || errno == EEXIST)) {
    // the directory was either successfuly created, or it already exists
    status = std::rename(object_so.c_str(), libname.c_str());
    if (status == 0) return true;
  }

  std::remove(object_so.c_str());
  return true;
}
#endif

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
