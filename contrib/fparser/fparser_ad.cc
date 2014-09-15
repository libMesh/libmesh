#include "fparser_ad.hh"
#include "extrasrc/fptypes.hh"
#include <stdlib.h>

using namespace FUNCTIONPARSERTYPES;

#include <iostream>

#ifdef LIBMESH_HAVE_FPARSER_JIT
#  include <fstream>
#  include <cstdio>
#  include <cstdlib>
#  include <dlfcn.h>
#  include "lib/sha1.h"
#endif

template<typename Value_t>
FunctionParserADBase<Value_t>::FunctionParserADBase() :
    FunctionParserBase<Value_t>(),
    mData(this->getParserData()),
    compiledFunction(NULL),
    mSilenceErrors(false)
{
}

template<typename Value_t>
FunctionParserADBase<Value_t>::FunctionParserADBase(const FunctionParserADBase& cpy) :
    FunctionParserBase<Value_t>(cpy),
    mData(this->getParserData()),
    compiledFunction(cpy.compiledFunction),
    mSilenceErrors(cpy.mSilenceErrors)
{
}

template<typename Value_t>
int FunctionParserADBase<Value_t>::OpcodeSize(unsigned op)
{
  if (int(op) >= VarBegin || op == cImmed) {
    return -1;
  }
  else {
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
    int this_size = OpcodeSize((*ip).first);
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

  // current opcode
  unsigned op = orig.back().first;

  // size of the current end opcode
  int op_size = OpcodeSize(op);

  // variable or immediate value
  if (op_size == -1)
  {
    if (op == cImmed)
      outer.push_back(OpcodeDataPair(cImmed, Value_t(0)));
    else
    {
      if (op == mVarOp)
        outer.push_back(OpcodeDataPair(cImmed, Value_t(1)));
      else
        outer.push_back(OpcodeDataPair(cImmed, Value_t(0)));
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
      outer.push_back(OpcodeDataPair(op, 0));
      return outer;

    case cMul:
      // a*db + da*b
      prog_da = DiffFunction(prog_a);
      prog_db = DiffFunction(prog_b);
      outer = prog_a;
      outer.insert(outer.end(), prog_db.begin(), prog_db.end());
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.push_back(OpcodeDataPair(cAdd, 0));
      return outer;

    case cSqr:
      // a*da*2
      prog_da = DiffFunction(prog_a);
      outer = prog_a;
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.push_back(OpcodeDataPair(cImmed, Value_t(2)));
      outer.push_back(OpcodeDataPair(cMul, 0));
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
      outer.push_back(OpcodeDataPair(cDiv, 0));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.insert(outer.end(), prog_db.begin(), prog_db.end());
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.push_back(OpcodeDataPair(cSqr, 0));
      outer.push_back(OpcodeDataPair(cDiv, 0));
      outer.push_back(OpcodeDataPair(cSub, 0));
      return outer;

    case cLog2:
    case cLog10:
      // *1/ln(base)
    case cLog:
      // da/a
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cDiv, 0));
      if (op == cLog2) {
        outer.push_back(OpcodeDataPair(cImmed, std::log(Value_t(2))));
        outer.push_back(OpcodeDataPair(cDiv, 0));
      }
      else if (op == cLog10) {
        outer.push_back(OpcodeDataPair(cImmed, std::log(Value_t(10))));
        outer.push_back(OpcodeDataPair(cDiv, 0));
      }
      return outer;

    case cExp:
    case cExp2:
      // da*exp(a)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(op, 0));
      outer.push_back(OpcodeDataPair(cMul, 0));
      if (op == cExp2) {
        outer.push_back(OpcodeDataPair(cImmed, std::log(Value_t(2))));
        outer.push_back(OpcodeDataPair(cMul, 0));
      }
      return outer;

    case cNeg:
      // -da
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodeDataPair(cNeg, 0));
      return outer;

    case cInv:
      // -da/a^2
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cSqr, 0));
      outer.push_back(OpcodeDataPair(cDiv, 0));
      outer.push_back(OpcodeDataPair(cNeg, 0));
      return outer;

    case cPow:
      // a**b * ( db*log(a) + b*da/a)
      prog_da = DiffFunction(prog_a);
      prog_db = DiffFunction(prog_b);
      outer = prog_a;
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.push_back(OpcodeDataPair(cPow, 0));
      outer.insert(outer.end(), prog_db.begin(), prog_db.end());
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cLog, 0));
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cDiv, 0));
      outer.push_back(OpcodeDataPair(cAdd, 0));
      outer.push_back(OpcodeDataPair(cMul, 0));
      return outer;

    case cSin :
      // da*cos(a)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cCos, 0));
      outer.push_back(OpcodeDataPair(cMul, 0));
      return outer;

    case cCos :
      // -da*sin(a)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cSin, 0));
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.push_back(OpcodeDataPair(cNeg, 0));
      return outer;

    case cTan :
      // (tan(a)^2+1)*da
      prog_da = DiffFunction(prog_a);
      outer = prog_a;
      outer.push_back(OpcodeDataPair(cTan, 0));
      outer.push_back(OpcodeDataPair(cSqr, 0));
      outer.push_back(OpcodeDataPair(cImmed, Value_t(1)));
      outer.push_back(OpcodeDataPair(cAdd, 0));
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.push_back(OpcodeDataPair(cMul, 0));
      return outer;

    case cAtan :
      // da/(a^2+1)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cSqr, 0));
      outer.push_back(OpcodeDataPair(cImmed, Value_t(1)));
      outer.push_back(OpcodeDataPair(cAdd, 0));
      outer.push_back(OpcodeDataPair(cDiv, 0));
      return outer;

    case cAsin:
      // da/sqrt(1-a^2)
    case cAcos:
      // -da/sqrt(1-a^2)
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodeDataPair(cImmed, Value_t(1)));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cSqr, 0));
      outer.push_back(OpcodeDataPair(cSub, 0));
      outer.push_back(OpcodeDataPair(cSqrt, 0));
      outer.push_back(OpcodeDataPair(cDiv, 0));
      if (op == cAcos)
        outer.push_back(OpcodeDataPair(cNeg, 0));
      return outer;

    case cSqrt :
      // da/(2*sqrt(a))
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodeDataPair(cImmed, Value_t(2)));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cSqrt, 0));
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.push_back(OpcodeDataPair(cDiv, 0));
      return outer;

    case cRSqrt :
      // -da/(2*a^(3/2)))
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodeDataPair(cImmed, Value_t(2)));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cImmed, Value_t(3)));
      outer.push_back(OpcodeDataPair(cPow, 0));
      outer.push_back(OpcodeDataPair(cSqrt, 0));
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.push_back(OpcodeDataPair(cDiv, 0));
      outer.push_back(OpcodeDataPair(cNeg, 0));
      return outer;

    case cCbrt :
      // da/(3*a^(2/3))
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.push_back(OpcodeDataPair(cImmed, Value_t(3)));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cCbrt, 0));
      outer.push_back(OpcodeDataPair(cSqr, 0));
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.push_back(OpcodeDataPair(cDiv, 0));
      return outer;

    case cAbs:
      // da*a/|a|
      prog_da = DiffFunction(prog_a);
      outer = prog_da;
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.push_back(OpcodeDataPair(cDup, 0));
      outer.push_back(OpcodeDataPair(cAbs, 0));
      outer.push_back(OpcodeDataPair(cDiv, 0));
      outer.push_back(OpcodeDataPair(cMul, 0));
      return outer;

    case cHypot:
      // (a*da+b*db)/hypot(a,b)
      prog_da = DiffFunction(prog_a);
      prog_db = DiffFunction(prog_b);
      outer = prog_a;
      outer.insert(outer.end(), prog_da.begin(), prog_da.end());
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.insert(outer.end(), prog_db.begin(), prog_db.end());
      outer.push_back(OpcodeDataPair(cMul, 0));
      outer.push_back(OpcodeDataPair(cAdd, 0));
      outer.insert(outer.end(), prog_a.begin(), prog_a.end());
      outer.insert(outer.end(), prog_b.begin(), prog_b.end());
      outer.push_back(OpcodeDataPair(cHypot, 0));
      outer.push_back(OpcodeDataPair(cDiv, 0));
      return outer;

    case cInt:
    case cFloor:
    case cCeil:
    case cTrunc:
      // no idea why anyone would like to diff those
      // (I'll pretend the discontinuities don't exist -
      // FP could not deal with them in a useful way in any case)
      outer.push_back(OpcodeDataPair(cImmed, Value_t(0)));
      return outer;
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
    int op = diff[i].first;
    int op_size = OpcodeSize(op);
    stack_size -= op_size > 0 ? op_size-1 : op_size;
    // mData->mStackSize is unsigned, stack_size should be an int
    // since we subtract from it, so cast mStackSize to int for
    // comparison to avoid compiler warnings.
    if (stack_size > static_cast<int>(mData->mStackSize))
      mData->mStackSize = stack_size;

    mData->mByteCode.push_back(op);
    if (diff[i].first == cImmed)
      mData->mImmed.push_back(diff[i].second);
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
      orig.push_back(OpcodeDataPair(op, Immed[nImmed++]));
    else if (op == cDup)
    {
      // substitute full code for cDup opcodes
      Interval arg = GetArgument(orig);
      orig.insert(orig.end(), arg.first, arg.second);
    }
    else if (op == cJump)
      throw UnsupportedOpcodeException;
    else if (op == cSinCos)
    {
      // this instruction puts two values on the stack!
      Interval arg = GetArgument(orig);
      DiffProgramFragment sub(arg.first, arg.second);
      orig.push_back(OpcodeDataPair(cSin, 0.0));
      orig.insert(orig.end(), sub.begin(), sub.end());
      orig.push_back(OpcodeDataPair(cCos, 0.0));
    }
    else if (op == cCsc)
    {
      orig.push_back(OpcodeDataPair(cSin, 0.0));
      orig.push_back(OpcodeDataPair(cInv, 0.0));
    }
    else if (op == cSec)
    {
      orig.push_back(OpcodeDataPair(cSin, 0.0));
      orig.push_back(OpcodeDataPair(cInv, 0.0));
    }
    else if (op == cCot)
    {
      orig.push_back(OpcodeDataPair(cTan, 0.0));
      orig.push_back(OpcodeDataPair(cInv, 0.0));
    }
#ifdef FP_SUPPORT_OPTIMIZER
    else if (op == cNop)
      continue;
#endif
    else
      orig.push_back(OpcodeDataPair(op, 0.0));
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

#ifdef LIBMESH_HAVE_FPARSER_JIT
template<typename Value_t>
Value_t FunctionParserADBase<Value_t>::Eval(const Value_t* Vars)
{
  return FunctionParserBase<Value_t>::Eval(Vars);
}

template<>
double FunctionParserADBase<double>::Eval(const double* Vars)
{
  if (compiledFunction == NULL)
    return FunctionParserBase<double>::Eval(Vars);
  else
    return (*compiledFunction)(Vars, &(mData->mImmed[0]));
}

// currently JIT is just implemented for double types (this could easily be extended)
template<>
bool FunctionParserADBase<double>::JITCompile(bool cacheFunction)
{
  DiffProgramFragment code = Expand();

  // generate a sha1 hash of the current program
  SHA1 *sha1 = new SHA1();
  char result[41];
  sha1->addBytes((char*) &(mData->mByteCode[0]), mData->mByteCode.size() * sizeof(unsigned));
  unsigned char* digest = sha1->getDigest();
  for (unsigned int i = 0; i<20; ++i)
    sprintf(&(result[i*2]), "%02x", digest[i]);
  delete sha1;

  // function name
  std::string fnname = "f_";
  fnname += result;

  // cache file name
  std::string jitprefix = "./.jit";
  std::string libname = jitprefix + result + ".so";

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
  char ccname[] = "tmp_jit_XXXXXX.cc";
  if (mkstemps(ccname, 3) == -1)
  {
    std::cerr << "Error creating JIT tmp file " << ccname << ".\n";
    return false;
  }
  char object[] = "tmp_jit_XXXXXX.so";
  if (mkstemps(object, 3) == -1)
  {
    std::cerr << "Error creating JIT tmp file " << object << ".\n";
    std::remove(ccname);
    return false;
  }

  int status;

  // save source
  std::ofstream ccfile;
  ccfile.open(ccname);
  ccfile << "#include <cmath>\n";
  ccfile << "extern \"C\" double " << fnname << "(double *params, double *immed) {\ndouble r, s[" << mData->mStackSize << "] ;\n";
  int nImmed = 0, sp = 0;
  for (unsigned int i = 0; i < code.size(); ++i)
  {
    int op = code[i].first;
    int op_size = OpcodeSize(op);
    if ( op>= VarBegin)
    {
      ccfile << "s[" << sp << "] = params[" << (op-VarBegin) << "];\n";
    }
    else
      switch (op)
      {
        case cImmed:
          ccfile << "s[" << sp << "] = immed[" << nImmed << "]; // " << code[i].second << ',' << mData->mImmed[nImmed] << '\n';
          nImmed++;
          break;
        case cAdd:
          ccfile << "s[" << (sp-2) << "] += s[" << (sp-1) << "];\n";
          break;
        case cSub: // Check if this is the right way around
          ccfile << "s[" << (sp-2) << "] -= s[" << (sp-1) << "];\n";
          break;
        case cRSub: // Check if this is the right way around
          ccfile << "s[" << (sp-2) << "] = s[" << (sp-1) << "] - s[" << (sp-2) << "];\n";
          break;
        case cMul:
          ccfile << "s[" << (sp-2) << "] *= s[" << (sp-1) << "];\n";
          break;
        case cDiv: // Check if this is the right way around
          ccfile << "s[" << (sp-1) << "] /= s[" << sp << "];\n";
          break;
        case cRDiv: // Check if this is the right way around
          ccfile << "s[" << (sp-1) << "] = s[" << sp << "] / s[" << (sp-1) << "];\n";
          break;

        case cSin:
          ccfile << "s[" << (sp-1) << "] = std::sin(s[" << (sp-1) << "]);\n";
          break;
        case cCos:
          ccfile << "s[" << (sp-1) << "] = std::cos(s[" << (sp-1) << "]);\n";
          break;
        case cTan:
          ccfile << "s[" << (sp-1) << "] = std::tan(s[" << (sp-1) << "]);\n";
          break;

        case cAbs:
          ccfile << "s[" << (sp-1) << "] = std::abs(s[" << (sp-1) << "]);\n";
          break;

        case cLog:
          ccfile << "s[" << (sp-1) << "] = std::log(s[" << (sp-1) << "]);\n";
          break;
        case cLog2:
          //ccfile << "s[" << (sp-1) << "] = std::log2(s[" << (sp-1) << "]);\n"; // C++11
          ccfile << "s[" << (sp-1) << "] = std::log(s[" << (sp-1) << "])/log(2.0);\n"; // C++11
          break;
        case cLog10:
          ccfile << "s[" << (sp-1) << "] = std::log10(s[" << (sp-1) << "]);\n";
          break;

        case cSqr:
          ccfile << "s[" << (sp-1) << "] *= s[" << (sp-1) << "];\n";
          break;
        case cSqrt:
          ccfile << "s[" << (sp-1) << "] = std::pow(s[" << (sp-1) << "], 0.5);\n";
          break;
        case cRSqrt:
          ccfile << "s[" << (sp-1) << "] = std::pow(s[" << (sp-1) << "], -0.5);\n";
          break;
        case cPow:
          ccfile << "s[" << (sp-2) << "] = std::pow(s[" << (sp-2) << "], s[" << (sp-1) << "]);\n";
          break;
        case cExp:
          ccfile << "s[" << (sp-1) << "] = std::exp(s[" << (sp-1) << "]);\n";
          break;
        case cExp2:
          ccfile << "s[" << (sp-1) << "] = std::pow(2.0, s[" << (sp-1) << "]);\n";
          break;

        default:
          std::cerr << "Opcode not supported by JIT.\n";
          ccfile.close();
          return false;
      }
    sp -= op_size > 0 ? op_size-1 : op_size;
  }
  ccfile << "return s[0]; }\n";
  ccfile.close();

  // run compiler
  std::string command = FPARSER_JIT_COMPILER" -O2 -shared -rdynamic -fPIC ";
  command += ccname;
  command += " -o ";
  command += object;
  status = system(command.c_str());
  std::remove(ccname);
  if (status != 0) {
    std::cerr << "JIT compile failed.\n";
    return false;
  }

  // load compiled object
  lib = dlopen(object, RTLD_NOW);
  if (lib == NULL) {
    std::cerr << "JIT object load failed.\n";
    std::remove(object);
    return false;
  }

  // fetch function pointer
  *(void **) (&compiledFunction) = dlsym(lib, fnname.c_str());
  char * error;
  if ((error = dlerror()) != NULL)  {
    std::cerr << "Error binding JIT compiled function\n" << error << '\n';
    compiledFunction = NULL;
    std::remove(object);
    return false;
  }

  // rename successfully compiled obj to cache file
  if (cacheFunction)
  {
    status = std::rename(object, libname.c_str());
    if (status != 0)
      std::remove(object);
  }
  else
    std::remove(object);

  return true;
}
#endif


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
