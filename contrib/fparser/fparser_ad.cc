#include "fparser_ad.hh"
#include "extrasrc/fpaux.hh"
#include "extrasrc/fptypes.hh"
#include "extrasrc/fpaux.hh"
#include <stdlib.h>
#include "Faddeeva.hh"

using namespace FUNCTIONPARSERTYPES;

#include "fpoptimizer/optimize.hh"
#include "fpoptimizer/codetree.hh"
using namespace FPoptimizer_CodeTree;

#include <iostream>
#include <fstream>
#include <cstdio>
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>

#if LIBMESH_HAVE_FPARSER_JIT
#  include <dlfcn.h>
#endif

#include "lib/sha1.h"

/**
 * The internals of the automatic differentiation algorithm are encapsulated in this class
 * and hidden from the public interface, as the installed FParser version does not have access
 * to the CodeTree<> class template. (see pimpl idiom).
 * The ADImplementation has full access to its FParser owner object through the parser pointer.
 */
template<typename Value_t>
class ADImplementation
{
  FunctionParserADBase<Value_t> * parser;
  typename FunctionParserADBase<Value_t>::Data * mData;
public:
  ADImplementation(FunctionParserADBase<Value_t> * _parser) :
      parser(_parser),
      UnsupportedOpcodeException(),
      RefuseToTakeCrazyDerivativeException() {}
  int AutoDiff(unsigned int, typename FunctionParserADBase<Value_t>::Data * mData, bool autoOptimize);

private:
  /**
   * CodeTreeAD is a helper class derived from CodeTree<> for the purpose of adding overloaded
   * operators, to express the elementary derivatives in a compact way.
   */
  class CodeTreeAD : public CodeTree<Value_t> {
  public:
    CodeTreeAD() : CodeTree<Value_t>() {}
    CodeTreeAD(const CodeTree<Value_t> & tree) : CodeTree<Value_t>(tree) {}
    CodeTreeAD(const Value_t & val) : CodeTree<Value_t>(CodeTreeImmed(val)) {}

    CodeTreeAD operator* (const CodeTreeAD & a) { return ADImplementation<Value_t>::MakeTree(cMul, *this, a); }
    CodeTreeAD operator/ (const CodeTreeAD & a) { return ADImplementation<Value_t>::MakeTree(cDiv, *this, a); }
    CodeTreeAD operator+ (const CodeTreeAD & a) { return ADImplementation<Value_t>::MakeTree(cAdd, *this, a); }
    CodeTreeAD operator- (const CodeTreeAD & a) { return ADImplementation<Value_t>::MakeTree(cSub, *this, a); }
    CodeTreeAD operator- () { return ADImplementation<Value_t>::MakeTree(cNeg, *this); }
  };

  // recursive subtree diferentiation - the heart of the AD algorithm
  CodeTreeAD D(const CodeTreeAD & func);

  // "named constructors" to build trees with parameters in a compact way
  static CodeTreeAD MakeTree(OPCODE op, const CodeTreeAD & param1);
  static CodeTreeAD MakeTree(OPCODE op, const CodeTreeAD & param1, const CodeTreeAD & param2);
  static CodeTreeAD MakeTree(OPCODE op, const CodeTreeAD & param1, const CodeTreeAD & param2, const CodeTreeAD & param3);

  // variable index we are differentiating w.r.t.
  unsigned int var;

  // Exceptions
  class UnsupportedOpcode : public std::exception {
    virtual const char* what() const throw() { return "Unsupported opcode"; }
  } UnsupportedOpcodeException;
  class RefuseToTakeCrazyDerivative : public std::exception {
    virtual const char* what() const throw() { return "The derivative of this expression would be undefined at a countable number of points."; }
  } RefuseToTakeCrazyDerivativeException;
};

template<typename Value_t>
FunctionParserADBase<Value_t>::FunctionParserADBase() :
    FunctionParserBase<Value_t>(),
    compiledFunction(NULL),
    mFPlog(this->mData->mFuncPtrs.size()),
    mADFlags(ADJITCache),
    mRegisteredDerivatives(),
    ad(new ADImplementation<Value_t>(this))
{
  this->AddFunction("plog", fp_plog, 2);
  mFErf = this->mData->mFuncPtrs.size();
  this->AddFunction("erf", fp_erf, 1);
}

template<typename Value_t>
FunctionParserADBase<Value_t>::FunctionParserADBase(const FunctionParserADBase& cpy) :
    FunctionParserBase<Value_t>(cpy),
    compiledFunction(cpy.compiledFunction),
    mFPlog(cpy.mFPlog),
    mFErf(cpy.mFErf),
    mADFlags(cpy.mADFlags),
    mRegisteredDerivatives(cpy.mRegisteredDerivatives),
    ad(new ADImplementation<Value_t>(this))
{
  pImmed = this->mData->mImmed.empty() ? NULL : &(this->mData->mImmed[0]);
}

template<typename Value_t>
Value_t FunctionParserADBase<Value_t>::fp_plog(const Value_t * params)
{
  const Value_t x = params[0];
  const Value_t a = params[1];
  return x < a ? fp_log(a) + (x-a)/a - (x-a)*(x-a)/(Value_t(2)*a*a) + (x-a)*(x-a)*(x-a)/(Value_t(3)*a*a*a) : fp_log(x);
}

template <typename Value_t>
Value_t FunctionParserADBase<Value_t>::fp_erf(const Value_t * params)
{
  return Faddeeva::erf(params[0]);
}

template<>
std::complex<float>
FunctionParserADBase<std::complex<float> >::fp_erf(const std::complex<float> * params)
{
  std::complex<double> result =
    Faddeeva::erf(std::complex<double>(params[0].real(),
                                       params[0].imag()));
  return std::complex<float>(result.real(), result.imag());
}

template<>
std::complex<long double>
FunctionParserADBase<std::complex<long double> >::fp_erf(const std::complex<long double> * params)
{
  std::complex<double> result =
    Faddeeva::erf(std::complex<double>(params[0].real(),
                                       params[0].imag()));
  return std::complex<long double>(result.real(), result.imag());
}

template<typename Value_t>
FunctionParserADBase<Value_t>::~FunctionParserADBase()
{
  delete ad;
}

template<typename Value_t>
bool FunctionParserADBase<Value_t>::AddVariable(const std::string & var_name)
{
  this->CopyOnWrite();

  // append new variable to variables string
  const std::string & vars = this->mData->mVariablesString;
  return this->ParseVariables(vars == "" ? var_name : vars + "," + var_name);
}

template<typename Value_t>
bool FunctionParserADBase<Value_t>::isZero()
{
  // determine if the program is a single cImmed 0
  return (this->mData->mByteCode.size() == 1  && this->mData->mImmed.size() == 1 &&
          this->mData->mByteCode[0] == cImmed && this->mData->mImmed[0] == Value_t(0));
}

template<typename Value_t>
void FunctionParserADBase<Value_t>::setZero()
{
  this->CopyOnWrite();

  // set program to a single cImmed 0
  this->mData->mByteCode.resize(1);
  this->mData->mImmed.resize(1);
  this->mData->mByteCode[0] = cImmed;
  this->mData->mImmed[0] = Value_t(0);
}


// this is a namespaced function because we cannot easily export CodeTree in the
// public interface of the FunctionParserADBase class in its installed state in libMesh
// as the codetree.hh header is not installed (part of FPoptimizer)

// quick subtree builder helper functions
template<typename Value_t>
typename ADImplementation<Value_t>::CodeTreeAD ADImplementation<Value_t>::MakeTree(OPCODE op, const CodeTreeAD & param1)
{
  CodeTreeAD tree = CodeTreeAD(CodeTreeOp<Value_t>(op));
  tree.AddParam(param1);
  tree.Rehash();
  return tree;
}
template<typename Value_t>
typename ADImplementation<Value_t>::CodeTreeAD ADImplementation<Value_t>::MakeTree(OPCODE op, const CodeTreeAD & param1, const CodeTreeAD & param2)
{
  CodeTreeAD tree = CodeTreeAD(CodeTreeOp<Value_t>(op));
  tree.AddParam(param1);
  tree.AddParam(param2);
  tree.Rehash();
  return tree;
}
template<typename Value_t>
typename ADImplementation<Value_t>::CodeTreeAD ADImplementation<Value_t>::MakeTree(OPCODE op, const CodeTreeAD & param1, const CodeTreeAD & param2, const CodeTreeAD & param3)
{
  CodeTreeAD tree = CodeTreeAD(CodeTreeOp<Value_t>(op));
  tree.AddParam(param1);
  tree.AddParam(param2);
  tree.AddParam(param3);
  tree.Rehash();
  return tree;
}


// return the derivative of func and put it into diff
template<typename Value_t>
typename ADImplementation<Value_t>::CodeTreeAD ADImplementation<Value_t>::D(const CodeTreeAD & func)
{
  // derivative of a constant number is 0
  if (func.IsImmed())
    return CodeTreeAD(0);

  // derivative of a variable is 1 for the variable we are diffing w.r.t. and 0 otherwise
  if (func.IsVar())
  {
    if (func.GetVar() == var) {
      // dx/dx = 1
      return CodeTreeAD(1);
    } else {
      // check if this derivative is registered
      typename std::vector<typename FunctionParserADBase<Value_t>::VariableDerivative>::const_iterator it;
      for (it = this->parser->mRegisteredDerivatives.begin();
           it != this->parser->mRegisteredDerivatives.end(); ++it)
        if (it->var == func.GetVar() && it->dependence == var)
          return CodeTreeAD(CodeTreeVar<Value_t>(it->derivative));

      // otherwise return zero
      return CodeTreeAD(0);
    }
  }

  // derivative being built for regular opcodes
  CodeTreeAD diff;

  // get opcode and parameter list
  OPCODE op = func.GetOpcode();
  const std::vector<CodeTree<Value_t> > & param_plain = func.GetParams();

  // augment the params to use the CodeTreeAD class
  std::vector<CodeTreeAD> param;
  unsigned int i, j, nparam = param_plain.size();
  for (i = 0; i < nparam; ++i)
    param.push_back(CodeTreeAD(param_plain[i]));

  // short hand for the first three parameters
  CodeTreeAD a = nparam > 0 ? param[0] : CodeTreeAD(),
             b = nparam > 1 ? param[1] : CodeTreeAD(),
             c = nparam > 2 ? param[2] : CodeTreeAD();

  switch (op)
  {
    //
    // these opcodes can take an arbitrary number of parameters
    //

    case cAdd:
    case cSub:
      diff.SetOpcode(op);
      for (i = 0; i < nparam; ++i)
        diff.AddParam(D(param[i]));
      diff.Rehash();
      break;

    case cMul:
      diff.SetOpcode(cAdd);
      for (i = 0; i < nparam; ++i)
      {
        CodeTree<Value_t> sub;
        sub.SetOpcode(cMul);
        for (j = 0; j < nparam; ++j)
          sub.AddParam(i==j ? D(param[j]) : param[j]);
        sub.Rehash();
        diff.AddParam(sub);
      }
      diff.Rehash();
      break;

    //
    // these opcodes can take a fixed number of parameters
    //

    case cDiv:
      // da/b - db/b^2
      return D(a)/b - D(b) / MakeTree(cPow, b, CodeTreeAD(2));

    case cNeg:
      return -D(a);
    case cInv:
      return -D(a) / MakeTree(cPow, a, CodeTreeAD(2));

    case cSin:
      return D(a) * MakeTree(cCos, a);
    case cCos:
      return -D(a) * MakeTree(cSin, a);

    case cSinh:
      return D(a) * MakeTree(cCosh, a);
    case cCosh:
      return D(a) * MakeTree(cSinh, a); // no -

    case cAsin:
      // da/sqrt(1-a^2)
      return D(a) / MakeTree(cSqrt, CodeTreeAD(1) - MakeTree(cSqr, a));
    case cAcos:
      // -da/sqrt(1-a^2)
      return -D(a) / MakeTree(cSqrt, CodeTreeAD(1) - MakeTree(cSqr, a));
    case cAtan :
      // da/(a^2+1)
      return D(a) / (MakeTree(cSqr, a) + CodeTreeAD(1));
    case cAtan2 :
      // (b*da-a*db)/(a^2+b^2)
      return (b * D(a) - a * D(b)) / (MakeTree(cSqr, a) + MakeTree(cSqr, b));

    case cPow:
      if (b.IsImmed())
      {
        Value_t exponent = b.GetImmed();
        if (exponent == Value_t(1))
          return D(a);
        if (exponent == Value_t(0))
          return CodeTreeAD(0);
        return MakeTree(cPow, a, CodeTreeAD(exponent - Value_t(1))) * b * D(a);
      }
      return MakeTree(cPow, a, b) * (D(b) * MakeTree(cLog, a) + b * D(a)/a);
    case cLog:
      return D(a)/a;
    case cExp:
      return D(a) * MakeTree(cExp, a);

    //
    // the derivatives of these could be undefined at a finite number of points (one)
    //

    case cIf:
      // we diff the two branches, but not the condition
      diff.SetOpcode(cIf);
      diff.AddParam(a);
      diff.AddParam(D(b));
      diff.AddParam(D(c));
      diff.Rehash();
      break;

    case cAbs:
      // da*a/|a|
      return D(a)*a/MakeTree(cAbs, a);

    case cMax:
      diff.SetOpcode(cIf);
      diff.AddParam(MakeTree(cLess, a, b));
      diff.AddParam(D(b));
      diff.AddParam(D(a));
      diff.Rehash();
      break;
    case cMin:
      diff.SetOpcode(cIf);
      diff.AddParam(MakeTree(cLess, a, b));
      diff.AddParam(D(a));
      diff.AddParam(D(b));
      diff.Rehash();
      break;

    //
    // the derivatives of these are undefined in a countable number of points, we refuse to take them
    // TODO: we could wrap those in an if statement, that returns 1/0 at the undefined points!
    //

    case cFloor:
    case cCeil:
    case cTrunc:
      throw RefuseToTakeCrazyDerivativeException;

    // we gloss over the single undifferentiable point here

    case cEqual:
    case cNEqual:
    case cLess:
    case cLessOrEq:
    case cGreater:
    case cGreaterOrEq:
      return CodeTreeAD(0);

    // these opcodes will never appear in the tree when building with keep_powi == false
    // they will be replaced by cPow, cMul, cDiv etc.:
    // cLog10, cLog2, cHypot, cExp2, cRSqrt, cCbrt, cInv, cCot, cCsc, cSec
    // cSqr, cSqrt, cInt, cLog2by, cSinCos, cSinhCosh, cRSub, cTan, cTanh

    case cFCall:
      if (func.GetFuncNo() == this->parser->mFPlog)
      {
        // a<b ? D(a) * [(1/b - 1/(b*b)*(a-b) + 1/b^3 * (a-b)^2) : 1/a]
        return D(a) * MakeTree(cIf,
          MakeTree(cLess, a, b),
          MakeTree(cInv, b) - MakeTree(cInv, b*b) * (a-b)
            + MakeTree(cPow, b, CodeTreeAD(-3)) * MakeTree(cSqr, a-b),
            MakeTree(cInv, a)
        );
      }
      else if (func.GetFuncNo() == this->parser->mFErf)
      {
        return D(a) * CodeTreeAD(Value_t(2) / fp_sqrt<Value_t>(M_PI)) * MakeTree(cExp, -(a*a));
      }
      // fall through to undefined

    case cPCall:
    default:
      throw UnsupportedOpcodeException;
  }

  return diff;
}

template<typename Value_t>
int FunctionParserADBase<Value_t>::AutoDiff(const std::string& var_name)
{
  this->CopyOnWrite();
  const bool cached = mADFlags & ADCacheDerivatives;

  try
  {
    unsigned int var_number = LookUpVarOpcode(var_name);

    // should and can we load a cached derivative?
    std::string cachename;
    const std::string jitdir = ".jitcache";
    if (cached)
    {
      // generate a sha1 hash of the Value type size, byte code, immediate list, and registered derivatives
      SHA1 *sha1 = new SHA1();
      char result[41]; // 40 sha1 chars plus null
      size_t value_t_size = sizeof(Value_t);
      sha1->addBytes(reinterpret_cast<const char *>(&value_t_size), sizeof(value_t_size));
      sha1->addBytes(reinterpret_cast<const char *>(&this->mData->mByteCode[0]), this->mData->mByteCode.size() * sizeof(unsigned));
      if (!this->mData->mImmed.empty())
        sha1->addBytes(reinterpret_cast<const char *>(&this->mData->mImmed[0]), this->mData->mImmed.size() * sizeof(Value_t));
      for (const auto & reg : this->mRegisteredDerivatives)
        sha1->addBytes(reinterpret_cast<const char *>(&reg), sizeof(VariableDerivative));

      unsigned char* digest = sha1->getDigest();
      for (unsigned int i = 0; i<20; ++i)
        sprintf(&(result[i*2]), "%02x", digest[i]);
      free(digest);
      delete sha1;

      // variable number
      char var_number_string[sizeof(unsigned int) * 2 + 1]; // two hex chars per byte plus null
      sprintf(var_number_string, "%x", var_number - VarBegin);

      // function name
      cachename = jitdir + "/d_";
      cachename += result;
      cachename += "_";
      cachename += var_number_string;

      // try to open cache file
      std::ifstream istr;
      istr.open(cachename.c_str(), std::ios::in | std::ios::binary);
      if (istr)
      {
        Unserialize(istr);
        // only claim success if the stream is not in a bad state
        if (istr.good()) return -1;
      }
    }

    // immediately optimize the derivative tree representation?
    const bool autoOptimize = mADFlags & ADAutoOptimize;

    // build derivative
    int result = ad->AutoDiff(var_number, this->mData, autoOptimize);

    // save the derivative if cacheing is enabled and derivative was successfully taken
    if (cached && result == -1)
    {
      // create cache directory
      if (mkdir(jitdir.c_str(), 0700) == 0 || errno == EEXIST)
      {
        // save to a temporary name and rename only when the file is fully written
        char cachetmpname[] = "./tmp_adc_XXXXXX";
        int cachetmpfile = mkstemp(cachetmpname);
        if (cachetmpfile == -1)
          std::cerr << "Error creating AD cache tmp file " << cachetmpname << ".\n";
        else
        {
          close(cachetmpfile);
          std::ofstream ostr;
          ostr.open(cachetmpname, std::ios::out | std::ios::binary);
          if (ostr)
          {
            Serialize(ostr);
            ostr.close();

            /**
             * MPI runs on asynchronous networked filesystems require a two-step
             * renaming. The link call will not do anything if the file cachename
             * has already been created by a different rank.
             */
            link(cachetmpname, cachename.c_str());
            std::remove(cachetmpname);
          }
        }
      }
    }

    return result;
  }
  catch(std::exception &e)
  {
    static bool printed_error = false;
    const bool silence_errors = mADFlags & ADSilenceErrors;
    if (!printed_error && !silence_errors)
    {
      std::cerr << "AutoDiff exception: " << e.what() << " (this message will only be shown once per process)"<< std::endl;
      printed_error = true;
    }
    setZero();
    return 0;
  }
}

template<typename Value_t>
unsigned int
FunctionParserADBase<Value_t>::LookUpVarOpcode(const std::string & var_name)
{
  // get c string and length of var argument
  const unsigned len = unsigned(var_name.size());
  const char* name = var_name.c_str();

  typename FUNCTIONPARSERTYPES::NamePtrsMap<Value_t> & NamePtrs = this->mData->mNamePtrs;
  typename FUNCTIONPARSERTYPES::NamePtrsMap<Value_t>::iterator vi;
  for (vi = NamePtrs.begin(); vi != NamePtrs.end(); ++vi)
    if (len == vi->first.nameLength && std::memcmp(name, vi->first.name, len) == 0)
      return vi->second.index;

  throw UnknownVariableException;
}

template<typename Value_t>
void
FunctionParserADBase<Value_t>::RegisterDerivative(const std::string & a, const std::string & b, const std::string & c)
{
  VariableDerivative rule;
  rule.var        = LookUpVarOpcode(a);
  rule.dependence = LookUpVarOpcode(b);
  rule.derivative = LookUpVarOpcode(c);
  mRegisteredDerivatives.push_back(rule);
}

template<typename Value_t>
int ADImplementation<Value_t>::AutoDiff(unsigned int _var, typename FunctionParserADBase<Value_t>::Data * mData, bool autoOptimize)
{
  CodeTreeAD orig;
  orig.GenerateFrom(*mData, false);

  // store variable we are diffing for as a member
  var = _var;

  // start recursing the code tree
  CodeTree<Value_t> diff = D(orig);
#ifndef FP_DUMMY_OPTIMIZER
  if (autoOptimize)
    FPoptimizer_Optimize::ApplyGrammars(diff);
#endif

  std::vector<unsigned> byteCode;
  std::vector<Value_t> immed;
  size_t stacktop_max = 0;

  diff.SynthesizeByteCode(byteCode, immed, stacktop_max);

  if(mData->mStackSize != stacktop_max)
  {
    mData->mStackSize = unsigned(stacktop_max); // Note: Ignoring GCC warning here.
  #if !defined(FP_USE_THREAD_SAFE_EVAL) && \
  !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
    mData->mStack.resize(stacktop_max);
  #endif
  }

  mData->mByteCode.swap(byteCode);
  mData->mImmed.swap(immed);

  return -1;
}

template<typename Value_t>
void FunctionParserADBase<Value_t>::Optimize()
{
  FunctionParserBase<Value_t>::Optimize();
  if (compiledFunction)
    JITCompile();
}

template<typename Value_t>
bool FunctionParserADBase<Value_t>::JITCompile()
{
  // JIT compile attempted for an unsupported value type
  return false;
}

#if LIBMESH_HAVE_FPARSER_JIT
// JIT compile for supported types
template<>
bool FunctionParserADBase<double>::JITCompile() { return JITCompileHelper("double"); }
template<>
bool FunctionParserADBase<float>::JITCompile() { return JITCompileHelper("float"); }
template<>
bool FunctionParserADBase<long double>::JITCompile() { return JITCompileHelper("long double"); }

template<typename Value_t>
Value_t FunctionParserADBase<Value_t>::Eval(const Value_t* Vars)
{
  if (compiledFunction == NULL)
    return FunctionParserBase<Value_t>::Eval(Vars);
  else
    return (*compiledFunction)(Vars, pImmed, Epsilon<Value_t>::value);
}

template<typename Value_t>
bool FunctionParserADBase<Value_t>::JITCompileHelper(const std::string & Value_t_name)
{
  // use the file cache for compiled functions?
  const bool cacheFunction = mADFlags & ADJITCache;

  // set compiled function pointer to zero to avoid stale values if JIT compilation fails
  compiledFunction = NULL;

  // get a pointer to the mImmed values
  pImmed = this->mData->mImmed.empty() ? NULL : &(this->mData->mImmed[0]);

  // get a reference to the stored bytecode
  const std::vector<unsigned>& ByteCode = this->mData->mByteCode;

  // drop out if the ByteCode is empty
  if (ByteCode.empty())
    return false;

  // generate a sha1 hash of the current program and the Value type name
  SHA1 *sha1 = new SHA1();
  char result[41];
  sha1->addBytes(reinterpret_cast<const char *>(&ByteCode[0]), ByteCode.size() * sizeof(unsigned));
  sha1->addBytes(Value_t_name.c_str(), Value_t_name.size());
  unsigned char* digest = sha1->getDigest();
  for (unsigned int i = 0; i<20; ++i)
    sprintf(&(result[i*2]), "%02x", digest[i]);
  free(digest);
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
  // opening the cached file did not work. (re)build it.

  // tmp filenames
  char ccname[] = "./tmp_jit_XXXXXX";
  int ccfile = mkstemp(ccname);
  if (ccfile == -1)
  {
    std::cerr << "Error creating JIT tmp file " << ccname << ".\n";
    return false;
  }
  close(ccfile); // close file. We reopen this using an ofstream below.

  char objectname[] = "./tmp_jit_XXXXXX";
  int objectfile = mkstemp(objectname);
  if (objectfile == -1)
  {
    std::cerr << "Error creating JIT tmp file " << objectname << ".\n";
    std::remove(ccname);
    return false;
  }
  close(objectfile); // close file. Will be reopened by the compiler.

  int status;
  std::vector<bool> jumpTarget(ByteCode.size());
  for (unsigned int i = 0; i < ByteCode.size(); ++i) jumpTarget[i] = false;

  // determine all jump targets in the current program
  unsigned long ip;
  for (unsigned int i = 0; i < ByteCode.size(); ++i)
    switch (ByteCode[i])
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
  std::ofstream ccout;
  ccout.open(ccname);
  ccout << "#define _USE_MATH_DEFINES\n";
  ccout << "#include <cmath>\n";
  ccout << "extern \"C\" " << Value_t_name << ' '
         << fnname << "(const " << Value_t_name << " *params, const " << Value_t_name << " *immed, const " << Value_t_name << " eps) {\n";
  ccout << Value_t_name << " r, s[" << this->mData->mStackSize << "];\n";
  int nImmed = 0, sp = -1, op;
  for (unsigned int i = 0; i < ByteCode.size(); ++i)
  {
    // place a label?
    if (jumpTarget[i])
    {
      ccout << "l" << i << ":\n";
      if (stackAtTarget[i] > -2)
        sp = stackAtTarget[i];
    }

    // translate opcode into C code
    switch (op = ByteCode[i])
    {
      case cImmed:
        ++sp; ccout << "s[" << sp << "] = immed[" << nImmed++ << "];\n"; break;
      case cAdd:
        --sp; ccout << "s[" << sp << "] += s[" << (sp+1) << "];\n"; break;
      case cSub:
        --sp; ccout << "s[" << sp << "] -= s[" << (sp+1) << "];\n"; break;
      case cRSub:
        --sp; ccout << "s[" << sp << "] = s[" << (sp+1) << "] - s[" << sp << "];\n"; break;
      case cMul:
        --sp; ccout << "s[" << sp << "] *= s[" << (sp+1) << "];\n"; break;
      case cDiv:
        --sp; ccout << "s[" << sp << "] /= s[" << (sp+1) << "];\n"; break;
      case cMod:
        --sp; ccout << "s[" << sp << "] = std::fmod(s[" << sp << "], s[" << (sp+1) << "]);\n"; break;
      case cRDiv:
        --sp; ccout << "s[" << sp << "] = s[" << (sp+1) << "] / s[" << sp << "];\n"; break;

      case cSin:
        ccout << "s[" << sp << "] = std::sin(s[" << sp << "]);\n"; break;
      case cCos:
        ccout << "s[" << sp << "] = std::cos(s[" << sp << "]);\n"; break;
      case cTan:
        ccout << "s[" << sp << "] = std::tan(s[" << sp << "]);\n"; break;
      case cSinh:
        ccout << "s[" << sp << "] = std::sinh(s[" << sp << "]);\n"; break;
      case cCosh:
        ccout << "s[" << sp << "] = std::cosh(s[" << sp << "]);\n"; break;
      case cTanh:
        ccout << "s[" << sp << "] = std::tanh(s[" << sp << "]);\n"; break;
      // TODO: div by zero -> this->mData->mEvalErrorType=1; return Value_t(0);
      case cCsc:
        ccout << "s[" << sp << "] = 1.0/std::sin(s[" << sp << "]);\n"; break;
      case cSec:
        ccout << "s[" << sp << "] = 1.0/std::cos(s[" << sp << "]);\n"; break;
      case cCot:
        ccout << "s[" << sp << "] = 1.0/std::tan(s[" << sp << "]);\n"; break;
      case cSinCos:
        ccout << "s[" << (sp+1) << "] = std::cos(s[" << sp << "]);\n";
        ccout << "s[" << sp << "] = std::sin(s[" << sp << "]);\n";
        ++sp;
        break;
      case cSinhCosh:
        ccout << "s[" << (sp+1) << "] = std::cosh(s[" << sp << "]);\n";
        ccout << "s[" << sp << "] = std::sinh(s[" << sp << "]);\n";
        ++sp;
        break;
      case cAsin:
        ccout << "s[" << sp << "] = std::asin(s[" << sp << "]);\n"; break;
      case cAcos:
        ccout << "s[" << sp << "] = std::acos(s[" << sp << "]);\n"; break;
      case cAsinh:
        ccout << "s[" << sp << "] = std::asinh(s[" << sp << "]);\n"; break;
      case cAcosh:
        ccout << "s[" << sp << "] = std::acosh(s[" << sp << "]);\n"; break;
      case cAtan:
        ccout << "s[" << sp << "] = std::atan(s[" << sp << "]);\n"; break;
      case cAtanh:
        ccout << "s[" << sp << "] = std::atanh(s[" << sp << "]);\n"; break;
      case cAtan2:
        --sp; ccout << "s[" << sp << "] = std::atan2(s[" << sp << "], s[" << (sp+1) << "]);\n"; break;
      case cHypot:
        --sp; ccout << "s[" << sp << "] = std::sqrt(s[" << sp << "]*s[" << sp << "] + s[" << (sp+1) << "]*s[" << (sp+1) << "]);\n"; break;

      case cAbs:
        ccout << "s[" << sp << "] = std::abs(s[" << sp << "]);\n"; break;
      case cMax:
        --sp; ccout << "s[" << sp << "] = s[" << sp << "] > s[" << (sp+1) << "] ? s[" << sp << "] : s[" << (sp+1) << "];\n"; break;
      case cMin:
        --sp; ccout << "s[" << sp << "] = s[" << sp << "] < s[" << (sp+1) << "] ? s[" << sp << "] : s[" << (sp+1) << "];\n"; break;
      case cTrunc:
        ccout << "s[" << sp << "] = s[" << sp << "] < 0 ? std::ceil(s[" << sp << "]) : std::floor(s[" << sp << "]);\n"; break;
      case cCeil:
        ccout << "s[" << sp << "] = std::ceil(s[" << sp << "]);\n"; break;
      case cFloor:
        ccout << "s[" << sp << "] = std::floor(s[" << sp << "]);\n"; break;
      case cInt:
        ccout << "s[" << sp << "] = s[" << sp << "] < 0 ? std::ceil(s[" << sp << "] - 0.5) : std::floor(s[" << sp << "] + 0.5);\n"; break;

      case cEqual:
        //--sp; ccout << "s[" << sp << "] = s[" << sp << "] == s[" << (sp+1) << "];\n"; break;
        --sp; ccout << "s[" << sp << "] = std::abs(s[" << sp << "] - s[" << (sp+1) << "]) <= eps;\n"; break;
      case cNEqual:
        //--sp; ccout << "s[" << sp << "] = s[" << sp << "] != s[" << (sp+1) << "];\n"; break;
        --sp; ccout << "s[" << sp << "] = std::abs(s[" << sp << "] - s[" << (sp+1) << "]) > eps;\n"; break;
      case cLess:
        --sp; ccout << "s[" << sp << "] = s[" << sp << "] < (s[" << (sp+1) << "] - eps);\n"; break;
      case cLessOrEq:
        --sp; ccout << "s[" << sp << "] = s[" << sp << "] <= (s[" << (sp+1) << "] + eps);\n"; break;
      case cGreater:
        --sp; ccout << "s[" << sp << "] = (s[" << sp << "] - eps) > s[" << (sp+1) << "];\n"; break;
      case cGreaterOrEq:
        --sp; ccout << "s[" << sp << "] = (s[" << sp << "] + eps) >= s[" << (sp+1) << "];\n"; break;
      case cNot:
        ccout << "s[" << sp << "] = std::abs(s[" << sp << "]) < 0.5;\n"; break;
      case cNotNot:
        ccout << "s[" << sp << "] = std::abs(s[" << sp << "]) >= 0.5;\n"; break;
      case cAbsNot:
        ccout << "s[" << sp << "] = s[" << sp << "] < 0.5;\n"; break;
      case cAbsNotNot:
        ccout << "s[" << sp << "] = s[" << sp << "] >= 0.5;\n"; break;
      case cOr:
        --sp; ccout << "s[" << sp << "] = (std::abs(s[" << sp << "]) >= 0.5) || (std::abs(s[" << (sp+1) << "]) >= 0.5);\n"; break;
      case cAbsOr:
        --sp; ccout << "s[" << sp << "] = (s[" << sp << "] >= 0.5) || (s[" << (sp+1) << "] >= 0.5);\n"; break;
      case cAnd:
        --sp; ccout << "s[" << sp << "] = (std::abs(s[" << sp << "]) >= 0.5) && (std::abs(s[" << (sp+1) << "]) >= 0.5);\n"; break;
      case cAbsAnd:
        --sp; ccout << "s[" << sp << "] = (s[" << sp << "] >= 0.5) && (s[" << (sp+1) << "] >= 0.5);\n"; break;

      case cLog:
        ccout << "s[" << sp << "] = std::log(s[" << sp << "]);\n"; break;
      case cLog2:
#ifdef FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS
        ccout << "s[" << (sp-1) << "] = std::log2(s[" << (sp-1) << "]);\n";
#else
        ccout << "s[" << sp << "] = std::log(s[" << sp << "])/log(2.0);\n";
#endif
        break;
      case cLog10:
        ccout << "s[" << sp << "] = std::log10(s[" << sp << "]);\n"; break;

      case cNeg:
        ccout << "s[" << sp << "] = -s[" << sp << "];\n"; break;
      case cInv:
        ccout << "s[" << sp << "] = 1.0/s[" << sp << "];\n"; break;
      case cDeg:
        ccout << "s[" << sp << "] *= 180.0/M_PI;\n"; break;
      case cRad:
        ccout << "s[" << sp << "] /= 180.0/M_PI;\n"; break;

      case cFetch:
        ++sp; ccout << "s[" << sp << "] = s[" << ByteCode[++i] << "];\n"; break;
      case cDup:
        ++sp; ccout << "s[" << sp << "] = s[" << (sp-1) << "];\n"; break;

      case cFCall:
      {
        unsigned function = ByteCode[++i];
        if (function == mFPlog)
        {
          // --sp; ccout << "s[" << sp << "] = s[" << sp << "] < s[" << (sp+1) << "] ? std::log(s[" << (sp+1) << "]) + (s[" << sp << "] - s[" << (sp+1) << "]) / s[" << (sp+1) << "] : std::log(s[" << sp << "]);\n";
          // --sp; ccout << "s[" << sp << "] = s[" << sp << "] < s[" << (sp+1) << "] ? std::log(s[" << (sp+1) << "]) - 1.5 + 2.0/s[" << (sp+1) << "] * s[" << sp << "] - 0.5/(s[" << (sp+1) << "]*s[" << (sp+1) << "]) * s[" << sp << "]*s[" << sp << "] : std::log(s[" << sp << "]);\n";
          --sp; ccout << "s[" << sp << "] = s[" << sp << "] < s[" << (sp+1) << "] ? std::log(s[" << (sp+1) << "])  +  (s[" << sp << "]-s[" << (sp+1) << "])/s[" << (sp+1) << "]  -  std::pow((s[" << sp << "]-s[" << (sp+1) << "])/s[" << (sp+1) << "],2.0)/2.0  +  std::pow((s[" << sp << "]-s[" << (sp+1) << "])/s[" << (sp+1) << "],3.0)/3.0 : std::log(s[" << sp << "]);\n";
        }
        else if (function == mFErf)
        {
#if LIBMESH_HAVE_CXX11_ERF
          ccout << "s[" << sp << "] = std::erf(s[" << sp << "]);\n";
#else
          std::cerr << "Libmesh is not compiled with c++11 so std::erf is not supported by JIT.\n";
          ccout.close();
          std::remove(ccname);
          return false;
#endif
        }
        else
        {
          std::cerr << "Function call not supported by JIT.\n";
          ccout.close();
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
        ccout << "s[" << dst << "] = s[" << src << "];\n";
        sp = dst;
        break;
      }
      case cLog2by:
        --sp; ccout << "s[" << sp << "] = std::log(s[" << sp << "])/log(2.0) * s[" << (sp+1) << "];\n"; break;
      case cNop:
        break;
#endif

      case cSqr:
        ccout << "s[" << sp << "] *= s[" << sp << "];\n"; break;
      case cSqrt:
        ccout << "s[" << sp << "] = std::sqrt(s[" << sp << "]);\n"; break;
      case cRSqrt:
        ccout << "s[" << sp << "] = std::pow(s[" << sp << "], " << Value_t_name << "(-0.5));\n"; break;
      case cPow:
        --sp; ccout << "s[" << sp << "] = std::pow(s[" << sp << "], s[" << (sp+1) << "]);\n"; break;
      case cExp:
        ccout << "s[" << sp << "] = std::exp(s[" << sp << "]);\n"; break;
      case cExp2:
        ccout << "s[" << sp << "] = std::pow(" << Value_t_name << "(2.0), s[" << sp << "]);\n"; break;
      case cCbrt:
#ifdef FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS
        ccout << "s[" << sp << "] = std::cbrt(s[" << sp << "]);\n"; break;
#else
        ccout << "s[" << sp << "] = s[" << sp << "] == 0 ? 0 : (s[" << sp << "] > 0 ? std::exp(std::log(s[" << sp << "])/3.0) : -std::exp(std::log(-s[" << sp << "])/3.0));\n"; break;
#endif

      case cJump:
      case cIf:
      case cAbsIf:
      {
        unsigned long ip = ByteCode[++i] + 1;

        if (op == cIf)
          ccout << "if (std::abs(s[" << sp-- << "]) < 0.5) ";
        if (op == cAbsIf)
          ccout << "if (s[" << sp-- << "] < 0.5) ";

        if (ip >= ByteCode.size())
          ccout << "return s[" << sp << "];\n";
        else
        {
          ccout << "goto l" << ip << ";\n";
          stackAtTarget[ip] = sp;
        }

        ++i;
        break;
      }

      default:
        if ( op>= VarBegin)
        {
          // load variable
          ++sp; ccout << "s[" << sp << "] = params[" << (op-VarBegin) << "];\n";
        }
        else
        {
          std::cerr << "Opcode not supported by JIT.\n";
          ccout.close();
          std::remove(ccname);
          return false;
        }
    }
  }
  ccout << "return s[" << sp << "]; }\n";
  ccout.close();

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
#if defined(__GNUC__) && defined(__APPLE__) && !defined(__INTEL_COMPILER)
  // gcc on OSX does neither need nor accept the  -rdynamic switch
  std::string command = FPARSER_JIT_COMPILER" -O2 -shared -fPIC ";
#else
  std::string command = FPARSER_JIT_COMPILER" -O2 -shared -rdynamic -fPIC ";
#endif
  command += ccname_cc + " -o " + objectname;
  status = system(command.c_str());
  std::remove(ccname_cc.c_str());
  if (status != 0) {
    std::cerr << "JIT compile failed.\n";
    return false;
  }

  // add a .so extension to the object (needed by dlopen on mac)
  std::string object_so = objectname;
  object_so += ".so";
  status = std::rename(objectname, object_so.c_str());
  if (status != 0)
  {
    std::cerr << "Unable to rename JIT compiled function object\n";
    std::remove(objectname);
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
  const char * error = dlerror();
  if (error != NULL)  {
    std::cerr << "Error binding JIT compiled function\n" << error << '\n';
    compiledFunction = NULL;
    std::remove(object_so.c_str());
    return false;
  }

  // clear evalerror (this will not get set again by the JIT code)
  this->mData->mEvalErrorType = 0;

  // hard link successfully compiled obj to final cache file
  if (cacheFunction && (mkdir(jitdir.c_str(), 0700) == 0 || errno == EEXIST)) {
    // the directory was either successfully created, or it already exists
    link(object_so.c_str(), libname.c_str());
  }

  // remove temporary object
  std::remove(object_so.c_str());
  return true;
}
#endif

template<typename Value_t>
void FunctionParserADBase<Value_t>::Serialize(std::ostream & ostr)
{
  // write version header
  const int version = 1;
  ostr.write(reinterpret_cast<const char *>(&version), sizeof(version));

  // write bytecode buffer
  const size_t byte_code_size = this->mData->mByteCode.size();
  ostr.write(reinterpret_cast<const char *>(&byte_code_size), sizeof(byte_code_size));
  if (byte_code_size > 0  )
    ostr.write(reinterpret_cast<const char *>(&this->mData->mByteCode[0]), byte_code_size * sizeof(unsigned));

  // write immediates
  const size_t immed_size = this->mData->mImmed.size();
  ostr.write(reinterpret_cast<const char *>(&immed_size), sizeof(immed_size));
  if (immed_size > 0)
    ostr.write(reinterpret_cast<const char *>(&this->mData->mImmed[0]), immed_size * sizeof(Value_t));

  // write stacktop max
  ostr.write(reinterpret_cast<const char *>(&this->mData->mStackSize), sizeof(this->mData->mStackSize));
}

template<typename Value_t>
void FunctionParserADBase<Value_t>::Unserialize(std::istream & istr)
{
  // read version header
  int version;
  istr.read(reinterpret_cast<char *>(&version), sizeof(version));
  if (version != 1) throw UnknownSerializationVersionException;

  // read bytecode buffer
  std::vector<unsigned> byteCode;
  size_t byte_code_size;
  istr.read(reinterpret_cast<char *>(&byte_code_size), sizeof(byte_code_size));
  if (byte_code_size > 0)
  {
    byteCode.resize(byte_code_size);
    istr.read(reinterpret_cast<char *>(&byteCode[0]), byte_code_size * sizeof(unsigned));
  }

  // read immediates
  std::vector<Value_t> immed;
  size_t immed_size;
  istr.read(reinterpret_cast<char *>(&immed_size), sizeof(immed_size));
  if (immed_size > 0)
  {
    immed.resize(immed_size);
    istr.read(reinterpret_cast<char *>(&immed[0]), immed_size * sizeof(Value_t));
  }

  // read stacktop
  unsigned stacktop_max;
  istr.read(reinterpret_cast<char *>(&stacktop_max), sizeof(unsigned));

  if(this->mData->mStackSize != stacktop_max)
  {
    this->mData->mStackSize = stacktop_max;
  #if !defined(FP_USE_THREAD_SAFE_EVAL) && \
  !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
    this->mData->mStack.resize(stacktop_max);
  #endif
  }

  this->mData->mByteCode.swap(byteCode);
  this->mData->mImmed.swap(immed);
}

#define FUNCTIONPARSERAD_INSTANTIATE_CLASS(type) \
    template class FunctionParserADBase< type >; \
    template class ADImplementation< type >;

#ifndef FP_DISABLE_DOUBLE_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(double)
#endif

#ifdef FP_SUPPORT_FLOAT_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(float)
#endif

#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(long double)
#endif

#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(std::complex<double>)
#endif

#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(std::complex<float>)
#endif

#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
FUNCTIONPARSERAD_INSTANTIATE_CLASS(std::complex<long double>)
#endif
