// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_PARSED_FUNCTION_PROGRAM_H
#define LIBMESH_PARSED_FUNCTION_PROGRAM_H

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_device.h"

#include <vector>

template <typename Value_t>
class FunctionParserBase;
template <typename Value_t>
class FunctionParserADBase;

namespace libMesh
{

enum class ParsedFunctionOpcode : unsigned int
{
  cAbs,
  cAcos,
  cAcosh,
  cArg,
  cAsin,
  cAsinh,
  cAtan,
  cAtan2,
  cAtanh,
  cCbrt,
  cCeil,
  cConj,
  cCos,
  cCosh,
  cCot,
  cCsc,
  cExp,
  cExp2,
  cFloor,
  cHypot,
  cIf,
  cImag,
  cInt,
  cLog,
  cLog10,
  cLog2,
  cMax,
  cMin,
  cPolar,
  cPow,
  cReal,
  cSec,
  cSin,
  cSinh,
  cSqrt,
  cTan,
  cTanh,
  cTrunc,
  cImmed,
  cJump,
  cNeg,
  cAdd,
  cSub,
  cMul,
  cDiv,
  cMod,
  cEqual,
  cNEqual,
  cLess,
  cLessOrEq,
  cGreater,
  cGreaterOrEq,
  cNot,
  cAnd,
  cOr,
  cNotNot,
  cDeg,
  cRad,
  cFCall,
  cPCall,
  cPopNMov,
  cLog2by,
  cNop,
  cSinCos,
  cSinhCosh,
  cAbsAnd,
  cAbsOr,
  cAbsNot,
  cAbsNotNot,
  cAbsIf,
  cDup,
  cFetch,
  cInv,
  cSqr,
  cRDiv,
  cRSub,
  cRSqrt,
  VarBegin
};

LIBMESH_DEVICE_INLINE constexpr unsigned int
parsed_function_var_begin()
{
  return static_cast<unsigned int>(ParsedFunctionOpcode::VarBegin);
}

LIBMESH_DEVICE_INLINE constexpr bool
parsed_function_is_var_opcode(const unsigned int opcode)
{
  return opcode >= parsed_function_var_begin();
}

template <typename Scalar>
struct ParsedFunctionProgram
{
  std::vector<unsigned int> bytecode;
  std::vector<Scalar> immediates;
  unsigned int stack_size = 0;
  unsigned int n_variables = 0;
  Scalar epsilon = 0;

  bool empty() const { return bytecode.empty(); }
};

template <typename Scalar>
struct ParsedFunctionProgramBundle
{
  ParsedFunctionProgram<Scalar> value;
  ParsedFunctionProgram<Scalar> dx;
#if LIBMESH_DIM > 1
  ParsedFunctionProgram<Scalar> dy;
#endif
#if LIBMESH_DIM > 2
  ParsedFunctionProgram<Scalar> dz;
#endif
  ParsedFunctionProgram<Scalar> dt;
};

template <typename Scalar>
struct ParsedFEMFunctionProgramBundle
{
  ParsedFunctionProgram<Scalar> value;
  ParsedFunctionProgram<Scalar> dx;
#if LIBMESH_DIM > 1
  ParsedFunctionProgram<Scalar> dy;
#endif
#if LIBMESH_DIM > 2
  ParsedFunctionProgram<Scalar> dz;
#endif
  ParsedFunctionProgram<Scalar> dt;
  std::vector<unsigned int> value_variable_numbers;
  std::vector<ParsedFunctionProgram<Scalar>> value_variable_derivatives;
  bool uses_field_gradients = false;
  bool uses_field_hessians = false;
  bool uses_normals = false;
  bool uses_additional_variables = false;

  bool supports_kokkos_value_goal() const
  {
    return !uses_field_gradients &&
           !uses_field_hessians &&
           !uses_normals &&
           !uses_additional_variables &&
           value_variable_numbers.size() == value_variable_derivatives.size();
  }
};

template <typename Scalar>
ParsedFunctionProgram<Scalar>
build_parsed_function_program(const FunctionParserADBase<Scalar> & parser);

} // namespace libMesh

#endif // LIBMESH_PARSED_FUNCTION_PROGRAM_H
