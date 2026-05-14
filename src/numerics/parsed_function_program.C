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

#include "libmesh/parsed_function_program.h"

#include "fparser_ad.hh"
#include "extrasrc/fptypes.hh"

namespace
{

template <typename Scalar>
void
validate_kokkos_program_opcode(const unsigned int opcode)
{
  using libMesh::ParsedFunctionOpcode;

  if (libMesh::parsed_function_is_var_opcode(opcode))
    return;

  switch (static_cast<ParsedFunctionOpcode>(opcode))
    {
    case ParsedFunctionOpcode::cAbs:
    case ParsedFunctionOpcode::cAcos:
    case ParsedFunctionOpcode::cAcosh:
    case ParsedFunctionOpcode::cAsin:
    case ParsedFunctionOpcode::cAsinh:
    case ParsedFunctionOpcode::cAtan:
    case ParsedFunctionOpcode::cAtan2:
    case ParsedFunctionOpcode::cAtanh:
    case ParsedFunctionOpcode::cCbrt:
    case ParsedFunctionOpcode::cCeil:
    case ParsedFunctionOpcode::cCos:
    case ParsedFunctionOpcode::cCosh:
    case ParsedFunctionOpcode::cCot:
    case ParsedFunctionOpcode::cCsc:
    case ParsedFunctionOpcode::cExp:
    case ParsedFunctionOpcode::cExp2:
    case ParsedFunctionOpcode::cFloor:
    case ParsedFunctionOpcode::cHypot:
    case ParsedFunctionOpcode::cIf:
    case ParsedFunctionOpcode::cInt:
    case ParsedFunctionOpcode::cLog:
    case ParsedFunctionOpcode::cLog10:
    case ParsedFunctionOpcode::cLog2:
    case ParsedFunctionOpcode::cMax:
    case ParsedFunctionOpcode::cMin:
    case ParsedFunctionOpcode::cPow:
    case ParsedFunctionOpcode::cSec:
    case ParsedFunctionOpcode::cSin:
    case ParsedFunctionOpcode::cSinh:
    case ParsedFunctionOpcode::cSqrt:
    case ParsedFunctionOpcode::cTan:
    case ParsedFunctionOpcode::cTanh:
    case ParsedFunctionOpcode::cTrunc:
    case ParsedFunctionOpcode::cImmed:
    case ParsedFunctionOpcode::cJump:
    case ParsedFunctionOpcode::cNeg:
    case ParsedFunctionOpcode::cAdd:
    case ParsedFunctionOpcode::cSub:
    case ParsedFunctionOpcode::cMul:
    case ParsedFunctionOpcode::cDiv:
    case ParsedFunctionOpcode::cMod:
    case ParsedFunctionOpcode::cEqual:
    case ParsedFunctionOpcode::cNEqual:
    case ParsedFunctionOpcode::cLess:
    case ParsedFunctionOpcode::cLessOrEq:
    case ParsedFunctionOpcode::cGreater:
    case ParsedFunctionOpcode::cGreaterOrEq:
    case ParsedFunctionOpcode::cNot:
    case ParsedFunctionOpcode::cAnd:
    case ParsedFunctionOpcode::cOr:
    case ParsedFunctionOpcode::cNotNot:
    case ParsedFunctionOpcode::cDeg:
    case ParsedFunctionOpcode::cRad:
    case ParsedFunctionOpcode::cPopNMov:
    case ParsedFunctionOpcode::cLog2by:
    case ParsedFunctionOpcode::cNop:
    case ParsedFunctionOpcode::cSinCos:
    case ParsedFunctionOpcode::cSinhCosh:
    case ParsedFunctionOpcode::cAbsAnd:
    case ParsedFunctionOpcode::cAbsOr:
    case ParsedFunctionOpcode::cAbsNot:
    case ParsedFunctionOpcode::cAbsNotNot:
    case ParsedFunctionOpcode::cAbsIf:
    case ParsedFunctionOpcode::cDup:
    case ParsedFunctionOpcode::cFetch:
    case ParsedFunctionOpcode::cInv:
    case ParsedFunctionOpcode::cSqr:
    case ParsedFunctionOpcode::cRDiv:
    case ParsedFunctionOpcode::cRSub:
    case ParsedFunctionOpcode::cRSqrt:
      return;

    case ParsedFunctionOpcode::cArg:
    case ParsedFunctionOpcode::cConj:
    case ParsedFunctionOpcode::cImag:
    case ParsedFunctionOpcode::cPolar:
    case ParsedFunctionOpcode::cReal:
      libmesh_error_msg("Kokkos parsed-function export does not support complex-valued fparser opcodes");

    case ParsedFunctionOpcode::cFCall:
    case ParsedFunctionOpcode::cPCall:
      libmesh_error_msg("Kokkos parsed-function export does not support user-defined or nested parser calls");

    case ParsedFunctionOpcode::VarBegin:
      return;
    }

  libmesh_error_msg("Kokkos parsed-function export encountered an unknown opcode " << opcode);
}

} // anonymous namespace

namespace libMesh
{

template <typename Scalar>
ParsedFunctionProgram<Scalar>
build_parsed_function_program(const FunctionParserADBase<Scalar> & parser)
{
  ParsedFunctionProgram<Scalar> program;
  const auto * data = parser.parser_data();
  libmesh_assert(data);

  program.bytecode.assign(data->mByteCode.begin(), data->mByteCode.end());
  program.immediates.assign(data->mImmed.begin(), data->mImmed.end());
  program.stack_size = data->mStackSize;
  program.n_variables = data->mVariablesAmount;
  program.epsilon = FunctionParserBase<Scalar>::epsilon();

  for (const auto opcode : program.bytecode)
    validate_kokkos_program_opcode<Scalar>(opcode);

  return program;
}

template ParsedFunctionProgram<Real>
build_parsed_function_program(const FunctionParserADBase<Real> & parser);

} // namespace libMesh
