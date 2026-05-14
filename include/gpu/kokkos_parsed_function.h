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

#ifndef LIBMESH_KOKKOS_PARSED_FUNCTION_H
#define LIBMESH_KOKKOS_PARSED_FUNCTION_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_KOKKOS

#include "libmesh/libmesh_device.h"
#include "libmesh/parsed_function_program.h"
#include "libmesh/point.h"
#include "libmesh/vector_value.h"

#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__

#include <cmath>
#include <string>
#include <vector>

namespace libMesh::Kokkos
{
namespace detail
{

template <typename Scalar>
struct DeviceParsedFunctionProgram
{
  ::Kokkos::View<unsigned int *> bytecode;
  ::Kokkos::View<Scalar *> immediates;
  unsigned int stack_size = 0;
  unsigned int n_variables = 0;
  Scalar epsilon = 0;

  bool empty() const { return bytecode.extent(0) == 0; }
};

template <typename Scalar>
inline void
validate_program_stack(const DeviceParsedFunctionProgram<Scalar> & program,
                       const char * program_name,
                       const unsigned int max_stack)
{
  libmesh_error_msg_if(program.stack_size > max_stack,
                       "KokkosParsedFunction requires a larger MaxStack bound for " <<
                         program_name << " bytecode");
}

template <typename Scalar>
inline void
validate_coordinate_program_variables(const DeviceParsedFunctionProgram<Scalar> & program,
                                      const char * program_name)
{
  libmesh_error_msg_if(program.n_variables > LIBMESH_DIM + 1,
                       "KokkosParsedFunction currently supports only x/y/z/t variables in " <<
                         program_name << " bytecode");
}

template <typename T>
inline ::Kokkos::View<T *>
upload_scalar_buffer(const std::vector<T> & values,
                     const std::string & label)
{
  ::Kokkos::View<T *> d(label, values.size());
  auto h = ::Kokkos::create_mirror_view(d);

  for (std::size_t i = 0; i < values.size(); ++i)
    h(i) = values[i];

  ::Kokkos::deep_copy(d, h);
  return d;
}

template <typename Scalar>
inline DeviceParsedFunctionProgram<Scalar>
make_device_program(const libMesh::ParsedFunctionProgram<Scalar> & program,
                    const std::string & label)
{
  DeviceParsedFunctionProgram<Scalar> d_program;
  d_program.bytecode = upload_scalar_buffer(program.bytecode, label + "_bytecode");
  d_program.immediates = upload_scalar_buffer(program.immediates, label + "_immediates");
  d_program.stack_size = program.stack_size;
  d_program.n_variables = program.n_variables;
  d_program.epsilon = program.epsilon;
  return d_program;
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_abs(const Scalar x)
{
  using std::abs;
  return abs(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_floor(const Scalar x)
{
  using std::floor;
  return floor(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_ceil(const Scalar x)
{
  using std::ceil;
  return ceil(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_log(const Scalar x)
{
  using std::log;
  return log(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_log10(const Scalar x)
{
  using std::log10;
  return log10(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_log2(const Scalar x)
{
  using std::log2;
  return log2(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_sin(const Scalar x)
{
  using std::sin;
  return sin(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_cos(const Scalar x)
{
  using std::cos;
  return cos(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_tan(const Scalar x)
{
  using std::tan;
  return tan(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_sinh(const Scalar x)
{
  using std::sinh;
  return sinh(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_cosh(const Scalar x)
{
  using std::cosh;
  return cosh(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_tanh(const Scalar x)
{
  using std::tanh;
  return tanh(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_exp(const Scalar x)
{
  using std::exp;
  return exp(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_exp2(const Scalar x)
{
  using std::exp2;
  return exp2(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_sqrt(const Scalar x)
{
  using std::sqrt;
  return sqrt(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_pow(const Scalar x,
       const Scalar y)
{
  using std::pow;
  return pow(x, y);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_hypot(const Scalar x,
         const Scalar y)
{
  using std::hypot;
  return hypot(x, y);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_cbrt(const Scalar x)
{
  using std::cbrt;
  return cbrt(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_asin(const Scalar x)
{
  using std::asin;
  return asin(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_acos(const Scalar x)
{
  using std::acos;
  return acos(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_atan(const Scalar x)
{
  using std::atan;
  return atan(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_atan2(const Scalar y,
         const Scalar x)
{
  using std::atan2;
  return atan2(y, x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_asinh(const Scalar x)
{
  using std::asinh;
  return asinh(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_acosh(const Scalar x)
{
  using std::acosh;
  return acosh(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_atanh(const Scalar x)
{
  using std::atanh;
  return atanh(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_mod(const Scalar x,
       const Scalar y)
{
  using std::fmod;
  return fmod(x, y);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_trunc(const Scalar x)
{
  return x < Scalar(0) ? pf_ceil(x) : pf_floor(x);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE Scalar
pf_int(const Scalar x)
{
  return x < Scalar(0) ? pf_ceil(x - Scalar(0.5)) : pf_floor(x + Scalar(0.5));
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE bool
pf_equal(const Scalar x,
         const Scalar y,
         const Scalar epsilon)
{
  return pf_abs(x - y) <= epsilon;
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE bool
pf_nequal(const Scalar x,
          const Scalar y,
          const Scalar epsilon)
{
  return pf_abs(x - y) > epsilon;
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE bool
pf_less(const Scalar x,
        const Scalar y,
        const Scalar epsilon)
{
  return x < y - epsilon;
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE bool
pf_less_or_eq(const Scalar x,
              const Scalar y,
              const Scalar epsilon)
{
  return x <= y + epsilon;
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE bool
pf_truth(const Scalar x)
{
  return pf_abs(x) >= Scalar(0.5);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE bool
pf_abs_truth(const Scalar x)
{
  return x >= Scalar(0.5);
}

template <typename Scalar, unsigned int MaxStack>
LIBMESH_DEVICE_INLINE Scalar
eval_parsed_function_program(const DeviceParsedFunctionProgram<Scalar> & program,
                             const Scalar * vars)
{
  if (program.empty())
    return 0;

  Scalar stack[MaxStack];
  unsigned int dp = 0;
  int sp = -1;

  for (unsigned int ip = 0; ip < program.bytecode.extent(0); ++ip)
    {
      const unsigned int opcode = program.bytecode(ip);

      if (libMesh::parsed_function_is_var_opcode(opcode))
        {
          stack[++sp] = vars[opcode - libMesh::parsed_function_var_begin()];
          continue;
        }

      switch (static_cast<libMesh::ParsedFunctionOpcode>(opcode))
        {
        case libMesh::ParsedFunctionOpcode::cAbs: stack[sp] = pf_abs(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cAcos: stack[sp] = pf_acos(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cAcosh: stack[sp] = pf_acosh(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cAsin: stack[sp] = pf_asin(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cAsinh: stack[sp] = pf_asinh(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cAtan: stack[sp] = pf_atan(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cAtan2: stack[sp - 1] = pf_atan2(stack[sp - 1], stack[sp]); --sp; break;
        case libMesh::ParsedFunctionOpcode::cAtanh: stack[sp] = pf_atanh(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cCbrt: stack[sp] = pf_cbrt(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cCeil: stack[sp] = pf_ceil(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cCos: stack[sp] = pf_cos(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cCosh: stack[sp] = pf_cosh(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cCot: stack[sp] = Scalar(1) / pf_tan(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cCsc: stack[sp] = Scalar(1) / pf_sin(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cExp: stack[sp] = pf_exp(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cExp2: stack[sp] = pf_exp2(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cFloor: stack[sp] = pf_floor(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cHypot: stack[sp - 1] = pf_hypot(stack[sp - 1], stack[sp]); --sp; break;

        case libMesh::ParsedFunctionOpcode::cIf:
          if (pf_truth(stack[sp--]))
            ip += 2;
          else
            {
              const unsigned int jump_ip = program.bytecode(ip + 1);
              const unsigned int jump_dp = program.bytecode(ip + 2);
              ip = jump_ip;
              dp = jump_dp;
            }
          break;

        case libMesh::ParsedFunctionOpcode::cInt: stack[sp] = pf_int(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cLog: stack[sp] = pf_log(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cLog10: stack[sp] = pf_log10(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cLog2: stack[sp] = pf_log2(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cMax: stack[sp - 1] = stack[sp - 1] > stack[sp] ? stack[sp - 1] : stack[sp]; --sp; break;
        case libMesh::ParsedFunctionOpcode::cMin: stack[sp - 1] = stack[sp - 1] < stack[sp] ? stack[sp - 1] : stack[sp]; --sp; break;
        case libMesh::ParsedFunctionOpcode::cPow: stack[sp - 1] = pf_pow(stack[sp - 1], stack[sp]); --sp; break;
        case libMesh::ParsedFunctionOpcode::cSec: stack[sp] = Scalar(1) / pf_cos(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cSin: stack[sp] = pf_sin(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cSinh: stack[sp] = pf_sinh(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cSqrt: stack[sp] = pf_sqrt(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cTan: stack[sp] = pf_tan(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cTanh: stack[sp] = pf_tanh(stack[sp]); break;
        case libMesh::ParsedFunctionOpcode::cTrunc: stack[sp] = pf_trunc(stack[sp]); break;

        case libMesh::ParsedFunctionOpcode::cImmed: stack[++sp] = program.immediates(dp++); break;
        case libMesh::ParsedFunctionOpcode::cJump:
          ip = program.bytecode(ip + 1);
          dp = program.bytecode(ip + 2);
          break;

        case libMesh::ParsedFunctionOpcode::cNeg: stack[sp] = -stack[sp]; break;
        case libMesh::ParsedFunctionOpcode::cAdd: stack[sp - 1] += stack[sp]; --sp; break;
        case libMesh::ParsedFunctionOpcode::cSub: stack[sp - 1] -= stack[sp]; --sp; break;
        case libMesh::ParsedFunctionOpcode::cMul: stack[sp - 1] *= stack[sp]; --sp; break;
        case libMesh::ParsedFunctionOpcode::cDiv: stack[sp - 1] /= stack[sp]; --sp; break;
        case libMesh::ParsedFunctionOpcode::cMod: stack[sp - 1] = pf_mod(stack[sp - 1], stack[sp]); --sp; break;
        case libMesh::ParsedFunctionOpcode::cEqual: stack[sp - 1] = Scalar(pf_equal(stack[sp - 1], stack[sp], program.epsilon)); --sp; break;
        case libMesh::ParsedFunctionOpcode::cNEqual: stack[sp - 1] = Scalar(pf_nequal(stack[sp - 1], stack[sp], program.epsilon)); --sp; break;
        case libMesh::ParsedFunctionOpcode::cLess: stack[sp - 1] = Scalar(pf_less(stack[sp - 1], stack[sp], program.epsilon)); --sp; break;
        case libMesh::ParsedFunctionOpcode::cLessOrEq: stack[sp - 1] = Scalar(pf_less_or_eq(stack[sp - 1], stack[sp], program.epsilon)); --sp; break;
        case libMesh::ParsedFunctionOpcode::cGreater: stack[sp - 1] = Scalar(pf_less(stack[sp], stack[sp - 1], program.epsilon)); --sp; break;
        case libMesh::ParsedFunctionOpcode::cGreaterOrEq: stack[sp - 1] = Scalar(pf_less_or_eq(stack[sp], stack[sp - 1], program.epsilon)); --sp; break;
        case libMesh::ParsedFunctionOpcode::cNot: stack[sp] = Scalar(!pf_truth(stack[sp])); break;
        case libMesh::ParsedFunctionOpcode::cAnd: stack[sp - 1] = Scalar(pf_truth(stack[sp - 1]) && pf_truth(stack[sp])); --sp; break;
        case libMesh::ParsedFunctionOpcode::cOr: stack[sp - 1] = Scalar(pf_truth(stack[sp - 1]) || pf_truth(stack[sp])); --sp; break;
        case libMesh::ParsedFunctionOpcode::cNotNot: stack[sp] = Scalar(pf_truth(stack[sp])); break;

        case libMesh::ParsedFunctionOpcode::cDeg: stack[sp] = stack[sp] * Scalar(180.) / libMesh::pi; break;
        case libMesh::ParsedFunctionOpcode::cRad: stack[sp] = stack[sp] * libMesh::pi / Scalar(180.); break;

        case libMesh::ParsedFunctionOpcode::cPopNMov:
          {
            const unsigned int target = program.bytecode(++ip);
            const unsigned int source = program.bytecode(++ip);
            stack[target] = stack[source];
            sp = static_cast<int>(target);
            break;
          }

        case libMesh::ParsedFunctionOpcode::cLog2by:
          stack[sp - 1] = pf_log2(stack[sp - 1]) * stack[sp];
          --sp;
          break;

        case libMesh::ParsedFunctionOpcode::cNop:
          break;

        case libMesh::ParsedFunctionOpcode::cSinCos:
          stack[sp + 1] = pf_cos(stack[sp]);
          stack[sp] = pf_sin(stack[sp]);
          ++sp;
          break;

        case libMesh::ParsedFunctionOpcode::cSinhCosh:
          stack[sp + 1] = pf_cosh(stack[sp]);
          stack[sp] = pf_sinh(stack[sp]);
          ++sp;
          break;

        case libMesh::ParsedFunctionOpcode::cAbsNot: stack[sp] = Scalar(!pf_abs_truth(stack[sp])); break;
        case libMesh::ParsedFunctionOpcode::cAbsNotNot: stack[sp] = Scalar(pf_abs_truth(stack[sp])); break;
        case libMesh::ParsedFunctionOpcode::cAbsAnd: stack[sp - 1] = Scalar(pf_abs_truth(stack[sp - 1]) && pf_abs_truth(stack[sp])); --sp; break;
        case libMesh::ParsedFunctionOpcode::cAbsOr: stack[sp - 1] = Scalar(pf_abs_truth(stack[sp - 1]) || pf_abs_truth(stack[sp])); --sp; break;

        case libMesh::ParsedFunctionOpcode::cAbsIf:
          if (pf_abs_truth(stack[sp--]))
            ip += 2;
          else
            {
              const unsigned int jump_ip = program.bytecode(ip + 1);
              const unsigned int jump_dp = program.bytecode(ip + 2);
              ip = jump_ip;
              dp = jump_dp;
            }
          break;

        case libMesh::ParsedFunctionOpcode::cDup: stack[sp + 1] = stack[sp]; ++sp; break;

        case libMesh::ParsedFunctionOpcode::cFetch:
          {
            const unsigned int stack_offset = program.bytecode(++ip);
            stack[sp + 1] = stack[stack_offset];
            ++sp;
            break;
          }

        case libMesh::ParsedFunctionOpcode::cInv: stack[sp] = Scalar(1) / stack[sp]; break;
        case libMesh::ParsedFunctionOpcode::cSqr: stack[sp] = stack[sp] * stack[sp]; break;
        case libMesh::ParsedFunctionOpcode::cRDiv: stack[sp - 1] = stack[sp] / stack[sp - 1]; --sp; break;
        case libMesh::ParsedFunctionOpcode::cRSub: stack[sp - 1] = stack[sp] - stack[sp - 1]; --sp; break;
        case libMesh::ParsedFunctionOpcode::cRSqrt: stack[sp] = Scalar(1) / pf_sqrt(stack[sp]); break;

        default:
          return Scalar(0);
        }
    }

  return stack[sp];
}

template <typename Scalar, unsigned int MaxStack>
LIBMESH_DEVICE_INLINE Scalar
eval_coordinate_parsed_function_program(const DeviceParsedFunctionProgram<Scalar> & program,
                                        const Point & p,
                                        const Real time)
{
  Scalar vars[LIBMESH_DIM + 1];
  vars[0] = p(0);
#if LIBMESH_DIM > 1
  vars[1] = p(1);
#endif
#if LIBMESH_DIM > 2
  vars[2] = p(2);
#endif
  vars[LIBMESH_DIM] = time;
  return eval_parsed_function_program<Scalar, MaxStack>(program, vars);
}

} // namespace detail

template <typename Scalar = libMesh::Number, unsigned int MaxStack = 64>
class KokkosParsedFunction;

template <typename Scalar = libMesh::Number, unsigned int MaxStack = 64>
class KokkosParsedScalarProgram
{
public:
  LIBMESH_DEVICE_INLINE
  KokkosParsedScalarProgram() = default;

  explicit KokkosParsedScalarProgram(const libMesh::ParsedFunctionProgram<Scalar> & program,
                                     const std::string & label)
    : _program(detail::make_device_program(program, label))
  {
    detail::validate_program_stack(_program, label.c_str(), MaxStack);
  }

  LIBMESH_DEVICE_INLINE
  unsigned int n_variables() const
  {
    return _program.n_variables;
  }

  template <typename VariableStorage>
  LIBMESH_DEVICE_INLINE
  Scalar operator()(const VariableStorage & vars) const
  {
    return detail::eval_parsed_function_program<Scalar, MaxStack>(_program, vars);
  }

private:
  detail::DeviceParsedFunctionProgram<Scalar> _program;
};

template <typename Scalar = libMesh::Number, unsigned int MaxStack = 64>
class KokkosParsedGradient
{
public:
  LIBMESH_DEVICE_INLINE
  KokkosParsedGradient() = default;

  LIBMESH_DEVICE_INLINE
  explicit KokkosParsedGradient(const KokkosParsedFunction<Scalar, MaxStack> & func)
    : _func(func)
  {
  }

  LIBMESH_DEVICE_INLINE
  Gradient operator()(const Point & p) const;

private:
  KokkosParsedFunction<Scalar, MaxStack> _func;
};

template <typename Scalar, unsigned int MaxStack>
class KokkosParsedFunction
{
public:
  LIBMESH_DEVICE_INLINE
  KokkosParsedFunction() = default;

  explicit KokkosParsedFunction(const libMesh::ParsedFunctionProgramBundle<Scalar> & program_bundle,
                                const Real time = 0.)
    : _value(detail::make_device_program(program_bundle.value, "parsed_function_value")),
      _dx(detail::make_device_program(program_bundle.dx, "parsed_function_dx")),
#if LIBMESH_DIM > 1
      _dy(detail::make_device_program(program_bundle.dy, "parsed_function_dy")),
#endif
#if LIBMESH_DIM > 2
      _dz(detail::make_device_program(program_bundle.dz, "parsed_function_dz")),
#endif
      _dt(detail::make_device_program(program_bundle.dt, "parsed_function_dt")),
      _time(time)
  {
    detail::validate_program_stack(_value, "value", MaxStack);
    detail::validate_coordinate_program_variables(_value, "value");
    detail::validate_program_stack(_dx, "dx", MaxStack);
    detail::validate_coordinate_program_variables(_dx, "dx");
#if LIBMESH_DIM > 1
    detail::validate_program_stack(_dy, "dy", MaxStack);
    detail::validate_coordinate_program_variables(_dy, "dy");
#endif
#if LIBMESH_DIM > 2
    detail::validate_program_stack(_dz, "dz", MaxStack);
    detail::validate_coordinate_program_variables(_dz, "dz");
#endif
    detail::validate_program_stack(_dt, "dt", MaxStack);
    detail::validate_coordinate_program_variables(_dt, "dt");
  }

  KokkosParsedFunction
  with_time(const Real time) const
  {
    auto copy = *this;
    copy._time = time;
    return copy;
  }

  LIBMESH_DEVICE_INLINE
  Scalar operator()(const Point & p) const
  {
    return detail::eval_coordinate_parsed_function_program<Scalar, MaxStack>(_value, p, _time);
  }

  LIBMESH_DEVICE_INLINE
  Scalar time_derivative(const Point & p) const
  {
    return detail::eval_coordinate_parsed_function_program<Scalar, MaxStack>(_dt, p, _time);
  }

  LIBMESH_DEVICE_INLINE
  Gradient gradient(const Point & p) const
  {
    Gradient g;
    g(0) = detail::eval_coordinate_parsed_function_program<Scalar, MaxStack>(_dx, p, _time);
#if LIBMESH_DIM > 1
    g(1) = detail::eval_coordinate_parsed_function_program<Scalar, MaxStack>(_dy, p, _time);
#endif
#if LIBMESH_DIM > 2
    g(2) = detail::eval_coordinate_parsed_function_program<Scalar, MaxStack>(_dz, p, _time);
#endif
    return g;
  }

  LIBMESH_DEVICE_INLINE
  KokkosParsedGradient<Scalar, MaxStack> gradient_function() const
  {
    return KokkosParsedGradient<Scalar, MaxStack>(*this);
  }

private:
  detail::DeviceParsedFunctionProgram<Scalar> _value;
  detail::DeviceParsedFunctionProgram<Scalar> _dx;
#if LIBMESH_DIM > 1
  detail::DeviceParsedFunctionProgram<Scalar> _dy;
#endif
#if LIBMESH_DIM > 2
  detail::DeviceParsedFunctionProgram<Scalar> _dz;
#endif
  detail::DeviceParsedFunctionProgram<Scalar> _dt;
  Real _time = 0.;

  friend class KokkosParsedGradient<Scalar, MaxStack>;
};

template <typename Scalar, unsigned int MaxStack>
LIBMESH_DEVICE_INLINE
Gradient
KokkosParsedGradient<Scalar, MaxStack>::operator()(const Point & p) const
{
  return _func.gradient(p);
}

template <typename Scalar = libMesh::Number,
          unsigned int MaxStack = 64,
          unsigned int MaxFieldVariables = 16>
class KokkosParsedFEMFunction
{
public:
  LIBMESH_DEVICE_INLINE
  KokkosParsedFEMFunction() = default;

  explicit KokkosParsedFEMFunction(const libMesh::ParsedFEMFunctionProgramBundle<Scalar> & program_bundle,
                                   const Real time = 0.)
    : _value(program_bundle.value, "parsed_fem_function_value"),
      _dx(program_bundle.dx, "parsed_fem_function_dx"),
#if LIBMESH_DIM > 1
      _dy(program_bundle.dy, "parsed_fem_function_dy"),
#endif
#if LIBMESH_DIM > 2
      _dz(program_bundle.dz, "parsed_fem_function_dz"),
#endif
      _dt(program_bundle.dt, "parsed_fem_function_dt"),
      _n_field_variables(cast_int<unsigned int>(program_bundle.value_variable_numbers.size())),
      _time(time)
  {
    libmesh_error_msg_if(!program_bundle.supports_kokkos_value_goal(),
                         "KokkosParsedFEMFunction currently supports only value-based ParsedFEMFunction expressions");
    libmesh_error_msg_if(_n_field_variables > MaxFieldVariables,
                         "KokkosParsedFEMFunction exceeds MaxFieldVariables");

    for (unsigned int i = 0; i != _n_field_variables; ++i)
      {
        _field_variable_numbers[i] = program_bundle.value_variable_numbers[i];
        _field_value_derivatives[i] =
          KokkosParsedScalarProgram<Scalar, MaxStack>(
            program_bundle.value_variable_derivatives[i],
            "parsed_fem_function_dvalue_" + std::to_string(i));
      }
  }

  LIBMESH_DEVICE_INLINE
  KokkosParsedFEMFunction
  with_time(const Real time) const
  {
    auto copy = *this;
    copy._time = time;
    return copy;
  }

  LIBMESH_DEVICE_INLINE
  unsigned int n_field_variables() const
  {
    return _n_field_variables;
  }

  LIBMESH_DEVICE_INLINE
  Real time() const
  {
    return _time;
  }

  LIBMESH_DEVICE_INLINE
  unsigned int field_variable_number(const unsigned int i) const
  {
    return _field_variable_numbers[i];
  }

  template <typename VariableStorage>
  LIBMESH_DEVICE_INLINE
  Scalar value(const VariableStorage & vars) const
  {
    return _value(vars);
  }

  template <typename VariableStorage>
  LIBMESH_DEVICE_INLINE
  Gradient gradient(const VariableStorage & vars,
                    const Gradient * field_gradients) const
  {
    Gradient g;
    g(0) = _dx(vars);
#if LIBMESH_DIM > 1
    g(1) = _dy(vars);
#endif
#if LIBMESH_DIM > 2
    g(2) = _dz(vars);
#endif

    for (unsigned int i = 0; i != _n_field_variables; ++i)
      g.add_scaled(field_gradients[i], _field_value_derivatives[i](vars));

    return g;
  }

private:
  KokkosParsedScalarProgram<Scalar, MaxStack> _value;
  KokkosParsedScalarProgram<Scalar, MaxStack> _dx;
#if LIBMESH_DIM > 1
  KokkosParsedScalarProgram<Scalar, MaxStack> _dy;
#endif
#if LIBMESH_DIM > 2
  KokkosParsedScalarProgram<Scalar, MaxStack> _dz;
#endif
  KokkosParsedScalarProgram<Scalar, MaxStack> _dt;
  KokkosParsedScalarProgram<Scalar, MaxStack> _field_value_derivatives[MaxFieldVariables];
  unsigned int _field_variable_numbers[MaxFieldVariables] = {};
  unsigned int _n_field_variables = 0;
  Real _time = 0.;
};

} // namespace libMesh::Kokkos

#endif // LIBMESH_HAVE_KOKKOS

#endif // LIBMESH_KOKKOS_PARSED_FUNCTION_H
