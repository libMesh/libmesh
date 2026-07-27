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


// Note: no include guards!  We want to be able to #include this
// header multiple times.

// TODO: icpc options
#if defined(__clang__) || defined(__NVCOMPILER)
#  pragma float_control(push)
#  pragma float_control(precise, on)
#  pragma clang fp contract(off) reassociate(off)
#elif defined(__GNUC__)
#  pragma GCC push_options
#  pragma GCC optimize("-fno-unsafe-math-optimizations")
#  pragma GCC optimize("-ffp-contract=off")
#endif

// If we start using these headers for C rather than just C++ then we
// might add #pragma STDC lines here
