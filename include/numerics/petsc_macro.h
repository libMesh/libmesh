// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PETSC_MACRO_H
#define LIBMESH_PETSC_MACRO_H

// Local includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// A convenient macro for comparing PETSc versions.  Returns 1 if the
// current PETSc version is < major.minor.subminor and zero otherwise.
//
// This macro does not require petscversion.h to be included for it to work correctly.
// It instead relies on the PETSc version numbers detected during configure.  Note that if
// LIBMESH_HAVE_PETSC is not defined, none of the LIBMESH_DETECTED_PETSC_VERSION_* variables will
// be defined either.
#define PETSC_VERSION_LESS_THAN(major,minor,subminor)                   \
  ((LIBMESH_DETECTED_PETSC_VERSION_MAJOR < (major) ||                   \
    (LIBMESH_DETECTED_PETSC_VERSION_MAJOR == (major) && (LIBMESH_DETECTED_PETSC_VERSION_MINOR < (minor) || \
                                                         (LIBMESH_DETECTED_PETSC_VERSION_MINOR == (minor) && \
                                                          LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR < (subminor))))))
#define PETSC_VERSION_EQUALS(major,minor,subminor)              \
  ((LIBMESH_DETECTED_PETSC_VERSION_MAJOR == (major)) &&         \
   (LIBMESH_DETECTED_PETSC_VERSION_MINOR == (minor)) &&         \
   (LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR == (subminor)))

// We used to have workarounds for missing extern "C" in old PETSc
// versions.  We no longer support PETSc versions so old, but we do
// still support libMesh applications old enough to have used these
// macros.
#define EXTERN_C_FOR_PETSC_BEGIN
#define EXTERN_C_FOR_PETSC_END

// Petsc include files
// Wrapped to avoid triggering our more paranoid warnings
#include <libmesh/ignore_warnings.h>
#ifdef I
# define LIBMESH_SAW_I
#endif
#include <petsc.h>
#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif
#include <libmesh/restore_warnings.h>

// These macros are still around for backwards compatibility, but we
// no longer use them within the library.
#define LibMeshVecDestroy(x)         VecDestroy(x)
#define LibMeshVecScatterDestroy(x)  VecScatterDestroy(x)
#define LibMeshMatDestroy(x)         MatDestroy(x)
#define LibMeshISDestroy(x)          ISDestroy(x)
#define LibMeshKSPDestroy(x)         KSPDestroy(x)
#define LibMeshSNESDestroy(x)        SNESDestroy(x)
#define LibMeshPetscViewerDestroy(x) PetscViewerDestroy(x)
#define LibMeshPCDestroy(x)          PCDestroy(x)

// PETSc devs temporarily considered adding VecScatterCreateWithData, but
// it was dropped in PETSc-130e142e39 and never made it into any release.
// We will keep the ifdef for backwards compatibility in case anyone wrote
// code directly using it, but that should be pretty unlikely.
#define LibMeshVecScatterCreate(xin,ix,yin,iy,newctx)  VecScatterCreate(xin,ix,yin,iy,newctx)
#define ISCreateLibMesh(comm,n,idx,mode,is) ISCreateGeneral((comm),(n),(idx),(mode),(is))

// As of release 3.8.0, MatGetSubMatrix was renamed to MatCreateSubMatrix.
#if PETSC_VERSION_LESS_THAN(3,8,0)
# define LibMeshCreateSubMatrix MatGetSubMatrix
#else
# define LibMeshCreateSubMatrix MatCreateSubMatrix
#endif

// As of release 3.17.0, SETERRQ was made to use __VA_ARGS__ and
// SETERRQ1, SETERRQ2, etc. were deprecated
#if PETSC_VERSION_LESS_THAN(3,17,0)
# define LIBMESH_SETERRQ1 SETERRQ1
# define LIBMESH_SETERRQ2 SETERRQ2
# define LIBMESH_SETERRQ3 SETERRQ3
#else
# define LIBMESH_SETERRQ1 SETERRQ
# define LIBMESH_SETERRQ2 SETERRQ
# define LIBMESH_SETERRQ3 SETERRQ
#endif

// As of 3.18, %D is no longer supported in format strings, but the
// replacement PetscInt_FMT didn't get added until 3.7.2
#if PETSC_VERSION_LESS_THAN(3,8,0)
# define LIBMESH_PETSCINT_FMT "D"
#else
# define LIBMESH_PETSCINT_FMT PetscInt_FMT
#endif

// As of 3.19, PETSC_NULL is deprecated
#if PETSC_VERSION_LESS_THAN(3,19,0)
# define LIBMESH_PETSC_NULLPTR PETSC_NULL
# define LIBMESH_PETSC_SUCCESS static_cast<PetscErrorCode>(0)
#else
# define LIBMESH_PETSC_NULLPTR PETSC_NULLPTR
# define LIBMESH_PETSC_SUCCESS PETSC_SUCCESS
#endif

// If we're using quad precision, we need to disambiguate std
// operations on PetscScalar

#if LIBMESH_DEFAULT_QUADRUPLE_PRECISION
# include <boost/multiprecision/float128.hpp>

namespace std
{
inline
std::ostream & operator<< (std::ostream & os, const PetscScalar in)
{
  os << (boost::multiprecision::float128(in));
  return os;
}

/*

// Whether these are necessary or not depends on gcc version!?

#define LIBMESH_PETSCSCALAR_UNARY(funcname) \
inline PetscScalar funcname \
  (const PetscScalar in) \
{ \
  return boost::multiprecision::funcname \
    (boost::multiprecision::float128(in)).backend().value(); \
}

LIBMESH_PETSCSCALAR_UNARY(sqrt)
LIBMESH_PETSCSCALAR_UNARY(exp)
LIBMESH_PETSCSCALAR_UNARY(log)
LIBMESH_PETSCSCALAR_UNARY(log10)
LIBMESH_PETSCSCALAR_UNARY(sin)
LIBMESH_PETSCSCALAR_UNARY(cos)
LIBMESH_PETSCSCALAR_UNARY(tan)
LIBMESH_PETSCSCALAR_UNARY(asin)
LIBMESH_PETSCSCALAR_UNARY(acos)
LIBMESH_PETSCSCALAR_UNARY(atan)
LIBMESH_PETSCSCALAR_UNARY(sinh)
LIBMESH_PETSCSCALAR_UNARY(cosh)
LIBMESH_PETSCSCALAR_UNARY(tanh)
LIBMESH_PETSCSCALAR_UNARY(abs)
LIBMESH_PETSCSCALAR_UNARY(fabs)
LIBMESH_PETSCSCALAR_UNARY(ceil)
LIBMESH_PETSCSCALAR_UNARY(floor)
*/

} // namespace std

// Helper functions for boost float128 compatibility
namespace libMesh
{
template <typename T>
PetscScalar PS(T val)
{
  return val.backend().value();
}

template <typename T>
PetscScalar * pPS(T * ptr)
{
  return &(ptr->backend().value());
}

template <typename T>
const PetscScalar * pPS(const T * ptr)
{
  return &(ptr->backend().value());
}

template <typename T>
PetscReal * pPR(T * ptr)
{
  return &(ptr->backend().value());
}

template <typename T>
const PetscReal * pPR(const T * ptr)
{
  return &(ptr->backend().value());
}
} // namespace libMesh

#else

namespace libMesh
{
template <typename T>
PetscScalar PS(T val)
{
  return val;
}

template <typename T>
PetscScalar * pPS(T * ptr)
{
  return ptr;
}

template <typename T>
const PetscScalar * pPS(const T * ptr)
{
  return ptr;
}

template <typename T>
PetscReal * pPR(T * ptr)
{
  return ptr;
}

template <typename T>
const PetscReal * pPR(const T * ptr)
{
  return ptr;
}
} // namespace libMesh

#endif // LIBMESH_ENABLE_QUADRUPLE_PRECISION

#else // LIBMESH_HAVE_PETSC

#define PETSC_VERSION_LESS_THAN(major,minor,subminor) 1
#define PETSC_VERSION_EQUALS(major,minor,revision) 0

#endif // LIBMESH_HAVE_PETSC

// The PETSC_VERSION_RELEASE constant was introduced just prior to 2.3.0 (ca. Apr 22 2005),
// so fall back to using PETSC_VERSION_LESS_THAN in case it doesn't exist.
#ifdef LIBMESH_DETECTED_PETSC_VERSION_RELEASE

#define PETSC_RELEASE_LESS_THAN(major, minor, subminor) \
  (PETSC_VERSION_LESS_THAN(major, minor, subminor) && LIBMESH_DETECTED_PETSC_VERSION_RELEASE)
#define PETSC_RELEASE_EQUALS(major, minor, subminor) \
  (PETSC_VERSION_EQUALS(major, minor, subminor) && LIBMESH_DETECTED_PETSC_VERSION_RELEASE)

#else

#define PETSC_RELEASE_LESS_THAN(major, minor, subminor) \
  (PETSC_VERSION_LESS_THAN(major, minor, subminor))
#define PETSC_RELEASE_EQUALS(major, minor, subminor) (PETSC_VERSION_EQUALS(major, minor, subminor))

#endif

#define PETSC_RELEASE_LESS_EQUALS(major, minor, subminor) \
  (PETSC_RELEASE_LESS_THAN(major, minor, subminor) || PETSC_RELEASE_EQUALS(major, minor, subminor))

#define PETSC_RELEASE_GREATER_EQUALS(major, minor, subminor) \
  (0 == PETSC_RELEASE_LESS_THAN(major, minor, subminor))

#define PETSC_RELEASE_GREATER_THAN(major, minor, subminor) \
  (0 == PETSC_RELEASE_LESS_EQUALS(major, minor, subminor))

#endif // LIBMESH_PETSC_MACRO_H
