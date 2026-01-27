// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local includes
#include "libmesh/libmesh_exceptions.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/perf_log.h"
#include "libmesh/print_trace.h"

// C/C++ includes
#include <cstdlib>

#ifdef LIBMESH_ENABLE_EXCEPTIONS
#include <exception>
#include <optional>
#endif

#ifdef LIBMESH_HAVE_OPENMP
#include <omp.h>
#endif

#include "signal.h"


// floating-point exceptions
#ifdef LIBMESH_HAVE_FENV_H
#  include <fenv.h>
#endif
#ifdef LIBMESH_HAVE_XMMINTRIN_H
#  include <xmmintrin.h>
#endif


#if defined(LIBMESH_HAVE_MPI)
# include "libmesh/ignore_warnings.h"
# include <mpi.h>
# include "libmesh/restore_warnings.h"
#endif // #if defined(LIBMESH_HAVE_MPI)

#if defined(LIBMESH_HAVE_PETSC)
# include "libmesh/petsc_solver_exception.h"
# include <petsc.h>
# include <petscerror.h>
# include "libmesh/petscdmlibmesh.h"
# if defined(LIBMESH_HAVE_SLEPC)
// Ignore unused variable warnings from SLEPc
#  include "libmesh/ignore_warnings.h"
#  include "libmesh/slepc_macro.h"
#  include <slepc.h>
#  include "libmesh/restore_warnings.h"
# endif // #if defined(LIBMESH_HAVE_SLEPC)
#endif // #if defined(LIBMESH_HAVE_PETSC)

#ifdef LIBMESH_HAVE_NETGEN
// We need the nglib namespace, because it's used everywhere in nglib
// and we don't get binary compatibility without it.
//
// We need the nglib namespace *here*, because somehow nobody ever
// figured out to just put it in nglib.h?
namespace nglib {
#include "netgen/nglib/nglib.h"
}
#endif

// If we're using MPI and VTK has been detected, we need to do some
// MPI initialize/finalize stuff for VTK.
#if defined(LIBMESH_HAVE_MPI) && defined(LIBMESH_HAVE_VTK)
#include "libmesh/ignore_warnings.h"
# include "vtkMPIController.h"
#include "libmesh/restore_warnings.h"
#endif

#include <mutex>

// --------------------------------------------------------
// Local anonymous namespace to hold miscellaneous bits
namespace {

/**
 * Floating point exception handler -- courtesy of Cody Permann & MOOSE team
 */
#if LIBMESH_HAVE_DECL_SIGACTION
void libmesh_handleFPE(int /*signo*/, siginfo_t * info, void * /*context*/)
{
  libMesh::err << std::endl;
  libMesh::err << "Floating point exception signaled (";
  switch (info->si_code)
    {
    case FPE_INTDIV: libMesh::err << "integer divide by zero"; break;
    case FPE_INTOVF: libMesh::err << "integer overflow"; break;
    case FPE_FLTDIV: libMesh::err << "floating point divide by zero"; break;
    case FPE_FLTOVF: libMesh::err << "floating point overflow"; break;
    case FPE_FLTUND: libMesh::err << "floating point underflow"; break;
    case FPE_FLTRES: libMesh::err << "floating point inexact result"; break;
    case FPE_FLTINV: libMesh::err << "invalid floating point operation"; break;
    case FPE_FLTSUB: libMesh::err << "subscript out of range"; break;
    default:         libMesh::err << "unrecognized"; break;
    }
  libMesh::err << ")!" << std::endl;

  libmesh_error_msg("\nTo track this down, compile in debug mode, then in gdb do:\n" \
                    << "  break libmesh_handleFPE\n"                    \
                    << "  run ...\n"                                    \
                    << "  bt");
}


void libmesh_handleSEGV(int /*signo*/, siginfo_t * info, void * /*context*/)
{
  libMesh::err << std::endl;
  libMesh::err << "Segmentation fault exception signaled (";
  switch (info->si_code)
    {
    case SEGV_MAPERR: libMesh::err << "Address not mapped"; break;
    case SEGV_ACCERR: libMesh::err << "Invalid permissions"; break;
    default:         libMesh::err << "unrecognized"; break;
    }
  libMesh::err << ")!" << std::endl;

  libmesh_error_msg("\nTo track this down, compile in debug mode, then in gdb do:\n" \
                    << "  break libmesh_handleSEGV\n"                    \
                    << "  run ...\n"                                    \
                    << "  bt");
}
#endif
} // anonymous namespace



namespace libMesh
{

// ------------------------------------------------------------
// libMesh functions

void libmesh_terminate_handler()
{
  bool quiet = false;

#ifdef LIBMESH_ENABLE_EXCEPTIONS
  // If we have an active exception, it may have an error message that
  // we should print, or it may have a type that tells us not to print
  // anything.
  std::optional<std::string> exception_message;
  std::exception_ptr ex = std::current_exception();
  if (ex)
    {
      try
        {
          std::rethrow_exception(ex);
        }
      // Capture the exception message to be used later.
      catch (const std::exception & std_ex)
        {
          exception_message = std_ex.what();
        }
      // We arrived here via TerminationException (likely from
      // libmesh_terminate()), which implies that a useful
      // error message has already been emitted.
      catch (const TerminationException &)
        {
          quiet = true;
        }
      // We're just trying to detect exception types here, not
      // actually rethrow
      catch (...)
        {
        }
    }
#endif

  if (!quiet)
    {
      libMesh::err << "libMesh terminating";
#ifdef LIBMESH_ENABLE_EXCEPTIONS
      if (exception_message)
        libMesh::err << ":\n" << *exception_message;
#endif
      libMesh::err << std::endl;

      // If this got called then we're probably crashing; let's print a
      // stack trace.  The trace files that are ultimately written depend on:
      // 1.) Who throws the exception.
      // 2.) Whether the C++ runtime unwinds the stack before the
      //     terminate_handler is called (this is implementation defined).
      //
      // The various cases are summarized in the table below:
      //
      //                        | libmesh exception | other exception
      //                        -------------------------------------
      // stack unwinds          |        A          |       B
      // stack does not unwind  |        C          |       D
      //
      // Case A: There will be two stack traces in the file: one "useful"
      //         one, and one nearly empty one due to stack unwinding.
      // Case B: You will get one nearly empty stack trace (not great, Bob!)
      // Case C: You will get two nearly identical stack traces, ignore one of them.
      // Case D: You will get one useful stack trace.
      //
      // Cases A and B (where the stack unwinds when an exception leaves
      // main) appear to be non-existent in practice.  I don't have a
      // definitive list, but the stack does not unwind for GCC on either
      // Mac or Linux.  I think there's good reasons for this behavior too:
      // it's much easier to get a stack trace when the stack doesn't
      // unwind, for example.
      libMesh::write_traceout();

      // We may care about performance data pre-crash; it would be sad to
      // throw that away.
      LibMeshInit::perf_log().print_log();
    }

  libmesh_abort();
}


/**
 * Toggle floating point exceptions -- courtesy of Cody Permann & MOOSE team
 */
void enableFPE(bool on)
{
#if !defined(LIBMESH_HAVE_FEENABLEEXCEPT) && defined(LIBMESH_HAVE_XMMINTRIN_H)
  static int flags = 0;
#endif

  if (on)
    {
#ifdef LIBMESH_HAVE_FEENABLEEXCEPT
      feenableexcept(FE_DIVBYZERO | FE_INVALID);
#elif  LIBMESH_HAVE_XMMINTRIN_H
      flags = _MM_GET_EXCEPTION_MASK();           // store the flags
      _MM_SET_EXCEPTION_MASK(flags & ~_MM_MASK_INVALID);
#endif

#if LIBMESH_HAVE_DECL_SIGACTION
      struct sigaction new_action, old_action;

      // Set up the structure to specify the new action.
      new_action.sa_sigaction = libmesh_handleFPE;
      sigemptyset (&new_action.sa_mask);
      new_action.sa_flags = SA_SIGINFO;

      sigaction (SIGFPE, nullptr, &old_action);
      if (old_action.sa_handler != SIG_IGN)
        sigaction (SIGFPE, &new_action, nullptr);
#endif
    }
  else
    {
#ifdef LIBMESH_HAVE_FEDISABLEEXCEPT
      fedisableexcept(FE_DIVBYZERO | FE_INVALID);
#elif  LIBMESH_HAVE_XMMINTRIN_H
      _MM_SET_EXCEPTION_MASK(flags);
#endif
      signal(SIGFPE, SIG_DFL);
    }
}


// Enable handling of SIGSEGV by libMesh
// (potentially instead of PETSc)
void enableSEGV(bool on)
{
#if LIBMESH_HAVE_DECL_SIGACTION
  static struct sigaction old_action;
  static bool was_on = false;

  if (on)
    {
      struct sigaction new_action;
      was_on = true;

      // Set up the structure to specify the new action.
      new_action.sa_sigaction = libmesh_handleSEGV;
      sigemptyset (&new_action.sa_mask);
      new_action.sa_flags = SA_SIGINFO;

      sigaction (SIGSEGV, &new_action, &old_action);
    }
  else if (was_on)
    {
      was_on = false;
      sigaction (SIGSEGV, &old_action, nullptr);
    }
#else
  libmesh_error_msg("System call sigaction not supported.");
#endif
}

} // namespace libMesh
