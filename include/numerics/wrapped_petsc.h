// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_WRAPPED_PETSC_H
#define LIBMESH_WRAPPED_PETSC_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// C++ includes
#include <utility> // std::swap, std::move

namespace libMesh
{

// Template class which wraps PETSc objects and specializes the
// destructor to call the appropriate "XXXDestroy()" routine.
template <typename T>
struct WrappedPetsc
{
  /**
   * Default constructor. This should mimic the way that we normally
   * write e.g.
   * KSP ksp;
   * and then proceed to use ksp in different PETSc routines. That is,
   * obj is not initialized to any particular value.
   */
  WrappedPetsc() : obj() {}

  /**
   * Constructor which initializes obj to a specific passed-in value.
   * This mimics code in which we explicitly create a PETSc object via
   * IS is = NULL;
   * Technically one could pass any pointer value in here, but it usually
   * only makes sense to pass nullptr.
   */
  WrappedPetsc(T obj_in) : obj(obj_in) {}

  /**
   * Destructor. Just calls destroy().
   */
  ~WrappedPetsc()
  {
    destroy();
  }

  /**
   * Calls destroy() _and_ sets the managed object to nullptr. As far
   * as I can tell, setting obj to nullptr is not done by the various
   * XXXDestroy() routines of PETSc, so we also don't do this in the
   * wrapping class's destructor, however, there are situations where
   * it is sometimes useful to both call the relevant XXXDestroy()
   * function and reset the pointer, hence the need for this function.
   */
  void reset_to_zero()
  {
    destroy();
    obj = nullptr;
  }

  /**
   * Copy constructor and copy assignment operator. These are deleted
   * since I don't think we can safely shallow copy PETSc objects like
   * KSP and Vec, which are internally reference-counted pointers and
   * probably don't do the right thing if they are shallow-copied.
   */
  WrappedPetsc(const WrappedPetsc & other) = delete;
  WrappedPetsc & operator= (const WrappedPetsc &) = delete;

  /**
   * Move constructor. We could almost default this, but we need to
   * set other.obj to nullptr so that when it is subsequently
   * Destroy()ed it's just a no-op rather than messing up the
   * reference count or trying to double-free memory.
   */
  WrappedPetsc(WrappedPetsc && other) noexcept
    : obj(other.obj)
  {
    other.obj = nullptr;
  }

  /**
   * Move-assignment operator. Use move-construct-and-swap idiom
   * instead of defaulting since we want to make sure our move
   * constructor leaves the passed-in object in a Destroy()able state.
   */
  WrappedPetsc & operator= (WrappedPetsc && other) noexcept
  {
    WrappedPetsc tmp(std::move(other));
    std::swap(tmp, *this);
    return *this;
  }

  /**
   * \returns pointer to the managed object
   * This is used to mimic code such as:
   * KSP ksp;
   * KSPCreate(comm, &ksp);
   * Since taking the address of the wrapping object doesn't make
   * sense in this context.
   */
  T * get() { return &obj; }

  /**
   * User-defined conversion function. We provide non-const access to
   * the underlying T object even when the "this" object is considered
   * const, since PETSc APIs which are "logically const" typically
   * still take non-const parameters.
   */
  operator T() const { return obj; }

  /**
   * The "dereferencing" operator. Returns a reference to the managed
   * object. This is needed for some situations in which the
   * user-defined conversion operator doesn't work, for example with
   * C-style casts:
   * KSP ksp;
   * ...
   * PetscObjectSetOptionsPrefix((PetscObject)(*ksp), "balance_");
   */
  T & operator*() { return obj; }

  /**
   * User-defined conversion to bool. This is intended to mimic code like:
   * IS is = nullptr;
   * ...
   * if (!is)
   *   ...
   * Note that this comparison is thus concerned with obj itself and
   * not &obj.
   */
  operator bool() const { return obj != nullptr; }

  /**
   * Must be specialized to call the appropriate XXXDestroy() routine
   * in order for a WrappedPetsc<T> object to be instantiated.
   *
   * We could try to do extra error checking in destroy() as shown
   * below, but note that:
   * 1.) destroy() is called from destructors, sometimes during
   *     stack unwinding. If there's an error code returned from
   *     XXXDestroy(), then our only option is to immediately terminate
   *     the program, which would then kill any chance of recovering from
   *     the exception.
   * 2.) It's not always safe to call non-Destroy() functions on PETSc objects
   *     which are about to be Destroy()ed, that is, we would have to check
   *     for nullptr, etc. which would lead to more complexity, and more code.
   *
   * One possible approach for extra error checking with immediate
   * abort on error:
   *
   * MPI_Comm comm;
   * PetscObjectGetComm((PetscObject)(&obj), &comm);
   * PetscErrorCode ierr = Type ## Destroy(&obj);
   * CHKERRABORT(comm, ierr);
   */
  void destroy();

private:
  T obj;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_PETSC

#endif
