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

#include <libmesh/libmesh_config.h>

#ifdef LIBMESH_HAVE_PETSC

#include <libmesh/petsc_mffd_matrix.h>
#include <libmesh/petsc_matrix_shell_matrix.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

class PetscMFFDMatrixTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(PetscMFFDMatrixTest);

  CPPUNIT_TEST(testAssignWithContext);
  CPPUNIT_TEST(testAssignWithoutContext);

  CPPUNIT_TEST_SUITE_END();

public:

  void setUp() {}
  void tearDown() {}

  void testAssignWithContext()
  {
    LOG_UNIT_TEST;

    // A real, owned Mat for the wrapper to adopt. init() attaches its own
    // context (pointing at owner) as a side effect.
    PetscMatrixShellMatrix<Number> owner(*TestCommWorld);
    owner.init(TestCommWorld->size(), TestCommWorld->size(), 1, 1);

    PetscMFFDMatrix<Number> mffd(*TestCommWorld);
    mffd.assign(owner.mat(), /*set_context=*/true);

    CPPUNIT_ASSERT(mffd.mat() == owner.mat());
    CPPUNIT_ASSERT(PetscMatrixBase<Number>::get_context(owner.mat(), *TestCommWorld) == &mffd);
  }

  void testAssignWithoutContext()
  {
    LOG_UNIT_TEST;

    PetscMatrixShellMatrix<Number> owner(*TestCommWorld);
    owner.init(TestCommWorld->size(), TestCommWorld->size(), 1, 1);

    {
      PetscMFFDMatrix<Number> mffd(*TestCommWorld);
      mffd.assign(owner.mat(), /*set_context=*/false);
      CPPUNIT_ASSERT(mffd.mat() == owner.mat());
    }

    // mffd is now destroyed. Since assign() above didn't attach a context,
    // owner's own context (set by its own init()) must still be intact --
    // this is the dangling-context scenario set_context=false exists to avoid.
    CPPUNIT_ASSERT(PetscMatrixBase<Number>::get_context(owner.mat(), *TestCommWorld) == &owner);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(PetscMFFDMatrixTest);

#endif // LIBMESH_HAVE_PETSC
