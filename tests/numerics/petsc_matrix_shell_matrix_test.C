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

#include <libmesh/petsc_matrix_shell_matrix.h>
#include <libmesh/petsc_vector.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

using namespace libMesh;

namespace
{
// A trivial MATOP_MULT: y = 2*x. Just enough to prove a PetscMatrixShellMatrix
// is a genuinely usable MatMult operand, e.g. as a matrix-free SNES Jacobian.
PetscErrorCode
petsc_matrix_shell_matrix_test_mult(Mat, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  LibmeshPetscCallQ(VecCopy(x, y));
  LibmeshPetscCallQ(VecScale(y, 2.));
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}
}

class PetscMatrixShellMatrixTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(PetscMatrixShellMatrixTest);

  CPPUNIT_TEST(testInitSizes);
  CPPUNIT_TEST(testGetContext);
  CPPUNIT_TEST(testMult);

  CPPUNIT_TEST_SUITE_END();

public:

  void setUp() {}
  void tearDown() {}

  void testInitSizes()
  {
    LOG_UNIT_TEST;

    PetscMatrixShellMatrix<Number> mat(*TestCommWorld);
    const numeric_index_type m = TestCommWorld->size();
    mat.init(m, m, /*m_l=*/1, /*n_l=*/1);

    CPPUNIT_ASSERT_EQUAL(mat.m(), m);
    CPPUNIT_ASSERT_EQUAL(mat.n(), m);
    CPPUNIT_ASSERT_EQUAL(mat.local_m(), numeric_index_type(1));
    CPPUNIT_ASSERT_EQUAL(mat.local_n(), numeric_index_type(1));

    MatType type;
    LibmeshPetscCall2((*TestCommWorld), MatGetType(mat.mat(), &type));
    CPPUNIT_ASSERT(std::string(type) == MATSHELL);
  }

  void testGetContext()
  {
    LOG_UNIT_TEST;

    PetscMatrixShellMatrix<Number> mat(*TestCommWorld);
    const numeric_index_type m = TestCommWorld->size();
    mat.init(m, m, 1, 1);

    // init() composes libmesh's own context onto the shell Mat, letting
    // get_context() recover the wrapping object later given only the raw Mat.
    CPPUNIT_ASSERT(PetscMatrixBase<Number>::get_context(mat.mat(), *TestCommWorld) == &mat);
  }

  void testMult()
  {
    LOG_UNIT_TEST;

    PetscMatrixShellMatrix<Number> mat(*TestCommWorld);
    const numeric_index_type m = TestCommWorld->size();
    mat.init(m, m, 1, 1);

    LibmeshPetscCall2((*TestCommWorld),
                      MatShellSetOperation(mat.mat(), MATOP_MULT,
                                           (void (*)(void))petsc_matrix_shell_matrix_test_mult));

    PetscVector<Number> x(*TestCommWorld, m, 1);
    PetscVector<Number> y(*TestCommWorld, m, 1);
    x.set(TestCommWorld->rank(), 3.);
    x.close();

    LibmeshPetscCall2((*TestCommWorld), MatMult(mat.mat(), x.vec(), y.vec()));

    LIBMESH_ASSERT_NUMBERS_EQUAL(Number(6.), y(TestCommWorld->rank()), TOLERANCE);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(PetscMatrixShellMatrixTest);

#endif // LIBMESH_HAVE_PETSC
