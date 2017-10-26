#include <libmesh/petsc_vector.h>

#ifdef LIBMESH_HAVE_PETSC

#include "numeric_vector_test.h"

#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

using namespace libMesh;

class PetscVectorTest : public NumericVectorTest<PetscVector<Number>> {
public:
  CPPUNIT_TEST_SUITE( PetscVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST( testGetArray );

  CPPUNIT_TEST_SUITE_END();

  void testGetArray()
  {
    unsigned int block_size  = 2;

    // a different size on each processor.
    unsigned int local_size  = block_size;
    unsigned int global_size = 0;

    for (libMesh::processor_id_type p=0; p<my_comm->size(); p++)
      global_size += (block_size + static_cast<unsigned int>(p));

    PetscVector<Number> v(*my_comm, global_size, local_size);

    PetscScalar * values = v.get_array();

    for (unsigned int i=0; i<local_size; i++)
      values[i] = i;

    v.restore_array();

    v.close();

    // Check the values through the interface
    for (unsigned int i=0; i<local_size; i++)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(std::abs(v(my_comm->rank()*2 + i)), i, TOLERANCE*TOLERANCE);

    // Check that we can see the same thing with get_array_read
    const PetscScalar * read_only_values = v.get_array_read();

    for (unsigned int i=0; i<local_size; i++)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(std::abs(read_only_values[i]), i, TOLERANCE*TOLERANCE);

    v.restore_array();

    // Test getting a read only array after getting a writable array
    values = v.get_array();
    read_only_values = v.get_array_read();
    CPPUNIT_ASSERT_EQUAL((intptr_t)read_only_values, (intptr_t)values);

    v.restore_array();
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( PetscVectorTest );

#endif // #ifdef LIBMESH_HAVE_PETSC
