#include <libmesh/petsc_vector.h>

#ifdef LIBMESH_HAVE_PETSC

#include "numeric_vector_test.h"


using namespace libMesh;

class PetscVectorTest : public NumericVectorTest<PetscVector<Number>> {
public:
  PetscVectorTest() :
    NumericVectorTest<PetscVector<Number>>() {
    if (unitlog->summarized_logs_enabled())
      this->libmesh_suite_name = "NumericVectorTest";
    else
      this->libmesh_suite_name = "PetscVectorTest";
  }

  CPPUNIT_TEST_SUITE( PetscVectorTest );

  NUMERICVECTORTEST

  CPPUNIT_TEST( testGetArray );

  CPPUNIT_TEST( testPetscOperations );

  CPPUNIT_TEST_SUITE_END();

  void testGetArray()
  {
    LOG_UNIT_TEST;

    unsigned int min_block_size  = 2;

    // a different size on each processor.
    unsigned int my_p = my_comm->rank();
    unsigned int local_size  = (min_block_size + my_p);
    unsigned int global_size = 0;
    unsigned int my_offset = 0;

    for (libMesh::processor_id_type p=0; p<my_comm->size(); p++)
      {
        const unsigned int p_size =
          (min_block_size + static_cast<unsigned int>(p));
        global_size += p_size;
        if (p < my_p)
          my_offset += p_size;
      }

    PetscVector<Number> v(*my_comm, global_size, local_size);

    PetscScalar * values = v.get_array();

    for (unsigned int i=0; i<local_size; i++)
      values[i] = i;

    v.restore_array();

    v.close();

    // Check the values through the interface
    for (unsigned int i=0; i<local_size; i++)
      LIBMESH_ASSERT_FP_EQUAL(i, std::abs(v(my_offset + i)), TOLERANCE*TOLERANCE);

    // Check that we can see the same thing with get_array_read
    const PetscScalar * read_only_values = v.get_array_read();

    for (unsigned int i=0; i<local_size; i++)
      LIBMESH_ASSERT_FP_EQUAL(i, std::abs(read_only_values[i]), TOLERANCE*TOLERANCE);

    v.restore_array();

    // Test getting a read only array after getting a writable array
    values = v.get_array();
    read_only_values = v.get_array_read();
    CPPUNIT_ASSERT_EQUAL(read_only_values, const_cast<const PetscScalar *>(values));

    v.restore_array();

    // Test to make sure we can get arrays after other operators
    for (unsigned int i = 0; i < local_size; i++)
      v.set(my_offset + i, i * 2.0);
    v.close();
    for (unsigned int i = 0; i < local_size; i++)
      LIBMESH_ASSERT_FP_EQUAL(i * 2.0, std::abs(v(my_offset + i)), TOLERANCE * TOLERANCE);
    values = v.get_array();
    read_only_values = v.get_array_read();
    for (unsigned int i = 0; i < local_size; i++)
    {
      LIBMESH_ASSERT_FP_EQUAL(i * 2.0, std::abs(values[i]), TOLERANCE * TOLERANCE);
      LIBMESH_ASSERT_FP_EQUAL(i * 2.0, std::abs(read_only_values[i]), TOLERANCE * TOLERANCE);
    }
    v.restore_array();
  }

  void testPetscOperations()
  {
    PetscVector<Number> v1(*my_comm, global_size, local_size);
    auto v2 = v1.clone();

    const libMesh::dof_id_type
      first = v1.first_local_index(),
      last  = v1.last_local_index();

    for (libMesh::dof_id_type n=first; n != last; n++)
    {
      v1.set (n, static_cast<libMesh::Number>(n+1));
      v2->set (n, static_cast<libMesh::Number>(2*(n+1)));
    }
    v1.close();
    v2->close();

    auto v_working_ptr = v1.clone();
    auto & v_working = *v_working_ptr;

    v_working.pointwise_mult(v1, *v2);

    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v_working(n)),
                              libMesh::Real((n+1)*2*(n+1)),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);

    v_working.pointwise_divide(v1, *v2);

    for (libMesh::dof_id_type n=first; n != last; n++)
      LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(v_working(n)),
                              libMesh::Real(0.5),
                              libMesh::TOLERANCE*libMesh::TOLERANCE);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( PetscVectorTest );

#endif // #ifdef LIBMESH_HAVE_PETSC
