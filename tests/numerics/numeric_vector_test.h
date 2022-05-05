#ifndef __numeric_vector_test_h__
#define __numeric_vector_test_h__

// test includes
#include "test_comm.h"

// libMesh includes
#include <libmesh/parallel.h>

#include "libmesh_cppunit.h"

#include <memory>

#define NUMERICVECTORTEST                       \
  CPPUNIT_TEST( testLocalize );                 \
  CPPUNIT_TEST( testLocalizeBase );             \
  CPPUNIT_TEST( testLocalizeIndices );          \
  CPPUNIT_TEST( testLocalizeIndicesBase );      \
  CPPUNIT_TEST( testLocalizeToOne );            \
  CPPUNIT_TEST( testLocalizeToOneBase );


template <class DerivedClass>
class NumericVectorTest : public CppUnit::TestCase {

protected:
  libMesh::Parallel::Communicator *my_comm;

  std::string libmesh_suite_name;

public:
  void setUp()
  {
    // By default we'll use the whole communicator in parallel;
    // Serial-only NumericVector subclasses will need to override
    // this.
    my_comm = TestCommWorld;
  }

  void tearDown()
  {}

  template <class Base, class Derived>
  void Localize(bool to_one=false)
  {
    const libMesh::processor_id_type root_pid = 0;
    unsigned int block_size  = 10;

    // a different size on each processor.
    unsigned int local_size  = block_size +
      static_cast<unsigned int>(my_comm->rank());
    unsigned int global_size = 0;

    for (libMesh::processor_id_type p=0; p<my_comm->size(); p++)
      global_size += (block_size + static_cast<unsigned int>(p));

    {
      auto v_ptr = std::make_unique<Derived>(*my_comm, global_size, local_size);
      Base & v = *v_ptr;
      std::vector<libMesh::Number> l(global_size);

      const libMesh::dof_id_type
        first = v.first_local_index(),
        last  = v.last_local_index();

      for (libMesh::dof_id_type n=first; n != last; n++)
        v.set (n, static_cast<libMesh::Number>(n));
      v.close();
      for (libMesh::dof_id_type n=first; n != last; n++)
        v.add (n, static_cast<libMesh::Number>(n));
      v.close();

      if (!to_one)
        v.localize(l);
      else
        v.localize_to_one(l,root_pid);

      if (!to_one || my_comm->rank() == root_pid)
        // Yes I really mean v.size()
        for (libMesh::dof_id_type i=0; i<v.size(); i++)
          LIBMESH_ASSERT_FP_EQUAL(2.*libMesh::libmesh_real(i),
                                  libMesh::libmesh_real(l[i]),
                                  libMesh::TOLERANCE*libMesh::TOLERANCE);

      for (libMesh::dof_id_type n=first; n != last; n++)
      {
        const auto value = static_cast<libMesh::Number>(n);
        v.insert (&value, std::vector<libMesh::numeric_index_type>({n}));
      }
      v.close();

      if (!to_one)
        v.localize(l);
      else
        v.localize_to_one(l,root_pid);

      if (!to_one || my_comm->rank() == root_pid)
        // Yes I really mean v.size()
        for (libMesh::dof_id_type i=0; i<v.size(); i++)
          LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(i),
                                  libMesh::libmesh_real(l[i]),
                                  libMesh::TOLERANCE*libMesh::TOLERANCE);
    }
  }


  template <class Base, class Derived>
  void LocalizeIndices()
  {
    unsigned int block_size  = 10;

    // a different size on each processor.
    unsigned int local_size  = block_size +
      static_cast<unsigned int>(my_comm->rank());
    unsigned int global_size = 0;

    for (libMesh::processor_id_type p=0; p<my_comm->size(); p++)
      global_size += (block_size + static_cast<unsigned int>(p));

    {
      auto v_ptr = std::make_unique<Derived>(*my_comm, global_size, local_size);
      Base & v = *v_ptr;

      // Let's try pulling the same number of entries from each processor
      std::vector<libMesh::Number> values(block_size * my_comm->size());
      std::vector<libMesh::dof_id_type> indices;
      indices.reserve(block_size * my_comm->size());

      const libMesh::dof_id_type
        first = v.first_local_index(),
        last  = v.last_local_index();

      for (libMesh::dof_id_type n=first; n != last; n++)
        v.set (n, static_cast<libMesh::Number>(n));
      v.close();

      libMesh::dof_id_type end_index = 0;
      for (libMesh::processor_id_type p=0; p<my_comm->size(); p++)
        {
          end_index += block_size + p;
          for (unsigned int j = 0; j != block_size; ++j)
            indices.push_back(end_index-j-1);
        }

      v.localize(values, indices);

      end_index = 0;
      for (libMesh::processor_id_type p=0; p<my_comm->size(); p++)
        {
          end_index += block_size + p;
          for (unsigned int j = 0; j != block_size; ++j)
            LIBMESH_ASSERT_FP_EQUAL(libMesh::libmesh_real(values[p*block_size+j]),
                                    libMesh::libmesh_real(end_index-j-1),
                                    libMesh::TOLERANCE*libMesh::TOLERANCE);
        }
    }
  }

  void testLocalize()
  {
    LOG_UNIT_TEST;

    Localize<DerivedClass,DerivedClass>();
  }

  void testLocalizeBase()
  {
    LOG_UNIT_TEST;

    Localize<libMesh::NumericVector<libMesh::Number>,DerivedClass>();
  }

  void testLocalizeToOne()
  {
    LOG_UNIT_TEST;

    Localize<DerivedClass,DerivedClass >(true);
  }

  void testLocalizeToOneBase()
  {
    LOG_UNIT_TEST;

    Localize<libMesh::NumericVector<libMesh::Number>,DerivedClass>(true);
  }

  void testLocalizeIndices()
  {
    LOG_UNIT_TEST;

    LocalizeIndices<DerivedClass,DerivedClass >();
  }

  void testLocalizeIndicesBase()
  {
    LOG_UNIT_TEST;

    LocalizeIndices<libMesh::NumericVector<libMesh::Number>,DerivedClass>();
  }
};

#endif // #ifdef __numeric_vector_test_h__
