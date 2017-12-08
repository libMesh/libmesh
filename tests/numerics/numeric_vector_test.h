#ifndef __numeric_vector_test_h__
#define __numeric_vector_test_h__

// test includes
#include "test_comm.h"

// libMesh includes
#include <libmesh/parallel.h>
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#define NUMERICVECTORTEST                       \
  CPPUNIT_TEST( testLocalize );                 \
  CPPUNIT_TEST( testLocalizeBase );             \
  CPPUNIT_TEST( testLocalizeToOne );            \
  CPPUNIT_TEST( testLocalizeToOneBase );

#ifndef LIBMESH_HAVE_CXX14_MAKE_UNIQUE
using libMesh::make_unique;
#endif

template <class DerivedClass>
class NumericVectorTest : public CppUnit::TestCase {

protected:
  libMesh::Parallel::Communicator *my_comm;

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
      auto v_ptr = libmesh_make_unique<Derived>(*my_comm, global_size, local_size);
      Base & v = *v_ptr;
      std::vector<libMesh::Number> l(global_size);

      const libMesh::dof_id_type
        first = v.first_local_index(),
        last  = v.last_local_index();

      for (libMesh::dof_id_type n=first; n != last; n++)
        v.set (n, static_cast<libMesh::Number>(n));
      v.close();

      if (!to_one)
        v.localize(l);
      else
        v.localize_to_one(l,root_pid);

      if (!to_one || my_comm->rank() == root_pid)
        // Yes I really mean v.size()
        for (libMesh::dof_id_type i=0; i<v.size(); i++)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(libMesh::libmesh_real(i),
                                       libMesh::libmesh_real(l[i]),
                                       libMesh::TOLERANCE*libMesh::TOLERANCE);
    }
  }

  void testLocalize()
  {
    Localize<DerivedClass,DerivedClass>();
  }

  void testLocalizeBase()
  {
    Localize<libMesh::NumericVector<libMesh::Number>,DerivedClass>();
  }

  void testLocalizeToOne()
  {
    Localize<DerivedClass,DerivedClass >(true);
  }

  void testLocalizeToOneBase()
  {
    Localize<libMesh::NumericVector<libMesh::Number>,DerivedClass>(true);
  }
};

#endif // #ifdef __numeric_vector_test_h__
