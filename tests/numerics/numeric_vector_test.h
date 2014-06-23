#ifndef __numeric_vector_test_h__
#define __numeric_vector_test_h__

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

using namespace libMesh;

template <class DerivedClass>
class NumericVectorTest : public CppUnit::TestCase {

public:
  void setUp()
  {}

  void tearDown()
  {}

  template <class Base, class Derived>
  void Localize(bool to_one=false)
  {
    libMesh::Parallel::Communicator CommTest(libMesh::GLOBAL_COMM_WORLD);

    const processor_id_type root_pid = 0;
    unsigned int block_size  = 10;

    // a different size on each processor.
    unsigned int local_size  = block_size +
      static_cast<unsigned int>(libMesh::global_processor_id());
    unsigned int global_size = 0;

    for (processor_id_type p=0; p<libMesh::global_n_processors(); p++)
      global_size += (block_size + static_cast<unsigned int>(p));

    {
      Base & v = *(new Derived(CommTest, global_size, local_size));
      std::vector<Number> l(global_size);

      const dof_id_type
        first = v.first_local_index(),
        last  = v.last_local_index();

      for (dof_id_type n=first; n != last; n++)
        v.set (n, static_cast<Number>(n));
      v.close();

      if(!to_one)
        v.localize(l);
      else
        v.localize_to_one(l,root_pid);

      if(!to_one || libMesh::global_processor_id() == root_pid)
        //Yes I really mean v.size()
        for (dof_id_type i=0; i<v.size(); i++)
          CPPUNIT_ASSERT_DOUBLES_EQUAL( libmesh_real(i) , libmesh_real(l[i]) , TOLERANCE*TOLERANCE );

      delete &v;
    }
  }

  void testLocalize()
  {
    Localize<DerivedClass,DerivedClass>();
  }

  void testLocalizeBase()
  {
    Localize<NumericVector<Number>,DerivedClass>();
  }

  void testLocalizeToOne()
  {
    Localize<DerivedClass,DerivedClass >(true);
  }

  void testLocalizeToOneBase()
  {
    Localize<NumericVector<Number>,DerivedClass>(true);
  }
};

#endif // #ifdef __numeric_vector_test_h__
