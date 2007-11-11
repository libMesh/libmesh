#ifndef __dof_object_test_h__
#define __dof_object_test_h__

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#define DOFOBJECTTEST \
  CPPUNIT_TEST( testValidId );

template <class DerivedClass>
class DofObjectTest {

private:
  DerivedClass * instance;
  
public:
  void setUp(DerivedClass * derived_instance)
  {
    instance=derived_instance;
  }

  void testValidId()
  {
    DofObject aobject(*instance);

    aobject.set_id(1);
    CPPUNIT_ASSERT( aobject.valid_id() );

    aobject.set_id(DofObject::invalid_id);
    CPPUNIT_ASSERT( !aobject.valid_id() );
  }
};

#endif // #ifdef __dof_object_test_h__
