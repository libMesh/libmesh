#ifndef __dof_object_test_h__
#define __dof_object_test_h__

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#define DOFOBJECTTEST \
  CPPUNIT_TEST( testSetId ); \
  CPPUNIT_TEST( testValidId ); \
  CPPUNIT_TEST( testInvalidateId ); \
  CPPUNIT_TEST( testSetProcId ); \
  CPPUNIT_TEST( testValidProcId ); \
  CPPUNIT_TEST( testInvalidateProcId ); \
  CPPUNIT_TEST( testSetNSystems ); \
  CPPUNIT_TEST( testSetNVariableGroups );

using namespace libMesh;

template <class DerivedClass>
class DofObjectTest {

private:
  DerivedClass * instance;
  
public:
  void setUp(DerivedClass * derived_instance)
  {
    instance=derived_instance;
  }

  void testSetId()
  {
    DofObject aobject(*instance);

    aobject.set_id(1);
    CPPUNIT_ASSERT_EQUAL( (unsigned int)1 , aobject.id() );
  }
  
  void testValidId()
  {
    DofObject aobject(*instance);

    aobject.set_id(1);
    CPPUNIT_ASSERT( aobject.valid_id() );

    aobject.set_id(DofObject::invalid_id);
    CPPUNIT_ASSERT( !aobject.valid_id() );
  }

  void testInvalidateId()
  {
    DofObject aobject(*instance);

    aobject.set_id(1);
    aobject.invalidate_id();

    CPPUNIT_ASSERT( !aobject.valid_id() );
  }

  void testSetProcId()
  {
    DofObject aobject(*instance);

    aobject.processor_id(libMesh::processor_id());
    CPPUNIT_ASSERT_EQUAL( (short unsigned int)libMesh::processor_id() , aobject.processor_id() );
  }

  void testValidProcId()
  {
    DofObject aobject(*instance);

    aobject.processor_id(libMesh::processor_id());
    CPPUNIT_ASSERT(aobject.valid_processor_id());
    
    aobject.processor_id(DofObject::invalid_processor_id);
    CPPUNIT_ASSERT(!aobject.valid_processor_id());
  }

  void testInvalidateProcId()
  {
    DofObject aobject(*instance);

    aobject.processor_id(libMesh::processor_id());
    aobject.invalidate_processor_id();

    CPPUNIT_ASSERT( !aobject.valid_processor_id() );
  }

  void testSetNSystems()
  {
    DofObject aobject(*instance);    
    
    aobject.set_n_systems (10);
    
    CPPUNIT_ASSERT_EQUAL( (unsigned int) 10, aobject.n_systems() );
  }

  void testSetNVariableGroups()
  {
    DofObject aobject(*instance);    
    
    aobject.set_n_systems (2);
    
    std::vector<unsigned int> nvpg;
    
    nvpg.push_back(10);
    nvpg.push_back(20);
    nvpg.push_back(30);

    aobject.set_n_vars_per_group (0, nvpg);
    aobject.set_n_vars_per_group (1, nvpg);
    
    for (unsigned int s=0; s<2; s++)
      {
	CPPUNIT_ASSERT_EQUAL( (unsigned int) 60, aobject.n_vars(s) );
	CPPUNIT_ASSERT_EQUAL( (unsigned int) 3,  aobject.n_var_groups(s) );

	for (unsigned int vg=0; vg<3; vg++)
	  CPPUNIT_ASSERT_EQUAL( nvpg[vg], aobject.n_vars(s,vg) );
      }
  }
};

#endif // #ifdef __dof_object_test_h__
