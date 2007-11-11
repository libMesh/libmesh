#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <point.h>

#include "../numerics/type_vector_test.h"

#define POINTTEST TYPEVECTORTEST

template <class DerivedClass>
class PointTestBase : public TypeVectorTestBase<DerivedClass> { 
};
