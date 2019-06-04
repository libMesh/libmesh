#include <libmesh/point.h>

#include "../numerics/type_vector_test.h"
#include "libmesh_cppunit.h"

#define POINTTEST TYPEVECTORTEST

using namespace libMesh;

template <class DerivedClass>
class PointTestBase : public TypeVectorTestBase<DerivedClass> {
public:
  virtual void setUp()
  {
    TypeVectorTestBase<DerivedClass>::setUp();
  }

  virtual void tearDown()
  {
    TypeVectorTestBase<DerivedClass>::tearDown();
  }
};
