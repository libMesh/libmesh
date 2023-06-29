#include "../geom/elem_test.h"

#include "libmesh/enum_to_string.h"
#include "libmesh/exodusII_io.h"

using namespace libMesh;

template <ElemType elem_type>
class ExodusTest : public PerElemTest<elem_type> {

public:

  void test_read_gold()
  {
  }

  void test_write()
  {
    LOG_UNIT_TEST;

    ExodusII_IO exii(*this->_mesh);
    exii.write("write_exodus_"+Utility::enum_to_string(elem_type)+".e");
  }
};

#define EXODUSTEST                              \
  CPPUNIT_TEST( test_read_gold );               \
  CPPUNIT_TEST( test_write );

#define INSTANTIATE_EXODUSTEST(elemtype)                        \
  class ExodusTest_##elemtype : public ExodusTest<elemtype> {   \
  public:                                                       \
  ExodusTest_##elemtype() :                                     \
    ExodusTest<elemtype>() {                                    \
    if (unitlog->summarized_logs_enabled())                     \
      this->libmesh_suite_name = "ExodusTest";                  \
    else                                                        \
      this->libmesh_suite_name = "ExodusTest_" #elemtype;       \
  }                                                             \
  CPPUNIT_TEST_SUITE( ExodusTest_##elemtype );                  \
  EXODUSTEST;                                                   \
  CPPUNIT_TEST_SUITE_END();                                     \
  };                                                            \
                                                                \
  CPPUNIT_TEST_SUITE_REGISTRATION( ExodusTest_##elemtype )

INSTANTIATE_EXODUSTEST(EDGE2);
INSTANTIATE_EXODUSTEST(EDGE3);
INSTANTIATE_EXODUSTEST(EDGE4);

#if LIBMESH_DIM > 1
INSTANTIATE_EXODUSTEST(TRI3);
INSTANTIATE_EXODUSTEST(TRI6);
INSTANTIATE_EXODUSTEST(TRI7);

INSTANTIATE_EXODUSTEST(QUAD4);
INSTANTIATE_EXODUSTEST(QUAD8);
INSTANTIATE_EXODUSTEST(QUAD9);
#endif // LIBMESH_DIM > 1

#if LIBMESH_DIM > 2
INSTANTIATE_EXODUSTEST(TET4);
INSTANTIATE_EXODUSTEST(TET10);
INSTANTIATE_EXODUSTEST(TET14);

INSTANTIATE_EXODUSTEST(HEX8);
INSTANTIATE_EXODUSTEST(HEX20);
INSTANTIATE_EXODUSTEST(HEX27);

INSTANTIATE_EXODUSTEST(PRISM6);
INSTANTIATE_EXODUSTEST(PRISM15);
INSTANTIATE_EXODUSTEST(PRISM18);
INSTANTIATE_EXODUSTEST(PRISM20);
INSTANTIATE_EXODUSTEST(PRISM21);

INSTANTIATE_EXODUSTEST(PYRAMID5);
INSTANTIATE_EXODUSTEST(PYRAMID13);
INSTANTIATE_EXODUSTEST(PYRAMID14);
INSTANTIATE_EXODUSTEST(PYRAMID18);
#endif // LIBMESH_DIM > 2
