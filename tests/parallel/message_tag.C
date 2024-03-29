
// libMesh includes
// Using *shims* here to test backwards compatibility via those
#include <libmesh/communicator.h>
#include <libmesh/message_tag.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class MessageTagTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MessageTagTest );

  CPPUNIT_TEST( testGetUniqueTagAuto );
  CPPUNIT_TEST( testGetUniqueTagManual );

  CPPUNIT_TEST_SUITE_END();

private:
  std::vector<std::string> _number;

public:
  void setUp()
  {}

  void tearDown()
  {}



  void testGetUniqueTagAuto()
  {
    LOG_UNIT_TEST;

    // We need to explicitly duplicate the communicator first, because
    // the original might already have tags used by other unit tests

    Parallel::Communicator newcomm;

    newcomm.duplicate(*TestCommWorld);

    const int n_vals = 5;
    const int n_vals_in_scope = 3;
    std::vector<int> vals(n_vals);

    {
      std::vector<Parallel::MessageTag> tags(n_vals_in_scope);
      for (int i=0; i != n_vals_in_scope; ++i)
        {
          tags[i] = newcomm.get_unique_tag();
          vals[i] = tags[i].value();
          for (int j=0; j != i; ++j)
            {
              CPPUNIT_ASSERT(vals[i] != vals[j]);
            }
        }
    }

    // Even after we go out of scope those values should be used up
    for (int i=n_vals_in_scope; i != n_vals; ++i)
      {
        Parallel::MessageTag another_tag = newcomm.get_unique_tag();
        vals[i] = another_tag.value();
        for (int j=0; j != i; ++j)
          {
            CPPUNIT_ASSERT(vals[i] != vals[j]);
          }
      }
  }



  void testGetUniqueTagManual()
  {
    LOG_UNIT_TEST;

    // Here we'll use the standard communicator, because even if it
    // used these tags in other contexts it should have freed them for
    // reuse later.

    const int requests[] = {2, 4, 6, 8, 8, 6, 8, 123, 3141, 3142};

    for (const int i : requests)
      {
        Parallel::MessageTag manual_tag =
          TestCommWorld->get_unique_tag(i);
        CPPUNIT_ASSERT_EQUAL(i, manual_tag.value());
      }
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( MessageTagTest );
