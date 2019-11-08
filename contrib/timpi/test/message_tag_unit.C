#include <timpi/communicator.h>
#include <timpi/message_tag.h>
#include <timpi/timpi_init.h>

#define TIMPI_UNIT_ASSERT(expr) \
  if (!(expr)) \
    timpi_error();

using namespace TIMPI;

Communicator *TestCommWorld;

  void testGetUniqueTagAuto()
  {
    // We need to duplicate the communicator first, because the
    // original might already have tags used by other unit tests

    TIMPI::Communicator newcomm;

    TestCommWorld->duplicate(newcomm);

    const int n_vals = 5;
    const int n_vals_in_scope = 3;
    std::vector<int> vals(n_vals);

    {
      std::vector<TIMPI::MessageTag> tags(n_vals_in_scope);
      for (int i=0; i != n_vals_in_scope; ++i)
        {
          tags[i] = newcomm.get_unique_tag();
          vals[i] = tags[i].value();
          for (int j=0; j != i; ++j)
            {
              TIMPI_UNIT_ASSERT(vals[i] != vals[j]);
            }
        }
    }

    // Even after we go out of scope those values should be used up
    for (int i=n_vals_in_scope; i != n_vals; ++i)
      {
        TIMPI::MessageTag another_tag = newcomm.get_unique_tag();
        vals[i] = another_tag.value();
        for (int j=0; j != i; ++j)
          {
            TIMPI_UNIT_ASSERT(vals[i] != vals[j]);
          }
      }
  }



  void testGetUniqueTagManual()
  {
    // Here we'll use the standard communicator, because even if it
    // used these tags in other contexts it should have freed them for
    // reuse later.

    const int requests[] = {2, 4, 6, 8, 8, 6, 8, 123, 3141, 3142};

    for (const int i : requests)
      {
        TIMPI::MessageTag manual_tag =
          TestCommWorld->get_unique_tag(i);
        TIMPI_UNIT_ASSERT(i == manual_tag.value());
      }
  }

int main(int argc, const char * const * argv)
{
  TIMPI::TIMPIInit init(argc, argv);
  TestCommWorld = &init.comm();

  testGetUniqueTagAuto();
  testGetUniqueTagManual();

  return 0;
}
