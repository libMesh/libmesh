#include "test_comm.h"

// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/parallel.h>
#include <libmesh/null_output_iterator.h>

#include <vector>

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

namespace libMesh {
namespace Parallel {
template <typename T>
class Packing<std::basic_string<T> > {
public:

  static const unsigned int size_bytes = 4;

  typedef T buffer_type;

  static unsigned int
  get_string_len (typename std::vector<T>::const_iterator in)
  {
    unsigned int string_len = reinterpret_cast<const unsigned char &>(in[size_bytes-1]);
    for (signed int i=size_bytes-2; i >= 0; --i)
      {
        string_len *= 256;
        string_len += reinterpret_cast<const unsigned char &>(in[i]);
      }
    return string_len;
  }


  static unsigned int
  packed_size (typename std::vector<T>::const_iterator in)
  {
    return get_string_len(in) + size_bytes;
  }

  static unsigned int packable_size
  (const std::basic_string<T> & s,
   const void *)
  {
    return s.size() + size_bytes;
  }


  template <typename Iter>
  static void pack (const std::basic_string<T> & b, Iter data_out,
                    const void *)
  {
    unsigned int string_len = b.size();
    for (unsigned int i=0; i != size_bytes; ++i)
      {
        *data_out++ = (string_len % 256);
        string_len /= 256;
      }
    std::copy(b.begin(), b.end(), data_out);
  }

  static std::basic_string<T>
  unpack (typename std::vector<T>::const_iterator in, void *)
  {
    unsigned int string_len = get_string_len(in);

    std::ostringstream oss;
    for (unsigned int i = 0; i < string_len; ++i)
      oss << reinterpret_cast<const unsigned char &>(in[i+size_bytes]);

    in += size_bytes + string_len;

    //  std::cout << oss.str() << std::endl;
    return std::string(oss.str());
  }

};


} // namespace Parallel

} // namespace libMesh



using namespace libMesh;

class PackedRangeTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( PackedRangeTest );

  CPPUNIT_TEST( testNullAllGather );
  CPPUNIT_TEST( testNullSendReceive );
  CPPUNIT_TEST( testContainerSendReceive );
  //  CPPUNIT_TEST( testAdapterSendReceive );
  //  CPPUNIT_TEST( testPointerAdapterSendReceive );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}



  void testNullAllGather()
  {
    std::vector<processor_id_type> vals;

    std::vector<std::string> send(1);
    if (TestCommWorld->rank() == 0)
      send[0].assign("Hello");
    else
      send[0].assign("Goodbye");

    TestCommWorld->allgather_packed_range
      ((void *)(NULL), send.begin(), send.end(),
       libMesh::null_output_iterator<std::string>());
  }


  void testNullSendReceive()
  {
    std::vector<processor_id_type> vals;

    std::vector<std::string> send(1);
    const unsigned int my_rank = TestCommWorld->rank();
    const unsigned int dest_rank =
      (my_rank + 1) % TestCommWorld->size();
    const unsigned int source_rank =
      (my_rank + TestCommWorld->size() - 1) % TestCommWorld->size();

    {
      std::ostringstream os;
      os << my_rank;
      send[0] = os.str();
    }

    TestCommWorld->send_receive_packed_range
      (dest_rank, (void *)(NULL), send.begin(), send.end(),
       source_rank, (void *)(NULL),
       libMesh::null_output_iterator<std::string>(),
       (std::string*)NULL);
  }


  void testContainerSendReceive()
  {
    std::vector<processor_id_type> vals;

    std::vector<std::string> send(1), recv;

    const unsigned int my_rank = TestCommWorld->rank();
    const unsigned int dest_rank =
      (my_rank + 1) % TestCommWorld->size();
    const unsigned int source_rank =
      (my_rank + TestCommWorld->size() - 1) % TestCommWorld->size();

    {
      std::ostringstream os;
      os << my_rank;
      send[0] = os.str();
    }

    TestCommWorld->send_receive_packed_range
      (dest_rank, (void *)(NULL), send.begin(), send.end(),
       source_rank, (void *)(NULL),
       std::back_inserter(recv),
       (std::string*)NULL);

    CPPUNIT_ASSERT_EQUAL(recv.size(), std::size_t(1));

    std::string check;
    {
      std::ostringstream os;
      os << source_rank;
      check = os.str();
    }

    CPPUNIT_ASSERT_EQUAL(recv[0], check);
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( PackedRangeTest );
