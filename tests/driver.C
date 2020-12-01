// Ignore overloaded-virtual warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <libmesh/restore_warnings.h>

// libMesh includes
#include <libmesh/libmesh.h>
#include "test_comm.h"

#ifdef LIBMESH_HAVE_CXX11_REGEX

// C++ includes
#include <regex>

// Add Tests to runner that match user-provided regex.
void add_matching_tests_to_runner(CppUnit::Test * test,
                                  const std::string & r_str,
                                  const std::regex & r,
                                  CppUnit::TextUi::TestRunner & runner,
                                  CppUnit::TestSuite & rejects)
{
  // If we running all tests we just add the "All Tests" test and then return
  if (test->getName() == "All Tests" && r_str == "All Tests")
  {
    libMesh::out << test->getName() << std::endl;
    runner.addTest(test);
    return;
  }

  if (test->getChildTestCount() == 0)
  {
    // Add the test to the runner
    if (std::regex_search(test->getName(), r))
    {
      libMesh::out << test->getName() << std::endl;
      runner.addTest(test);
    }
    // Add the test to the rejects it can be cleaned up later
    else
      rejects.addTest(test);
  }

  // Call this recursively on each of our children, if any.
  for (int i = 0; i < test->getChildTestCount(); i++)
    add_matching_tests_to_runner(test->getChildTestAt(i), r_str, r, runner, rejects);
}

#endif


int main(int argc, char ** argv)
{
  // Initialize the library.  This is necessary because the library
  // may depend on a number of other libraries (i.e. MPI  and Petsc)
  // that require initialization before use.
  libMesh::LibMeshInit init(argc, argv);
  TestCommWorld = &init.comm();

  // We can now run all tests that match a regular expression, for
  // example, "--re PartitionerTest" will match all the Partitioner
  // unit tests.  If the user does not specify a regex, we run all the
  // tests returned by makeTest().

  // An example regex_string that would _exactly_ match a _single_ test is:
  // "PartitionerTest_HilbertSFCPartitioner_ReplicatedMesh::testPartition2"
  // On the other hand, the regex "HilbertSFC" would match all of the
  // following tests:
  //
  // PartitionerTest_HilbertSFCPartitioner_ReplicatedMesh
  // PartitionerTest_HilbertSFCPartitioner_ReplicatedMesh::testPartitionEmpty
  // PartitionerTest_HilbertSFCPartitioner_ReplicatedMesh::testPartition1
  // PartitionerTest_HilbertSFCPartitioner_ReplicatedMesh::testPartition2
  // PartitionerTest_HilbertSFCPartitioner_ReplicatedMesh::testPartitionNProc
  //
  // If the user does not provide a a regex, the default re is "All Tests",
  // which runs all the unit tests.

  // Read command line argument specifying the regular expression.
  std::string regex_string = "All Tests";
  regex_string = libMesh::command_line_next("--re", regex_string);

  // Recursively add tests matching the regex tothe runner object.
  CppUnit::TextUi::TestRunner runner;

  // The Cppunit registry object that knows about all the tests.
  CppUnit::TestFactoryRegistry & registry = CppUnit::TestFactoryRegistry::getRegistry();

  // A test suite container for holding tests not matching the regex. When main's
  // scope ends, this class's destructor will delete the rejected tests
  CppUnit::TestSuite rejects("rejects");

#ifdef LIBMESH_HAVE_CXX11_REGEX
  // Make regex object from user's input.
  std::regex the_regex(regex_string);

  // Add all tests which match the re to the runner object.
  libMesh::out << "Will run the following tests:" << std::endl;
  add_matching_tests_to_runner(registry.makeTest(), regex_string, the_regex, runner, rejects);
#else
  // If no C++11 <regex> just run all the tests.
  runner.addTest(registry.makeTest());
#endif

  // Finally, run the matching tests.
  if (runner.run())
    return 0;

  return 1;
}

libMesh::Parallel::Communicator *TestCommWorld;
