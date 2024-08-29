// Ignore overloaded-virtual warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/TestResult.h>
#include <libmesh/restore_warnings.h>

// libMesh includes
#include <libmesh/libmesh.h>

#include "libmesh_cppunit.h"
#include "test_comm.h"

#ifdef LIBMESH_HAVE_CXX11_REGEX

// C++ includes
#include <regex>

// Add Tests to runner that match user-provided regex.
int add_matching_tests_to_runner(CppUnit::Test * test,
                                 const std::string & allow_r_str,
                                 const std::regex & allow_r,
                                 const std::string & deny_r_str,
                                 const std::regex & deny_r,
                                 CppUnit::TextUi::TestRunner & runner,
                                 CppUnit::TestSuite & rejects)
{
  int n_tests_added = 0;

  // If running all tests, just add the "All Tests" test and return
  if (test->getName() == "All Tests" && allow_r_str == "All Tests" &&
      deny_r_str == "^$")
  {
    libMesh::out << test->getName() << std::endl;
    runner.addTest(test);
    return -12345;
  }

  if (test->getChildTestCount() == 0)
  {
    // Add the test to the runner
    if ((allow_r_str == "All Tests" ||
         std::regex_search(test->getName(), allow_r)) &&
        !std::regex_search(test->getName(), deny_r))
    {
      libMesh::out << test->getName() << std::endl;
      n_tests_added ++;
      runner.addTest(test);
    }
    // Add the test to the rejects it can be cleaned up later
    else
      rejects.addTest(test);
  }

  // Call this recursively on each of our children, if any.
  for (int i = 0; i < test->getChildTestCount(); i++)
    n_tests_added +=
      add_matching_tests_to_runner(test->getChildTestAt(i), allow_r_str, allow_r,
                                   deny_r_str, deny_r, runner, rejects);

  return n_tests_added;
}

#endif


int main(int argc, char ** argv)
{
  // Initialize the library.  This is necessary because the library
  // may depend on a number of other libraries (i.e. MPI  and Petsc)
  // that require initialization before use.
  libMesh::LibMeshInit init(argc, argv);
  TestCommWorld = &init.comm();

  // See how long each of our tests are taking to run.  This should be
  // coarse-grained enough to enable even if we're not performance
  // logging inside the library itself.  We need to do this early, so
  // we can query unitlog when first initializing tests.
  libMesh::PerfLog driver_unitlog ("Unit Tests");
  unitlog = &driver_unitlog;

  // Print just logs summarized by test suite, not every test
  // individually
  if (!libMesh::on_command_line("--full-logs"))
    driver_unitlog.enable_summarized_logs();

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

  // We can also skip tests that match a regular expression, with e.g.
  // "--deny_re PartitionerTest" to skip all the Partitioner unit
  // tests (even if a "--re" option would have included them.)

  // Read command line argument specifying the allowlist regular expression.
  std::string allow_regex_string = "All Tests";
  allow_regex_string = libMesh::command_line_next("--re", allow_regex_string);

  // Read command line argument specifying the allowlist regular expression.
  std::string deny_regex_string = "^$";
  deny_regex_string = libMesh::command_line_next("--deny_re", deny_regex_string);

  // Recursively add tests matching the regex to the runner object.
  CppUnit::TextUi::TestRunner runner;

  // The Cppunit registry object that knows about all the tests.
  CppUnit::TestFactoryRegistry & registry = CppUnit::TestFactoryRegistry::getRegistry();

  // A test suite container for holding tests not matching the regex. When main's
  // scope ends, this class's destructor will delete the rejected tests
  CppUnit::TestSuite rejects("rejects");

#ifdef LIBMESH_HAVE_CXX11_REGEX
  // Make regex objects from user's input.
  const std::regex allow_regex(allow_regex_string);
  const std::regex deny_regex(deny_regex_string);

  // Add all tests which match the re to the runner object.
  libMesh::out << "Will run the following tests:" << std::endl;
  const int n_tests_added =
    add_matching_tests_to_runner(registry.makeTest(),
                                 allow_regex_string, allow_regex,
                                 deny_regex_string, deny_regex,
                                 runner, rejects);
  if (n_tests_added >= 0)
    libMesh::out << "--- Running " << n_tests_added << " tests in total." << std::endl;
#else
  // If no C++11 <regex> just run all the tests.
  runner.addTest(registry.makeTest());
#endif

  std::unique_ptr<CppUnit::TestResult> controller;
  std::unique_ptr<CppUnit::BriefTestProgressListener> listener;

  // Actually run all the requested tests, possibly with verbose
  // output of test names as they are run
  if (libMesh::on_command_line("--verbose"))
    {
      listener = std::make_unique<CppUnit::BriefTestProgressListener>();
      runner.eventManager().addListener(listener.get());
    }

  bool succeeded = runner.run();

  // Many users won't care at all about the PerfLog
#ifndef LIBMESH_ENABLE_PERFORMANCE_LOGGING
  if (!libMesh::on_command_line("--full-logs"))
    driver_unitlog.clear();
#endif

  // 1 for failure, 0 for success
  return !succeeded;
}

libMesh::Parallel::Communicator * TestCommWorld;

libMesh::PerfLog * unitlog;
