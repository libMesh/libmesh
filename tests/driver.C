// Ignore overloaded-virtual warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/TestPath.h>
#include <cppunit/TestResult.h>
#include <libmesh/restore_warnings.h>

// libMesh includes
#include <libmesh/libmesh.h>

#include "libmesh_cppunit.h"
#include "test_comm.h"

#ifdef LIBMESH_HAVE_CXX11_REGEX

// C++ includes
#include <regex>

// We want to be able to add a selection of tests from our full suite
// to our test runner.  But then the test runner will want to destruct
// those tests afterward, so we can't destruct the suite without it
// *also* trying to destruct all its tests (undefined behavior), and
// we can't not destruct the suite without the suite and its
// unselected tests being leaked (which would be fine since just exit
// afterward, except valgrind sees the leaks and screams).
//
// So, instead of adding selected tests, we add shims of selected
// tests that can be deleted without deleting the tests they shim.
class TestShim : public CppUnit::Test
{
public:
  TestShim (CppUnit::Test & shimmed_test) : _shimmed_test(shimmed_test) {}

  virtual void run(CppUnit::TestResult * result) override { _shimmed_test.run(result); }

  virtual int countTestCases () const override { return _shimmed_test.countTestCases(); }

  virtual int getChildTestCount () const override { return _shimmed_test.getChildTestCount(); }

  virtual std::string getName () const override { return _shimmed_test.getName(); }

  virtual bool findTestPath (const std::string & testName, CppUnit::TestPath & testPath) const override { return _shimmed_test.findTestPath(testName, testPath); }

  virtual bool findTestPath (const CppUnit::Test * test, CppUnit::TestPath & testPath) const override { return _shimmed_test.findTestPath(test, testPath); }

  virtual CppUnit::Test * findTest(const std::string & testName) const override { return _shimmed_test.findTest(testName); }

  virtual CppUnit::TestPath resolveTestPath(const std::string & testPath) const override { return _shimmed_test.resolveTestPath(testPath); }

protected:
  virtual CppUnit::Test * doGetChildTestAt(int index) const override { return _shimmed_test.getChildTestAt(index); }

private:
  CppUnit::Test & _shimmed_test;
};

// Add Tests to runner that match user-provided regex.
int add_matching_tests_to_runner(CppUnit::Test * test,
                                 const std::string & allow_r_str,
                                 const std::regex & allow_r,
                                 const std::string & deny_r_str,
                                 const std::regex & deny_r,
                                 CppUnit::TextUi::TestRunner & runner)
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

      // yes, explicit new; this is how CppUnit works
      runner.addTest(new TestShim(*test));
    }
  }

  // Call this recursively on each of our children, if any.
  for (int i = 0; i < test->getChildTestCount(); i++)
    n_tests_added +=
      add_matching_tests_to_runner(test->getChildTestAt(i), allow_r_str, allow_r,
                                   deny_r_str, deny_r, runner);

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

  // We might have to delete the test suite ourselves, after the
  // runner has deleted whatever subtests it has.
  std::unique_ptr<CppUnit::Test> owned_suite;

  // Recursively add tests matching the regex to the runner object.
  CppUnit::TextUi::TestRunner runner;

  // The Cppunit registry object that knows about all the tests, and
  // the test suite it creates.
  CppUnit::TestFactoryRegistry & registry = CppUnit::TestFactoryRegistry::getRegistry();
  CppUnit::Test * suite = registry.makeTest();

#ifdef LIBMESH_HAVE_CXX11_REGEX
  // Make regex objects from user's input.
  const std::regex allow_regex(allow_regex_string);
  const std::regex deny_regex(deny_regex_string);

  // Add all tests which match the re to the runner object.
  libMesh::out << "Will run the following tests:" << std::endl;
  const int n_tests_added =
    add_matching_tests_to_runner(suite,
                                 allow_regex_string, allow_regex,
                                 deny_regex_string, deny_regex,
                                 runner);
  if (n_tests_added >= 0)
    libMesh::out << "--- Running " << n_tests_added << " tests in total." << std::endl;

  // If we didn't add the whole suite to the runner, we need to clean
  // it up ourselves
  if (n_tests_added != -12345)
    owned_suite.reset(suite);
#else
  // If no C++11 <regex> just run all the tests.
  runner.addTest(suite);
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
