#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <libmesh/libmesh.h>

#include "test_comm.h"

int main( int argc, char **argv)
{
  // Initialize the library.  This is necessary because the library
  // may depend on a number of other libraries (i.e. MPI  and Petsc)
  // that require initialization before use.
  libMesh::LibMeshInit init(argc, argv);
  TestCommWorld = &init.comm();

  CppUnit::TextUi::TestRunner runner;
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest( registry.makeTest() );

  // If the tests all succeed, report success
  if (runner.run())
    return 0;

  // If any test fails report failure
  return 1;
}

libMesh::Parallel::Communicator *TestCommWorld;
