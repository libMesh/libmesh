// Standalone driver for the Kokkos numerics "oracle" tests.
//
// This is a .K translation unit so it gets compiled with the Kokkos
// toolchain via the custom rules in tests/Makefile.am, independently of
// the CppUnit-based unit_tests driver (these tests need Kokkos but not
// CppUnit).  The oracle suites compute results on the device through
// the gpu/ headers and compare them against references computed on the
// host with the owning libMesh types, along with a few host-only
// checks; run_unit_tests.sh runs one such driver per enabled METHOD,
// and a nonzero exit status marks the run as failed.

#include "libmesh/libmesh_config.h"
#include "numerics/kokkos_tensor_ops_oracle_runners.h"
#include "numerics/kokkos_vector_ops_oracle_runners.h"

#include <cstdio>

int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int total_fail = 0;
  total_fail += libMeshTest::KokkosVectorOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosTensorOracle::run_all_oracles();

  Kokkos::finalize();

  if (total_fail == 0)
    std::printf("ALL TESTS PASSED\n");
  else
    std::printf("%d TEST(S) FAILED\n", total_fail);

  return total_fail ? 1 : 0;
}
