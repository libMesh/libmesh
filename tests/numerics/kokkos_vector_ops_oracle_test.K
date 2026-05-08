#include "libmesh/libmesh_config.h"
#include "kokkos_vector_ops_oracle_runners.h"

int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  const int total_fail = libMeshTest::KokkosVectorOracle::run_all_oracles();

  Kokkos::finalize();

  if (total_fail == 0)
    std::printf("ALL TESTS PASSED\n");
  else
    std::printf("%d TEST(S) FAILED\n", total_fail);

  return total_fail ? 1 : 0;
}
