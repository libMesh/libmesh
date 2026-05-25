#include "libmesh/libmesh_config.h"
#include "fe/kokkos_fe_oracle_runners.h"
#include "numerics/kokkos_tensor_ops_oracle_runners.h"
#include "numerics/kokkos_vector_ops_oracle_runners.h"

#include <cstdio>

int
main(int argc, char ** argv)
{
  const int child_result = libMeshTest::KokkosFEOracle::dispatch_child_case(argc, argv);
  if (child_result >= 0)
    return child_result;

  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int total_fail = 0;
  total_fail += libMeshTest::KokkosVectorOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosTensorOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosFEOracle::run_all_oracles(argv[0]);

  Kokkos::finalize();

  if (total_fail == 0)
    std::printf("ALL TESTS PASSED\n");
  else
    std::printf("%d TEST(S) FAILED\n", total_fail);

  return total_fail ? 1 : 0;
}
