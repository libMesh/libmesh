#ifndef LIBMESH_TESTS_FE_KOKKOS_FE_ORACLE_RUNNERS_H
#define LIBMESH_TESTS_FE_KOKKOS_FE_ORACLE_RUNNERS_H

#include "kokkos_fe_contract_test.K"
#include "kokkos_fe_invariant_test.K"
#include "kokkos_fe_map_oracle_test.K"
#include "kokkos_fe_permuted_map_oracle_test.K"
#include "kokkos_fe_reconstruction_oracle_test.K"
#include "kokkos_fe_shape_oracle_test.K"
#include "kokkos_fe_side_trace_oracle_test.K"
#include "kokkos_fe_types_oracle_test.K"
#include "kokkos_quadrature_oracle_test.K"

namespace libMeshTest
{
namespace KokkosFEOracle
{

inline int
dispatch_child_case(int argc, char ** argv)
{
  return libMeshTest::KokkosFEContractOracle::dispatch_child_case(argc, argv);
}

inline int
run_all_oracles(const char * argv0)
{
  int total_fail = 0;

  total_fail += libMeshTest::KokkosFETypesOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosFEShapeOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosFEMapOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosFEInvariantOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosFEPermutedMapOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosFEReconstructionOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosFESideTraceOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosFEQuadratureOracle::run_all_oracles();
  total_fail += libMeshTest::KokkosFEContractOracle::run_all_oracles(argv0);

  return total_fail;
}

} // namespace KokkosFEOracle
} // namespace libMeshTest

#endif
