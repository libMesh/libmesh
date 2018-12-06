
#include "partitioner_test.h"

// If we don't have SFC this should fall back on Linear so we'll test
// heedless of configuration
#include <libmesh/hilbert_sfc_partitioner.h>

// This test is temporarily disabled since it fails sporadically on
// macOS. For more information, see:
// https://github.com/libMesh/libmesh/issues/1967
// INSTANTIATE_PARTITIONER_TEST(HilbertSFCPartitioner,ReplicatedMesh);
