
#include "partitioner_test.h"

// If we don't have SFC this should fall back on Linear so we'll test
// heedless of configuration
#include <libmesh/morton_sfc_partitioner.h>

INSTANTIATE_PARTITIONER_TEST(MortonSFCPartitioner,ReplicatedMesh);
