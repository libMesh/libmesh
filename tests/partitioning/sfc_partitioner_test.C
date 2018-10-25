
#include "partitioner_test.h"

// If we don't have SFC this should fall back on Linear so we'll test
// heedless of configuration
#include <libmesh/sfc_partitioner.h>

INSTANTIATE_PARTITIONER_TEST(SFCPartitioner,ReplicatedMesh);
