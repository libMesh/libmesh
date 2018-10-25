
#include "partitioner_test.h"

// If we don't have METIS this should fall back on SFC or Linear so
// we'll test heedless of configuration
#include <libmesh/metis_partitioner.h>

INSTANTIATE_PARTITIONER_TEST(MetisPartitioner,ReplicatedMesh);
