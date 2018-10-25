
#include "partitioner_test.h"

// If we don't have ParMETIS this should fall back on Metis or SFC or
// Linear so we'll test heedless of configuration
#include <libmesh/parmetis_partitioner.h>

INSTANTIATE_PARTITIONER_TEST(ParmetisPartitioner,ReplicatedMesh);
INSTANTIATE_PARTITIONER_TEST(ParmetisPartitioner,DistributedMesh);
