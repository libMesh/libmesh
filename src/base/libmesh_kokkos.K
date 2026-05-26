#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_KOKKOS

// Avoid conflicting complex operators between CUDA and PETSc
#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__

namespace libMesh
{

bool
libmesh_kokkos_initialize_if_needed(int & argc, char ** & argv)
{
  if (Kokkos::is_initialized())
    return false;

  Kokkos::initialize(argc, argv);
  return true;
}

void
libmesh_kokkos_finalize_if_needed(const bool should_finalize)
{
  if (should_finalize && Kokkos::is_initialized())
    Kokkos::finalize();
}

} // namespace libMesh

#endif
