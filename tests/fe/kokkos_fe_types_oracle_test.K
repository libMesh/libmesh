// GPU kernel oracle tests for libMesh::Kokkos FE type helpers.
// The test suite covers:
//   A. get_side_topology() against libMesh side topology, with the 1D
//      NODEELEM -> EDGE2 surrogate used by the Kokkos path.
//   B. class_from_topology() against a class inferred from libMesh
//      first-order LAGRANGE FE spaces.
//   C. n_dofs(FEShapeKey) for Kokkos-supported exact LAGRANGE keys against
//      libMesh::FEInterface::n_dofs().
//   D. n_dofs(FEShapeKey) for Kokkos-supported MONOMIAL keys against
//      libMesh::FEInterface::n_dofs().
//   E. supports_shape()/supports_grad_shape()/supports_n_dofs() agree on the
//      current Kokkos evaluator support boundary.
//
// Returns 0 on success, non-zero on failure.

#include "gpu/kokkos_fe_types.h"

#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/libmesh.h"

#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

// Avoid conflicting complex operators between CUDA and PETSc
#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__

#include "kokkos_fe_oracle_test_utils.h"

#include <cstdio>
#include <memory>
#include <string>

using kokkos_test_utils::build_reference_elem;

namespace
{

struct side_topology_case
{
  libMesh::ElemType parent;
  libMesh::ElemType expected;
};

struct class_from_topology_case
{
  libMesh::ElemType    topo;
  libMesh::FEElemClass expected;
};

struct n_dof_case
{
  libMesh::Kokkos::FEShapeKey key;
  unsigned int                expected;
};

struct support_case
{
  libMesh::Kokkos::FEShapeKey key;
  bool                        expected;
};

static libMesh::ElemType
host_side_topology_oracle(libMesh::ElemType parent_type)
{
  auto elem = build_reference_elem(parent_type);

  if (elem->dim() == 1)
    return libMesh::EDGE2;

  auto first_side = elem->build_side_ptr(0);
  const libMesh::ElemType side_topo = first_side->type();

  for (unsigned int s = 1; s < elem->n_sides(); ++s)
  {
    auto side = elem->build_side_ptr(s);
    if (side->type() != side_topo)
      return libMesh::INVALID_ELEM;
  }

  return side_topo;
}

static libMesh::FEElemClass
host_class_from_topology_oracle(libMesh::ElemType topo)
{
  auto elem = build_reference_elem(topo);

  if (elem->dim() == 1)
    return libMesh::FEElemClass::EDGE;

  const libMesh::FEType fe_type(libMesh::FIRST, libMesh::LAGRANGE);
  const unsigned int ndofs = libMesh::FEInterface::n_dofs(fe_type, 0, elem.get());

  switch (elem->dim())
  {
    case 2:
      switch (ndofs)
      {
        case 3: return libMesh::FEElemClass::TRI;
        case 4: return libMesh::FEElemClass::QUAD;
        default: return libMesh::FEElemClass::N_CLASSES;
      }

    case 3:
      switch (ndofs)
      {
        case 4: return libMesh::FEElemClass::TET;
        case 8: return libMesh::FEElemClass::HEX;
        case 6: return libMesh::FEElemClass::PRISM;
        case 5: return libMesh::FEElemClass::PYRAMID;
        default: return libMesh::FEElemClass::N_CLASSES;
      }

    default:
      return libMesh::FEElemClass::N_CLASSES;
  }
}

static unsigned int
host_n_dofs_oracle(libMesh::Kokkos::FEShapeKey key)
{
  auto elem = build_reference_elem(key.elem_type);
  return libMesh::FEInterface::n_dofs(libMesh::FEType(key.order, key.family), 0, elem.get());
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Test 1: get_side_topology() against libMesh side topology.
// ---------------------------------------------------------------------------
static int
test_side_topology()
{
  static const libMesh::ElemType parents[] = {
    libMesh::EDGE2,
    libMesh::EDGE3,
    libMesh::EDGE4,
    libMesh::TRI3,
    libMesh::TRI7,
    libMesh::QUAD4,
    libMesh::TRI6,
    libMesh::QUAD8,
    libMesh::QUAD9,
    libMesh::TET4,
    libMesh::HEX8,
    libMesh::TET10,
    libMesh::TET14,
    libMesh::HEX20,
    libMesh::HEX27
  };
  constexpr int n_cases = sizeof(parents) / sizeof(parents[0]);

  Kokkos::View<side_topology_case *> d_cases(std::string("side_cases"), n_cases);
  {
    auto h = Kokkos::create_mirror_view(d_cases);
    for (int i = 0; i < n_cases; ++i)
    {
      h(i).parent = parents[i];
      h(i).expected = host_side_topology_oracle(parents[i]);
    }
    Kokkos::deep_copy(d_cases, h);
  }

  Kokkos::View<int> d_fail(std::string("side_fail"));
  Kokkos::deep_copy(d_fail, 0);

  Kokkos::parallel_for(
    n_cases,
    KOKKOS_LAMBDA(int i) {
      using namespace libMesh::Kokkos;
      if (get_side_topology(d_cases(i).parent) != d_cases(i).expected)
        Kokkos::atomic_add(&d_fail(), 1);
    });
  Kokkos::fence();

  int fail = 0;
  Kokkos::deep_copy(fail, d_fail);
  return fail;
}

// ---------------------------------------------------------------------------
// Test 2: class_from_topology() against libMesh FE oracle.
// ---------------------------------------------------------------------------
static int
test_class_from_topology()
{
  static const libMesh::ElemType topos[] = {
    libMesh::EDGE2,
    libMesh::EDGE3,
    libMesh::EDGE4,
    libMesh::TRI3,
    libMesh::TRI6,
    libMesh::TRI7,
    libMesh::QUAD4,
    libMesh::QUAD8,
    libMesh::QUAD9,
    libMesh::TET4,
    libMesh::TET10,
    libMesh::TET14,
    libMesh::HEX8,
    libMesh::HEX20,
    libMesh::HEX27,
    libMesh::PRISM6,
    libMesh::PRISM15,
    libMesh::PRISM18,
    libMesh::PRISM20,
    libMesh::PRISM21,
    libMesh::PYRAMID5,
    libMesh::PYRAMID13,
    libMesh::PYRAMID14,
    libMesh::PYRAMID18
  };
  constexpr int n_cases = sizeof(topos) / sizeof(topos[0]);

  Kokkos::View<class_from_topology_case *> d_cases(std::string("class_cases"), n_cases);
  {
    auto h = Kokkos::create_mirror_view(d_cases);
    for (int i = 0; i < n_cases; ++i)
    {
      h(i).topo = topos[i];
      h(i).expected = host_class_from_topology_oracle(topos[i]);
    }
    Kokkos::deep_copy(d_cases, h);
  }

  Kokkos::View<int> d_fail(std::string("class_fail"));
  Kokkos::deep_copy(d_fail, 0);

  Kokkos::parallel_for(
    n_cases,
    KOKKOS_LAMBDA(int i) {
      using namespace libMesh::Kokkos;
      if (class_from_topology(d_cases(i).topo) != d_cases(i).expected)
        Kokkos::atomic_add(&d_fail(), 1);
    });
  Kokkos::fence();

  int fail = 0;
  Kokkos::deep_copy(fail, d_fail);
  return fail;
}

// ---------------------------------------------------------------------------
// Test 3: n_dofs() for Kokkos-supported exact LAGRANGE keys against
// libMesh FEInterface.
// ---------------------------------------------------------------------------
static int
test_n_dofs_lagrange()
{
  using libMesh::Kokkos::FEShapeKey;

  static const FEShapeKey keys[] = {
    { libMesh::LAGRANGE, libMesh::EDGE2,     libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::EDGE3,     libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::EDGE3,     libMesh::SECOND   },
    { libMesh::LAGRANGE, libMesh::EDGE4,     libMesh::FIRST    },

    { libMesh::LAGRANGE, libMesh::TRI3,      libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::TRI6,      libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::TRI6,      libMesh::SECOND   },
    { libMesh::LAGRANGE, libMesh::TRI7,      libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::TRI7,      libMesh::SECOND   },
    { libMesh::LAGRANGE, libMesh::QUAD4,     libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::QUAD8,     libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::QUAD8,     libMesh::SECOND   },
    { libMesh::LAGRANGE, libMesh::QUAD9,     libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::QUAD9,     libMesh::SECOND   },

    { libMesh::LAGRANGE, libMesh::TET4,      libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::TET10,     libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::TET10,     libMesh::SECOND   },
    { libMesh::LAGRANGE, libMesh::TET14,     libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::TET14,     libMesh::SECOND   },
    { libMesh::LAGRANGE, libMesh::HEX8,      libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::HEX20,     libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::HEX20,     libMesh::SECOND   },
    { libMesh::LAGRANGE, libMesh::HEX27,     libMesh::FIRST    },
    { libMesh::LAGRANGE, libMesh::HEX27,     libMesh::SECOND   }
  };
  constexpr int n_cases = sizeof(keys) / sizeof(keys[0]);

  Kokkos::View<n_dof_case *> d_cases(std::string("lagrange_cases"), n_cases);
  {
    auto h = Kokkos::create_mirror_view(d_cases);
    for (int i = 0; i < n_cases; ++i)
    {
      h(i).key = keys[i];
      h(i).expected = host_n_dofs_oracle(keys[i]);
    }
    Kokkos::deep_copy(d_cases, h);
  }

  Kokkos::View<int> d_fail(std::string("lagrange_fail"));
  Kokkos::deep_copy(d_fail, 0);

  Kokkos::parallel_for(
    n_cases,
    KOKKOS_LAMBDA(int i) {
      using namespace libMesh::Kokkos;
      if (n_dofs(d_cases(i).key) != d_cases(i).expected)
        Kokkos::atomic_add(&d_fail(), 1);
    });
  Kokkos::fence();

  int fail = 0;
  Kokkos::deep_copy(fail, d_fail);
  return fail;
}

// ---------------------------------------------------------------------------
// Test 4: n_dofs() for Kokkos-supported MONOMIAL keys against libMesh
// FEInterface.
// ---------------------------------------------------------------------------
static int
test_n_dofs_monomial()
{
  using libMesh::Kokkos::FEShapeKey;

  static const libMesh::Order orders[] = {
    libMesh::CONSTANT,
    libMesh::FIRST,
    libMesh::SECOND,
    libMesh::THIRD,
    libMesh::FOURTH,
    libMesh::FIFTH
  };
  static const libMesh::ElemType higher_dim_topos[] = {
    libMesh::TRI7,
    libMesh::QUAD9,
    libMesh::TET14,
    libMesh::HEX27,
    libMesh::PRISM21
  };

  constexpr int n_cases = sizeof(orders) / sizeof(orders[0]) *
                          (2 + sizeof(higher_dim_topos) / sizeof(higher_dim_topos[0]));

  Kokkos::View<n_dof_case *> d_cases(std::string("monomial_cases"), n_cases);
  {
    auto h = Kokkos::create_mirror_view(d_cases);
    int i = 0;
    for (const auto order : orders)
    {
      const libMesh::ElemType edge_topo =
        (order <= libMesh::THIRD) ? libMesh::EDGE4 : libMesh::EDGE3;
      const libMesh::ElemType pyramid_topo =
        (order <= libMesh::THIRD) ? libMesh::PYRAMID18 : libMesh::PYRAMID14;

      h(i).key = { libMesh::MONOMIAL, edge_topo, order };
      h(i).expected = host_n_dofs_oracle(h(i).key);
      ++i;

      h(i).key = { libMesh::MONOMIAL, pyramid_topo, order };
      h(i).expected = host_n_dofs_oracle(h(i).key);
      ++i;

      for (const auto topo : higher_dim_topos)
      {
        h(i).key = { libMesh::MONOMIAL, topo, order };
        h(i).expected = host_n_dofs_oracle(h(i).key);
        ++i;
      }
    }
    Kokkos::deep_copy(d_cases, h);
  }

  Kokkos::View<int> d_fail(std::string("monomial_fail"));
  Kokkos::deep_copy(d_fail, 0);

  Kokkos::parallel_for(
    n_cases,
    KOKKOS_LAMBDA(int i) {
      using namespace libMesh::Kokkos;
      if (n_dofs(d_cases(i).key) != d_cases(i).expected)
        Kokkos::atomic_add(&d_fail(), 1);
    });
  Kokkos::fence();

  int fail = 0;
  Kokkos::deep_copy(fail, d_fail);
  return fail;
}

// ---------------------------------------------------------------------------
// Test 5: support predicates agree on the Kokkos evaluator boundary.
// ---------------------------------------------------------------------------
static int
test_support_contract()
{
  using libMesh::Kokkos::FEShapeKey;

  static const support_case cases[] = {
    { { libMesh::LAGRANGE, libMesh::EDGE2,     libMesh::FIRST    }, true  },
    { { libMesh::LAGRANGE, libMesh::EDGE3,     libMesh::SECOND   }, true  },
    { { libMesh::LAGRANGE, libMesh::TRI7,      libMesh::SECOND   }, true  },
    { { libMesh::LAGRANGE, libMesh::QUAD9,     libMesh::SECOND   }, true  },
    { { libMesh::LAGRANGE, libMesh::TET14,     libMesh::SECOND   }, true  },
    { { libMesh::LAGRANGE, libMesh::HEX27,     libMesh::SECOND   }, true  },
    { { libMesh::LAGRANGE, libMesh::NODEELEM,  libMesh::CONSTANT }, false },
    { { libMesh::LAGRANGE, libMesh::NODEELEM,  libMesh::FIRST    }, false },
    { { libMesh::LAGRANGE, libMesh::EDGE4,     libMesh::THIRD    }, false },
    { { libMesh::LAGRANGE, libMesh::TRI7,      libMesh::THIRD    }, false },
    { { libMesh::LAGRANGE, libMesh::TET14,     libMesh::THIRD    }, false },
    { { libMesh::LAGRANGE, libMesh::PRISM6,    libMesh::FIRST    }, false },
    { { libMesh::LAGRANGE, libMesh::PRISM15,   libMesh::SECOND   }, false },
    { { libMesh::LAGRANGE, libMesh::PYRAMID5,  libMesh::FIRST    }, false },
    { { libMesh::LAGRANGE, libMesh::PYRAMID14, libMesh::SECOND   }, false },
    { { libMesh::LAGRANGE, libMesh::PYRAMID18, libMesh::THIRD    }, false },
    { { libMesh::LAGRANGE, libMesh::EDGE2,     libMesh::INVALID_ORDER }, false },
    { { libMesh::MONOMIAL, libMesh::EDGE4,     libMesh::THIRD    }, true  },
    { { libMesh::MONOMIAL, libMesh::TRI7,      libMesh::FIFTH    }, true  },
    { { libMesh::MONOMIAL, libMesh::QUAD9,     libMesh::FIFTH    }, true  },
    { { libMesh::MONOMIAL, libMesh::TET14,     libMesh::FIFTH    }, true  },
    { { libMesh::MONOMIAL, libMesh::PRISM21,   libMesh::FIFTH    }, true  },
    { { libMesh::MONOMIAL, libMesh::PYRAMID14, libMesh::FIFTH    }, true  },
    { { libMesh::MONOMIAL, libMesh::NODEELEM,  libMesh::CONSTANT }, false },
    { { libMesh::MONOMIAL, libMesh::EDGE4,     libMesh::FOURTH   }, false },
    { { libMesh::MONOMIAL, libMesh::PYRAMID18, libMesh::FOURTH   }, false },
    { { libMesh::MONOMIAL, libMesh::TRI7,      libMesh::SIXTH    }, false },
    { { libMesh::MONOMIAL, libMesh::HEX27,     libMesh::SIXTH    }, false },
    { { libMesh::MONOMIAL, libMesh::EDGE2,     libMesh::INVALID_ORDER }, false }
  };
  constexpr int n_cases = sizeof(cases) / sizeof(cases[0]);

  Kokkos::View<support_case *> d_cases(std::string("support_cases"), n_cases);
  {
    auto h = Kokkos::create_mirror_view(d_cases);
    for (int i = 0; i < n_cases; ++i)
      h(i) = cases[i];
    Kokkos::deep_copy(d_cases, h);
  }

  Kokkos::View<int> d_fail(std::string("support_fail"));
  Kokkos::deep_copy(d_fail, 0);

  Kokkos::parallel_for(
    n_cases,
    KOKKOS_LAMBDA(int i) {
      using namespace libMesh::Kokkos;

      const bool shape_supported = supports_shape(d_cases(i).key);
      const bool grad_supported = supports_grad_shape(d_cases(i).key);
      const bool ndofs_supported = supports_n_dofs(d_cases(i).key);

      if (shape_supported != d_cases(i).expected ||
          grad_supported != d_cases(i).expected ||
          ndofs_supported != d_cases(i).expected ||
          shape_supported != grad_supported ||
          shape_supported != ndofs_supported)
        Kokkos::atomic_add(&d_fail(), 1);
    });
  Kokkos::fence();

  int fail = 0;
  Kokkos::deep_copy(fail, d_fail);
  return fail;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  libMesh::LibMeshInit init(argc, argv);

  int total_fail = 0;

  {
    const int f = test_side_topology();
    std::printf("[side_topology_oracle]     %s  (%d failures)\n", f ? "FAIL" : "PASS", f);
    total_fail += f;
  }
  {
    const int f = test_class_from_topology();
    std::printf("[class_topology_oracle]    %s  (%d failures)\n", f ? "FAIL" : "PASS", f);
    total_fail += f;
  }
  {
    const int f = test_n_dofs_lagrange();
    std::printf("[lagrange_ndofs_oracle]    %s  (%d failures)\n", f ? "FAIL" : "PASS", f);
    total_fail += f;
  }
  {
    const int f = test_n_dofs_monomial();
    std::printf("[monomial_ndofs_oracle]    %s  (%d failures)\n", f ? "FAIL" : "PASS", f);
    total_fail += f;
  }
  {
    const int f = test_support_contract();
    std::printf("[support_contract_oracle]  %s  (%d failures)\n", f ? "FAIL" : "PASS", f);
    total_fail += f;
  }

  Kokkos::finalize();

  if (total_fail == 0)
    std::printf("ALL TESTS PASSED\n");
  else
    std::printf("%d TEST(S) FAILED\n", total_fail);

  return total_fail ? 1 : 0;
}
