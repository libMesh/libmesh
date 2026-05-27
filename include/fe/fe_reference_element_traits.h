// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

#ifndef LIBMESH_FE_REFERENCE_ELEMENT_TRAITS_H
#define LIBMESH_FE_REFERENCE_ELEMENT_TRAITS_H

#include "libmesh/enum_elem_type.h"
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_device.h"
#include "libmesh/point.h"

namespace libMesh
{

template <unsigned int N>
struct ReferenceElementVector
{
  unsigned int values[N];

  LIBMESH_DEVICE_INLINE constexpr unsigned int operator[](unsigned int i) const
  { return values[i]; }
};

template <unsigned int Rows, unsigned int Cols>
struct ReferenceElementTable
{
  unsigned int values[Rows][Cols];

  LIBMESH_DEVICE_INLINE constexpr const unsigned int * operator[](unsigned int i) const
  { return values[i]; }

  LIBMESH_DEVICE_INLINE constexpr unsigned int operator()(unsigned int i, unsigned int j) const
  { return values[i][j]; }
};


LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<2>
edge2_side_node_counts()
{
  return {{1, 1}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<2>
edge3_side_node_counts()
{
  return {{1, 1}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<2>
edge4_side_node_counts()
{
  return {{1, 1}};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<3>
tri3_side_node_counts()
{
  return {{2, 2, 2}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<3>
tri6_side_node_counts()
{
  return {{3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<3>
tri7_side_node_counts()
{
  return {{3, 3, 3}};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<4>
quad4_side_node_counts()
{
  return {{2, 2, 2, 2}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<4>
quad8_side_node_counts()
{
  return {{3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<4>
quad9_side_node_counts()
{
  return {{3, 3, 3, 3}};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<4>
tet4_side_node_counts()
{
  return {{3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<4>
tet10_side_node_counts()
{
  return {{6, 6, 6, 6}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<4>
tet14_side_node_counts()
{
  return {{7, 7, 7, 7}};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<6>
hex8_side_node_counts()
{
  return {{4, 4, 4, 4, 4, 4}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<6>
hex20_side_node_counts()
{
  return {{8, 8, 8, 8, 8, 8}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<6>
hex27_side_node_counts()
{
  return {{9, 9, 9, 9, 9, 9}};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<5>
prism6_side_node_counts()
{
  return {{3, 4, 4, 4, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<5>
prism15_side_node_counts()
{
  return {{6, 8, 8, 8, 6}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<5>
prism18_side_node_counts()
{
  return {{6, 9, 9, 9, 6}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<5>
prism20_side_node_counts()
{
  return {{7, 9, 9, 9, 7}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<5>
prism21_side_node_counts()
{
  return {{7, 9, 9, 9, 7}};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<5>
pyramid5_side_node_counts()
{
  return {{3, 3, 3, 3, 4}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<5>
pyramid13_side_node_counts()
{
  return {{6, 6, 6, 6, 8}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<5>
pyramid14_side_node_counts()
{
  return {{6, 6, 6, 6, 9}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<5>
pyramid18_side_node_counts()
{
  return {{7, 7, 7, 7, 9}};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<6>
tet10_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<6>
tet14_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<12>
hex20_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<12>
hex27_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<9>
prism15_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<9>
prism18_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<9>
prism20_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<9>
prism21_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<8>
pyramid13_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<8>
pyramid14_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3, 3, 3}};
}
LIBMESH_DEVICE_INLINE constexpr ReferenceElementVector<8>
pyramid18_edge_node_counts()
{
  return {{3, 3, 3, 3, 3, 3, 3, 3}};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<5, 4>
prism6_side_nodes()
{
  return {{
    {0, 2, 1, 99},
    {0, 1, 4, 3},
    {1, 2, 5, 4},
    {2, 0, 3, 5},
    {3, 4, 5, 99}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<5, 8>
prism15_side_nodes()
{
  return {{
    {0, 2, 1, 8, 7, 6, 99, 99},
    {0, 1, 4, 3, 6, 10, 12, 9},
    {1, 2, 5, 4, 7, 11, 13, 10},
    {2, 0, 3, 5, 8, 9, 14, 11},
    {3, 4, 5, 12, 13, 14, 99, 99}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<5, 9>
prism18_side_nodes()
{
  return {{
    {0, 2, 1, 8, 7, 6, 99, 99, 99},
    {0, 1, 4, 3, 6, 10, 12, 9, 15},
    {1, 2, 5, 4, 7, 11, 13, 10, 16},
    {2, 0, 3, 5, 8, 9, 14, 11, 17},
    {3, 4, 5, 12, 13, 14, 99, 99, 99}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<5, 9>
prism20_side_nodes()
{
  return {{
    {0, 2, 1, 8, 7, 6, 18, 99, 99},
    {0, 1, 4, 3, 6, 10, 12, 9, 15},
    {1, 2, 5, 4, 7, 11, 13, 10, 16},
    {2, 0, 3, 5, 8, 9, 14, 11, 17},
    {3, 4, 5, 12, 13, 14, 19, 99, 99}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<5, 9>
prism21_side_nodes()
{
  return {{
    {0, 2, 1, 8, 7, 6, 18, 99, 99},
    {0, 1, 4, 3, 6, 10, 12, 9, 15},
    {1, 2, 5, 4, 7, 11, 13, 10, 16},
    {2, 0, 3, 5, 8, 9, 14, 11, 17},
    {3, 4, 5, 12, 13, 14, 19, 99, 99}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<5, 4>
pyramid5_side_nodes()
{
  return {{
    {0, 1, 4, 99},
    {1, 2, 4, 99},
    {2, 3, 4, 99},
    {3, 0, 4, 99},
    {0, 3, 2, 1}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<5, 8>
pyramid13_side_nodes()
{
  return {{
    {0, 1, 4, 5, 10, 9, 99, 99},
    {1, 2, 4, 6, 11, 10, 99, 99},
    {2, 3, 4, 7, 12, 11, 99, 99},
    {3, 0, 4, 8, 9, 12, 99, 99},
    {0, 3, 2, 1, 8, 7, 6, 5}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<5, 9>
pyramid14_side_nodes()
{
  return {{
    {0, 1, 4, 5, 10, 9, 99, 99, 99},
    {1, 2, 4, 6, 11, 10, 99, 99, 99},
    {2, 3, 4, 7, 12, 11, 99, 99, 99},
    {3, 0, 4, 8, 9, 12, 99, 99, 99},
    {0, 3, 2, 1, 8, 7, 6, 5, 13}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<5, 9>
pyramid18_side_nodes()
{
  return {{
    {0, 1, 4, 5, 10, 9, 14, 99, 99},
    {1, 2, 4, 6, 11, 10, 15, 99, 99},
    {2, 3, 4, 7, 12, 11, 16, 99, 99},
    {3, 0, 4, 8, 9, 12, 17, 99, 99},
    {0, 3, 2, 1, 8, 7, 6, 5, 13}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<3, 2>
tri3_side_nodes()
{
  return {{
    {0, 1},
    {1, 2},
    {2, 0}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<3, 3>
tri6_side_nodes()
{
  return {{
    {0, 1, 3},
    {1, 2, 4},
    {2, 0, 5}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<3, 3>
tri7_side_nodes()
{
  return {{
    {0, 1, 3},
    {1, 2, 4},
    {2, 0, 5}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<4, 2>
quad4_side_nodes()
{
  return {{
    {0, 1},
    {1, 2},
    {2, 3},
    {3, 0}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<4, 3>
quad8_side_nodes()
{
  return {{
    {0, 1, 4},
    {1, 2, 5},
    {2, 3, 6},
    {3, 0, 7}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<4, 3>
quad9_side_nodes()
{
  return {{
    {0, 1, 4},
    {1, 2, 5},
    {2, 3, 6},
    {3, 0, 7}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<4, 3>
tet4_side_nodes()
{
  return {{
    {0, 2, 1},
    {0, 1, 3},
    {1, 2, 3},
    {2, 0, 3}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<4, 6>
tet10_side_nodes()
{
  return {{
    {0, 2, 1, 6, 5, 4},
    {0, 1, 3, 4, 8, 7},
    {1, 2, 3, 5, 9, 8},
    {2, 0, 3, 6, 7, 9}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<4, 7>
tet14_side_nodes()
{
  return {{
    {0, 2, 1, 6, 5, 4, 10},
    {0, 1, 3, 4, 8, 7, 11},
    {1, 2, 3, 5, 9, 8, 12},
    {2, 0, 3, 6, 7, 9, 13}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<6, 4>
hex8_side_nodes()
{
  return {{
    {0, 3, 2, 1},
    {0, 1, 5, 4},
    {1, 2, 6, 5},
    {2, 3, 7, 6},
    {3, 0, 4, 7},
    {4, 5, 6, 7}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<6, 8>
hex20_side_nodes()
{
  return {{
    {0, 3, 2, 1, 11, 10,  9,  8},
    {0, 1, 5, 4,  8, 13, 16, 12},
    {1, 2, 6, 5,  9, 14, 17, 13},
    {2, 3, 7, 6, 10, 15, 18, 14},
    {3, 0, 4, 7, 11, 12, 19, 15},
    {4, 5, 6, 7, 16, 17, 18, 19}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<6, 9>
hex27_side_nodes()
{
  return {{
    {0, 3, 2, 1, 11, 10,  9,  8, 20},
    {0, 1, 5, 4,  8, 13, 16, 12, 21},
    {1, 2, 6, 5,  9, 14, 17, 13, 22},
    {2, 3, 7, 6, 10, 15, 18, 14, 23},
    {3, 0, 4, 7, 11, 12, 19, 15, 24},
    {4, 5, 6, 7, 16, 17, 18, 19, 25}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<2, 1>
edge2_side_nodes()
{
  return {{
    {0},
    {1}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<2, 1>
edge3_side_nodes()
{
  return {{
    {0},
    {1}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<2, 1>
edge4_side_nodes()
{
  return {{
    {0},
    {1}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<6, 3>
tet10_edge_nodes()
{
  return {{
    {0, 1, 4},
    {1, 2, 5},
    {0, 2, 6},
    {0, 3, 7},
    {1, 3, 8},
    {2, 3, 9}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<6, 3>
tet14_edge_nodes()
{
  return {{
    {0, 1, 4},
    {1, 2, 5},
    {0, 2, 6},
    {0, 3, 7},
    {1, 3, 8},
    {2, 3, 9}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<12, 3>
hex20_edge_nodes()
{
  return {{
    {0, 1, 8},
    {1, 2, 9},
    {2, 3, 10},
    {0, 3, 11},
    {0, 4, 12},
    {1, 5, 13},
    {2, 6, 14},
    {3, 7, 15},
    {4, 5, 16},
    {5, 6, 17},
    {6, 7, 18},
    {4, 7, 19}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<12, 3>
hex27_edge_nodes()
{
  return {{
    {0, 1, 8},
    {1, 2, 9},
    {2, 3, 10},
    {0, 3, 11},
    {0, 4, 12},
    {1, 5, 13},
    {2, 6, 14},
    {3, 7, 15},
    {4, 5, 16},
    {5, 6, 17},
    {6, 7, 18},
    {4, 7, 19}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<9, 3>
prism15_edge_nodes()
{
  return {{
    {0, 1, 6},
    {1, 2, 7},
    {0, 2, 8},
    {0, 3, 9},
    {1, 4, 10},
    {2, 5, 11},
    {3, 4, 12},
    {4, 5, 13},
    {3, 5, 14}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<9, 3>
prism18_edge_nodes()
{
  return {{
    {0, 1, 6},
    {1, 2, 7},
    {0, 2, 8},
    {0, 3, 9},
    {1, 4, 10},
    {2, 5, 11},
    {3, 4, 12},
    {4, 5, 13},
    {3, 5, 14}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<9, 3>
prism20_edge_nodes()
{
  return {{
    {0, 1, 6},
    {1, 2, 7},
    {0, 2, 8},
    {0, 3, 9},
    {1, 4, 10},
    {2, 5, 11},
    {3, 4, 12},
    {4, 5, 13},
    {3, 5, 14}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<9, 3>
prism21_edge_nodes()
{
  return {{
    {0, 1, 6},
    {1, 2, 7},
    {0, 2, 8},
    {0, 3, 9},
    {1, 4, 10},
    {2, 5, 11},
    {3, 4, 12},
    {4, 5, 13},
    {3, 5, 14}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<8, 3>
pyramid13_edge_nodes()
{
  return {{
    {0, 1, 5},
    {1, 2, 6},
    {2, 3, 7},
    {0, 3, 8},
    {0, 4, 9},
    {1, 4, 10},
    {2, 4, 11},
    {3, 4, 12}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<8, 3>
pyramid14_edge_nodes()
{
  return {{
    {0, 1, 5},
    {1, 2, 6},
    {2, 3, 7},
    {0, 3, 8},
    {0, 4, 9},
    {1, 4, 10},
    {2, 4, 11},
    {3, 4, 12}
  }};
}

LIBMESH_DEVICE_INLINE constexpr ReferenceElementTable<8, 3>
pyramid18_edge_nodes()
{
  return {{
    {0, 1, 5},
    {1, 2, 6},
    {2, 3, 7},
    {0, 3, 8},
    {0, 4, 9},
    {1, 4, 10},
    {2, 4, 11},
    {3, 4, 12}
  }};
}

LIBMESH_DEVICE_INLINE bool
requires_side_specific_topology(ElemType parent)
{
  switch (parent)
  {
    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      return true;
    default:
      return false;
  }
}

LIBMESH_DEVICE_INLINE ElemType
side_topology_or_invalid(ElemType parent,
                         unsigned int side)
{
  switch (parent)
  {
    case PRISM6:
      switch (side)
      {
        case 0:
        case 4:
          return TRI3;
        case 1:
        case 2:
        case 3:
          return QUAD4;
        default:
          return INVALID_ELEM;
      }

    case PRISM15:
      switch (side)
      {
        case 0:
        case 4:
          return TRI6;
        case 1:
        case 2:
        case 3:
          return QUAD8;
        default:
          return INVALID_ELEM;
      }

    case PRISM18:
      switch (side)
      {
        case 0:
        case 4:
          return TRI6;
        case 1:
        case 2:
        case 3:
          return QUAD9;
        default:
          return INVALID_ELEM;
      }

    case PRISM20:
    case PRISM21:
      switch (side)
      {
        case 0:
        case 4:
          return TRI7;
        case 1:
        case 2:
        case 3:
          return QUAD9;
        default:
          return INVALID_ELEM;
      }

    case PYRAMID5:
      switch (side)
      {
        case 0:
        case 1:
        case 2:
        case 3:
          return TRI3;
        case 4:
          return QUAD4;
        default:
          return INVALID_ELEM;
      }

    case PYRAMID13:
      switch (side)
      {
        case 0:
        case 1:
        case 2:
        case 3:
          return TRI6;
        case 4:
          return QUAD8;
        default:
          return INVALID_ELEM;
      }

    case PYRAMID14:
      switch (side)
      {
        case 0:
        case 1:
        case 2:
        case 3:
          return TRI6;
        case 4:
          return QUAD9;
        default:
          return INVALID_ELEM;
      }

    case PYRAMID18:
      switch (side)
      {
        case 0:
        case 1:
        case 2:
        case 3:
          return TRI7;
        case 4:
          return QUAD9;
        default:
          return INVALID_ELEM;
      }

    default:
      return INVALID_ELEM;
  }
}

LIBMESH_DEVICE_INLINE unsigned int
side_node_count_or_zero(ElemType parent,
                        unsigned int side)
{
  switch (parent)
  {
    case EDGE2:
      return side < 2 ? edge2_side_node_counts()[side] : 0;
    case EDGE3:
      return side < 2 ? edge3_side_node_counts()[side] : 0;
    case EDGE4:
      return side < 2 ? edge4_side_node_counts()[side] : 0;
    case TRI3:
    case TRISHELL3:
      return side < 3 ? tri3_side_node_counts()[side] : 0;
    case TRI6:
      return side < 3 ? tri6_side_node_counts()[side] : 0;
    case TRI7:
      return side < 3 ? tri7_side_node_counts()[side] : 0;
    case QUAD4:
    case QUADSHELL4:
      return side < 4 ? quad4_side_node_counts()[side] : 0;
    case QUAD8:
    case QUADSHELL8:
      return side < 4 ? quad8_side_node_counts()[side] : 0;
    case QUAD9:
    case QUADSHELL9:
      return side < 4 ? quad9_side_node_counts()[side] : 0;
    case TET4:
      return side < 4 ? tet4_side_node_counts()[side] : 0;
    case TET10:
      return side < 4 ? tet10_side_node_counts()[side] : 0;
    case TET14:
      return side < 4 ? tet14_side_node_counts()[side] : 0;
    case HEX8:
      return side < 6 ? hex8_side_node_counts()[side] : 0;
    case HEX20:
      return side < 6 ? hex20_side_node_counts()[side] : 0;
    case HEX27:
      return side < 6 ? hex27_side_node_counts()[side] : 0;
    case PRISM6:
      return side < 5 ? prism6_side_node_counts()[side] : 0;
    case PRISM15:
      return side < 5 ? prism15_side_node_counts()[side] : 0;
    case PRISM18:
      return side < 5 ? prism18_side_node_counts()[side] : 0;
    case PRISM20:
      return side < 5 ? prism20_side_node_counts()[side] : 0;
    case PRISM21:
      return side < 5 ? prism21_side_node_counts()[side] : 0;
    case PYRAMID5:
      return side < 5 ? pyramid5_side_node_counts()[side] : 0;
    case PYRAMID13:
      return side < 5 ? pyramid13_side_node_counts()[side] : 0;
    case PYRAMID14:
      return side < 5 ? pyramid14_side_node_counts()[side] : 0;
    case PYRAMID18:
      return side < 5 ? pyramid18_side_node_counts()[side] : 0;
    default:
      return 0;
  }
}

LIBMESH_DEVICE_INLINE unsigned int
edge_node_count_or_zero(ElemType parent,
                        unsigned int edge)
{
  switch (parent)
  {
    case TET10:
      return edge < 6 ? tet10_edge_node_counts()[edge] : 0;
    case TET14:
      return edge < 6 ? tet14_edge_node_counts()[edge] : 0;
    case HEX20:
      return edge < 12 ? hex20_edge_node_counts()[edge] : 0;
    case HEX27:
      return edge < 12 ? hex27_edge_node_counts()[edge] : 0;
    case PRISM15:
      return edge < 9 ? prism15_edge_node_counts()[edge] : 0;
    case PRISM18:
      return edge < 9 ? prism18_edge_node_counts()[edge] : 0;
    case PRISM20:
      return edge < 9 ? prism20_edge_node_counts()[edge] : 0;
    case PRISM21:
      return edge < 9 ? prism21_edge_node_counts()[edge] : 0;
    case PYRAMID13:
      return edge < 8 ? pyramid13_edge_node_counts()[edge] : 0;
    case PYRAMID14:
      return edge < 8 ? pyramid14_edge_node_counts()[edge] : 0;
    case PYRAMID18:
      return edge < 8 ? pyramid18_edge_node_counts()[edge] : 0;
    default:
      return 0;
  }
}

LIBMESH_DEVICE_INLINE bool
try_local_side_node(ElemType parent,
                    unsigned int side,
                    unsigned int side_node,
                    unsigned int & node)
{
  const unsigned int count = side_node_count_or_zero(parent, side);
  if (!count || side_node >= count)
    return false;

  switch (parent)
  {
    case EDGE2:
      node = edge2_side_nodes()(side, side_node);
      return true;
    case EDGE3:
      node = edge3_side_nodes()(side, side_node);
      return true;
    case EDGE4:
      node = edge4_side_nodes()(side, side_node);
      return true;
    case TRI3:
    case TRISHELL3:
      node = tri3_side_nodes()(side, side_node);
      return true;
    case TRI6:
      node = tri6_side_nodes()(side, side_node);
      return true;
    case TRI7:
      node = tri7_side_nodes()(side, side_node);
      return true;
    case QUAD4:
    case QUADSHELL4:
      node = quad4_side_nodes()(side, side_node);
      return true;
    case QUAD8:
    case QUADSHELL8:
      node = quad8_side_nodes()(side, side_node);
      return true;
    case QUAD9:
    case QUADSHELL9:
      node = quad9_side_nodes()(side, side_node);
      return true;
    case TET4:
      node = tet4_side_nodes()(side, side_node);
      return true;
    case TET10:
      node = tet10_side_nodes()(side, side_node);
      return true;
    case TET14:
      node = tet14_side_nodes()(side, side_node);
      return true;
    case HEX8:
      node = hex8_side_nodes()(side, side_node);
      return true;
    case HEX20:
      node = hex20_side_nodes()(side, side_node);
      return true;
    case HEX27:
      node = hex27_side_nodes()(side, side_node);
      return true;
    case PRISM6:
      node = prism6_side_nodes()(side, side_node);
      return true;
    case PRISM15:
      node = prism15_side_nodes()(side, side_node);
      return true;
    case PRISM18:
      node = prism18_side_nodes()(side, side_node);
      return true;
    case PRISM20:
      node = prism20_side_nodes()(side, side_node);
      return true;
    case PRISM21:
      node = prism21_side_nodes()(side, side_node);
      return true;
    case PYRAMID5:
      node = pyramid5_side_nodes()(side, side_node);
      return true;
    case PYRAMID13:
      node = pyramid13_side_nodes()(side, side_node);
      return true;
    case PYRAMID14:
      node = pyramid14_side_nodes()(side, side_node);
      return true;
    case PYRAMID18:
      node = pyramid18_side_nodes()(side, side_node);
      return true;
    default:
      return false;
  }
}

LIBMESH_DEVICE_INLINE bool
try_local_edge_node(ElemType parent,
                    unsigned int edge,
                    unsigned int edge_node,
                    unsigned int & node)
{
  const unsigned int count = edge_node_count_or_zero(parent, edge);
  if (!count || edge_node >= count)
    return false;

  switch (parent)
  {
    case TET10:
      node = tet10_edge_nodes()(edge, edge_node);
      return true;
    case TET14:
      node = tet14_edge_nodes()(edge, edge_node);
      return true;
    case HEX20:
      node = hex20_edge_nodes()(edge, edge_node);
      return true;
    case HEX27:
      node = hex27_edge_nodes()(edge, edge_node);
      return true;
    case PRISM15:
      node = prism15_edge_nodes()(edge, edge_node);
      return true;
    case PRISM18:
      node = prism18_edge_nodes()(edge, edge_node);
      return true;
    case PRISM20:
      node = prism20_edge_nodes()(edge, edge_node);
      return true;
    case PRISM21:
      node = prism21_edge_nodes()(edge, edge_node);
      return true;
    case PYRAMID13:
      node = pyramid13_edge_nodes()(edge, edge_node);
      return true;
    case PYRAMID14:
      node = pyramid14_edge_nodes()(edge, edge_node);
      return true;
    case PYRAMID18:
      node = pyramid18_edge_nodes()(edge, edge_node);
      return true;
    default:
      return false;
  }
}

LIBMESH_DEVICE_INLINE bool
try_reference_node(ElemType type,
                   unsigned int node,
                   Point & pt)
{
  switch (type)
  {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      switch (node)
      {
        case 0:
          pt = Point(-1.0);
          return true;
        case 1:
          pt = Point(1.0);
          return true;
        case 2:
          if (type == EDGE3)
          {
            pt = Point(0.0);
            return true;
          }
          if (type == EDGE4)
          {
            pt = Point(-1. / 3.);
            return true;
          }
          return false;
        case 3:
          if (type == EDGE4)
          {
            pt = Point(1. / 3.);
            return true;
          }
          return false;
        default:
          return false;
      }

    case TRI3:
    case TRI6:
    case TRI7:
      switch (node)
      {
        case 0:
          pt = Point(0.0, 0.0);
          return true;
        case 1:
          pt = Point(1.0, 0.0);
          return true;
        case 2:
          pt = Point(0.0, 1.0);
          return true;
        case 3:
          pt = Point(0.5, 0.0);
          return true;
        case 4:
          pt = Point(0.5, 0.5);
          return true;
        case 5:
          pt = Point(0.0, 0.5);
          return true;
        case 6:
          if (type == TRI7)
          {
            pt = Point(1. / 3., 1. / 3.);
            return true;
          }
          return false;
        default:
          return false;
      }

    case QUAD4:
    case QUAD8:
    case QUAD9:
      switch (node)
      {
        case 0:
          pt = Point(-1.0, -1.0);
          return true;
        case 1:
          pt = Point(1.0, -1.0);
          return true;
        case 2:
          pt = Point(1.0, 1.0);
          return true;
        case 3:
          pt = Point(-1.0, 1.0);
          return true;
        case 4:
          pt = Point(0.0, -1.0);
          return true;
        case 5:
          pt = Point(1.0, 0.0);
          return true;
        case 6:
          pt = Point(0.0, 1.0);
          return true;
        case 7:
          pt = Point(-1.0, 0.0);
          return true;
        case 8:
          if (type == QUAD9)
          {
            pt = Point(0.0, 0.0);
            return true;
          }
          return false;
        default:
          return false;
      }

    case TET4:
    case TET10:
    case TET14:
      switch (node)
      {
        case 0:
          pt = Point(0.0, 0.0, 0.0);
          return true;
        case 1:
          pt = Point(1.0, 0.0, 0.0);
          return true;
        case 2:
          pt = Point(0.0, 1.0, 0.0);
          return true;
        case 3:
          pt = Point(0.0, 0.0, 1.0);
          return true;
        case 4:
          pt = Point(0.5, 0.0, 0.0);
          return true;
        case 5:
          pt = Point(0.5, 0.5, 0.0);
          return true;
        case 6:
          pt = Point(0.0, 0.5, 0.0);
          return true;
        case 7:
          pt = Point(0.0, 0.0, 0.5);
          return true;
        case 8:
          pt = Point(0.5, 0.0, 0.5);
          return true;
        case 9:
          pt = Point(0.0, 0.5, 0.5);
          return true;
        case 10:
          if (type == TET14)
          {
            pt = Point(1. / 3., 1. / 3., 0.0);
            return true;
          }
          return false;
        case 11:
          if (type == TET14)
          {
            pt = Point(1. / 3., 0.0, 1. / 3.);
            return true;
          }
          return false;
        case 12:
          if (type == TET14)
          {
            pt = Point(1. / 3., 1. / 3., 1. / 3.);
            return true;
          }
          return false;
        case 13:
          if (type == TET14)
          {
            pt = Point(0.0, 1. / 3., 1. / 3.);
            return true;
          }
          return false;
        default:
          return false;
      }

    case HEX8:
    case HEX20:
    case HEX27:
      switch (node)
      {
        case 0:
          pt = Point(-1.0, -1.0, -1.0);
          return true;
        case 1:
          pt = Point(1.0, -1.0, -1.0);
          return true;
        case 2:
          pt = Point(1.0, 1.0, -1.0);
          return true;
        case 3:
          pt = Point(-1.0, 1.0, -1.0);
          return true;
        case 4:
          pt = Point(-1.0, -1.0, 1.0);
          return true;
        case 5:
          pt = Point(1.0, -1.0, 1.0);
          return true;
        case 6:
          pt = Point(1.0, 1.0, 1.0);
          return true;
        case 7:
          pt = Point(-1.0, 1.0, 1.0);
          return true;
        case 8:
          pt = Point(0.0, -1.0, -1.0);
          return true;
        case 9:
          pt = Point(1.0, 0.0, -1.0);
          return true;
        case 10:
          pt = Point(0.0, 1.0, -1.0);
          return true;
        case 11:
          pt = Point(-1.0, 0.0, -1.0);
          return true;
        case 12:
          pt = Point(-1.0, -1.0, 0.0);
          return true;
        case 13:
          pt = Point(1.0, -1.0, 0.0);
          return true;
        case 14:
          pt = Point(1.0, 1.0, 0.0);
          return true;
        case 15:
          pt = Point(-1.0, 1.0, 0.0);
          return true;
        case 16:
          pt = Point(0.0, -1.0, 1.0);
          return true;
        case 17:
          pt = Point(1.0, 0.0, 1.0);
          return true;
        case 18:
          pt = Point(0.0, 1.0, 1.0);
          return true;
        case 19:
          pt = Point(-1.0, 0.0, 1.0);
          return true;
        case 20:
          if (type == HEX27)
          {
            pt = Point(0.0, 0.0, -1.0);
            return true;
          }
          return false;
        case 21:
          if (type == HEX27)
          {
            pt = Point(0.0, -1.0, 0.0);
            return true;
          }
          return false;
        case 22:
          if (type == HEX27)
          {
            pt = Point(1.0, 0.0, 0.0);
            return true;
          }
          return false;
        case 23:
          if (type == HEX27)
          {
            pt = Point(0.0, 1.0, 0.0);
            return true;
          }
          return false;
        case 24:
          if (type == HEX27)
          {
            pt = Point(-1.0, 0.0, 0.0);
            return true;
          }
          return false;
        case 25:
          if (type == HEX27)
          {
            pt = Point(0.0, 0.0, 1.0);
            return true;
          }
          return false;
        case 26:
          if (type == HEX27)
          {
            pt = Point(0.0, 0.0, 0.0);
            return true;
          }
          return false;
        default:
          return false;
      }

    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      switch (node)
      {
        case 0:
          pt = Point(-1.0, -1.0, 0.0);
          return true;
        case 1:
          pt = Point(1.0, -1.0, 0.0);
          return true;
        case 2:
          pt = Point(1.0, 1.0, 0.0);
          return true;
        case 3:
          pt = Point(-1.0, 1.0, 0.0);
          return true;
        case 4:
          pt = Point(0.0, 0.0, 1.0);
          return true;
        case 5:
          pt = Point(0.0, -1.0, 0.0);
          return true;
        case 6:
          pt = Point(1.0, 0.0, 0.0);
          return true;
        case 7:
          pt = Point(0.0, 1.0, 0.0);
          return true;
        case 8:
          pt = Point(-1.0, 0.0, 0.0);
          return true;
        case 9:
          if (type == PYRAMID13 || type == PYRAMID14 || type == PYRAMID18)
          {
            pt = Point(-0.5, -0.5, 0.5);
            return true;
          }
          return false;
        case 10:
          if (type == PYRAMID13 || type == PYRAMID14 || type == PYRAMID18)
          {
            pt = Point(0.5, -0.5, 0.5);
            return true;
          }
          return false;
        case 11:
          if (type == PYRAMID13 || type == PYRAMID14 || type == PYRAMID18)
          {
            pt = Point(0.5, 0.5, 0.5);
            return true;
          }
          return false;
        case 12:
          if (type == PYRAMID13 || type == PYRAMID14 || type == PYRAMID18)
          {
            pt = Point(-0.5, 0.5, 0.5);
            return true;
          }
          return false;
        case 13:
          if (type == PYRAMID14 || type == PYRAMID18)
          {
            pt = Point(0.0, 0.0, 0.0);
            return true;
          }
          return false;
        case 14:
          if (type == PYRAMID18)
          {
            pt = Point(-2. / 3., 0.0, 1. / 3.);
            return true;
          }
          return false;
        case 15:
          if (type == PYRAMID18)
          {
            pt = Point(0.0, 2. / 3., 1. / 3.);
            return true;
          }
          return false;
        case 16:
          if (type == PYRAMID18)
          {
            pt = Point(2. / 3., 0.0, 1. / 3.);
            return true;
          }
          return false;
        case 17:
          if (type == PYRAMID18)
          {
            pt = Point(0.0, -2. / 3., 1. / 3.);
            return true;
          }
          return false;
        default:
          return false;
      }

    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      switch (node)
      {
        case 0:
          pt = Point(0.0, 0.0, -1.0);
          return true;
        case 1:
          pt = Point(1.0, 0.0, -1.0);
          return true;
        case 2:
          pt = Point(0.0, 1.0, -1.0);
          return true;
        case 3:
          pt = Point(0.0, 0.0, 1.0);
          return true;
        case 4:
          pt = Point(1.0, 0.0, 1.0);
          return true;
        case 5:
          pt = Point(0.0, 1.0, 1.0);
          return true;
        case 6:
          if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.5, 0.0, -1.0);
            return true;
          }
          return false;
        case 7:
          if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.5, 0.5, -1.0);
            return true;
          }
          return false;
        case 8:
          if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.0, 0.5, -1.0);
            return true;
          }
          return false;
        case 9:
          if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.0, 0.0, 0.0);
            return true;
          }
          return false;
        case 10:
          if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(1.0, 0.0, 0.0);
            return true;
          }
          return false;
        case 11:
          if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.0, 1.0, 0.0);
            return true;
          }
          return false;
        case 12:
          if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.5, 0.0, 1.0);
            return true;
          }
          return false;
        case 13:
          if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.5, 0.5, 1.0);
            return true;
          }
          return false;
        case 14:
          if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.0, 0.5, 1.0);
            return true;
          }
          return false;
        case 15:
          if (type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.5, 0.0, 0.0);
            return true;
          }
          return false;
        case 16:
          if (type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.5, 0.5, 0.0);
            return true;
          }
          return false;
        case 17:
          if (type == PRISM18 || type == PRISM20 || type == PRISM21)
          {
            pt = Point(0.0, 0.5, 0.0);
            return true;
          }
          return false;
        case 18:
          if (type == PRISM20 || type == PRISM21)
          {
            pt = Point(1. / 3., 1. / 3., -1.0);
            return true;
          }
          return false;
        case 19:
          if (type == PRISM20 || type == PRISM21)
          {
            pt = Point(1. / 3., 1. / 3., 1.0);
            return true;
          }
          return false;
        case 20:
          if (type == PRISM21)
          {
            pt = Point(1. / 3., 1. / 3., 0.0);
            return true;
          }
          return false;
        default:
          return false;
      }

    default:
      return false;
  }
}

LIBMESH_DEVICE_INLINE bool
try_refspace_node(ElemType type,
                  unsigned int node,
                  Point & pt)
{
  switch (type)
  {
    case NODEELEM:
      if (!node)
      {
        pt = Point(0.0, 0.0, 0.0);
        return true;
      }
      return false;

    case TRISHELL3:
      return try_reference_node(TRI3, node, pt);

    case QUADSHELL4:
      return try_reference_node(QUAD4, node, pt);

    case QUADSHELL8:
      return try_reference_node(QUAD8, node, pt);

    case QUADSHELL9:
      return try_reference_node(QUAD9, node, pt);

    default:
      return try_reference_node(type, node, pt);
  }
}

LIBMESH_DEVICE_INLINE bool
try_reference_side_node(ElemType parent,
                        unsigned int side,
                        unsigned int side_node,
                        Point & pt)
{
  unsigned int node = libMesh::invalid_uint;
  if (!try_local_side_node(parent, side, side_node, node))
    return false;

  return try_reference_node(parent, node, pt);
}

} // namespace libMesh

#endif // LIBMESH_FE_REFERENCE_ELEMENT_TRAITS_H
