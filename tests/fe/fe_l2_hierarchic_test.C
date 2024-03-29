
#include "fe_test.h"

INSTANTIATE_FETEST(FIRST, L2_HIERARCHIC, EDGE2);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, EDGE2);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, EDGE3);
INSTANTIATE_FETEST(THIRD, L2_HIERARCHIC, EDGE3);
INSTANTIATE_FETEST(FOURTH, L2_HIERARCHIC, EDGE3);

#if LIBMESH_DIM > 1
INSTANTIATE_FETEST(FIRST, L2_HIERARCHIC, TRI3);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, TRI3);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, TRI6);
INSTANTIATE_FETEST(THIRD, L2_HIERARCHIC, TRI6);
INSTANTIATE_FETEST(FOURTH, L2_HIERARCHIC, TRI6);

INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, TRI7);
INSTANTIATE_FETEST(THIRD, L2_HIERARCHIC, TRI7);
INSTANTIATE_FETEST(FOURTH, L2_HIERARCHIC, TRI7);

INSTANTIATE_FETEST(FIRST, L2_HIERARCHIC, QUAD4);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, QUAD4);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, QUAD9);
INSTANTIATE_FETEST(THIRD, L2_HIERARCHIC, QUAD9);
INSTANTIATE_FETEST(FOURTH, L2_HIERARCHIC, QUAD9);
#endif

#if LIBMESH_DIM > 2
INSTANTIATE_FETEST(FIRST, L2_HIERARCHIC, TET10);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, TET4);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, TET14);
INSTANTIATE_FETEST(THIRD, L2_HIERARCHIC, TET4);
INSTANTIATE_FETEST(FOURTH, L2_HIERARCHIC, TET10);

INSTANTIATE_FETEST(FIRST, L2_HIERARCHIC, PRISM6);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, PRISM6);
INSTANTIATE_FETEST(FOURTH, L2_HIERARCHIC, PRISM6);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, PRISM15);
INSTANTIATE_FETEST(FIRST, L2_HIERARCHIC, PRISM18);
INSTANTIATE_FETEST(THIRD, L2_HIERARCHIC, PRISM20);
INSTANTIATE_FETEST(FOURTH, L2_HIERARCHIC, PRISM20);
INSTANTIATE_FETEST(THIRD, L2_HIERARCHIC, PRISM21);
INSTANTIATE_FETEST(FOURTH, L2_HIERARCHIC, PRISM21);

INSTANTIATE_FETEST(FIRST, L2_HIERARCHIC, HEX8);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, HEX8);
INSTANTIATE_FETEST(SECOND, L2_HIERARCHIC, HEX27);
INSTANTIATE_FETEST(THIRD, L2_HIERARCHIC, HEX27);
INSTANTIATE_FETEST(FOURTH, L2_HIERARCHIC, HEX27);
#endif
