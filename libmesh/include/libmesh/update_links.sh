#!/bin/sh

set -x

ln -s ../base/*.h \
      ../enums/*.h \
      ../error_estimation/*.h \
      ../fe/*.h \
      ../geom/*.h \
      ../mesh/*.h \
      ../numerics/*.h \
      ../parallel/*.h \
      ../partitioning/*.h \
      ../physics/*.h \
      ../quadrature/*.h \
      ../reduced_basis/*.h \
      ../solvers/*.h \
      ../systems/*.h \
      ../utils/*.h \
      . 2>/dev/null