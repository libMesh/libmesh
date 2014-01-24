#!/bin/sh
#set -x

# This shell just tests the group renaming.

set -e
echo ""

# Create the input cdl file
rm -f tst_grp_rename.cdl
cat >tst_grp_rename.cdl <<EOF
netcdf tst_grp_rename {

types:
  int(*) vlen_t;

dimensions:
  d2 = 2;

variables:
  vlen_t v1(d2);

group: inner {

  types:
    compound c_t { int f1; float f2; };

  dimensions:
    d3 = 3;

  variables:
    c_t vc(d3);

  group: inner_inner {
    dimensions:
      d3 = 4;
  }
}
}
EOF

# Create the reference cdl file
rm -f ref_grp_rename.cdl
cat >ref_grp_rename.cdl <<EOF
netcdf tst_grp_rename {
types:
  int(*) vlen_t ;
dimensions:
	d2 = 2 ;
variables:
	vlen_t v1(d2) ;
data:

 v1 = {}, {} ;

group: renamed {
  types:
    compound c_t {
      int f1 ;
      float f2 ;
    }; // c_t
  dimensions:
  	d3 = 3 ;
  variables:
  	c_t vc(d3) ;
  data:

   vc = {0, 0}, {0, 0}, {0, 0} ;

  group: inner_renamed {
    dimensions:
    	d3 = 4 ;
    } // group inner_renamed
  } // group renamed
}
EOF

echo "*** Running group_rename test"

FAIL=0

# Create ref_tst_group_rename.nc
rm -f tst_grp_rename.nc
../ncgen/ncgen -k3 ./tst_grp_rename.cdl

# Try to rename 2nd level group
if ! ./renamegroup tst_grp_rename.nc "/inner/inner_inner" "inner_renamed" ; then
  echo "***FAIL: attempt to rename /inner/inner_inner failed"
  FAIL=1
fi

# Try to 1st level group
if ! ./renamegroup tst_grp_rename.nc "/inner" "renamed" ; then
  echo "***FAIL: attempt to rename /inner failed"
  FAIL=1
fi

# Dump the final .nc and compare with the reference
rm -f tst_grp_rename.dmp
../ncdump/ncdump tst_grp_rename.nc > ./tst_grp_rename.dmp

if ! diff ref_grp_rename.cdl tst_grp_rename.dmp ; then
  echo "***FAIL: output and reference output differ"
  FAIL=1
fi

# Finally, try to rename root group; should fail
if ./renamegroup tst_grp_rename.nc "/" "rootgroup" ; then
  echo "***FAIL: attempt to rename root group should not have succeeded"
  FAIL=1
else
  echo "***XFAIL : attempt to rename root group failed as expected"
fi

rm -f tst_grp_rename.cdl tst_grp_rename.nc ref_grp_rename.nc

exit $FAIL
