#include "config.h"
#include <nc_tests.h>
#include "err_macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define FILECLASSIC "tst_dimsize_classic.nc"
#define FILE64OFFSET "tst_dimsize_64offset.nc"
#define FILE64DATA "tst_dimsize_64data.nc"

#define DIMMAXCLASSIC (NC_MAX_INT - 3)
#define DIMMAX64OFFSET (NC_MAX_UINT - 3)

#ifdef ENABLE_CDF5
#define DIMMAX64DATA (NC_MAX_UINT64 - 3)
#endif

/*
Test that at least the meta-data works
for dimension sizes X modes.
NC_CLASSIC => NC_INT_MAX - 3
NC_64BIT_OFFSET => NC_UINT_MAX - 3
NC_64BIT_DATA => NC_UINT64_MAX - 3
Note that this will not test the last case when
|size_t| == 4.
Also, leave the files around so we can test with ncdump.
*/

int
main(int argc, char **argv)
{
    int ncid;
    size_t dimsize;
    int dimid;
    int stat = NC_NOERR;
    printf("\n*** Testing Max Dimension Sizes\n");

    printf("\n|size_t|=%lu\n",(unsigned long)sizeof(size_t));

    printf("\n*** Writing Max Dimension Size (%d) For NC_CLASSIC\n",DIMMAXCLASSIC);
    if ((stat=nc_create(FILECLASSIC, NC_CLOBBER, &ncid))) ERRSTAT(stat);
    dimsize = DIMMAXCLASSIC;
    if ((stat=nc_def_dim(ncid, "testdim", dimsize, &dimid))) ERRSTAT(stat);
    if ((stat=nc_close(ncid))) ERRSTAT(stat);

    printf("\n*** Reading Max Dimension Size For NC_CLASSIC\n");
    if ((stat=nc_open(FILECLASSIC, NC_NOCLOBBER, &ncid))) ERRSTAT(stat);
    if ((stat=nc_inq_dimid(ncid, "testdim", &dimid))) ERRSTAT(stat);
    if ((stat=nc_inq_dimlen(ncid, dimid, &dimsize))) ERRSTAT(stat);
    if(dimsize != DIMMAXCLASSIC) ERR;
    if ((stat=nc_close(ncid))) ERRSTAT(stat);

    printf("\n*** Writing Max Dimension Size (%u) For NC_64BIT_OFFSET\n",DIMMAX64OFFSET);
    if ((stat=nc_create(FILE64OFFSET, NC_CLOBBER | NC_64BIT_OFFSET, &ncid))) ERRSTAT(stat);
    dimsize = DIMMAX64OFFSET;
    if ((stat=nc_def_dim(ncid, "testdim", dimsize, &dimid))) ERRSTAT(stat);
    if ((stat=nc_enddef(ncid))) ERRSTAT(stat);
    if ((stat=nc_close(ncid))) ERRSTAT(stat);

    printf("\n*** Reading Max Dimension Size For NC_64BIT_OFFSET\n");
    if ((stat=nc_open(FILE64OFFSET, NC_NOCLOBBER|NC_64BIT_OFFSET, &ncid))) ERRSTAT(stat);
    if ((stat=nc_inq_dimid(ncid, "testdim", &dimid))) ERRSTAT(stat);
    if ((stat=nc_inq_dimlen(ncid, dimid, &dimsize))) ERRSTAT(stat);
    if(dimsize != DIMMAX64OFFSET) ERR;
    if ((stat=nc_close(ncid))) ERRSTAT(stat);

#ifdef ENABLE_CDF5
    if(sizeof(size_t) == 8) {
      printf("\n*** Writing Max Dimension Size (%llu) For NC_64BIT_DATA\n",DIMMAX64DATA);
        if ((stat=nc_create(FILE64DATA, NC_CLOBBER | NC_64BIT_DATA, &ncid))) ERRSTAT(stat);
        dimsize = (size_t)DIMMAX64DATA;
        if ((stat=nc_def_dim(ncid, "testdim", dimsize, &dimid))) ERRSTAT(stat);
        if ((stat=nc_close(ncid))) ERRSTAT(stat);

	printf("\n*** Reading Max Dimension Size For NC_64BIT_DATA\n");
	if ((stat=nc_open(FILE64DATA, NC_NOCLOBBER|NC_64BIT_DATA, &ncid))) ERRSTAT(stat);
	if ((stat=nc_inq_dimid(ncid, "testdim", &dimid))) ERRSTAT(stat);
	if ((stat=nc_inq_dimlen(ncid, dimid, &dimsize))) ERRSTAT(stat);
	if(dimsize != DIMMAX64DATA) ERR;
	if ((stat=nc_close(ncid))) ERRSTAT(stat);
    }
#endif /* ENABLE_CDF5 */

    SUMMARIZE_ERR;
    FINAL_RESULTS;
}
