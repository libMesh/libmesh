/*! Testing for proper read of little-endian variables in an hdf4 file.
 *
 * Added to debug issue NCF-332. Based on code submitted by
 * https://github.com/Unidata/netcdf-c/issues/113.
 */

#include <stdio.h>
#include <config.h>
#include <unistd.h>
#include <nc_tests.h>
#include "err_macros.h"
#include <hdf5.h>
#include <H5DSpublic.h>
#include "mfhdf.h"

#define DIM1 5
#define DIM0 5
#define RANK 2
#define FILENAME "tst_h4_lendian.h4"
#define SDSNAME "data"

int read_hdf_file(int dtype) {

  int ncid = 0;
  int le_int16_varid = 0;
  int retval = 0;
  int ed = 0;

  printf("\to Reading hdf4 file with a little-endian datatype %d\n",dtype);

  printf("\t\to Opening file....\t\t\t\t\t");
  retval = nc_open(FILENAME, NC_NETCDF4 | NC_NOWRITE, &ncid);
  if(retval) {printf("Failure [%d]\n",retval); return retval;}
  else {printf("Success\n");}

  printf("\t\to Getting varid....\t\t\t\t\t");
  retval = nc_inq_varid(ncid,SDSNAME,&le_int16_varid);
  if(retval) {printf("Failure [%d]\n",retval); return retval;}
  else {printf("Success\n");}

  printf("\t\to Querying endianness of the variable....\t\t");
  retval = nc_inq_var_endian(ncid,le_int16_varid,&ed);
  if(retval) {printf("Failure [%d]\n",retval); return retval;}
  else {printf("Success\n");}

  printf("\t\to Checking that endianness is NC_ENDIAN_LITTLE....\t");
  if (ed == NC_ENDIAN_LITTLE) printf("Success\n");
  else {printf("Failure [%d]\n\n",ed); nc_close(ncid); return -1;}

  printf("\t\to Closing file....\t\t\t\t\t");
  retval = nc_close(ncid);
  if(retval) {printf("Failure [%d]\n\n",retval); return retval;}
  else {printf("Success\n\n");}

  return 0;
}

int create_hdf_file(int dtype) {

    int32 sd_id, sds_id, istat, sd_index;
    int32 dims[2], start[2], edges[2], rank;
    int16 array_data[DIM0][DIM1];
    intn i, j, count;

    start[0] = 0;
    start[1] = 0;
    edges[0] = DIM1;
    edges[1] = DIM0;

    /* populate data array */
    count = 0;
    for (j = 0; j < DIM0; j++)
      {
        for (i = 0; i < DIM1; i++)
          array_data[j][i] = count++;
      }

    printf("\to Creating hdf4 file with little-endian datatype %d....\t",dtype);

    sd_id = SDstart(FILENAME, DFACC_CREATE);
    /* sds_id = SDcreate(sd_id, SDSNAME, DFNT_LITEND|dtype, RANK, edges); */
    sds_id = SDcreate(sd_id, SDSNAME, dtype, RANK, edges);

    istat = SDendaccess(sds_id);
    if(istat) {printf("Failure %d\n", istat); SDend(sd_id); return istat;}

    istat = SDend(sd_id);
    if(istat) {printf("Failure %d\n", istat); SDend(sd_id); return istat;}

    sd_id = SDstart(FILENAME, DFACC_WRITE);

    sd_index = 0;
    sds_id = SDselect(sd_id, sd_index);

    istat = SDwritedata(sds_id, start, NULL, edges, (VOIDP)array_data);
    if(istat) {printf("Failure %d\n", istat); SDend(sd_id); return istat;}

    istat = SDendaccess(sds_id);
    if(istat) {printf("Failure %d\n", istat); SDend(sd_id); return istat;}

    istat = SDend(sd_id);
    if(istat) {printf("Failure %d\n", istat); return istat;}

    printf("Success\n");
    return 0;
}


int test_read_write(int dtype) {

  int res = 0;

  res = create_hdf_file(dtype);
  if(res) {unlink(FILENAME); return res;}

  res = read_hdf_file(dtype);

  unlink(FILENAME);
  return res;
}

/*! Standard main function.
 *
 */
int main() {

  int res = 0;

  printf("Test reading from an hdf4 file with a little-endian datatype.\n");

  /* True Positives. */
  res = test_read_write(DFNT_LINT8);
  res = test_read_write(DFNT_LUINT8);
  res = test_read_write(DFNT_LINT16);
  res = test_read_write(DFNT_LUINT16);
  res = test_read_write(DFNT_LINT32);
  res = test_read_write(DFNT_LUINT32);
  res = test_read_write(DFNT_LFLOAT32);
  res = test_read_write(DFNT_LFLOAT64);

  /* True Negatives. */
  printf("\t**** Testing for True Negatives. THESE SHOULD FAIL.****\n\n");
  res = test_read_write(DFNT_INT8);
  if(!res) {printf("Should have failed. Error!\n"); return -1;}

  res = test_read_write(DFNT_UINT8);
  if(!res) {printf("Should have failed. Error!\n"); return -1;}

  res = test_read_write(DFNT_INT16);
  if(!res) {printf("Should have failed. Error!\n"); return -1;}

  res = test_read_write(DFNT_UINT16);
  if(!res) {printf("Should have failed. Error!\n"); return -1;}

  res = test_read_write(DFNT_INT32);
  if(!res) {printf("Should have failed. Error!\n"); return -1;}

  res = test_read_write(DFNT_UINT32);
  if(!res) {printf("Should have failed. Error!\n"); return -1;}

  res = test_read_write(DFNT_FLOAT32);
  if(!res) {printf("Should have failed. Error!\n"); return -1;}

  res = test_read_write(DFNT_FLOAT64);
  if(!res) {printf("Should have failed. Error!\n"); return -1;}

  printf("Finished.\n");
  return 0;
}
