/*
 * Example illustrates the use of SZIP compression in HDF5
 */

#include <stdio.h>
#include "hdf5.h"

#define NX 500
#define NY 600
#define CH_NX 100
#define CH_NY 25

static void initialize(void);
static int compare(void);

static float buf[NX][NY];
static float buf_r[NX][NY];

static const char* filename = "test.h5";

static int
writeszip()
{
   hid_t file;
   hid_t dataset32;
   hid_t properties;
   hid_t lcpl_id, dapl_id;
   hid_t data_space;
   hsize_t dims[2], chunk_size[2];
   unsigned szip_options_mask;
   unsigned szip_pixels_per_block;

  /*
   * Create a new file using read/write access, default file 
   * creation properties, and default file access properties.
   */
   file = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

   /* Describe the size of the array. */
   dims[0] = NX;
   dims[1] = NY;
   data_space = H5Screate_simple (2, dims, NULL);

  /*
   * Set the dataset creation property list to specify that
   * the raw data is to be partitioned into 100x100 element
   * chunks and that each chunk is to be compressed.
   */
   chunk_size[0] = CH_NX;
   chunk_size[1] = CH_NY;
   properties = H5Pcreate (H5P_DATASET_CREATE);
   H5Pset_chunk (properties, 2, chunk_size);

  /* 
   * Set parameters for SZIP compression; check the description of
   * the H5Pset_szip function in the HDF5 Reference Manual for more 
   * information.
   */
   szip_options_mask=H5_SZIP_NN_OPTION_MASK;
   szip_pixels_per_block=32;

   H5Pset_szip (properties, szip_options_mask, szip_pixels_per_block);

  /*
   * Create a new dataset within the file.  The datatype
   * and data space describe the data on disk, which may
   * be different from the format used in the application's
   * memory.
   */

   lcpl_id = H5Pcreate (H5P_LINK_CREATE);
   dapl_id = H5Pcreate (H5P_DATASET_ACCESS);

   dataset32 = H5Dcreate (file, "datasetF32", H5T_NATIVE_FLOAT, data_space, lcpl_id, properties, dapl_id);

  /*
   * Write the array to the file.  The datatype and dataspace
   * describe the format of the data in the `buf' buffer.
   * The raw data is translated to the format required on disk, 
   * as defined above.  We use default raw data transfer properties.
   */

   H5Dwrite (dataset32, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, buf);

   H5Sclose (data_space);
   H5Pclose(lcpl_id);
   H5Pclose(dapl_id);
   H5Pclose (properties);
   H5Dclose (dataset32);
   H5Fclose (file);

   return 1;
}

static int
readszip()
{
    hid_t file;
    hid_t dataset32;
    hid_t properties;
    int errcnt = 0;

    file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    properties = H5Pcreate(H5P_DATASET_ACCESS);
    dataset32 = H5Dopen(file, "datasetF32", properties);

    /*
     * Read the array.  This is similar to writing data,
     * except the data flows in the opposite direction.
     * Note: Decompression is automatic.
     */

    H5Dread(dataset32, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf_r);

    errcnt = compare();

    H5Pclose(properties);
    H5Dclose (dataset32);
    H5Fclose (file);

    return (errcnt==0 ? 1 : 0);
}

static int
compare(void)
{
    int i,j;
    int errs = 0;
    
    /* Do comparison */
    for (i=0; i < NX; i++) {
	for (j=0; j < NY; j++) {
            if(buf[i][j] != buf_r[i][j]) {
		errs++;
	        printf("mismatch: [%d][%d]: write = %f read=%f\n",
		        i,j,buf[i][j],buf_r[i][j]);
	}
     }
   }
   return errs;
}

static void
initialize(void)
{
   int i, j;
   /* Initialize data buffer with some bogus data. */
   for(i=0; i < NX; i++) {
     for(j=0; j < NY; j++) {
       buf[i][j] = (float)(i + j);
     }
   }
}

int
main(int argc, char** argv)
{
    int extfile = 0;
    if(argc > 1) {
	filename = argv[1];
	extfile = 1;
    }

    initialize();
    if(!extfile) {
	if(!writeszip()) {
	    fprintf(stderr,"writeszip failed.\n");
	    goto fail;
	}
    }

    if(!readszip()) {
	fprintf(stderr,"openfile failed.\n");
	goto fail;
    }
    fprintf(stderr,"***PASS\n");
    return 0;
fail:
    fprintf(stderr,"***FAIL\n");
    return 1;
}
