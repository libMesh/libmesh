/* This is part of the netCDF package. Copyright 2006-2018 University
   Corporation for Atmospheric Research/Unidata. See COPYRIGHT file
   for conditions of use.

   Test the netCDF-4 attribute code.

   Ed Hartnett
*/

#include <config.h>
#include <nc_tests.h>
#include "err_macros.h"
#include "hdf5internal.h"

/* The data file we will create. */
#define FILE_NAME "tst_atts.nc"

/* Names of attributes. */
#define OLD_NAME "Constantinople"
#define OLD_NAME_2 "Constantinopolis"
#define NEW_NAME "Istanbul"

/* Contents of attributes. */
#define CONTENTS "Lots of people!"
#define CONTENTS_2 "Lots of people!!" /* 1 longer than CONTENTS */
#define CONTENTS_3 "Lots 0f pe0ple!"  /* same len as CONTENTS */
#define VAR_NAME "Earth"

/**
WARNING: following should match lists in libsrc4/nc4file.c
*/

/**
 * @internal Define the names of attributes to ignore added by the
 * HDF5 dimension scale; these attached to variables. They cannot be
 * modified thru the netcdf-4 API.
 */
static const char* NC_RESERVED_VARATT_LIST[] = {
   NC_ATT_REFERENCE_LIST,
   NC_ATT_CLASS,
   NC_ATT_DIMENSION_LIST,
   NC_ATT_NAME,
   NC_ATT_COORDINATES,
   NC_DIMID_ATT_NAME,
   NULL
};

/**
 * @internal Define the names of attributes to ignore because they are
 * "hidden" global attributes. They can be read, but not modified thru
 * the netcdf-4 API.
 */
static const char* NC_RESERVED_ATT_LIST[] = {
   NC_ATT_FORMAT,
   NC3_STRICT_ATT_NAME,
   NCPROPS,
   ISNETCDF4ATT,
   SUPERBLOCKATT,
   NULL
};

/**
 * @internal Define the subset of the reserved list that is readable
 * by name only
*/
static const char* NC_RESERVED_SPECIAL_LIST[] = {
   ISNETCDF4ATT,
   SUPERBLOCKATT,
   NCPROPS,
   NULL
};

int
main(int argc, char **argv)
{
   printf("\n*** Testing netCDF-4 attributes.\n");
   printf("*** testing attribute renaming for read-only file...");
   {
      int ncid;

      /* Create a file with an att. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLOBBER, &ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS),
                          CONTENTS)) ERR;
      if (nc_close(ncid)) ERR;

      /* Reopen the file read-only. */
      if (nc_open(FILE_NAME, NC_NOWRITE, &ncid)) ERR;

      /* Try to rename the att, but it won't work. */
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, NEW_NAME) != NC_EPERM) ERR;

      /* Try to create another att, it also won't work. */
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME_2, strlen(CONTENTS),
                          CONTENTS) != NC_EPERM) ERR;

      /* Try to delete the att. More failure ensues. */
      if (nc_del_att(ncid, NC_GLOBAL, OLD_NAME) != NC_EPERM) ERR;

      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("*** testing deleting atts...");
   {
      int ncid;
      int natts;

      /* Create a file with two atts. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLOBBER, &ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS),
                          CONTENTS)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME_2, 0, NULL)) ERR;

      /* These will not work. */
      if (nc_del_att(ncid + TEST_VAL_42, NC_GLOBAL, OLD_NAME) != NC_EBADID) ERR;
      if (nc_del_att(ncid, TEST_VAL_42, OLD_NAME) != NC_ENOTVAR) ERR;
      if (nc_del_att(ncid, NC_GLOBAL, NULL) != NC_EINVAL) ERR;
      if (nc_del_att(ncid, NC_GLOBAL, NEW_NAME) != NC_ENOTATT) ERR;

      /* End define mode. It redef will be called automatically. */
      if (nc_enddef(ncid)) ERR;

      /* Delete the attribute. */
      if (nc_del_att(ncid, NC_GLOBAL, OLD_NAME)) ERR;
      if (nc_close(ncid)) ERR;

      /* Reopen the file. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_inq_natts(ncid, &natts)) ERR;
      if (natts != 1) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("*** testing deleting atts classic model...");
   {
      int ncid;
      int natts;

      /* Create a file with an att. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLASSIC_MODEL, &ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS),
                          CONTENTS)) ERR;

      /* End define mode. */
      if (nc_enddef(ncid)) ERR;

      /* This will not work. */
      if (nc_del_att(ncid, NC_GLOBAL, OLD_NAME) != NC_ENOTINDEFINE) ERR;

      /* Delete the attribute. Redef is needed since this is a classic
       * model file. */
      if (nc_redef(ncid)) ERR;
      if (nc_del_att(ncid, NC_GLOBAL, OLD_NAME)) ERR;
      if (nc_close(ncid)) ERR;

      /* Reopen the file. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_inq_natts(ncid, &natts)) ERR;
      if (natts) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("*** testing over-writing atts classic model...");
   {
      int ncid;
      int natts;
      char *data_in;

      if (!(data_in = malloc(strlen(CONTENTS) + 1))) ERR;

      /* Create a file with an att. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLASSIC_MODEL, &ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS),
                          CONTENTS)) ERR;

      /* End define mode. */
      if (nc_enddef(ncid)) ERR;

      /* Try and write a new att. Won't work. */
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME_2, strlen(CONTENTS_2),
                          CONTENTS_2) != NC_ENOTINDEFINE) ERR;

      /* This will not work. Overwriting att must be same length or
       * shorter if not in define mode. */
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS_2),
                          CONTENTS_2) != NC_ENOTINDEFINE) ERR;

      /* Now overwrite the att. */
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS_3),
                          CONTENTS_3)) ERR;

      /* Delete the attribute. Redef is needed since this is a classic
       * model file. This should work but does not. */
      if (nc_redef(ncid)) ERR;
      if (nc_del_att(ncid, NC_GLOBAL, OLD_NAME)) ERR;
      if (nc_close(ncid)) ERR;

      /* Reopen the file. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_inq_natts(ncid, &natts)) ERR;
      /* If delete worked, natts would be 0. */
      /* if (natts != 0) ERR; */
      if (natts != 1) ERR;

      /* Get the attribute. */
      if (nc_get_att_text(ncid, NC_GLOBAL, OLD_NAME, data_in)) ERR;
      /* if (strncmp(CONTENTS_3, data_in, strlen(CONTENTS))) ERR; */
      free(data_in);
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("*** testing attribute renaming for a global attribute...");
   {
      int ncid, attid;
      char *data_in;
      char too_long_name[NC_MAX_NAME + 2];
      char name_in[NC_MAX_NAME + 1];

      /* Set up a name that is too long for netCDF. */
      memset(too_long_name, 'a', NC_MAX_NAME + 1);
      too_long_name[NC_MAX_NAME + 1] = 0;

      if (!(data_in = malloc(strlen(CONTENTS) + 1))) ERR;

      /* Create a file with an att. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLOBBER, &ncid)) ERR;
      if (nc_def_var(ncid, VAR_NAME, NC_INT, 0, NULL, NULL)) ERR;

      /* These will not work. */
      if (nc_put_att_text(ncid + TEST_VAL_42, NC_GLOBAL, OLD_NAME,
                          strlen(CONTENTS), CONTENTS) != NC_EBADID) ERR;
      if (nc_put_att_text(ncid, TEST_VAL_42, OLD_NAME, strlen(CONTENTS),
                          CONTENTS) != NC_ENOTVAR) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, NULL, strlen(CONTENTS),
                          CONTENTS) != NC_EBADNAME) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, BAD_NAME, strlen(CONTENTS),
                          CONTENTS) != NC_EBADNAME) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, too_long_name, strlen(CONTENTS),
                          CONTENTS) != NC_EBADNAME) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS),
                          NULL) != NC_EINVAL) ERR;
      {
         /* Check that the NC_GLOBAL reserved words are rejected. */
         const char** reserved = NC_RESERVED_ATT_LIST;
         for ( ; *reserved; reserved++)
         {
            if (nc_put_att_text(ncid, NC_GLOBAL, *reserved, strlen(CONTENTS),
                                CONTENTS) != NC_ENAMEINUSE) ERR;
         }
      }
      {
         /* Check that the variable reserved words are rejected. */
         const char** reserved = NC_RESERVED_VARATT_LIST;
         for ( ; *reserved; reserved++)
         {
           if (nc_put_att_text(ncid, 0, *reserved, strlen(CONTENTS),
                                CONTENTS) != NC_ENAMEINUSE) ERR;
         }
      }
      {
         /* Check that the read-only reserved words are rejected. */
         const char** reserved = NC_RESERVED_SPECIAL_LIST;
         for ( ; *reserved; reserved++)
         {
            if (nc_put_att_text(ncid, NC_GLOBAL, *reserved, strlen(CONTENTS),
                                CONTENTS) != NC_ENAMEINUSE) ERR;
         }
      }

      /* Write the attribute at last. */
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS),
                          CONTENTS)) ERR;

      /* Write another with different name. */
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME_2, strlen(CONTENTS),
                          CONTENTS)) ERR;

      /* These will not work. */
      if (nc_rename_att(ncid + TEST_VAL_42, NC_GLOBAL, OLD_NAME, NEW_NAME) != NC_EBADID) ERR;
      if (nc_rename_att(ncid, TEST_VAL_42, OLD_NAME, NEW_NAME) != NC_ENOTVAR) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, NULL) != NC_EINVAL) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, NULL, NEW_NAME) != NC_EINVAL) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, NULL, NULL) != NC_EINVAL) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, BAD_NAME) != NC_EBADNAME) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, too_long_name) != NC_EMAXNAME) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, OLD_NAME_2) != NC_ENAMEINUSE) ERR;

      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS),
                          NULL) != NC_EINVAL) ERR;
      {
         /* Check that the NC_GLOBAL reserved words are rejected. */
         const char** reserved = NC_RESERVED_ATT_LIST;
         for ( ; *reserved; reserved++)
         {
            if (nc_put_att_text(ncid, NC_GLOBAL, *reserved, strlen(CONTENTS),
                                CONTENTS) != NC_ENAMEINUSE) ERR;
         }
      }
      {
         /* Check that the variable reserved words are rejected. */
         const char** reserved = NC_RESERVED_VARATT_LIST;
         for ( ; *reserved; reserved++)
         {
            if (nc_put_att_text(ncid, 0, *reserved, strlen(CONTENTS),
                                CONTENTS) != NC_ENAMEINUSE) ERR;
         }
      }
      {
         /* Check that the read-only reserved words are rejected. */
         const char** reserved = NC_RESERVED_SPECIAL_LIST;
         for ( ; *reserved; reserved++)
         {
            if (nc_put_att_text(ncid, NC_GLOBAL, *reserved, strlen(CONTENTS),
                                CONTENTS) != NC_ENAMEINUSE) ERR;
         }
      }
      
      /* Write the attribute at last. */
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME, strlen(CONTENTS),
                          CONTENTS)) ERR;
      
      /* Write another with different name. */
      if (nc_put_att_text(ncid, NC_GLOBAL, OLD_NAME_2, strlen(CONTENTS),
                          CONTENTS)) ERR;

      /* These will not work. */
      if (nc_rename_att(ncid + TEST_VAL_42, NC_GLOBAL, OLD_NAME, NEW_NAME) != NC_EBADID) ERR;
      if (nc_rename_att(ncid, TEST_VAL_42, OLD_NAME, NEW_NAME) != NC_ENOTVAR) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, NULL) != NC_EINVAL) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, NULL, NEW_NAME) != NC_EINVAL) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, NULL, NULL) != NC_EINVAL) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, BAD_NAME) != NC_EBADNAME) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, too_long_name) != NC_EMAXNAME) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, OLD_NAME_2) != NC_ENAMEINUSE) ERR;
      
      /* Rename the att. */
      if (nc_rename_att(ncid, NC_GLOBAL, OLD_NAME, NEW_NAME)) ERR;

      /* These will not work. */
      if (nc_inq_attid(ncid + TEST_VAL_42, NC_GLOBAL, NEW_NAME, &attid) != NC_EBADID) ERR;
      if (nc_inq_attid(ncid, TEST_VAL_42, NEW_NAME, &attid) != NC_ENOTVAR) ERR;
      if (nc_inq_attid(ncid, NC_GLOBAL, NULL, &attid) != NC_EBADNAME) ERR;

      /* Check the file. */
      if (nc_inq_attid(ncid, NC_GLOBAL, NEW_NAME, &attid)) ERR;
      if (attid != 0) ERR;

      /* This also works. */
      if (nc_inq_attid(ncid, NC_GLOBAL, NEW_NAME, NULL)) ERR;

      /* These won't work. */
      if (nc_inq_attname(ncid + TEST_VAL_42, NC_GLOBAL, attid, name_in) != NC_EBADID) ERR;
      if (nc_inq_attname(ncid, TEST_VAL_42, attid, name_in) != NC_ENOTVAR) ERR;
      if (nc_inq_attname(ncid, NC_GLOBAL, -1, name_in) != NC_ENOTATT) ERR;

      /* Get the name from the ID. */
      if (nc_inq_attname(ncid, NC_GLOBAL, attid, name_in)) ERR;
      if (strcmp(name_in, NEW_NAME)) ERR;

      /* Also works but does little. */
      if (nc_inq_attname(ncid, NC_GLOBAL, attid, NULL)) ERR;

      /* These will not work. */
      if (nc_get_att_text(ncid + TEST_VAL_42, NC_GLOBAL, NEW_NAME, data_in) != NC_EBADID) ERR;
      if (nc_get_att_text(ncid, TEST_VAL_42, NEW_NAME, data_in) != NC_ENOTVAR) ERR;
      if (nc_get_att_text(ncid, NC_GLOBAL, NULL, data_in) != NC_EBADNAME) ERR;

      /* Get the attribute at last. */
      if (nc_get_att_text(ncid, NC_GLOBAL, NEW_NAME, data_in)) ERR;
      if (strncmp(CONTENTS, data_in, strlen(CONTENTS))) ERR;

      /* This also works. */
      if (nc_get_att_text(ncid, NC_GLOBAL, NEW_NAME, NULL)) ERR;
      if (nc_close(ncid)) ERR;

      /* Reopen the file and check again. */
      if (nc_open(FILE_NAME, NC_NOWRITE, &ncid)) ERR;

      if (nc_inq_attid(ncid, NC_GLOBAL, NEW_NAME, &attid)) ERR;
      if (attid != 0) ERR;
      if (nc_get_att_text(ncid, NC_GLOBAL, NEW_NAME, data_in)) ERR;
      if (strncmp(CONTENTS, data_in, strlen(CONTENTS))) ERR;
      if (nc_close(ncid)) ERR;

      free(data_in);
   }
   SUMMARIZE_ERR;
   printf("*** testing attribute renaming for a variable attribute...");
   {
#define OLD_NAME1 "Constantinople"
#define NEW_NAME1 "Istanbul____________"
#define CONTENTS1 "Lots of people!"

      int ncid, attid, varid;
      char *data_in;

      if (!(data_in = malloc(strlen(CONTENTS1) + 1))) ERR;

      /* Create a file with an att. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLOBBER, &ncid)) ERR;
      if (nc_def_var(ncid, VAR_NAME, NC_INT, 0, NULL, &varid)) ERR;
      if (nc_put_att_text(ncid, varid, OLD_NAME1, strlen(CONTENTS1),
                          CONTENTS1)) ERR;

      /* Rename the att. */
      if (nc_rename_att(ncid, varid, OLD_NAME1, NEW_NAME1)) ERR;

      /* Check the file. */
      if (nc_inq_attid(ncid, varid, NEW_NAME1, &attid)) ERR;
      if (attid != 0) ERR;
      if (nc_get_att_text(ncid, varid, NEW_NAME1, data_in)) ERR;
      if (strncmp(CONTENTS1, data_in, strlen(CONTENTS1))) ERR;
      if (nc_close(ncid)) ERR;

      /* Reopen the file and check again. */
      if (nc_open(FILE_NAME, NC_NOWRITE, &ncid)) ERR;

      if (nc_inq_attid(ncid, varid, NEW_NAME1, &attid)) ERR;
      if (attid != 0) ERR;
      if (nc_get_att_text(ncid, varid, NEW_NAME1, data_in)) ERR;
      if (strncmp(CONTENTS1, data_in, strlen(CONTENTS1))) ERR;
      if (nc_close(ncid)) ERR;

      free(data_in);
   }
   SUMMARIZE_ERR;
   FINAL_RESULTS;
}
