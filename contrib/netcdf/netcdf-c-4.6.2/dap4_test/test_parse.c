/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

/**
This provides a simple dap4  metadata -> xml printer.
Used to test the parser
*/

#include "test_common.h"

int
main(int argc, char** argv)
{
    int ret = NC_NOERR;

    setup(TDMR_PARSE,argc,argv);

    if((ret = NCD4_parse(metadata))) goto done;
    ret = NCD4_print(metadata,output);
    ncbytesnull(output);
    if(ret == NC_NOERR) {
        fprintf(stdout,"%s\n",ncbytescontents(output));
	fflush(stdout);
    }
done:
    return cleanup(ret);
}
