#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "netcdf.h"
#include "netcdf_mem.h"
#include "ncbytes.h"
#include "nc_tests.h"

#ifdef USE_NETCDF4
#include <hdf5.h>
extern int H5Eprint1(FILE * stream);
#endif

static int
readfile(const char* path, NC_memio* memio)
{
    int status = NC_NOERR;
    FILE* f = NULL;
    size_t filesize = 0;
    size_t count = 0;
    char* memory = NULL;
    char* p = NULL;

    /* Open the file for reading */
#ifdef _MSC_VER
    f = fopen(path,"rb");
#else
    f = fopen(path,"r");
#endif
    if(f == NULL)
	{status = errno; goto done;}
    /* get current filesize */
    if(fseek(f,0,SEEK_END) < 0)
	{status = errno; goto done;}
    filesize = (size_t)ftell(f);
    /* allocate memory */
    memory = malloc((size_t)filesize);
    if(memory == NULL)
	{status = NC_ENOMEM; goto done;}
    /* move pointer back to beginning of file */
    rewind(f);
    count = filesize;
    p = memory;
    while(count > 0) {
        size_t actual;
        actual = fread(p,1,count,f);
	if(actual == 0 || ferror(f))
	    {status = NC_EIO; goto done;}	 
	count -= actual;
	p += actual;
    }
    if(memio) {
	memio->size = (size_t)filesize;
	memio->memory = memory;
    }    
done:
    if(status != NC_NOERR && memory != NULL)
	free(memory);
    if(f != NULL) fclose(f);
    return status;    
}

int
main(int argc, char** argv)
{
    int retval = NC_NOERR;
    int ncid;
    NC_memio mem;
    char* path = "f03tst_open_mem.nc";

    if(argc > 1)
	path = argv[1];

    if((retval=readfile(path,&mem)))
	goto exit;

    if((retval = nc_open_mem("mem", NC_INMEMORY|NC_NETCDF4, mem.size, mem.memory, &ncid)))
	goto exit;
    if((retval = nc_close(ncid)))
	goto exit;
    return 0;
exit:
    fprintf(stderr,"retval=%d\n",retval);
    return 1;
}
