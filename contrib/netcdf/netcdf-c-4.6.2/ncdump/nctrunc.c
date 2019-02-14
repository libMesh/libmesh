/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define BUFLEN 100000

int
main(int argc, char** argv)
{
    unsigned char buffer[BUFLEN];
    size_t count, red, avail, trunc;
    unsigned char* p = buffer;
    long i;
    FILE* input = stdin;

    if(argc > 1) {
	input = fopen(argv[1],"r");
    }

    /* Read the whole file */
    p = buffer;
    red = 0;
    avail = BUFLEN;
    for(;;) {
	count = fread(p,1,avail,input);
	p += count;
	red += count;
	avail -= count;
	if(feof(input)) break;
    }
    trunc = red;
    for(i=(red-1);i>=0;i--) {
	if(buffer[i] != '\0') {
	    trunc = i + 1;
	    break;	    
	}
    }
    p = buffer;
    avail = trunc;
    for(;;) {
	count = fwrite(p,1,avail,stdout);
	if(count == 0) break;
	p += count;
	avail -= count;
	if(avail == 0) break;	
    }
    if(avail > 0)
        exit(1);
    else
	exit(0);
}
