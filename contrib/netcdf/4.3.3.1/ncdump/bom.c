/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* BOM Sequences */
static char* U8   = "\xEF\xBB\xBF";    /* UTF-8 */
static char* BE32 = "\x00\x00\xFE\xFF"; /* UTF-32; big-endian */
static char* LE32 = "\xFF\xFE";       /* UTF-32; little-endian */
static char* BE16 = "\xFE\xFF";       /* UTF-16; big-endian */
static char* LE16 = "\xFF\xFE";       /* UTF-16; little-endian */

int
main(int argc, char** argv)
{
    char* bom = U8;
    size_t bomlen = 3;
    if(argc > 1 && strlen(argv[1]) > 0) {
	char* which = argv[1];
	switch (which[0]) {
	case '1': bom = BE16; bomlen = 2; break;
	case '3': bom = BE32; bomlen = 2; break;
	default: break;
	}
    }
    fwrite(bom,1,bomlen,stdout);
    exit(0);
}
