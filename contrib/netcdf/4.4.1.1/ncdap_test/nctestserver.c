#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "ncdispatch.h"

/**
Given a partial suffix path,
try to find a server for which
a request to that server + path
returns some kind of result.
This indicates that the server is up
and running.
Return the complete url for the server
plus the path.
*/

int
main(int argc, char** argv)
{
    const char* url = NULL;
    const char* path = NULL;
    char* serverlist[64];
    int nservers = 0;

    if(argc == 1)
	path = "";
    else if(argc >= 2)
	path = argv[1];
    if(argc >= 3) {
	int i;
	for(i=2;i<argc;i++,nservers++)
	    serverlist[i-2] = argv[i];
        serverlist[nservers] = NULL;
    }
    url = NC_findtestserver(path,(nservers==0?(const char**)NULL:(const char**)serverlist));
    if(url == NULL)
	url = "";
    printf("%s",url);
    fflush(stdout);
    exit(0);
}
