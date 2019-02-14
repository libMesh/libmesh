#define XSTRINGIFY(s) #s
#define STRINGIFY(s) XSTRINGIFY(s)

static const char*
gettopsrcdir(void)
{
    const char* topsrcdir = NULL;
#ifdef TOPSRCDIR
    topsrcdir = STRINGIFY(TOPSRCDIR);
#else
    static char tsd[4096];
    extern char *getcwd(char *buf, size_t size);
    tsd[0] = '\0';
    getcwd(tsd,sizeof(tsd));
    if(strlen(tsd) > 0) {
        strcat(tsd,"/..");
        topsrcdir = tsd;
    }
#endif
    if(topsrcdir == NULL) {
        fprintf(stderr,"*** FAIL: $abs_top_srcdir not defined\n");
        exit(1);
    }    
    fprintf(stderr,"topsrcdir=%s\n",topsrcdir);
    return topsrcdir;
}
