/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncgen/genlib.c,v 1.57 2010/04/04 19:39:47 dmh Exp $
 *********************************************************************/

#include "includes.h"

/* invoke netcdf calls (or generate C or Fortran code) to create netcdf
 * from in-memory structure.
*/
void
define_netcdf(void)
{

    /* Execute exactly one of these */
#ifdef ENABLE_C
    if (l_flag == L_C) genc_netcdf(); else /* create C code to create netcdf */
#endif
#ifdef ENABLE_F77
    if (l_flag == L_F77) genf77_netcdf(); else /* create Fortran code */
#endif
#ifdef ENABLE_JAVA
    if(l_flag == L_JAVA) genjava_netcdf(); else
#endif
/* Binary is the default */
#ifdef ENABLE_BINARY
    genbin_netcdf(); /* create netcdf */
#else
    derror("No language specified");
#endif
    close_netcdf();
    cleanup();
}

void
close_netcdf(void)
{
#ifdef ENABLE_C
    if (l_flag == L_C) genc_close(); else /* create C code to close netcdf */
#endif
#ifdef ENABLE_F77
    if (l_flag == L_F77) genf77_close(); else
#endif
#ifdef ENABLE_JAVA
    if (l_flag == L_JAVA) genjava_close(); else
#endif
#ifdef ENABLE_BINARY
    if (l_flag == L_BINARY) genbin_close();
#endif
}

/**
Return a string representing
the fully qualified name of the symbol.
Symbol must be top level
Caller must free.
*/
void
topfqn(Symbol* sym)
{
#ifdef USE_NETCDF4
    char* fqn;
    char* fqnname;
    char* parentfqn;
    Symbol* parent;
#endif

    if(sym->fqn != NULL)
	return; /* already defined */

#ifdef USE_NETCDF4
    if(!usingclassic) {
        parent = sym->container;
        /* Recursively compute parent fqn */
        if(parent == NULL) { /* implies this is the rootgroup */
            assert(sym->grp.is_root);
            sym->fqn = estrdup("");
            return;
        } else if(parent->fqn == NULL) {
            topfqn(parent);
        }
        parentfqn = parent->fqn;

        fqnname = fqnescape(sym->name);
        fqn = (char*)ecalloc(strlen(fqnname) + strlen(parentfqn) + 1 + 1);
        strcpy(fqn,parentfqn);
        strcat(fqn,"/");
        strcat(fqn,fqnname);
        sym->fqn = fqn;
    } else
#endif /*USE_NETCDF4*/
    {
	sym->fqn = strdup(sym->name);
    }
}

/**
Return a string representing
the fully qualified name of a nested symbol
(i.e. field or econst).
Caller must free.
*/
void
nestedfqn(Symbol* sym)
{
    char* fqn;
    char* fqnname;
    Symbol* parent;

    if(sym->fqn != NULL)
	return; /* already defined */

    /* Parent must be a type */
    parent = sym->container;
    assert (parent->objectclass == NC_TYPE);

    assert(parent->fqn != NULL);

    fqnname = fqnescape(sym->name);
    fqn = (char*)ecalloc(strlen(fqnname) + strlen(parent->fqn) + 1 + 1);
    strcpy(fqn,parent->fqn);
    strcat(fqn,".");
    strcat(fqn,fqnname);
    sym->fqn = fqn;
}

/**
Return a string representing
the fully qualified name of an attribute.
Caller must free.
*/
void
attfqn(Symbol* sym)
{
    char* fqn;
    char* fqnname;
    char* parentfqn;
    Symbol* parent;

    if(sym->fqn != NULL)
	return; /* already defined */

    assert (sym->objectclass == NC_ATT);

    parent = sym->container;
    if(parent == NULL)
	parentfqn = "";
    else
	parentfqn = parent->fqn;

    fqnname = fqnescape(sym->name);
    fqn = (char*)ecalloc(strlen(fqnname) + strlen(parentfqn) + 1 + 1);
    strcpy(fqn,parentfqn);
    strcat(fqn,"_");
    strcat(fqn,fqnname);
    sym->fqn = fqn;
}

#if 0
/* Result is pool alloc'd*/
char*
cprefixed(List* prefix, char* suffix, char* separator)
{
    int slen;
    int plen;
    int i;
    char* result;

    ASSERT(suffix != NULL);
    plen = prefixlen(prefix);
    if(prefix == NULL || plen == 0) return codify(suffix);
    /* plen > 0*/
    slen = 0;
    for(i=0;i<plen;i++) {
	Symbol* sym = (Symbol*)listget(prefix,i);
	slen += (strlen(sym->name)+strlen(separator));
    }
    slen += strlen(suffix);
    slen++; /* for null terminator*/
    result = poolalloc(slen);
    result[0] = '\0';
    /* Leave off the root*/
    i = (rootgroup == (Symbol*)listget(prefix,0))?1:0;
    for(;i<plen;i++) {
	Symbol* sym = (Symbol*)listget(prefix,i);
        strcat(result,sym->name); /* append "<prefix[i]/>"*/
	strcat(result,separator);
    }
    strcat(result,suffix); /* append "<suffix>"*/
    return result;
}
#endif /*0*/
