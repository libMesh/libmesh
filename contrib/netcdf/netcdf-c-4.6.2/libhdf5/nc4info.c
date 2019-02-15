/**
 * @file
 * @internal Add provenance info for netcdf-4 files.
 *
 * Copyright 2010, UCAR/Unidata See netcdf/COPYRIGHT file for copying
 * and redistribution conditions.
 * @author Dennis Heimbigner
 */

#include "config.h"
#include "nc4internal.h"
#include "hdf5internal.h"
#include "nclist.h"
#include "ncbytes.h"

/* Various Constants */
#define NCPROPS_MAX_NAME 1024 /* max key name size */
#define NCPROPS_MAX_VALUE 1024 /* max value size */
#define HDF5_MAX_NAME 1024 /**< HDF5 max name. */

#define ESCAPECHARS "\\=|,"

/** @internal Check NetCDF return code. */
#define NCHECK(expr) {if((expr)!=NC_NOERR) {goto done;}}

/** @internal Check HDF5 return code. */
#define HCHECK(expr) {if((expr)<0) {ncstat = NC_EHDFERR; goto done;}}

static int globalpropinitialized = 0;
struct NCPROPINFO globalpropinfo; /**< Global property info. */

/* Forward */
static int properties_parse(const char* text0, NClist* pairs);

/**
 * @internal Initialize default provenance info
 * This will only be used for newly created files
 * or for opened files that do not contain an _NCProperties
 * attribute.
 *
 * @return ::NC_NOERR No error.
 * @author Dennis Heimbigner
 */
int
NC4_provenance_init(void)
{
    int stat = NC_NOERR;
    int i;
    NClist* other = NULL;
    char* name = NULL;
    char* value = NULL;
    unsigned major,minor,release;

    if(globalpropinitialized)
	return stat;

    /* Build _NCProperties info */

    /* Initialize globalpropinfo */
    memset((void*)&globalpropinfo,0,sizeof(globalpropinfo));
    globalpropinfo.version = NCPROPS_VERSION;
    globalpropinfo.properties = nclistnew();
    if(globalpropinfo.properties == NULL)
	{stat = NC_ENOMEM; goto done;}

    /* Insert primary library version as first entry */
    if((name = strdup(NCPNCLIB2)) == NULL)
	{stat = NC_ENOMEM; goto done;}
    nclistpush(globalpropinfo.properties,name);
    name = NULL; /* Avoid multiple free() */

    if((value = strdup(PACKAGE_VERSION)) == NULL)
	{stat = NC_ENOMEM; goto done;}
    nclistpush(globalpropinfo.properties,value);
    value = NULL;
    
    /* Insert the HDF5 as underlying storage format library */
    if((name = strdup(NCPHDF5LIB2)) == NULL)
	{stat = NC_ENOMEM; goto done;}
    nclistpush(globalpropinfo.properties,name);
    name = NULL;

    stat = NC4_hdf5get_libversion(&major,&minor,&release);
    if(stat) goto done;
    {
	char sversion[64];
        snprintf(sversion,sizeof(sversion),"%1u.%1u.%1u",major,minor,release);
	if((value = strdup(sversion)) == NULL)
	    {stat = NC_ENOMEM; goto done;}
    }
    nclistpush(globalpropinfo.properties,value);
    value = NULL;

    /* Add any extra fields */
    /*Parse them into an NClist */
    other = nclistnew();
    if(other == NULL) {stat = NC_ENOMEM; goto done;}
#ifdef NCPROPERTIES
    stat = properties_parse(NCPROPERTIES_EXTRA,other);
    if(stat) goto done;
#endif
    /* merge into the properties list */
    for(i=0;i<nclistlength(other);i++)
	nclistpush(globalpropinfo.properties,strdup(nclistget(other,i)));
    nclistfreeall(other);
    other = NULL;

done:
    if(name != NULL) free(name);
    if(value != NULL) free(value);
    if(other != NULL)
	nclistfreeall(other);    
    if(stat && globalpropinfo.properties != NULL) {
	nclistfreeall(globalpropinfo.properties);
        globalpropinfo.properties = NULL;
    }
    if(stat == NC_NOERR)
        globalpropinitialized = 1; /* avoid repeating it */
    return stat;
}

/* Locate a specific character and return its pointer
   or EOS if not found
   take \ escapes into account */
static char*
locate(char* p, char tag)
{
    char* next;
    int c;
    assert(p != NULL);
    for(next = p;(c = *next);next++) {
	if(c == tag)
	    return next;
	else if(c == '\\' && next[1] != '\0')
	    next++; /* skip escaped char */
    }
    return next; /* not found */
}

/**
 * @internal finalize default provenance info
 *
 * @return ::NC_NOERR No error.
 * @author Dennis Heimbigner
 */
int
NC4_provenance_finalize(void)
{
    nclistfreeall(globalpropinfo.properties);
    return NC_NOERR;
}

/**
 * @internal Parse file properties.
 *
 * @param text0 Text properties.
 * @param pairs list of parsed (key,value) pairs
 *
 * @return ::NC_NOERR No error.
 * @author Dennis Heimbigner
 */
static int
properties_parse(const char* text0, NClist* pairs)
{
    int ret = NC_NOERR;
    char* p;
    char* q;
    char* text = NULL;

    if(text0 == NULL || strlen(text0) == 0)
	goto done;

    text = strdup(text0);
    if(text == NULL) return NC_ENOMEM;

    /* For back compatibility with version 1, translate '|' -> ',' */
    for(p=text;*p;p++) {
	if(*p == NCPROPSSEP1)
	    *p = NCPROPSSEP2;
    }

    /* Walk and fill in ncinfo */
    p = text;
    while(*p) {
	char* name = p;
	char* value = NULL;
	char* next = NULL;

	/* Delimit whole (key,value) pair */
	q = locate(p,NCPROPSSEP2);
	if(*q != '\0') /* Never go beyond the final nul term */
  	    *q++ = '\0';
	next = q;
	/* split key and value */
	q = locate(p,'=');      
 	name = p;
        *q++ = '\0';
	value = q;
	/* Set up p for next iteration */
	p = next;
	nclistpush(pairs,strdup(name));
	nclistpush(pairs,strdup(value));
    }
done:
    if(text) free(text);
    return ret;
}


/* Utility to transfer a string to a buffer with escaping */
static void
escapify(NCbytes* buffer, const char* s)
{
    const char* p;
    for(p=s;*p;p++) {
	if(strchr(ESCAPECHARS,*p) != NULL)
	    ncbytesappend(buffer,'\\');
	ncbytesappend(buffer,*p);
    }
}

/* Utility to copy contents of the dfalt into an NCPROPINFO object */
static int
propinfo_default(struct NCPROPINFO* dst, const struct NCPROPINFO* dfalt)
{
    int i;
    if(dst->properties == NULL) {
	dst->properties = nclistnew();
	if(dst->properties == NULL) return NC_ENOMEM;
    }
    dst->version = dfalt->version;
    for(i=0;i<nclistlength(dfalt->properties);i++) {
	char* s = nclistget(dfalt->properties,i);
	s = strdup(s);
	if(s == NULL) return NC_ENOMEM;
        nclistpush(dst->properties,s);
    }
    return NC_NOERR;
}

/**
 * @internal Build _NCProperties attribute value.
 *
 * Convert a NCPROPINFO instance to a single string.
 *
 * @param info Properties info.
 * @param propdatap Pointer that gets properties string.
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_EINVAL failed.
 * @author Dennis Heimbigner
 */
int
NC4_buildpropinfo(struct NCPROPINFO* info, char** propdatap)
{
    int stat = NC_NOERR;
    int i;
    NCbytes* buffer = NULL;
    char sversion[64];

    if(info == NULL || info->version == 0 || propdatap == NULL)
      {stat = NC_EINVAL; goto done;}

    *propdatap = NULL;

    buffer = ncbytesnew();
    if(!buffer) {stat = NC_ENOMEM; goto done;}

    /* start with version */
    ncbytescat(buffer,NCPVERSION);
    ncbytesappend(buffer,'=');
    snprintf(sversion,sizeof(sversion),"%d",info->version);
    ncbytescat(buffer,sversion);

    for(i=0;i<nclistlength(info->properties);i+=2) {
	char* value, *name;
	name = nclistget(info->properties,i);
	if(name == NULL) continue;
	value = nclistget(info->properties,i+1);
        ncbytesappend(buffer,NCPROPSSEP2); /* terminate last entry */
	escapify(buffer,name);
        ncbytesappend(buffer,'=');
	escapify(buffer,value);
    }
    /* Force null termination */
    ncbytesnull(buffer);
    *propdatap = ncbytesextract(buffer);

done:
    if(buffer != NULL) ncbytesfree(buffer);
    return stat;
}

#if 0
/**
 * @internal Write the properties attribute to file.
 *
 * @param h5 Pointer to HDF5 file info struct.
 *
 * @return ::NC_NOERR No error.
 * @author Dennis Heimbigner
 */
int
NC4_put_ncproperties(NC_FILE_INFO_T* file)
{
    int ncstat = NC_NOERR;
    char* text = NULL;

    /* Get root group */
    grp = ((NC_HDF5_GRP_INFO_T *)(h5->root_grp->format_grp_info))->hdf_grpid;
    /* See if the NCPROPS attribute exists */
    if(H5Aexists(grp,NCPROPS) <= 0) { /* Does not exist */
      ncstat = NC4_buildpropinfo(&h5->fileinfo->propattr,&text);
      if(text == NULL || ncstat != NC_NOERR) {
        goto done;
      }
      /* Create a datatype to refer to. */
      HCHECK((atype = H5Tcopy(H5T_C_S1)));
      HCHECK((H5Tset_cset(atype, H5T_CSET_ASCII)));
      HCHECK((H5Tset_size(atype, strlen(text)+1))); /*keep nul term */
      HCHECK((aspace = H5Screate(H5S_SCALAR)));
      HCHECK((attid = H5Acreate(grp, NCPROPS, atype, aspace, H5P_DEFAULT)));
      HCHECK((H5Awrite(attid, atype, text)));
    }
 done:
    if(text != NULL) {
      free(text);
      text = NULL;
    }

    if(attid >= 0) HCHECK((H5Aclose(attid)));
    if(aspace >= 0) HCHECK((H5Sclose(aspace)));
    if(atype >= 0) HCHECK((H5Tclose(atype)));
    return ncstat;
}
#endif

/**
 * @internal
 *
 * Construct the provenance information for a newly created file
 * using dfalt as the default.
 * Note that creation of the _NCProperties attribute is deferred
 * to the sync_netcdf4_file function.
 *
 * @param file Pointer to file object.
 * @param dfalt
 *
 * @return ::NC_NOERR No error.
 * [Note: other errors are reported via LOG()]
 * @author Dennis Heimbigner
 */
int
NC4_set_provenance(NC_FILE_INFO_T* file, const struct NCPROPINFO* dfalt)
{
    int ncstat = NC_NOERR;
    struct NCPROVENANCE* provenance = NULL;
    int superblock = -1;

    assert(file->provenance == NULL);
    provenance = calloc(1,sizeof(struct NCPROVENANCE));
    if(provenance == NULL) {ncstat = NC_ENOMEM; goto done;}

    /* Initialize from the default */
    provenance->propattr.version = globalpropinfo.version;
    /* Get the superblock number */
    if((ncstat = NC4_hdf5get_superblock(file,&superblock)))
	goto done;
    provenance->superblockversion = superblock;

    /* Capture properties */
    provenance->propattr.properties = nclistnew();
    if(provenance->propattr.properties == NULL)
	{ncstat = NC_ENOMEM; goto done;}
    /* add in the dfalt values */
    if(dfalt != NULL) {
	int i;
	for(i=0;i<nclistlength(dfalt->properties);i++) {
	    char* prop = nclistget(dfalt->properties,i);
	    if(prop != NULL) {
		prop = strdup(prop);
		if(prop == NULL) {ncstat = NC_ENOMEM; goto done;}
	        nclistpush(provenance->propattr.properties,prop);
	    }
	}
    }

done:
    if(ncstat) {
	LOG((0,"Could not create _NCProperties attribute"));
	(void)NC4_free_provenance(provenance);
    } else
        file->provenance = provenance;
    return NC_NOERR;
}

/**
 * @internal
 *
 * Construct the provenance information for a newly opened file
 * Using the specified _NCProperties value. If NULL, then
 * initialize using dfalt.
 *
 * @param file Pointer to file object.
 * @param propstring The contents of _NCProperties
 * @param dfalt
 *
 * @return ::NC_NOERR No error.
 * @return ::NC_ENOMEM
 * @return ::NC_EINVAL
 * @author Dennis Heimbigner
 */
int
NC4_get_provenance(NC_FILE_INFO_T* file, const char* propstring, const struct NCPROPINFO* dfalt)
{
    int ncstat = NC_NOERR;
    struct NCPROVENANCE* provenance;
    char *name = NULL;
    char *value = NULL;
    int v = 0;
    int superblock = -1;

    assert(file->provenance == NULL);
    if((file->provenance = calloc(1,sizeof(struct NCPROVENANCE))) == NULL)
	{ncstat = NC_ENOMEM; goto done;}
    provenance = file->provenance;
    if((provenance->propattr.properties = nclistnew()) == NULL)
	{ncstat = NC_ENOMEM; goto done;}

    /* Set the superblock */
    if((ncstat = NC4_hdf5get_superblock(file,&superblock)))
	goto done;
    provenance->superblockversion = superblock;

    if(propstring == NULL) {
	/* Use dfalt */
	if((ncstat=propinfo_default(&provenance->propattr,dfalt)))
	    goto done;
    } else {
	NClist* list = provenance->propattr.properties;
        if((ncstat=properties_parse(propstring,list)))
	    goto done;
	/* Check the version and remove from properties list*/
        if(nclistlength(list) < 2)
	    {ncstat = NC_EINVAL; goto done;} /* bad _NCProperties attribute */
	/* Extract the purported version=... */
	name = nclistremove(list,0);
        value = nclistremove(list,0);
	if(strcmp(name,NCPVERSION) == 0) {
            if(sscanf(value,"%d",&v) != 1)
	        {ncstat = NC_EINVAL; goto done;} /* illegal version */
            if(v <= 0 || v > NCPROPS_VERSION)
	        {ncstat = NC_EINVAL; goto done;} /* unknown version */
	    provenance->propattr.version = v;
	} else
	    {ncstat = NC_EINVAL; goto done;} /* bad _NCProperties attribute */
#if 0
        /* Now, rebuild from version 1 to version 2 if necessary */
        if(provenance->propattr.version == 1) {
	    int i;
	    for(i=0;i<nclistlength(list);i+=2) {
	        char* newname = NULL;
	        name = nclistget(list,i);
	        if(name == NULL) continue; /* ignore */
	        if(strcmp(name,NCPNCLIB1) == 0)
		    newname = NCPNCLIB2; /* change name */
	        else if(strcmp(name,NCPHDF5LIB1) == 0)
		    newname = NCPHDF5LIB2;
		else continue; /* ignore */
		/* Do any rename */
	        nclistset(list,i,strdup(newname));
	        if(name) free(name);
	    }
        }
#endif
    }
done:
    if(name != NULL) free(name);
    if(value != NULL) free(value);
    return ncstat;
}

/**
 * @internal
 *
 * Free the NCPROVENANCE object
 * @param prov Pointer to provenance object
 *
 * @return ::NC_NOERR No error.
 * @author Dennis Heimbigner
 */
int
NC4_free_provenance(struct NCPROVENANCE* prov)
{
    if(prov == NULL) return NC_NOERR;
    if(prov->propattr.properties != NULL)
	nclistfreeall(prov->propattr.properties);
    prov->propattr.properties = NULL;
    free(prov);
    return NC_NOERR;
}

/* HDF5 Specific attribute read/write of _NCProperties */
int
NC4_read_ncproperties(NC_FILE_INFO_T* h5)
{
    int retval = NC_NOERR;
    hid_t hdf5grpid = -1;
    hid_t attid = -1;
    hid_t aspace = -1;
    hid_t atype = -1;
    hid_t ntype = -1;
    char* text = NULL;
    H5T_class_t t_class;
    hsize_t size;

    hdf5grpid = ((NC_HDF5_GRP_INFO_T *)(h5->root_grp->format_grp_info))->hdf_grpid;

    if(H5Aexists(hdf5grpid,NCPROPS) <= 0) { /* Does not exist */
	/* File did not contain a _NCProperties attribute */		
        retval=NC4_get_provenance(h5,NULL,&globalpropinfo);
        goto done;
    }

    /* NCPROPS Attribute exists, make sure it is legitimate */
    attid = H5Aopen_name(hdf5grpid, NCPROPS);
    assert(attid > 0);
    aspace = H5Aget_space(attid);
    atype = H5Aget_type(attid);
    /* Verify atype and size */
    t_class = H5Tget_class(atype);
    if(t_class != H5T_STRING)
	{retval = NC_EINVAL; goto done;}
    size = H5Tget_size(atype);
    if(size == 0)
	{retval = NC_EINVAL; goto done;}
    text = (char*)malloc(1+(size_t)size);
    if(text == NULL)
	{retval = NC_ENOMEM; goto done;}
    if((ntype = H5Tget_native_type(atype, H5T_DIR_DEFAULT)) < 0)
	{retval = NC_EHDFERR; goto done;}
    if((H5Aread(attid, ntype, text)) < 0)
	{retval = NC_EHDFERR; goto done;}
    /* Make sure its null terminated */
    text[(size_t)size] = '\0';
    /* Process the _NCProperties value */
    if((retval = NC4_get_provenance(h5, text, &globalpropinfo)))
	goto done;

done:
    if(text != NULL) free(text);
    /* Close out the HDF5 objects */
    if(attid > 0 && H5Aclose(attid) < 0) retval = NC_EHDFERR;
    if(aspace > 0 && H5Sclose(aspace) < 0) retval = NC_EHDFERR;
    if(atype > 0 && H5Tclose(atype) < 0) retval = NC_EHDFERR;
    if(ntype > 0 && H5Tclose(ntype) < 0) retval = NC_EHDFERR;

    /* For certain errors, actually fail, else log that attribute was invalid and ignore */
    if(retval != NC_ENOMEM && retval != NC_EHDFERR) {
	LOG((0,"Invalid _NCProperties attribute"));
	retval = NC_NOERR;
    }
    return retval;
}

int
NC4_write_ncproperties(NC_FILE_INFO_T* h5)
{
    int retval = NC_NOERR;
    hid_t hdf5grpid = -1;
    hid_t attid = -1;
    hid_t aspace = -1;
    hid_t atype = -1;
    char* text = NULL;
    size_t len = 0;

    /* If the file is read-only, return an error. */
    if (h5->no_write)
      {retval = NC_EPERM; goto done;}

    hdf5grpid = ((NC_HDF5_GRP_INFO_T *)(h5->root_grp->format_grp_info))->hdf_grpid;

    if(H5Aexists(hdf5grpid,NCPROPS) > 0) /* Already exists, no overwrite */
	goto done;

    /* Build the attribute string */
    if((retval = NC4_buildpropinfo(&h5->provenance->propattr,&text)))
	goto done;

    /* Build the HDF5 string type */
    if ((atype = H5Tcopy(H5T_C_S1)) < 0)
	{retval = NC_EHDFERR; goto done;}
    if (H5Tset_strpad(atype, H5T_STR_NULLTERM) < 0)
	{retval = NC_EHDFERR; goto done;}
    if(H5Tset_cset(atype, H5T_CSET_ASCII) < 0)
	{retval = NC_EHDFERR; goto done;}

    /* Create NCPROPS attribute */

   len = strlen(text);
   if(H5Tset_size(atype, len) < 0)
      {retval = NC_EFILEMETA; goto done;}
   if((aspace = H5Screate(H5S_SCALAR)) < 0)
      {retval = NC_EFILEMETA; goto done;}
   if ((attid = H5Acreate(hdf5grpid, NCPROPS, atype, aspace, H5P_DEFAULT)) < 0)
      {retval = NC_EFILEMETA; goto done;}
   if (H5Awrite(attid, atype, text) < 0)
      {retval = NC_EFILEMETA; goto done;}

done:
    if(text != NULL) free(text);
    /* Close out the HDF5 objects */
    if(attid > 0 && H5Aclose(attid) < 0) retval = NC_EHDFERR;
    if(aspace > 0 && H5Sclose(aspace) < 0) retval = NC_EHDFERR;
    if(atype > 0 && H5Tclose(atype) < 0) retval = NC_EHDFERR;

    /* For certain errors, actually fail, else log that attribute was invalid and ignore */
    switch (retval) {
    case NC_ENOMEM:
    case NC_EHDFERR:
    case NC_EPERM:
    case NC_EFILEMETA:
    case NC_NOERR:
	break;
    default:
	LOG((0,"Invalid _NCProperties attribute"));
	retval = NC_NOERR;
	break;
    }
    return retval;
}

/* Debugging */

void
ncprintpropinfo(struct NCPROPINFO* info)
{
    int i;
    fprintf(stderr,"[%p] version=%d\n",info,info->version);
    for(i=0;i<nclistlength(info->properties);i+=2) {
	char* name = nclistget(info->properties,i);
	char* value = nclistget(info->properties,i+1);
	fprintf(stderr,"\t[%d] name=|%s| value=|%s|\n",i,name,value);
    }    
}

void
ncprintprovenance(struct NCPROVENANCE* prov)
{
    fprintf(stderr,"[%p] superblockversion=%d\n",prov,prov->superblockversion);
    ncprintpropinfo(&prov->propattr);
}

