/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#include "config.h"
#include "ocinternal.h"
#include "occompile.h"
#include "ocdebug.h"

static const unsigned int MAX_UINT = 0xffffffff;

static OCerror mergedas1(OCnode* dds, OCnode* das);
static OCerror mergedods1(OCnode* dds, OCnode* das);
static OCerror mergeother1(OCnode* root, OCnode* das);
static char* pathtostring(OClist* path, char* separator);
static void computefullname(OCnode* node);

/* Process ocnodes to fix various semantic issues*/
void
occomputesemantics(OClist* ocnodes)
{
    unsigned int i,j;
    OCASSERT((ocnodes != NULL));
    for(i=0;i<oclistlength(ocnodes);i++) {
	OCnode* node = (OCnode*)oclistget(ocnodes,i);
	/* set the container for dims*/
	if(node->octype == OC_Dimension && node->dim.array != NULL) {
	    node->container = node->dim.array->container;
	}
    }
    /* Fill in array.sizes */
    for(i=0;i<oclistlength(ocnodes);i++) {
	OCnode* node = (OCnode*)oclistget(ocnodes,i);
	if(node->array.rank > 0) {
	    node->array.sizes = (size_t*)malloc(node->array.rank*sizeof(size_t));
	    for(j=0;j<node->array.rank;j++) {
		OCnode* dim = (OCnode*)oclistget(node->array.dimensions,j);
		node->array.sizes[j] = dim->dim.declsize;		
	    }
	}
    }
}

void
occomputefullnames(OCnode* root)
{
    unsigned int i;
    if(root->name != NULL) computefullname(root);
    if(root->subnodes != NULL) { /* recurse*/
        for(i=0;i<oclistlength(root->subnodes);i++) {
	    OCnode* node = (OCnode*)oclistget(root->subnodes,i);
	    occomputefullnames(node);
	}
    }
}

static void
computefullname(OCnode* node)
{
    char* tmp;
    char* fullname;
    OClist* path;

    OCASSERT((node->name != NULL));
    if(node->fullname != NULL)
	return;
    path = oclistnew();
    occollectpathtonode(node,path);
    tmp = pathtostring(path,PATHSEPARATOR);
    if(tmp == NULL) {
        fullname = nulldup(node->name);
    } else {
        fullname = tmp;
    }
    node->fullname = fullname;
    oclistfree(path);
}

/* Convert path to a string */
static char*
pathtostring(OClist* path, char* separator)
{
    int slen,i,len;
    char* pathname;
    if(path == NULL) return NULL;
    len = oclistlength(path);
    if(len == 0) return NULL;
    for(i=0,slen=0;i<len;i++) {
	OCnode* node = (OCnode*)oclistget(path,(size_t)i);
	if(node->container == NULL || node->name == NULL) continue;
	slen += strlen(node->name);
    }
    slen += ((len-1)*strlen(separator));
    slen += 1;   /* for null terminator*/
    pathname = (char*)ocmalloc((size_t)slen);
    MEMCHECK(pathname,NULL);
    pathname[0] = '\0';
    for(i=0;i<len;i++) {
	OCnode* node = (OCnode*)oclistget(path,(size_t)i);
	if(node->container == NULL || node->name == NULL) continue;
	if(strlen(pathname) > 0) strcat(pathname,separator);
        strcat(pathname,node->name);
    }
    return pathname;
}

/* Collect the set of nodes ending in "node"*/
void
occollectpathtonode(OCnode* node, OClist* path)
{
    if(node == NULL) return;
    occollectpathtonode(node->container,path);
    oclistpush(path,(void*)node);
}

OCnode*
ocnode_new(char* name, OCtype ptype, OCnode* root)
{
    OCnode* cdf = (OCnode*)ocmalloc(sizeof(OCnode));
    MEMCHECK(cdf,(OCnode*)NULL);
    memset((void*)cdf,0,sizeof(OCnode));
    cdf->header.magic = OCMAGIC;
    cdf->header.occlass = OC_Node;
    cdf->name = (name?nulldup(name):NULL);
    cdf->octype = ptype;
    cdf->array.dimensions = NULL;
    cdf->root = root;
    return cdf;
}

static OCattribute*
makeattribute(char* name, OCtype ptype, OClist* values)
{
    OCattribute* att = (OCattribute*)ocmalloc(sizeof(OCattribute)); /* ocmalloc zeros*/
    MEMCHECK(att,(OCattribute*)NULL);
    att->name = nulldup(name);
    att->etype = ptype;
    att->nvalues = oclistlength(values);
    att->values = NULL;
    if(att->nvalues > 0) {
	int i;
        att->values = (char**)ocmalloc(sizeof(char*)*att->nvalues);
        for(i=0;i<att->nvalues;i++)
	    att->values[i] = nulldup((char*)oclistget(values,(size_t)i));
    }
    return att;
}

void
ocroot_free(OCnode* root)
{
    OCtree* tree;
    OCstate* state;
    int i;

    if(root == NULL || root->tree == NULL) return;

    tree = root->tree;
    state = tree->state;

    /* Free up the OCDATA instance, if any */
    if(tree->data.data != NULL)
	ocdata_free(state,tree->data.data);

    for(i=0;i<oclistlength(state->trees);i++) {
	OCnode* node = (OCnode*)oclistget(state->trees,(size_t)i);
	if(root == node)
	    oclistremove(state->trees,(size_t)i);
    }
    /* Note: it is ok if state->trees does not contain this root */    
    octree_free(tree);
}

void
octree_free(OCtree* tree)
{
    if(tree == NULL) return;
    ocnodes_free(tree->nodes);
    ocfree(tree->constraint);
    ocfree(tree->text);
    if(tree->data.xdrs != NULL) {
        xxdr_free(tree->data.xdrs);
    }
    ocfree(tree->data.filename); /* may be null */
    if(tree->data.file != NULL) fclose(tree->data.file);
    ocfree(tree->data.memory);
    ocfree(tree);
}

void
ocnodes_free(OClist* nodes)
{
    unsigned int i,j;
    for(i=0;i<oclistlength(nodes);i++) {
	OCnode* node = (OCnode*)oclistget(nodes,i);
        ocfree(node->name);
        ocfree(node->fullname);
        while(oclistlength(node->att.values) > 0) {
	    char* value = (char*)oclistpop(node->att.values);
	    ocfree(value);
        }
        while(oclistlength(node->attributes) > 0) {
            OCattribute* attr = (OCattribute*)oclistpop(node->attributes);
	    ocfree(attr->name);
#if 0
	    /* If the attribute type is string, then we need to free them*/
all values are strings now
	    if(attr->etype == OC_String || attr->etype == OC_URL)
#endif
	    {
		char** strings = (char**)attr->values;
		for(j=0;j<attr->nvalues;j++) {ocfree(*strings); strings++;}
	    }
	    ocfree(attr->values);
	    ocfree(attr);
        }
        if(node->array.dimensions != NULL) oclistfree(node->array.dimensions);
        if(node->subnodes != NULL) oclistfree(node->subnodes);
        if(node->att.values != NULL) oclistfree(node->att.values);
        if(node->attributes != NULL) oclistfree(node->attributes);
	if(node->array.sizes != NULL) free(node->array.sizes);
        ocfree(node);
    }
    oclistfree(nodes);
}

/*
In order to be as compatible as possible with libdap,
we try to use the same algorithm for DAS->DDS matching.
As described there, the algorithm is as follows.
    If the [attribute] name contains one or
    more field separators then look for a [DDS]variable whose
    name matches exactly. If the name contains no field separators then
    the look first in the top level [of the DDS] and then in all subsequent
    levels and return the first occurrence found. In general, this
    searches constructor types in the order in which they appear
    in the DDS, but there is no requirement that it do so.

    Note: If a dataset contains two constructor types which have field names
    that are the same (say point.x and pair.x) one should use fully qualified
    names to get each of those variables.

*/

OCerror
ocddsdasmerge(OCstate* state, OCnode* dasroot, OCnode* ddsroot)
{
    OCerror stat = OC_NOERR;
    OClist* dasglobals = oclistnew();
    OClist* dodsglobals = oclistnew(); /* top-level DODS_XXX {...} */
    OClist* dasnodes = oclistnew();
    OClist* varnodes = oclistnew();
    OClist* ddsnodes;
    unsigned int i,j;

    if(dasroot->tree == NULL || dasroot->tree->dxdclass != OCDAS)
	{stat = OCTHROW(OC_EINVAL); goto done;}
    if(ddsroot->tree == NULL || (ddsroot->tree->dxdclass != OCDDS
        && ddsroot->tree->dxdclass != OCDATADDS))
	{stat = OCTHROW(OC_EINVAL); goto done;}

    ddsnodes = ddsroot->tree->nodes;

    /* 1. collect all the relevant DAS nodes;
          namely those that contain at least one
          attribute value, not including the leaf attributes.
          Simultaneously look for potential ambiguities
          if found; complain but continue: result are indeterminate.
          also collect globals separately*/
    for(i=0;i<oclistlength(dasroot->tree->nodes);i++) {
	OCnode* das = (OCnode*)oclistget(dasroot->tree->nodes,i);
	int hasattributes = 0;
	if(das->octype == OC_Attribute) continue; /* ignore these for now*/
	if(das->name == NULL || das->att.isglobal) {
	    oclistpush(dasglobals,(void*)das);
	    continue;
	}
	if(das->att.isdods) {
	    oclistpush(dodsglobals,(void*)das);
	    continue;
	}
	for(j=0;j<oclistlength(das->subnodes);j++) {
	    OCnode* subnode = (OCnode*)oclistget(das->subnodes,j);
	    if(subnode->octype == OC_Attribute) {hasattributes = 1; break;}
	}
	if(hasattributes) {
	    /* Look for previously collected nodes with same name*/
            for(j=0;j<oclistlength(dasnodes);j++) {
	        OCnode* das2 = (OCnode*)oclistget(dasnodes,j);
		if(das->name == NULL || das2->name == NULL) continue;
		if(strcmp(das->name,das2->name)==0) {
		    oclog(OCLOGWARN,"oc_mergedas: potentially ambiguous DAS name: %s",das->name);
		}
	    }
	    oclistpush(dasnodes,(void*)das);
	}
    }

    /* 2. collect all the leaf DDS nodes (of type OC_Atomic)*/
    for(i=0;i<oclistlength(ddsnodes);i++) {
	OCnode* dds = (OCnode*)oclistget(ddsnodes,i);
	if(dds->octype == OC_Atomic) oclistpush(varnodes,(void*)dds);
    }

    /* 3. For each das node, locate matching DDS node(s) and attach
          attributes to the DDS node(s).
          Match means:
          1. DAS->fullname :: DDS->fullname
          2. DAS->name :: DDS->fullname (support DAS names with embedded '.')
          3. DAS->name :: DDS->name
    */
    for(i=0;i<oclistlength(dasnodes);i++) {
	OCnode* das = (OCnode*)oclistget(dasnodes,i);
        for(j=0;j<oclistlength(varnodes);j++) {
	    OCnode* dds = (OCnode*)oclistget(varnodes,j);
	    if(strcmp(das->fullname,dds->fullname)==0
	       || strcmp(das->name,dds->fullname)==0
	       || strcmp(das->name,dds->name)==0) {
		mergedas1(dds,das);
		/* remove from dasnodes list*/
		oclistset(dasnodes,i,(void*)NULL);
	    }
	}
    }

    /* 4. Assign globals*/
    for(i=0;i<oclistlength(dasglobals);i++) {
	OCnode* das = (OCnode*)oclistget(dasglobals,i);
	if(das == NULL) continue;
	mergedas1(ddsroot,das);
    }
    /* 5. Assign DODS_*/
    for(i=0;i<oclistlength(dodsglobals);i++) {
	OCnode* das = (OCnode*)oclistget(dodsglobals,i);
	if(das == NULL) continue;
	mergedods1(ddsroot,das);
    }
    /* 6. Assign other orphan attributes, which means
	  construct their full name and assign as a global attribute. */
    for(i=0;i<oclistlength(dasnodes);i++) {
	OCnode* das = (OCnode*)oclistget(dasnodes,i);
	if(das == NULL) continue;
	mergeother1(ddsroot, das);
    }

done:
    /* cleanup*/
    oclistfree(dasglobals);
    oclistfree(dodsglobals);
    oclistfree(dasnodes);
    oclistfree(varnodes);
    return OCTHROW(stat);
}

static OCerror
mergedas1(OCnode* dds, OCnode* das)
{
    unsigned int i;
    OCerror stat = OC_NOERR;
    if(das == NULL) return OC_NOERR; /* nothing to do */
    if(dds->attributes == NULL) dds->attributes = oclistnew();
    /* assign the simple attributes in the das set to this dds node*/
    for(i=0;i<oclistlength(das->subnodes);i++) {
	OCnode* attnode = (OCnode*)oclistget(das->subnodes,i);
	if(attnode->octype == OC_Attribute) {
	    OCattribute* att = makeattribute(attnode->name,
						attnode->etype,
						attnode->att.values);
            oclistpush(dds->attributes,(void*)att);
	}
    }
    return OCTHROW(stat);
}

static OCerror
mergedods1(OCnode* dds, OCnode* dods)
{
    unsigned int i;
    OCerror stat = OC_NOERR;
    if(dods == NULL) return OC_NOERR; /* nothing to do */
    OCASSERT(dods->octype == OC_Attributeset);
    if(dds->attributes == NULL) dds->attributes = oclistnew();
    /* assign the simple attributes in the das set to this dds node
       with renaming to tag as DODS_
    */
    for(i=0;i<oclistlength(dods->subnodes);i++) {
	OCnode* attnode = (OCnode*)oclistget(dods->subnodes,i);
	if(attnode->octype == OC_Attribute) {
	    OCattribute* att;
	    /* prefix the attribute name with the name of the attribute
               set plus "."
            */
	    size_t len =   strlen(attnode->name)
                         + strlen(dods->name)
			 + strlen(".")
			 + 1; /*null*/
	    char* newname = (char*)malloc(len);
	    if(newname == NULL) return OC_ENOMEM;
	    strcpy(newname,dods->name);
	    strcat(newname,".");
	    strcat(newname,attnode->name);
	    att = makeattribute(newname,attnode->etype,attnode->att.values);
	    free(newname);
            oclistpush(dds->attributes,(void*)att);
	}
    }
    return OCTHROW(stat);
}

static OCerror
mergeother1(OCnode* root, OCnode* das)
{
    OCerror stat = OC_NOERR;
    OCattribute* att = NULL;

    OCASSERT(root != NULL);
    if(root->attributes == NULL) root->attributes = oclistnew();

    if(das->octype == OC_Attribute) {
        /* compute the full name of this attribute */
        computefullname(das);
        /* create attribute */
        att = makeattribute(das->fullname,das->etype,das->att.values);	
        oclistpush(root->attributes,(void*)att);
    } else if(das->octype == OC_Attributeset) {
	int i;
	/* Recurse */
        for(i=0;i<oclistlength(das->subnodes);i++) {
	    OCnode* sub = (OCnode*)oclistget(das->subnodes,i);
	    if(sub == NULL) continue;
	    mergeother1(root,sub);
	}	
    } else
	stat = OC_EDAS;
    return OCTHROW(stat);
}

static void
ocuncorrelate(OCnode* root)
{
    OCtree* tree = root->tree;
    unsigned int i;
    if(tree == NULL) return;
    for(i=0;i<oclistlength(tree->nodes);i++) {
	OCnode* node = (OCnode*)oclistget(tree->nodes,i);
	node->datadds = NULL;
    }        
}

static OCerror
occorrelater(OCnode* dds, OCnode* dxd)
{
    int i,j;
    OCerror ocstat = OC_NOERR;

    if(dds->octype != dxd->octype) {
	OCTHROWCHK((ocstat = OC_EINVAL)); goto fail;
    }
    if(dxd->name != NULL && dxd->name != NULL
       && strcmp(dxd->name,dds->name) != 0) {
	OCTHROWCHK((ocstat = OC_EINVAL)); goto fail;
    } else if(dxd->name != dxd->name) { /* test NULL==NULL */
	OCTHROWCHK((ocstat = OC_EINVAL)); goto fail;
    }

    if(dxd->array.rank != dds->array.rank) {
	OCTHROWCHK((ocstat = OC_EINVAL)); goto fail;
    }

    dds->datadds = dxd;

    switch (dds->octype) {
    case OC_Dataset:
    case OC_Structure:
    case OC_Grid:
    case OC_Sequence:
	/* Remember: there may be fewer datadds fields than dds fields */
	for(i=0;i<oclistlength(dxd->subnodes);i++) {
	    OCnode* dxd1 = (OCnode*)oclistget(dxd->subnodes,(size_t)i);
	    for(j=0;j<oclistlength(dds->subnodes);j++) {
		OCnode* dds1 = (OCnode*)oclistget(dds->subnodes,(size_t)j);
		if(strcmp(dxd1->name,dds1->name) == 0) {
		    ocstat = occorrelater(dds1,dxd1);
		    if(ocstat != OC_NOERR) {OCTHROWCHK(ocstat); goto fail;}
		    break;
		}
	    }
	}
	break;
    case OC_Dimension:
    case OC_Atomic:
	break;
    default: OCPANIC1("unexpected node type: %d",dds->octype);
    }
    /* Correlate the dimensions */
    if(dds->array.rank > 0) {
	for(i=0;i<oclistlength(dxd->subnodes);i++) {
	    OCnode* ddsdim = (OCnode*)oclistget(dds->array.dimensions,(size_t)i);
	    OCnode* dxddim = (OCnode*)oclistget(dxd->array.dimensions,(size_t)i);
	    ocstat = occorrelater(ddsdim,dxddim);
	    if(!ocstat) goto fail;	    
	}	
    }

fail:
    return OCTHROW(ocstat);

}

OCerror
occorrelate(OCnode* dds, OCnode* dxd)
{
    if(dds == NULL || dxd == NULL) return OCTHROW(OC_EINVAL);
    ocuncorrelate(dds);
    return occorrelater(dds,dxd);
}

/*
Mark cacheable those atomic String/URL typed nodes
that are contained only in structures with rank > 0.
*/
void
ocmarkcacheable(OCstate* state, OCnode* ddsroot)
{
    int i,j;
#if 0
    int ok;
#endif
    OClist* treenodes = ddsroot->tree->nodes;
    OClist* path = oclistnew();
    for(i=0;i<oclistlength(treenodes);i++) {
        OCnode* node = (OCnode*)oclistget(treenodes,(size_t)i);
	if(node->octype != OC_Atomic) continue;
	if(node->etype != OC_String && node->etype != OC_URL) continue;
	/* collect node path */
        oclistclear(path);
        occollectpathtonode(node,path);	
#if 0
        ok = 1;
#endif
	for(j=1;j<oclistlength(path)-1;j++) {/* skip top level dataset and node itself*/
            OCnode* pathnode = (OCnode*)oclistget(path,(size_t)j);
	    if(pathnode->octype != OC_Structure
		|| pathnode->array.rank > 0) {
#if 0
	    ok=0;
#endif
	    break;
	    }
	}	
#if 0
	if(ok) {
   	    node->cache.cacheable = 1;
	    node->cache.valid = 0;
	}
#endif
    }
    oclistfree(path);
}

#if 0

OCerror
ocddsdasmerge(OCstate* state, OCnode* ddsroot, OCnode* dasroot)
{
    int i,j;
    OCerror stat = OC_NOERR;
    OClist* globals = oclistnew();
    if(dasroot == NULL) return OCTHROW(stat);
    /* Start by looking for global attributes*/
    for(i=0;i<oclistlength(dasroot->subnodes);i++) {
	OCnode* node = (OCnode*)oclistget(dasroot->subnodes,i);
	if(node->att.isglobal) {
	    for(j=0;j<oclistlength(node->subnodes);j++) {
		OCnode* attnode = (OCnode*)oclistget(node->subnodes,j);
		Attribute* att = makeattribute(attnode->name,
						attnode->etype,
						attnode->att.values);
		oclistpush(globals,(void*)att);
	    }
	}
    }
    ddsroot->attributes = globals;
    /* Now try to match subnode names with attribute set names*/
    for(i=0;i<oclistlength(dasroot->subnodes);i++) {
	OCnode* das = (OCnode*)oclistget(dasroot->subnodes,i);
	int match = 0;
        if(das->att.isglobal) continue;
        if(das->octype == OC_Attributeset) {
            for(j=0;j<oclistlength(ddsroot->subnodes) && !match;j++) {
	        OCnode* dds = (OCnode*)oclistget(ddsroot->subnodes,j);
	        if(strcmp(das->name,dds->name) == 0) {
		    match = 1;
	            stat = mergedas1(dds,das);
	            if(stat != OC_NOERR) break;
		}
	    }
	}
        if(!match) {marklostattribute(das);}
    }
    if(stat == OC_NOERR) ddsroot->attributed = 1;
    return OCTHROW(stat);
}

/* Merge das attributes into the dds node*/

static int
mergedas1(OCnode* dds, OCnode* das)
{
    int i,j;
    int stat = OC_NOERR;
    if(dds->attributes == NULL) dds->attributes = oclistnew();
    /* assign the simple attributes in the das set to this dds node*/
    for(i=0;i<oclistlength(das->subnodes);i++) {
	OCnode* attnode = (OCnode*)oclistget(das->subnodes,i);
	if(attnode->octype == OC_Attribute) {
	    Attribute* att = makeattribute(attnode->name,
						attnode->etype,
						attnode->att.values);
            oclistpush(dds->attributes,(void*)att);
	}
    }
    /* Try to merge any enclosed sets with subnodes of dds*/
    for(i=0;i<oclistlength(das->subnodes);i++) {
	OCnode* dasnode = (OCnode*)oclistget(das->subnodes,i);
	int match = 0;
        if(dasnode->octype == OC_Attribute) continue; /* already dealt with above*/
        for(j=0;j<oclistlength(dds->subnodes) && !match;j++) {
	    OCnode* ddsnode = (OCnode*)oclistget(dds->subnodes,j);
	    if(strcmp(dasnode->name,ddsnode->name) == 0) {
	        match = 1;
	        stat = mergedas1(ddsnode,dasnode);
	        if(stat != OC_NOERR) break;
	    }
	}
        if(!match) {marklostattribute(dasnode);}
    }
    return OCTHROW(stat);
}

void*
oclinearize(OCtype etype, unsigned int nstrings, char** strings)
{
    int i;
    size_t typesize;
    char* memp;
    char* memory;

    if(nstrings == 0) return NULL;
    typesize = octypesize(etype);
    memory = (char*)ocmalloc(nstrings*typesize);
    MEMCHECK(memory,NULL);
    memp = memory;
    for(i=0;i<nstrings;i++) {
	char* value = strings[i];
        converttype(etype,value,memp);
	memp += typesize;
    }
    return memory;
}

static int
converttype(OCtype etype, char* value, char* memory)
{
    long iv;
    unsigned long uiv;
    double dv;
    char c[1];
    int outofrange = 0;
#ifdef HAVE_LONG_LONG_INT
    long long llv;
    unsigned long long ullv;
#endif

    switch (etype) {
    case OC_Char:
	if(sscanf(value,"%c",c) != 1) goto fail;
	*((char*)memory) = c[0];
	break;
    case OC_Byte:
	if(sscanf(value,"%ld",&iv) != 1) goto fail;
        else if(iv > OC_BYTE_MAX || iv < OC_BYTE_MIN) {iv = OC_BYTE_MAX; outofrange = 1;}
	*((signed char*)memory) = (signed char)iv;
	break;
    case OC_UByte:
	if(sscanf(value,"%lu",&uiv) != 1) goto fail;
        else if(uiv > OC_UBYTE_MAX) {uiv = OC_UBYTE_MAX; outofrange = 1;}
	*((unsigned char*)memory) = (unsigned char)uiv;
	break;
    case OC_Int16:
	if(sscanf(value,"%ld",&iv) != 1) goto fail;
        else if(iv > OC_INT16_MAX || iv < OC_INT16_MIN) {iv = OC_INT16_MAX; outofrange = 1;}
	*((signed short*)memory) = (signed short)iv;
	break;
    case OC_UInt16:
	if(sscanf(value,"%lu",&uiv) != 1) goto fail;
        else if(uiv > OC_UINT16_MAX) {uiv = OC_UINT16_MAX; outofrange = 1;}
	*((unsigned short*)memory) = (unsigned short)uiv;
	break;
    case OC_Int32:
	if(sscanf(value,"%ld",&iv) != 1) goto fail;
        else if(iv > OC_INT32_MAX || iv < OC_INT32_MIN) {iv = OC_INT32_MAX; outofrange = 1;}
	*((signed int*)memory) = (signed int)iv;
	break;
    case OC_UInt32:
	if(sscanf(value,"%lu",&uiv) != 1) goto fail;
        else if(uiv > OC_UINT32_MAX) {uiv = OC_UINT32_MAX; outofrange = 1;}
	*((unsigned char*)memory) = (unsigned int)uiv;
	break;
#ifdef HAVE_LONG_LONG_INT
    case OC_Int64:
	if(sscanf(value,"%lld",&llv) != 1) goto fail;
        /*else if(iv > OC_INT64_MAX || iv < OC_INT64_MIN) goto fail;*/
	*((signed long long*)memory) = (signed long long)llv;
	break;
    case OC_UInt64:
	if(sscanf(value,"%llu",&ullv) != 1) goto fail;
	*((unsigned long long*)memory) = (unsigned long long)ullv;
	break;
#endif
    case OC_Float32:
	if(sscanf(value,"%lf",&dv) != 1) goto fail;
	*((float*)memory) = (float)dv;
	break;
    case OC_Float64:
	if(sscanf(value,"%lf",&dv) != 1) goto fail;
	*((double*)memory) = (double)dv;
	break;
    case OC_String: case OC_URL:
	*((char**)memory) = nulldup(value);
	break;
    default:
	goto fail;
    }
    if(outofrange)
        oc_log(LOGWARN,"converttype range failure: %d: %s",etype,value);
    return 1;
fail:
    oc_log(LOGERR,"converttype bad value: %d: %s",etype,value);
    return 0;
}
#endif /*0*/
