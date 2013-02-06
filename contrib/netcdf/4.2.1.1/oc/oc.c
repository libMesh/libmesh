/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#include "config.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "oc.h"
#include "ocinternal.h"
#include "occontent.h"
#include "ocdebug.h"
#include "ocdump.h"
#include "oclog.h"
#include "occlientparams.h"
#include "ochttp.h"

#undef TRACK

/**************************************************/

static int ocinitialized = 0;

/**************************************************/
/* Track legal ids */

#ifdef OC_FASTCONSISTENCY

#define ocverify(object) ((object) != NULL && (*(object) == OCMAGIC)?1:0)

#define ocassign(object) ((OCobject)(object))
#define ocassignall(list)

#else /*!OC_FASTCONSISTENCY*/

static OClist* ocmap = NULL;

static int
ocverify(unsigned long object)
{
    unsigned int i;
    void** map = (void**)oclistcontents(ocmap);
    unsigned int len = oclistlength(ocmap);
#ifdef TRACK
fprintf(stderr,"verify: %lu\n",object); fflush(stderr);
#endif
    if(object > 0) {
        for(i=0;i<len;i++,map++) {
	    if(*map == (void*)object) {
		return 1;
	    }
	}
    }
    fprintf(stderr,"illegal object id: %lu\n",(unsigned long)object);
    fflush(stderr);
    return 0;
}

static OCobject
ocassign(void* object)
{
#ifdef TRACK
fprintf(stderr,"assign: %lu\n",(unsigned long)object); fflush(stderr);
#endif
    oclistpush(ocmap,(ocelem)object);
    return (OCobject)object;
}

static void
ocassignall(OClist* list)
{
    unsigned int i;
    if(list != NULL)
	for(i=0;i<oclistlength(list);i++) {
            void* object = (void*)oclistget(list,i);
#ifdef TRACK
fprintf(stderr,"assign: %lu\n",(unsigned long)object); fflush(stderr);
#endif
            oclistpush(ocmap,(ocelem)object);
	}
}
#endif /*!OC_FASTCONSISTENCY*/

#define OCVERIFYX(T,s,x,r) if(!ocverify(x)) {return (r);}
#define OCVERIFY(T,s,x) OCVERIFYX(T,s,x,OC_EINVAL)

#define OCDEREF(T,s,x) (s)=(T)(x)

/**************************************************/

static int
oc_initialize(void)
{
    int status = OC_NOERR;
#ifndef OC_FASTCONSISTENCY
    ocmap = oclistnew();    
    oclistsetalloc(ocmap,1024);
#endif
    status = ocinternalinitialize();
    ocinitialized = 1;
    return status;
}

/**************************************************/

OCerror
oc_open(const char* url, OCconnection* connp)
{
    OCerror ocerr;
    OCstate* state;
    if(!ocinitialized) oc_initialize();
    ocerr = ocopen(&state,url);
    if(ocerr == OC_NOERR && connp) {
	*connp = (OCconnection)ocassign(state);
    }
    return ocerr;
}

OCerror
oc_close(OCconnection conn)
{
    OCstate* state;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    occlose(state);
    return OC_NOERR;
}

/* Release/reclaim the tree of objects associated with a given root */
OCerror
oc_root_free(OCconnection conn, OCobject root0)
{
    OCnode* root;
    OCVERIFY(OCnode*,root,root0);
    OCDEREF(OCnode*,root,root0);

    ocfreeroot(root);
    return OC_NOERR;
}

/* Return the # of  OCobjects associated with a tree with specified root */
unsigned int
oc_inq_nobjects(OCconnection conn, OCobject root0)
{
    OCnode* root;
    OClist* nodes;
    unsigned int nobjects;
    OCVERIFYX(OCnode*,root,root0,-1);
    OCDEREF(OCnode*,root,root0);

    if(root == NULL) return 0;
    root = root->root;
    if(root == NULL) return 0;
    nodes = root->tree->nodes;
    nobjects = oclistlength(nodes);
    return nobjects;
}

/* Return all the OCobjects associated with a tree with specified root */
OCobject*
oc_inq_objects(OCconnection conn, OCobject root0)
{
    unsigned int i;
    OCnode* root;
    OClist* nodes;
    OCobject* objects = NULL;
    unsigned int nobjects;
    OCVERIFYX(OCnode*,root,root0,OCNULL);
    OCDEREF(OCnode*,root,root0);

    if(root == NULL) return NULL;
    root = root->root;
    if(root == NULL) return NULL;
    nodes = root->tree->nodes;
    nobjects = oclistlength(nodes);
    if(nodes != NULL && nobjects > 0) {
	size_t len = sizeof(OCobject)*(1+nobjects);
        objects = (OCobject*)ocmalloc(len);
	for(i=0;i<oclistlength(nodes);i++) {
	    objects[i] = (OCobject)oclistget(nodes,i);
	}
	objects[nobjects] = OCNULL; /* null terminate */
    }
    return objects;
}

/* Return the text of the DDS or DAS as received from the server */
/* Return NULL if no fetch was ever made. */
const char*
oc_inq_text(OCconnection conn, OCobject root0)
{
    OCnode* root;
    OCVERIFYX(OCnode*,root,root0,NULL);
    OCDEREF(OCnode*,root,root0);

    if(root == NULL) return NULL;
    root = root->root;
    if(root == NULL) return NULL;
    return root->tree->text;
}

OCerror
oc_inq_object(OCconnection conn,
	  OCobject node0,
	  char** namep,
	  OCtype* objecttypep,
	  OCtype* primitivetypep, /* if objecttype == OC_Primitive */
	  OCobject* containerp,
	  unsigned int* rankp,
	  unsigned int* subnodesp,
	  unsigned int* nattrp)
{
    OCnode* node;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(namep) *namep = nulldup(node->name);
    if(objecttypep) *objecttypep = node->octype;
    if(primitivetypep) *primitivetypep = node->etype;
    if(rankp) *rankp = node->array.rank;
    if(containerp) *containerp = (OCobject)node->container;    
    if(subnodesp) *subnodesp = oclistlength(node->subnodes);
    if(nattrp) {
        if(node->octype == OC_Attribute) {
            *nattrp = oclistlength(node->att.values);
        } else {
            *nattrp = oclistlength(node->attributes);
	}
    }
    return OC_NOERR;
}

/* Useful accessor functions */
OCerror
oc_inq_name(OCconnection conn, OCobject node0, char** namep)
{
    OCstate* state;
    OCnode* node;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(state == NULL || node == NULL) return OC_EINVAL;
    if(namep) *namep = nulldup(node->name);
    return OC_NOERR;
}

OCerror
oc_inq_nsubnodes(OCconnection conn, OCobject node0, unsigned int* nsubnodesp)
{
    OCnode* node;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(nsubnodesp) *nsubnodesp = oclistlength(node->subnodes);
    return OC_NOERR;
}

OCerror
oc_inq_primtype(OCconnection conn, OCobject node0, OCtype* typep)
{
    OCnode* node;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(typep) *typep = node->etype;
    return OC_NOERR;
}

/* Alias for oc_inq_class */
OCerror
oc_inq_type(OCconnection conn, OCobject node0, OCtype* typep)
{
    return oc_inq_class(conn,node0,typep);
}

OCerror
oc_inq_class(OCconnection conn, OCobject node0, OCtype* typep)
{
    OCnode* node;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(typep) *typep = node->octype;
    return OC_NOERR;
}

OCerror
oc_inq_rank(OCconnection conn, OCobject node0, unsigned int* rankp)
{
    OCnode* node;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(rankp) *rankp = node->array.rank;
    return OC_NOERR;
}

OCerror
oc_inq_nattr(OCconnection conn, OCobject node0, unsigned int* nattrp)
{
    OCnode* node;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(nattrp) {
        if(node->octype == OC_Attribute) {
            *nattrp = oclistlength(node->att.values);
        } else {
            *nattrp = oclistlength(node->attributes);
	}
    }
    return OC_NOERR;
}

OCerror
oc_inq_root(OCconnection conn, OCobject node0, OCobject* rootp)
{
    OCnode* node;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(rootp) *rootp = (OCobject)node->root;
    return OC_NOERR;
}

OCerror
oc_inq_container(OCconnection conn, OCobject node0, OCobject* containerp)
{
    OCnode* node;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(containerp) *containerp = (OCobject)node->container;
    return OC_NOERR;
}

/* Return the subnode objects, if any, for a given object */
/* Caller must free returned list */
OCerror
oc_inq_subnodes(OCconnection conn, OCobject node0, OCobject** subnodesp)
{
    OCnode* node;
    OCobject* subnodes = NULL;
    unsigned int len;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    len = oclistlength(node->subnodes);
    if(len > 0) {
	unsigned int i;
	subnodes = (OCobject*)occalloc(sizeof(OCobject),len+1);
	for(i=0;i<len;i++) {
	    OCnode* ocnode = (OCnode*)oclistget(node->subnodes,i);
	    subnodes[i] = (OCobject)ocnode;
	}	
        subnodes[len] = OCNULL; /* NULL terminate */
    }
    if(subnodesp) *subnodesp = subnodes;
    return OC_NOERR;
}

OCerror
oc_inq_ith(OCconnection conn,
           OCobject node0, unsigned int index, OCobject* subnodeidp)
{
    OCnode* node;
    OCobject subnodeid = OCNULL;
    unsigned int nsubnodes;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    nsubnodes = oclistlength(node->subnodes);
    if(nsubnodes > 0 &&  index < nsubnodes) {
        subnodeid = (OCobject)oclistget(node->subnodes,index);
    } else
	return OC_EINVAL;
    if(subnodeidp) *subnodeidp = subnodeid;
    return OC_NOERR;
}

/* Return the dimension objects, if any, for a given object */
/* Caller must free returned dimids */
OCerror
oc_inq_dimset(OCconnection conn, OCobject node0, OCobject** dimids)
{
    OCnode* node;
    OCobject* dims = NULL;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(node->array.rank > 0) {
	unsigned int i;
	dims = (OCobject*)occalloc(sizeof(OCobject),node->array.rank+1);
	for(i=0;i<node->array.rank;i++) {
	    OCnode* dim = (OCnode*)oclistget(node->array.dimensions,i);
	    dims[i] = (OCobject)dim;
	}	
	dims[node->array.rank] = OCNULL;
    }
    if(dimids) *dimids = dims;
    return OC_NOERR;
}


OCerror
oc_inq_ithdim(OCconnection conn, OCobject node0, unsigned int index, OCobject* dimidp)
{
    OCnode* node;
    OCobject dimid = OCNULL;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    if(node->array.rank > 0 && index < node->array.rank) {
        dimid = (OCobject)oclistget(node->array.dimensions,index);
    } else
	return OC_EINVAL;
    if(dimidp) *dimidp = dimid;
    return OC_NOERR;
}

OCerror
oc_inq_dim(OCconnection conn, OCobject node0, size_t* sizep, char** namep)
{
    OCnode* dim;
    OCVERIFY(OCnode*,dim,node0);
    OCDEREF(OCnode*,dim,node0);

    if(dim->octype != OC_Dimension) return OC_EINVAL;
    if(sizep) *sizep = dim->dim.declsize;
    if(namep) *namep = nulldup(dim->name);
    return OC_NOERR;
}

/* Obtain info about the ith attribute attached to a given DDS node*/

/* This procedure returns the value as the original DAS string */
OCerror
oc_inq_attrstrings(OCconnection conn, OCobject node0, unsigned int i,
			   char** namep, OCtype* octypep,
			   unsigned int* nvaluesp, char*** stringsp)
{
    OCnode* node;
    OCattribute* attr;
    unsigned int nattrs;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    nattrs = oclistlength(node->attributes);
    if(i >= nattrs) return OC_EINVAL;
    attr = (OCattribute*)oclistget(node->attributes,i);
    if(namep) *namep = strdup(attr->name);
    if(octypep) *octypep = attr->etype;
    if(nvaluesp) *nvaluesp = attr->nvalues;
    if(stringsp) {
	if(attr->nvalues > 0) {
	    size_t space = attr->nvalues * sizeof(char*);
	    char** strings = ocmalloc(space);
	    for(i=0;i<attr->nvalues;i++)
	        strings[i] = nulldup(attr->values[i]);
	    *stringsp = strings;
	} else
	    *stringsp = NULL;
    }
    return OC_NOERR;    
}

/* This procedure returns the value as the default binary value 
   corresponding to the string value
*/

OCerror
oc_inq_attr(OCconnection conn, OCobject node0, unsigned int i,
	    char** namep, OCtype* octypep, unsigned int* nvaluesp, void** valuesp)
{
    OCnode* node;
    OCattribute* attr;
    unsigned int nattrs;
    OCVERIFY(OCnode*,node,node0);
    OCDEREF(OCnode*,node,node0);

    nattrs = oclistlength(node->attributes);
    if(i >= nattrs) return OC_EINVAL;
    attr = (OCattribute*)oclistget(node->attributes,i);
    if(namep) *namep = strdup(attr->name);
    if(octypep) *octypep = attr->etype;
    if(nvaluesp) *nvaluesp = attr->nvalues;
    if(valuesp && attr->nvalues > 0) {
	void* memory = NULL;
        memory = oclinearize(attr->etype,attr->nvalues,attr->values);
	*valuesp = memory;
    }
    return OC_NOERR;    
}

/* Convenience function */
void
oc_attr_reclaim(OCtype etype, unsigned int nvalues, void* values)
{
    if(nvalues == 0 || values == NULL) return;
    if(etype == OC_String || etype == OC_URL) {
        unsigned int i;
	char** strings = (char**)values;
	for(i=0;i<nvalues;i++) {ocfree(strings[i]);}	
    }
    ocfree(values);    
}

OCerror
oc_inq_dasattr_nvalues(OCconnection conn, OCobject node0,
			unsigned int* nvaluesp)
{
    OCnode* attr;
    OCVERIFY(OCnode*,attr,node0);
    OCDEREF(OCnode*,attr,node0);
    if(attr->octype != OC_Attribute) return OC_EINVAL;
    if(nvaluesp) *nvaluesp = oclistlength(attr->att.values);
    return OC_NOERR;
}

OCerror
oc_inq_dasattr(OCconnection conn, OCobject node0, unsigned int i,
               OCtype* primtypep, char** valuep)
{
    OCnode* attr;
    unsigned int nvalues;
    OCVERIFY(OCnode*,attr,node0);
    OCDEREF(OCnode*,attr,node0);

    if(attr->octype != OC_Attribute) return OC_EINVAL;
    nvalues = oclistlength(attr->att.values);
    if(i >= nvalues) return OC_EINVAL;
    if(valuep) *valuep = nulldup((char*)oclistget(attr->att.values,i));
    if(primtypep) *primtypep = attr->etype;
    return OC_NOERR;
}

/**************************************************/
/* Fetch and parse a given class of DXD the server specified
   at open time, and using a specified set of constraints.
   Return the root node of the parsed tree of objects.
*/
OCerror oc_fetch(OCconnection conn, const char* constraint,
                 OCdxd dxdkind, OCobject* rootp)
{
    return oc_fetchf(conn,constraint,dxdkind,0,rootp);
}

OCerror oc_fetchf(OCconnection conn, const char* constraint,
                 OCdxd dxdkind, OCflags flags, OCobject* rootp)
{
    OCstate* state;
    OCerror ocerr = OC_NOERR;
    OCnode* root;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);

    ocerr = ocfetchf(state,constraint,dxdkind,flags,&root);
    if(ocerr) return ocerr;

    ocassignall(root->tree->nodes);
    if(rootp) *rootp = (OCobject)ocassign(root);
    return ocerr;
}


OCerror
oc_data_root(OCconnection conn, OCobject root0, OCdata content0)
{
    OCstate* state;
    OCnode* root;
    OCcontent* content;
    OCerror ocerr = OC_NOERR;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCnode*,root,root0);
    OCDEREF(OCnode*,root,root0);
    OCVERIFY(OCcontent*,content,content0);
    OCDEREF(OCcontent*,content,content0);

    if(root->tree == NULL) {OCTHROWCHK((ocerr=OC_EINVAL)); goto fail;}
    ocerr = ocrootdata(state,root,content);

fail:
    return ocerr;
}

OCdata
oc_data_new(OCconnection conn)
{
    OCstate* state;
    OCVERIFYX(OCstate*,state,conn,OCNULL);
    OCDEREF(OCstate*,state,conn);

    return (OCdata)ocassign(ocnewcontent(state));
}

OCerror
oc_data_free(OCconnection conn, OCdata content0)
{
    OCstate* state;
    OCcontent* content;
    if(content0 == OCNULL) return OC_NOERR;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCcontent*,content,content0);
    OCDEREF(OCcontent*,content,content0);

    ocfreecontent(state,content);
    return OC_NOERR;
}

OCerror
oc_data_ith(OCconnection conn, OCdata parentdata, size_t index, OCdata subdata)
{
    OCstate* state;
    OCcontent* parent;
    OCcontent* child;
    OCerror ocerr = OC_NOERR;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCcontent*,parent,parentdata);
    OCDEREF(OCcontent*,parent,parentdata);
    OCVERIFY(OCcontent*,child,subdata);
    OCDEREF(OCcontent*,child,subdata);
    ocerr = ocdataith(state,parent,index,child);
    return ocerr;
}

OCerror
oc_data_get(OCconnection conn, OCdata currentcontent,
            void* memory, size_t memsize, size_t start, size_t count)
{
    OCstate* state;
    OCcontent* current;
    OCerror ocerr = OC_NOERR;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCcontent*,current,currentcontent);
    OCDEREF(OCcontent*,current,currentcontent);

    ocerr = ocgetcontent(state,current,memory,memsize,start,count);

    if(ocerr == OC_EDATADDS) ocdataddsmsg(state,current->tree);

    return ocerr;
}

OCerror
oc_data_count(OCconnection conn, OCdata content0, size_t* sizep)
{
    OCstate* state;
    OCcontent* current;
    OCerror ocerr = OC_NOERR;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCcontent*,current,content0);
    OCDEREF(OCcontent*,current,content0);
    ocerr = ocdatacount(state, current, sizep);
    return ocerr;
}

OCerror
oc_data_index(OCconnection conn, OCdata content0, size_t* sizep)
{
    OCcontent* current;
    OCerror ocerr = OC_NOERR;
    OCVERIFY(OCcontent*,current,content0);
    OCDEREF(OCcontent*,current,content0);

    if(sizep)
	*sizep = (current->cache.valid ? current->cache.index : 0);
    return ocerr;
}

OCerror
oc_data_object(OCconnection conn, OCdata content0, OCobject* op)
{
    OCcontent* current;
    OCerror ocerr = OC_NOERR;
    OCVERIFY(OCcontent*,current,content0);
    OCDEREF(OCcontent*,current,content0);

    if(op) *op = (OCobject)current->node;
    return ocerr;
}

OCerror
oc_data_mode(OCconnection conn, OCdata content0, OCmode* modep)
{
    OCcontent* current;
    OCerror ocerr = OC_NOERR;
    OCVERIFY(OCcontent*,current,content0);
    OCDEREF(OCcontent*,current,content0);

    if(modep) *modep = current->mode;
    return ocerr;
}

/**************************************************/
/* OCtype management */
size_t oc_typesize(OCtype etype)
{
    return octypesize(etype);
}

const char*
oc_typetostring(OCtype octype)
{
    return octypetoddsstring(octype);
}

OCerror
oc_typeprint(OCtype etype, char* buf, size_t bufsize, void* value)
{
    return octypeprint(etype,buf,bufsize,value);
}

/**************************************************/
/* The oc_logXXX procedures are define in oclog.c */

/**************************************************/
/* Miscellaneous */

const char*
oc_errstring(int err)
{
    return ocerrstring(err);
}

/* Get clientparameters from the URL */
const char*
oc_clientparam_get(OCconnection conn, const char* param)
{
    OCstate* state;
    OCVERIFYX(OCstate*,state,conn,NULL);
    OCDEREF(OCstate*,state,conn);

    return ocparamlookup(state,param);
}

#ifdef OCIGNORE
/* Delete client parameter
   return value:
	OC_NOERR => defined; deletion performed
	OC_EINVAL => not already defined
*/
OCerror
oc_clientparam_delete(OCconnection conn, const char* param)
{
    OCstate* state;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);

    return ocparamdelete(state->clientparams,param);
}

/* Insert client parameter
   return value:
	OC_NOERR => not already define; insertion performed
	OC_EINVAL => already defined
*/
OCerror
oc_clientparam_insert(OCconnection conn, const char* param, const char* value)
{
    OCstate* state;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);

    state->clientparams = dapparaminsert(state->clientparams,param,value);
    return OC_NOERR;
}

/* Replace client parameter
   return value:
	OC_NOERR => already define; replacement performed
	OC_EINVAL => not already defined
*/
OCerror
oc_clientparam_replace(OCconnection conn, const char* param, const char* value)
{
    OCstate* state;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);

    return dapparamreplace(state->clientparams,param,value);
}
#endif

OCerror
oc_dd(OCconnection conn, OCobject root0, int level)
{
    OCstate* state;
    OCnode* root;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCnode*,root,root0);
    OCDEREF(OCnode*,root,root0);

    ocdd(state,root,1,level);
    return OC_NOERR;
}

OCerror
oc_ddnode(OCconnection conn, OCobject root0)
{
    OCnode* root;
    OCVERIFY(OCnode*,root,root0);
    OCDEREF(OCnode*,root,root0);

    ocdumpnode(root);
    return OC_NOERR;
}

/* Merge a specified DAS into a specified DDS or DATADDS */
OCerror
oc_attach_das(OCconnection conn, OCobject dasroot, OCobject ddsroot)
{
    OCstate* state;
    OCnode* das;
    OCnode* dds;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCnode*,das,dasroot);
    OCDEREF(OCnode*,das,dasroot);
    OCVERIFY(OCnode*,dds,ddsroot);
    OCDEREF(OCnode*,dds,ddsroot);

    return ocddsdasmerge(state,das,dds);
}

#if 0
/*
I suppressed this operation because I realized that it
appears to be impossible to implement correctly in general.
It also exposes a flaw in the protocol.
If I have a dds with two identical Grids, G1 and G2,
and I ask for the projection ?G1.temp,G2.temp,
then I get back a DATADDS with duplicated temp fields,
which means the DATADDS is illegal.  The problem is that
the protocol throws away important scoping information
about the fact that each temp field is part of a different
grid.
*/
/* Connect a specified DATADDS tree to a specified DDS tree */
OCerror
oc_attach_datadds(OCconnection conn, OCobject dataddsroot, OCobject ddsroot)
{
    OCstate* state;
    OCnode* dxd;
    OCnode* dds;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCnode*,dxd,dataddsroot);
    OCDEREF(OCnode*,dxd,dataddsroot);
    OCVERIFY(OCnode*,dds,ddsroot);
    OCDEREF(OCnode*,dds,ddsroot);

    /* get the true roots */
    dds = dds->root;
    dxd = dxd->root;
    if(dds == NULL || dxd == NULL) return OC_EBADID;    
    /* correlate the DATADDS to the DDS */
    return occorrelate(dxd,dds);
}

/* Return the attached DATADDS object for a given DDS object */
OCerror oc_inq_datadds(OCconnection conn, OCobject dds0, OCobject* dataddsp)
{
    OCstate* state;
    OCnode* dds;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    OCVERIFY(OCnode*,dds,dds0);
    OCDEREF(OCnode*,dds,dds0);

    if(dataddsp) *dataddsp = (OCobject)dds->datadds;
    return OC_NOERR;
}
#endif

/**************************************************/

OCerror
oc_svcerrordata(OCconnection conn, char** codep,
                               char** msgp, long* httpp)
{
    OCstate* state;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    return ocsvcerrordata(state,codep,msgp,httpp);
}


/**************************************************/
/* Experimental: this is useful for the netcdf
   DRNO project.
*/

/* New 10/31/2009: return the size (in bytes)
   of the fetched datadds.
*/

OCerror
oc_raw_xdrsize(OCconnection conn, OCobject root0, size_t* sizep)
{
    OCnode* root;
    OCerror ocerr = OC_NOERR;
    OCVERIFY(OCnode*,root,root0);
    OCDEREF(OCnode*,root,root0);

    if(sizep == NULL) goto done;
    if(root->tree == NULL || root->tree->dxdclass != OCDATADDS)
	{OCTHROWCHK((ocerr=OC_EINVAL)); goto done;}
    if(sizep) *sizep = root->tree->data.datasize;

done:
    return ocerr;
}

/* Resend a url as a head request to check the Last-Modified time */
OCerror
oc_update_lastmodified_data(OCconnection conn)
{
    OCstate* state;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    return ocupdatelastmodifieddata(state);
}

long
oc_get_lastmodified_data(OCconnection conn)
{
    OCstate* state;
    OCVERIFY(OCstate*,state,conn);
    OCDEREF(OCstate*,state,conn);
    return state->datalastmodified;
}

int
oc_dumpnode(OCconnection conn, OCobject root0)
{
    OCnode* root;
    OCerror ocerr = OC_NOERR;
    OCVERIFY(OCnode*,root,root0);
    OCDEREF(OCnode*,root,root0);
    ocdumpnode(root);
    return ocerr;
}

OCerror
oc_ping(const char* url)
{
    return ocping(url);
}
    
