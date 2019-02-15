/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "d4includes.h"
#include <stdarg.h>
#include <assert.h>
#include "ezxml.h"

/**
 * Implement the Dap4 Parser Using a DOM Parser
 *
 * This code creates in internal representation of the netcdf-4 metadata
 * to avoid having to make so many calls into the netcdf library.
 */

/***************************************************/

/*
Map an xml node name to an interpretation.
We use NCD4_NULL when we directly interpret
the tag and do not need to search for it.
E.g Dataset, Dim, econst, etc.
If the sort is NCD4_NULL, then that means it is
irrelevant because that keyword will never be
searched for in this table.
*/
static const struct KEYWORDINFO {
    char* tag; /* The xml tag e.g. <tag...> */
    NCD4sort sort; /* What kind of node are we building */
    nc_type subsort; /* discriminator */
    char* aliasfor; /* Some names are aliases for others */
} keywordmap[] = {
{"Attribute", NCD4_ATTR,NC_NAT,NULL},
{"Byte", NCD4_VAR,NC_BYTE,"Int8"},
{"Char", NCD4_VAR,NC_CHAR,NULL},
{"Dataset", NCD4_NULL,NC_NAT,NULL},
{"Dim", NCD4_NULL,NC_NAT,NULL},
{"Dimension", NCD4_DIM,NC_NAT,NULL},
{"Enum", NCD4_VAR,NC_ENUM,NULL},
{"Enumconst", NCD4_NULL,NC_NAT,NULL},
{"Enumeration", NCD4_TYPE,NC_ENUM,NULL},
{"Float32", NCD4_VAR,NC_FLOAT,NULL},
{"Float64", NCD4_VAR,NC_DOUBLE,NULL},
{"Group", NCD4_GROUP,NC_NAT,NULL},
{"Int16", NCD4_VAR,NC_SHORT,NULL},
{"Int32", NCD4_VAR,NC_INT,NULL},
{"Int64", NCD4_VAR,NC_INT64,NULL},
{"Int8", NCD4_VAR,NC_BYTE,NULL},
{"Map", NCD4_NULL,NC_NAT,NULL},
{"Opaque", NCD4_VAR,NC_OPAQUE,NULL},
{"OtherXML", NCD4_XML,NC_NAT,NULL},
{"Sequence", NCD4_VAR,NC_SEQ,NULL},
{"String", NCD4_VAR,NC_STRING,NULL},
{"Structure", NCD4_VAR,NC_STRUCT,NULL},
{"UByte", NCD4_VAR,NC_UBYTE,"UInt8"},
{"UInt16", NCD4_VAR,NC_USHORT,NULL},
{"UInt32", NCD4_VAR,NC_UINT,NULL},
{"UInt64", NCD4_VAR,NC_UINT64,NULL},
{"UInt8", NCD4_VAR,NC_UBYTE,NULL},
{"URL", NCD4_VAR,NC_STRING,"String"},
};
typedef struct KEYWORDINFO KEYWORDINFO;

static const struct ATOMICTYPEINFO {
    char* name; nc_type type; size_t size;
} atomictypeinfo[] = {
/* Keep in sorted order for binary search */
{"Byte",NC_BYTE,sizeof(char)},
{"Char",NC_CHAR,sizeof(char)},
{"Float32",NC_FLOAT,sizeof(float)},
{"Float64",NC_DOUBLE,sizeof(double)},
{"Int16",NC_SHORT,sizeof(short)},
{"Int32",NC_INT,sizeof(int)},
{"Int64",NC_INT64,sizeof(long long)},
{"Int8",NC_BYTE,sizeof(char)},
{"String",NC_STRING,sizeof(char*)},
{"UByte",NC_UBYTE,sizeof(unsigned char)},
{"UInt16",NC_USHORT,sizeof(unsigned short)},
{"UInt32",NC_UINT,sizeof(unsigned int)},
{"UInt64",NC_UINT64,sizeof(unsigned long long)},
{"UInt8",NC_UBYTE,sizeof(unsigned char)},
{NULL,NC_NAT,0}
};

/***************************************************/

#ifdef D4DEBUG
static void setname(NCD4node* node, const char* name)
{
    nullfree(node->name);
    (node)->name = strdup(name); \
    fprintf(stderr,"setname: node=%lx name=%s\n",(unsigned long)(node),(node)->name); \
}
#define SETNAME(node,src) setname((node),src)
#else
#define SETNAME(node,src) do {nullfree((node)->name); (node)->name = strdup(src);} while(0);
#endif

/***************************************************/

extern const char** ezxml_all_attr(ezxml_t xml, int* countp);

/* Forwards */

static int addOrigType(NCD4parser*, NCD4node* src, NCD4node* dst, const char* tag);
static int defineAtomicTypes(NCD4parser*);
static void classify(NCD4node* container, NCD4node* node);
static int convertString(union ATOMICS*, NCD4node* type, const char* s);
static int downConvert(union ATOMICS*, NCD4node* type);
static int fillgroup(NCD4parser*, NCD4node* group, ezxml_t xml);
static NCD4node* getOpaque(NCD4parser*, ezxml_t varxml, NCD4node* group);
static int getValueStrings(NCD4parser*, NCD4node*, ezxml_t xattr, NClist*);
static int isReserved(const char* name);
static const KEYWORDINFO* keyword(const char* name);
static NCD4node* lookupAtomictype(NCD4parser*, const char* name);
static NCD4node* lookFor(NClist* elems, const char* name, NCD4sort sort);
static NCD4node* lookupFQN(NCD4parser*, const char* sfqn, NCD4sort);
static int lookupFQNList(NCD4parser*, NClist* fqn, NCD4sort sort, NCD4node** result);
static NCD4node* makeAnonDim(NCD4parser*, const char* sizestr);
static int makeNode(NCD4parser*, NCD4node* parent, ezxml_t, NCD4sort, nc_type, NCD4node**);
static int parseAtomicVar(NCD4parser*, NCD4node* container, ezxml_t xml, NCD4node**);
static int parseAttributes(NCD4parser*, NCD4node* container, ezxml_t xml);
static int parseDimensions(NCD4parser*, NCD4node* group, ezxml_t xml);
static int parseDimRefs(NCD4parser*, NCD4node* var, ezxml_t xml);
static int parseEconsts(NCD4parser*, NCD4node* en, ezxml_t xml);
static int parseEnumerations(NCD4parser*, NCD4node* group, ezxml_t dom);
static int parseFields(NCD4parser*, NCD4node* container, ezxml_t xml);
static int parseError(NCD4parser*, ezxml_t errxml);
static int parseGroups(NCD4parser*, NCD4node* group, ezxml_t dom);
static int parseMaps(NCD4parser*, NCD4node* var, ezxml_t xml);
static int parseMetaData(NCD4parser*, NCD4node* node, ezxml_t xml);
static int parseStructure(NCD4parser*, NCD4node* container, ezxml_t dom, NCD4node**);
static int parseSequence(NCD4parser*, NCD4node* container, ezxml_t dom,NCD4node**);
static int parseLL(const char* text, long long*);
static int parseULL(const char* text, unsigned long long*);
static int parseVariables(NCD4parser*, NCD4node* group, ezxml_t xml);
static int parseVariable(NCD4parser*, NCD4node* group, ezxml_t xml, NCD4node**);
static void reclaimParser(NCD4parser* parser);
static void record(NCD4parser*, NCD4node* node);
static int splitOrigType(NCD4parser*, const char* fqn, NCD4node* var);
static void track(NCD4parser*, NCD4node* node);
static int traverse(NCD4parser*, ezxml_t dom);
#ifndef FIXEDOPAQUE
static int defineBytestringType(NCD4parser*);
#endif

/***************************************************/
/* API */

int
NCD4_parse(NCD4meta* metadata)
{
    int ret = NC_NOERR;
    NCD4parser* parser = NULL;
    int ilen;
    ezxml_t dom = NULL;

    /* Create and fill in the parser state */
    parser = (NCD4parser*)calloc(1,sizeof(NCD4parser));
    if(parser == NULL) {ret=NC_ENOMEM; goto done;}
    parser->metadata = metadata;
    ilen = strlen(parser->metadata->serial.dmr);
    dom = ezxml_parse_str(parser->metadata->serial.dmr,ilen);
    if(dom == NULL) {ret=NC_ENOMEM; goto done;}
    parser->types = nclistnew();
    parser->dims = nclistnew();
    parser->vars = nclistnew();
#ifdef D4DEBUG
    parser->debuglevel = 1;
#endif

    /*Walk the DOM tree */
    ret = traverse(parser,dom);

done:
    if(dom != NULL)
	ezxml_free(dom);
    reclaimParser(parser);
    return THROW(ret);
}

static void
reclaimParser(NCD4parser* parser)
{
    int i,len;
    if(parser == NULL) return;
    nclistfree(parser->types);
    nclistfree(parser->dims);
    nclistfree(parser->vars);
    /* Reclaim unused atomic type nodes */
    len = nclistlength(parser->atomictypes);    
    for(i=0;i<len;i++) {
	if(parser->used[i])
	    reclaimNode((NCD4node*)nclistget(parser->atomictypes,i));
    }
    nclistfree(parser->atomictypes);
    nullfree(parser->used);
    free (parser);
}

/**************************************************/

/* Recursively walk the DOM tree to create the metadata */
static int
traverse(NCD4parser* parser, ezxml_t dom)
{
    int ret = NC_NOERR;

    /* See if we have an <Error> or <Dataset> */
    if(strcmp(dom->name,"Error")==0) {
	ret=parseError(parser,dom);
	ret=NC_EDMR;
	goto done;
    } else if(strcmp(dom->name,"Dataset")==0) {
	const char* xattr = NULL;
        if((ret=makeNode(parser,NULL,NULL,NCD4_GROUP,NC_NULL,&parser->metadata->root))) goto done;
        parser->metadata->root->group.isdataset = 1;
        parser->metadata->root->meta.id = parser->metadata->ncid;
        parser->metadata->groupbyid = nclistnew();
        SETNAME(parser->metadata->root,"/");
	xattr = ezxml_attr(dom,"name");
	if(xattr != NULL) parser->metadata->root->group.datasetname = strdup(xattr);
	xattr = ezxml_attr(dom,"dapVersion");
	if(xattr != NULL) parser->metadata->root->group.dapversion = strdup(xattr);
	xattr = ezxml_attr(dom,"dmrVersion");
	if(xattr != NULL) parser->metadata->root->group.dmrversion = strdup(xattr);
        /* fill in the atomic types */
        if((ret=defineAtomicTypes(parser))) goto done;
        /* Recursively walk the tree */
        if((ret = fillgroup(parser,parser->metadata->root,dom))) goto done;
    } else
	FAIL(NC_EINVAL,"Unexpected dom root name: %s",dom->name);
done:
    return THROW(ret);
}

static int
fillgroup(NCD4parser* parser, NCD4node* group, ezxml_t xml)
{
    int ret = NC_NOERR;

    /* Extract Dimensions */
    if((ret = parseDimensions(parser,group,xml))) goto done;
    /* Extract Enum types */
    if((ret = parseEnumerations(parser,group,xml))) goto done;
    /* Extract variables */
    if((ret = parseVariables(parser,group,xml))) goto done;
    /* Extract subgroups*/
    if((ret = parseGroups(parser,group,xml))) goto done;
    /* Parse group level attributes */
    if((ret = parseAttributes(parser,group,xml))) goto done;
done:
    return THROW(ret);
}

static int
parseDimensions(NCD4parser* parser, NCD4node* group, ezxml_t xml)
{
    int ret = NC_NOERR;
    ezxml_t x;
    for(x=ezxml_child(xml, "Dimension");x != NULL;x = ezxml_next(x)) {
	NCD4node* dimnode = NULL;
	unsigned long long size;
	const char* sizestr;
	const char* unlimstr;
	sizestr = ezxml_attr(x,"size");
	if(sizestr == NULL)
	    FAIL(NC_EDIMSIZE,"Dimension has no size");
	unlimstr = ezxml_attr(x,UCARTAGUNLIM);
	if((ret = parseULL(sizestr,&size))) goto done;
	if((ret=makeNode(parser,group,x,NCD4_DIM,NC_NULL,&dimnode))) goto done;
	dimnode->dim.size = (long long)size;
	dimnode->dim.isunlimited = (unlimstr != NULL);
	/* Process attributes */
	if((ret = parseAttributes(parser,dimnode,x))) goto done;
	classify(group,dimnode);
    }
done:
    return THROW(ret);
}

static int
parseEnumerations(NCD4parser* parser, NCD4node* group, ezxml_t xml)
{
    int ret = NC_NOERR;
    ezxml_t x;

    for(x=ezxml_child(xml, "Enumeration");x != NULL;x = ezxml_next(x)) {
	NCD4node* node = NULL;
	NCD4node* basetype = NULL;
	const char* fqn = ezxml_attr(x,"basetype");
	basetype = lookupFQN(parser,fqn,NCD4_TYPE);
	if(basetype == NULL) {
	    FAIL(NC_EBADTYPE,"Enumeration has unknown type: ",fqn);
	}
	if((ret=makeNode(parser,group,x,NCD4_TYPE,NC_ENUM,&node))) goto done;
	node->basetype = basetype;
	if((ret=parseEconsts(parser,node,x))) goto done;
	if(nclistlength(node->en.econsts) == 0)
	    FAIL(NC_EINVAL,"Enumeration has no values");
	classify(group,node);
	/* Finally, see if this type has UCARTAGORIGTYPE xml attribute */
	if(parser->metadata->controller->controls.translation == NCD4_TRANSNC4) {
	    const char* typetag = ezxml_attr(x,UCARTAGORIGTYPE);
	    if(typetag != NULL) {
	    }
	}
    }
done:
    return THROW(ret);
}

static int
parseEconsts(NCD4parser* parser, NCD4node* en, ezxml_t xml)
{
    int ret = NC_NOERR;
    ezxml_t x;
    NClist* econsts = nclistnew();

    for(x=ezxml_child(xml, "EnumConst");x != NULL;x = ezxml_next(x)) {
        NCD4node* ec = NULL;
	const char* name;
	const char* svalue;
	name = ezxml_attr(x,"name");
	if(name == NULL) FAIL(NC_EBADNAME,"Enum const with no name");
	if((ret=makeNode(parser,en,x,NCD4_ECONST,NC_NULL,&ec))) goto done	;
	svalue = ezxml_attr(x,"value");
	if(svalue == NULL)
	    FAIL(NC_EINVAL,"Enumeration Constant has no value");
	if((ret=convertString(&ec->en.ecvalue,en->basetype,svalue)))
	    FAIL(NC_EINVAL,"Non-numeric Enumeration Constant: %s->%s",ec->name,svalue);
	PUSH(econsts,ec);
    }
    en->en.econsts = econsts;
done:
    return THROW(ret);
}

static int
parseVariables(NCD4parser* parser, NCD4node* group, ezxml_t xml)
{
    int ret = NC_NOERR;
    ezxml_t x;
    for(x=xml->child;x != NULL;x=x->ordered) {
	NCD4node* node = NULL;
	const KEYWORDINFO* info = keyword(x->name);
	if(info == NULL)
	    FAIL(NC_ETRANSLATION,"Unexpected node type: %s",x->name);
	/* Check if we need to process this node */
	if(!ISVAR(info->sort)) continue; /* Handle elsewhere */
	node = NULL;
	ret = parseVariable(parser,group,x,&node);
        if(ret != NC_NOERR || node == NULL) break;
    }
done:
    return THROW(ret);
}

static int
parseVariable(NCD4parser* parser, NCD4node* container, ezxml_t xml, NCD4node** nodep)
{
    int ret = NC_NOERR;
    NCD4node* node = NULL;
    const KEYWORDINFO* info = keyword(xml->name);

    switch (info->subsort) {
    case NC_STRUCT:
	ret = parseStructure(parser,container,xml,&node);
	break;
    case NC_SEQ:
	ret = parseSequence(parser,container,xml,&node);
	break;
    default:
	ret = parseAtomicVar(parser,container,xml,&node);
    }
    *nodep = node;

    return THROW(ret);
}

static int
parseMetaData(NCD4parser* parser, NCD4node* container, ezxml_t xml)
{
    int ret = NC_NOERR;
    /* Process dimrefs */
    if((ret=parseDimRefs(parser,container,xml))) goto done;
    /* Process attributes */
    if((ret = parseAttributes(parser,container,xml))) goto done;
    /* Process maps */
    if((ret = parseMaps(parser,container,xml))) goto done;
done:
    return THROW(ret);
}

static int
parseStructure(NCD4parser* parser, NCD4node* container, ezxml_t xml, NCD4node** nodep)
{
    int ret = NC_NOERR;
    NCD4node* var = NULL;
    NCD4node* type = NULL;
    NCD4node* group = NULL;
    char* fqnname = NULL;

    group = NCD4_groupFor(container); /* default: put type in the same group as var */

    /* Make the structure as a variable with same name as structure; will be fixed later */
    if((ret=makeNode(parser,container,xml,NCD4_VAR,NC_STRUCT,&var))) goto done;
    classify(container,var);

    /* Make the structure as a type with (for now) partial fqn name from the variable */
    if((ret=makeNode(parser,group,xml,NCD4_TYPE,NC_STRUCT,&type))) goto done;
    classify(group,type);
    /* Set the basetype */
    var->basetype = type;
    /* Now change the struct typename */
    fqnname = NCD4_makeName(var,"_");
    if(fqnname == NULL)
	FAIL(NC_ENOMEM,"Out of memory");
    SETNAME(type,fqnname);

    /* Parse Fields into the type */
    if((ret = parseFields(parser,type,xml))) goto done;

    /* Parse attributes, dims, and maps into the var */
    if((ret = parseMetaData(parser,var,xml))) goto done;

    record(parser,var);

    /* See if this var has UCARTAGORIGTYPE attribute */
    if(parser->metadata->controller->controls.translation == NCD4_TRANSNC4) {
	const char* typetag = ezxml_attr(xml,UCARTAGORIGTYPE);
	if(typetag != NULL) {
	    /* yes, place it on the type */
	    if((ret=addOrigType(parser,var,type,typetag))) goto done;
 	}
    }

    if(nodep) *nodep = var;

done:
    nullfree(fqnname);
    return THROW(ret);
}

static int
parseFields(NCD4parser* parser, NCD4node* container, ezxml_t xml)
{
    int ret = NC_NOERR;
    ezxml_t x;
    for(x=xml->child;x != NULL;x=x->ordered) {
	NCD4node* node = NULL;
        const KEYWORDINFO* info = keyword(x->name);
	if(!ISVAR(info->sort)) continue; /* not a field */
	ret = parseVariable(parser,container,x,&node);
	if(ret) goto done;
    }
done:
    return THROW(ret);
}

/*
Specialized version of parseFields that is used
to attach a singleton field to a vlentype
*/
static int
parseVlenField(NCD4parser* parser, NCD4node* container, ezxml_t xml, NCD4node** fieldp)
{
    int ret = NC_NOERR;
    NCD4node* field = NULL;
    ezxml_t x;
    for(x=xml->child;x != NULL;x=x->ordered) {
        const KEYWORDINFO* info = keyword(x->name);
	if(!ISVAR(info->sort)) continue; /* not a field */
	if(field != NULL)
	    {ret = NC_EBADTYPE; goto done;}
	if((ret = parseVariable(parser,container,x,&field)))
	    goto done;
    }
    if(fieldp) *fieldp = field;
done:
    return THROW(ret);
}

static int
parseSequence(NCD4parser* parser, NCD4node* container, ezxml_t xml, NCD4node** nodep)
{
    int ret = NC_NOERR;
    NCD4node* var = NULL;
    NCD4node* structtype = NULL;
    NCD4node* vlentype = NULL;
    NCD4node* group = NULL;
    char name[NC_MAX_NAME];
    char* fqnname = NULL;
    int usevlen = 0;

    group = NCD4_groupFor(container);

    /* Convert a sequence variable into two or three things:
	1. a compound type representing the fields of the sequence.
	2. a vlen type whose basetype is #1
	3. a variable whose basetype is #2.
	If we can infer that the sequence was riginally produced
	from a netcdf-4 vlen, then we can avoid createing #1.
	Naming is as follows. Assume the var name is V
	and the NCD4_makeName of the var is V1..._Vn.
	1. var name is V.
	2. vlen type is V1..._VN_t
	3. compound type (if any) is V1..._VN_cmpd (Note, d4meta will append _t to this)
     */

    /* Determine if we need to build a structure type or can go straight to a vlen
       Test:  UCARTAGVLEN xml attribute is set
    */
    if(parser->metadata->controller->controls.translation == NCD4_TRANSNC4) {
	const char* vlentag = ezxml_attr(xml,UCARTAGVLEN);
	if(vlentag != NULL)
	    usevlen = 1;
    } else
	usevlen = 0;

    /* make secondary names from the var fqn name */
    if(usevlen) {
	/* Parse the singleton field and then use it to fix up the var */
	if((ret=parseVlenField(parser,container,xml,&var)))
	    goto done;
	/* compute a partial fqn */
        fqnname = NCD4_makeName(var,"_");
        if(fqnname == NULL)
	    {ret = NC_ENOMEM; goto done;}
	/* Now, create the vlen type using the field's basetype */
        if((ret=makeNode(parser,group,xml,NCD4_TYPE,NC_SEQ,&vlentype))) goto done;
        classify(group,vlentype);
	vlentype->basetype = var->basetype;
	/* Use name <fqnname>_t */
	strncpy(name,fqnname,sizeof(name));
	strncat(name,"_t", sizeof(name) - strlen(name) - 1);
        SETNAME(vlentype,name);
        /* Set the basetype */
        var->basetype = vlentype;
    } else {
	/* Start by creating the var node; will be fixed up later */
	if((ret=makeNode(parser,container,xml,NCD4_VAR,NC_SEQ,&var))) goto done;
	classify(container,var);
        fqnname = NCD4_makeName(var,"_");
        if(fqnname == NULL)
	    {ret = NC_ENOMEM; goto done;}
        if((ret=makeNode(parser,group,xml,NCD4_TYPE,NC_STRUCT,&structtype))) goto done;
        classify(group,structtype);
	/* Use name <fqnname>_base */
	strncpy(name,fqnname,sizeof(name));
	strncat(name,"_base", sizeof(name) - strlen(name) - 1);
        SETNAME(structtype,name);
        /* Parse Fields into type */
        if((ret = parseFields(parser,structtype,xml))) goto done;
	/* Create a seq type whose basetype is the compound type */
        if((ret=makeNode(parser,group,xml,NCD4_TYPE,NC_SEQ,&vlentype))) goto done;
        classify(group,vlentype);
	/* Use name <xname>_t */
	strncpy(name,fqnname,sizeof(name));
	strncat(name,"_t", sizeof(name) - strlen(name) - 1);
        SETNAME(vlentype,name);
	vlentype->basetype = structtype;
        /* Set the basetype */
        var->basetype = vlentype;
    }

    /* Parse attributes, dims, and maps into var*/
    if((ret = parseMetaData(parser,var,xml))) goto done;

    record(parser,var);

    /* See if this var has UCARTAGORIGTYPE attribute */
    if(parser->metadata->controller->controls.translation == NCD4_TRANSNC4) {
	const char* typetag = ezxml_attr(xml,UCARTAGORIGTYPE);
	if(typetag != NULL) {
	    /* yes, place it on the type */
	    if((ret=addOrigType(parser,var,vlentype,typetag))) goto done;
 	}
    }
    if(nodep) *nodep = var;

done:
    if(fqnname) free(fqnname);
    return THROW(ret);
}

static int
parseGroups(NCD4parser* parser, NCD4node* parent, ezxml_t xml)
{
    int ret = NC_NOERR;
    ezxml_t x;
    for(x=ezxml_child(xml, "Group");x != NULL;x = ezxml_next(x)) {
	NCD4node* group = NULL;
	const char* name = ezxml_attr(x,"name");
	if(name == NULL) FAIL(NC_EBADNAME,"Group has no name");
	if((ret=makeNode(parser,parent,x,NCD4_GROUP,NC_NULL,&group))) goto done;
	group->group.varbyid = nclistnew();
        if((ret = fillgroup(parser,group,x))) goto done;
        /* Parse group attributes */
        if((ret = parseAttributes(parser,group,x))) goto done;
	PUSH(parent->groups,group);
    }
done:
    return THROW(ret);
}

static int
parseAtomicVar(NCD4parser* parser, NCD4node* container, ezxml_t xml, NCD4node** nodep)
{
    int ret = NC_NOERR;
    NCD4node* node = NULL;
    NCD4node* base = NULL;
    const char* typename;
    const KEYWORDINFO* info;
    NCD4node* group;

    /* Check for aliases */
    for(typename=xml->name;;) {
	info = keyword(typename);
	if(info->aliasfor == NULL) break;
	typename = info->aliasfor;
    }
    group = NCD4_groupFor(container);
    /* Locate its basetype; handle opaque and enum separately */
    if(info->subsort == NC_ENUM) {
        const char* enumfqn = ezxml_attr(xml,"enum");
	if(enumfqn == NULL)
	    base = NULL;
	else
	    base = lookupFQN(parser,enumfqn,NCD4_TYPE);
    } else if(info->subsort == NC_OPAQUE) {
	/* See if the xml references an opaque type name */
	base = getOpaque(parser,xml,group);
    } else {
	base = lookupFQN(parser,info->tag,NCD4_TYPE);
    }
    if(base == NULL || !ISTYPE(base->sort)) {
	FAIL(NC_EBADTYPE,"Unexpected variable type: %s",info->tag);
    }
    if((ret=makeNode(parser,container,xml,NCD4_VAR,base->subsort,&node))) goto done;
    classify(container,node);
    node->basetype = base;
    /* Parse attributes, dims, and maps */
    if((ret = parseMetaData(parser,node,xml))) goto done;
    /* See if this var has UCARTAGORIGTYPE attribute */
    if(parser->metadata->controller->controls.translation == NCD4_TRANSNC4) {
	const char* typetag = ezxml_attr(xml,UCARTAGORIGTYPE);
	if(typetag != NULL) {
	    /* yes, place it on the type */
	    if((ret=addOrigType(parser,node,node,typetag))) goto done;
 	}
    }
    if(nodep) *nodep = node;
done:
    return THROW(ret);
}

static int
parseDimRefs(NCD4parser* parser, NCD4node* var, ezxml_t xml)
{
    int ret = NC_NOERR;
    ezxml_t x;
    for(x=ezxml_child(xml, "Dim");x!= NULL;x=ezxml_next(x)) {
	NCD4node* dim = NULL;
	const char* fqn;

	fqn = ezxml_attr(x,"name");
	if(fqn != NULL) {
   	    dim = lookupFQN(parser,fqn,NCD4_DIM);
	    if(dim == NULL) {
	        FAIL(NC_EBADDIM,"Cannot locate dim with name: %s",fqn);
	    }
	} else {
	    const char* sizestr = ezxml_attr(x,"size");
	    if(sizestr == NULL) {
	        FAIL(NC_EBADDIM,"Dimension reference has no name and no size");
	    }
	    /* Make or reuse anonymous dimension in root group */
	    dim = makeAnonDim(parser,sizestr);
	    if(dim == NULL)
		FAIL(NC_EBADDIM,"Cannot create anonymous dimension for size: %s",sizestr);
	}
	PUSH(var->dims,dim);
    }
done:
    return THROW(ret);
}

static int
parseMaps(NCD4parser* parser, NCD4node* var, ezxml_t xml)
{
    int ret = NC_NOERR;
    ezxml_t x;

    for(x=ezxml_child(xml, "Map");x!= NULL;x=ezxml_next(x)) {
	NCD4node* mapref = NULL;
	const char* fqn;
	fqn = ezxml_attr(x,"name");
	if(fqn == NULL)
	    FAIL(NC_ENOTVAR,"<Map> has no name attribute");
        mapref = lookupFQN(parser,fqn,NCD4_VAR);
	if(mapref == NULL)
	    FAIL(NC_ENOTVAR,"<Map> name does not refer to a variable: %s",fqn);
	PUSH(var->maps,mapref);
    }
done:
    return THROW(ret);
}

static int
parseAttributes(NCD4parser* parser, NCD4node* container, ezxml_t xml)
{
    int ret = NC_NOERR;
    ezxml_t x;
    NClist* values = NULL;

    /* First, transfer any reserved xml attributes */
    {
	int count = 0;
	const char** all = NULL;
	all = ezxml_all_attr(xml,&count);
	if(all != NULL && count > 0) {
	    const char** p;
	    container->xmlattributes = nclistnew();
	    for(p=all;*p;p+=2) {
		if(isReserved(*p)) {
		    nclistpush(container->xmlattributes,strdup(p[0]));
		    nclistpush(container->xmlattributes,strdup(p[1]));
		}
	    }
	}
    }

    for(x=ezxml_child(xml, "Attribute");x!= NULL;x=ezxml_next(x)) {
	const char* name = ezxml_attr(x,"name");
	const char* type = ezxml_attr(x,"type");
	NCD4node* attr = NULL;
	NCD4node* basetype;

	if(name == NULL) FAIL(NC_EBADNAME,"Missing <Attribute> name");
#ifdef HYRAXHACK
	/* Hyrax specifies type="container" for container types */
	if(strcmp(type,"container")==0
	   || strcmp(type,"Container")==0)
	    type = NULL;
#endif
	if(type == NULL) {
	    /* <Attribute> containers not supported; ignore */
	    continue;
	}

	if((ret=makeNode(parser,container,x,NCD4_ATTR,NC_NULL,&attr))) goto done;
	basetype = lookupFQN(parser,type,NCD4_TYPE);
	if(basetype == NULL)
	    FAIL(NC_EBADTYPE,"Unknown <Attribute> type: ",type);
	if(basetype->subsort == NC_NAT && basetype->subsort != NC_ENUM)
	    FAIL(NC_EBADTYPE,"<Attribute> type must be atomic or enum: ",type);
	attr->basetype = basetype;
	values = nclistnew();
	if((ret=getValueStrings(parser,basetype,x,values))) {
	    FAIL(NC_EINVAL,"Malformed attribute: %s",name);
	}
	attr->attr.values = values; values = NULL;
	PUSH(container->attributes,attr);
    }
done:
    if(ret != NC_NOERR) {
        nclistfreeall(values);
    }
    return THROW(ret);
}

static int
parseError(NCD4parser* parser, ezxml_t errxml)
{
    const char* shttpcode = ezxml_attr(errxml,"httpcode");
    ezxml_t x;
    if(shttpcode == NULL) shttpcode = "400";
    if(sscanf(shttpcode,"%d",&parser->metadata->error.httpcode) != 1)
        nclog(NCLOGERR,"Malformed <ERROR> response");
    x=ezxml_child(errxml, "Message");
    if(x != NULL) {
	const char* txt = ezxml_txt(x);
	parser->metadata->error.message = (txt == NULL ? NULL : strdup(txt));
    }
    x=ezxml_child(errxml, "Context");
    if(x != NULL) {
	const char* txt = ezxml_txt(x);
	parser->metadata->error.context = (txt == NULL ? NULL : strdup(txt));
    }
    x=ezxml_child(errxml, "OtherInformation");
    if(x != NULL) {
	const char* txt = ezxml_txt(x);
	parser->metadata->error.otherinfo = (txt == NULL ? NULL : strdup(txt));
    }
    return THROW(NC_NOERR);
}

/*
Find or create an opaque type
*/
static NCD4node*
getOpaque(NCD4parser* parser, ezxml_t varxml, NCD4node* group)
{
    int i, ret = NC_NOERR;
    long long len;
    NCD4node* opaquetype = NULL;
    const char* xattr;

#ifndef FIXEDOPAQUE
    len = 0;
#else
    len = parser->metadata->controller->controls.opaquesize;
#endif
    if(parser->metadata->controller->controls.translation == NCD4_TRANSNC4) {
        /* See if this var has UCARTAGOPAQUE attribute */
        xattr = ezxml_attr(varxml,UCARTAGOPAQUE);
        if(xattr != NULL) {
	    long long tmp = 0;
            if((ret = parseLL(xattr,&tmp)) || (tmp < 0))
	        FAIL(NC_EINVAL,"Illegal opaque len: %s",xattr);
	    len = tmp;
        }
    }
#ifndef FIXEDOPAQUE
    if(len == 0) {
        /* Need to use _bytestring */
	if((ret=defineBytestringType(parser)))
  	    goto done;
	assert(parser->metadata->_bytestring != NULL);
	opaquetype = parser->metadata->_bytestring;
    } else
#endif
    { /*(len > 0) || FIXEDOPAQUE */
        /* Try to locate existing opaque type with this length */
        for(i=0;i<nclistlength(parser->types); i++) {
	    NCD4node* op = (NCD4node*)nclistget(parser->types,i);
	    if(op->subsort != NC_OPAQUE) continue;
	    if(op->opaque.size == len) {opaquetype = op; break;}
	}
        if(opaquetype == NULL) {/* create it */
	    char name[NC_MAX_NAME+1];
	    /* Make name be "opaqueN" */
	    snprintf(name,NC_MAX_NAME,"opaque%lld_t",len);
	    /* Opaque types are always created in the current group */
	    if((ret=makeNode(parser,group,NULL,NCD4_TYPE,NC_OPAQUE,&opaquetype)))
	        goto done;
  	    SETNAME(opaquetype,name);
	    opaquetype->opaque.size = len;
	    if(opaquetype != NULL)
	        record(parser,opaquetype);
	}
    }
done:
    return opaquetype;
}

/* get all value strings */
static int
getValueStrings(NCD4parser* parser, NCD4node* type, ezxml_t xattr, NClist* svalues)
{
    const char* s;
    /* See first if we have a "value" xml attribute */
    s = ezxml_attr(xattr,"value");
    if(s != NULL)
	PUSH(svalues,strdup(s));
    else {/* look for <Value> subnodes */
	ezxml_t x;
        for(x=ezxml_child(xattr, "Value");x != NULL;x = ezxml_next(x)) {
	    char* es;
	    char* ds;
	    /* We assume that either their is a single xml attribute called "value",
               or there is a single chunk of text containing possibly multiple values.
	    */
	    s = ezxml_attr(x,"value");
	    if(s == NULL) {/* See if there is a text part. */
		s = x->txt;
		if(s == NULL) s = "";
	    }
	    /* Need to de-escape the string */
	    es = NCD4_entityescape(s);
	    ds = NCD4_deescape(es);
	    nclistpush(svalues,ds);
	    nullfree(es);
	}
    }
    return THROW(NC_NOERR);
}

/***************************************************/
/* Utilities */

NCD4node*
NCD4_groupFor(NCD4node* node)
{
    while(node->sort != NCD4_GROUP) node = node->container;
    return node;
}

/* Determine is a name is reserved */
static int
isReserved(const char* name)
{
    if(name == NULL) return 0;
    return (name[0] == RESERVECHAR);
}

/* If a node has the UCARTAGORIGTYPE attribute,
   then capture that annotation. */
static int
addOrigType(NCD4parser* parser, NCD4node* src, NCD4node* dst, const char* oldname)
{
    int ret = NC_NOERR;

    if(dst == NULL) dst = src;
    /* Record the original type in the destination*/
    if((ret=splitOrigType(parser,oldname,dst))) goto done;
done:
    return THROW(ret);
}

static int
splitOrigType(NCD4parser* parser, const char* fqn, NCD4node* type)
{
    int ret = NC_NOERR;
    NClist* pieces = nclistnew();
    NCD4node* group = NULL;
    char* name = NULL;

    if((ret=NCD4_parseFQN(fqn,pieces))) goto done;
    /* It should be the case that the pieces are {/group}+/name */
    name = (char*)nclistpop(pieces);
    if((ret = lookupFQNList(parser,pieces,NCD4_GROUP,&group))) goto done;
    if(group == NULL) {
	FAIL(NC_ENOGRP,"Non-existent group in FQN: ",fqn);
    }
    type->nc4.orig.name = strdup(name+1); /* plus 1 to skip the leading separator */
    type->nc4.orig.group = group;

done:
    return THROW(ret);
}

/* Locate an attribute.
   If not found, then *attrp will be null
*/
NCD4node*
NCD4_findAttr(NCD4node* container, const char* attrname)
{
    int i;
    /* Look directly under this xml for <Attribute> */
    for(i=0;i<nclistlength(container->attributes);i++) {
	NCD4node* attr = (NCD4node*)nclistget(container->attributes,i);
	if(strcmp(attr->name,attrname)!=0) continue;
	return attr;
    }
    return NULL;
}

/*
Parse a simple string of digits into an unsigned long long
Return the value.
*/

static int
parseULL(const char* text, unsigned long long* ullp)
{
    extern int errno;
    char* endptr;
    unsigned long long uint64 = 0;

    errno = 0; endptr = NULL;
#ifdef HAVE_STRTOULL
    uint64 = strtoull(text,&endptr,10);
    if(errno == ERANGE)
	return THROW(NC_ERANGE);
#else /*!(defined HAVE_STRTOLL && defined HAVE_STRTOULL)*/
    sscanf((char*)text, "%llu", &uint64);
    /* Have no useful way to detect out of range */
#endif /*!(defined HAVE_STRTOLL && defined HAVE_STRTOULL)*/
    if(ullp) *ullp = uint64;
    return THROW(NC_NOERR);
}

/*
Parse a simple string of digits into an signed long long
Return the value.
*/

static int
parseLL(const char* text, long long* llp)
{
    extern int errno;
    char* endptr;
    long long int64 = 0;

    errno = 0; endptr = NULL;
#ifdef HAVE_STRTOLL
    int64 = strtoll(text,&endptr,10);
    if(errno == ERANGE)
	return THROW(NC_ERANGE);
#else /*!(defined HAVE_STRTOLL && defined HAVE_STRTOLL)*/
    sscanf((char*)text, "%lld", &int64);
    /* Have no useful way to detect out of range */
#endif /*!(defined HAVE_STRTOLL && defined HAVE_STRTOLL)*/
    if(llp) *llp = int64;
    return THROW(NC_NOERR);
}

/*
Convert a sequence of fqn names into a specific node.
WARNING: This is highly specialized in that it assumes
that the final object is one of: dimension, type, or var.
This means that e.g. groups, attributes, econsts, cannot
be found by this procedure.
*/
static int
lookupFQNList(NCD4parser* parser, NClist* fqn, NCD4sort sort, NCD4node** result)
{
    int ret = NC_NOERR;
    int i,nsteps;
    NCD4node* current;
    char* name = NULL;
    NCD4node* node = NULL;

    /* Step 1: walk thru groups until can go no further */
    current = parser->metadata->root;
    nsteps = nclistlength(fqn);
    for(i=1;i<nsteps;i++) { /* start at 1 to side-step root name */
	assert(ISGROUP(current->sort));
	name = (char*)nclistget(fqn,i);
        /* See if we can find a matching subgroup */
	node = lookFor(current->group.elements,name,NCD4_GROUP);
	if(node == NULL)
	    break; /* reached the end of the group part of the fqn */
	current = node;
    }
    /* Invariant:
	1. i == nsteps => node != null => last node was a group:
                                          it must be our target
	2. i == (nsteps-1) => non-group node at the end; disambiguate
	3. i < (nsteps - 1) => need a compound var to continue
    */
    if(i == nsteps) {
	if(sort != NCD4_GROUP) goto sortfail;
	goto done;
    }
    if(i == (nsteps - 1)) {
	assert (node == NULL);
        node = lookFor(current->group.elements,name,sort);
	if(node == NULL) goto sortfail;
	goto done;
    }
    assert (i < (nsteps - 1)); /* case 3 */
    /* We have steps to take, so node better be a compound var */
    node = lookFor(current->group.elements,name,NCD4_VAR);
    if(node == NULL || !ISCMPD(node->basetype->subsort))
	goto fail;
    /* So we are at a compound variable, so walk its fields recursively */
    /* Use the type to do the walk */
    current = node->basetype;
    assert (i < (nsteps - 1));
    i++; /* skip variable name */
    for(;;i++) {
	int j;
	name = (char*)nclistget(fqn,i);
	assert(ISTYPE(current->sort) && ISCMPD(current->subsort));
	for(node=NULL,j=0;j<nclistlength(current->vars);j++) {
	    NCD4node* field = (NCD4node*)nclistget(current->vars,j);
	    if(strcmp(field->name,name)==0)
		{node = field; break;}
	}
	if(node == NULL)
	    goto sortfail; /* no match, so failed */
	if(i == (nsteps - 1))
	    break;
	if(!ISCMPD(node->basetype->subsort))
	    goto fail; /* more steps, but no compound field, so failed */
	current = node->basetype;
    }
done:
    if(result) *result = node;
    return THROW(ret);
fail:
    ret = NC_EINVAL;
    goto done;
sortfail:
    ret = NC_EBADID;
    goto done;
}

static NCD4node*
lookFor(NClist* elems, const char* name, NCD4sort sort)
{
    int n,i;
    if(elems == NULL || nclistlength(elems) == 0) return NULL;
    n = nclistlength(elems);
    for(i=0;i<n;i++) {
	NCD4node* node = (NCD4node*)nclistget(elems,i);
	if(strcmp(node->name,name) == 0 && (sort == node->sort))
	    return node;
    }
    return NULL;
}

void
NCD4_printElems(NCD4node* group)
{
    int n,i;
    NClist* elems;
    elems = group->group.elements;
    if(elems == NULL || nclistlength(elems) == 0) return;
    n = nclistlength(elems);
    for(i=0;i<n;i++) {
	NCD4node* node = (NCD4node*)nclistget(elems,i);
	fprintf(stderr,"name=%s sort=%d subsort=%d\n",
		node->name,node->sort,node->subsort);
    }
    fflush(stderr);
}

static NCD4node*
lookupFQN(NCD4parser* parser, const char* sfqn, NCD4sort sort)
{
    int ret = NC_NOERR;
    NClist* fqnlist = nclistnew();
    NCD4node* match = NULL;

    /* Short circuit atomic types */
    if(NCD4_TYPE == sort) {
        match = lookupAtomictype(parser,(sfqn[0]=='/'?sfqn+1:sfqn));
        if(match != NULL)
	    goto done;
    }
    if((ret=NCD4_parseFQN(sfqn,fqnlist))) goto done;
    if((ret=lookupFQNList(parser,fqnlist,sort,&match))) goto done;
done:
    nclistfreeall(fqnlist);
    return (ret == NC_NOERR ? match : NULL);
}

static const KEYWORDINFO*
keyword(const char* name)
{
    int n = sizeof(keywordmap)/sizeof(KEYWORDINFO);
    int L = 0;
    int R = (n - 1);
    int m, cmp;
    const struct KEYWORDINFO* p;
    for(;;) {
	if(L > R) break;
        m = (L + R) / 2;
	p = &keywordmap[m];
	cmp = strcasecmp(p->tag,name);
	if(cmp == 0) return p;
	if(cmp < 0)
	    L = (m + 1);
	else /*cmp > 0*/
	    R = (m - 1);
    }
    return NULL;
}

#ifndef FIXEDOPAQUE
static int
defineBytestringType(NCD4parser* parser)
{
    int ret = NC_NOERR;
    NCD4node* bstring = NULL;
    if(parser->metadata->_bytestring == NULL) {
        /* Construct a single global opaque type for mapping DAP opaque type */
        ret = makeNode(parser,parser->metadata->root,NULL,NCD4_TYPE,NC_OPAQUE,&bstring);
        if(ret != NC_NOERR) goto done;
        SETNAME(bstring,"_bytestring");
	bstring->opaque.size = 0;
	bstring->basetype = lookupAtomictype(parser,"UInt8");
        PUSH(parser->metadata->root->types,bstring);
	parser->metadata->_bytestring = bstring;
    } else
	bstring = parser->metadata->_bytestring;
done:
    return THROW(ret);
}
#endif

static int
defineAtomicTypes(NCD4parser* parser)
{
    int ret = NC_NOERR;
    NCD4node* node;
    const struct ATOMICTYPEINFO* ati;

    parser->atomictypes = nclistnew();
    if(parser->atomictypes == NULL)
	return THROW(NC_ENOMEM);
    for(ati=atomictypeinfo;ati->name;ati++) {
        if((ret=makeNode(parser,parser->metadata->root,NULL,NCD4_TYPE,ati->type,&node))) goto done;
	SETNAME(node,ati->name);
        node->container = parser->metadata->root;
	record(parser,node);
	PUSH(parser->atomictypes,node);
    }
    parser->used = (char*)calloc(1,nclistlength(parser->atomictypes));
    if(parser->used == NULL) {ret = NC_ENOMEM; goto done;}

done:
    return THROW(ret);
}

/* Binary search the set of set of atomictypes */
static NCD4node*
lookupAtomictype(NCD4parser* parser, const char* name)
{
    int n = nclistlength(parser->atomictypes);
    int L = 0;
    int R = (n - 1);
    int m, cmp;
    NCD4node* p;

    for(;;) {
	if(L > R) break;
        m = (L + R) / 2;
	p = (NCD4node*)nclistget(parser->atomictypes,m);
	cmp = strcasecmp(p->name,name);
	if(cmp == 0) return p;
	if(cmp < 0)
	    L = (m + 1);
	else /*cmp > 0*/
	    R = (m - 1);
    }
    return NULL;
}

/**************************************************/

static int
makeNode(NCD4parser* parser, NCD4node* parent, ezxml_t xml, NCD4sort sort, nc_type subsort, NCD4node** nodep)
{
    int ret = NC_NOERR;
    NCD4node* node = (NCD4node*)calloc(1,sizeof(NCD4node));

    if(node == NULL) return THROW(NC_ENOMEM);
    node->sort = sort;
    node->subsort = subsort;
    node->container = parent;
    /* Set node name, if it exists */
    if(xml != NULL) {
        const char* name = ezxml_attr(xml,"name");
        if(name != NULL) {
	    if(strlen(name) > NC_MAX_NAME) {
	        nclog(NCLOGERR,"Name too long: %s",name);
	    }
	    SETNAME(node,name);
	}
    }
    if(parent != NULL) {
	if(parent->sort == NCD4_GROUP)
	    PUSH(parent->group.elements,node);
    }
    track(parser,node);
    if(nodep) *nodep = node;
    return THROW(ret);
}

static NCD4node*
makeAnonDim(NCD4parser* parser, const char* sizestr)
{
    long long size = 0;
    int ret;
    char name[NC_MAX_NAME+1];
    NCD4node* dim = NULL;
    NCD4node* root = parser->metadata->root;

    ret = parseLL(sizestr,&size);
    if(ret) return NULL;
    snprintf(name,NC_MAX_NAME,"/_Anonymous%lld",size);
    /* See if it exists already */
    dim = lookupFQN(parser,name,NCD4_DIM);
    if(dim == NULL) {/* create it */
	if((ret=makeNode(parser,root,NULL,NCD4_DIM,NC_NULL,&dim))) goto done;
	SETNAME(dim,name+1); /* leave out the '/' separator */
	dim->dim.size = (long long)size;
	dim->dim.isanonymous = 1;
	PUSH(root->dims,dim);
    }
done:
    return (ret?NULL:dim);
}

/*
Classify inserts the node into the proper container list
based on the node's sort.
*/
static void
classify(NCD4node* container, NCD4node* node)
{
    if(ISGROUP(container->sort))
	nclistpush(container->group.elements,node);
    switch (node->sort) {
    case NCD4_GROUP:
	PUSH(container->groups,node);
        break;
    case NCD4_DIM:
	PUSH(container->dims,node);
        break;
    case NCD4_TYPE:
	PUSH(container->types,node);
        break;
    case NCD4_VAR:
	PUSH(container->vars,node);
        break;
    case NCD4_ATTR: case NCD4_XML:
	PUSH(container->attributes,node);
        break;
    default: break;
    }
}

/*
Classify inserts the node into the proper parser global list
based on the node's sort.
*/
static void
record(NCD4parser* parser, NCD4node* node)
{
    switch (node->sort) {
    case NCD4_GROUP:
	PUSH(parser->groups,node);
        break;
    case NCD4_DIM:
	PUSH(parser->dims,node);
        break;
    case NCD4_TYPE:
	PUSH(parser->types,node);
        break;
    case NCD4_VAR:
	PUSH(parser->vars,node);
        break;
    default: break;
    }
}

/*
Undo a classify and record
for a field node.
Used by buildSequenceType.
*/
#if 0
static void
forget(NCD4parser* parser, NCD4node* var)
{
    int i;
    NCD4node* container = var->container;
    assert(ISVAR(var->sort) && ISTYPE(container->sort) && ISCMPD(container->subsort));
    /* Unrecord: remove from the parser lists */
    for(i=0;i<parser->vars;i++) {
	NCD4node* test = nclistget(parser->vars,i);
	if(test == var) {
	    nclistremove(parser->vars,i);
	    break;
	}
    }
    /* Unclassify: remove from the container var list */
    for(i=0;i<container->vars;i++) {
	NCD4node* test = nclistget(container->vars,i);
	if(test == var) {
	    nclistremove(container->vars,i);
	    break;
	}
    }
}
#endif

static void
track(NCD4parser* parser, NCD4node* node)
{
#ifdef D4DEBUG
    fprintf(stderr,"track: node=%lx sort=%d subsort=%d",(unsigned long)node,node->sort,node->subsort);
    if(node->name != NULL)
        fprintf(stderr," name=%s\n",node->name);
    fprintf(stderr,"\n");
#endif
    PUSH(parser->metadata->allnodes,node);
#ifdef D4DEBUG
    fprintf(stderr,"track: |allnodes|=%ld\n",nclistlength(parser->metadata->allnodes));
    fflush(stderr);
#endif
}

/**************************************************/

static int
convertString(union ATOMICS* converter, NCD4node* type, const char* s)
{
    switch (type->subsort) {
    case NC_BYTE:
    case NC_SHORT:
    case NC_INT:
    case NC_INT64:
	if(sscanf(s,"%lld",converter->i64) != 1) return THROW(NC_ERANGE);
	break;
    case NC_UBYTE:
    case NC_USHORT:
    case NC_UINT:
    case NC_UINT64:
	if(sscanf(s,"%llu",converter->u64) != 1) return THROW(NC_ERANGE);
	break;
    case NC_FLOAT:
    case NC_DOUBLE:
	if(sscanf(s,"%lf",converter->f64) != 1) return THROW(NC_ERANGE);
	break;
    case NC_STRING:
	converter->s[0]= strdup(s);
	break;
    }/*switch*/
    return downConvert(converter,type);
}

static int
downConvert(union ATOMICS* converter, NCD4node* type)
{
    unsigned long long u64 = converter->u64[0];
    long long i64 = converter->i64[0];
    double f64 = converter->f64[0];
    char* s = converter->s[0];
    switch (type->subsort) {
    case NC_BYTE:
	converter->i8[0] = (char)i64;
	break;
    case NC_UBYTE:
	converter->u8[0] = (unsigned char)u64;
	break;
    case NC_SHORT:
	converter->i16[0] = (short)i64;
	break;
    case NC_USHORT:
	converter->u16[0] = (unsigned short)u64;
	break;
    case NC_INT:
	converter->i32[0] = (int)i64;
	break;
    case NC_UINT:
	converter->u32[0] = (unsigned int)u64;
	break;
    case NC_INT64:
	converter->i64[0] = i64;
	break;
    case NC_UINT64:
	converter->u64[0]= u64;
	break;
    case NC_FLOAT:
	converter->f32[0] = (float)f64;
	break;
    case NC_DOUBLE:
	converter->f64[0] = f64;
	break;
    case NC_STRING:
	converter->s[0]= s;
	break;
    }/*switch*/
    return THROW(NC_NOERR);
}

#if 0
/* Try to remove excess text from a value set */
static int
valueParse(NCD4node* type, const char* values0, NClist* vlist)
{
    char* values;
    char* p;
    char* q;
    char* s;
    char* line;
    ptrdiff_t len;

    if(values0 == NULL || (len=strlen(values0)) == 0)
	return THROW(NC_NOERR);
    values = strdup(values0);
    /* Compress the text by removing sequences of blanks and newlines:
       note that this will fail for string typed values that might have
       embedded blanks, so in that case, we assume each string is on a separate line.
       For NC_CHAR, we treat like strings, except we use the last char in the line
       as the value. This is all heuristic.
    */
    switch (type->subsort) {
    case NC_STRING:
	p = values;
	for(;;) {
	    if(*p == '\0') break;
	    line = p;
	    /* Start by looking for \n or \r\n */
            for(;*p;p++) {if(*p == '\n') break;}
	    q = p - 1;
	    *p++ = '\0';
	    if(*q == '\r') {*q = '\0';}
            nclistpush(vlist,strdup(line));
	}
	break;
    case NC_CHAR:
	p = values;
	for(;*p;) {
	    char c[2];
	    line = p;
	    /* Start by looking for \n or \r\n */
            for(;*p;p++) {if(*p == '\n') break;}
	    q = p;
	    *p++ = '\0';
	    q--;
	    if(*q == '\r') {*q = '\0';}
	    len = strlen(line);
	    if(len > 0) {
		c[0] = *q;
		c[1] = '\0';
		nclistpush(vlist,strdup(c));
	    }
	}
	break;
    default:
	p = values;
        for(;*p;p++) {if(*p <= ' ') *p = '\n';}
	line = values;
	for(p=line;*p;p++) {if(*p == '\n') break;}
	for(line=values;*line;) {
	    size_t size = strlen(line);
	    if(size > 0)
	        nclistpush(vlist,strdup(line));
	    line += (size+1); /* skip terminating nul */
	}
	break;
    }
    free(values);
    return THROW(NC_NOERR);
}
#endif
