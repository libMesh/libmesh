/*********************************************************************
 *   Copyright 2009, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/
/* $Id: semantics.c,v 1.4 2010/05/24 19:59:58 dmh Exp $ */
/* $Header: /upc/share/CVS/netcdf-3/ncgen/semantics.c,v 1.4 2010/05/24 19:59:58 dmh Exp $ */

#include        "includes.h"
#include        "dump.h"
#include        "offsets.h"

/* Forward*/
static void computefqns(void);
static void filltypecodes(void);
static void processenums(void);
static void processeconstrefs(void);
static void processtypes(void);
static void processtypesizes(void);
static void processvars(void);
static void processattributes(void);
static void processunlimiteddims(void);
static void processeconstrefs(void);
static void processeconstrefsR(Datalist*);

static List* ecsearchgrp(Symbol* grp, List* candidates);
static List* findecmatches(char* ident);
static void fixeconstref(NCConstant* con);
static void inferattributetype(Symbol* asym);
static void validateNIL(Symbol* sym);
static void checkconsistency(void);
static int tagvlentypes(Symbol* tsym);

static void computefqns(void);

static Symbol* uniquetreelocate(Symbol* refsym, Symbol* root);
static Symbol* checkeconst(Symbol* en, const char* refname);


List* vlenconstants;  /* List<Constant*>;*/
			  /* ptr to vlen instances across all datalists*/

/* Post-parse semantic checks and actions*/
void
processsemantics(void)
{
    /* Fill in the fqn for every defining symbol */
    computefqns();
    /* Process each type and sort by dependency order*/
    processtypes();
    /* Make sure all typecodes are set if basetype is set*/
    filltypecodes();
    /* Process each type to compute its size*/
    processtypesizes();
    /* Process each var to fill in missing fields, etc*/
    processvars();
    /* Process attributes to connect to corresponding variable*/
    processattributes();
    /* Fix up enum constant values*/
    processenums();
    /* Fix up enum constant references*/
    processeconstrefs();
    /* Compute the unlimited dimension sizes */
    processunlimiteddims();
    /* check internal consistency*/
    checkconsistency();
}

/*
Given a reference symbol, produce the corresponding
definition symbol; return NULL if there is no definition
Note that this is somewhat complicated to conform to
various scoping rules, namely:
1. look into parent hierarchy for un-prefixed dimension names.
2. look in whole group tree for un-prefixed type names;
   search is depth first. MODIFIED 5/26/2009: Search is as follows:
   a. search parent hierarchy for matching type names.
   b. search whole tree for unique matching type name
   c. complain and require prefixed name.
3. look in the same group as ref for un-prefixed variable names.
4. ditto for group references
5. look in whole group tree for un-prefixed enum constants;
   result must be unique
*/

Symbol*
locate(Symbol* refsym)
{
    Symbol* sym = NULL;
    switch (refsym->objectclass) {
    case NC_DIM:
	if(refsym->is_prefixed) {
	    /* locate exact dimension specified*/
	    sym = lookup(NC_DIM,refsym);
	} else { /* Search for matching dimension in all parent groups*/
	    Symbol* parent = lookupgroup(refsym->prefix);/*get group for refsym*/
	    while(parent != NULL) {
		/* search this parent for matching name and type*/
		sym = lookupingroup(NC_DIM,refsym->name,parent);
		if(sym != NULL) break;
		parent = parent->container;
	    }
	}		
	break;
    case NC_TYPE:
	if(refsym->is_prefixed) {
	    /* locate exact type specified*/
	    sym = lookup(NC_TYPE,refsym);
	} else {
	    Symbol* parent;
	    int i; /* Search for matching type in all groups (except...)*/
	    /* Short circuit test for primitive types*/
	    for(i=NC_NAT;i<=NC_STRING;i++) {
		Symbol* prim = basetypefor(i);
		if(prim == NULL) continue;
	        if(strcmp(refsym->name,prim->name)==0) {
		    sym = prim;
		    break;
		}
	    }
	    if(sym == NULL) {
	        /* Added 5/26/09: look in parent hierarchy first */
	        parent = lookupgroup(refsym->prefix);/*get group for refsym*/
	        while(parent != NULL) {
		    /* search this parent for matching name and type*/
		    sym = lookupingroup(NC_TYPE,refsym->name,parent);
		    if(sym != NULL) break;
		    parent = parent->container;
		}
	    }
	    if(sym == NULL) {
	        sym = uniquetreelocate(refsym,rootgroup); /* want unique */
	    }
	}		
	break;
    case NC_VAR:
	if(refsym->is_prefixed) {
	    /* locate exact variable specified*/
	    sym = lookup(NC_VAR,refsym);
	} else {
	    Symbol* parent = lookupgroup(refsym->prefix);/*get group for refsym*/
   	    /* search this parent for matching name and type*/
	    sym = lookupingroup(NC_VAR,refsym->name,parent);
	}		
        break;
    case NC_GRP:
	if(refsym->is_prefixed) {
	    /* locate exact group specified*/
	    sym = lookup(NC_GRP,refsym);
	} else {
 	    Symbol* parent = lookupgroup(refsym->prefix);/*get group for refsym*/
   	    /* search this parent for matching name and type*/
	    sym = lookupingroup(NC_GRP,refsym->name,parent);
	}		
	break;

    default: PANIC1("locate: bad refsym type: %d",refsym->objectclass);
    }
    if(debug > 1) {
	char* ncname;
	if(refsym->objectclass == NC_TYPE)
	    ncname = ncclassname(refsym->subclass);
	else
	    ncname = ncclassname(refsym->objectclass);
	fdebug("locate: %s: %s -> %s\n",
		ncname,fullname(refsym),(sym?fullname(sym):"NULL"));
    }   
    return sym;
}

/*
Search for an object in all groups using preorder depth-first traversal.
Return NULL if symbol is not unique or not found at all.
*/
static Symbol*
uniquetreelocate(Symbol* refsym, Symbol* root)
{
    int i;
    Symbol* sym = NULL;
    /* search the root for matching name and major type*/
    sym = lookupingroup(refsym->objectclass,refsym->name,root);
    if(sym == NULL) {
	for(i=0;i<listlength(root->subnodes);i++) {
	    Symbol* grp = (Symbol*)listget(root->subnodes,i);
	    if(grp->objectclass == NC_GRP && !grp->ref.is_ref) {
		Symbol* nextsym = uniquetreelocate(refsym,grp);
		if(nextsym != NULL) {
		    if(sym != NULL) return NULL; /* not unique */	
		    sym = nextsym;
		}
	    }
	}
    }
    return sym;
}


/*
Compute the fqn for every top-level definition symbol
*/
static void
computefqns(void)
{
    int i,j;
    /* Groups first */
    for(i=0;i<listlength(grpdefs);i++) {
        Symbol* sym = (Symbol*)listget(grpdefs,i);
	topfqn(sym);
    }
    /* Dimensions */
    for(i=0;i<listlength(dimdefs);i++) {
        Symbol* sym = (Symbol*)listget(dimdefs,i);
	topfqn(sym);
    }
    /* types */
    for(i=0;i<listlength(typdefs);i++) {
        Symbol* sym = (Symbol*)listget(typdefs,i);
	topfqn(sym);
    }
    /* variables */
    for(i=0;i<listlength(vardefs);i++) {
        Symbol* sym = (Symbol*)listget(vardefs,i);
	topfqn(sym);
    }
    /* fill in the fqn names of econsts */
    for(i=0;i<listlength(typdefs);i++) {
        Symbol* sym = (Symbol*)listget(typdefs,i);
	if(sym->subclass == NC_ENUM) {
	    for(j=0;j<listlength(sym->subnodes);j++) {
		Symbol* econ = (Symbol*)listget(sym->subnodes,j);
		nestedfqn(econ);
	    }
	}
    }
    /* fill in the fqn names of fields */
    for(i=0;i<listlength(typdefs);i++) {
        Symbol* sym = (Symbol*)listget(typdefs,i);
	if(sym->subclass == NC_COMPOUND) {
	    for(j=0;j<listlength(sym->subnodes);j++) {
		Symbol* field = (Symbol*)listget(sym->subnodes,j);
		nestedfqn(field);
	    }
	}
    }
    /* fill in the fqn names of attributes */
    for(i=0;i<listlength(gattdefs);i++) {
        Symbol* sym = (Symbol*)listget(gattdefs,i);
        attfqn(sym);
    }
    for(i=0;i<listlength(attdefs);i++) {
        Symbol* sym = (Symbol*)listget(attdefs,i);
        attfqn(sym);
    }
}

/* 1. Do a topological sort of the types based on dependency*/
/*    so that the least dependent are first in the typdefs list*/
/* 2. fill in type typecodes*/
/* 3. mark types that use vlen*/
static void
processtypes(void)
{
    int i,j,keep,added;
    List* sorted = listnew(); /* hold re-ordered type set*/
    /* Prime the walk by capturing the set*/
    /*     of types that are dependent on primitive types*/
    /*     e.g. uint vlen(*) or primitive types*/
    for(i=0;i<listlength(typdefs);i++) {
        Symbol* sym = (Symbol*)listget(typdefs,i);
	keep=0;
	switch (sym->subclass) {
	case NC_PRIM: /*ignore pre-defined primitive types*/
	    sym->touched=1;
	    break;
	case NC_OPAQUE:
	case NC_ENUM:
	    keep=1;
	    break;
        case NC_VLEN: /* keep if its basetype is primitive*/
	    if(sym->typ.basetype->subclass == NC_PRIM) keep=1;
	    break;	    	
	case NC_COMPOUND: /* keep if all fields are primitive*/
	    keep=1; /*assume all fields are primitive*/
	    for(j=0;j<listlength(sym->subnodes);j++) {
		Symbol* field = (Symbol*)listget(sym->subnodes,j);
		ASSERT(field->subclass == NC_FIELD);
		if(field->typ.basetype->subclass != NC_PRIM) {keep=0;break;}
	    }	  
	    break;
	default: break;/* ignore*/
	}
	if(keep) {
	    sym->touched = 1;
	    listpush(sorted,(void*)sym);
	}
    }	
    /* 2. repeated walk to collect level i types*/
    do {
        added=0;
        for(i=0;i<listlength(typdefs);i++) {
	    Symbol* sym = (Symbol*)listget(typdefs,i);
	    if(sym->touched) continue; /* ignore already processed types*/
	    keep=0; /* assume not addable yet.*/
	    switch (sym->subclass) {
	    case NC_PRIM: 
	    case NC_OPAQUE:
	    case NC_ENUM:
		PANIC("type re-touched"); /* should never happen*/
	        break;
            case NC_VLEN: /* keep if its basetype is already processed*/
	        if(sym->typ.basetype->touched) keep=1;
	        break;	    	
	    case NC_COMPOUND: /* keep if all fields are processed*/
	        keep=1; /*assume all fields are touched*/
	        for(j=0;j<listlength(sym->subnodes);j++) {
		    Symbol* field = (Symbol*)listget(sym->subnodes,j);
		    ASSERT(field->subclass == NC_FIELD);
		    if(!field->typ.basetype->touched) {keep=1;break;}
	        }	  
	        break;
	    default: break;				
	    }
	    if(keep) {
		listpush(sorted,(void*)sym);
		sym->touched = 1;
		added++;
	    }	    
	}
    } while(added > 0);
    /* Any untouched type => circular dependency*/
    for(i=0;i<listlength(typdefs);i++) {
	Symbol* tsym = (Symbol*)listget(typdefs,i);
	if(tsym->touched) continue;
	semerror(tsym->lineno,"Circular type dependency for type: %s",fullname(tsym));
    }
    listfree(typdefs);
    typdefs = sorted;
    /* fill in type typecodes*/
    for(i=0;i<listlength(typdefs);i++) {
        Symbol* sym = (Symbol*)listget(typdefs,i);
	if(sym->typ.basetype != NULL && sym->typ.typecode == NC_NAT)
	    sym->typ.typecode = sym->typ.basetype->typ.typecode;
    }
    /* Identify types containing vlens */
    for(i=0;i<listlength(typdefs);i++) {
        Symbol* tsym = (Symbol*)listget(typdefs,i);
	tagvlentypes(tsym);
    }
}

/* Recursively check for vlens*/
static int
tagvlentypes(Symbol* tsym)
{
    int tagged = 0;
    int j;
    switch (tsym->subclass) {
        case NC_VLEN: 
	    tagged = 1;
	    tagvlentypes(tsym->typ.basetype);
	    break;	    	
	case NC_COMPOUND: /* keep if all fields are primitive*/
	    for(j=0;j<listlength(tsym->subnodes);j++) {
		Symbol* field = (Symbol*)listget(tsym->subnodes,j);
		ASSERT(field->subclass == NC_FIELD);
		if(tagvlentypes(field->typ.basetype)) tagged = 1;
	    }	  
	    break;
	default: break;/* ignore*/
    }
    if(tagged) tsym->typ.hasvlen = 1;
    return tagged;
}

/* Make sure all typecodes are set if basetype is set*/
static void
filltypecodes(void)
{
    Symbol* sym;
    for(sym=symlist;sym != NULL;sym = sym->next) {    
	if(sym->typ.basetype != NULL && sym->typ.typecode == NC_NAT)
	    sym->typ.typecode = sym->typ.basetype->typ.typecode;
    }
}

static void
processenums(void)
{
    int i,j;
    List* enumids = listnew();
    for(i=0;i<listlength(typdefs);i++) {
	Symbol* sym = (Symbol*)listget(typdefs,i);
	ASSERT(sym->objectclass == NC_TYPE);
	if(sym->subclass != NC_ENUM) continue;
	for(j=0;j<listlength(sym->subnodes);j++) {
	    Symbol* esym = (Symbol*)listget(sym->subnodes,j);
	    ASSERT(esym->subclass == NC_ECONST);
	    listpush(enumids,(void*)esym);
	}
    }	    
    /* Convert enum values to match enum type*/
    for(i=0;i<listlength(typdefs);i++) {
	Symbol* tsym = (Symbol*)listget(typdefs,i);
	ASSERT(tsym->objectclass == NC_TYPE);
	if(tsym->subclass != NC_ENUM) continue;
	for(j=0;j<listlength(tsym->subnodes);j++) {
	    Symbol* esym = (Symbol*)listget(tsym->subnodes,j);
	    NCConstant newec;
	    ASSERT(esym->subclass == NC_ECONST);
	    newec.nctype = esym->typ.typecode;
	    convert1(&esym->typ.econst,&newec);
	    esym->typ.econst = newec;
	}	
    }
}

/* Walk all data lists looking for econst refs
   and convert to point to actual definition
*/
static void
processeconstrefs(void)
{
    int i;
    /* locate all the datalist and walk them recursively */
    for(i=0;i<listlength(attdefs);i++) {
	Symbol* att = (Symbol*)listget(attdefs,i);
	if(att->data != NULL && listlength(att->data) > 0)
	    processeconstrefsR(att->data);
    }
    for(i=0;i<listlength(vardefs);i++) {
	Symbol* var = (Symbol*)listget(vardefs,i);
	if(var->data != NULL && listlength(var->data) > 0)
	    processeconstrefsR(var->data);
    }
}

/* Recursive helper for processeconstrefs */
static void
processeconstrefsR(Datalist* data)
{
    NCConstant* con;
    int i;
    for(i=0,con=data->data;i<data->alloc;i++,con++) {
	if(con->nctype == NC_COMPOUND) {
	    /* Iterate over the sublists */
	    processeconstrefsR(con->value.compoundv);
	} else if(con->nctype == NC_ECONST) {
	    fixeconstref(con);
	}
    }
}

static void
fixeconstref(NCConstant* con)
{
    Symbol* match = NULL;
    Symbol* parent = NULL;
    Symbol* refsym = con->value.enumv;
    List* grpmatches;

    /* Locate all possible matching enum constant definitions */
    List* candidates = findecmatches(refsym->name);
    if(candidates == NULL) {
	semerror(con->lineno,"Undefined enum or enum constant reference: %s",refsym->name);
	return;
    }
    /* One hopes that 99% of the time, the match is unique */
    if(listlength(candidates) == 1) {
	con->value.enumv = (Symbol*)listget(candidates,0);
	goto done;
    }
    /* If this ref has a specified group prefix, then find that group
       and search only within it for matches to the candidates */
    if(refsym->is_prefixed && refsym->prefix != NULL) {
	parent = lookupgroup(refsym->prefix);
	if(parent == NULL) {
	    semerror(con->lineno,"Undefined group reference: ",fullname(refsym));
	    goto done;
	}
	/* Search this group only for matches */
	grpmatches = ecsearchgrp(parent,candidates);
	switch (listlength(grpmatches)) {
	case 0:
	    semerror(con->lineno,"Undefined enum or enum constant reference: ",refsym->name);
	    listfree(grpmatches);
	    goto done;
	case 1:
	    break;
	default:
	    semerror(con->lineno,"Ambiguous enum constant reference: %s", fullname(refsym));
	}
	con->value.enumv = listget(grpmatches,0);
	listfree(grpmatches);
	goto done;
    }
    /* Sigh, we have to search up the tree to see if any of our candidates are there */
    parent = refsym->container;
    assert(parent == NULL || parent->objectclass == NC_GRP);
    while(parent != NULL && match == NULL) {
	grpmatches = ecsearchgrp(parent,candidates);
	switch (listlength(grpmatches)) {
	case 0: break;
	case 1: match = listget(grpmatches,0); break;
	default:
	    semerror(con->lineno,"Ambiguous enum constant reference: %s", fullname(refsym));
	    match = listget(grpmatches,0);
	    break;
	}
	listfree(grpmatches);
    }
    if(match != NULL) {
	con->value.enumv = match;
	goto done;
    }
    /* Not unique and not in the parent tree, so complains and pick the first candidate */
    semerror(con->lineno,"Ambiguous enum constant reference: %s", fullname(refsym));
    con->value.enumv = (Symbol*)listget(candidates,0);
done:
    listfree(candidates);
}

/*
Locate enums whose name is a prefix of ident
and contains the suffix as an enum const
and capture that enum constant.
*/
static List*
findecmatches(char* ident)
{
    List* matches = listnew();
    int i;

    for(i=0;i<listlength(typdefs);i++) {
	int len;
	Symbol* ec;
	Symbol* en = (Symbol*)listget(typdefs,i);
	if(en->subclass != NC_ENUM)
	    continue;
        /* First, assume that the ident is the econst name only */
	ec = checkeconst(en,ident);
	if(ec != NULL)
	    listpush(matches,ec);
	/* Second, do the prefix check */	
	len = strlen(en->name);
	if(strncmp(ident,en->name,len) == 0) {
		Symbol *ec;
		/* Find the matching ec constant, if any */
	    if(*(ident+len) != '.') continue;
	    ec = checkeconst(en,ident+len+1); /* +1 for the dot */
	    if(ec != NULL)
		listpush(matches,ec);
	}
    }
    if(listlength(matches) == 0) {
	listfree(matches);
        matches = NULL;
    }
    return matches;    
}

static List*
ecsearchgrp(Symbol* grp, List* candidates)
{
    List* matches = listnew();
    int i,j;
    /* do the intersection of grp subnodes and candidates */
    for(i=0;i<listlength(grp->subnodes);i++) {
	Symbol* sub= (Symbol*)listget(grp->subnodes,i);
	if(sub->subclass != NC_ENUM)
	    continue;
	for(j=0;j<listlength(candidates);j++) {
	    Symbol* ec = (Symbol*)listget(candidates,j);
	    if(ec->container == sub)
		listpush(matches,ec);
	}
    }
    if(listlength(matches) == 0) {
        listfree(matches);
	matches = NULL;
    }
    return matches;
}

static Symbol*
checkeconst(Symbol* en, const char* refname)
{
    int i;
    for(i=0;i<listlength(en->subnodes);i++) {
	Symbol* ec = (Symbol*)listget(en->subnodes,i);
	if(strcmp(ec->name,refname) == 0)
	    return ec;
    }
    return NULL;
}


/* Compute type sizes and compound offsets*/
void
computesize(Symbol* tsym)
{
    int i;
    int offset = 0;
    unsigned long totaldimsize;
    if(tsym->touched) return;
    tsym->touched=1;
    switch (tsym->subclass) {
        case NC_VLEN: /* actually two sizes for vlen*/
	    computesize(tsym->typ.basetype); /* first size*/
	    tsym->typ.size = ncsize(tsym->typ.typecode);
	    tsym->typ.alignment = nctypealignment(tsym->typ.typecode);
	    tsym->typ.nelems = 1; /* always a single compound datalist */
	    break;
	case NC_PRIM:
	    tsym->typ.size = ncsize(tsym->typ.typecode);
	    tsym->typ.alignment = nctypealignment(tsym->typ.typecode);
	    tsym->typ.nelems = 1;
	    break;
	case NC_OPAQUE:
	    /* size and alignment already assigned*/
	    tsym->typ.nelems = 1;
	    break;
	case NC_ENUM:
	    computesize(tsym->typ.basetype); /* first size*/
	    tsym->typ.size = tsym->typ.basetype->typ.size;
	    tsym->typ.alignment = tsym->typ.basetype->typ.alignment;
	    tsym->typ.nelems = 1;
	    break;
	case NC_COMPOUND: /* keep if all fields are primitive*/
	    /* First, compute recursively, the size and alignment of fields*/
	    for(i=0;i<listlength(tsym->subnodes);i++) {
		Symbol* field = (Symbol*)listget(tsym->subnodes,i);
		ASSERT(field->subclass == NC_FIELD);
		computesize(field);
		/* alignment of struct is same as alignment of first field*/
		if(i==0) tsym->typ.alignment = field->typ.alignment;
	    }	  
	    /* now compute the size of the compound based on*/
	    /* what user specified*/
	    offset = 0;
	    for(i=0;i<listlength(tsym->subnodes);i++) {
		Symbol* field = (Symbol*)listget(tsym->subnodes,i);
		/* only support 'c' alignment for now*/
		int alignment = field->typ.alignment;
		offset += getpadding(offset,alignment);
		field->typ.offset = offset;
		offset += field->typ.size;
	    }
	    tsym->typ.size = offset;
	    break;
        case NC_FIELD: /* Compute size assume no unlimited dimensions*/
	    if(tsym->typ.dimset.ndims > 0) {
	        computesize(tsym->typ.basetype);
	        totaldimsize = crossproduct(&tsym->typ.dimset,0,0);
	        tsym->typ.size = tsym->typ.basetype->typ.size * totaldimsize;
	        tsym->typ.alignment = tsym->typ.basetype->typ.alignment;
	        tsym->typ.nelems = 1;
	    } else {
	        tsym->typ.size = tsym->typ.basetype->typ.size;
	        tsym->typ.alignment = tsym->typ.basetype->typ.alignment;
	        tsym->typ.nelems = tsym->typ.basetype->typ.nelems;
	    }
	    break;
	default:
	    PANIC1("computesize: unexpected type class: %d",tsym->subclass);
	    break;
    }
}

void
processvars(void)
{
    int i,j;
    for(i=0;i<listlength(vardefs);i++) {
	Symbol* vsym = (Symbol*)listget(vardefs,i);
	Symbol* basetype = vsym->typ.basetype;
	/* fill in the typecode*/
	vsym->typ.typecode = basetype->typ.typecode;
	/* validate uses of NIL */
        validateNIL(vsym);
	for(j=0;j<vsym->typ.dimset.ndims;j++) {
	    /* validate the dimensions*/
            /* UNLIMITED must only be in first place if using classic */
	    if(vsym->typ.dimset.dimsyms[j]->dim.declsize == NC_UNLIMITED) {
	        if(usingclassic && j != 0)
		    semerror(vsym->lineno,"Variable: %s: UNLIMITED must be in first dimension only",fullname(vsym));
	    }
	}	
    }
}

static void
processtypesizes(void)
{
    int i;
    /* use touch flag to avoid circularity*/
    for(i=0;i<listlength(typdefs);i++) {
	Symbol* tsym = (Symbol*)listget(typdefs,i);
	tsym->touched = 0;
    }
    for(i=0;i<listlength(typdefs);i++) {
	Symbol* tsym = (Symbol*)listget(typdefs,i);
	computesize(tsym); /* this will recurse*/
    }
}

static void
processattributes(void)
{
    int i,j;
    /* process global attributes*/
    for(i=0;i<listlength(gattdefs);i++) {
	Symbol* asym = (Symbol*)listget(gattdefs,i);
	if(asym->typ.basetype == NULL) inferattributetype(asym);
        /* fill in the typecode*/
	asym->typ.typecode = asym->typ.basetype->typ.typecode;
	if(asym->data->length == 0) {
	    /* If the attribute has a zero length, then default it;
               note that it must be of type NC_CHAR */
	    if(asym->typ.typecode != NC_CHAR)
	        semerror(asym->lineno,"Empty datalist can only be assigned to attributes of type char",fullname(asym));
	    asym->data = builddatalist(1);
	    emptystringconst(asym->lineno,&asym->data->data[asym->data->length]);
	}
	validateNIL(asym);
    }
    /* process per variable attributes*/
    for(i=0;i<listlength(attdefs);i++) {
	Symbol* asym = (Symbol*)listget(attdefs,i);
	/* If no basetype is specified, then try to infer it;
           the exception is _Fillvalue, whose type is that of the
           containing variable.
        */
        if(strcmp(asym->name,specialname(_FILLVALUE_FLAG)) == 0) {
	    /* This is _Fillvalue */
	    asym->typ.basetype = asym->att.var->typ.basetype; /* its basetype is same as its var*/
	    /* put the datalist into the specials structure */
	    if(asym->data == NULL) {
		/* Generate a default fill value */
	        asym->data = getfiller(asym->typ.basetype);
	    }
	    asym->att.var->var.special._Fillvalue = asym->data;
	} else if(asym->typ.basetype == NULL) {
	    inferattributetype(asym);
	}
	/* fill in the typecode*/
	asym->typ.typecode = asym->typ.basetype->typ.typecode;
	if(asym->data->length == 0) {
	    /* If the attribute has a zero length, and is char type, then default it */
	    if(asym->typ.typecode != NC_CHAR)
	        semerror(asym->lineno,"Empty datalist can only be assigned to attributes of type char",fullname(asym));
	    asym->data = builddatalist(1);
	    emptystringconst(asym->lineno,&asym->data->data[asym->data->length]);
	}
	validateNIL(asym);
    }
    /* collect per-variable attributes per variable*/
    for(i=0;i<listlength(vardefs);i++) {
	Symbol* vsym = (Symbol*)listget(vardefs,i);
	List* list = listnew();
        for(j=0;j<listlength(attdefs);j++) {
	    Symbol* asym = (Symbol*)listget(attdefs,j);
	    ASSERT(asym->att.var != NULL);
	    if(asym->att.var != vsym) continue;	    
            listpush(list,(void*)asym);
	}
	vsym->var.attributes = list;
    }
}

/*
 Look at the first primitive value of the
 attribute's datalist to infer the type of the attribute.
 There is a potential ambiguity when that value is a string.
 Is the attribute type NC_CHAR or NC_STRING?
 The answer is we always assume it is NC_CHAR in order to
 be back compatible with ncgen.
*/

static nc_type
inferattributetype1(Datasrc* src)
{
    nc_type result = NC_NAT;
    /* Recurse down any enclosing compound markers to find first non-fill "primitive"*/
    while(result == NC_NAT && srcmore(src)) {
	if(issublist(src)) {
	    srcpush(src);
	    result = inferattributetype1(src);
	    srcpop(src);
	} else {	
	    NCConstant* con = srcnext(src);
	    if(isprimplus(con->nctype)) result = con->nctype;
	    /* else keep looking*/
	}
    }
    return result;
}

static void
inferattributetype(Symbol* asym)
{
    Datalist* datalist;
    Datasrc* src;
    nc_type nctype;
    ASSERT(asym->data != NULL);
    datalist = asym->data;
    if(datalist->length == 0) {
        /* Default for zero length attributes */
	asym->typ.basetype = basetypefor(NC_CHAR);
	return;
    }
    src = datalist2src(datalist);
    nctype = inferattributetype1(src);    
    freedatasrc(src);
    /* get the corresponding primitive type built-in symbol*/
    /* special case for string*/
    if(nctype == NC_STRING)
        asym->typ.basetype = basetypefor(NC_CHAR);
    else if(usingclassic) {
        /* If we are in classic mode, then restrict the inferred type
           to the classic types */
	switch (nctype) {
	case NC_UBYTE:
	    nctype = NC_SHORT;
	    break;	
	case NC_USHORT:
	case NC_UINT:
	case NC_INT64:
	case NC_UINT64:
	case NC_OPAQUE:
	case NC_ENUM:
	    nctype = NC_INT;
	    break;
	default: /* leave as is */
	    break;
	}
	asym->typ.basetype = basetypefor(nctype);
    } else
	asym->typ.basetype = basetypefor(nctype);
}

#ifdef USE_NETCDF4
/* recursive helper for validataNIL */
static void
validateNILr(Datalist* src)
{
    int i;
    for(i=0;i<src->length;i++) {
	NCConstant* con = datalistith(src,i);
	if(isnilconst(con))
            semerror(con->lineno,"NIL data can only be assigned to variables or attributes of type string");
	else if(islistconst(con)) /* recurse */
	    validateNILr(con->value.compoundv);
    }
}
#endif

static void
validateNIL(Symbol* sym)
{
#ifdef USE_NETCDF4
    Datalist* datalist = sym->data;

    if(sym->data == NULL || datalist->length == 0) return;
    if(sym->typ.typecode == NC_STRING) return;
    validateNILr(sym->data);
#endif
}


/* Find name within group structure*/
Symbol*
lookupgroup(List* prefix)
{
#ifdef USE_NETCDF4
    if(prefix == NULL || listlength(prefix) == 0)
	return rootgroup;
    else
	return (Symbol*)listtop(prefix);
#else
    return rootgroup;
#endif
}

/* Find name within given group*/
Symbol*
lookupingroup(nc_class objectclass, char* name, Symbol* grp)
{
    int i;
    if(name == NULL) return NULL;
    if(grp == NULL) grp = rootgroup;
dumpgroup(grp);
    for(i=0;i<listlength(grp->subnodes);i++) {
	Symbol* sym = (Symbol*)listget(grp->subnodes,i);
	if(sym->ref.is_ref) continue;
	if(sym->objectclass != objectclass) continue;
	if(strcmp(sym->name,name)!=0) continue;
	return sym;
    }
    return NULL;
}

/* Find symbol within group structure*/
Symbol*
lookup(nc_class objectclass, Symbol* pattern)
{
    Symbol* grp;
    if(pattern == NULL) return NULL;
    grp = lookupgroup(pattern->prefix);
    if(grp == NULL) return NULL;
    return lookupingroup(objectclass,pattern->name,grp);
}


/* return internal size for values of specified netCDF type */
size_t
nctypesize(
     nc_type type)			/* netCDF type code */
{
    switch (type) {
      case NC_BYTE: return sizeof(char);
      case NC_CHAR: return sizeof(char);
      case NC_SHORT: return sizeof(short);
      case NC_INT: return sizeof(int);
      case NC_FLOAT: return sizeof(float);
      case NC_DOUBLE: return sizeof(double);
      case NC_UBYTE: return sizeof(unsigned char);
      case NC_USHORT: return sizeof(unsigned short);
      case NC_UINT: return sizeof(unsigned int);
      case NC_INT64: return sizeof(long long);
      case NC_UINT64: return sizeof(unsigned long long);
      case NC_STRING: return sizeof(char*);
      default:
	PANIC("nctypesize: bad type code");
    }
    return 0;
}

static int
sqContains(List* seq, Symbol* sym)
{
    int i;
    if(seq == NULL) return 0;
    for(i=0;i<listlength(seq);i++) {
        Symbol* sub = (Symbol*)listget(seq,i);
	if(sub == sym) return 1;
    }
    return 0;
}

static void
checkconsistency(void)
{
    int i;
    for(i=0;i<listlength(grpdefs);i++) {
	Symbol* sym = (Symbol*)listget(grpdefs,i);
	if(sym == rootgroup) {
	    if(sym->container != NULL)
	        PANIC("rootgroup has a container");
	} else if(sym->container == NULL && sym != rootgroup)
	    PANIC1("symbol with no container: %s",sym->name);
	else if(sym->container->ref.is_ref != 0)
	    PANIC1("group with reference container: %s",sym->name);
	else if(sym != rootgroup && !sqContains(sym->container->subnodes,sym))
	    PANIC1("group not in container: %s",sym->name);
	if(sym->subnodes == NULL)
	    PANIC1("group with null subnodes: %s",sym->name);
    }
    for(i=0;i<listlength(typdefs);i++) {
	Symbol* sym = (Symbol*)listget(typdefs,i);
        if(!sqContains(sym->container->subnodes,sym))
	    PANIC1("type not in container: %s",sym->name);
    }
    for(i=0;i<listlength(dimdefs);i++) {
	Symbol* sym = (Symbol*)listget(dimdefs,i);
        if(!sqContains(sym->container->subnodes,sym))
	    PANIC1("dimension not in container: %s",sym->name);
    }
    for(i=0;i<listlength(vardefs);i++) {
	Symbol* sym = (Symbol*)listget(vardefs,i);
        if(!sqContains(sym->container->subnodes,sym))
	    PANIC1("variable not in container: %s",sym->name);
	if(!(isprimplus(sym->typ.typecode)
	     || sqContains(typdefs,sym->typ.basetype)))
	    PANIC1("variable with undefined type: %s",sym->name);
    }
}

static void
computeunlimitedsizes(Dimset* dimset, int dimindex, Datalist* data, int ischar)
{
    int i;
    size_t xproduct, unlimsize;
    int nextunlim,lastunlim;
    Symbol* thisunlim = dimset->dimsyms[dimindex];
    size_t length;
    
    ASSERT(thisunlim->dim.isunlimited);
    nextunlim = findunlimited(dimset,dimindex+1);
    lastunlim = (nextunlim == dimset->ndims);

    xproduct = crossproduct(dimset,dimindex+1,nextunlim);

    if(!lastunlim) {
	/* Compute candiate size of this unlimited */
        length = data->length;
	unlimsize = length / xproduct;
	if(length % xproduct != 0)
	    unlimsize++; /* => fill requires at some point */
#ifdef GENDEBUG2
fprintf(stderr,"unlimsize: dim=%s declsize=%lu xproduct=%lu newsize=%lu\n",
thisunlim->name,
(unsigned long)thisunlim->dim.declsize,
(unsigned long)xproduct,
(unsigned long)unlimsize);
#endif
	if(thisunlim->dim.declsize < unlimsize) /* want max length of the unlimited*/
            thisunlim->dim.declsize = unlimsize;
        /*!lastunlim => data is list of sublists, recurse on each sublist*/
	for(i=0;i<data->length;i++) {
	    NCConstant* con = data->data+i;
	    ASSERT(con->nctype == NC_COMPOUND);
	    computeunlimitedsizes(dimset,nextunlim,con->value.compoundv,ischar);
	}
    } else {			/* lastunlim */
	if(ischar) {
	    /* Char case requires special computations;
	       compute total number of characters */
	    length = 0;
	    for(i=0;i<data->length;i++) {
		NCConstant* con = &data->data[i];
		switch (con->nctype) {
	        case NC_CHAR: case NC_BYTE: case NC_UBYTE:
		    length++;
		    break;
		case NC_STRING:
		    length += con->value.stringv.len;
	            break;
		case NC_COMPOUND:
		    semwarn(datalistline(data),"Expected character constant, found {...}");
		    break;
		default:
		    semwarn(datalistline(data),"Illegal character constant: %d",con->nctype);
	        }
	    }
	} else { /* Data list should be a list of simple non-char constants */
   	    length = data->length;
	}
	unlimsize = length / xproduct;
	if(length % xproduct != 0)
	    unlimsize++; /* => fill requires at some point */
#ifdef GENDEBUG2
fprintf(stderr,"unlimsize: dim=%s declsize=%lu xproduct=%lu newsize=%lu\n",
thisunlim->name,
(unsigned long)thisunlim->dim.declsize,
(unsigned long)xproduct,
(unsigned long)unlimsize);
#endif
	if(thisunlim->dim.declsize < unlimsize) /* want max length of the unlimited*/
            thisunlim->dim.declsize = unlimsize;
    }
}

static void
processunlimiteddims(void)
{
    int i;
    /* Set all unlimited dims to size 0; */
    for(i=0;i<listlength(dimdefs);i++) {
	Symbol* dim = (Symbol*)listget(dimdefs,i);
	if(dim->dim.isunlimited)
	    dim->dim.declsize = 0;
    }
    /* Walk all variables */
    for(i=0;i<listlength(vardefs);i++) {
	Symbol* var = (Symbol*)listget(vardefs,i);
	int first,ischar;
	Dimset* dimset = &var->typ.dimset;
	if(dimset->ndims == 0) continue; /* ignore scalars */
	if(var->data == NULL) continue; /* no data list to walk */
	ischar = (var->typ.basetype->typ.typecode == NC_CHAR);
	first = findunlimited(dimset,0);
	if(first == dimset->ndims) continue; /* no unlimited dims */
	if(first == 0) {
	    computeunlimitedsizes(dimset,first,var->data,ischar);
	} else {
	    for(i=0;i<var->data->length;i++) {
	        NCConstant* con = var->data->data+i;
	        if(con->nctype != NC_COMPOUND)
		    semerror(con->lineno,"UNLIMITED dimension (other than first) must be enclosed in {}");
		else
	            computeunlimitedsizes(dimset,first,con->value.compoundv,ischar);
	    }
	}
    }
#ifdef GENDEBUG1
    /* print unlimited dim size */
    if(listlength(dimdefs) == 0)
        fprintf(stderr,"unlimited: no unlimited dimensions\n");
    else for(i=0;i<listlength(dimdefs);i++) {
	Symbol* dim = (Symbol*)listget(dimdefs,i);
	if(dim->dim.isunlimited)
	    fprintf(stderr,"unlimited: %s = %lu\n",
		    dim->name,
	            (unsigned long)dim->dim.declsize);
    }
#endif
}
