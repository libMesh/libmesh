/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Id: ncgen.y,v 1.42 2010/05/18 21:32:46 dmh Exp $
 *********************************************************************/

/* yacc source for "ncgen", a netCDL parser and netCDF generator */

%error-verbose

%{
/*
static char SccsId[] = "$Id: ncgen.y,v 1.42 2010/05/18 21:32:46 dmh Exp $";
*/
#include        "includes.h"
#include        "offsets.h"

/* Following are in ncdump (for now)*/
/* Need some (unused) definitions to get it to compile */
#define ncatt_t void*
#define ncvar_t void
#include        "nctime.h"

/* parser controls */
#define YY_NO_INPUT 1

/* True if string a equals string b*/
#define STREQ(a, b)     (*(a) == *(b) && strcmp((a), (b)) == 0)
#define VLENSIZE  (sizeof(nc_vlen_t))
#define MAXFLOATDIM 4294967295.0

/* mnemonic */
typedef enum Attrkind {ATTRVAR, ATTRGLOBAL, DONTKNOW} Attrkind;

typedef nc_vlen_t vlen_t;

/* We retain the old representation of the symbol list
   as a linked list.
*/
Symbol* symlist;

/* Track rootgroup separately*/
Symbol* rootgroup;

/* Track the group sequence */
static List* groupstack;

/* Provide a separate sequence for accumulating values
   during the parse.
*/
static List* stack;

/* track homogeneity of types for data lists*/
static nc_type consttype;

/* Misc. */
static int stackbase;
static int stacklen;
static int count;
static int opaqueid; /* counter for opaque constants*/
static int arrayuid; /* counter for pseudo-array types*/

char* primtypenames[PRIMNO] = {
"nat",
"byte", "char", "short",
"int", "float", "double",
"ubyte", "ushort", "uint",
"int64", "uint64",
"string"
};

/*Defined in ncgen.l*/
extern int lineno;              /* line number for error messages */
extern Bytebuffer* lextext;           /* name or string with escapes removed */

extern double double_val;       /* last double value read */
extern float float_val;         /* last float value read */
extern long long int64_val;         /* last int64 value read */
extern int int32_val;             /* last int32 value read */
extern short int16_val;         /* last short value read */
extern unsigned long long uint64_val;         /* last int64 value read */
extern unsigned int uint32_val;             /* last int32 value read */
extern unsigned short uint16_val;         /* last short value read */
extern char char_val;           /* last char value read */
extern signed char byte_val;    /* last byte value read */
extern unsigned char ubyte_val;    /* last byte value read */

/* Track definitions of dims, types, attributes, and vars*/
List* grpdefs;
List* dimdefs;
List* attdefs; /* variable-specific attributes*/
List* gattdefs; /* global attributes only*/
List* xattdefs; /* unknown attributes*/
List* typdefs;
List* vardefs;
List* condefs; /* non-dimension constants used in type defs*/
List* tmp;

/* Forward */
static Constant makeconstdata(nc_type);
static Constant evaluate(Symbol* fcn, Datalist* arglist);
static Constant makeenumconst(Symbol*);
static void addtogroup(Symbol*);
static Symbol* currentgroup(void);
static Symbol* createrootgroup(void);
static Symbol* creategroup(Symbol*);
static int dupobjectcheck(nc_class,Symbol*);
static void setpathcurrent(Symbol* sym);
static Symbol* makeattribute(Symbol*,Symbol*,Symbol*,Datalist*,Attrkind);
static Symbol* makeprimitivetype(nc_type i);
static Symbol* makespecial(int tag, Symbol* vsym, Symbol* tsym, void* data, int isconst);
static int containsfills(Datalist* list);
static void datalistextend(Datalist* dl, Constant* con);
static void vercheck(int ncid);

int yylex(void);

#ifndef NO_STDARG
static void yyerror(const char *fmt, ...);
#else
static void yyerror(fmt,va_alist) const char* fmt; va_dcl;
#endif

/* Extern */
extern int lex_init(void);

%}

/* DECLARATIONS */

%union {
Symbol* sym;
unsigned long  size; /* allow for zero size to indicate e.g. UNLIMITED*/
long           mark; /* track indices into the sequence*/
int            nctype; /* for tracking attribute list type*/
Datalist*      datalist;
Constant       constant;
}

%token <sym>
        NC_UNLIMITED_K /* keyword for unbounded record dimension */
        CHAR_K      /* keyword for char datatype */
        BYTE_K      /* keyword for byte datatype */
        SHORT_K     /* keyword for short datatype */
        INT_K       /* keyword for int datatype */
        FLOAT_K     /* keyword for float datatype */
        DOUBLE_K    /* keyword for double datatype */
        UBYTE_K      /* keyword for unsigned byte datatype */
        USHORT_K     /* keyword for unsigned short datatype */
        UINT_K       /* keyword for unsigned int datatype */
        INT64_K       /* keyword for long long datatype */
        UINT64_K       /* keyword for unsigned long long datatype */
        IDENT       /* name for a dimension, variable, or attribute */
        TERMSTRING  /* terminal string */
        CHAR_CONST  /* char constant (not ever generated by ncgen.l) */
        BYTE_CONST  /* byte constant */
        SHORT_CONST /* short constant */
        INT_CONST   /* int constant */
        INT64_CONST   /* long long constant */
        UBYTE_CONST  /* unsigned byte constant */
        USHORT_CONST /* unsigned short constant */
        UINT_CONST   /* unsigned long long constant */
        UINT64_CONST   /* unsigned int constant */
        FLOAT_CONST /* float constant */
        DOUBLE_CONST /* double constant */
        DIMENSIONS  /* keyword starting dimensions section, if any */
        VARIABLES   /* keyword starting variables section, if any */
        NETCDF      /* keyword declaring netcdf name */
        DATA        /* keyword starting data section, if any */
        TYPES
	COMPOUND
        ENUM
        OPAQUE
        OPAQUESTRING    /* 0x<even number of hexdigits> */
        GROUP
	PATH            /* / or (/IDENT)+ */
	FILLMARKER	/* "_" as opposed to the attribute */
        _FILLVALUE
        _FORMAT
        _STORAGE
        _CHUNKSIZES
        _DEFLATELEVEL
        _SHUFFLE
        _ENDIANNESS
        _NOFILL
        _FLETCHER32
	DATASETID	

%type <sym> ident typename primtype dimd varspec
	    attrdecl enumid path dimref fielddim fieldspec
%type <sym> typeref
%type <sym> varref
%type <sym> type_var_ref
%type <mark> enumidlist fieldlist fields varlist dimspec dimlist field
	     fielddimspec fielddimlist
%type <constant> dataitem constdata constint conststring constbool
%type <constant> simpleconstant function 
%type <datalist> datalist intlist datalist1 datalist0 arglist


%start  ncdesc /* start symbol for grammar */

%%

/* RULES */

ncdesc: NETCDF
	DATASETID
        rootgroup
        {if (error_count > 0) YYABORT;}
        ;

rootgroup: '{'
           groupbody
           subgrouplist
           '}';

/* 2/3/08 - Allow group body with only attributes. (H/T John Storrs). */
groupbody:
		attrdecllist
                typesection     /* Type definitions */
                dimsection      /* dimension declarations */
                vasection       /* variable and attribute declarations */
                datasection     /* data for variables within the group */
                ;

subgrouplist: /*empty*/ | subgrouplist namedgroup;

namedgroup: GROUP ident '{'
            {
		Symbol* id = $2;
                markcdf4("Group specification");
		if(creategroup(id) == NULL) 
                    yyerror("duplicate group declaration within parent group for %s",
                                id->name);
            }
            groupbody
            subgrouplist
            {listpop(groupstack);}
            '}'
	    attrdecllist
	    ;
      
typesection:    /* empty */
                | TYPES {}
		| TYPES typedecls
			{markcdf4("Type specification");}
                ;

typedecls: type_or_attr_decl | typedecls type_or_attr_decl ;

typename: ident
	    { /* Use when defining a type */
              $1->objectclass = NC_TYPE;
              if(dupobjectcheck(NC_TYPE,$1))
                    yyerror("duplicate type declaration for %s",
                            $1->name);
              listpush(typdefs,(void*)$1);
	    }
	  ;

type_or_attr_decl: typedecl {} | attrdecl ';' {} ;

typedecl:
	  enumdecl optsemicolon
	| compounddecl optsemicolon
	| vlendecl optsemicolon
	| opaquedecl optsemicolon
	;

optsemicolon: /*empty*/ | ';' ;


enumdecl: primtype ENUM typename
          '{' enumidlist '}'
              {
		int i;
                addtogroup($3); /* sets prefix*/
                $3->objectclass=NC_TYPE;
                $3->subclass=NC_ENUM;
                $3->typ.basetype=$1;
                $3->typ.size = $1->typ.size;
                $3->typ.alignment = $1->typ.alignment;
                stackbase=$5;
                stacklen=listlength(stack);
                $3->subnodes = listnew();
                /* Variety of field fixups*/
		/* 1. add in the enum values*/
		/* 2. make this type be their container*/
		/* 3. make constant names visible in the group*/
		/* 4. set field basetype to be same as enum basetype*/
                for(i=stackbase;i<stacklen;i++) {
                   Symbol* eid = (Symbol*)listget(stack,i);
		   assert(eid->subclass == NC_ECONST);
		   addtogroup(eid);
                   listpush($3->subnodes,(void*)eid);
                   eid->container = $3;
		   eid->typ.basetype = $3->typ.basetype;
                }               
                listsetlength(stack,stackbase);/* remove stack nodes*/
              }
          ;

enumidlist:   enumid
		{$$=listlength(stack); listpush(stack,(void*)$1);}
	    | enumidlist ',' enumid
		{
		    int i;
		    $$=$1;
		    /* check for duplicates*/
		    stackbase=$1;
		    stacklen=listlength(stack);
		    for(i=stackbase;i<stacklen;i++) {
		      Symbol* elem = (Symbol*)listget(stack,i);
		      if(strcmp($3->name,elem->name)==0)
  	                yyerror("duplicate enum declaration for %s",
        	                 elem->name);
		    }    	    
		    listpush(stack,(void*)$3);
		}
	    ;

enumid: ident '=' constdata
        {
            $1->objectclass=NC_TYPE;
            $1->subclass=NC_ECONST;
            $1->typ.econst=$3;
	    $$=$1;
        }
        ;

opaquedecl: OPAQUE '(' INT_CONST ')' typename
                {
		    vercheck(NC_OPAQUE);
                    addtogroup($5); /*sets prefix*/
                    $5->objectclass=NC_TYPE;
                    $5->subclass=NC_OPAQUE;
                    $5->typ.typecode=NC_OPAQUE;
                    $5->typ.size=int32_val;
                    $5->typ.alignment=nctypealignment(NC_OPAQUE);
                }
            ;

vlendecl: typeref '(' '*' ')' typename
                {
                    Symbol* basetype = $1;
		    vercheck(NC_VLEN);
                    addtogroup($5); /*sets prefix*/
                    $5->objectclass=NC_TYPE;
                    $5->subclass=NC_VLEN;
                    $5->typ.basetype=basetype;
                    $5->typ.typecode=NC_VLEN;
                    $5->typ.size=VLENSIZE;
                    $5->typ.alignment=nctypealignment(NC_VLEN);
                }
          ;

compounddecl: COMPOUND typename '{' fields '}'
          {
	    int i,j;
	    vercheck(NC_COMPOUND);
            addtogroup($2);
	    /* check for duplicate field names*/
	    stackbase=$4;
	    stacklen=listlength(stack);
	    for(i=stackbase;i<stacklen;i++) {
	      Symbol* elem1 = (Symbol*)listget(stack,i);
	      for(j=i+1;j<stacklen;j++) {
	          Symbol* elem2 = (Symbol*)listget(stack,j);
	          if(strcmp(elem1->name,elem2->name)==0) {
	            yyerror("duplicate field declaration for %s",elem1->name);
		  }
	      }
	    }
	    $2->objectclass=NC_TYPE;
            $2->subclass=NC_COMPOUND;
            $2->typ.basetype=NULL;
            $2->typ.typecode=NC_COMPOUND;
	    $2->subnodes = listnew();
	    /* Add in the fields*/
	    for(i=stackbase;i<stacklen;i++) {
	        Symbol* fsym = (Symbol*)listget(stack,i);
		fsym->container = $2;
 	        listpush($2->subnodes,(void*)fsym);
	    }    	    
	    listsetlength(stack,stackbase);/* remove stack nodes*/
          }
            ;


fields:   field ';' {$$=$1;}
	   | fields field ';' {$$=$1;}
	   ;

field: typeref fieldlist
        {
	    int i;
	    $$=$2;
	    stackbase=$2;
	    stacklen=listlength(stack);
	    /* process each field in the fieldlist*/
            for(i=stackbase;i<stacklen;i++) {
                Symbol* f = (Symbol*)listget(stack,i);
		f->typ.basetype = $1;
            }
        }
        ;

primtype:         CHAR_K  { $$ = primsymbols[NC_CHAR]; }
                | BYTE_K  { $$ = primsymbols[NC_BYTE]; }
                | SHORT_K { $$ = primsymbols[NC_SHORT]; }
                | INT_K   { $$ = primsymbols[NC_INT]; }
                | FLOAT_K { $$ = primsymbols[NC_FLOAT]; }
                | DOUBLE_K{ $$ = primsymbols[NC_DOUBLE]; }
                | UBYTE_K  { vercheck(NC_UBYTE); $$ = primsymbols[NC_UBYTE]; }
                | USHORT_K { vercheck(NC_USHORT); $$ = primsymbols[NC_USHORT]; }
                | UINT_K   { vercheck(NC_UINT); $$ = primsymbols[NC_UINT]; }
                | INT64_K   { vercheck(NC_INT64); $$ = primsymbols[NC_INT64]; }
                | UINT64_K   { vercheck(NC_UINT64); $$ = primsymbols[NC_UINT64]; }
                ;

dimsection:     /* empty */
                | DIMENSIONS {}
		| DIMENSIONS dimdecls {}
                ;

dimdecls:       dim_or_attr_decl ';'
                | dimdecls dim_or_attr_decl ';'
                ;

dim_or_attr_decl: dimdeclist {} | attrdecl {} ;

dimdeclist:     dimdecl
                | dimdeclist ',' dimdecl
                ;

dimdecl:
	  dimd '=' UINT_CONST
              {
		$1->dim.declsize = (size_t)uint32_val;
#ifdef DEBUG1
fprintf(stderr,"dimension: %s = %lu\n",$1->name,(unsigned long)$1->dim.declsize);
#endif
	      }
	| dimd '=' INT_CONST
              {
		if(int32_val <= 0) {
		    derror("dimension size must be positive");
		    YYABORT;
		}
		$1->dim.declsize = (size_t)int32_val;
#ifdef DEBUG1
fprintf(stderr,"dimension: %s = %lu\n",$1->name,(unsigned long)$1->dim.declsize);
#endif
	      }
        | dimd '=' DOUBLE_CONST
                   { /* for rare case where 2^31 < dimsize < 2^32 */
                       if (double_val <= 0)
                         yyerror("dimension length must be positive");
                       if (double_val > MAXFLOATDIM)
                         yyerror("dimension too large");
                       if (double_val - (size_t) double_val > 0)
                         yyerror("dimension length must be an integer");
                       $1->dim.declsize = (size_t)double_val;
#ifdef DEBUG1
fprintf(stderr,"dimension: %s = %lu\n",$1->name,(unsigned long)$1->dim.declsize);
#endif
                   }
        | dimd '=' NC_UNLIMITED_K
                   {
		        $1->dim.declsize = NC_UNLIMITED;
		        $1->dim.isunlimited = 1;
#ifdef DEBUG1
fprintf(stderr,"dimension: %s = UNLIMITED\n",$1->name);
#endif
		   }
                ;

dimd:           ident
                   { 
                     $1->objectclass=NC_DIM;
                     if(dupobjectcheck(NC_DIM,$1))
                        yyerror( "Duplicate dimension declaration for %s",
                                $1->name);
		     addtogroup($1);
		     $$=$1;
		     listpush(dimdefs,(void*)$1);
                   }
                ;

vasection:      /* empty */
                | VARIABLES {}
                | VARIABLES vadecls {}
                ;

vadecls:        vadecl_or_attr ';'
                | vadecls vadecl_or_attr ';'
                ;

vadecl_or_attr: vardecl {} | attrdecl {} ;

vardecl:        typeref varlist
		{
		    int i;
		    stackbase=$2;
		    stacklen=listlength(stack);
		    /* process each variable in the varlist*/
	            for(i=stackbase;i<stacklen;i++) {
	                Symbol* sym = (Symbol*)listget(stack,i);
			sym->objectclass = NC_VAR;
		        if(dupobjectcheck(NC_VAR,sym)) {
                            yyerror("Duplicate variable declaration for %s",
                                    sym->name);
			} else {
		  	    sym->typ.basetype = $1;
	                    addtogroup(sym);
		            listpush(vardefs,(void*)sym);
			}
		    }
		    listsetlength(stack,stackbase);/* remove stack nodes*/
		}
                ;

varlist:      varspec
	        {$$=listlength(stack);
                 listpush(stack,(void*)$1);
		}
            | varlist ',' varspec
	        {$$=$1; listpush(stack,(void*)$3);}
            ;

varspec:        ident dimspec
                    {
		    int i;
		    Dimset dimset;
		    stacklen=listlength(stack);
		    stackbase=$2;
		    count = stacklen - stackbase;
		    if(count >= NC_MAX_VAR_DIMS) {
			yyerror("%s has too many dimensions",$1->name);
			count = NC_MAX_VAR_DIMS - 1;
			stacklen = stackbase + count;
		    }
  	            dimset.ndims = count;
		    /* extract the actual dimensions*/
		    if(dimset.ndims > 0) {
		        for(i=0;i<count;i++) {
			    Symbol* dsym = (Symbol*)listget(stack,stackbase+i);
			    dimset.dimsyms[i] = dsym;
			}
			$1->typ.dimset = dimset;
		    }
		    $1->typ.basetype = NULL; /* not yet known*/
                    $1->objectclass=NC_VAR;
		    listsetlength(stack,stackbase);/* remove stack nodes*/
		    }
                ;

dimspec:        /* empty */ {$$=listlength(stack);}
                | '(' dimlist ')' {$$=$2;}
                ;

dimlist:        dimref {$$=listlength(stack); listpush(stack,(void*)$1);}
                | dimlist ',' dimref
		    {$$=$1; listpush(stack,(void*)$3);}
                ;

dimref: path
            {Symbol* dimsym = $1;
		dimsym->objectclass = NC_DIM;
		/* Find the actual dimension*/
		dimsym = locate(dimsym);
		if(dimsym == NULL) {
		    derror("Undefined or forward referenced dimension: %s",$1->name);
		    YYABORT;
		}
		$$=dimsym;
	    }
	;

fieldlist:
	  fieldspec
	    {$$=listlength(stack);
             listpush(stack,(void*)$1);
	    }
	| fieldlist ',' fieldspec
	    {$$=$1; listpush(stack,(void*)$3);}
        ;

fieldspec:
	ident fielddimspec
	    {
		int i;
		Dimset dimset;
		stackbase=$2;
		stacklen=listlength(stack);
		count = stacklen - stackbase;
		if(count >= NC_MAX_VAR_DIMS) {
		    yyerror("%s has too many dimensions",$1->name);
		    count = NC_MAX_VAR_DIMS - 1;
		    stacklen = stackbase + count;
		}
  	        dimset.ndims = count;
		if(count > 0) {
		    /* extract the actual dimensions*/
		    for(i=0;i<count;i++) {
		        Symbol* dsym = (Symbol*)listget(stack,stackbase+i);
		        dimset.dimsyms[i] = dsym;
		    }
		    $1->typ.dimset = dimset;
		}
		$1->typ.basetype = NULL; /* not yet known*/
                $1->objectclass=NC_TYPE;
                $1->subclass=NC_FIELD;
		listsetlength(stack,stackbase);/* remove stack nodes*/
		$$ = $1;
	    }
	;

fielddimspec:        /* empty */ {$$=listlength(stack);}
                | '(' fielddimlist ')' {$$=$2;}
                ;

fielddimlist:
	  fielddim {$$=listlength(stack); listpush(stack,(void*)$1);}
	| fielddimlist ',' fielddim
	    {$$=$1; listpush(stack,(void*)$3);}
        ;

fielddim:
	  UINT_CONST
	    {  /* Anonymous integer dimension.
	         Can only occur in type definitions*/
	     char anon[32];
	     sprintf(anon,"const%u",uint32_val);
	     $$ = install(anon);
	     $$->objectclass = NC_DIM;
	     $$->dim.isconstant = 1;
	     $$->dim.declsize = uint32_val;
	    }
	| INT_CONST
	    {  /* Anonymous integer dimension.
	         Can only occur in type definitions*/
	     char anon[32];
	     if(int32_val <= 0) {
		derror("field dimension must be positive");
		YYABORT;
	     }
	     sprintf(anon,"const%d",int32_val);
	     $$ = install(anon);
	     $$->objectclass = NC_DIM;
	     $$->dim.isconstant = 1;
	     $$->dim.declsize = int32_val;
	    }
	;


/* Use this when referencing defined objects */

varref:
	type_var_ref
	    {Symbol* vsym = $1;
		if(vsym->objectclass != NC_VAR) {
		    derror("Undefined or forward referenced variable: %s",vsym->name);
		    YYABORT;
		}
		$$=vsym;
	    }
	  ;

typeref:
	type_var_ref	  
	    {Symbol* tsym = $1;
		if(tsym->objectclass != NC_TYPE) {
		    derror("Undefined or forward referenced type: %s",tsym->name);
		    YYABORT;
		}
		$$=tsym;
	    }
	;

type_var_ref: 
	path 
	    {Symbol* tvsym = $1; Symbol* sym;
		/* disambiguate*/
		tvsym->objectclass = NC_VAR;
		sym = locate(tvsym);
		if(sym == NULL) {
		    tvsym->objectclass = NC_TYPE;
		    sym = locate(tvsym);
		    if(tvsym == NULL) {
		        derror("Undefined or forward referenced name: %s",$1->name);
		        YYABORT;
		    } else tvsym = sym;
		} else tvsym = sym;
		if(tvsym == NULL) {
		    derror("Undefined name: %s",$1->name);
		    YYABORT;
		}
		$$=tvsym;
	    }
	| primtype {$$=$1;}
	;

/* Use this for all attribute decls */
/* Global vs var-specific will be separated in makeattribute */

/* Watch out; this is left recursive */
attrdecllist: /*empty*/ {} | attrdecl ';' attrdecllist {} ;

attrdecl:
	  ':' ident '=' datalist
	    { $$=makeattribute($2,NULL,NULL,$4,ATTRGLOBAL);}
	| typeref type_var_ref ':' ident '=' datalist
	    {Symbol* tsym = $1; Symbol* vsym = $2; Symbol* asym = $4;
		if(vsym->objectclass == NC_VAR) {
		    $$=makeattribute(asym,vsym,tsym,$6,ATTRVAR);
		} else {
		    derror("Doubly typed attribute: %s",asym->name);
		    YYABORT;
		}
	    }
	| type_var_ref ':' ident '=' datalist
	    {Symbol* sym = $1; Symbol* asym = $3;
		if(sym->objectclass == NC_VAR) {
		    $$=makeattribute(asym,sym,NULL,$5,ATTRVAR);
		} else if(sym->objectclass == NC_TYPE) {
		    $$=makeattribute(asym,NULL,sym,$5,ATTRGLOBAL);
		} else {
		    derror("Attribute prefix not a variable or type: %s",asym->name);
		    YYABORT;
		}
	    }
	| type_var_ref ':' _FILLVALUE '=' datalist
	    {$$ = makespecial(_FILLVALUE_FLAG,$1,NULL,(void*)$5,0);}
	| typeref type_var_ref ':' _FILLVALUE '=' datalist
	    {$$ = makespecial(_FILLVALUE_FLAG,$2,$1,(void*)$6,0);}
	| type_var_ref ':' _STORAGE '=' conststring
	    {$$ = makespecial(_STORAGE_FLAG,$1,NULL,(void*)&$5,1);}
	| type_var_ref ':' _CHUNKSIZES '=' intlist
	    {$$ = makespecial(_CHUNKSIZES_FLAG,$1,NULL,(void*)$5,0);}
	| type_var_ref ':' _FLETCHER32 '=' constbool
	    {$$ = makespecial(_FLETCHER32_FLAG,$1,NULL,(void*)&$5,1);}
	| type_var_ref ':' _DEFLATELEVEL '=' constint
	    {$$ = makespecial(_DEFLATE_FLAG,$1,NULL,(void*)&$5,1);}
	| type_var_ref ':' _SHUFFLE '=' constbool
	    {$$ = makespecial(_SHUFFLE_FLAG,$1,NULL,(void*)&$5,1);}
	| type_var_ref ':' _ENDIANNESS '=' conststring
	    {$$ = makespecial(_ENDIAN_FLAG,$1,NULL,(void*)&$5,1);}
	| type_var_ref ':' _NOFILL '=' constbool
	    {$$ = makespecial(_NOFILL_FLAG,$1,NULL,(void*)&$5,1);}
	| ':' _FORMAT '=' conststring
	    {$$ = makespecial(_FORMAT_FLAG,NULL,NULL,(void*)&$4,1);}
	;

path:
	  ident
	    {
	        $$=$1;
                $1->is_ref=1;
                setpathcurrent($1);
	    }
	| PATH
	    {
	        $$=$1;
                $1->is_ref=1;
                $1->is_prefixed=1;
	        /* path is set in ncgen.l*/
	    }
	;

datasection:    /* empty */
                | DATA {}
                | DATA datadecls {}
                ;

datadecls:      datadecl ';'
                | datadecls datadecl ';'
                ;

datadecl:       varref '=' datalist
                   {$1->data = $3;}
                ;
datalist:
	  datalist0 {$$ = $1;}
	| datalist1 {$$ = $1;}
	;

datalist0:
	/*empty*/ {$$ = builddatalist(0);}
	;

datalist1: /* Must have at least 1 element */
	  dataitem {$$ = builddatalist(0); datalistextend($$,&($1));}
	| datalist ',' dataitem
	    {datalistextend($1,&($3)); $$=$1;}
	;

dataitem:
	  constdata {$$=$1;}
	| '{' datalist '}' {$$=builddatasublist($2);}
	;

constdata:
	  simpleconstant      {$$=$1;}
	| OPAQUESTRING	{$$=makeconstdata(NC_OPAQUE);}
	| FILLMARKER	{$$=makeconstdata(NC_FILLVALUE);}
	| path		{$$=makeenumconst($1);}
	| function
	;

function:
	ident '(' arglist ')' {$$=evaluate($1,$3);}
	;

arglist:
	  simpleconstant
	    {$$ = builddatalist(0); datalistextend($$,&($1));}
	| arglist ',' simpleconstant
	    {datalistextend($1,&($3)); $$=$1;}
	;

simpleconstant:
	  CHAR_CONST	{$$=makeconstdata(NC_CHAR);} /* never used apparently*/
	| BYTE_CONST	{$$=makeconstdata(NC_BYTE);}
	| SHORT_CONST	{$$=makeconstdata(NC_SHORT);}
	| INT_CONST	{$$=makeconstdata(NC_INT);}
	| INT64_CONST	{$$=makeconstdata(NC_INT64);}
	| UBYTE_CONST	{$$=makeconstdata(NC_UBYTE);}
	| USHORT_CONST	{$$=makeconstdata(NC_USHORT);}
	| UINT_CONST	{$$=makeconstdata(NC_UINT);}
	| UINT64_CONST	{$$=makeconstdata(NC_UINT64);}
	| FLOAT_CONST	{$$=makeconstdata(NC_FLOAT);}
	| DOUBLE_CONST	{$$=makeconstdata(NC_DOUBLE);}
	| TERMSTRING	{$$=makeconstdata(NC_STRING);}
	;

intlist:
	  constint {$$ = builddatalist(0); datalistextend($$,&($1));}
	| intlist ',' constint {$$=$1; datalistextend($1,&($3));}
	;

constint:
	  INT_CONST
		{$$=makeconstdata(NC_INT);}
	| UINT_CONST
		{$$=makeconstdata(NC_UINT);}
	| INT64_CONST
		{$$=makeconstdata(NC_INT64);}
	| UINT64_CONST
		{$$=makeconstdata(NC_UINT64);}
	;

conststring:
	TERMSTRING	{$$=makeconstdata(NC_STRING);}
	;

constbool:
	  conststring {$$=$1;}
	| constint {$$=$1;}

/* End OF RULES */

/* Push all idents thru here*/
ident:
	IDENT {$$=$1;}
	;

%%

#ifndef NO_STDARG
static void
yyerror(const char *fmt, ...)
#else
static void
yyerror(fmt,va_alist) const char* fmt; va_dcl
#endif
{
    va_list argv;
    vastart(argv,fmt);
    (void)fprintf(stderr,"%s: %s line %d: ", progname, cdlname, lineno);
    vderror(fmt,argv);
}

/* undefine yywrap macro, in case we are using bison instead of yacc */
#ifdef yywrap
#undef yywrap
#endif

static int
ncgwrap(void)                    /* returns 1 on EOF if no more input */
{
    return  1;
}

/* get lexical input routine generated by lex  */
#include "ncgenyy.c"

/* Really should init our data within this file */
void
parse_init(void)
{
    int i;
    opaqueid = 0;
    arrayuid = 0;
    symlist = NULL;
    stack = listnew();
    groupstack = listnew();
    consttype = NC_NAT;
    grpdefs = listnew();
    dimdefs = listnew();
    attdefs = listnew();
    gattdefs = listnew();
    xattdefs = listnew();
    typdefs = listnew();
    vardefs = listnew();
    condefs = listnew();
    tmp = listnew();
    createrootgroup();
    /* Create the primitive types */
    for(i=NC_NAT+1;i<=NC_STRING;i++) {
        primsymbols[i] = makeprimitivetype(i);
    }
    lex_init();
}

static Symbol*
makeprimitivetype(nc_type nctype)
{
    Symbol* sym = install(primtypenames[nctype]);
    sym->objectclass=NC_TYPE;
    sym->subclass=NC_PRIM;
    sym->ncid = nctype;
    sym->typ.typecode = nctype;
    sym->typ.size = ncsize(nctype);
    sym->typ.nelems = 1;
    sym->typ.alignment = nctypealignment(nctype);
    /* Make the basetype circular so we can always ask for it */
    sym->typ.basetype = sym;
    sym->prefix = listnew();
    return sym;
}

/* Symbol table operations for ncgen tool */
/* install sname in symbol table even if it is already there */
Symbol*
install(const char *sname)
{
    Symbol* sp;
    sp = (Symbol*) emalloc (sizeof (struct Symbol));
    memset((void*)sp,0,sizeof(struct Symbol));    
    sp->name = nulldup(sname);
    sp->next = symlist;
    sp->lineno = lineno;
    sp->location = currentgroup();
    sp->container = currentgroup();
    symlist = sp;
    return sp;
}


static Symbol*
currentgroup(void)
{
    if(listlength(groupstack) == 0) return rootgroup;
    return (Symbol*)listtop(groupstack);
}

static Symbol*
createrootgroup(void)
{
    Symbol* gsym = install(ROOTGROUPNAME);
    gsym->objectclass = NC_GRP;
    gsym->container = NULL;
    gsym->subnodes = listnew();
    gsym->grp.is_root = 1;
    gsym->prefix = listnew();
    listpush(grpdefs,(void*)gsym);
    rootgroup = gsym;
    return gsym;
}

static Symbol*
creategroup(Symbol * gsym)
{
    /* See if this group already exists in currentgroup */
    gsym->objectclass = NC_GRP;
    if(dupobjectcheck(NC_GRP,gsym)) {
        derror("Duplicate group name in same scope: %s",gsym->name);
	return NULL;
    }
    addtogroup(gsym);
    gsym->subnodes = listnew();
    listpush(groupstack,(void*)gsym);
    listpush(grpdefs,(void*)gsym);
    return gsym;
}

static Constant
makeconstdata(nc_type nctype)
{
    Constant con = nullconstant;
    consttype = nctype;
    con.nctype = nctype;
    con.lineno = lineno;
    con.filled = 0;
    switch (nctype) {
	case NC_CHAR: con.value.charv = char_val; break;
        case NC_BYTE: con.value.int8v = byte_val; break;
        case NC_SHORT: con.value.int16v = int16_val; break;
        case NC_INT: con.value.int32v = int32_val; break;
        case NC_FLOAT:
	    con.value.floatv = float_val;
	    break;
        case NC_DOUBLE:
	    con.value.doublev = double_val;
	    break;
        case NC_STRING: { /* convert to a set of chars*/
	    size_t len;
	    len = bbLength(lextext);
	    con.value.stringv.len = len;
	    con.value.stringv.stringv = bbDup(lextext);
	    bbClear(lextext);	    
	    }
	    break;

	/* Allow these constants even in netcdf-3 */
        case NC_UBYTE: con.value.uint8v = ubyte_val; break;
        case NC_USHORT: con.value.uint16v = uint16_val; break;
        case NC_UINT: con.value.uint32v = uint32_val; break;
        case NC_INT64: con.value.int64v = int64_val; break;
        case NC_UINT64: con.value.uint64v = uint64_val; break;

#ifdef USE_NETCDF4
	case NC_OPAQUE: {
	    char* s;
	    int len;
	    len = bbLength(lextext);
	    s = (char*)emalloc(len+1);
	    strncpy(s,bbContents(lextext),len);
	    s[len] = '\0';
	    con.value.opaquev.stringv = s;
	    con.value.opaquev.len = len;
	    } break;
#endif

	case NC_FILLVALUE:
	    break; /* no associated value*/
	    
	default:
	    yyerror("Data constant: unexpected NC type: %s",
		    nctypename(nctype));
	    con.value.stringv.stringv = NULL;    
	    con.value.stringv.len = 0;
    }
    return con;
}

static Constant
makeenumconst(Symbol* econst)
{
    Constant con;
    markcdf4("Enum type");
    consttype = NC_ENUM;
    con.nctype = NC_ECONST;
    con.lineno = lineno;
    con.filled = 0;
    /* fix up econst to be a ref to an econst*/
    econst->objectclass = NC_TYPE;
    econst->subclass = NC_ECONST;
    {
	Symbol* defsym;
	defsym = locate(econst);
	if(defsym == NULL)
	    derror("Undefined or forward referenced enum constant: %s",econst->name);
	econst = defsym;
    }
    con.value.enumv = econst;
    return con;
}

static void
addtogroup(Symbol* sym)
{
    Symbol* grp = currentgroup();
    sym->container = grp;
    listpush(grp->subnodes,(void*)sym);
    setpathcurrent(sym);
}

/* Check for duplicate name of given type within current group*/
static int
dupobjectcheck(nc_class objectclass, Symbol* pattern)
{
    int i;
    Symbol* grp;
    if(pattern == NULL) return 0;
    grp = pattern->container;
    if(grp == NULL || grp->subnodes == NULL) return 0;
    for(i=0;i<listlength(grp->subnodes);i++) {
	Symbol* sym = (Symbol*)listget(grp->subnodes,i);
	if(!sym->is_ref && sym->objectclass == objectclass
	   && strcmp(sym->name,pattern->name)==0) return 1;
    }
    return 0;
}

static void
setpathcurrent(Symbol* sym)
{
    sym->is_prefixed = 0;
    sym->prefix = prefixdup(groupstack);
}

/* Convert an nc_type code to the corresponding Symbol*/
Symbol*
basetypefor(nc_type nctype)
{
    return primsymbols[nctype];
}

char*
specialname(int flag)
{
    switch (flag) {
    case _FILLVALUE_FLAG: return "_FillValue";
    case _FORMAT_FLAG: return "_Format";
    case _STORAGE_FLAG: return "_Storage";
    case _CHUNKSIZES_FLAG: return "_ChunkSizes";
    case _FLETCHER32_FLAG: return "_Fletcher32";
    case _DEFLATE_FLAG: return "_DeflateLevel";
    case _SHUFFLE_FLAG: return "_Shuffle";
    case _ENDIAN_FLAG: return "_Endianness";
    case _NOFILL_FLAG: return "_NoFill";
    default: break;
    }
    return "<unknown>";
}

static int
truefalse(Constant* con, int tag)
{
    if(con->nctype == NC_STRING) {
	char* sdata = con->value.stringv.stringv;
	if(strncmp(sdata,"false",NC_MAX_NAME) == 0
           || strncmp(sdata,"0",NC_MAX_NAME) == 0)
	    return 0;
	else if(strncmp(sdata,"true",NC_MAX_NAME) == 0
           || strncmp(sdata,"1",NC_MAX_NAME) == 0)
	    return 1;
	else goto fail;
    } else if(con->value.int32v < 0 || con->value.int32v > 1)
	goto fail;
    return con->value.int32v;

fail:
    derror("%s: illegal value",specialname(tag));
    return 0;
}

/* Since this may be affected by the _Format attribute, which
   may come last, capture all the special info and sort it out
   in semantics.
*/
static Symbol*
makespecial(int tag, Symbol* vsym, Symbol* tsym, void* data, int isconst)
{
    Symbol* attr = NULL;
    Datalist* list;
    Constant* con;
    Specialdata* special = (Specialdata*)malloc(sizeof(Specialdata));
    Constant iconst;
    int tf = 0;
    char* sdata = NULL;
    int idata =  -1;
    
    specials_flag += (tag == _FILLVALUE_FLAG ? 0 : 1);

    if(isconst) {
	con = (Constant*)data;
	list = builddatalist(1);
        dlappend(list,(Constant*)data);
    } else {
        list = (Datalist*)data;
        con = (Constant*)list->data;
    }

    if(tag == _FORMAT && vsym != NULL) {
	derror("_Format: must be global attribute");
	vsym = NULL;
    }

    switch (tag) {
    case _FLETCHER32_FLAG:
    case _SHUFFLE_FLAG:
    case _NOFILL_FLAG:
	iconst.nctype = (con->nctype == NC_STRING?NC_STRING:NC_INT);
	convert1(con,&iconst);
	tf = truefalse(&iconst,tag);
	break;
    case _FORMAT_FLAG:
    case _STORAGE_FLAG:
    case _ENDIAN_FLAG:
	iconst.nctype = NC_STRING;
	convert1(con,&iconst);
	if(iconst.nctype == NC_STRING)
	    sdata = iconst.value.stringv.stringv;
	else
	    derror("%s: illegal value",specialname(tag));
	break;
    case _DEFLATE_FLAG:
	iconst.nctype = NC_INT;
	convert1(con,&iconst);
	if(iconst.nctype == NC_INT)
	    idata = iconst.value.int32v;
	else
	    derror("%s: illegal value",specialname(tag));
	break;
    case _CHUNKSIZES_FLAG:
    case _FILLVALUE_FLAG:
	/* Handle below */
	break;
    default: PANIC1("unexpected special tag: %d",tag);
    }
    
    if(vsym != NULL) special = &vsym->var.special;
    if(tag == _FORMAT_FLAG) {
	struct Kvalues* kvalue;
	int found;
	found = 0;
	/* Use the table in main.c */
        for(kvalue=legalkinds;kvalue->name;kvalue++) {
	    if(strcmp(sdata,kvalue->name) == 0) {
		format_flag = kvalue->k_flag;
		found = 1;
	        break;
	    }
	}
	if(!found)
	    derror("_Format: illegal value: %s",sdata);
    } else if(tag == _FILLVALUE_FLAG) {
	special->_Fillvalue = list;
	/* fillvalue must be a single value*/
	if(list->length != 1)
	    derror("_FillValue: must be a single (possibly compound) value",
			vsym->name);
        /* check that the attribute value contains no fill values*/
        if(containsfills(list)) {
	    derror("Attribute data may not contain fill values (i.e. _ )");
        }
	/* _FillValue is also a real attribute*/
	if(vsym->objectclass != NC_VAR) {
	    derror("_FillValue attribute not associated with variable: %s",vsym->name);
	}
	if(tsym  == NULL) tsym = vsym->typ.basetype;
	else if(vsym->typ.basetype != tsym) {
	    derror("_FillValue attribute type does not match variable type: %s",vsym->name);
	}
	attr=makeattribute(install("_FillValue"),vsym,tsym,list,ATTRVAR);
    } else switch (tag) {
        case _STORAGE_FLAG:
            if(strcmp(sdata,"contiguous") == 0)
                special->_Storage = NC_CONTIGUOUS;
            else if(strcmp(sdata,"chunked") == 0)
                special->_Storage = NC_CHUNKED;
            else
                derror("_Storage: illegal value: %s",sdata);
            special->flags |= _STORAGE_FLAG;
            break;
        case _FLETCHER32_FLAG:
            special->_Fletcher32 = tf;
            special->flags |= _FLETCHER32_FLAG;
            break;
        case _DEFLATE_FLAG:
            special->_DeflateLevel = idata;
            special->flags |= _DEFLATE_FLAG;
            break;
        case _SHUFFLE_FLAG:
	    special->_Shuffle = tf;
            special->flags |= _SHUFFLE_FLAG;
            break;
        case _ENDIAN_FLAG:
            if(strcmp(sdata,"little") == 0)
                special->_Endianness = 1;
            else if(strcmp(sdata,"big") == 0)
                special->_Endianness = 2;
            else
                derror("_Endianness: illegal value: %s",sdata);
            special->flags |= _ENDIAN_FLAG;
            break;
        case _NOFILL_FLAG:
            special->_Fill = (1 - tf); /* negate */
            special->flags |= _NOFILL_FLAG;
            break;
        case _CHUNKSIZES_FLAG: {
	    int i;
            special->nchunks = list->length;
            special->_ChunkSizes = (size_t*)emalloc(sizeof(size_t)*special->nchunks);
            for(i=0;i<special->nchunks;i++) {
		iconst.nctype = NC_INT;
		convert1(&list->data[i],&iconst);
	        if(iconst.nctype == NC_INT) {
		    special->_ChunkSizes[i] = (size_t)iconst.value.int32v;
	        } else
		    derror("%s: illegal value",specialname(tag));
            }
            special->flags |= _CHUNKSIZES_FLAG;
	    /* Chunksizes => storage == chunked */
            special->flags |= _STORAGE_FLAG;
            special->_Storage = NC_CHUNKED;
            } break;
        default: PANIC1("makespecial: illegal token: %d",tag);
     }

    return attr;
}

static Symbol*
makeattribute(Symbol* asym,
		Symbol* vsym,
		Symbol* tsym,
		Datalist* data,
		Attrkind kind) /* global var or unknown*/
{
    asym->objectclass = NC_ATT;
    asym->data = data;
    addtogroup(asym);
    switch (kind) {
    case ATTRVAR:
        asym->att.var = vsym;
        asym->typ.basetype = tsym;
        listpush(attdefs,(void*)asym);
	break;
    case ATTRGLOBAL:
        asym->att.var = NULL; /* NULL => NC_GLOBAL*/
        asym->typ.basetype = tsym;
        listpush(gattdefs,(void*)asym);
	break;
    default: PANIC1("unexpected attribute type: %d",kind);
    }
    /* finally; check that the attribute value contains no fill values*/
    if(containsfills(data)) {
	derror("Attribute data may not contain fill values (i.e. _ ): %s",asym->name);
    }
    return asym;
}

static int
containsfills(Datalist* list)
{
    if(list != NULL) {
        int i;
        Constant* con = list->data;
        for(i=0;i<list->length;i++,con++) {
	    if(con->nctype == NC_COMPOUND) {
	        if(containsfills(con->value.compoundv)) return 1;	
	    } else if(con->nctype == NC_FILLVALUE) return 1;	
	}
    }
    return 0;
}

static void
datalistextend(Datalist* dl, Constant* con)
{
    dlappend(dl,con);
}

static void
vercheck(int ncid)
{
    char* tmsg = NULL;
    switch (ncid) {
    case NC_UBYTE: tmsg = "netCDF4 type: UBYTE"; break;
    case NC_USHORT: tmsg = "netCDF4 type: USHORT"; break;
    case NC_UINT: tmsg = "netCDF4 type: UINT"; break;
    case NC_INT64: tmsg = "netCDF4 type: INT64"; break;
    case NC_UINT64: tmsg = "netCDF4 type: UINT64"; break;
    case NC_STRING: tmsg = "netCDF4 type: STRING"; break;
    case NC_VLEN: tmsg = "netCDF4 type: VLEN"; break;
    case NC_OPAQUE: tmsg = "netCDF4 type: OPAQUE"; break;
    case NC_ENUM: tmsg = "netCDF4 type: ENUM"; break;
    case NC_COMPOUND: tmsg = "netCDF4 type: COMPOUND"; break;
    default: break;
    }
    if(tmsg != NULL) markcdf4(tmsg);
}

/*
Since the arguments are all simple constants,
we can evaluate the function immediately
and return its value.
Note that currently, only a single value can
be returned.
*/

static Constant
evaluate(Symbol* fcn, Datalist* arglist)
{
    Constant result;

    /* prepare the result */
    result.lineno = fcn->lineno;
    result.filled = 0;

    if(strcasecmp(fcn->name,"time") == 0) {
        char* timekind = NULL;
        char* timevalue = NULL;
        result.nctype = NC_DOUBLE;
        result.value.doublev = 0;
	/* int time([string],string) */
	switch (arglist->length) {
	case 2:
	    if(arglist->data[1].nctype != NC_STRING) {
	        derror("Expected function signature: time([string,]string)");
	        goto done;
	    }
	    /* fall thru */
	case 1:
	    if(arglist->data[0].nctype != NC_STRING) {
	        derror("Expected function signature: time([string,]string)");
	        goto done;
	    }
	    break;
	case 0:
	default: 
	    derror("Expected function signature: time([string,]string)");
	    goto done;
	}
	if(arglist->length == 2) {
	    timekind = arglist->data[0].value.stringv.stringv;
            timevalue = arglist->data[1].value.stringv.stringv;
	} else
            timevalue = arglist->data[0].value.stringv.stringv;
	if(timekind == NULL) { /* use cd time as the default */
            cdCompTime comptime;
	    CdTime cdtime;
	    cdCalenType timetype = cdStandard;
	    cdChar2Comp(timetype,timevalue,&comptime);
	    /* convert comptime to cdTime */
	    cdtime.year = comptime.year;	    
	    cdtime.month = comptime.month;
	    cdtime.day = comptime.day;    
	    cdtime.hour = comptime.hour;
	    cdtime.baseYear = 1970;
	    cdtime.timeType = CdChron;
	    /* convert to double value */
	    Cdh2e(&cdtime,&result.value.doublev);
        } else {
	    derror("Time conversion '%s' not supported",timekind);
	    goto done;
	}
    } else {	/* Unknown function */
	derror("Unknown function name: %s",fcn->name);
	goto done;
    }

done:
    return result;
}

