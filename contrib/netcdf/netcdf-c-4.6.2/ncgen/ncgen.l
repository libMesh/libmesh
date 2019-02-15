%{
/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Id: ncgen.l,v 1.24 2009/09/25 18:22:35 dmh Exp $
 *********************************************************************/

/* Problems:
1. We assume the input is true utf8.
   Unfortunately, we may actually get iso-latin-8859-1.
   This means that there will be ambiguity about the characters
   in the range 128-255 because they will look like n-byte unicode
   when they are 1-byte 8859 characters. Because of our encoding,
   8859 characters above 128 will be handles as n-byte utf8 and so
   will probably not lex correctly.
   Solution: assume utf8 and note in the documentation that
   ISO8859 is specifically unsupported.
2. The netcdf function NC_check_name in string.c must be modified to
   conform to the use of UTF8.
3. We actually have three tests for UTF8 of increasing correctness
   (in the sense that the least correct will allow some sequences that
    are technically illegal UTF8).
   The tests are derived from the table at
	    http://www.w3.org/2005/03/23-lex-U
   We include lexical definitions for all three, but use the second version.
4. Single character constants enclosed in '...' cannot be
   utf-8, so we assume they are by default encoded using the 1-byte
   subset of utf-8. It turns out that this subset is in fact
   equivalent to US-Ascii (7-bit).
   We could use ISO-8859-1, but that conflicts with UTF-8 above value 127.
*/

/* lex specification for tokens for ncgen */

/* Fill value used by ncdump from version 2.4 and later.  Should match
   definition of FILL_STRING in ../ncdump/vardata.h */

#include "ncgen.h"
#include "ncgeny.h"

#define FILL_STRING "_"
#define XDR_INT32_MIN (-2147483647-1)
#define XDR_INT32_MAX 2147483647
#define XDR_INT64_MIN  (-9223372036854775807LL-1)
#define XDR_INT64_MAX  (9223372036854775807LL)

#undef DEBUG
#ifdef DEBUG
static int MIN_BYTE = NC_MIN_BYTE;
static int MIN_SHORT = NC_MIN_SHORT;
static int MIN_INT = NC_MIN_INT;
static int MAX_BYTE = NC_MAX_BYTE;
static int MAX_SHORT = NC_MAX_SHORT;
static int MAX_INT = NC_MAX_INT;
static int MAX_UBYTE = NC_MAX_UBYTE;
static int MAX_USHORT = NC_MAX_USHORT;
static unsigned int MAX_UINT = NC_MAX_UINT;

#undef NC_MIN_BYTE
#undef NC_MIN_SHORT
#undef NC_MIN_INT
#undef NC_MAX_BYTE
#undef NC_MAX_SHORT
#undef NC_MAX_INT
#undef NC_MAX_UBYTE
#undef NC_MAX_USHORT
#undef NC_MAX_UINT

#define NC_MIN_BYTE MIN_BYTE
#define NC_MIN_SHORT MIN_SHORT
#define NC_MIN_INT MIN_INT
#define NC_MAX_BYTE MAX_BYTE
#define NC_MAX_SHORT MAX_SHORT
#define NC_MAX_INT MAX_INT
#define NC_MAX_UBYTE MAX_UBYTE
#define NC_MAX_USHORT MAX_USHORT
#define NC_MAX_UINT MAX_UINT
#endif

#define TAGCHARS "BbSsLlUu"

#define tstdecimal(ch) ((ch) >= '0' && (ch) <= '9')

/*Mnemonics*/
#define ISIDENT 1

/* Define a fake constant indicating that
   no tag was specified */
#define NC_NOTAG (-1)

char errstr[100];		/* for short error messages */

int lineno;                    /* line number for error messages */
Bytebuffer* lextext;           /* name or string with escapes removed */

#define YY_BREAK                /* defining as nothing eliminates unreachable
				   statement warnings from flex output,
                                   but make sure every action ends with
                                   "return" or "break"! */

int specialconstants;		/* 1 if nan, nanf, infinity, etc is used */
double double_val;             /* last double value read */
float float_val;               /* last float value read */
long long int64_val;           /* last int64 value read */
int int32_val;                 /* last int32 value read */
short int16_val;               /* last short value read */
unsigned long long uint64_val; /* last int64 value read */
unsigned int uint32_val;       /* last int32 value read */
unsigned short uint16_val;     /* last short value read */
char char_val;                 /* last char value read */
signed char byte_val;                 /* last byte value read */
unsigned char ubyte_val;       /* last byte value read */

static Symbol* makepath(char* text);
static int lexdebug(int);
static unsigned long long parseULL(char* text, int*);
static nc_type downconvert(unsigned long long uint64, int*, int, int);
static int tagmatch(nc_type nct, int tag);
static int nct2lexeme(nc_type nct);
static int collecttag(char* text, char** stagp);

struct Specialtoken specials[] = {
{"_FillValue",_FILLVALUE,_FILLVALUE_FLAG},
{"_Format",_FORMAT,_FORMAT_FLAG},
{"_Storage",_STORAGE,_STORAGE_FLAG},
{"_ChunkSizes",_CHUNKSIZES,_CHUNKSIZES_FLAG},
{"_Fletcher32",_FLETCHER32,_FLETCHER32_FLAG},
{"_DeflateLevel",_DEFLATELEVEL,_DEFLATE_FLAG},
{"_Shuffle",_SHUFFLE,_SHUFFLE_FLAG},
{"_Endianness",_ENDIANNESS,_ENDIAN_FLAG},
{"_NoFill",_NOFILL,_NOFILL_FLAG},
{"_NCProperties",_NCPROPS,_NCPROPS_FLAG},
{"_IsNetcdf4",_ISNETCDF4,_ISNETCDF4_FLAG},
{"_SuperblockVersion",_SUPERBLOCK,_SUPERBLOCK_FLAG},
{"_Filter",_FILTER,_FILTER_FLAG},
{NULL,0} /* null terminate */
};

%}
%x ST_C_COMMENT
%x TEXT
%p 6000

/* The most correct (validating) version of UTF8 character set
   (Taken from: http://www.w3.org/2005/03/23-lex-U)

The lines of the expression cover the UTF8 characters as follows:
1. non-overlong 2-byte
2. excluding overlongs
3. straight 3-byte
4. excluding surrogates
5. straight 3-byte
6. planes 1-3
7. planes 4-15
8. plane 16

UTF8   ([\xC2-\xDF][\x80-\xBF])                       \
     | (\xE0[\xA0-\xBF][\x80-\xBF])                   \
     | ([\xE1-\xEC][\x80-\xBF][\x80-\xBF])            \
     | (\xED[\x80-\x9F][\x80-\xBF])                   \
     | ([\xEE-\xEF][\x80-\xBF][\x80-\xBF])            \
     | (\xF0[\x90-\xBF][\x80-\xBF][\x80-\xBF])        \
     | ([\xF1-\xF3][\x80-\xBF][\x80-\xBF][\x80-\xBF]) \
     | (\xF4[\x80-\x8F][\x80-\xBF][\x80-\xBF])        \

*/

/* Wish there was some way to ifdef lex files */

/*The most relaxed version of UTF8 (not used)
UTF8 ([\xC0-\xD6].)|([\xE0-\xEF]..)|([\xF0-\xF7]...)
*/

/*The partially relaxed version of UTF8, and the one used here */
UTF8 ([\xC0-\xD6][\x80-\xBF])|([\xE0-\xEF][\x80-\xBF][\x80-\xBF])|([\xF0-\xF7][\x80-\xBF][\x80-\xBF][\x80-\xBF])

/* The old definition of ID
ID ([A-Za-z_]|{UTF8})([A-Z.@#\[\]a-z_0-9+-]|{UTF8})*
*/

/* Don't permit control characters or '/' in names, but other special
   chars OK if escaped.  Note that to preserve backwards
   compatibility, none of the characters _.@+- should be escaped, as
   they were previously permitted in names without escaping.  */

idescaped \\[ !"#$%&'()*,:;<=>?\[\\\]^`{|}~]
numescaped \\[0-9]

/* New definition to conform to a subset of string.c */
ID ([a-zA-Z_]|{UTF8}|{numescaped})([a-zA-Z0-9_.@+-]|{UTF8}|{idescaped})*

escaped		\\.

/* Capture a datasetidentifier */
/* DATASETID ([a-zA-Z0-9!#$%&*:;<=>?/^|~_.@+-]|{UTF8})* */
DATASETID [^{][^{]*


/* Note: this definition of string will work for utf8 as well,
   although it is a very relaxed definition
*/
nonquotes	([^"\\]|{escaped})*
exp		([eE][+-]?[0-9]+)

OPAQUESTRING	(0[xX][0-9A-Fa-f][0-9A-Fa-f]*)

PATH    ([/]|([/]{ID})([/]{ID})*)

XUNUMBER     {OPAQUESTRING}([Ss]|[Ll]|[Ll][Ll])?
NUMBER       [+-]?[0-9][0-9]*[Uu]?([BbSs]|[Ll]|[Ll][Ll])?
DBLNUMBER    [+-]?[0-9]*\.[0-9]*{exp}?[LlDd]?|[+-]?[0-9]*{exp}[LlDd]?
FLTNUMBER    [+-]?[0-9]*\.[0-9]*{exp}?[Ff]|[+-]?[0-9]*{exp}[Ff]

SPECIAL "_FillValue"|"_Format"|"_Storage"|"_ChunkSizes"|"_Fletcher32"|"_DeflateLevel"|"_Shuffle"|"_Endianness"|"_NoFill"|"_NCProperties"|"_IsNetcdf4"|"_SuperblockVersion"|"_Filter"

USASCII   [\x01-\x7F]

%%
[ \r\t\f]+	{ /* whitespace */
		  break;
		}

\/\/.*		        { /* comment */
                          break;
                        }

\"{nonquotes}\"		{int len; char* s = NULL;
			 /* In netcdf4, this will be used in a variety
                            of places, so only remove escapes */
/*
if(yyleng > MAXTRST) {
yyerror("string too long, truncated\n");
yytext[MAXTRST-1] = '\0';
}
*/
			len = unescape((char *)yytext+1,yyleng-2,!ISIDENT,&s);
			if(len < 0) {
			    sprintf(errstr,"Illegal character: %s",yytext);
			    yyerror(errstr);
			}
			bbClear(lextext);
			bbAppendn(lextext,s,len);
			bbNull(lextext);
			if(s) efree(s);
		 	return lexdebug(TERMSTRING);
		        }

{OPAQUESTRING}	{ /* drop leading 0x; pad to even number of chars */
		char* p = yytext+2;
		int len = yyleng - 2;
		bbClear(lextext);
	        bbAppendn(lextext,p,len);
	        if((len % 2) == 1) bbAppend(lextext,'0');
		bbNull(lextext);
		/* convert all chars to lower case */
		for(p=bbContents(lextext);(int)*p;p++) *p = tolower(*p);
		return lexdebug(OPAQUESTRING);
		}

compound|struct|structure {return lexdebug(COMPOUND);}
enum				{return lexdebug(ENUM);}
opaque				{return lexdebug(OPAQUE_);}

float|real		{return lexdebug(FLOAT_K);}
char			{return lexdebug(CHAR_K);}
byte			{return lexdebug(BYTE_K);}
ubyte			{return lexdebug(UBYTE_K);}
short			{return lexdebug(SHORT_K);}
ushort			{return lexdebug(USHORT_K);}
long|int|integer	{return lexdebug(INT_K);}
ulong|uint|uinteger	{return lexdebug(UINT_K);}
int64			{return lexdebug(INT64_K);}
uint64			{return lexdebug(UINT64_K);}
double			{return lexdebug(DOUBLE_K);}
string                  {return lexdebug(STRING_K);}

unlimited|UNLIMITED	{int32_val = -1;
			 return lexdebug(NC_UNLIMITED_K);}

types:			{return lexdebug(TYPES);}
dimensions:		{return lexdebug(DIMENSIONS);}
variables:		{return lexdebug(VARIABLES);}
data:			{return lexdebug(DATA);}
group:			{return lexdebug(GROUP);}

(netcdf|NETCDF|netCDF)  {BEGIN(TEXT);return lexdebug(NETCDF);}

DoubleInf|-?Infinity { /* missing value (pre-2.4 backward compatibility) */
                if (yytext[0] == '-') {
		    double_val = NEGNC_INFINITE;
                } else {
		    double_val = NC_INFINITE;
                }
		specialconstants = 1;
		return lexdebug(DOUBLE_CONST);
		}
NaN|nan { /* missing value (pre-2.4 backward compatibility) */
		double_val = NAN;
		specialconstants = 1;
		return lexdebug(DOUBLE_CONST);
		}

FloatInf|-?Infinityf|-?Inff {/* missing value (pre-2.4 backward compatibility)*/
                if (yytext[0] == '-') {
		    float_val = NEGNC_INFINITEF;
                } else {
		    float_val = NC_INFINITEF;
                }
		specialconstants = 1;
		return lexdebug(FLOAT_CONST);
		}
NaNf|nanf { /* missing value (pre-2.4 backward compatibility) */
		float_val = NANF;
		specialconstants = 1;
		return lexdebug(FLOAT_CONST);
		}

NIL|nil|Nil {
#ifdef USE_NETCDF4
		if(l_flag == L_C || l_flag == L_BINARY)
		    return lexdebug(NIL);
		yyerror("NIL only allowed for netcdf-4 and for -lc or -lb");
#else
		yyerror("NIL only allowed for netcdf-4 and for -lc or -lb");
#endif
		}

{PATH}		{
		bbClear(lextext);
		bbAppendn(lextext,(char*)yytext,yyleng+1); /* include null */
	        bbNull(lextext);
		yylval.sym = makepath(bbContents(lextext));
		return lexdebug(PATH);
		}


{SPECIAL}      {struct Specialtoken* st;
		bbClear(lextext);
		bbAppendn(lextext,(char*)yytext,yyleng+1); /* include null */
		bbNull(lextext);
		for(st=specials;st->name;st++) {
		    if(strcmp(bbContents(lextext),st->name)==0) {return lexdebug(st->token);}
		}
		return 0;
		}

<TEXT>{DATASETID} {
		    int c;
		    char* p; char* q;
		    /* copy the trimmed name */
		    bbClear(lextext);
		    bbAppendn(lextext,(char*)yytext,yyleng+1); /* include null */
		    bbNull(lextext);
		    p = bbContents(lextext);
		    q = p;
		    while((c=*p++)) {if(c > ' ') *q++ = c;}
		    *q = '\0';
		    if(datasetname == NULL)
		        datasetname = bbDup(lextext);
		    BEGIN(INITIAL);
		    return lexdebug(DATASETID);
		}

{ID}		{ char* id = NULL; int len;
		    len = strlen(yytext);
		    len = unescape(yytext,len,ISIDENT,&id);
		    if(NCSTREQ(id, FILL_STRING)) {
			efree(id);
		        return lexdebug(FILLMARKER);
		    }
		    yylval.sym = install(id);
		    efree(id);
		    return lexdebug(IDENT);
		}

{NUMBER}    {
		/*
		  We need to try to see what size of integer ((u)int).
		  Technically, the user should specify, but...
		  If out of any integer range, then complain
		*/
		    int slen = strlen(ncgtext);
		    char* stag = NULL;
		    int tag = NC_NAT;
		    int isneg = 0;
		    int c = ncgtext[0];
		    int fail = 0;
		    nc_type nct = 0;
		    char* pos = NULL;
		    int hasU = 0;

		    /* capture the tag string */
		    tag = collecttag(ncgtext,&stag);
		    if(tag == NC_NAT) {
			sprintf(errstr,"Illegal integer suffix: %s",stag);
			yyerror(errstr);
			goto done;
		    }
		    /* drop the tag from the input text */
		    ncgtext[slen - strlen(stag)] = '\0';
		    hasU = isuinttype(tag);
		    if(!tstdecimal(c)) {
			pos = ncgtext+1;
			isneg = (c == '-');
		    } else
		        pos = ncgtext;
		    if(isneg && hasU) {
			sprintf(errstr,"Unsigned integer cannot be signed: %s",ncgtext);
			yyerror(errstr);
			goto done;
		    }
		    uint64_val = parseULL(pos,&fail);
		    if(fail) {
			sprintf(errstr,"integer constant out of range: %s",ncgtext);
			yyerror(errstr);
			goto done;
		    }
		    /* Down convert to smallest possible range */
		    nct = downconvert(uint64_val,&tag,isneg,hasU);
		    switch (k_flag) {
		    case NC_FORMAT_64BIT_DATA:
		    case NC_FORMAT_NETCDF4:
			    return lexdebug(nct2lexeme(nct));
		    case NC_FORMAT_CLASSIC:
		    case NC_FORMAT_64BIT_OFFSET:
		    case NC_FORMAT_NETCDF4_CLASSIC:
			    if(nct > NC_INT) {
				sprintf(errstr,"Illegal integer constant for classic format: %s",ncgtext);
				yyerror(errstr);
				goto done;
			    }
		    }

		    if(!tagmatch(nct,tag))  {
			semwarn(lineno,"Warning: Integer out of range for tag: %s; tag treated as changed.",ncgtext);
		    }
		    return lexdebug(nct2lexeme(nct));
done: return 0;
	    }

{XUNUMBER}    {
		int c;
		int token = 0;
		int slen = strlen(yytext);
		char* stag = NULL;
	        int tag = NC_NAT;
		char* hex = yytext+2; /* point to first true hex digit */
		int xlen = (slen - 3);  /* true hex length */

		yytext[slen-1] = '\0';
	        /* capture the tag string */
		tag = collecttag(yytext,&stag);
		if(tag == NC_NAT) {
		    sprintf(errstr,"Illegal integer suffix: %s",stag);
		    yyerror(errstr);
		    goto done;
		}
		yytext[slen - strlen(stag)] = '\0';
	        if(xlen > 16) { /* truncate hi order digits */
		    hex += (xlen - 16);
		}
		/* convert to an unsigned long long */
		uint64_val = 0;
		while((c=*hex++)) {
		    unsigned int hexdigit = (c <= '9'?(c-'0'):(c-'a')+0xa);
		    uint64_val = ((uint64_val << 4) | hexdigit);
		}
		switch (tag) {
		case NC_USHORT:
		    uint16_val = (unsigned short)uint64_val;
		    token = USHORT_CONST;
		    break;
		case NC_UINT:
		    token = UINT_CONST;
		    break;
		case NC_UINT64:
		    token = UINT64_CONST;
		    break;
		default: /* should never happen */
		    if (sscanf((char*)yytext, "%i", &uint32_val) != 1) {
		        sprintf(errstr,"bad unsigned int constant: %s",(char*)yytext);
		        yyerror(errstr);
		    }
		    token = UINT_CONST;
		}
		return lexdebug(token);
	    }
{DBLNUMBER}	{
		if (sscanf((char*)yytext, "%le", &double_val) != 1) {
		    sprintf(errstr,"bad long or double constant: %s",(char*)yytext);
		    yyerror(errstr);
		}
                return lexdebug(DOUBLE_CONST);
                }
{FLTNUMBER}  {
		if (sscanf((char*)yytext, "%e", &float_val) != 1) {
		    sprintf(errstr,"bad float constant: %s",(char*)yytext);
		    yyerror(errstr);
		}
                return lexdebug(FLOAT_CONST);
                }
\'[^\\]\'          {
	        (void) sscanf((char*)&yytext[1],"%c",&byte_val);
		return lexdebug(BYTE_CONST);
                }
\'\\[0-7][0-7][0-7]\'  {
		int oct = unescapeoct(&yytext[2]);
		if(oct < 0) {
		    sprintf(errstr,"bad octal character constant: %s",(char*)yytext);
		    yyerror(errstr);
		}
	        byte_val = (unsigned int)oct;
		return lexdebug(BYTE_CONST);
                }
\'\\[xX][0-9a-fA-F][0-9a-fA-F]\'  {
		int hex = unescapehex(&yytext[3]);
		if(byte_val < 0) {
		    sprintf(errstr,"bad hex character constant: %s",(char*)yytext);
		    yyerror(errstr);
		}
		byte_val = (unsigned int)hex;
		return lexdebug(BYTE_CONST);
                }
\'\\.\'        {
	       switch ((char)yytext[2]) {
	          case 'a': byte_val = '\007'; break; /* not everyone under-
						       * stands '\a' yet */
     	          case 'b': byte_val = '\b'; break;
		  case 'f': byte_val = '\f'; break;
		  case 'n': byte_val = '\n'; break;
		  case 'r': byte_val = '\r'; break;
		  case 't': byte_val = '\t'; break;
		  case 'v': byte_val = '\v'; break;
		  case '\\': byte_val = '\\'; break;
		  case '?': byte_val = '\177'; break;
		  case '\'': byte_val = '\''; break;
		  default: byte_val = (char)yytext[2];
	           }
		return lexdebug(BYTE_CONST);
                }

\n		{
		lineno++ ;
                break;
		}

"/""*"	{/*initial*/
	    BEGIN(ST_C_COMMENT);
	    break;
	}

<ST_C_COMMENT>([^*]|"*"[^/])*	{/* continuation */
				     break;
				}

<ST_C_COMMENT>"*/"	{/* final */
			    BEGIN(INITIAL);
			    break;
			}

<ST_C_COMMENT><<EOF>>	{/* final, error */
			    fprintf(stderr,"unterminated /**/ comment");
			    BEGIN(INITIAL);
			    break;
			}

.		{/* Note: this next rule will not work for UTF8 characters */
		return lexdebug(yytext[0]) ;
		}
%%
static int
lexdebug(int token)
{
    if(debug >= 2)
    {
	char* text = yytext;
	text[yyleng] = 0;
        fprintf(stderr,"Token=%d |%s| line=%d\n",token,text,lineno);
    }
    return token;
}

int
lex_init(void)
{
    lineno = 1;
    lextext = bbNew();
    if(0) unput(0); /* keep -Wall quiet */
    return 0;
}

static Symbol*
makepath(char* text0)
{
    /* Create a reference symbol.
       Convert path to a sequence of symbols.
       Use last name as symbol name (with root group reference ('/') as exception).
    */
    Symbol* refsym = NULL;
    /* walk the path converting to a sequence of symbols */
    if(strcmp(text0,"/")==0) {
	/* special case of root reference */
	refsym = rootgroup;
    } else {
        List* prefix = listnew();
	/* split the text into IDENT chunks, convert to symbols */
        Symbol* container = rootgroup;
	char *ident, *p, *match;
        char* text = estrdup(text0);
	int lastident;

	ident=text+1; p=ident; /* skip leading '/' */
	lastident = 0;
	do {
	    match = esc_strchr(p,'/',0);
	    lastident = (*match == '\0');
	    *match='\0';
	    (void)unescape(p,strlen(p),ISIDENT,&ident);
	    refsym = lookupingroup(NC_GRP,ident,container);
	    if(!lastident) {
	        if(refsym == NULL) {
	            sprintf(errstr,"Undefined or forward referenced group: %s",ident);
	            yyerror(errstr);
		    refsym = rootgroup;
	        } else
	             listpush(prefix,(void*)refsym);
	    } else {/* lasiident is true */
		refsym = install(ident); /* create as symbol */
	        refsym->objectclass = NC_GRP;/* tentative */
		refsym->ref.is_ref = 1;
		refsym->container = container;
		refsym->subnodes = listnew();
	    }
   	    container = refsym;
	    p = (lastident?match:match+1);
	    if(ident) efree(ident);
	} while(!lastident);
        refsym->prefix = prefix;
	efree(text);
    }
    return refsym;
}

/*
Parse a simple string of digits into an unsigned long long
Return the value.
*/
static unsigned long long
parseULL(char* text, int* failp)
{
    extern int errno;
    char* endptr;
    unsigned long long uint64 = 0;

    errno = 0; endptr = NULL;
    assert(tstdecimal(text[0]));
#ifdef HAVE_STRTOULL
    uint64 = strtoull(text,&endptr,10);
    if(errno == ERANGE) {
	if(failp) *failp = ERANGE;
	return 0;
    }
#else /*!(defined HAVE_STRTOLL && defined HAVE_STRTOULL)*/
    sscanf((char*)text, "%llu", &uint64);
    /* Have no useful way to detect out of range */
#endif /*!(defined HAVE_STRTOLL && defined HAVE_STRTOULL)*/
    return uint64;
}


/**
Given the raw bits, the sign char, the tag, and hasU
fill in the appropriate *_val field
and return the type.
Note that we cannot return unsigned types if running pure netcdf classic.
The rule is to pick the smallest enclosing type.

The rule used here is that the tag (the suffix, if any)
always takes precedence and the value is modified to conform
if possible, otherwise out-of-range is signalled.
For historical reasons (ncgen3), values that fit as unsigned
are acceptable for the signed tag and conversion is attempted;
e.g. 65535s; is legal and is return as a negative short.
*/
static nc_type
downconvert(unsigned long long uint64, int* tagp, int isneg, int hasU)
{
    nc_type nct = NC_NAT;
    int tag = *tagp;
    int bit63set = (uint64 >> 63);
    long long int64 = *((long long*)&uint64);

    if(isneg && hasU) {
	return (*tagp = NC_NAT);
    }
    /* To simplify the code, we look for special case of NC_UINT64
       constants that will not fit into an NC_INT64 constant.
     */
    if(tag == NC_UINT64 && bit63set) {
        uint64_val = uint64;
	return tag;
    }
    /* At this point we need deal only with int64 value */
    /* Apply the isneg */
    if(isneg)
	int64 = - int64;

    if(tag == NC_NOTAG) {
        /* If we have no other info, then assume NC_(U)INT(64) */
	if(int64 >= NC_MIN_INT && int64 <= NC_MAX_INT) {
	    nct = (tag = NC_INT);
	    int32_val = (signed int)int64;
	} else if(int64 >= 0 && int64 <= NC_MAX_UINT) {
	        nct = (tag = NC_UINT);
	        uint32_val = (unsigned int)int64;
	} else if(int64 < 0) {
		nct = (tag = NC_INT64);
	        int64_val = (signed long long)int64;
	} else {
	        nct = (tag = NC_UINT64);
	        uint64_val = (unsigned long long)int64;
	}
        goto done;
    }
    if(isuinttype(tag) && int64 < 0)
	goto outofrange;
    switch (tag) {
    case NC_UBYTE:
	    if(int64 <= NC_MAX_UBYTE) {
	        nct = NC_UBYTE;
	        ubyte_val = (unsigned char)int64;
	    } else
		goto outofrange;
	    break;
    case NC_USHORT:
	    if(int64 <= NC_MAX_USHORT) {
	        nct = NC_USHORT;
	        uint16_val = (unsigned short)int64;
	    } else
	       goto outofrange;
	    break;
    case NC_UINT:
	    if(int64 <= NC_MAX_UINT) {
	        nct = NC_UINT;
	        uint32_val = (unsigned int)int64;
	    } else
		goto outofrange;
	    break;
    case NC_UINT64:
            if(int64 <= NC_MAX_UINT64) {
	        nct = NC_UINT64;
	        uint64_val = uint64;
	    } else
		goto outofrange;
	    break;
    case NC_INT64:
	    nct = NC_INT64;
	    int64_val = int64;
	    break;
    case NC_BYTE:
	    nct = NC_BYTE;
	    byte_val = (signed char)int64;
	    break;
    case NC_SHORT:
	    nct = NC_SHORT;
	    int16_val = (signed short)int64;
	    break;
    case NC_INT:
	    nct = NC_INT;
	    int32_val = (signed int)int64;
	    break;
    default:
	    goto outofrange;
    }

done:
    *tagp = tag;
    return nct;
outofrange:
    yyerror("Value out of range");
    return NC_NAT;
}

static int
nct2lexeme(nc_type nct)
{
    switch(nct) {
    case NC_BYTE: return BYTE_CONST;
    case NC_CHAR: return CHAR_CONST;
    case NC_SHORT: return SHORT_CONST;
    case NC_INT: return INT_CONST;
    case NC_UBYTE: return UBYTE_CONST;
    case NC_USHORT: return USHORT_CONST;
    case NC_UINT: return UINT_CONST;
    case NC_INT64: return INT64_CONST;
    case NC_UINT64: return UINT64_CONST;
    }
    return 0;
}

static int
tagmatch(nc_type nct, int tag)
{
    if(tag == NC_NAT || tag ==  NC_NOTAG)
	return 1;
    return nct == tag;
}

/* capture the tag string */
static int
collecttag(char* text, char** stagp)
{
    char* stag0;
#define MAXTAGLEN 3
    char stag[MAXTAGLEN+1];
    int slen = strlen(text);
    int staglen;
    int tag = NC_NAT;
    int hasU = 0;

    for(stag0 = text+(slen-1);stag0 > 0;stag0--) {
	if(strchr(TAGCHARS,*stag0) == NULL) {stag0++; break;}
    }
    if(stagp) *stagp = stag0;
    staglen = strlen(stag0);
    if(staglen == 0)
	return NC_NOTAG;
    if(staglen > MAXTAGLEN)
	return tag;
    strncpy(stag,stag0,sizeof(stag));
    stag[MAXTAGLEN] = '\0';
    if(stag[0] == 'U' || stag[0] == 'u') {
	hasU = 1;
        memmove(stag,stag+1,MAXTAGLEN);
	staglen--;
    } else if(stag[staglen-1] == 'U' || stag[staglen-1] == 'u') {
	hasU = 1;
	staglen--;
	stag[staglen] = '\0';
    }
    if(strlen(stag) == 0 && hasU) {
	tag = NC_UINT;
    } else if(strlen(stag) == 1) {
	switch (stag[0]) {
	case 'B': case 'b': tag = (hasU ? NC_UBYTE : NC_BYTE); break;
	case 'S': case 's': tag = (hasU ? NC_USHORT : NC_SHORT); break;
	case 'L': case 'l': tag = (hasU ? NC_UINT : NC_INT); break;
	default: break;
	}
    } else if(strcasecmp(stag,"ll") == 0) {
	tag = (hasU ? NC_UINT64 : NC_INT64);
    }
    if(tag == NC_NAT) {
	if(strlen(stag) > 0)
	    return tag;
	tag = NC_NAT;
    }
    return tag;
}
