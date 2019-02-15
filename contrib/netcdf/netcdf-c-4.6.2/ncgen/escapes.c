/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/ncgen/escapes.c,v 1.5 2010/04/04 19:39:44 dmh Exp $
 *********************************************************************/

#include "includes.h"
#include "ncutf8.h"

#define HEXCHARS "0123456789abcdefABCDEF"
#define OCTCHARS "01234567"

/* Forward*/
static void initcodify(void);
static char* ccodify(const char*);
static char* f77codify(const char*);
static char* jcodify(const char*);

/*
Given a character c, fill s with the character suitably escaped.
E.g. c = '\t' => s="\t"
Caller must ensure enough space.
Watch out for embedded NULs.
Currently passes unicode thru unchanged.
Returns s as its result.
*/

char*
escapifychar(unsigned int c, char* s0, int quote)
{
    char* s = s0;
    if(c == '\\') {
	*s++ = '\\'; *s++='\\';
    } else if(c == quote) {
	*s++ = '\\'; *s++=(char)quote;
    } else if(c >= ' ' && c != '\177') {
	*s++ = (char)c;
    } else if((c & 0x80) != 0) {/* Unicode */
	*s++ = (char)c;
    } else {
        switch (c) {
	case '\b': strcpy(s,"\\b"); s+=2; break;
	case '\f': strcpy(s,"\\f"); s+=2; break;
	case '\n': strcpy(s,"\\n"); s+=2; break;
	case '\r': strcpy(s,"\\r"); s+=2; break;
	case '\t': strcpy(s,"\\t"); s+=2; break;
	case '\v': strcpy(s,"\\v"); s+=2; break;
	default: {
	    unsigned int oct1 = (c & 007);
	    unsigned int oct2 = ((c >> 3) & 007);
	    unsigned int oct3 = ((c >> 6) & 003);
	    *s++ = '\\';
	    *s++ = oct3 + '0';
	    *s++ = oct2 + '0';
	    *s++ = oct1 + '0';
	} break;
	}
    }
    *s = '\0';
    return s0;
}

/* Return a pool string that is s0 with all characters*/
/* ecaped that require it.  The resulting string is not*/
/* surrounded by quotes.*/
/* Since the string might actually contain nulls, specify the length.*/

char*
escapify(char* s0, int quote, size_t len)
{
    int i;
    char* result;
    result = poolalloc(1+4*len); /* overkill to support maximal expansion*/
    result[0] = '\0';
    for(i=0;i<len;i++) {
	char tmp[8];
	escapifychar((unsigned int)s0[i],tmp,quote);
        strcat(result,tmp);
    }
    return result;        
}

char*
escapifyname(char* s0)
{
    return escapify(s0,'"',strlen(s0));
}

void
cquotestring(Bytebuffer* databuf, char quote)
{
    char* escaped = escapify(bbContents(databuf),'"',bbLength(databuf));
    bbClear(databuf);
    bbAppend(databuf,quote);
    bbCat(databuf,escaped);
    bbAppend(databuf,quote);
}

/*
 * Replace special chars in name so it can be used in C and Fortran
 * variable names without causing syntax errors.  Here we just replace
 * each "-" in a name with "_MINUS_", each "." with "_PERIOD_", etc.
 * For bytes with high bit set, from UTF-8 encoding of Unicode, just
 * replace with "_xHH", where each H is the appropriate hex digit.
 * However, if the utf flag is set, then just pass utf characters as is.
 * If a name begins with a number N, such as "4LFTX", replace with
 * "DIGIT_N_", such as "DIGIT_4_LFTX".
 * Note that apparently, FORTRAN will not allow a leading underscore,
 * so remove if we are doing fortran.
 *
 * It is required that codify be idempotent:
 * i.e. codify(codify(s)) == codify(s)
 *
 * Returned name is pool alloc'd so is transient
 */

static int init = 0;
static char* repls[256];	/* replacement string for each char */
static int lens[256];	/* lengths of replacement strings */
static struct {
	char c;
	char *s;
} ctable[] = {
	{' ', "_SPACE_"},
	{'!', "_EXCLAMATION_"},
	{'"', "_QUOTATION_"},
	{'#', "_HASH_"},
	{'$', "_DOLLAR_"},
	{'%', "_PERCENT_"},
	{'&', "_AMPERSAND_"},
	{'\'', "_APOSTROPHE_"},
	{'(', "_LEFTPAREN_"},
	{')', "_RIGHTPAREN_"},
	{'*', "_ASTERISK_"},
	{'+', "_PLUS_"},
	{',', "_COMMA_"},
	{'-', "_MINUS_"},
	{'.', "_PERIOD_"},
	{':', "_COLON_"},
	{';', "_SEMICOLON_"},
	{'<', "_LESSTHAN_"},
	{'=', "_EQUALS_"},
	{'>', "_GREATERTHAN_"},
	{'?', "_QUESTION_"},
	{'@', "_ATSIGN_"},
	{'[', "_LEFTBRACKET_"},
	{'\\', "_BACKSLASH_"},
	{']', "_RIGHTBRACKET_"},
	{'^', "_CIRCUMFLEX_"},
	{'`', "_BACKQUOTE_"},
	{'{', "_LEFTCURLY_"},
	{'|', "_VERTICALBAR_"},
	{'}', "_RIGHTCURLY_"},
	{'~', "_TILDE_"},
 	{'/', "_SLASH_"},
};
static int idtlen;
static int hexlen;
static Bytebuffer* newname;

static void
initcodify(void)
{
    int nctable = (sizeof(ctable))/(sizeof(ctable[0]));
    int i;
    char *rp;

    newname = bbNew();
    idtlen = strlen("DIGIT_n_"); /* initial digit template */
    hexlen = strlen("_XHH"); /* template for hex of non-ASCII bytes */
    for(i = 0; i < 128; i++) {
        rp = ecalloc(2);
        rp[0] = i;
        rp[1] = '\0';
        repls[i] = rp;
    }
    for(i=0; i < nctable; i++) {
        size_t j = ctable[i].c;
        efree(repls[j]);
        repls[j] = ctable[i].s;
    }
    for(i = 128; i < 256; i++) {
        rp = ecalloc(hexlen+1);
        snprintf(rp, hexlen+1, "_X%2.2X", i); /* need to include null*/
        rp[hexlen] = '\0';
        repls[i] = rp;
    }
    for(i = 0; i < 256; i++) {
        lens[i] = strlen(repls[i]);
    }
    init = 1;               /* only do this initialization once */
}

/*
Convert a name to a 
form suitable for use in a
language file.
Conversion depends on l_flag.
*/
char*
codify(const char *name0)
{
    /* If the name is rooted, then elide
       the leading '/'.
    */
    if(name0[0] == '/')
	name0++;
    switch (l_flag) {
    case L_BINARY:
	return pooldup(name0);
    case L_C:
        return ccodify(name0);
    case L_F77:
        return f77codify(name0);
    case L_JAVA:
        return jcodify(name0);
    default:
        assert(0); /*no such language*/
    }
    return NULL;
}

static char*
ccodify(const char *name0)
{
    const unsigned char *cp;
    unsigned int c;
    char* name;

    if(init == 0) initcodify();
    bbClear(newname);
    cp = (const unsigned char*) name0;
    if('0' <= *cp && *cp <= '9') { /* handle initial digit, if any */
	char tmp[16];
	snprintf(tmp,sizeof(tmp),"DIGIT_%c_", *cp);
	bbCat(newname,tmp);
	cp++;
    }
    while((c=*cp++)) { /* copy name to newname, replacing special chars */
	ASSERT(c <= 256);
	bbCat(newname,repls[c]);
    }
    /* Remove leading _, if any */
    name = bbContents(newname);
    if(bbGet(newname,0) == '_') name++;
    return pooldup(name);
}

char*
cescapifychar(unsigned int c, int quote)
{
    char* s = poolalloc(4+1);
    escapifychar(c,s,quote);
    return s;
}

/**************************************************/
/* CML String Escapes */
/**************************************************/

/*
Given a character c, fill s with the character suitably escaped
for use with xml.
Caller must ensure enough space
Currently does not handle unicode
Returns s as it result.
*/

static char hexdigits[] = "0123456789ABCDEF";

static char printescapable[] = "\"&<>";
static char* printescape[] = {"quot", "amp", "lt", "gt"};

static void
xescapifychar(unsigned int c, int quote, Bytebuffer* s)
{
    if(c >= ' ' && c < '\177') {
	char* p;
	char** q;
	for(p=printescapable,q=printescape;*p;p++,q++) {if(c==*p) break;}
	if(*p) {
	    bbAppend(s,'&');
	    bbCat(s,*q);
	    bbAppend(s,';');
	} else
	    bbAppend(s,(char)c);
    } else {
	/* Do hex escape */
	    unsigned int hex1 = (c & 0x0f);
	    unsigned int hex2 = ((c >> 4) & 0x0f);
	    bbCat(s,"&#");
	    bbAppend(s,hexdigits[hex2]);
	    bbAppend(s,hexdigits[hex1]);
	    bbAppend(s,';');
    }
}

/* Return a pool string that is s0 with all characters
   ecaped that require it.  The resulting string is not
   surrounded by quotes.
   Since the string might actually contain nulls, specify the length.
*/

char*
xescapify(char* s0, int quote, size_t len)
{
    int i;
    char* result;
    Bytebuffer* escaped = bbNew();
    for(i=0;i<len;i++) {
	xescapifychar((unsigned int)s0[i],quote,escaped);
    }
    result = pooldup(bbContents(escaped));
    bbFree(escaped);
    return result;        
}

/**************************************************/
/* Java String Escapes */
/**************************************************/

/*
Given a utf16 character c,
fill s with the characters needed
to suitably escape c for use with Java.
*/

static void
jescapifychar(unsigned short c, int quote, Bytebuffer* s)
{
    /* Separate out ascii from UTF16 */
    if(c <= '\177') {
	/* Separate printables from controls */
        if(c >= ' ' && c < '\177') {
	    if (c == quote) {
		bbAppend(s,'\\');
	    }
	    bbAppend(s,(char)c);
	} else switch (c) {
	    case '\t': bbCat(s,"\\t"); break;
	    case '\b': bbCat(s,"\\b"); break;
	    case '\n': bbCat(s,"\\n"); break;
	    case '\r': bbCat(s,"\\r"); break;
	    case '\f': bbCat(s,"\\f"); break;
	    default:
		{ /* Do hex escape */
		int hex1 = (c & 0x0f);
		int hex2 = ((c >> 4) & 0x0f);
		int hex3 = ((c >> 8) & 0x0f);
		int hex4 = ((c >> 12) & 0x0f);
		bbAppend(s,'\\');
		bbAppend(s,'u');
		bbAppend(s,hexdigits[hex4]);
		bbAppend(s,hexdigits[hex3]);
		bbAppend(s,hexdigits[hex2]);
		bbAppend(s,hexdigits[hex1]);
		} break;
	}
    } else { /* Do \uxxxx escapes */
	/* Do hex escape */
	int hex1 = (c & 0x0f);
	int hex2 = ((c >> 4) & 0x0f);
	int hex3 = ((c >> 8) & 0x0f);
	int hex4 = ((c >> 12) & 0x0f);
	bbAppend(s,'\\');
	bbAppend(s,'u');
	bbAppend(s,hexdigits[hex4]);
	bbAppend(s,hexdigits[hex3]);
	bbAppend(s,hexdigits[hex2]);
	bbAppend(s,hexdigits[hex1]);
    }
}

/* Return a pool string that is s0 with all characters
   ecaped that require it.  The resulting string is not
   surrounded by quotes.
   Since the string might actually contain nulls, specify the length.
*/

char*
jescapify(char* s0, int quote, size_t len)
{
    int stat = NC_NOERR;
    int i;
    char* result = NULL;
    unsigned char* s8 = NULL;
    unsigned short* s16 = NULL; /* for storing the utf16 string */
    Bytebuffer* escaped = bbNew();
    size_t len16;

    s8 = (unsigned char*)s0;
    stat = nc_utf8_to_utf16(s8, &s16, &len16);
    if(stat != NC_NOERR) {
	derror("Cannot convert UTF8 string to UTF16: %s",s0);
	if(s16) efree(s16);
	return NULL;
    }
    for(i=0;i<len16;i++) {
	jescapifychar(s16[i],quote,escaped);
    }
    efree(s16);    
    result = pooldup(bbContents(escaped));
    bbFree(escaped);
    return result;        
}

char*
jescapifyname(char* s0)
{
    return jescapify(s0,'"',strlen(s0));
}

/*
Convert a java name that might possibly
contain utf8 characters to one that is
acceptable to the Java compiler.
Basically this means convert the printables
using ccodify (above) equivalent and then escape
all the utf chars.
*/
static char*
jcodify (const char *name)
{
    return ccodify(name);
}

/**************************************************/
/* FORTRAN does escapes differently than e.g. C */

char*
f77escapifychar(unsigned int c, char* s0)
{
    char* s = s0;
    s0[0] = '\0';
    if(c == '\'') {
	*s++ = '\''; *s++='\'';
    } else if(c >= ' ' && c < '\177') {
	*s++ = (char)c;
    } else {
	char tmp[32];
	nprintf(tmp,sizeof(tmp),"//char(%u)",c);	
	strcat(s,tmp);
	s += strlen(tmp);
    }
    *s = '\0';
    return s0;
}

void
f77quotestring(Bytebuffer* databuf)
{
    int i;
    int lastcharescaped;
    unsigned int slen = bbLength(databuf);
    unsigned char* s;

    /* Handle the empty string case */
    if(slen == 0) {
	bbCat(databuf,"char(0)");
	return;
    }

    s = (unsigned char*)ecalloc(slen+1);
    memcpy((void*)s,bbContents(databuf),slen);
    s[slen] = '\0';
    bbClear(databuf);

    lastcharescaped = 0;    
    for(i=0;i<slen;i++) {
	char tmp[32];
	unsigned int c = s[i];
	int thischarescaped = (c < ' ' || c >= '\177');
	if(i > 0) {
            if(!lastcharescaped && thischarescaped) bbAppend(databuf,'\'');
	    else if(lastcharescaped && !thischarescaped) bbCat(databuf,"//'");
	} else if(!thischarescaped)
	    bbAppend(databuf,'\'');
	f77escapifychar(c,tmp);
	if(i == 0 && thischarescaped)
            bbCat(databuf,tmp+2);
	else
            bbCat(databuf,tmp);
	lastcharescaped = thischarescaped;
    }
    if(!lastcharescaped) bbAppend(databuf,'\'');
}

static char*
f77codify(const char* s0)
{
    Bytebuffer* buf = bbNew();
    char* name;
    bbCat(buf,s0);
    f77quotestring(buf);
    name = bbDup(buf);
    bbFree(buf);
    return name;
}

/**************************************************/
/* Escape Fqn segment names by replacing
   '/' and '.' by alternate representation.
*/

char*
fqnescape(const char* s)
{
    const char* p;
    char* q;
    int c;
    int l = strlen(s);

/*
1234567
_SLASH_
_DOT__
*/
    char* newname = poolalloc(l*7+1);
    *newname = '\0';
    for(q=newname,p=s;(c=*p++);) {
#if 0
        if(c == '/' || c == '.') {
	    /* Do hex escape */
	    int hex1 = (c & 0x0f);
	    int hex2 = ((c >> 4) & 0x0f);
	    *q++ = 'X';
            *q++ = hexdigits[hex1];
            *q++ = hexdigits[hex2];
        } else
	    *q++ = c;
#else
        if(c == '/') {
	    strcat(q,"_SLASH_");
	    q += 7;
        } else if(c == '.') {
	    strcat(q,"_DOT_");
	    q += 5;
	} else {
	    *q++ = c;
	    *q = '\0';
	}
#endif
    }
    return newname;
}

/**************************************************/

/*
 * Given a pointer to a string of the form
 * 'xdd', return the corresponding hex byte
 */

int
unescapehex(const char* s)
{
    int b;
    int c1 = s[0];
    int c2 = s[1];
    if(strchr(HEXCHARS,c1) == NULL
       || strchr(HEXCHARS,c2) == NULL)
	return -1;
    b = 0;
    if(c1 < 'a') c1 = (c1 - 'A') + 'a';/* lowercase */
    if(c1 <= '9') b = (c1 - '0') << 4;
    else b = ((c1 - 'a')+10) << 4;
    if(c2 < 'a') c2 |= (c2 - 'A') + 'a';/* lowercase */
    if(c2 <= '9') b = (c2 - '0');
    else b |= ((c2 - 'a')+10);
    return b;
}

/*
 * Given a pointer to a string of the form
 * 'ddd', return the corresponding 
 * unsigned octal byte
 */

int
unescapeoct(const char* s)
{
    int b;
    int c1 = s[0];
    int c2 = s[1];
    int c3 = s[2];
    if(strchr(OCTCHARS,c1) == NULL
       || strchr(OCTCHARS,c2) == NULL
       || strchr(OCTCHARS,c3) == NULL)
	return -1;
    b = (c1 - '0') << 6;
    b |= (c2 - '0') << 3;
    b |= (c3 - '0');
    return b;
}

/*
 * "Un-escapes" valid escape sequences in yystring (read by lex) into the
 * appropriate unescaped characters.  For example, the two character
 * sequence "\t" in yystring would be converted into a single tab character.
 * On return, termstring is nul terminated.
 * Watch out for embedded nuls and utf-8 characters.
 * Return # of characters written.
 * Note that the escape handling for identifiers is different
 * than for string constants
 */

int
unescape(
     const char *yytext, /* text to unescape */
     int yyleng, /* length of yytext */
     int isident, /* Is this an identifier? */
     char** sp /* Return the unescaped version of yytext */ 
     )
{
    char* s = NULL; /* unescaped string */
    const char *t, *tend;
    char* p;
    int b;

    s = (char*)emalloc(yyleng+1);
    memcpy(s,yytext,yyleng);
    s[yyleng] = '\0';

    /* translate "\" escapes, e.g. "\t" to tab character  */
    t = yytext;
    tend = t + yyleng;
    p = s;
    while(*t && t < tend) {
	if (*t == '\\') {
	    t++;
	    switch (*t) {
	      case 'a':
		*p++ = ('\007'); t++; /* will use '\a' when STDC */
		break;
	      case 'b':
		*p++ = ('\b'); t++;
		break;
	      case 'f':
		*p++ = ('\f'); t++;
		break;
	      case 'n':
		*p++ = ('\n'); t++;
		break;
	      case 'r':
		*p++ = ('\r'); t++;
		break;
	      case 't':
		*p++ = ('\t'); t++;
		break;
	      case 'v':
		*p++ = ('\v'); t++;
		break;
	      case '\\':
		*p++ = ('\\'); t++;
		break;
	      case '?':
		*p++ = ('\177'); t++;
		break;
	      case '\'':
		*p++ = ('\''); t++;
		break;
	      case '\"':
		*p++ = ('\"'); t++;
		break;
	      /* Hex or oct constants not allowed in identifiers */
	      case 'x':
		if(!isident) {
		    /* t now points to hex */
		    b = unescapehex(t);
		    t += 2;
		} else
		    b = *t++;
	        *p++ = ((char)b);
		break;
	      case '0': case '1': case '2': case '3':
	      case '4': case '5': case '6': case '7':
		if(!isident) {
		    /* t now points to octal */
		    b = unescapeoct(t);
		    if(b < 0) {
		        derror("Bad octal constant: %s",yytext);
		        b = 0;
		    }
		    t += 3;
		} else
		    b = *t++;
 		*p++ = ((char)b);
		break;
	      default:
		*p++ = (*t); t++;
		break;
	    }
	} else {
	    *p++ = (*t); t++;
	}
    }
    *p = '\0';
    if(sp) *sp = s;
    return (p-s);
}


		
static int
ishex(int c)
{
    return ((c >= 'a' && c <= 'z')
	    || (c >= 'A' && c <= 'Z')
	    || (c >= '0' && c <= '9'));
}

static int
isoct(int c)
{
    return ((c >= '0' && c <= '7'));
}

/**
Do equivalent of strchr, but taking escapes into account
Set octhex to true if the string might contain \ddd or \xdd.
Note that escaped versions of the stop character are not
considered to match.
WARNING: if char is not found, then return pointer to the
trailing nul char
*/

char*
esc_strchr(char* s, int stopc, int octhex)
{
    char* p;
    int c;

    for(p=s;(c=*p);) {
	if(c == '\\') {
	    if(p[1] == '\0' && c == stopc)
	        return p; /* special case of '\\' at end of string */
	    p++;
	    c = *p;
/* FiX: allow utf8 digits? */
	    if(octhex && (c == 'x' || c == 'X')
	       && ishex(p[1]) && ishex(p[2]))
		    p += 3;
	    else if(octhex && isoct(c) && isoct(p[1]) && isoct(p[2]))
		    p += 3;
	    else /* treat like just '\<c>' */
	        p++;
	} else if(stopc == c)
	    break;
	else
	    p++;
    }
    return p; /* not found */
}
