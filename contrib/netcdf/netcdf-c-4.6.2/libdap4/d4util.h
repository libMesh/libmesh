/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#ifndef D4UTIL_H
#define D4UTIL_H 1

#ifdef HAVE_MEMMOVE
#define d4memmove(dst,src,n) memmove(dst,src,n)
#else
#define d4memmove(dst,src,n) localmemmove(dst,src,n)
#endif

/* This is intended to be big enough to work as
   an offset/position/size for a file or a memory block.
*/
typedef unsigned long long d4size_t;

/* Define an counted memory marker */
typedef struct D4blob {d4size_t size; void* memory;} D4blob;

/**************************************************/

/* signature: void swapinline16(void* ip) */
#define swapinline16(ip) \
{ \
    union {char b[2]; unsigned short i;} u; \
    char* src = (char*)(ip); \
    u.b[0] = src[1]; \
    u.b[1] = src[0]; \
    *((unsigned short*)ip) = u.i; \
}

/* signature: void swapinline32(void* ip) */
#define swapinline32(ip) \
{ \
    union {char b[4]; unsigned int i;} u; \
    char* src = (char*)(ip); \
    u.b[0] = src[3]; \
    u.b[1] = src[2]; \
    u.b[2] = src[1]; \
    u.b[3] = src[0]; \
    *((unsigned int*)ip) = u.i; \
}

/* signature: void swapinline64(void* ip) */
#define swapinline64(ip) \
{ \
    union {char b[8]; unsigned long long i;} u; \
    char* src = (char*)(ip); \
    u.b[0] = src[7]; \
    u.b[1] = src[6]; \
    u.b[2] = src[5]; \
    u.b[3] = src[4]; \
    u.b[4] = src[3]; \
    u.b[5] = src[2]; \
    u.b[6] = src[1]; \
    u.b[7] = src[0]; \
    *((unsigned long long*)ip) = u.i; \
}

/***************************************************/
/* Define the NCD4node.data.flags */

#define HASNIL   (0) /* no flags set */
#define HASSEQ   (1) /* transitively contains sequence(s)*/
#define HASSTR   (2) /* transitively contains strings */
#define HASOPFIX (4) /* transitively contains fixed size opaques */
#define HASOPVAR (8) /* transitively contains variable size opaques */
#define LEAFSEQ (16) /* mark leaf sequences */
#define HASANY   (HASNIL|HASSEQ|HASSTR|HASOPTFIX|HASOPVAR)
/***************************************************/

extern int ncd4__testurl(const char* parth, char** basename);

#endif /*D4UTIL_H*/
