/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef D4BYTES_H
#define D4BYTES_H 1

typedef struct D4bytes {
  size_t alloc;
  size_t used;
  void* memory;
} D4bytes;

#if defined(_CPLUSPLUS_) || defined(__CPLUSPLUS__) || defined(__CPLUSPLUS)
extern "C" {
#endif

extern D4bytes* d4bytesnew(void);
extern void d4bytesfree(D4bytes*);
extern void* d4bytesalloc(D4bytes*,size_t);
extern void* d4byteszero(D4bytes*,size_t);
extern D4bytes* d4bytesconcat(D4bytes*,D4bytes*);

#define d4byteslength(d4) ((d4)->used)
#define d4bytesmemory(d4) ((d4)->memory)
#define d4bytesreset(d4) {(d4)->used = 0;}

#if defined(_CPLUSPLUS_) || defined(__CPLUSPLUS__) || defined(__CPLUSPLUS)
}
#endif

#endif /*D4BYTES_H*/
