/*********************************************************************
  *   Copyright 2016, UCAR/Unidata
  *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
  *********************************************************************/

#ifndef D4HTTP_H
#define D4HTTP_H 1

extern int curlopen(CURL** curlp);
extern void curlclose(CURL*);

extern ncerror ncd4_fetchurl(CURL*, const char*, NCbytes*, long*, struct credentials*);
extern ncerror ncd4_fetchurl_file(CURL*, const char*, FILE*, d4size_t*, long*);

extern long ncd4_fetchhttpcode(CURL* curl);

extern ncerror ncd4_fetchlastmodified(CURL* curl, char* url, long* filetime);

extern ncerror ncd4_curlopen(CURL** curlp);
extern void ncd4_curlclose(CURL* curlp);

extern ncerror ncd4_ping(const char* url);

#endif /*D4HTTP_H*/
