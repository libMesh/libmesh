/*********************************************************************
 *   Copyright 2010, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header$
 *********************************************************************/

#ifndef OCLOG_H
#define OCLOG_H

#define OCENVFLAG "OCLOGFILE"

/* Suggested tag values */
#define OCLOGNOTE 0
#define OCLOGWARN 1
#define OCLOGERR 2
#define OCLOGDBG 3

extern void ocloginit(void);
extern int ocsetlogging(int tf);
extern int oclogopen(const char* file);
extern void oclogclose(void);

/* The tag value is an arbitrary integer */
extern void oclog(int tag, const char* fmt, ...);
extern void oclogtext(int tag, const char* text);
extern void oclogtextn(int tag, const char* text, size_t count);

/* Provide printable names for tags */
extern void oclogsettags(char** tagset, char* dfalt);

#endif /*OCLOG_H*/
