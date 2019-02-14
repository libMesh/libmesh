/*********************************************************************
 *   Copyright 2007, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   Russ Rew
 *********************************************************************/
#ifndef _INDENT_H
#define _INDENT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Handle nested group indentation */
extern void indent_init();	/* initialize indent to zero */
extern void indent_out();	/* output current indent */
extern void indent_more();	/* increment current indent */
extern void indent_less();	/* decrement current indent */
extern int  indent_get();	/* return current indent */

#ifdef __cplusplus
}
#endif

#endif /* _INDENT_H */
