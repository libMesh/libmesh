/*
 * NOTICE and LICENSE for Tecplot Input/Output Library (TecIO) - OpenFOAM
 *
 * Copyright (C) 1988-2009 Tecplot, Inc.  All rights reserved worldwide.
 *
 * Tecplot hereby grants OpenCFD limited authority to distribute without
 * alteration the source code to the Tecplot Input/Output library, known 
 * as TecIO, as part of its distribution of OpenFOAM and the 
 * OpenFOAM_to_Tecplot converter.  Users of this converter are also hereby
 * granted access to the TecIO source code, and may redistribute it for the
 * purpose of maintaining the converter.  However, no authority is granted
 * to alter the TecIO source code in any form or manner.
 *
 * This limited grant of distribution does not supersede Tecplot, Inc.'s 
 * copyright in TecIO.  Contact Tecplot, Inc. for further information.
 * 
 * Tecplot, Inc.
 * 3535 Factoria Blvd, Ste. 550
 * Bellevue, WA 98006, USA
 * Phone: +1 425 653 1200
 * http://www.tecplot.com/
 *
 */
#ifndef Q_MSG_H
#define Q_MSG_H
/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2008 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
#if defined EXTERN
#undef EXTERN
#endif
#if defined Q_MSGMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

#define MAX_STATUS_LINE_MSG_LEN 255

#include "TranslatedString.h"

EXTERN Boolean_t WrapString(const char  *OldString,
                            char       **NewString);
EXTERN void Warning(tecplot::strutil::TranslatedString Format,
                    ...); /* zero or more arguments */
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
EXTERN void ErrMsg(tecplot::strutil::TranslatedString Format,
                   ...); /* zero or more arguments */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE
#endif
#if !defined ENGINE
#if defined MOTIF
#endif
#endif
#if !defined ENGINE
#endif
#if defined Q_MSGMODULE
#else
#endif
#endif // TECPLOTKERNEL

#endif // Q_MSG_H
