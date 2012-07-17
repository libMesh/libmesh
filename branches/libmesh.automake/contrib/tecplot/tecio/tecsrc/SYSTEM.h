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
/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2008 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/
#if defined EXTERN
#undef EXTERN
#endif
#if defined SYSTEMMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
EXTERN int       OpenFileListGetCount(void);
EXTERN char     *GetLongFileName(const char *FileName);
EXTERN Boolean_t VerifyToOverwriteFile(const char *FName);
EXTERN Boolean_t IsValidDirectory(const char *FName);
EXTERN Boolean_t FileExists(const char *F,
                            Boolean_t   ShowErr);
EXTERN Boolean_t IsOkFNameChar(unsigned char ch);
EXTERN void ErrFName(const char *FName);
EXTERN Boolean_t IsValidFileName(const char *FileName,
                                 Boolean_t   IsReading,
                                 Boolean_t   ShowError);
EXTERN Boolean_t ResizeFile(FILE   *File,
                            Int64_t Length);
EXTERN Boolean_t Close_File(FILE     **F,
                            Boolean_t  ShowErr);
EXTERN Boolean_t Open_File(FILE       **F,
                           const char *FName,
                           Boolean_t  IsReading,
                           Boolean_t  IsAppending,
                           Boolean_t  ForceOpen,
                           Boolean_t  ShowErr,
                           Boolean_t  IsAscii);

