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
 ******     Copyright (C) 1988-2008 Tecplot, Inc.          *******
 *******                                                  ********
 *****************************************************************
 *****************************************************************
*/
#if !defined FILESTREAM_h
#define FILESTREAM_h

#if defined EXTERN
#  undef EXTERN
#endif
#if defined FILESTREAMMODULE
#  define EXTERN
#else
#  define EXTERN extern
#endif

typedef struct
{
    FILE      *File;
    Boolean_t  IsByteOrderNative;
} FileStream_s;

/**
 * Creates a structure for associating an open file stream with its byte
 * order. The byte order can changed at any time.
 *
 * @param File
 *     Handle to a file which can be NULL.
 * @param IsByteOrderNative
 *     TRUE if the file's byte order is native, FALSE if foreign.
 *
 * @return
 *     An allocated structure associating an open file to its byte order.
 */
EXTERN FileStream_s *FileStreamAlloc(FILE      *File,
                                     Boolean_t  IsByteOrderNative);

/**
 * Deallocates the structure associating the file stream with the byte order.
 * This function does NOT close the file.
 *
 * @param FileStream
 *     Pointer to an open file stream or a pointer to NULL.
 */
EXTERN void FileStreamDealloc(FileStream_s **FileStream);

#endif
