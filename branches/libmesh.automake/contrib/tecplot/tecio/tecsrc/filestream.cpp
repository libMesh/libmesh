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
#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
 *****************************************************************
 *****************************************************************
 *******                                                  ********
 ******     Copyright (C) 1988-2008 Tecplot, Inc.          *******
 *******                                                  ********
 *****************************************************************
 *****************************************************************
*/

#define FILESTREAMMODULE

#include "GLOBAL.h"
#include "TASSERT.h"
#include "ALLOC.h"
#include "SYSTEM.h"
#include "FILESTREAM.h"

/**
 */
FileStream_s *FileStreamAlloc(FILE      *File,
                              Boolean_t  IsByteOrderNative)
{
    REQUIRE(VALID_REF(File) || File == NULL);

    FileStream_s *Result = ALLOC_ITEM(FileStream_s, "FileStream");
    if (Result != NULL)
    {
        Result->File              = File;
        Result->IsByteOrderNative = IsByteOrderNative;
    }

    ENSURE(VALID_REF(Result) || Result == NULL);
    return Result;
}

/**
 */
void FileStreamDealloc(FileStream_s **FileStream)
{
    REQUIRE(VALID_REF(FileStream));
    REQUIRE(VALID_REF(*FileStream) || *FileStream == NULL);

    if (*FileStream != NULL)
    {
        FREE_ITEM(*FileStream, "FileStream");
        *FileStream = NULL;
    }

    ENSURE(*FileStream == NULL);
}
