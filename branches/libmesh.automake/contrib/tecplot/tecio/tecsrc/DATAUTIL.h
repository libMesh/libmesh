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
 * DATAUTIL.h : COPYRIGHT (C)1987-2002 Tecplot, Inc.
 *                 ALL RIGHTS RESERVED
 *
 * NOTE:  THIS MODULE NOW IS PART OF THE TECPLOT SOURCE
 *        ONLY EDIT THIS IN THE MAIN TECPLOT SOURCE DIRECTORY.
 *
 *
 */
#ifndef DATAUTIL_H
#define DATAUTIL_H
#define DATAUTIL_VERSION 61

#if defined MAKEARCHIVE
extern void InitInputSpecs(void);
#endif


/*
 *
 * Read a binary tecplot datafile.
 *
 * @param GetHeaderInfoOnly
 *   Return only the header info from the datafile.
 *
 * @param FName
 *  Name of the file to read.
 *
 * @param IVersion
 *  Returns version of the input file.
 *
 * @param DataSetTitle
 *  Allocates space for and returns dataset title.
 *
 * @param NumZones
 *  Returns the number of zones.
 *
 * @param NumVars
 *  Returns the number of variables.
 *
 * @param VarNames
 *  Allocates space for and returns the var names.
 *
 * @param ZoneNames
 *  Allocates space for and returns the zone names.
 *
 * @param NumPtsI, NumPtsJ, NumPtsK
 *  Zone dimensions loaded into LgIndex_t arrays.
 *
 * @param ZoneNames
 *  Zone types loaded into ZoneType_e array.
 *
 * @param UserRec
 *  Allocates space for and returns the user records.
 *
 * @param RawDataspaceAllocated
 *  Only used if GetHeaderInfoOnly is FALSE. TRUE = calling program has alloced space for
 *  the raw data. FALSE= let ReadTec allocate space for the raw data.
 *
 * @param NodeMap
 *  Finite Element connectivity information. ReadTec
 *  will allocate the space for you if RawDataspaceAllocated is FALSE.
 *
 * @param VDataBase
 *  Raw field data loaded into double arrays.  ReadTec
 *  will allocate the space for you if RawDataspaceAllocated is
 *  FALSE.  If RawDataspaceAllocated is TRUE then ReadTec will
 *  only load the arrays that have non NULL addresses.
 *
 */
LIBFUNCTION Boolean_t STDCALL ReadTec(Boolean_t       GetHeaderInfoOnly,
                                      char           *FName,
                                      short          *IVersion,
                                      char          **DataSetTitle,
                                      EntIndex_t     *NumZones,
                                      EntIndex_t     *NumVars,
                                      StringList_pa  *VarNames,
                                      StringList_pa  *ZoneNames,
                                      LgIndex_t     **NumPtsI,
                                      LgIndex_t     **NumPtsJ,
                                      LgIndex_t     **NumPtsK,
                                      ZoneType_e    **ZoneType,
                                      StringList_pa  *UserRec,
                                      Boolean_t       RawDataspaceAllocated,
                                      NodeMap_t    ***NodeMap,
                                      double       ***VDataBase);

LIBFUNCTION void * STDCALL TecAlloc(size_t size);

LIBFUNCTION void STDCALL TecFree(void *ptr);


#endif /* !DATAUTIL_H */
