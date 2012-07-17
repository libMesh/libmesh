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
#if defined EXTERN
#undef EXTERN
#endif
#if defined DATAIOMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN Boolean_t OpenBinaryFileAndCheckMagicNumber(FileStream_s **FileStream,
                                                   char          *FName,
                                                   FileOffset_t   StartOffset,
                                                   short         *IVersion);

EXTERN Boolean_t ReadDataFileHeader(FileStream_s    *FileStream,
                                    short            IVersion,
                                    Boolean_t        ShowDataIOStatus,
                                    EntIndex_t      *NumZones,
                                    EntIndex_t      *NumVars,
                                    SmInteger_t     *NumCustomLabelSets,
                                    char           **DataSetTitle,
                                    Text_s         **BaseText,
                                    Geom_s         **BaseGeom,
                                    StringList_pa  **CustomLabelBase,
                                    StringList_pa   *UserRec,
                                    AuxData_pa      *DataSetAuxData,
                                    Set_pa         **IsVarCellCentered,
                                    Boolean_t       *HasText,
                                    Boolean_t       *HasGeoms,
                                    ArrayList_pa    *ZoneSpecList,
                                    StringList_pa   *VarNames,
                                    ArrayList_pa    *VarAuxDataList, /*<AuxData_pa>[NumVars]*/
                                    Set_pa          *IsRawFNAvailable, /* classic data only */
                                    LgIndex_t      **FNNumBndryConns,  /* classic data only */
                                    DataFileType_e  *FileType);


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
