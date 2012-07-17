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
#if defined TEXTMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

#define _TEXT_H_INCLUDED

/* These macros for checking CoordSys_e and Units of text objects (i.e., those associated with the text tool). */
#define VALID_TEXT_COORDSYS(sys)  (((sys)==CoordSys_Frame)||((sys)==CoordSys_Grid)||((sys)==CoordSys_Grid3D))
#define VALID_TEXT_UNITS(units)  (((units)==Units_Grid)||((units)==Units_Frame)||((units)==Units_Point))
#define VALID_TEXT_COORDSYS_AND_UNITS(pos_sys, size_units) \
           ( VALID_TEXT_COORDSYS((pos_sys)) && \
             VALID_TEXT_UNITS((size_units)) && \
             ! ((pos_sys) == CoordSys_Frame && (size_units) == Units_Grid) )

/* This is for any type of font in Tecplot. */
#define VALID_FONT_SIZEUNITS(units)  (((units)==Units_Grid)||((units)==Units_Frame)||((units)==Units_Point)||(units)==Units_AxisPercentage)

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if 0 /* contract template */
#endif
#if 0 /* contract template */
#endif
#endif /* TECPLOTKERNEL */
