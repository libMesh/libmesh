***********************************************
**                   README                  **
***********************************************

To build the TecIO library and/or the pltview utility 
simply run the Runmake script in this directory.

If customization is needed it will most likely be done
in GLOBAL.h (to identify machine as 64 bit) and/or in
dataio4.c.  Just look for CRAY in dataio4.c and you
will find most of the critical areas.  Note that the
existing code defined by CRAY is quite old and has
not been in use for some time.

Each example has its own Makefile. You may have to adjust
the variables at the top of the Makefile for your platform.


ReadTec()

The ReadTec() is included in the tecio library but is
not supported by Tecplot, Inc.  ReadTec is used 
to read Tecplot binary data files (all versions at or 
older than the Tecplot version providing the tecio 
library). See tecsrc/DATAUTIL.h for more information.

The pltview example app gives an example of  using ReadTec
to read just the header from a file as well as loading all
field data from a file./*
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
