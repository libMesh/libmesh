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
#if defined GEOMMODULE
#define EXTERN
#else
#define EXTERN extern
#endif


/* * macros for checking CoordSys_e * */
#define VALID_RECTANGLE_COORDSYS(sys) \
          (((sys)==CoordSys_Frame) || \
           ((sys)==CoordSys_Grid))
#define VALID_SQUARE_COORDSYS(sys) VALID_RECTANGLE_COORDSYS((sys))
#define VALID_ELLIPSE_COORDSYS(sys) VALID_RECTANGLE_COORDSYS((sys))
#define VALID_CIRCLE_COORDSYS(sys) VALID_ELLIPSE_COORDSYS((sys))
#define VALID_IMAGE_COORDSYS(sys) VALID_RECTANGLE_COORDSYS((sys))
#define VALID_LINESEG_COORDSYS(sys) \
          (((sys)==CoordSys_Frame) || \
           ((sys)==CoordSys_Grid)  || \
           ((sys)==CoordSys_Grid3D))
#define VALID_GEOM_COORDSYS(sys) \
          (((sys)==CoordSys_Frame) || \
           ((sys)==CoordSys_Grid)  || \
           ((sys)==CoordSys_Grid3D))

#define VALID_GEOM_TYPE(geomtype) \
          ( VALID_ENUM((geomtype),GeomType_e) && \
            (geomtype)!=GeomType_LineSegs3D )

#define VALID_GEOM_FIELD_DATA_TYPE(datatype) \
          ( ( (datatype) == FieldDataType_Float ) || \
            ( (datatype) == FieldDataType_Double ) )

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
