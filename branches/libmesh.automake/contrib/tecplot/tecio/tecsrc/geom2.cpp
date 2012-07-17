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
****** Copyright (C) 1988-2008 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#define GEOM2MODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ALLOC.h"

#include "GEOM.h"
#include "TEXT.h"
#include "STRUTIL.h"
#include "GEOM2.h"

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "DATASET0.h"


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# if 0 /* 3D geometry arrowheads are not drawn at this time. */
#endif
# if 0 /* 3D geometry arrowheads are not drawn at this time. */
# endif
#if 0
#endif
#             ifndef NO_ASSERTS
#             endif
#             ifndef NO_ASSERTS
#             endif
#endif /* TECPLOTKERNEL */


FieldDataType_e GetGeomFieldDataType(Geom_s const* Geom)
{
    FieldDataType_e Result;

    REQUIRE(VALID_REF(Geom));
    REQUIRE(VALID_REF(Geom->GeomData.Generic.V1Base));

    Result = Geom->DataType;

    ENSURE(VALID_GEOM_FIELD_DATA_TYPE(Result));
    /*
     * Check that the geom's field data arrays (if they exist)
     * have the same type as the geometry.
     */
    ENSURE(IMPLICATION(VALID_REF(Geom->GeomData.Generic.V1Base), Result == GetFieldDataType(Geom->GeomData.Generic.V1Base)));
    ENSURE(IMPLICATION(VALID_REF(Geom->GeomData.Generic.V2Base), Result == GetFieldDataType(Geom->GeomData.Generic.V2Base)));
    ENSURE(IMPLICATION(VALID_REF(Geom->GeomData.Generic.V3Base), Result == GetFieldDataType(Geom->GeomData.Generic.V3Base)));

    return Result;
}
