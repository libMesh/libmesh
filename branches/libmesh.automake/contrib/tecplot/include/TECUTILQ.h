#if defined EXTERN
#undef EXTERN
#endif
#if defined TECUTILQMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

#if defined ADDON

/* CORE SOURCE CODE REMOVED */

/* 
BEGINMSWINSECTION
*/
#ifdef MSWIN 
/* Add new Windows only functions here */
#endif
/* 
ENDMSWINSECTION 
*/

/* 
BEGINMOTIFSECTION 
*/

#ifdef MOTIF
/* Add new Motif only functions here */
LINKTOADDON void STDCALL TecUtilInterfaceGetMotifHandles(XtAppContext *AppContext,
                                                         Widget       *MainWidget);
#endif
/* 
ENDMOTIFSECTION
*/

LINKTOADDON Boolean_t STDCALL TecUtilBlankingCheckDataPoint(EntIndex_t Zone,
                                                            LgIndex_t  PointIndex);
 
LINKTOADDON Boolean_t STDCALL TecUtilBlankingCheckFECell(EntIndex_t Zone,
                                                         LgIndex_t  CellIndex);
 
LINKTOADDON Boolean_t STDCALL TecUtilBlankingCheckIJKCell(EntIndex_t  Zone,
                                                          IJKPlanes_e ZonePlane,
                                                          LgIndex_t   CellIndex);

LINKTOADDON int           STDCALL TecUtilLockGetCount(void);
LINKTOADDON Boolean_t     STDCALL TecUtilLockIsOn(void);

LINKTOADDON LgIndex_t     STDCALL TecUtilGetTecplotVersion(void);
LINKTOADDON char *        STDCALL TecUtilTecplotGetHomeDirectory(void);

LINKTOADDON int           STDCALL TecUtilFrameGetCount(void);
LINKTOADDON FrameMode_e   STDCALL TecUtilFrameGetMode(void);

LINKTOADDON int           STDCALL TecUtilPickListGetCount (void);
LINKTOADDON PickObjects_e STDCALL TecUtilPickListGetType (int PickListIndex);
LINKTOADDON char *        STDCALL TecUtilPickListGetFrameName (int PickListIndex);
LINKTOADDON void          STDCALL TecUtilAxisGetRange(char   Axis,
                                                      short  AxisNum,
                                                      double *AxisMin,
                                                      double *AxisMax);
LINKTOADDON char          STDCALL TecUtilPickListGetAxisKind (int PickListIndex);
LINKTOADDON int           STDCALL TecUtilPickListGetAxisNumber (int PickListIndex);
LINKTOADDON EntIndex_t    STDCALL TecUtilPickListGetZoneNumber (int PickListIndex);
LINKTOADDON void          STDCALL TecUtilPickListGetZoneIndices (int PickListIndex,
                                                                LgIndex_t *IIndex,
                                                                LgIndex_t *JIndex,
                                                                LgIndex_t *KIndex);
LINKTOADDON EntIndex_t    STDCALL TecUtilPickListGetXYMapNumber (int PickListIndex);
LINKTOADDON LgIndex_t     STDCALL TecUtilPickListGetXYMapIndex (int PickListIndex);
LINKTOADDON Text_ID       STDCALL TecUtilPickListGetText (int PickListIndex);
LINKTOADDON Geom_ID       STDCALL TecUtilPickListGetGeom (int PickListIndex);
LINKTOADDON void          STDCALL TecUtilPickListGetGeomInfo (int PickListIndex,
                                                                SmInteger_t *PolylineNum,
                                                                LgIndex_t   *PointIndex);

LINKTOADDON void STDCALL TecUtilVarGetMinMax(EntIndex_t Var,
                                             double     *VarMin,
                                             double     *VarMax);


LINKTOADDON void STDCALL TecUtilDataFECellGetNodes(EntIndex_t Zone,
                                                 int        Face,
                                                 LgIndex_t  CellIndex,
                                                 LgIndex_t *I1,
                                                 LgIndex_t *I2,
                                                 LgIndex_t *I3,
                                                 LgIndex_t *I4);

LINKTOADDON void STDCALL TecUtilDataIJKCellGetIndices(EntIndex_t  Zone,
                                                      IJKPlanes_e Plane,
                                                      LgIndex_t   CellIndex,
                                                      LgIndex_t  *I1,
                                                      LgIndex_t  *I2,
                                                      LgIndex_t  *I3,
                                                      LgIndex_t  *I4);

/* General query */
LINKTOADDON ArbParam_t STDCALL TecUtilFieldStyleGetArbValue(EntIndex_t  Zone,
                                                            const char *S1,
                                                            const char *S2,
                                                            const char *S3);

LINKTOADDON double STDCALL TecUtilFieldStyleGetDoubleValue( EntIndex_t Zone,
                                                            const char *S1,
                                                            const char *S2,
                                                            const char *S3);

LINKTOADDON ArbParam_t STDCALL TecUtilXYMapStyleGetArbValue(EntIndex_t XYMap,
                                                            const char *S1,
                                                            const char *S2,
                                                            const char *S3);

LINKTOADDON double STDCALL TecUtilXYMapStyleGetDoubleValue( EntIndex_t XYMap,
                                                            const char *S1,
                                                            const char *S2,
                                                            const char *S3);



/*
 * Query for Unique ID's
 */
LINKTOADDON UniqueID_t STDCALL TecUtilFrameGetUniqueID(void);
LINKTOADDON UniqueID_t STDCALL TecUtilDataSetGetUniqueID(void);
LINKTOADDON UniqueID_t STDCALL TecUtilZoneGetUniqueID(EntIndex_t Zone);
LINKTOADDON UniqueID_t STDCALL TecUtilVarGetUniqueID(EntIndex_t Var);


/*
 * Get var or zone number from unique id.
 */
LINKTOADDON EntIndex_t STDCALL TecUtilVarGetNumByUniqueID(UniqueID_t UniqueID);
LINKTOADDON EntIndex_t STDCALL TecUtilZoneGetNumByUniqueID(UniqueID_t UniqueID);

/*
 * Get the variable number within the dataset.
 *
 * Var may be one of 'X','Y','Z','U','V','W','C','S'
 *
 * This returns -1 on error.
 *
 */
LINKTOADDON EntIndex_t STDCALL TecUtilVarGetNumByAssignment(char Var);

/*
 * Get the variable number within the dataset.
 *
 * Varname is a "simple" variable name
 *
 * This returns -1 on error.
 *
 */
LINKTOADDON EntIndex_t STDCALL TecUtilVarGetNumByName(const char *VarName);



/*
 * Get the pointer to the raw field data in tecplot.
 * Use this function with extreme caution.  
 *
 * DO NOT ASSUME THAT RAW DATA INTERNAL TO TECPLOT
 * REMAINS IN THE SAME LOCATION AT ALL TIMES.  ALWAYS
 * CALL THIS FUNCTION AFTER ANY EVENT WHERE
 * TECPLOT ITSELF MAY MOVE/ALTER THE RAW DATA.
 *
 * ALSO MAKE SURE AND CALL Action_SetResetFlagsOnVarValueChange
 * (See ACTION2.h) * AFTER ANY FIELD VALUES HAVE CHANGED.
 *
 * Return variables:
 *
 *    DataPtr....... The address to raw data.
 *    FieldDataType. The data type of the raw data.
 *                   currently this is one of:
 *
 *                         FieldDataType_Float,
 *                         FieldDataType_Double,
 *                         FieldDataType_LongInt,
 *                         FieldDataType_ShortInt,
 *                         FieldDataType_Byte,
 *                         FieldDataType_Bit,
 *                         FieldDataType_Invalid,
 *
 * A return value of FALSE indicates invalid dataset, zone, or
 * var parameters.
 *
 */
LINKTOADDON void STDCALL TecUtilDataValueGetRawPtr(EntIndex_t        Zone, /* <-activex> */
                                                   EntIndex_t        Var,
                                                   void            **DataPtr,
                                                   FieldDataType_e  *FieldDataType);

LINKTOADDON void STDCALL TecUtilDataNodeGetRawPtr(EntIndex_t Zone,
                                                  NodeMap_t  **NodeMapPtr);

LINKTOADDON void STDCALL TecUtilDataFaceNbrGetRawPtr(EntIndex_t Zone,
                                                     LgIndex_t  **FNPtr);

/*
 * TecUtilZoneGetName, TecUtilVarGetName, TecUtilXYMapGetName return a copy
 * of the zone and variable names. The client is responsible
 * for releasing the string with TecUtilStringDealloc().
 */
LINKTOADDON Boolean_t STDCALL TecUtilZoneGetName(EntIndex_t   Zone,
                                                 char       **ZName);
LINKTOADDON Boolean_t STDCALL TecUtilVarGetName(EntIndex_t   VarNum,
                                                char       **VName);
LINKTOADDON Boolean_t STDCALL TecUtilXYMapGetName(EntIndex_t Map,
                                                  char       **Name);


/*
 * Do all that Action_SetIJK does but also assign globals:
 * XD,YD,ZD,VU,VV,VW,VB,VC,VS,NM,SD, and FS
 */
LINKTOADDON void STDCALL TecUtilZoneGetInfo(EntIndex_t    CurZone,
                                            LgIndex_t     *IMax,
                                            LgIndex_t     *JMax,
                                            LgIndex_t     *KMax,
                                            FieldData_pa  *XVar,
                                            FieldData_pa  *YVar,
                                            FieldData_pa  *ZVar,
                                            NodeMap_pa    *NMap,
                                            FieldData_pa  *UVar,
                                            FieldData_pa  *VVar,
                                            FieldData_pa  *WVar,
                                            FieldData_pa  *BVar,
                                            FieldData_pa  *CVar,
                                            FieldData_pa  *SVar);


/*
 * TecUtilDataSetGetInfo returns a copy of the DataSetTitle
 * which the addon must call TecUtilStringDealloc to free.
 */
LINKTOADDON Boolean_t  STDCALL TecUtilDataSetGetInfo(char       **DataSetTitle,
                                                     EntIndex_t  *NumZones,
                                                     EntIndex_t  *NumVars);


LINKTOADDON Boolean_t STDCALL TecUtilDataSetRequiresSaving(void);

LINKTOADDON void STDCALL TecUtilFrameGetPosAndSize(double *X,
                                                   double *Y,
                                                   double *Width,
                                                   double *Height);

LINKTOADDON Text_ID  STDCALL TecUtilTextGetBase(void);

LINKTOADDON Geom_ID  STDCALL TecUtilGeomGetBase(void);

/*
 * Use tecplot's "probe" capability to return field values
 * at a given X,Y [Z] location.
 *
 * X,Y,Z         - Position to probe at.
 *
 * ICell         - Address of LgIndex_t variable to return the index
 *                 of the cell the data point was found in.
 *                 If StartWithLocalCell is TRUE then this must be
 *                 preset to the I index of the starting cell to
 *                 look in.
 * 
 * JCell         - See ICell
 * KCell         - See ICell
 * 
 * Plane         - Plane of cell where datapoint located.  Ignore if
 *                 SearchVolume is TRUE.
 *
 * CurZone       - Zone of cell where datapoint located.
 *
 * StartWithLocalCell - Start search in area of cell *ICell,*JCell,*KCell.
 *
 * VValue        - Array of doubles the size of the number of variables in
 *                 the dataset.  Must be allocated by calling function.
 *
 *                 NOTE2: If the frame mode is 3D and SearchVolume is FALSE, then
 *                        only X and Y are used and they are assumed to be the
 *                        X,Y coordinate in the Eye-Coordinate system.  The 
 *                        returned probe will be the location on the first surface
 *                        encountered in the Eye-Coordinate system.
 *
 * SourceZones   - Set of zones to limit the search to.  Set to NULL to search
 *                 all zones.
 *
 * SearchVolume  - Set to TRUE if X,Y,Z represent datapoint inside of a 3D volume
 *                 zone.  Set to FALSE to use X,Y only and a point on the surface
 *                 of a zone is returned.
 *
 * GetZoneOnly   - Do minimal work necessary to update the CurZone variable only.
 *
 * GetNearestPoint Return values for the nearest grid point.  This will return FALSE
 *                 if the initial probe does not fall within a cell.  (Need to use 
 *                 CheckScatterPoints() to get nearest datapoint if this fails).
 *
 */
LINKTOADDON Boolean_t STDCALL TecUtilProbeAtPosition(double          X,
                                                     double          Y,
                                                     double          Z,
                                                     LgIndex_t      *ICell,
                                                     LgIndex_t      *JCell,
                                                     LgIndex_t      *KCell,
                                                     IJKPlanes_e    *Plane,
                                                     EntIndex_t     *CurZone,
                                                     Boolean_t       StartWithLocalCell,
                                                     double         *VValue_Array,
                                                     Set_pa          SourceZones,
                                                     Boolean_t       SearchVolume,
                                                     Boolean_t       GetZoneOnly,
                                                     Boolean_t       GetNearestPoint);
  

LINKTOADDON Boolean_t       STDCALL TecUtilZoneGetEnabled(Set_pa *EnabledZones);
LINKTOADDON Boolean_t       STDCALL TecUtilVarGetEnabled(Set_pa *EnabledVars);
LINKTOADDON Boolean_t       STDCALL TecUtilZoneGetActive(Set_pa *ActiveZones);
LINKTOADDON Boolean_t       STDCALL TecUtilXYMapGetActive(Set_pa *ActiveXYMaps);
LINKTOADDON void            STDCALL TecUtilXYMapGetAssignment(EntIndex_t            XYMap,
                                                              EntIndex_t           *Zone,
                                                              EntIndex_t           *XAxisVar,
                                                              EntIndex_t           *YAxisVar,
                                                              SmInteger_t          *XAxis,
                                                              SmInteger_t          *YAxis,
                                                              FunctionDependency_e *FunctionDependency);
LINKTOADDON Boolean_t       STDCALL TecUtilZoneIsFiniteElement(EntIndex_t Zone);
LINKTOADDON Boolean_t       STDCALL TecUtilZoneIsOrdered(EntIndex_t Zone);
LINKTOADDON ZoneType_e      STDCALL TecUtilZoneGetType(EntIndex_t Zone);
LINKTOADDON double          STDCALL TecUtilDataValueGetByZoneVar(EntIndex_t   Zone,
                                                                 EntIndex_t   Var,
                                                                 LgIndex_t    PointIndex);
LINKTOADDON FieldData_pa    STDCALL TecUtilDataValueGetRef(EntIndex_t Zone,
                                                           EntIndex_t Var);
LINKTOADDON NodeMap_pa      STDCALL TecUtilDataNodeGetRef(EntIndex_t Zone);
LINKTOADDON FaceNeighbor_pa STDCALL TecUtilDataFaceNbrGetRef(EntIndex_t Zone);
LINKTOADDON FieldDataType_e STDCALL TecUtilDataValueGetRefType(FieldData_pa FieldData);

LINKTOADDON Boolean_t STDCALL TecUtilImportGetLoaderInstr(char          **DataSetReaderName,
                                                          StringList_pa *DataSetLoaderInstructions);

LINKTOADDON Boolean_t STDCALL TecUtilDialogMessageBox(const char      *Message,
                                                      MessageBoxType_e MessageBoxType);
LINKTOADDON Boolean_t STDCALL TecUtilDialogGetIndexRange(LgIndex_t   MaxRangeValue,
                                                         LgIndex_t  *Min,
                                                         LgIndex_t  *Max,
                                                         LgIndex_t  *Skip);
LINKTOADDON Boolean_t STDCALL TecUtilDialogGetVariables(const char *Instructions,
                                                        const char *TextField1Label,
                                                        const char *TextField2Label,
                                                        const char *TextField3Label,
                                                        EntIndex_t *Var1,
                                                        EntIndex_t *Var2,
                                                        EntIndex_t *Var3);

/*
 * Set DefaultText to NULL if you want the default text
 * to be blank.
 *
 * Extra care her is taken to isolate tecplot internal 
 * method for allocating memory from addon program.
 * The AddOn calling function MUST call TecUtilStringDealloc
 * to free the resulting text.
 */
LINKTOADDON Boolean_t STDCALL TecUtilDialogGetSimpleText(const char  *Instructions,
                                                         const char  *DefaultText,
                                                         char       **Text);


LINKTOADDON void STDCALL TecUtilTextBoxGetPosition(Text_ID   T,
                                                   double    *X1,
                                                   double    *Y1,
                                                   double    *X2,
                                                   double    *Y2,
                                                   double    *X3,
                                                   double    *Y3,
                                                   double    *X4,
                                                   double    *Y4);

#endif

LINKTOADDON Boolean_t STDCALL TecUtilMacroFunctionExists(const char *FunctionName);
LINKTOADDON Boolean_t STDCALL TecUtilMacroIsBatchModeActive(void);


LINKTOADDON void STDCALL TecUtilInterfaceGetDotsPerInch(double *VDotsPerInch,
                                                        double *HDotsPerInch);

LINKTOADDON int STDCALL TecUtilInterfaceGetBaseFontSize(void);


LINKTOADDON double STDCALL TecUtilDataValueGetByRef(FieldData_pa FieldData,
                                                    LgIndex_t    PointIndex);
LINKTOADDON void   STDCALL TecUtilDataValueGetMinMaxByRef(FieldData_pa  FieldData,
                                                          double       *Min,
                                                          double       *Max);
LINKTOADDON LgIndex_t STDCALL TecUtilDataNodeGetByZone(EntIndex_t Zone,
                                                       LgIndex_t  Element,
                                                       LgIndex_t  Corner);

LINKTOADDON LgIndex_t STDCALL TecUtilDataNodeGetByRef(NodeMap_pa NodeMapPtr,
                                                      LgIndex_t  Element,
                                                      LgIndex_t  Corner);
LINKTOADDON LgIndex_t STDCALL TecUtilDataFaceNbrGetByZone(EntIndex_t Zone,
                                                          LgIndex_t  Element,
                                                          LgIndex_t  Face);
LINKTOADDON LgIndex_t STDCALL TecUtilDataFaceNbrGetByRef(FaceNeighbor_pa FNPtr,
                                                         LgIndex_t       Element,
                                                         LgIndex_t       Face);
LINKTOADDON EntIndex_t STDCALL TecUtilXYMapGetCount (void);

LINKTOADDON Boolean_t STDCALL TecUtilMacroIsRecordingActive (void);

LINKTOADDON LgIndex_t STDCALL TecUtilLimitGetValue(const char *LimitString);

LINKTOADDON Boolean_t STDCALL TecUtilDataSetIsAvailable (void);

LINKTOADDON Boolean_t STDCALL TecUtilVarIsEnabled (EntIndex_t Var);

LINKTOADDON Boolean_t STDCALL TecUtilZoneIsEnabled (EntIndex_t Zone);

LINKTOADDON Boolean_t STDCALL TecUtilZoneIsActive (EntIndex_t Zone);

LINKTOADDON Boolean_t STDCALL TecUtilXYMapIsActive(EntIndex_t XYMap);

LINKTOADDON Boolean_t STDCALL TecUtilGetTempFileName(char **TempFileName);

LINKTOADDON void STDCALL TecUtilColorMapGetBasicColorRGB(ColorIndex_t BasicColor,
                                                         ColorIndex_t *Red,
                                                         ColorIndex_t *Green,
                                                         ColorIndex_t *Blue );

LINKTOADDON LgIndex_t STDCALL TecUtilColorMapNumBasicColors(void);

LINKTOADDON Boolean_t STDCALL TecUtilAutoRedrawIsActive(void);

/* CORE SOURCE CODE REMOVED */

