#if defined EXTERN
#undef EXTERN
#endif
#if defined TECUTILMMODULE
#define EXTERN
#else
#define EXTERN extern
#endif


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
#endif
/* 
ENDMOTIFSECTION
*/

LINKTOADDON ArgList_pa STDCALL TecUtilArgListAlloc(void);
LINKTOADDON void       STDCALL TecUtilArgListClear(ArgList_pa ArgList);
LINKTOADDON Boolean_t  STDCALL TecUtilArgListAppendInt(ArgList_pa   ArgList,
                                                       const char  *Name,
                                                       LgIndex_t    Value);
LINKTOADDON Boolean_t  STDCALL TecUtilArgListAppendDouble(ArgList_pa ArgList,
                                                          const char *Name,
                                                          double     Value);
LINKTOADDON Boolean_t  STDCALL TecUtilArgListAppendString(ArgList_pa ArgList,
                                                          const char *Name,
                                                          const char *Value);
LINKTOADDON Boolean_t  STDCALL TecUtilArgListAppendArray(ArgList_pa ArgList,
                                                         const char *Name,
                                                         const void *Value);
LINKTOADDON void       STDCALL TecUtilArgListDealloc(ArgList_pa *ArgList);

LINKTOADDON Boolean_t STDCALL TecUtilColorMapRedistControlPts(void);
LINKTOADDON Boolean_t STDCALL TecUtilColorMapCopyStandard(ContourColorMap_e ColorMap);
LINKTOADDON Boolean_t STDCALL TecUtilColorMapResetToFactory(void);

LINKTOADDON Boolean_t STDCALL TecUtilRedraw(Boolean_t DoFullDrawing);/* <help> "Redraw the current frame." */
LINKTOADDON Boolean_t STDCALL TecUtilRedrawAll(Boolean_t DoFullDrawing);
LINKTOADDON Boolean_t STDCALL TecUtilDoubleBuffer(DoubleBufferAction_e DoubleBufferAction);
LINKTOADDON Boolean_t STDCALL TecUtilDrawGraphics(Boolean_t DoDrawing);
LINKTOADDON Boolean_t STDCALL TecUtilQuit(void);
LINKTOADDON Boolean_t STDCALL TecUtilFrameCreateNew(Boolean_t UseSuppliedFrameSize,
                                                    double    XPos,
                                                    double    YPos,
                                                    double    Width,
                                                    double    Height);

LINKTOADDON Boolean_t STDCALL TecUtilFramePopByName(const char *Name);
LINKTOADDON Boolean_t STDCALL TecUtilFramePushByName(const char *Name);
LINKTOADDON Boolean_t STDCALL TecUtilFramePopByUniqueID(UniqueID_t UniqueID);
LINKTOADDON Boolean_t STDCALL TecUtilFramePushByUniqueID(UniqueID_t UniqueID);


LINKTOADDON Boolean_t STDCALL TecUtilFramePushTop(void);
LINKTOADDON Boolean_t STDCALL TecUtilFramePush(int FrameNum);
LINKTOADDON Boolean_t STDCALL TecUtilFramePop(int FrameNum);
LINKTOADDON Boolean_t STDCALL TecUtilFramePopAtPosition(double X,
                                                        double Y);
LINKTOADDON Boolean_t STDCALL TecUtilFrameDeleteTop(void);
LINKTOADDON Boolean_t STDCALL TecUtilFrameFitAllToPaper(void);

LINKTOADDON ColorIndex_t STDCALL TecUtilFrameGetBackgroundColor(void);

LINKTOADDON Boolean_t STDCALL TecUtilFrameGetName(char **Name);
LINKTOADDON Boolean_t STDCALL TecUtilStyleSetBase(StyleBase_e StyleBase);
LINKTOADDON Boolean_t STDCALL TecUtilZoneDelete(Set_pa ZoneList);
LINKTOADDON Boolean_t STDCALL TecUtilReadDataSet(ReadDataOption_e  ReadDataOption,
                                                 Boolean_t         ResetStyle,
                                                 StringList_pa     FileNamesOrInstructions,
                                                 const char       *DataSetReader,
                                                 FrameMode_e       InitialFrameMode,        
                                                 Boolean_t         IncludeText,
                                                 Boolean_t         IncludeGeom,
                                                 Boolean_t         IncludeCustomLabels,
                                                 Boolean_t         IncludeData,
                                                 Boolean_t         CollapseZonesAndVars,
                                                 Set_pa            ZonesToRead,    
                                                 VarLoadMode_e     VarLoadMode,
                                                 Set_pa            VarPositionList,     
                                                 StringList_pa     VarNameList,     
                                                 LgIndex_t         ISkip,
                                                 LgIndex_t         JSkip,
                                                 LgIndex_t         KSkip);
LINKTOADDON Boolean_t STDCALL TecUtilWriteDataSet(const char       *FName,      
                                                   Boolean_t         IncludeText,
                                                   Boolean_t         IncludeGeom,
                                                   Boolean_t         IncludeCustomLabels,
                                                   Boolean_t         IncludeData,
                                                   Set_pa            ZonesToWrite,    
                                                   Set_pa            VarsToWrite,
                                                   Boolean_t         WriteBinary,
                                                   Boolean_t         UsePointFormat,
                                                   SmInteger_t       AsciiPrecision);
LINKTOADDON Boolean_t STDCALL TecUtilCreateRectangularZone(LgIndex_t       IMax,
                                                            LgIndex_t       JMax,
                                                            LgIndex_t       KMax,
                                                            double          XMin,
                                                            double          YMin,
                                                            double          ZMin,
                                                            double          XMax,
                                                            double          YMax,
                                                            double          ZMax,
                                                            FieldDataType_e FieldDataType);
LINKTOADDON Boolean_t STDCALL TecUtilCreateCircularZone(LgIndex_t       IMax,
                                                         LgIndex_t       JMax,
                                                         LgIndex_t       KMax,
                                                         double          XOrigin,
                                                         double          YOrigin,
                                                         double          Radius,
                                                         double          ZMin,
                                                         double          ZMax,
                                                         FieldDataType_e FieldDataType);
LINKTOADDON Boolean_t STDCALL TecUtilCreateSimpleXYZone(LgIndex_t       NumPts,
                                                        double         *XValues_Array,
                                                        double         *YValues_Array,
                                                        FieldDataType_e FieldDataType);

LINKTOADDON char * STDCALL TecUtilGetBasePath(const char *FName);


LINKTOADDON Boolean_t STDCALL TecUtilPublish(const char      *FName,
                                             Boolean_t        IncludeLayoutPackage,
                                             ImageSelection_e ImageSelection);
LINKTOADDON Boolean_t STDCALL TecUtilNewLayout(void);
LINKTOADDON Boolean_t STDCALL TecUtilOpenLayout(const char    *FName,
                                                StringList_pa AltInstructionList,
                                                Boolean_t     Append);
LINKTOADDON Boolean_t STDCALL TecUtilSaveLayout(const char     *FName,
                                                Boolean_t UseRelativePaths);
/*
BEGIN ARGLIST NOTES TecUtilSaveLayoutX
  TecUtilArgListAppendInt(ArgList,    SV_INCLUDEDATA,      <boolean>);
  TecUtilArgListAppendInt(ArgList,    SV_INCLUDEPREVIEW,   <boolean>);
  TecUtilArgListAppendInt(ArgList,    SV_USERELATIVEPATHS, <boolean>);
  TecUtilArgListAppendString(ArgList, SV_FNAME,            <string>);
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilSaveLayoutX(ArgList_pa ArgList);

LINKTOADDON Boolean_t STDCALL TecUtilReadStylesheet(const char *FName,
                                                    Boolean_t   IncludePlotStyle,
                                                    Boolean_t   IncludeText,
                                                    Boolean_t   IncludeGeom,
                                                    Boolean_t   IncludeStreamPositions,
                                                    Boolean_t   IncludeContourLevels,
                                                    Boolean_t   MergeStyle,    
                                                    Boolean_t   IncludeFrameSizeAndPosition);
LINKTOADDON Boolean_t STDCALL TecUtilWriteStylesheet(const char *FName,
                                                     Boolean_t   IncludePlotStyle,
                                                     Boolean_t   IncludeText,
                                                     Boolean_t   IncludeGeom,
                                                     Boolean_t   IncludeStreamPositions,
                                                     Boolean_t   IncludeContourLevels,
                                                     Boolean_t   IncludeFactoryDefaults);
LINKTOADDON Boolean_t STDCALL TecUtilReadColorMap(const char *FName);
LINKTOADDON Boolean_t STDCALL TecUtilRawColorMap(int           NumRawRGBValues,
                                                 ColorIndex_t *RawRValues_Array,
                                                 ColorIndex_t *RawGValues_Array,
                                                 ColorIndex_t *RawBValues_Array);
LINKTOADDON Boolean_t STDCALL TecUtilWriteColorMap(const char *FName);

LINKTOADDON Boolean_t STDCALL TecUtilExport(Boolean_t Append);
LINKTOADDON Boolean_t STDCALL TecUtilWorkViewFitSelectFrames(void);
LINKTOADDON Boolean_t STDCALL TecUtilWorkViewFitAllFrames(void);
LINKTOADDON Boolean_t STDCALL TecUtilWorkViewFitPaper(void);
LINKTOADDON Boolean_t STDCALL TecUtilWorkViewMaximize(void);
LINKTOADDON Boolean_t STDCALL TecUtilWorkViewLastView(void);
LINKTOADDON Boolean_t STDCALL TecUtilWorkViewZoom(double X1, 
                                                  double Y1, 
                                                  double X2, 
                                                  double Y2);
LINKTOADDON Boolean_t STDCALL TecUtilWorkViewTranslate(double X,
                                                       double Y);

LINKTOADDON Boolean_t STDCALL TecUtilViewPush(void);
LINKTOADDON Boolean_t STDCALL TecUtilViewPaste(void);
LINKTOADDON Boolean_t STDCALL TecUtilViewCopy(void);
LINKTOADDON Boolean_t STDCALL TecUtilViewLast(void);
LINKTOADDON Boolean_t STDCALL TecUtilViewZoom(double X1, 
                                              double Y1, 
                                              double X2, 
                                              double Y2);
LINKTOADDON Boolean_t STDCALL TecUtilViewTranslate(double X,
                                                   double Y);
LINKTOADDON Boolean_t STDCALL TecUtilViewCenter(void);
LINKTOADDON Boolean_t STDCALL TecUtilViewScale(double Scale);
LINKTOADDON Boolean_t STDCALL TecUtilViewAxisFit(char Axis,
                                                 short AxisNum);
LINKTOADDON Boolean_t STDCALL TecUtilViewDataFit(void);
LINKTOADDON Boolean_t STDCALL TecUtilViewFit(void);


LINKTOADDON Boolean_t STDCALL TecUtilReset3DAxes(void);
LINKTOADDON Boolean_t STDCALL TecUtilReset3DScaleFactors(void);
LINKTOADDON Boolean_t STDCALL TecUtilPrint(void);

LINKTOADDON Boolean_t STDCALL TecUtilPickAtPosition(double X, 
                                                    double Y, 
                                                    Boolean_t CollectingObjects, 
                                                    Boolean_t DiggingForObjects);

LINKTOADDON Boolean_t STDCALL TecUtilPickDeselectAll (void);
LINKTOADDON Boolean_t STDCALL TecUtilPickAddZones(Boolean_t CollectingObjects,
                                                  Set_pa    ZoneSet);
LINKTOADDON Boolean_t STDCALL TecUtilPickAddXYMaps(Boolean_t CollectingObjects,
                                                   Set_pa    XYMapsSet);
LINKTOADDON Boolean_t STDCALL TecUtilPickAddAll(PickObjects_e ObjectType);
LINKTOADDON Boolean_t STDCALL TecUtilPickAddAllInRect(double         X1,
                                                      double         Y1,
                                                      double         X2,
                                                      double         Y2,
                                                      PickObjects_e  ObjectType,
                                                      const char    *Filter);
LINKTOADDON Boolean_t STDCALL TecUtilPickEdit(const char *Action);
LINKTOADDON Boolean_t STDCALL TecUtilPickCut(void);
LINKTOADDON Boolean_t STDCALL TecUtilPickCopy(void);
LINKTOADDON Boolean_t STDCALL TecUtilPickClear(void);
LINKTOADDON Boolean_t STDCALL TecUtilPickPaste(void);
LINKTOADDON Boolean_t STDCALL TecUtilPickShift(double DXPaper, 
                                               double DYPaper, 
                                               PointerStyle_e PointerStyle);
LINKTOADDON Boolean_t STDCALL TecUtilPickMagnify(double MagFactor);
LINKTOADDON Boolean_t STDCALL TecUtilPickPush(void);
LINKTOADDON Boolean_t STDCALL TecUtilPickPop(void);
LINKTOADDON Boolean_t STDCALL TecUtilPickSetMouseMode(MouseButtonMode_e MouseMode);







LINKTOADDON Boolean_t STDCALL TecUtilXYMapCopy(EntIndex_t SourceMap,
                                                    EntIndex_t DestMap);
LINKTOADDON Boolean_t STDCALL TecUtilXYMapCreate(void);
LINKTOADDON Boolean_t STDCALL TecUtilXYMapDelete(Set_pa MapsToDelete);
LINKTOADDON Boolean_t STDCALL TecUtilXYMapShiftToTop(Set_pa MapsToShift);
LINKTOADDON Boolean_t STDCALL TecUtilXYMapShiftToBottom(Set_pa MapsToShift);
LINKTOADDON Boolean_t STDCALL TecUtilViewRotate(RotateAxis_e  RotateAxis,
                                                double        RotateAmountInDegrees);
LINKTOADDON Boolean_t STDCALL TecUtilViewRotate3D(RotateAxis_e           RotateAxis,
                                                  double                 RotateAmountInDegrees,
                                                  double                 VectorX,
                                                  double                 VectorY,
                                                  double                 VectorZ,
                                                  RotateOriginLocation_e RotateOriginLocation);
LINKTOADDON Boolean_t STDCALL TecUtilReset3DOrigin(void);
LINKTOADDON Boolean_t STDCALL TecUtilSet3DEyeDistance(double EyeDistance);
/*
BEGIN ARGLIST NOTES TecUtilReset3DOriginX
  TecUtilArgListAppendInt(ArgList,    SV_ORIGINRESETLOCATION, (LgIndex_t)<OriginResetLocation_e>);
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilReset3DOriginX(ArgList_pa ArgList);

LINKTOADDON Boolean_t STDCALL TecUtilResetVectorLength(void);

LINKTOADDON Boolean_t STDCALL TecUtilContourLevelAdd(int       NumEntries, 
                                                     double   *RawData_Array, 
                                                     Boolean_t ShowTrace);

LINKTOADDON Boolean_t STDCALL TecUtilContourLevelNew(int       NumEntries, 
                                                     double   *RawData_Array, 
                                                     Boolean_t ShowTrace);

LINKTOADDON Boolean_t STDCALL TecUtilContourLevelDeleteRange(double    RangeMin, 
                                                             double    RangeMax, 
                                                             Boolean_t ShowTrace);

LINKTOADDON Boolean_t STDCALL TecUtilContourLevelReset(int NumEntries);

LINKTOADDON Boolean_t STDCALL TecUtilContourLevelDelNearest(double    Level, 
                                                            Boolean_t ShowTrace);



LINKTOADDON Boolean_t STDCALL TecUtilContourLabelAdd(double    X, 
                                                     double    Y, 
                                                     double    Z, 
                                                     Boolean_t IsAligned);

LINKTOADDON Boolean_t STDCALL TecUtilContourLabelDeleteAll(void);
LINKTOADDON Boolean_t STDCALL TecUtilStreamtraceAdd(int           NumRakePoints,
                                                    Streamtrace_e StreamType,
                                                    StreamDir_e   Direction,
                                                    double        StartXPos,
                                                    double        StartYPos,
                                                    double        StartZPos,
                                                    double        AltStartXPos,
                                                    double        AltStartYPos,
                                                    double        AltStartZPos);
LINKTOADDON Boolean_t STDCALL TecUtilStreamtraceDeleteAll(void);
LINKTOADDON Boolean_t STDCALL TecUtilStreamtraceDeleteRange(int Start,
                                                            int End);
LINKTOADDON Boolean_t STDCALL TecUtilStreamtraceSetTermLine(int     NumPoints,
                                                            double *XTermLinePts_Array,
                                                            double *YTermLinePts_Array);
LINKTOADDON Boolean_t STDCALL TecUtilStreamtraceResetDelta(void);
LINKTOADDON int STDCALL TecUtilStreamtraceGetCount(void);
LINKTOADDON void STDCALL TecUtilStreamtraceGetPos(int    StreamNumber,
                                                  double *X,
                                                  double *Y,
                                                  double *Z);

LINKTOADDON Boolean_t STDCALL TecUtilDataValueSetByZoneVar(EntIndex_t Zone,
                                                             EntIndex_t Var,
                                                             LgIndex_t  PointIndex,
                                                             double     Value);
LINKTOADDON Boolean_t STDCALL TecUtilDataAlter(const char     *Equation,
                                               Set_pa          ZoneSet,
                                               LgIndex_t       IMin,
                                               LgIndex_t       IMax,
                                               LgIndex_t       ISkip,
                                               LgIndex_t       JMin,
                                               LgIndex_t       JMax,
                                               LgIndex_t       JSkip,
                                               LgIndex_t       KMin,
                                               LgIndex_t       KMax,
                                               LgIndex_t       KSkip,
                                               FieldDataType_e DestDataType);
LINKTOADDON Boolean_t STDCALL TecUtilSmooth(EntIndex_t          Zone,
                                            EntIndex_t          SmoothVar,
                                            LgIndex_t           NumSmoothPasses,
                                            double              SmoothWeight,
                                            BoundaryCondition_e SmoothBndryCond);
LINKTOADDON Boolean_t STDCALL TecUtilWriteCurveInfo(const char     *FName,
                                                    EntIndex_t      XYMapNum,
                                                    ProcessXYMode_e ProcessXYMode);
LINKTOADDON Boolean_t STDCALL TecUtilZoneCopy(EntIndex_t ZoneUsed,
                                                   LgIndex_t  IMin,
                                                   LgIndex_t  IMax,
                                                   LgIndex_t  ISkip,
                                                   LgIndex_t  JMin,
                                                   LgIndex_t  JMax,
                                                   LgIndex_t  JSkip,
                                                   LgIndex_t  KMin,
                                                   LgIndex_t  KMax,
                                                   LgIndex_t  KSkip);
LINKTOADDON Boolean_t STDCALL TecUtilCreateMirrorZones(Set_pa SourceZones,
                                                       char   MirrorVar);
LINKTOADDON Boolean_t STDCALL TecUtilCreateStreamZones(Boolean_t ConcatenateStreams);
LINKTOADDON Boolean_t STDCALL TecUtilCreateIsoZones(void);
LINKTOADDON Boolean_t STDCALL TecUtilCreateSliceZones(void);
LINKTOADDON Boolean_t STDCALL TecUtilCreateContourLineZones(void);
LINKTOADDON Boolean_t STDCALL TecUtilCreateFEBoundary(EntIndex_t SourceZone,
                                                      Boolean_t  RemoveBlankedSurfaces);

/* DEPRECATED */
LINKTOADDON Boolean_t STDCALL TecUtilCreateSliceZone(double OriginX,
                                                     double OriginY,
                                                     double OriginZ,
                                                     double NormalX,
                                                     double NormalY,
                                                     double NormalZ);

LINKTOADDON Boolean_t STDCALL TecUtilCreateSliceZoneFromPlane(SliceSource_e SliceSource,
                                                              double        OriginX,
                                                              double        OriginY,
                                                              double        OriginZ,
                                                              double        NormalX,
                                                              double        NormalY,
                                                              double        NormalZ);


LINKTOADDON Boolean_t STDCALL TecUtilExtractFromPolyline(double   *PolylineXPts_Array,
                                                         double   *PolylineYPts_Array,
                                                         double   *PolylineZPts_Array,
                                                         LgIndex_t NumPtsInPolyline,
                                                         Boolean_t ExtractThroughVolume,
                                                         Boolean_t ExtractOnlyPointsOnPolyline,
                                                         Boolean_t IncludeDistanceVariable,
                                                         LgIndex_t NumPtsToExtractAlongPolyline,
                                                         Boolean_t ExtractToFile,
                                                         const char     *ExtractFName);
LINKTOADDON Boolean_t STDCALL TecUtilExtractFromGeom(Boolean_t ExtractOnlyPointsOnPolyline,
                                                     Boolean_t IncludeDistanceVariable,
                                                     LgIndex_t NumPtsToExtractAlongPolyline,
                                                     Boolean_t ExtractToFile,
                                                     const char     *ExtractFName);
LINKTOADDON Boolean_t  STDCALL TecUtilPolarToRectangular(Set_pa ZoneSet);
LINKTOADDON Boolean_t  STDCALL TecUtilRotate2D(Set_pa ZoneSet,
                                               double RotateAmountInDegrees,
                                               double XOrigin,
                                               double YOrigin);
LINKTOADDON Boolean_t  STDCALL TecUtilDataRotate2D(Set_pa ZoneSet,
                                                   double RotateAmountInDegrees,
                                                   double XOrigin,
                                                   double YOrigin);
LINKTOADDON Boolean_t STDCALL TecUtilAverageCellCenterData(Set_pa ZoneSet,
                                                           Set_pa VarSet);
LINKTOADDON Boolean_t STDCALL TecUtilLinearInterpolate(Set_pa             SourceZones,
                                                       EntIndex_t         DestZone,
                                                       Set_pa             VarList,
                                                       double             LinearInterpConst,
                                                       LinearInterpMode_e LinearInterpMode);
LINKTOADDON Boolean_t STDCALL TecUtilInverseDistInterpolation(Set_pa          SourceZones,
                                                              EntIndex_t      DestZone,
                                                              Set_pa          VarList,
                                                              double          InvDistExponent,
                                                              double          InvDistMinRadius,
                                                              PtSelection_e   InterpPtSelection,
                                                              LgIndex_t       InterpNPoints);
LINKTOADDON Boolean_t STDCALL TecUtilKrig(Set_pa         SourceZones,
                                           EntIndex_t     DestZone,
                                           Set_pa         VarList,
                                           double         KrigRange,
                                           double         KrigZeroValue,
                                           Drift_e        KrigDrift,
                                           PtSelection_e  InterpPtSelection,
                                           LgIndex_t      InterpNPoints);
LINKTOADDON Boolean_t STDCALL TecUtilTriangulate(Set_pa     SourceZones,
                                                  Boolean_t  DoBoundary,
                                                  Set_pa     BoundaryZones,
                                                  Boolean_t  IncludeBoundaryPts,
                                                  LgIndex_t *NumCoincidentPts,
                                                  double     TriangleKeepFactor);

/*
BEGIN ARGLIST NOTES TecUtilAnimateZonesX
  TecUtilArgListAppendInt(ArgList,    SV_START,               <int>);
  TecUtilArgListAppendInt(ArgList,    SV_END,                 <int>);
  TecUtilArgListAppendInt(ArgList,    SV_SKIP,                <int>);
  TecUtilArgListAppendInt(ArgList,    SV_CREATEMOVIEFILE,     <boolean>);
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilAnimateZonesX(ArgList_pa ArgList);
/*
BEGIN ARGLIST NOTES TecUtilAnimateXYMapsX
  TecUtilArgListAppendInt(ArgList,    SV_START,               <int>);
  TecUtilArgListAppendInt(ArgList,    SV_END,                 <int>);
  TecUtilArgListAppendInt(ArgList,    SV_SKIP,                <int>);
  TecUtilArgListAppendInt(ArgList,    SV_CREATEMOVIEFILE,     <boolean>);
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilAnimateXYMapsX(ArgList_pa ArgList);
/*
BEGIN ARGLIST NOTES TecUtilAnimateContourLevelsX
  TecUtilArgListAppendInt(ArgList,    SV_START,               <int>);
  TecUtilArgListAppendInt(ArgList,    SV_END,                 <int>);
  TecUtilArgListAppendInt(ArgList,    SV_SKIP,                <int>);
  TecUtilArgListAppendInt(ArgList,    SV_CREATEMOVIEFILE,     <boolean>);
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilAnimateContourLevelsX(ArgList_pa ArgList);
/*
BEGIN ARGLIST NOTES TecUtilAnimateIJKPlanesX
  TecUtilArgListAppendInt(ArgList,    SV_START,               <int>);
  TecUtilArgListAppendInt(ArgList,    SV_END,                 <int>);
  TecUtilArgListAppendInt(ArgList,    SV_SKIP,                <int>);
  TecUtilArgListAppendInt(ArgList,    SV_PLANES,              <int>);
  TecUtilArgListAppendInt(ArgList,    SV_CREATEMOVIEFILE,     <boolean>);
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilAnimateIJKPlanesX(ArgList_pa ArgList);
/*
BEGIN ARGLIST NOTES TecUtilAnimateIJKBlankingX
  TecUtilArgListAppendDouble(ArgList, SV_IMINFRACT,          <double>);
  TecUtilArgListAppendDouble(ArgList, SV_JMINFRACT,          <double>);
  TecUtilArgListAppendDouble(ArgList, SV_KMINFRACT,          <double>);
  TecUtilArgListAppendDouble(ArgList, SV_IMAXFRACT,          <double>);
  TecUtilArgListAppendDouble(ArgList, SV_JMAXFRACT,          <double>);
  TecUtilArgListAppendDouble(ArgList, SV_KMAXFRACT,          <double>);
  TecUtilArgListAppendDouble(ArgList, SV_IMINFRACT2,         <double>);
  TecUtilArgListAppendDouble(ArgList, SV_JMINFRACT2,         <double>);
  TecUtilArgListAppendDouble(ArgList, SV_KMINFRACT2,         <double>);
  TecUtilArgListAppendDouble(ArgList, SV_IMAXFRACT2,         <double>);
  TecUtilArgListAppendDouble(ArgList, SV_JMAXFRACT2,         <double>);
  TecUtilArgListAppendDouble(ArgList, SV_KMAXFRACT2,         <double>);
  TecUtilArgListAppendInt(ArgList,    SV_NUMSTEPS,           <double>);
  TecUtilArgListAppendInt(ArgList,    SV_CREATEMOVIEFILE,    <boolean>);
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilAnimateIJKBlankingX(ArgList_pa ArgList);
/*
BEGIN ARGLIST NOTES TecUtilAnimateStreamX
  TecUtilArgListAppendInt(ArgList,    SV_NUMCYCLES,          <int>);
  TecUtilArgListAppendInt(ArgList,    SV_STEPSPERCYCLE,      <int>);
  TecUtilArgListAppendInt(ArgList,    SV_CREATEMOVIEFILE,    <boolean>);
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilAnimateStreamX(ArgList_pa ArgList);
/*
BEGIN ARGLIST NOTES TecUtilAnimateSlicesX
  TecUtilArgListAppendInt(ArgList,    SV_START,               <int>);
  TecUtilArgListAppendInt(ArgList,    SV_END,                 <int>);
  TecUtilArgListAppendInt(ArgList,    SV_NUMSLICES,           <int>);
  TecUtilArgListAppendInt(ArgList,    SV_CREATEMOVIEFILE,     <boolean>);
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilAnimateSlicesX(ArgList_pa ArgList);
LINKTOADDON Boolean_t STDCALL TecUtilAnimateZones(EntIndex_t  StartZone,
                                                  EntIndex_t  EndZone,
                                                  EntIndex_t  ZoneSkip,
                                                  Boolean_t   CreateMoviefile,
                                                  const char *MovieFName);
LINKTOADDON Boolean_t STDCALL TecUtilAnimateXYMaps(EntIndex_t  StartMap,
                                                   EntIndex_t  EndMap,
                                                   EntIndex_t  MapSkip,
                                                   Boolean_t   CreateMoviefile,
                                                   const char *MovieFName);
LINKTOADDON Boolean_t STDCALL TecUtilAnimateContourLevels(SmInteger_t StartLevel,
                                                          SmInteger_t EndLevel,
                                                          SmInteger_t LevelSkip,
                                                          Boolean_t   CreateMoviefile,
                                                          const char *MovieFName);
LINKTOADDON Boolean_t STDCALL TecUtilAnimateIJKPlanes(char       IJOrK,
                                                      LgIndex_t  StartIndex,
                                                      LgIndex_t  EndIndex,
                                                      LgIndex_t  IndexSkip,
                                                      Boolean_t  CreateMoviefile,
                                                      const char      *MovieFName);
LINKTOADDON Boolean_t STDCALL TecUtilAnimateIJKBlanking(double      StartIMinFract,
                                                        double      StartJMinFract,
                                                        double      StartKMinFract,
                                                        double      StartIMaxFract,
                                                        double      StartJMaxFract,
                                                        double      StartKMaxFract,
                                                        double      EndIMinFract,
                                                        double      EndJMinFract,
                                                        double      EndKMinFract,
                                                        double      EndIMaxFract,
                                                        double      EndJMaxFract,
                                                        double      EndKMaxFract,
                                                        int         NumSteps,
                                                        Boolean_t   CreateMoviefile,
                                                        const char *MovieFName);
LINKTOADDON Boolean_t STDCALL TecUtilAnimateStream(int         NumStepsPerCycle,
                                                   int         NumCycles,
                                                   Boolean_t   CreateMoviefile,
                                                   const char *MovieFName);
LINKTOADDON Boolean_t STDCALL TecUtilAnimateSlices(SmInteger_t  StartSlice,
                                                   SmInteger_t  EndSlice,
                                                   SmInteger_t  NumSlices,
                                                   Boolean_t    CreateMovieFile,
                                                   const char  *MovieFName);
LINKTOADDON Boolean_t STDCALL TecUtilDelay(LgIndex_t Seconds);
LINKTOADDON Boolean_t STDCALL TecUtilMacroRunFunction(const char *QuickMacroName,
                                                      const char *MacroParameters);
LINKTOADDON Boolean_t STDCALL TecUtilDataSetSetTitle(const char *DataSetTitle);
LINKTOADDON Boolean_t STDCALL TecUtilVarRename(EntIndex_t  VarNum,
                                               const char *VarName);
LINKTOADDON Boolean_t STDCALL TecUtilZoneRename(EntIndex_t  Zone,
                                                const char *ZoneName);
LINKTOADDON Boolean_t STDCALL TecUtilSystem(const char *Command,
                                            Boolean_t Wait);
LINKTOADDON Boolean_t STDCALL TecUtilMacroPanelAddTitle(const char *Title);
LINKTOADDON Boolean_t STDCALL TecUtilDialogLaunch(Dialog_e DialogToLaunch);
LINKTOADDON Boolean_t STDCALL TecUtilDialogDrop(Dialog_e DialogToDrop);
LINKTOADDON Boolean_t STDCALL TecUtilMacroRunFile(const char *FName);
LINKTOADDON Boolean_t STDCALL TecUtilAddOnLoad(const char   *LibName,
                                               const char   *InitFunctionName,
                                               LibraryType_e AddOnType);


LINKTOADDON Boolean_t STDCALL TecUtilExportIsRecording(void);
LINKTOADDON Boolean_t STDCALL TecUtilExportFinish(void);

LINKTOADDON Boolean_t STDCALL TecUtilExportStart(void);
LINKTOADDON Boolean_t STDCALL TecUtilExportNextFrame(void);
LINKTOADDON void STDCALL TecUtilExportCancel(void);

LINKTOADDON void STDCALL TecUtilThreeDViewGetViewerAngle( double *PsiAngle,
                                                  double *ThetaAngle,
                                                  double *AlphaAngle );
LINKTOADDON void STDCALL TecUtilThreeDViewGetViewerPos(double *XPos,
                                                       double *YPos,
                                                       double *ZPos);

LINKTOADDON void STDCALL TecUtilThreeDViewGetProjection(double *FieldOfView,
                                                        double *ViewWidth,
                                                        Boolean_t *IsInPerspective);

LINKTOADDON void STDCALL TecUtilConvert3DPositionToGrid(double XPosition,
                                                        double YPosition,
                                                        double ZPosition,
                                                        double *XGridPosition,
                                                        double *YGridPosition,
                                                        double *ZGridPosition);

LINKTOADDON double STDCALL TecUtil3DViewGetNearZPlane(void);
LINKTOADDON void STDCALL TecUtilSetupTransformations(void);
LINKTOADDON void STDCALL TecUtilFrameLightweightPopStart(void);
LINKTOADDON Boolean_t STDCALL TecUtilFrameLightweightPopNext(void);
LINKTOADDON void STDCALL TecUtilFrameLightweightPopEnd(void);

/* Custom Exporters */
LINKTOADDON Boolean_t STDCALL TecUtilImageRGBBitmapCreate(BitDumpRegion_e Region);
LINKTOADDON Boolean_t STDCALL TecUtilImageIndexedBitmapCreate(BitDumpRegion_e Region,
                                                              short *RedColorTable_Array,
                                                              short *GreenColorTable_Array,
                                                              short *BlueColorTable_Array);
LINKTOADDON void      STDCALL TecUtilImageBitmapDestroy(void);
LINKTOADDON Boolean_t STDCALL TecUtilImageGetDimensions(short *Width,
                                                        short *Height);
LINKTOADDON Boolean_t STDCALL TecUtilImageRGBGetScanLine( short ScanLine,
                                                          short *Red_Array,
                                                          short *Green_Array,
                                                          short *Blue_Array);
LINKTOADDON Boolean_t STDCALL TecUtilImageIndexedGetScanLine( short ScanLine,
                                                              short *RGBIndex_Array);

LINKTOADDON void STDCALL TecUtilImageGetColorTable(Byte_t *Red_Array,
                                                   Byte_t *Green_Array,
                                                   Byte_t *Blue_Array);

LINKTOADDON Boolean_t STDCALL TecUtilImageBitmapCreateX(ArgList_pa ArgList);









