#if defined EXTERN
#undef EXTERN
#endif
#if defined TECUTILSMODULE
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


LINKTOADDON SetValueReturnCode_e STDCALL TecUtilStyleSetLowLevel(Widget       TextFieldWidget,
                                                            double       DValue,
                                                            ArbParam_t   IValue,
                                                            ArbParam_t   SetOrOffset,
                                                            AssignOp_e   AssignModifier,
                                                            const char  *P1,
                                                            const char  *P2,
                                                            const char  *P3,
                                                            const char  *P4,
                                                            const char  *P5,
                                                            const char  *P6,
                                                            Boolean_t    OkToRecord);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilFrameSetMode (FrameMode_e NewFrameMode);
LINKTOADDON SetValueReturnCode_e STDCALL TecUtilFrameSetName(const char *Name);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetActive (Set_pa     ZoneSet,
                                                                     AssignOp_e AssignModifier);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetActive (Set_pa     XYMapSet,
                                                                 AssignOp_e AssignModifier);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilFieldSetLayer (const char *Layer,
                                                               Boolean_t   TurnOnFieldLayer);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYSetLayer (const char *Layer,
                                                            Boolean_t   TurnOnXYLayer);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilContourSetVariable (EntIndex_t NewVariable);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetMesh(const char *Attribute,
                                                        Set_pa      ZoneSet,
                                                        double      DValue,
                                                        ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetContour(const char *Attribute,
                                                           Set_pa      ZoneSet,
                                                           double      DValue,
                                                           ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetVector(const char *Attribute,
                                                          Set_pa      ZoneSet,
                                                          double      DValue,
                                                          ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetVectorIJKSkip(const char *Attribute,
                                                                 Set_pa      ZoneSet,
                                                                 LgIndex_t   Skip);
LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetScatter(const char *Attribute,
                                                           Set_pa      ZoneSet,
                                                           double      DValue,
                                                           ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetScatterIJKSkip(const char *Attribute,
                                                                  Set_pa      ZoneSet,
                                                                  LgIndex_t   Skip);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetScatterSymbolShape(const char *Attribute,
                                                                      Set_pa      ZoneSet,
                                                                      ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetShade(const char *Attribute,
                                                         Set_pa      ZoneSet,
                                                         double      DValue,
                                                         ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetBoundary(const char *Attribute,
                                                            Set_pa      ZoneSet,
                                                            double      DValue,
                                                            ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetVolumeMode(const char *Attribute,
                                                                  const char *SubAttribute,
                                                                  Set_pa      ZoneSet,
                                                                  ArbParam_t  IValue);


/* deprecated - use TecUtilZoneSetVolumeMode instead */
LINKTOADDON SetValueReturnCode_e STDCALL TecUtilZoneSetIJKMode(const char *Attribute,
                                                           const char *SubAttribute,
                                                           Set_pa      ZoneSet,
                                                           ArbParam_t  IValue);


LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetName(Set_pa      XYMapSet,
                                                             const char *NewName);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetAssignment(const char *Attribute,
                                                                   Set_pa      XYMapSet,
                                                                   double      DValue,
                                                                   ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetLine(const char *Attribute,
                                                              Set_pa      XYMapSet,
                                                              double      DValue,
                                                              ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetCurve(const char *Attribute,
                                                               Set_pa      XYMapSet,
                                                               double      DValue,
                                                               ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetSymbol(const char *Attribute,
                                                                Set_pa      XYMapSet,
                                                                double      DValue,
                                                                ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetSymbolShape(const char *Attribute,
                                                                    Set_pa      XYMapSet,
                                                                    ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetBarChart(const char *Attribute,
                                                                  Set_pa      XYMapSet,
                                                                  double      DValue,
                                                                  ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetErrorBar(const char *Attribute,
                                                                  Set_pa      XYMapSet,
                                                                  double      DValue,
                                                                  ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilXYMapSetIndices(const char *Attribute,
                                                                const char *SubAttribute,
                                                                Set_pa      XYMapSet,
                                                                ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilCurveSetExtendedSettings(EntIndex_t  XYMapNum,
                                                                         const char *Settings);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilPrintSetup(const char *Attribute,
                                                           const char *SubAttribute,
                                                           double      DValue,
                                                           ArbParam_t  IValue);

LINKTOADDON SetValueReturnCode_e STDCALL TecUtilExportSetup(const char *Attribute,
                                                            const char *SubAttribute,
                                                            double      DValue,
                                                            ArbParam_t  IValue);
LINKTOADDON SetValueReturnCode_e STDCALL TecUtilFrameSetLinking(const char  *Attribute,
                                                                ArbParam_t   IValue);
LINKTOADDON SetValueReturnCode_e STDCALL TecUtilColorMapSetBase(ContourColorMap_e BaseColorMap);
