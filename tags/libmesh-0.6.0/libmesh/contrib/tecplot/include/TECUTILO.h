#if defined EXTERN
#undef EXTERN
#endif
#if defined TECUTILOMODULE
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
LINKTOADDON Boolean_t STDCALL TecUtilInterfaceWinAddPreMsgFn(PreTranslateMessage_pf PreTranslateMessageProc);
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



LINKTOADDON Boolean_t STDCALL   TecUtilTimerAddCallback(UInt32_t               Interval,
                                                        ArbParam_t             ClientData,
                                                        AddOnTimerCallback_pf  Callback);

/*
 * Use Tecplot's file dialogs to get a filename(s) to read or write.
 */
LINKTOADDON Boolean_t STDCALL TecUtilDialogGetFileName (SelectFileOption_e   DialogOption,
                                                        char               **FileName,
                                                        const char          *FileTypeName,
                                                        const char          *DefaultFileName,
                                                        const char          *DefaultFilter);
LINKTOADDON Boolean_t STDCALL TecUtilDialogGetFileNames (SelectFileOption_e   DialogOption,
                                                         StringList_pa       *FileNames,
                                                         const char          *FileTypeName,
                                                         StringList_pa       DefaultFileNames,
                                                         const char          *DefaultFilter);
/*
 * Turn on or off the automatic updating of the sensitivities of the
 * mouse mode buttons in the sidebar.
 */
LINKTOADDON void STDCALL TecUtilSidebarAutoSensitivity(Boolean_t DoAuto);
                                                   
/*
 * Set the sensitivity for a mouse mode button in the tecplot sidebar.
 * This is usually called after calling TecUtilSidebarAutoSensitivity
 * with a value of FALSE.
 * USE CAUTION WHEN ACTIVATING A BUTTON - FOR EXAMPLE, DO NOT ACTIVATE THE
 * CONTOUR BUTTON WHEN NO CONTOURING VARIABLE HAS BEEN DEFINED!
 * When you wish to return control of the sidebar back to tecplot then
 * call TecUtilSidebarAutoSensitivity with TRUE.
 */
LINKTOADDON void STDCALL TecUtilSidebarSetSensitivity(MouseButtonMode_e MouseMode,
                                                      Boolean_t         IsSensitive);




/*
 * Probe callback function and supporting query functions.  User installs callback
 * and later uses query functions to get probed information.
 */
LINKTOADDON Boolean_t STDCALL TecUtilProbeInstallCallback(ProbeDestination_pf ProbeDestination,
                                                          const char         *InformationLineText);
LINKTOADDON void STDCALL TecUtilProbeAllowCOBs(void);
LINKTOADDON double STDCALL TecUtilProbeXYGetIndValue(void);
LINKTOADDON Boolean_t STDCALL TecUtilProbeXYGetDepValue(EntIndex_t MapNum,
                                                        double    *DepValue);
LINKTOADDON EntIndex_t STDCALL TecUtilProbeXYGetSourceMap(void);
LINKTOADDON LgIndex_t STDCALL TecUtilProbeGetPointIndex(void);
LINKTOADDON CZType_e STDCALL TecUtilProbeFieldGetCZType(void);
LINKTOADDON double STDCALL TecUtilProbeFieldGetValue(EntIndex_t VarNum);
LINKTOADDON EntIndex_t STDCALL TecUtilProbeFieldGetZone(void);
LINKTOADDON IJKPlanes_e STDCALL TecUtilProbeFieldGetPlane(void);
LINKTOADDON LgIndex_t  STDCALL TecUtilProbeFieldGetCell(void);



LINKTOADDON Boolean_t STDCALL TecUtilExtractInstallCallback(ExtractDestination_pf ExtractDestination,
                                                            const char           *InformationLineText);




/*
 *  TecUtilDataSetCreate
 *
 *  Input variables are:
 *
 *  VarNames................... String list of Names of the variables.
 *                              This array list strings must be allocated and
 *                              deallocated by the calling function using the
                                TecUtilStringList functions.
 *                              Tecplot makes its own copy of this list.
 *  ResetStyle ................ TRUE=Remove existing Text/Geoms etc.
 */
LINKTOADDON Boolean_t STDCALL TecUtilDataSetCreate(const char    *DataSetTitle,
                                                   StringList_pa VarNames,
                                                   Boolean_t     ResetStyle);





LINKTOADDON VarLoadMode_e STDCALL TecUtilDataSetGetVarLoadMode(void);


/*
 *
 *   ZoneType can be one of:
 *                            ZoneType_Ordered, 
 *                            ZoneType_FETriangle
 *                            ZoneType_FEQuad
 *                            ZoneType_FETetra 
 *                            ZoneType_FEBrick  
 *
 *   If ZoneType_Ordered, then IMax,JMax,KMax are simply the 
 *   dimensions of the zone.  If FE_xxxx then IMax is the number 
 *   of data points and JMax is the number of elements.
 *   KMax is not used if ZoneType is FE_xxxx.
 *
 *   VarDataType is an array of FieldDataType_e which currently 
 *   can be set to one of:       
 *
 *           FieldDataType_Float,
 *           FieldDataType_Double,
 *           FieldDataType_LongInt,
 *           FieldDataType_ShortInt,
 *           FieldDataType_Byte,
 *           FieldDataType_Bit
 *
 *   Or you may set VarDataType to be NULL, in which case
 *   the data type of the vars in the first zone of the
 *   dataset are used.  NOTE: TecUtilDataSetAddZone cannot be
 *   called if a dataset does not already exist.  Call
 *   TecUtilDataSetCreate first.
 */
LINKTOADDON Boolean_t STDCALL TecUtilDataSetAddZone(const char      *Name,
                                                    LgIndex_t        IMax,
                                                    LgIndex_t        JMax,
                                                    LgIndex_t        KMax,
                                                    ZoneType_e       Type,
                                                    FieldDataType_e *VarDataType_Array);

/*
BEGIN ARGLIST NOTES TecUtilAddZoneX
  TecUtilArgListAppendString(ArgList, SV_NAME,             <string>);
  TecUtilArgListAppendInt(ArgList,    SV_ZONETYPE,         (LgIndex_t)<ZoneType_e>);
  TecUtilArgListAppendInt(ArgList,    SV_ZONE,             <int>);
  TecUtilArgListAppendInt(ArgList,    SV_IMAX,             <int>);
  TecUtilArgListAppendInt(ArgList,    SV_JMAX,             <int>);
  TecUtilArgListAppendInt(ArgList,    SV_KMAX,             <int>);
  TecUtilArgListAppendArray(ArgList,  SV_VARDATATYPE       (FieldDataType_e *));
END ARGLIST NOTES
*/
LINKTOADDON Boolean_t STDCALL TecUtilDataSetAddZoneX(ArgList_pa ArgList);

/*
 * Reallocate the size of a zone.   An attempt is made to copy the
 * old values into the new but no guarentees are made.  
 *
 * If the zone is reduced in size and the zone was FE then a check
 * is made on the connectivity list.  Any index references outside
 * the allowed range are set to 0 (i.e. the first datapoint).
 *
 * If the zone is expanded then the calling function is responsible
 * for inserting the new values in the expanded area.
 *
 */
LINKTOADDON Boolean_t STDCALL TecUtilZoneRealloc(EntIndex_t Zone,
                                                 LgIndex_t  NewIMaxOrNumDataPoints,
                                                 LgIndex_t  NewJMaxOrNumElements,
                                                 LgIndex_t  NewKMax);


/*
 * Add a variable to the current dataset. 
 * The variable is initialized to zero.
 *
 *   FieldDataType is an array of FieldDataType_e which currently 
 *   can be set to one of:       
 *
 *           FieldDataType_Float,
 *           FieldDataType_Double,
 *           FieldDataType_LongInt,
 *           FieldDataType_ShortInt,
 *           FieldDataType_Byte,
 *           FieldDataType_Bit
 *
 *   Or you may set FieldDataType to be NULL, in which case
 *   Tecplot tries to figure out a good type for you?
 *
 */
LINKTOADDON Boolean_t STDCALL TecUtilDataSetAddVar(const char *VarName,
                                                   FieldDataType_e *FieldDataType_Array);
 

/************************************************************************
 * AddOn callbacks.                                                     *
 *                                                                      *
 * The functions provided below allow you to add or remove callbacks to *
 * tecplot for various events that may occur during a tecplot session.  *
 ************************************************************************/

/*
 * Include a function in the list of functions to call when tecplot is
 * Considering to quit.
 */
LINKTOADDON Boolean_t STDCALL TecUtilQuitAddQueryCallback(MopupQueryAddOnCallback_pf QuitQueryCallback);

LINKTOADDON Boolean_t STDCALL TecUtilStateChangeSetMode(StateChangeAddOnCallback_pf Callback,
                                                        StateChangeMode_e           Mode);

/*
 * Include a function in the list of functions to call when either a dataset
 * has changed or the style has changed in tecplot.  See ADDON.h for the
 * parameters for the callback function.
 */
LINKTOADDON Boolean_t STDCALL TecUtilStateChangeAddCallback(StateChangeAddOnCallback_pf StateChangeCallback);

/*
 * Inform tecplot of a state change.  Currently this must be called in the following
 * situations:
 *
 *       - Launch and dismiss of modal dialogs (Windows only).
 *       - After a variable has been added and subsequently modified.
 *       - After a variable has been modified.
 *       - After TecUtilDataSetAddZone has been called and the field data has
 *         been modified (Use StateChange_ZonesAdded).
 *       - After the node map has been altered.
 *
 * The CallData parameter is required for the following state changes:
 *
 *   StateChange_VarsAltered        (A set of variables)
 *   StateChange_VarsAdded          (A set of variables)
 *   StateChange_ZonesDeleted       (A set of zones)
 *   StateChange_ZonesAdded         (A set of zones)
 *   StateChange_NodeMapsAltered    (A set of zones)
 * 
 *  See GLOBAL.h for a list of the enumerated values.  They start with StateChange_
 */
LINKTOADDON void STDCALL TecUtilStateChanged(StateChange_e StateChange,
                                             ArbParam_t    CallData);

/*
 * Include a function in the list of functions to call when a special AddOn macro command 
 * is processed.  AddOnIDString must be a unique string used to determine the appropriate
 * macro callback to use when a $!ADDONCOMMAND command is processed.
 */
LINKTOADDON Boolean_t STDCALL TecUtilMacroAddCommandCallback(const char                  *AddOnIDString,
                                                             MacroCommandAddOnCallback_pf CallbackFunction);


/*
 * If tecplot is currently "Recording" then append 
 * a $!ADDONCOMMAND command to the file.
 * AddOnIDString is a unique string that is used to determine the
 * appropriate callback function to use when the macro command is later
 * processed.
 *
 */
LINKTOADDON Boolean_t STDCALL TecUtilMacroRecordAddOnCommand(const char *AddOnIDString,
                                                             const char *Command);

 /*
  * Like TecUtilMacroRecordAddOnCommand except you can take advantage of the
  * optional RAWDATA section that can hang off of an $!ADDONCOMMAND command.
  */

LINKTOADDON Boolean_t STDCALL TecUtilMacroRecordAddOnComRaw(const char *AddOnIDString,
                                                            const char *Command,
                                                            const char *RawData);

/*
 * Send any old thing you want to a macro record file.
 */
LINKTOADDON Boolean_t STDCALL TecUtilMacroRecordRawCommand(const char *Command);




/*
 * Send an event to the low level event dispatcher in tecplot.  Use this with
 * caution.
 */
LINKTOADDON void STDCALL TecUtilDispatchWorkAreaEvent(int       I,
                                                      int       J,
                                                      int       ButtonOrKey,
                                                      Event_e   Event,
                                                      Boolean_t IsShifted,
                                                      Boolean_t IsAlted,
                                                      Boolean_t IsControlled);

/*
 * Add a user defined menu to tecplot.
 *
 * Parameter                 Notes
 *
 * MenuLabel................ Text to put on the menu option.
 * Mnemonic................. Character to underline in MenuLabel.
 *                           Use '\0' if you don't want a mnemonic.
 * Callbackvalue............ Value you want passed to MenuOptionCallback.
 * MenuOptionCallback....... Function you create.  This is called
 *                           with Callbackvalue as its only parameter.
 * NewSubOption              No longer a parameter.
 *
 */

LINKTOADDON Boolean_t STDCALL TecUtilMenuAddOption(const char            *ParentPath,
                                                   const char            *MenuLabel,
                                                   char                   Mnemonic,
                                                   DynamicMenuCallback_pf MenuOptionCallback);


/*
 * Set the sensitivity of a menu option based on its label.
 */
LINKTOADDON Boolean_t STDCALL TecUtilMenuSetSensitivity(const char      *ParentPath,
                                                         const char      *MenuLabel,
                                                         Boolean_t  IsSensitive);

/****************************************************************************
 *          DataSet Converters and Loaders.                                 *
 *                                                                          *
 * A Converter is an addon that connects to tecplot in a way such that      *
 * tecplot uses it's own file/io dialogs to read the non-Tecplot data into  *
 * tecplot and the converter is a simple function that only knows how       *
 * to convert from it's own data format into tecplot's binary format.       *
 *                                                                          *
 * Converters are registered with tecplot by calling:                       *
 *                                                                          *
 *              TecUtilImportAddConverter                                   *
 *                                                                          *
 * An loader works like this:                                               *
 *     - TecUtilMenuAddOption is called to register a callback              *
 *       (most often used to bring up a dialog) to the loader addon         *
 *                                                                          *
 *     - The loader registers itself with tecplot by calling                *
 *       TecUtilImportAddLoader                                             *
 *                                                                          *
 *     When a request is made to use the loader function:                   *
 *                                                                          *
 *        - The loader creates a new dataset and loads the data into        *
 *          tecplot by hand.                                                *
 *                                                                          *
 *        - The loader then calls TecUtilImportSetLoaderInstr               *
 *          to assign the specific instructions needed to load the data     *
 *          into the current dataset.   These instructions will be included *
 *          in any layout files written out.                                *
 ****************************************************************************/
LINKTOADDON Boolean_t STDCALL TecUtilImportAddConverter(DataSetConverter_pf DataSetConverterFunction,
                                                          const char         *DataSetConverterName,
                                                          const char         *FNameFilter);


LINKTOADDON Boolean_t STDCALL TecUtilImportAddLoader(DataSetLoader_pf                    DataSetLoaderFunction,
                                                      const char                         *DataSetLoaderName,
                                                      DynamicMenuCallback_pf              LoaderSelectionCallback,
                                                      DataSetLoaderInstructionOverride_pf InstructionOverride);


LINKTOADDON Boolean_t STDCALL TecUtilImportSetLoaderInstr(const char    *DataSetLoaderName,
                                                          StringList_pa Instructions);


                                                                     
LINKTOADDON void STDCALL TecUtilImportWriteLoaderInstr(const char   *DataSetLoaderName,
                                                        StringList_pa Instructions);
/*
 * Register version information for an addon with tecplot.  This MUST only be called from the
 * addon initialization function.
 */
LINKTOADDON void STDCALL TecUtilAddOnRegisterInfo(const char *OfficialName,
                                                         const char *Version,
                                                         const char *Author);
LINKTOADDON Boolean_t STDCALL TecUtilAddOnGetRegisteredInfo(const char *OfficialName,
                                                    char **Version,
                                                    char **Author);


/* 
 *  Extended Curve Fit Registration
 */

LINKTOADDON Boolean_t STDCALL TecUtilCurveRegisterExtCrvFit(const char   *CurveFitName,
                                  GetXYDataPointsCallback_pf              GetXYDataPointsCallback,
                                  GetProbeValueCallback_pf                GetProbeValueCallback,
                                  GetCurveInfoStringCallback_pf           GetCurveInfoStringCallback,
                                  GetCurveSettingsCallback_pf             GetCurveSettingsCallback,                                  
                                  GetAbbreviatedSettingsStringCallback_pf GetAbbreviatedSettingsStringCallback);


/*
 * Use tecplot's error message capability.  This has the following
 * features:
 *
 *     1.  Automatic word wrapping.
 *     2.  Goes to batch.log file if in batch mode.
 *     3.  Goes to stderr if graphics are not initialized.
 */
LINKTOADDON void STDCALL TecUtilDialogErrMsg(const char *Message);

/*
 * NOTE: Offset starts at 1.
 */
LINKTOADDON void STDCALL TecUtilDataValueSetByRef(FieldData_pa FD,
                                              LgIndex_t    PointIndex,
                                              double       Value);

LINKTOADDON void STDCALL TecUtilDataNodeSetByZone(EntIndex_t Zone,
                                                       LgIndex_t  Element,
                                                       LgIndex_t  Corner,
                                                       LgIndex_t  Node);

LINKTOADDON void STDCALL TecUtilDataNodeSetByRef(NodeMap_pa  NM,
                                             LgIndex_t   Element,
                                             LgIndex_t   Corner,
                                             LgIndex_t   Node);

LINKTOADDON void STDCALL TecUtilLockOn(void);

LINKTOADDON void STDCALL TecUtilLockOff(void);

LINKTOADDON void STDCALL TecUtilDialogLaunchPercentDone(const char     *Label,
                                                         Boolean_t ShowTheScale);
LINKTOADDON void STDCALL TecUtilDialogSetPercentDoneText(const char *Text);

LINKTOADDON Boolean_t STDCALL TecUtilDialogCheckPercentDone(int PercentDone);

LINKTOADDON void STDCALL TecUtilDialogDropPercentDone(void);

LINKTOADDON Boolean_t STDCALL TecUtilMacroExecuteCommand(const char *Command);

LINKTOADDON void STDCALL TecUtilInterrupt(void);


LINKTOADDON void STDCALL TecUtilGeomDelete(Geom_ID GID);
LINKTOADDON void STDCALL TecUtilTextDelete(Text_ID TID);

LINKTOADDON Boolean_t STDCALL TecUtilPickGeom(Geom_ID GID);
LINKTOADDON Boolean_t STDCALL TecUtilPickText(Text_ID TID);

LINKTOADDON Boolean_t  STDCALL TecUtilGeomIsValid(Geom_ID GID);
LINKTOADDON Boolean_t  STDCALL TecUtilTextIsValid(Text_ID TID);
LINKTOADDON char *     STDCALL TecUtilStringAlloc(int   MaxLength, /* <-activex> */
                                                  const char *DebugInfo);
LINKTOADDON void       STDCALL TecUtilStringDealloc(char **S); /* <-activex> */


/* STRING LIST FUNCTIONS */
LINKTOADDON void          STDCALL TecUtilStringListClear(StringList_pa StringList);
LINKTOADDON void          STDCALL TecUtilStringListRemoveStrings(StringList_pa StringList,
                                                                 LgIndex_t     StringNumber,
                                                                 LgIndex_t     Count);
LINKTOADDON void          STDCALL TecUtilStringListRemoveString(StringList_pa StringList,
                                                                LgIndex_t     StringNumber);
LINKTOADDON void          STDCALL TecUtilStringListDealloc(StringList_pa *StringList);
LINKTOADDON StringList_pa STDCALL TecUtilStringListAlloc(void);
LINKTOADDON Boolean_t     STDCALL TecUtilStringListAppendString(StringList_pa StringList,
                                                                const char    *String);
LINKTOADDON LgIndex_t     STDCALL TecUtilStringListGetCount(StringList_pa StringList);
LINKTOADDON char *        STDCALL TecUtilStringListGetString(StringList_pa StringList,
                                                              LgIndex_t     StringNumber);
LINKTOADDON Boolean_t     STDCALL TecUtilStringListSetString(StringList_pa StringList,
                                                             LgIndex_t     StringNumber,
                                                             const char    *String);
LINKTOADDON Boolean_t     STDCALL TecUtilStringListInsertString(StringList_pa StringList,
                                                                LgIndex_t     StringNumber,
                                                                const char    *String);
LINKTOADDON StringList_pa STDCALL TecUtilStringListCopy(StringList_pa StringList);
LINKTOADDON Boolean_t     STDCALL TecUtilStringListAppend(StringList_pa Target,
                                                          StringList_pa Source);
LINKTOADDON char *        STDCALL TecUtilStringListToNLString(StringList_pa StringList);
LINKTOADDON StringList_pa STDCALL TecUtilStringListFromNLString(const char *String);


/* * SET FUNCTIONS * */

LINKTOADDON Set_pa     STDCALL TecUtilSetAlloc(Boolean_t ShowErr);
LINKTOADDON void       STDCALL TecUtilSetDealloc(Set_pa *Set);
LINKTOADDON Boolean_t  STDCALL TecUtilSetCopy(Set_pa    DstSet,
                                              Set_pa    SrcSet,
                                              Boolean_t ShowErr);
LINKTOADDON void       STDCALL TecUtilSetClear(Set_pa Set);
LINKTOADDON Boolean_t  STDCALL TecUtilSetAddMember(Set_pa     Set,
                                               SetIndex_t Member,
                                               Boolean_t  ShowErr);
LINKTOADDON void       STDCALL TecUtilSetRemoveMember(Set_pa     Set,
                                                    SetIndex_t Member);
LINKTOADDON Boolean_t  STDCALL TecUtilSetIsMember(Set_pa     Set,
                                            SetIndex_t Member);
LINKTOADDON Boolean_t  STDCALL TecUtilSetIsEmpty(Set_pa Set);
LINKTOADDON SetIndex_t STDCALL TecUtilSetGetMemberCount(Set_pa Set);
LINKTOADDON Boolean_t  STDCALL TecUtilSetIsEqual(Set_pa Set1,
                                                Set_pa Set2);
LINKTOADDON SetIndex_t STDCALL TecUtilSetGetMember(Set_pa     Set,
                                                   SetIndex_t Position);
LINKTOADDON SetIndex_t STDCALL TecUtilSetGetPosition(Set_pa     Set,
                                                     SetIndex_t Member);
LINKTOADDON SetIndex_t STDCALL TecUtilSetGetNextMember(Set_pa     Set,
                                                       SetIndex_t Member);
#define TecUtilSetForEachMember(Member, Set) \
            for (Member = TecUtilSetGetNextMember(Set, TECUTILSETNOTMEMBER); \
                 Member != TECUTILSETNOTMEMBER; \
                 Member = TecUtilSetGetNextMember(Set, Member))

LINKTOADDON double STDCALL TecUtilConvertXPosition(CoordSys_e  OldCoordSys,
                                                   CoordSys_e  NewCoordSys,
                                                   double      OldX);
LINKTOADDON double STDCALL TecUtilConvertXDimension(CoordSys_e  OldCoordSys,
                                                    CoordSys_e  NewCoordSys,
                                                    double      OldDimension);
LINKTOADDON double STDCALL TecUtilConvertYPosition(CoordSys_e  OldCoordSys,
                                                   CoordSys_e  NewCoordSys,
                                                   double      OldY);
LINKTOADDON double STDCALL TecUtilConvertYDimension(CoordSys_e  OldCoordSys,
                                                    CoordSys_e  NewCoordSys,
                                                    double      OldDimension);
LINKTOADDON double STDCALL TecUtilConvertUnits(Units_e OldUnits,
                                               Units_e NewUnits,
                                               double  OldSize);

/*
 *
 * WARNING:  ReadBinaryData can only be called if libtec.so was built WITHOUT
 *           tecplot's memory allocation checking turned on.
 *
 */
LINKTOADDON Boolean_t STDCALL TecUtilReadBinaryData(Boolean_t       GetHeaderInfoOnly, /* <-activex> */
                                                    const char     *FName,
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

LINKTOADDON LgIndex_t STDCALL TecUtilTecIni(const char      *Title,
                                             const char      *Variables,
                                             const char      *FName,
                                             const char      *ScratchDir,
                                             LgIndex_t *Debug,
                                             LgIndex_t *VIsDouble);

LINKTOADDON LgIndex_t STDCALL TecUtilTecZne(const char      *ZoneTitle,
                                             LgIndex_t *IMx,
                                             LgIndex_t *JMx,
                                             LgIndex_t *KMx,
                                             const char      *ZFormat,
                                             const char      *DupList);

LINKTOADDON LgIndex_t STDCALL TecUtilTecDat(LgIndex_t *N,
                                            void      *FieldData_Array,
                                            LgIndex_t *IsDouble);

LINKTOADDON LgIndex_t STDCALL TecUtilTecNod(LgIndex_t *NData_Array);

LINKTOADDON LgIndex_t STDCALL TecUtilTecEnd(void);

LINKTOADDON LgIndex_t STDCALL TecUtilTecLab(const char *S);

LINKTOADDON LgIndex_t STDCALL TecUtilTecUsr(const char *S);

LINKTOADDON LgIndex_t STDCALL TecUtilTecFil(LgIndex_t *F);

LINKTOADDON LgIndex_t STDCALL TecUtilTecTxt(double    *XPos,
                                            double    *YPos,
                                            LgIndex_t *PosCoordMode,
                                            LgIndex_t *AttachToZone,
                                            LgIndex_t *Zone,
                                            LgIndex_t *BFont,
                                            LgIndex_t *FontHeightUnits,
                                            double    *FontHeight,
                                            LgIndex_t *BoxType,
                                            double    *BoxMargin,
                                            double    *BoxLineThickness,
                                            LgIndex_t *BoxColor,
                                            LgIndex_t *BoxFillColor,
                                            double    *Angle,
                                            LgIndex_t *Anchor,
                                            double    *LineSpacing,
                                            LgIndex_t *TextColor,
                                            LgIndex_t *Scope,
                                            const char*Text,
                                            const char*MacroFunctionCommand);

LINKTOADDON LgIndex_t STDCALL TecUtilTecGeo(double    *XPos,
                                            double    *YPos,
                                            double    *ZPos,
                                            LgIndex_t *PosCoordMode,
                                            LgIndex_t *AttachToZone,
                                            LgIndex_t *Zone,
                                            LgIndex_t *Color,
                                            LgIndex_t *FillColor,
                                            LgIndex_t *IsFilled,
                                            LgIndex_t *GeomType,
                                            LgIndex_t *LinePattern,
                                            double    *PatternLength,
                                            double    *LineThickness,
                                            LgIndex_t *NumEllipsePts,
                                            LgIndex_t *ArrowheadStyle,
                                            LgIndex_t *ArrowheadAttachment,
                                            double    *ArrowheadSize,
                                            double    *ArrowheadAngle,
                                            LgIndex_t *Scope,
                                            LgIndex_t *NumSegments,
                                            LgIndex_t *NumSegPts,
                                            float     *XGeomData,
                                            float     *YGeomData,
                                            float     *ZGeomData,
                                            const char*MacroFunctionCommand);


/* Text and Geom Functions */

LINKTOADDON void     STDCALL TecUtilTextGetXYPos (Text_ID  TID,
                                                  double  *XPos,
                                                  double  *YPos);
LINKTOADDON CoordSys_e    STDCALL TecUtilTextGetPositionCoordSys (Text_ID TID);
LINKTOADDON EntIndex_t    STDCALL TecUtilTextGetZoneOrMap (Text_ID TID);
LINKTOADDON Boolean_t     STDCALL TecUtilTextIsAttached (Text_ID TID);
LINKTOADDON ColorIndex_t  STDCALL TecUtilTextGetColor (Text_ID TID);
LINKTOADDON Font_e        STDCALL TecUtilTextGetFont (Text_ID TID);
LINKTOADDON double        STDCALL TecUtilTextGetHeight (Text_ID TID);
LINKTOADDON Units_e       STDCALL TecUtilTextGetSizeUnits (Text_ID TID);
LINKTOADDON TextBox_e     STDCALL TecUtilTextBoxGetType (Text_ID TID);
LINKTOADDON double        STDCALL TecUtilTextBoxGetMargin (Text_ID TID);
LINKTOADDON double        STDCALL TecUtilTextBoxGetLineThickness (Text_ID TID);
LINKTOADDON ColorIndex_t  STDCALL TecUtilTextBoxGetColor (Text_ID TID);
LINKTOADDON ColorIndex_t  STDCALL TecUtilTextBoxGetFillColor (Text_ID TID);
LINKTOADDON double        STDCALL TecUtilTextGetAngle (Text_ID TID);
LINKTOADDON TextAnchor_e  STDCALL TecUtilTextGetAnchor (Text_ID TID);
LINKTOADDON double        STDCALL TecUtilTextGetLineSpacing (Text_ID TID);
LINKTOADDON Scope_e       STDCALL TecUtilTextGetScope (Text_ID TID);
/*
 * TecUtilTextGetMacroFunctionCmd and TecUtilTextGetString return a copy
 * of the appropriate strings which the addon must call TecUtilStringDealloc to free.
 */
LINKTOADDON Boolean_t     STDCALL TecUtilTextGetMacroFunctionCmd (Text_ID   TID,
                                                                  char    **MacroFunctionCmd);
LINKTOADDON Boolean_t     STDCALL TecUtilTextGetString (Text_ID   TID,
                                                        char    **TextString);
LINKTOADDON Text_ID       STDCALL TecUtilTextGetNext (Text_ID TID);
LINKTOADDON Text_ID       STDCALL TecUtilTextGetPrev (Text_ID TID);

LINKTOADDON void          STDCALL TecUtilGeomGetXYZAnchorPos(Geom_ID  GID,
                                                             double  *XPos,
                                                             double  *YPos,
                                                             double  *ZPos);
LINKTOADDON EntIndex_t    STDCALL TecUtilGeomGetZoneOrMap (Geom_ID GID);
LINKTOADDON Boolean_t     STDCALL TecUtilGeomIsAttached (Geom_ID GID);
LINKTOADDON ColorIndex_t  STDCALL TecUtilGeomGetColor (Geom_ID GID);
LINKTOADDON ColorIndex_t  STDCALL TecUtilGeomGetFillColor (Geom_ID GID);
LINKTOADDON Boolean_t     STDCALL TecUtilGeomGetIsFilled (Geom_ID GID);
LINKTOADDON GeomForm_e    STDCALL TecUtilGeomGetType (Geom_ID GID);
LINKTOADDON LinePattern_e STDCALL TecUtilGeomGetLinePattern (Geom_ID GID);
LINKTOADDON double        STDCALL TecUtilGeomGetPatternLength (Geom_ID GID);
LINKTOADDON double        STDCALL TecUtilGeomGetLineThickness (Geom_ID GID);
LINKTOADDON SmInteger_t   STDCALL TecUtilGeomEllipseGetNumPoints (Geom_ID GID);
LINKTOADDON ArrowheadStyle_e        STDCALL TecUtilGeomArrowheadGetStyle (Geom_ID GID);
LINKTOADDON ArrowheadAttachment_e   STDCALL TecUtilGeomArrowheadGetAttach (Geom_ID GID);
LINKTOADDON double        STDCALL TecUtilGeomArrowheadGetSize (Geom_ID GID);
LINKTOADDON double        STDCALL TecUtilGeomArrowheadGetAngle (Geom_ID GID);
LINKTOADDON Scope_e       STDCALL TecUtilGeomGetScope (Geom_ID GID);
LINKTOADDON CoordSys_e    STDCALL TecUtilGeomGetPositionCoordSys (Geom_ID GID);
/*
 * TecUtilGeomGetMacroFunctionCmd returns a copy of the macro function
 * command string which the addon must call TecUtilStringDealloc to free.
 */
LINKTOADDON Boolean_t STDCALL TecUtilGeomGetMacroFunctionCmd (Geom_ID   GID,
                                                              char    **MacroFunctionCmd);

LINKTOADDON Geom_ID       STDCALL TecUtilGeomGetNext (Geom_ID GID);
LINKTOADDON Geom_ID       STDCALL TecUtilGeomGetPrev (Geom_ID GID);

LINKTOADDON void       STDCALL TecUtilTextSetXYPos (Text_ID TID,
                                                    double XPos,
                                                    double YPos);
LINKTOADDON void       STDCALL TecUtilTextSetCoordSysAndUnits (Text_ID TID,
                                                               CoordSys_e PositionCoordSys,
                                                               Units_e    HeightUnits);
LINKTOADDON void       STDCALL TecUtilTextSetZoneOrMap (Text_ID TID,
                                                        EntIndex_t ZoneOrMap);
LINKTOADDON void       STDCALL TecUtilTextSetAttached (Text_ID TID,
                                                       Boolean_t Attached);
LINKTOADDON void       STDCALL TecUtilTextSetColor (Text_ID TID,
                                                    ColorIndex_t Color);
LINKTOADDON void       STDCALL TecUtilTextSetFont (Text_ID TID,
                                                   Font_e Font);
LINKTOADDON void       STDCALL TecUtilTextSetHeight (Text_ID TID,
                                                     double Height);
LINKTOADDON void       STDCALL TecUtilTextBoxSetType (Text_ID TID,
                                                      TextBox_e TextBoxType);
LINKTOADDON void       STDCALL TecUtilTextBoxSetMargin (Text_ID TID,
                                                        double Margin);
LINKTOADDON void       STDCALL TecUtilTextBoxSetLineThickness (Text_ID TID,
                                                               double LineThickness);
LINKTOADDON void       STDCALL TecUtilTextBoxSetColor (Text_ID TID,
                                                       ColorIndex_t BoxColor);
LINKTOADDON void       STDCALL TecUtilTextBoxSetFillColor (Text_ID TID,
                                                           ColorIndex_t BoxFillColor);
LINKTOADDON void       STDCALL TecUtilTextSetAngle (Text_ID TID,
                                                    double Angle);
LINKTOADDON void       STDCALL TecUtilTextSetAnchor (Text_ID TID,
                                                     TextAnchor_e Anchor);
LINKTOADDON void       STDCALL TecUtilTextSetLineSpacing (Text_ID TID,
                                                          double LineSpacing);
LINKTOADDON void       STDCALL TecUtilTextSetScope (Text_ID TID,
                                                    Scope_e Scope);
LINKTOADDON Boolean_t  STDCALL TecUtilTextSetMacroFunctionCmd (Text_ID TID,
                                                               const char *Command);
LINKTOADDON Boolean_t  STDCALL TecUtilTextSetString (Text_ID TID,
                                                     const char *TextString);

LINKTOADDON void       STDCALL TecUtilGeomSetXYZAnchorPos(Geom_ID GID,
                                                          double  XPos,
                                                          double  YPos,
                                                          double  ZPos);
LINKTOADDON void       STDCALL TecUtilGeomSetZoneOrMap (Geom_ID GID,
                                                        EntIndex_t ZoneOrMap);
LINKTOADDON void       STDCALL TecUtilGeomSetAttached (Geom_ID GID,
                                                       Boolean_t Attached);
LINKTOADDON void       STDCALL TecUtilGeomSetColor (Geom_ID GID,
                                                     ColorIndex_t Color);
LINKTOADDON void       STDCALL TecUtilGeomSetFillColor (Geom_ID GID,
                                                         ColorIndex_t FillColor);
LINKTOADDON void       STDCALL TecUtilGeomSetIsFilled (Geom_ID GID,
                                                       Boolean_t IsFilled);
LINKTOADDON void       STDCALL TecUtilGeomSetLinePattern (Geom_ID GID,
                                                          LinePattern_e LinePattern);
LINKTOADDON void       STDCALL TecUtilGeomSetPatternLength (Geom_ID GID,
                                                            double PatternLength);
LINKTOADDON void       STDCALL TecUtilGeomSetLineThickness (Geom_ID GID,
                                                            double LineThickness);
LINKTOADDON void       STDCALL TecUtilGeomEllipseSetNumPoints (Geom_ID GID,
                                                            SmInteger_t NumEllipsePts);
LINKTOADDON void       STDCALL TecUtilGeomArrowheadSetStyle (Geom_ID GID,
                                                             ArrowheadStyle_e ArrowheadStyle);
LINKTOADDON void       STDCALL TecUtilGeomArrowheadSetAttach (Geom_ID GID,
                                                              ArrowheadAttachment_e ArrowheadAttachment);
LINKTOADDON void       STDCALL TecUtilGeomArrowheadSetSize (Geom_ID  GID,
                                                            double   ArrowheadSize);
LINKTOADDON void       STDCALL TecUtilGeomArrowheadSetAngle (Geom_ID  GID,
                                                             double   ArrowheadAngle);
LINKTOADDON void       STDCALL TecUtilGeomSetScope (Geom_ID  GID,
                                                    Scope_e  Scope);
LINKTOADDON void       STDCALL TecUtilGeomSetPositionCoordSys (Geom_ID    GID,
                                                               CoordSys_e CoordSys);
LINKTOADDON Boolean_t  STDCALL TecUtilGeomSetMacroFunctionCmd (Geom_ID  GID,
                                                               const char    *Command);
LINKTOADDON void       STDCALL TecUtilDropOpeningBanner (void);


/* Geometry/text convenience functions */

LINKTOADDON Text_ID STDCALL TecUtilTextCreate(CoordSys_e  PositionCoordSys,
                                              double      PosX, 
                                              double      PosY,
                                              Units_e     HeightUnits,
                                              double      Height,
                                              const char *Text);

LINKTOADDON Geom_ID STDCALL TecUtilGeomSquareCreate(CoordSys_e PositionCoordSys,
                                                    double     CornerX,
                                                    double     CornerY,
                                                    double     Size);

LINKTOADDON Geom_ID STDCALL TecUtilGeomCircleCreate(CoordSys_e PositionCoordSys,
                                                    double     CenterX,
                                                    double     CenterY,
                                                    double     Radius);

LINKTOADDON Geom_ID STDCALL TecUtilGeomRectangleCreate(CoordSys_e PositionCoordSys,
                                                       double     CornerX, 
                                                       double     CornerY, 
                                                       double     Width, 
                                                       double     Height);

LINKTOADDON Geom_ID STDCALL TecUtilGeomEllipseCreate(CoordSys_e PositionCoordSys,
                                                     double     CenterX, 
                                                     double     CenterY, 
                                                     double     HAxis, 
                                                     double     VAxis);

LINKTOADDON Geom_ID STDCALL TecUtilGeom2DPolylineCreate(CoordSys_e PositionCoordSys,
                                                        double    *PtsX_Array,
                                                        double    *PtsY_Array,
                                                        LgIndex_t  NumPts);

LINKTOADDON Geom_ID STDCALL TecUtilGeom3DPolylineCreate(double   *PtsX_Array,
                                                        double   *PtsY_Array,
                                                        double   *PtsZ_Array,
                                                        LgIndex_t NumPts);
LINKTOADDON Geom_ID STDCALL TecUtilGeom2DMPolyCreate(CoordSys_e PositionCoordSys,
                                                     LgIndex_t  NumPolys,
                                                     LgIndex_t  *NumPointsInPolylines_Array);

LINKTOADDON Geom_ID STDCALL TecUtilGeom3DMPolyCreate(LgIndex_t NumPolys,
                                                     LgIndex_t *NumPointsInPolylines_Array);
LINKTOADDON Geom_ID STDCALL TecUtilGeomArcCreate(CoordSys_e PositionCoordSys,
                                                 double     CenterX, 
                                                 double     CenterY, 
                                                 double     Radius, 
                                                 double     StartAngle, 
                                                 double     EndAngle);

LINKTOADDON Geom_ID STDCALL TecUtilGeom2DLineSegmentCreate(CoordSys_e PositionCoordSys,
                                                           double     X1, 
                                                           double     Y1, 
                                                           double     X2, 
                                                           double     Y2);

LINKTOADDON Geom_ID STDCALL TecUtilGeom3DLineSegmentCreate(double X1, 
                                                           double Y1, 
                                                           double Z1, 
                                                           double X2, 
                                                           double Y2, 
                                                           double Z2);
LINKTOADDON LgIndex_t STDCALL TecUtilGeomMPolyGetPolylineCnt(Geom_ID GID);
LINKTOADDON LgIndex_t STDCALL TecUtilGeomPolyGetPointCount(Geom_ID GID);
LINKTOADDON LgIndex_t STDCALL TecUtilGeomMPolyGetPointCount(Geom_ID   GID,
                                                            LgIndex_t PolyNum);
LINKTOADDON void STDCALL TecUtilGeom2DMPolyGetPoint(Geom_ID   GID,
                                                    LgIndex_t PolyNum,
                                                    LgIndex_t PointIndex,
                                                    double    *X,
                                                    double    *Y);
LINKTOADDON void STDCALL TecUtilGeom2DPolylineGetPoint(Geom_ID   GID,
                                                       LgIndex_t PointIndex,
                                                       double    *X,
                                                       double    *Y);
LINKTOADDON void STDCALL TecUtilGeom2DMPolySetPoint(Geom_ID   GID,
                                                    LgIndex_t PolyNum,
                                                    LgIndex_t PointIndex,
                                                    double    X,
                                                    double    Y);
LINKTOADDON void STDCALL TecUtilGeom2DPolylineSetPoint(Geom_ID   GID,
                                                       LgIndex_t PointIndex,
                                                       double    X,
                                                       double    Y);
LINKTOADDON void STDCALL TecUtilGeom2DMPolySetPolyline(Geom_ID   GID,
                                                       LgIndex_t PolyNum,
                                                       double    *X_Array,
                                                       double    *Y_Array);
LINKTOADDON void STDCALL TecUtilGeom3DMPolyGetPoint(Geom_ID   GID,
                                                    LgIndex_t PolyNum,
                                                    LgIndex_t PointIndex,
                                                    double    *X,
                                                    double    *Y,
                                                    double    *Z);
LINKTOADDON void STDCALL TecUtilGeom3DPolylineGetPoint(Geom_ID   GID,
                                                       LgIndex_t PointIndex,
                                                       double    *X,
                                                       double    *Y,
                                                       double    *Z);
LINKTOADDON void STDCALL TecUtilGeom3DMPolySetPoint(Geom_ID   GID,
                                                    LgIndex_t PolyNum,
                                                    LgIndex_t PointIndex,
                                                    double    X,
                                                    double    Y,
                                                    double    Z);
LINKTOADDON void STDCALL TecUtilGeom3DPolylineSetPoint(Geom_ID   GID,
                                                       LgIndex_t PointIndex,
                                                       double    X,
                                                       double    Y,
                                                       double    Z);
LINKTOADDON void STDCALL TecUtilGeom3DMPolySetPolyline(Geom_ID   GID,
                                                       LgIndex_t PolyNum,
                                                       double    *X_Array,
                                                       double    *Y_Array,
                                                       double    *Z_Array);
LINKTOADDON double STDCALL TecUtilGeomCircleGetRadius(Geom_ID GID);
LINKTOADDON void STDCALL TecUtilGeomCircleSetRadius(Geom_ID GID,
                                                    double  Radius);
LINKTOADDON double STDCALL TecUtilGeomSquareGetSize(Geom_ID GID);
LINKTOADDON void STDCALL TecUtilGeomSquareSetSize(Geom_ID GID,
                                                  double  Size);
LINKTOADDON void STDCALL TecUtilGeomRectangleGetSize(Geom_ID GID,
                                                     double  *Width,
                                                     double  *Height);
LINKTOADDON void STDCALL TecUtilGeomRectangleSetSize(Geom_ID GID,
                                                     double  Width,
                                                     double  Height);
LINKTOADDON void STDCALL TecUtilGeomEllipseGetSize(Geom_ID GID,
                                                   double  *HAxis,
                                                   double  *VAxis);
LINKTOADDON void STDCALL TecUtilGeomEllipseSetSize(Geom_ID GID,
                                                   double  HAxis,
                                                   double  VAxis);


LINKTOADDON char * STDCALL TecUtilGetCurLayoutFName(void);

/* CORE SOURCE CODE REMOVED */


LINKTOADDON void STDCALL TecUtilHelp(const char *HelpFName,
                                     Boolean_t   GoToID,
                                     int         HelpID);



LINKTOADDON Boolean_t STDCALL TecUtilDataSetLockOn(const char *LockString);
LINKTOADDON Boolean_t STDCALL TecUtilDataSetLockOff(const char *LockString);
LINKTOADDON Boolean_t STDCALL TecUtilDataSetIsLocked(char **LockString);

#endif  /* ADDON */
