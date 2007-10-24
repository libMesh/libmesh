/* CORE SOURCE CODE REMOVED */

/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** (C) Copyright 1988-2001  by AMTEC ENGINEERING INC. *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#ifndef _GLOBAL_H
#define _GLOBAL_H

#if defined EXTERN
#undef EXTERN
#endif
#if defined Q_MAINMODULE && defined TECPLOTKERNEL
#define EXTERN
#else
#define EXTERN extern
#endif



/****************************************************************
 *                                                              *
 *                          MACROS                              *
 *                                                              *
 ****************************************************************/
#if defined TRUE
#undef TRUE
#endif
#if defined FALSE
#undef FALSE
#endif
#if defined MIN
#undef MIN
#endif
#if defined MAX
#undef MAX
#endif
#if defined ROUND
#undef ROUND
#endif
#if defined ROUND2
#undef ROUND2
#endif
#if defined TRUNC
#undef TRUNC
#endif

#define TRUE                  ((Boolean_t)1)
#define FALSE                 ((Boolean_t)0)

/****************************************************************
 *                                                              *
 *                           MACROS                             *
 *                                                              *
 ****************************************************************/
#define ABS(X)                ((X) >= 0 ? (X) : -(X) )
#define MAX(X,Y)              ((X) > (Y) ? (X) : (Y) )
#define MIN(X,Y)              ((X) < (Y) ? (X) : (Y) )
#define BESTSHOWCOLOR(X)      ((X) == White_C ? Black_C : White_C)
#define ROUND(X)              ((short)((X)+0.499))
#define ROUNDL(X)             ((LgIndex_t)((X)+0.499))
#define ROUND2(X)             ((X) >= 0 ? ((int)((X)+0.499)) : ((int)((X)-0.499)))
#define TRUNC(X)              ((short) (X))
#define RAD_TO_DEG(rad)       (180.*(rad)/PI)
#define DEG_TO_RAD(deg)       (PI*(deg)/180.)
#define CAPITAL(C)            ((islower(C) ? toupper(C) : (C)))
#define ISWHITESPACE(C)       ((C == ' ') || (C == '\t') || (C == '\n'))
#define ISSEPARATOR(C)        ((C == ' ') || (C == '\t') || (C == ','))
#define IJKINDEX(I,J,K)       ((I) + \
                               ((J)*CZData->NumIPts) + \
                               ((K)*CZData->NumIJPts))
#define IINDEX(N)             ((N) % CZData->NumIPts)
#define JINDEX(N)             (((N) % CZData->NumIJPts)/CZData->NumIPts)
#define KINDEX(N)             ((N)/CZData->NumIJPts)
#define SWITCH(A,B)           {double T = A; A = B; B = T;}
#define FPRINTFOK(x)          (Boolean_t)((x) > 0)
#define GRAPHICSARE3D(F)      ((F->FrameMode == Frame_ThreeD))

/* convenience macros for implication, P -> Q, and equivalence, P <-> Q. */
#define IMPLICATION(P,Q) (!(P) || (Q))
#define EQUIVALENCE(P,Q) ((P) == (Q))

/* suppress compiler warnings about unused parameters */
#ifndef UNUSED
# define UNUSED(param) (void)param
#endif 

/**
 * Reverses the byte order of the specified 2 byte buffer.
 *
 * param Buffer
 *     Pointer to the 2 bytes needing byte order reversal.
 */
#define REVERSE_2BYTES(Buffer) \
          { \
            char Byte0 = ((char *)(Buffer))[0]; \
            ((char *)(Buffer))[0] = ((char *)(Buffer))[1]; \
            ((char *)(Buffer))[1] = Byte0; \
          }

/**
 * Reverses the byte order of the specified 4 byte buffer.
 *
 * param Buffer
 *     Pointer to the 4 bytes needing byte order reversal.
 */
#define REVERSE_4BYTES(Buffer) \
          { \
            char Byte0 = ((char *)(Buffer))[0]; \
            char Byte1 = ((char *)(Buffer))[1]; \
            ((char *)(Buffer))[0] = ((char *)(Buffer))[3]; \
            ((char *)(Buffer))[3] = Byte0; \
            ((char *)(Buffer))[1] = ((char *)(Buffer))[2]; \
            ((char *)(Buffer))[2] = Byte1; \
          }

/**
 * Reverses the byte order of the specified 8 byte buffer.
 *
 * param Buffer
 *     Pointer to the 8 bytes needing byte order reversal.
 */
#define REVERSE_8BYTES(Buffer) \
          { \
            char Byte0 = ((char *)(Buffer))[0]; \
            char Byte1 = ((char *)(Buffer))[1]; \
            char Byte2 = ((char *)(Buffer))[2]; \
            char Byte3 = ((char *)(Buffer))[3]; \
            ((char *)(Buffer))[0] = ((char *)(Buffer))[7]; \
            ((char *)(Buffer))[7] = Byte0; \
            ((char *)(Buffer))[1] = ((char *)(Buffer))[6]; \
            ((char *)(Buffer))[6] = Byte1; \
            ((char *)(Buffer))[2] = ((char *)(Buffer))[5]; \
            ((char *)(Buffer))[5] = Byte2; \
            ((char *)(Buffer))[3] = ((char *)(Buffer))[4]; \
            ((char *)(Buffer))[4] = Byte3; \
          }


/****************************************************************
 *                                                              *
 *             ADDON MSWIN IMPORT/EXPORT DEFINITIONS            *
 *                                                              *
 ****************************************************************/
#if defined MSWIN
#  define STDCALL __stdcall
#else
#  define STDCALL
#endif /* MSWIN */

#if defined (__cplusplus)
# define EXTERNC extern "C"
#else
# define EXTERNC
#endif /* __cplusplus */

#if defined MSWIN && defined ADDON
# define LINKTOADDON EXTERNC _declspec ( dllimport )
#else
# define LINKTOADDON EXTERNC
#endif /* MSWIN && ADDON */

/* Note: Add-ons should never define AMTEC_INTERNAL_MAKELIBTEC */
#if defined MSWIN && !defined TECPLOTKERNEL && defined AMTEC_INTERNAL_MAKELIBTEC
# undef LINKTOADDON
# define LINKTOADDON EXTERNC _declspec ( dllexport )
#endif

/*
 *
 * Usage:
 * EXPORTFROMADDON void STDCALL InitTecAddOn(void) { ... }
 *
 */

#if defined MSWIN && defined ADDON
# define EXPORTFROMADDON EXTERNC _declspec ( dllexport )
#else
# define EXPORTFROMADDON EXTERNC
#endif /* MSWIN && ADDON */

#define EXPORTFROMDLL EXPORTFROMADDON 

/* CORE SOURCE CODE REMOVED */


/* CORE SOURCE CODE REMOVED */


/****************************************************************
 *                                                              *
 *                       HARD CONSTANTS                         *
 *                                                              *
 ****************************************************************/
/* CORE SOURCE CODE REMOVED */

#define MAXINDEX                (LgIndex_t)2147483646   /* int */
#define MAXZONEMAP               32700                  /* int */
#define LARGEDOUBLE              1.0e+150               /* double */
#define SMALLDOUBLE              1.0e-150               /* double */
#define LARGESTEXPONENT          150                    /* int */
#define SMALLESTEXPONENT         -150                   /* int */
#if defined VMS
#  define LARGESTDOUBLEEXPONENT  307  /* int */
#  define SMALLESTDOUBLEEXPONENT -307 /* int */
#  define LARGESTDOUBLE          1.0e+307 /* double */
#  define LARGEFLOAT             1.0e+37  /* float */
#  define SMALLFLOAT             1.0e-37  /* float */
#else
#  define LARGESTDOUBLEEXPONENT  308
#  define SMALLESTDOUBLEEXPONENT -307
#  define LARGESTDOUBLE          1.0e+308
#  define LARGEFLOAT             3.40282347E+38
#  define SMALLFLOAT             1.17549435E-38
/* Do not remove SMALLSTDOUBLE: needed for ActiveX library */
#  define SMALLSTDOUBLE          1.0e-307 /* double */
#endif
#define LARGELONG                MAXINDEX /* long */
#define LARGESHORT               32766    /* short */
#define ETX                      3        /* int */
#define LN10                     2.30258509299404568402 /* double */
#define PIOVER2                  1.57079632679489661923 /* double */
#define TWOPI                    6.28318530717958647692 /* double */
#if defined PI
#undef PI
#endif
#define PI                       3.14159265358979323846 /* double */
#define ANGLEEPSILON             1.0e-10                /* double */
#define LARGESTANGLE             (4*PI+ANGLEEPSILON)    /* double */
#define DEGPERRADIANS            57.295779513082323     /* double */
#define CMPERINCH                2.54                   /* double */
#define POINTSPERINCH            72.0                   /* double */
#define FONTMOVEMARK             192                    /* int */
#define FONTDECISIONMARK         128                    /* int */
#define FONTLINEMARK             64                     /* int */ 
#define BAD_SET_VALUE            ((SetIndex_t)-1)       /* int */
#define SOLID_TRANSLUCENCY       0                      /* int */
#define BAD_DISTANCE             (-1.0)                 /* double */
/* MIN_CIRCUMFERENTIAL_INDEX is the min J dimension for circular zones */
#define MIN_CIRCUMFERENTIAL_INDEX  4                    /* int */

/* CORE SOURCE CODE REMOVED */

/* CORE SOURCE CODE REMOVED */

#define CurBinaryFileVersion         75  /* int */
#define TecplotInterfaceVersion      90  /* int */            

#define    MaxNumZonesOrVars         MAXZONEMAP /* int */
#define    MaxNumPickObjects         1500 /* int */
#define    MaxXAxes                  5    /* int */
#define    MaxYAxes                  5    /* int */
#define    MaxGeoSegments            50   /* int */
#define    MaxPtsCircleOrEllipse     720  /* int */
#define    MaxFrames                 128  /* int */
#define    MaxCustomLabelSets        10   /* int */
#define    MaxCustomLabelsPerSet     5000 /* int */
#define    MaxFontMoves              20000  /* int */
#define    MaxColorMapOverrides      16     /* int */
#define    MaxValueBlankConstraints  8      /* int */
#define    MaxChrsDatasetTitle       256    /* int */
#define    MaxChrsZnTitle            64     /* int */
#define    MaxChrsVarName            64     /* int */
#define    MaxNumViews               16     /* int */
#define    MaxBasicSizes             5      /* int */
#define    MaxColorMapControlPoints  9      /* int */
#define    MaxRawColorMapEntries     800    /* int */
#define    MaxDataSetReaders         100    /* int */
#define    MaxExtendedCurveFits      100    /* int */
#define    MaxColorMapCycles         20     /* int */


/* Dimension Limits */

#define    MinPaperDimInWorkArea     0.5    /* double */
#define    MinFrameWidth             0.25   /* double */
#define    MinFrameHeight            0.25   /* double */
#define    MinAxisLength             0.1    /* double */


#define    BadEnumValue              255    /* int */


/* CORE SOURCE CODE REMOVED */

/* CORE SOURCE CODE REMOVED */

#define    Black_C           ((ColorIndex_t)0)
#define    Red_C             ((ColorIndex_t)1)
#define    Green_C           ((ColorIndex_t)2)
#define    Blue_C            ((ColorIndex_t)3)
#define    Cyan_C            ((ColorIndex_t)4)
#define    Yellow_C          ((ColorIndex_t)5)
#define    Purple_C          ((ColorIndex_t)6)
#define    White_C           ((ColorIndex_t)7)

#define    Custom1_C         ((ColorIndex_t)8)
#define    Custom2_C         ((ColorIndex_t)9)
#define    Custom3_C         ((ColorIndex_t)10)
#define    Custom4_C         ((ColorIndex_t)11)
#define    Custom5_C         ((ColorIndex_t)12)
#define    Custom6_C         ((ColorIndex_t)13)
#define    Custom7_C         ((ColorIndex_t)14)
#define    Custom8_C         ((ColorIndex_t)15)
#define    Custom9_C         ((ColorIndex_t)16)

#define    Custom10_C         ((ColorIndex_t)17)
#define    Custom11_C         ((ColorIndex_t)18)
#define    Custom12_C         ((ColorIndex_t)19)
#define    Custom13_C         ((ColorIndex_t)20)
#define    Custom14_C         ((ColorIndex_t)21)
#define    Custom15_C         ((ColorIndex_t)22)
#define    Custom16_C         ((ColorIndex_t)23)
#define    Custom17_C         ((ColorIndex_t)24)
#define    Custom18_C         ((ColorIndex_t)25)
#define    Custom19_C         ((ColorIndex_t)26)

#define    Custom20_C         ((ColorIndex_t)27)
#define    Custom21_C         ((ColorIndex_t)28)
#define    Custom22_C         ((ColorIndex_t)29)
#define    Custom23_C         ((ColorIndex_t)30)
#define    Custom24_C         ((ColorIndex_t)31)
#define    Custom25_C         ((ColorIndex_t)32)
#define    Custom26_C         ((ColorIndex_t)33)
#define    Custom27_C         ((ColorIndex_t)34)
#define    Custom28_C         ((ColorIndex_t)35)
#define    Custom29_C         ((ColorIndex_t)36)

#define    Custom30_C         ((ColorIndex_t)37)
#define    Custom31_C         ((ColorIndex_t)38)
#define    Custom32_C         ((ColorIndex_t)39)
#define    Custom33_C         ((ColorIndex_t)40)
#define    Custom34_C         ((ColorIndex_t)41)
#define    Custom35_C         ((ColorIndex_t)42)
#define    Custom36_C         ((ColorIndex_t)43)
#define    Custom37_C         ((ColorIndex_t)44)
#define    Custom38_C         ((ColorIndex_t)45)
#define    Custom39_C         ((ColorIndex_t)46)

#define    Custom40_C         ((ColorIndex_t)47)
#define    Custom41_C         ((ColorIndex_t)48)
#define    Custom42_C         ((ColorIndex_t)49)
#define    Custom43_C         ((ColorIndex_t)50)
#define    Custom44_C         ((ColorIndex_t)51)
#define    Custom45_C         ((ColorIndex_t)52)
#define    Custom46_C         ((ColorIndex_t)53)
#define    Custom47_C         ((ColorIndex_t)54)
#define    Custom48_C         ((ColorIndex_t)55)
#define    Custom49_C         ((ColorIndex_t)56)

#define    Custom50_C         ((ColorIndex_t)57)
#define    Custom51_C         ((ColorIndex_t)58)
#define    Custom52_C         ((ColorIndex_t)59)
#define    Custom53_C         ((ColorIndex_t)60)
#define    Custom54_C         ((ColorIndex_t)61)
#define    Custom55_C         ((ColorIndex_t)62)
#define    Custom56_C         ((ColorIndex_t)63)

#define    MultiColor_C      ((ColorIndex_t)(-1))
#define    NoColor_C         ((ColorIndex_t)(-2))
#define    InvalidColor_C    ((ColorIndex_t)(-255))

/* CORE SOURCE CODE REMOVED */

/****************************************************************
 *                                                              *
 *                          SIMPLE TYPEDEFS                     *
 *                                                              *
 ****************************************************************/



#if defined DECALPHA || defined LINUXALPHA || defined COMPAQALPHA
#define LONGIS64
#endif


#if defined CRAY
typedef    unsigned int    UInt64_t;
#endif

#if defined MSWIN
typedef    unsigned __int64 UInt64_t;
#endif

#if defined LONGIS64
typedef    unsigned int    UInt32_t;
typedef    int             LgInteger_t;
#else
typedef    unsigned long   UInt32_t;
typedef    long int        LgInteger_t;
#endif

typedef    unsigned short  UInt16_t;

typedef    int             LgIndex_t;
typedef    LgIndex_t       NodeMap_t;
typedef    LgIndex_t       ScreenDim_t;

/*
 *  The following type is used for passing
 *  arbituary integers or pointers in parameters
 *  to the setvalue function.
 *  NOTE the following:
 *
 *  Pointers are 8 bytes on the DEC ALPHA (same with long's)
 *  Just about everything is 8 bytes on CRAY.
 */

#if defined CRAY
typedef char*      ArbParam_t;
#else
#  if defined LONGIS64
   typedef long ArbParam_t;
#  else
   typedef LgIndex_t  ArbParam_t;
#  endif
#endif
typedef    ArbParam_t      UniqueID_t;

/* used to hold file offset and size values */
typedef long FileOffset_t;

/*
 *  SmInteger must be at least a short....
 */

typedef    unsigned char    Byte_t;
typedef    short            SmInteger_t;
typedef    SmInteger_t      ColorIndex_t;
typedef    SmInteger_t      EntIndex_t;

typedef    char             Boolean_t;
typedef    char            *ZoneName_t;
typedef    char            *VarName_t;
typedef    char            *LString_t;

typedef    LgIndex_t        HeapLength_t;
typedef    LgIndex_t        SegPtsArray_t[MaxGeoSegments];
typedef    double           BasicSize_t[MaxBasicSizes];
typedef    double          *VarList_t;

typedef    long             SetIndex_t;

typedef    unsigned long    SetData_t;
typedef    SetData_t       *SetData_pt;

/* CORE SOURCE CODE REMOVED */

typedef    char             SymbolChar_t[3];




/****************************************************************
 *                                                              *
 *                     ENUMERATED TYPEDEFS                      *
 *                                                              *
 ****************************************************************/


/* CORE SOURCE CODE REMOVED */

/* Used by some of the image exporters/animators */
typedef struct
  {
    Byte_t  R;
    Byte_t  G;
    Byte_t  B;
  } RGBTriple_s; 

typedef RGBTriple_s RGBPalette_t[256];

typedef enum
  {
    StateChange_VarsAltered,
    StateChange_VarsAdded,
    StateChange_ZonesDeleted,
    StateChange_ZonesAdded,
    StateChange_NodeMapsAltered,
    StateChange_FrameDeleted,
    StateChange_NewTopFrame,
    StateChange_Style,
    StateChange_DataSetReset,
    StateChange_NewLayout,
    StateChange_CompleteReset,
    StateChange_XYMapAssignment,
    StateChange_ContourLevels,
    StateChange_ModalDialogLaunch,
    StateChange_ModalDialogDismiss,
    StateChange_QuitTecplot,
    StateChange_ZoneName,
    StateChange_VarName,
    StateChange_XYMapName,
    StateChange_XYMapAddDeleteOrReorder,
    StateChange_View,
    StateChange_ColorMap,
    StateChange_ContourVar,
    StateChange_Streamtrace,
    StateChange_NewAxisVariables,
    StateChange_MouseModeUpdate,
    StateChange_PickListCleared,
    StateChange_PickListGroupSelect,
    StateChange_PickListSingleSelect,
    StateChange_PickListStyle,
    StateChange_DataSetFileName,
    StateChange_UnsuspendInterface, /* was StateChange_DrawGraphicsOn */
    StateChange_SuspendInterface, /* was StateChange_DrawGraphicsOff */
    StateChange_DataSetLockOn,
    StateChange_DataSetLockOff,
    StateChange_Text,
    StateChange_Geom,
    StateChange_DataSetTitle,
    StateChange_DrawingInterrupted,
    StateChange_PrintPreviewLaunch,
    StateChange_PrintPreviewDismiss,
    END_StateChange_e,
    StateChange_Invalid = BadEnumValue,
    /* Deprecated values */
    StateChange_DrawGraphicsOn = StateChange_UnsuspendInterface,
    StateChange_DrawGraphicsOff = StateChange_SuspendInterface
  } StateChange_e; /*<help> "StateChange_DrawGraphicsOn is deprecated. Use StateChange_UnsuspendInterface\n"*/
                   /*<help> "StateChange_DrawGraphicsOff is deprecated. Use StateChange_SuspendInterface"*/

typedef enum
  {
    StateChangeMode_v75,
    StateChangeMode_v80,
    END_StateChangeMode_e,
    StateChangeMode_Invalid = BadEnumValue
  } StateChangeMode_e;

typedef enum
{
  LayoutPackageObject_Image,
  LayoutPackageObject_Layout,
  LayoutPackageObject_Data,
  END_LayoutPackageObject_e,
  LayoutPackageObject_Invalid = BadEnumValue
} LayoutPackageObject_e;

typedef enum
  {
    VarLoadMode_ByName,
    VarLoadMode_ByPosition,
    END_VarLoadMode_e,
    VarLoadMode_Invalid = BadEnumValue
  } VarLoadMode_e;

typedef enum
  {
    ImageSelection_OnePerFrame,
    ImageSelection_WorkspaceOnly,
    END_ImageSelection_e,
    ImageSelection_Invalid = BadEnumValue
  } ImageSelection_e;

typedef enum
  {
    LibraryType_Unknown,
    LibraryType_V7Standard,
    LibraryType_V7ActiveX,
    END_LibraryType_e,
    LibraryType_Invalid = BadEnumValue
  } LibraryType_e; /* <help> "Add-on types" */


typedef enum
  {
    AssignOp_Equals,
    AssignOp_PlusEquals,
    AssignOp_MinusEquals,
    AssignOp_TimesEquals,
    AssignOp_DivideEquals,
    AssignOp_ConvertFromCm,
    AssignOp_ConvertFromIn,
    AssignOp_ConvertFromPt,
    AssignOp_ConvertFromPix,
    END_AssignOp_e,
    AssignOp_Invalid = BadEnumValue
  } AssignOp_e;

typedef enum
  {
    Dialog_ColorMap,
    Dialog_Equation,
    Dialog_MacroViewer,
    Dialog_PlotAttributes,
    Dialog_QuickEdit,
    Dialog_QuickMacroPanel,
    END_Dialog_e,
    Dialog_Invalid = BadEnumValue
  } Dialog_e; /* <help> "Tecplot dialog types" */


typedef enum
  {
    ProcessXYMode_Draw,
    ProcessXYMode_GetXYMinMax,
    ProcessXYMode_GetSinglePick,
    ProcessXYMode_CheckOnlyForGroupPick,
    ProcessXYMode_GetGroupPick,
    ProcessXYMode_GetFirstValidDataPoint,
    ProcessXYMode_GetNearestPoint,
    ProcessXYMode_GetDependentValue,
    ProcessXYMode_DisplayCurveCoef,
    ProcessXYMode_WriteCurveCoef,
    ProcessXYMode_WriteCurvePoints,
    ProcessXYMode_InsertLabels,
    END_ProcessXYMode_e,
    ProcessXYMode_Invalid = BadEnumValue
  } ProcessXYMode_e;


typedef enum
  {
    StyleBase_Factory,
    StyleBase_Config,
    END_StyleBase_e,
    StyleBase_Invalid = BadEnumValue
  } StyleBase_e;


typedef enum
  {
    ReadDataOption_NewData,
    ReadDataOption_AppendData,
    ReadDataOption_ReplaceData,
    END_ReadDataOption_e,
    ReadDataOption_Invalid = BadEnumValue
  } ReadDataOption_e;

typedef enum
  {
    NodeLabel_Index,
    NodeLabel_VarValue,
    NodeLabel_XAndYVarValue,
    END_NodeLabel_e,
    NodeLabel_Invalid = BadEnumValue
  } NodeLabel_e;


typedef enum
  {
    SubBoundaryEditOption_All,
    SubBoundaryEditOption_Add,
    SubBoundaryEditOption_Remove,
    SubBoundaryEditOption_AddOnly,
    END_SubBoundaryEditOption_e,
    SubBoundaryEditOption_Invalid = BadEnumValue
  } SubBoundaryEditOption_e;


typedef enum
  {
    PointerStyle_Undefined,
    PointerStyle_Normal,
    PointerStyle_Adjuster,
    PointerStyle_AllDirections,
    PointerStyle_Rotate,
    PointerStyle_Zoom,
    PointerStyle_Locate,
    PointerStyle_UpperLeftBracket,
    PointerStyle_UpperRightBracket,
    PointerStyle_LeftBracket,
    PointerStyle_LowerLeftBracket,
    PointerStyle_LowerRightBracket,
    PointerStyle_RightBracket,
    PointerStyle_BottomBracket,
    PointerStyle_TopBracket,
    PointerStyle_UpDown,
    PointerStyle_LeftRight,
    PointerStyle_Waiting,
    END_PointerStyle_e,
    PointerStyle_Invalid = BadEnumValue
  } PointerStyle_e;



typedef enum
  {
    NotEditing,
    ActivelyEditing,
    WasEditing,
    END_SingleEditState_e,
    EditingInvalid = BadEnumValue
  } SingleEditState_e;


typedef enum
  {
    SetValue_Ok,
    SetValue_DuplicateValue,
    SetValue_InvalidCommandOption,
    SetValue_NoAttachedDatasetError,
    SetValue_NoAttachedFrameError,
    SetValue_NotAllowedInConfigError,
    SetValue_ValueRangeError,
    SetValue_ValueSyntaxError,
    SetValue_AssignOpError,
    SetValue_InvalidVarOrZone,
    SetValue_InternalMemoryError,
    SetValue_ContextError1,
    SetValue_ContextError2,
    SetValue_OnlyAllowedInConfigError,
    END_SetValueReturnCode_e,
    SetValue_Invalid = BadEnumValue
  } SetValueReturnCode_e; 


typedef enum
  {
    ObjectAlign_LeftJustify,
    ObjectAlign_RightJustify,
    ObjectAlign_Center,
    ObjectAlign_Top,
    ObjectAlign_Bottom,
    END_ObjectAlign_e,
    ObjectAlign_Invalid = BadEnumValue
  } ObjectAlign_e;


/*
 * For 3D axis labels only.
 */
typedef enum
  {
    LabelAlignment_ByAngle,
    LabelAlignment_AlongAxis,
    LabelAlignment_PerpendicularToAxis,
    END_LabelAlignment_e,
    LabelAlignment_Invalid = BadEnumValue
  } LabelAlignment_e; /* <help> Label alignment for 3D axis labels only" */

typedef enum
  {
    View_Fit,
    View_DataFit,
    View_AxisFit,
    View_Scale,
    View_Center,
    View_Translate,
    View_Zoom,
    View_Last,
    View_Copy,
    View_Paste,
    View_Push,
    END_View_e,
    View_Invalid = BadEnumValue
  } View_e;



typedef enum
  {
    WorkspaceView_FitSelectedFrames,
    WorkspaceView_FitAllFrames,
    WorkspaceView_FitPaper,
    WorkspaceView_Maximize,
    WorkspaceView_LastView,
    WorkspaceView_Zoom,
    WorkspaceView_Translate,
    END_WorkspaceView_e,
    WorkspaceView_Invalid = BadEnumValue
  } WorkspaceView_e;


typedef enum
  {
    Arrowhead_Plain,
    Arrowhead_Filled,
    Arrowhead_Hollow,
    END_ArrowheadStyle_e,
    Arrowhead_Invalid = BadEnumValue
  } ArrowheadStyle_e;


typedef enum
  {
    ArrowheadAttach_None,
    ArrowheadAttach_AtBeginning,
    ArrowheadAttach_AtEnd,
    ArrowheadAttach_AtBothEnds,
    END_ArrowheadAttachment_e,
    ArrowheadAttach_Invalid = BadEnumValue
  } ArrowheadAttachment_e;



typedef enum
  {
    StatusInfo_Hover,
    StatusInfo_Identify,
    StatusInfo_Instruction,
    StatusInfo_Working,
    END_StatusInfo_e,
    StatusInfo_Invalid = BadEnumValue
  } StatusInfo_e;



typedef enum
  {
    Frame_Empty,
    Frame_ThreeD,
    Frame_TwoD,
    Frame_XY,
    Frame_Sketch,
    END_FrameMode_e,
    Frame_Invalid = BadEnumValue
  } FrameMode_e;



typedef enum
  {
    PickObject_None,
    PickObject_Frame,
    PickObject_Axis,
    PickObject_3DOrientationAxis,
    PickObject_Geom,
    PickObject_Text,
    PickObject_ContourLegend,
    PickObject_ContourLabel,
    PickObject_ScatterLegend,
    PickObject_XYLegend,
    PickObject_ReferenceVector,
    PickObject_ReferenceScatterSymbol,
    PickObject_StreamtracePosition,
    PickObject_StreamtraceTermLine,
    PickObject_Paper,
    PickObject_Zone,
    PickObject_XYMapping,
    PickObject_StreamtraceCOB,
    PickObject_SliceCOB,
    PickObject_IsoSurfaceCOB,
    END_PickObjects_e,
    PickObject_Invalid = BadEnumValue
  } PickObjects_e;



typedef enum
  {
    AxisSubObject_GridArea,
    AxisSubObject_AxisLine,
    AxisSubObject_Title,
    END_AxisSubObject_e,
    AxisSubObject_Invalid = BadEnumValue
  } AxisSubObject_e;


typedef enum
  {
    AltMouseButtonMode_Regen,
    AltMouseButtonMode_RevertToSelect,
    END_AltMouseButtonMode_e,
    AltMouseButtonMode_Invalid = BadEnumValue
  } AltMouseButtonMode_e;


typedef enum
  {
    Mouse_NoMode,
    Mouse_Select,
    Mouse_Adjust,
    Mouse_Zoom,
    Mouse_Translate,
    Mouse_Probe,
    Mouse_Text,
    Mouse_GeomPolyline,
    Mouse_GeomSquare,
    Mouse_GeomCircle,
    Mouse_GeomRectangle,
    Mouse_GeomEllipse,
    Mouse_GeomSpline,
    Mouse_CreateFrame,
    Mouse_RotateSpherical,
    Mouse_RotateRollerBall,
    Mouse_RotateTwist,
    Mouse_RotateXAxis,
    Mouse_RotateYAxis,
    Mouse_RotateZAxis,
    Mouse_ContourLabel,
    Mouse_ContourAdd,
    Mouse_ContourDelete,
    Mouse_StreamPoints,
    Mouse_StreamEndLine,
    Mouse_ExtractPoints,
    Mouse_ExtractLine,
    Mouse_CreateRectangularZone,
    Mouse_CreateCircularZone,
    Mouse_Slice,
    Mouse_User1,
    Mouse_User2,
    Mouse_User3,
    Mouse_User4,
    END_MouseButtonMode_e,
    Mouse_Invalid = BadEnumValue
  } MouseButtonMode_e;


typedef enum
  {
    DetailsButtonState_QuickEdit,
    DetailsButtonState_ObjectDetails,
    DetailsButtonState_ToolDetails,
    END_DetailsButtonState_e,
    DetailsButtonState_Invalid = BadEnumValue
  } DetailsButtonState_e;


typedef enum
  {
    Event_ButtonPress,
    Event_ButtonRelease,
    Event_ButtonDoublePress,
    Event_Motion,
    Event_Drag,
    Event_KeyPress,
    END_Event_e,
    Event_Invalid = BadEnumValue
  } Event_e;


typedef enum
  {
    ObjectDrawMode_DrawFirst,
    ObjectDrawMode_Move,
    ObjectDrawMode_Remove,
    ObjectDrawMode_Place,
    END_ObjectDrawMode_e,
    ObjectDrawMode_Invalid = BadEnumValue
  } ObjectDrawMode_e;


typedef enum
  {
    ThreeDViewChangeDrawLevel_Full,
    ThreeDViewChangeDrawLevel_Trace,
    END_ThreeDViewChangeDrawLevel_e,
    ThreeDViewChangeDrawLevel_Invalid = BadEnumValue
  } ThreeDViewChangeDrawLevel_e; /* <help> "ThreeDViewChangeDrawLevel is deprecated. Use PlotApproximateMode.\n"*/

typedef enum
  {
    NonCurrentFrameRedrawLevel_Full,
    NonCurrentFrameRedrawLevel_Trace,
    END_NonCurrentFrameRedrawLevel_e,
    NonCurrentFrameRedrawLevel_Invalid = BadEnumValue
  } NonCurrentFrameRedrawLevel_e; /* <help> "NonCurrentFrameRedrawLevel is deprecated. Use PlotApproximateMode.\n"*/


typedef enum
  {
    RotationMode_XYZAxis,
    RotationMode_Spherical,
    RotationMode_RollerBall,
    END_RotationMode_e,
    RotationMode_Invalid = BadEnumValue
  } RotationMode_e;




typedef enum
  {
    RotateAxis_X,
    RotateAxis_Y,
    RotateAxis_Z,
    RotateAxis_Psi,
    RotateAxis_Theta,
    RotateAxis_Alpha,
    RotateAxis_Twist,
    RotateAxis_VertRollerBall,
    RotateAxis_HorzRollerBall,
    RotateAxis_AboutVector,
/* CORE SOURCE CODE REMOVED */
    END_RotateAxis_e,
    RotateAxis_Invalid = BadEnumValue
  } RotateAxis_e;

typedef enum
  {
    RotateOriginLocation_DefinedOrigin,
    RotateOriginLocation_Viewer,
    END_RotateOriginLocation_e,
    RotateOriginLocation_Invalid = BadEnumValue
  } RotateOriginLocation_e;

/*
 * NOTE: This is only used with the $!Reset3DOrigin command.
 */
typedef enum
  {
    OriginResetLocation_DataCenter,
    OriginResetLocation_ViewCenter,
    END_OriginResetLocation_e,
    OriginResetLocation_Invalid = BadEnumValue
  } OriginResetLocation_e; 

/*
 * NOTE: This is only used with the $!CreateSliceZoneFromPlane command.
 */
typedef enum
  {
    SliceSource_SurfaceZones,
    SliceSource_VolumeZones,
    SliceSource_SurfacesOfVolumeZones,
    END_SliceSource_e,
    SliceSource_Invalid = BadEnumValue
  } SliceSource_e; 





typedef enum
  {
    Input_SmInteger,
    Input_Short,
    Input_Integer,
    Input_Float,
    Input_Double,
    END_Input_e,
    Input_Invalid = BadEnumValue
  } Input_e;



typedef enum
  {
    PtSelection_All,
    PtSelection_NearestN,
    PtSelection_OctantN,
    END_PtSelection_e,
    PtSelection_Invalid = BadEnumValue
  } PtSelection_e;



typedef enum
  {
    Drift_None,
    Drift_Linear,
    Drift_Quad,
    END_Drift_e,
    Drift_Invalid = BadEnumValue
  } Drift_e;



/* atpoint is simple boundary condition.
   atpointb2 is better boundary condition.
*/
typedef enum
  {
    DerivPos_atpoint,
    DerivPos_atpointb2,
    DerivPos_kphalf,
    DerivPos_jphalf,
    DerivPos_iphalf,
    END_DerivPos_e,
    DerivPos_Invalid = BadEnumValue
  } DerivPos_e; /*<help>"atpoint is the simple boundary condition\n"*/
                /*<help>"atpointb2 is a better boundary condition"*/



typedef enum
  {
    LinearInterpMode_DontChange,
    LinearInterpMode_SetToConst,
    END_LinearInterpMode_e,
    LinearInterpMode_Invalid = BadEnumValue
  } LinearInterpMode_e;



typedef enum
  {
    ConstraintOp2Mode_UseVar,
    ConstraintOp2Mode_UseConstant,
    END_ConstraintOp2Mode_e,
    ConstraintOp2Mode_Invalid = BadEnumValue
  } ConstraintOp2Mode_e;



typedef enum
  {
    ValueBlankCellMode_AllCorners,
    ValueBlankCellMode_AnyCorner,
    ValueBlankCellMode_PrimaryCorner,
    END_ValueBlankCellMode_e,
    ValueBlankCellMode_Invalid = BadEnumValue
  } ValueBlankCellMode_e;


/*
 * DEPRECATED: ValueBlankMode_e enumeration will not be supported after
 *             version 8. This API was retained for add-on developers
 *             using the TecUtilStyleSetLowLevel API.
 */
typedef enum
  {
    ValueBlankMode_AndRule,
    ValueBlankMode_OrRule,
    ValueBlankMode_CornerRule,
    END_ValueBlankMode_e,
    ValueBlankMode_Invalid = BadEnumValue
  } ValueBlankMode_e; /*<help>"DEPRECATED: ValueBlankMode_e will not be supported after version 8"*/


typedef enum
  {
    CellBlankedCond_NotBlanked,
    CellBlankedCond_PartiallyBlanked,
    CellBlankedCond_EntirelyBlanked,
    CellBlankedCond_Uncertain,
    END_CellBlankedCond_e,
    CellBlankedCond_Invalid = BadEnumValue
  } CellBlankedCond_e;


typedef enum
  {
    RelOp_LessThanOrEqual,
    RelOp_GreaterThanOrEqual,
    END_RelOp_e,
    RelOp_Invalid = BadEnumValue
  } RelOp_e;



typedef enum
  {
    IJKBlankMode_BlankInterior,
    IJKBlankMode_BlankExterior,
    END_IJKBlankMode_e,
    IJKBlankMode_Invalid = BadEnumValue
  } IJKBlankMode_e;


typedef enum
  {
    PlotApproximationMode_Default,
    PlotApproximationMode_NonCurrentAlwaysReduced,
    PlotApproximationMode_AllFramesAlwaysReduced,
    END_PlotApproximationMode_e,
    PlotApproximationMode_Invalid = BadEnumValue
  } PlotApproximationMode_e;


/*
 * NOTE: FillPat_e is deprecated.  It must be retained to maintain
 *       backward compatibility with the TecUtil layer however.
 *       This has been replaced by Translucency_e.
 */
typedef enum                     
  {                             
    Pattern_Solid,             
    Pattern_LowTranslucent,   
    Pattern_MedTranslucent,  
    Pattern_HighTranslucent, 
    END_FillPat_e,           
    Pattern_Invalid = BadEnumValue 
  } FillPat_e; /*<help>"DEPRECATED: Replaced by Translucency_e"*/             


typedef enum                     
  {                             
    Translucency_Solid,             
    Translucency_Low,   
    Translucency_Medium,  
    Translucency_High, 
    END_Translucency_e,           
    Translucency_Invalid = BadEnumValue 
  } Translucency_e;                    



typedef enum
  {
    SunRaster_OldFormat,
    SunRaster_Standard,
    SunRaster_ByteEncoded,
    END_SunRaster_e,
    SunRaster_Invalid = BadEnumValue
  } SunRaster_e;


typedef enum
  {
    BoundaryCondition_Fixed,
    BoundaryCondition_ZeroGradient,
    BoundaryCondition_Zero2nd,
    END_BoundaryCondition_e,
    BoundaryCondition_Invalid = BadEnumValue
  } BoundaryCondition_e;



/* Note:
 *   In 2D: AxisMode_Independent and AxisMode_XYZDependent are used;
 *   in 3D: AxisMode_Independent, AxisMode_XYZDependent, and AxisMode_XYDependent are used.
 */
typedef enum
  {
    AxisMode_Independent,
    AxisMode_XYZDependent,
    AxisMode_XYDependent,
    END_AxisMode_e,
    AxisMode_Invalid = BadEnumValue
  } AxisMode_e;/*<help>"In 2D AxisMode_Independent and AxisMode_XYZDependent are used\n"*/
               /*<help>"In 3D AxisMode_Independent, "*/
               /*<help>"AxisMode_XYZDependent, and AxisMode_XYDependent are used."*/

typedef enum
  {
    Quick_LineColor,
    Quick_FillColor,
    Quick_TextColor,
    END_QuickColorMode_e,
    Quick_Invalid = BadEnumValue
  } QuickColorMode_e;




typedef enum
  {
    LinePattern_Solid,
    LinePattern_Dashed,
    LinePattern_DashDot,
    LinePattern_Dotted,
    LinePattern_LongDash,
    LinePattern_DashDotDot,
    END_LinePattern_e,
    LinePattern_Invalid = BadEnumValue
  } LinePattern_e;



typedef enum
  {
    Join_Miter,
    Join_Round,
    Join_Bevel,
    END_LineJoin_e,
    Join_Invalid = BadEnumValue
  } LineJoin_e;



typedef enum
  {
    Cap_Flat,
    Cap_Round,
    Cap_Square,
    END_LineCap_e,
    Cap_Invalid = BadEnumValue
  } LineCap_e;



typedef enum
  {
    GeomForm_LineSegs,
    GeomForm_Rectangle,
    GeomForm_Square,      /* New!  Make sure and convert old style sheets!!! */
    GeomForm_Circle,
    GeomForm_Ellipse,
    GeomForm_LineSegs3D,
    END_GeomForm_e,
    GeomForm_Invalid = BadEnumValue
  } GeomForm_e;



/* Note: This replaces Element_e */
typedef enum
  {
    ZoneType_Ordered,
    ZoneType_FETriangle,
    ZoneType_FEQuad,
    ZoneType_FETetra,
    ZoneType_FEBrick,
    END_ZoneType_e,
    ZoneType_Invalid = BadEnumValue
  } ZoneType_e;

typedef enum
  {
    ZoneOrder_I,
    ZoneOrder_J,
    ZoneOrder_K,
    ZoneOrder_IJ,
    ZoneOrder_IK,
    ZoneOrder_JK,
    ZoneOrder_IJK,
    END_ZoneOrder_e,
    ZoneOrder_Invalid = BadEnumValue
  } ZoneOrder_e;

typedef enum
  {
    DataFormat_IJKBlock,
    DataFormat_IJKPoint,
    DataFormat_FEBlock,
    DataFormat_FEPoint,
    END_DataFormat_e,
    DataFormat_Invalid = BadEnumValue
  } DataFormat_e;



typedef enum
  {
    PD_HPGL,
    PD_HPGL2,
    PD_PS,
    PD_LASERG, /* deprecated */
    PD_EPS,
    PD_WINDOWS, /* Windows Print Driver */
    PD_WMF, /* Windows MetaFile (used from Export only) */
    END_PrinterDriver_e,
    PD_Invalid = BadEnumValue
  } PrinterDriver_e;



typedef enum
  {
    Image_None,
    Image_TIFF,
    Image_EPSI2,
    Image_FRAME,
    END_EPSPreviewImage_e,
    Image_Invalid = BadEnumValue
  } EPSPreviewImage_e;

typedef enum
  {
    TIFFByteOrder_Intel,
    TIFFByteOrder_Motorola,
    END_TIFFByteOrder_e,
    TIFFByteOrder_Invalid = BadEnumValue
  } TIFFByteOrder_e;

typedef enum
  {
    ExportFormat_RasterMetafile,
    ExportFormat_TIFF,
    ExportFormat_SGI,
    ExportFormat_SunRaster,
    ExportFormat_XWindows,
    ExportFormat_PSImage,
    ExportFormat_HPGL,
    ExportFormat_HPGL2,
    ExportFormat_PS,
    ExportFormat_EPS,
    ExportFormat_LaserGraphics, /* deprecated */
    ExportFormat_WindowsMetafile,
    ExportFormat_BMP,
    ExportFormat_PNG,
    ExportFormat_AVI,
    ExportFormat_Custom,  /* May be used in a future version */
    END_ExportFormat_e,
    ExportFormat_Invalid = BadEnumValue
  } ExportFormat_e;

typedef enum
  {
    AnimationDest_Screen,
    AnimationDest_AVI,
    AnimationDest_RM,
    END_AnimationDest_e,
    AnimationDest_Invalid = BadEnumValue
  } AnimationDest_e;



/* CORE SOURCE CODE REMOVED */

typedef enum
  {
    BitDumpRegion_CurrentFrame,
    BitDumpRegion_AllFrames,
    BitDumpRegion_WorkArea,
    END_BitDumpRegion_e,
    BitDumpRegion_Invalid = BadEnumValue
  } BitDumpRegion_e;


typedef enum
  {
    Paper_Letter,
    Paper_Double,
    Paper_A4,
    Paper_A3,
    Paper_Custom1,
    Paper_Custom2,
    END_PaperSize_e,
    Paper_Invalid = BadEnumValue
  } PaperSize_e;



typedef enum
  {
    PaperUnitSpacing_HalfCentimeter,
    PaperUnitSpacing_OneCentimeter,
    PaperUnitSpacing_TwoCentimeters,
    PaperUnitSpacing_QuarterInch,
    PaperUnitSpacing_HalfInch,
    PaperUnitSpacing_OneInch,
    PaperUnitSpacing_TenPoints,
    PaperUnitSpacing_TwentyFourPoints,
    PaperUnitSpacing_ThirtySixPoints,
    PaperUnitSpacing_FiftyPoints,
    PaperUnitSpacing_SeventyTwoPoints,
    PaperUnitSpacing_OneTenthInch,
    PaperUnitSpacing_OneTenthCentimeter,
    END_PaperUnitSpacing_e,
    PaperUnitSpacing_Invalid = BadEnumValue
  } PaperUnitSpacing_e;


typedef enum
  {
    Palette_Monochrome,
    Palette_PenPlotter,
    Palette_Color,
    END_Palette_e,
    Palette_Invalid = BadEnumValue
  } Palette_e;


typedef enum
  {
    PrintRenderType_Vector,
    PrintRenderType_Image,
    END_PrintRenderType_e,
    PrintRenderType_Invalid = BadEnumValue
  } PrintRenderType_e;


typedef enum
  {
    Units_Grid,
    Units_Frame,
    Units_Point,
    Units_Screen,
    END_Units_e,
    Units_Invalid = BadEnumValue
  } Units_e;


typedef enum
  {
    Scale_Linear,
    Scale_Log,
    END_CoordScale_e,
    Scale_Invalid = BadEnumValue
  } CoordScale_e;





typedef enum
  {
    CoordSys_Grid,
    CoordSys_Frame,
    CoordSys_FrameOffset,
    CoordSys_Paper,
    CoordSys_Screen,
    CoordSys_Hardcopy,
    END_CoordSys_e,
    CoordSys_Invalid = BadEnumValue
  } CoordSys_e;

/*
 *  NOTE:  CoordSys_FrameOffset always is stored in inches internally.
 *         in stylesheet this may be written in other units if
 *         appropriate suffix is added.
 *
 */



typedef enum
  {
    Scope_Global,
    Scope_Local,
    END_Scope_e,
    Scope_Invalid = BadEnumValue
  } Scope_e;



typedef enum
  {
    TextAnchor_Left,
    TextAnchor_Center,
    TextAnchor_Right,
    TextAnchor_MidLeft,
    TextAnchor_MidCenter,
    TextAnchor_MidRight,
    TextAnchor_HeadLeft,
    TextAnchor_HeadCenter,
    TextAnchor_HeadRight,
    TextAnchor_OnSide,
    END_TextAnchor_e,
    TextAnchor_Invalid = BadEnumValue
  } TextAnchor_e;



typedef enum
  {
    TextBox_None,
    TextBox_Filled,
    TextBox_Hollow,
    END_TextBox_e,
    TextBox_Invalid = BadEnumValue
  } TextBox_e;



typedef enum
  {
    GeomShape_Square,
    GeomShape_Del,
    GeomShape_Grad,
    GeomShape_RTri,
    GeomShape_LTri,
    GeomShape_Diamond,
    GeomShape_Circle,
    GeomShape_Cube,
    GeomShape_Sphere,
    GeomShape_Tetra,
    GeomShape_Pyramid,
    END_GeomShape_e,
    GeomShape_Invalid = BadEnumValue
  } GeomShape_e;


typedef enum
  {
    BasicSize_Tiny,
    BasicSize_Small,
    BasicSize_Medium,
    BasicSize_Large,
    BasicSize_Huge,
    END_BasicSize_e,
    BasicSize_Invalid = BadEnumValue
  } BasicSize_e;



/*
 * NOTE: LineForm_e is deprecated.  It must be retained to maintain
 *       backward compatibility with the TecUtil layer however.
 *       This has been replaced by CurveType_e.
 */
typedef enum
  {
    LineForm_LineSeg,
    LineForm_CurvFit,
    LineForm_EToRFit,
    LineForm_PowerFit,
    LineForm_Spline,
    LineForm_ParaSpline,
    END_LineForm_e,
    LineForm_Invalid = BadEnumValue
  } LineForm_e;


typedef enum
  {
    CurveType_LineSeg,
    CurveType_CurvFit,
    CurveType_EToRFit,
    CurveType_PowerFit,
    CurveType_Spline,
    CurveType_ParaSpline,
    CurveType_Extended,
    END_CurveType_e,
    CurveType_Invalid = BadEnumValue
  } CurveType_e;

typedef enum
  {
    Script_None,
    Script_Super,
    Script_Sub,
    END_Script_e,
    Script_Invalid = BadEnumValue
  } Script_e;


typedef enum
  {
    Font_Helvetica,
    Font_HelveticaBold,
    Font_Greek,
    Font_Math,
    Font_UserDefined,
    Font_Times,
    Font_TimesItalic,
    Font_TimesBold,
    Font_TimesItalicBold,
    Font_Courier,
    Font_CourierBold,
    END_Font_e,
    Font_Invalid = BadEnumValue
  } Font_e;

typedef enum
  {
    TwoDDrawOrder_ByZone,
    TwoDDrawOrder_ByLayer,
    END_TwoDDrawOrder_e,
    TwoDDrawOrder_Invalid = BadEnumValue
  } TwoDDrawOrder_e;

/*
 *
 * NOTE: Streamtrace_TwoDLine is new.  All 2D
 *       streamtraces are assigned this value.
 */
typedef enum
  {
    Streamtrace_SurfaceLine,
    Streamtrace_SurfaceRibbon,
    Streamtrace_VolumeLine,
    Streamtrace_VolumeRibbon,
    Streamtrace_VolumeRod,
    Streamtrace_TwoDLine,
    END_Streamtrace_e,
    Streamtrace_Invalid = BadEnumValue
  } Streamtrace_e;



typedef enum
  {
    StreamDir_Forward,
    StreamDir_Reverse,
    StreamDir_Both,
    END_StreamDir_e,
    StreamDir_Invalid = BadEnumValue
  } StreamDir_e;

typedef enum
  {
    IsoSurfaceSelection_AllContourLevels,
    IsoSurfaceSelection_OneSpecificValue,
    IsoSurfaceSelection_TwoSpecificValues,
    IsoSurfaceSelection_ThreeSpecificValues,
    END_IsoSurfaceSelection_e,
    IsoSurfaceSelection_Invalid = BadEnumValue
  } IsoSurfaceSelection_e;


typedef enum
  {
    FieldDataType_Reserved, /* never use */
    FieldDataType_Float,
    FieldDataType_Double,
    FieldDataType_LongInt,
    FieldDataType_ShortInt,
    FieldDataType_Byte,
    FieldDataType_Bit,
    FieldDataType_IJKFunction,   /* Not used yet */
    END_FieldDataType_e,
    FieldDataType_Invalid = BadEnumValue
  } FieldDataType_e;


typedef enum
  {
    Mesh_Wireframe,
    Mesh_Overlay,
    Mesh_HiddenLine,
    END_MeshPlotType_e,
    Mesh_Invalid = BadEnumValue
  } MeshPlotType_e;




typedef enum
  {
    Contour_Lines,
    Contour_Flood,
    Contour_Overlay,
    Contour_AverageCell,
    Contour_CornerCell,
    END_ContourPlotType_e,
    Contour_Invalid = BadEnumValue
  } ContourPlotType_e;




typedef enum
  {
    Vector_TailAtPoint,
    Vector_HeadAtPoint,
    Vector_MidAtPoint,
    Vector_HeadOnly,
    END_VectorPlotType_e,
    Vector_Invalid = BadEnumValue
  } VectorPlotType_e;



/*
 * NOTE: ShadePlotType_e is deprecated.  It must be retained to maintain
 *       backward compatibility with the TecUtil layer however.
 *       This has been replaced by LightingEffect_e.
 */
typedef enum
  {
    Shade_SolidColor,
    Shade_Paneled,
    Shade_Gouraud,
    Shade_ColoredPaneled,
    Shade_ColoredGouraud,
    END_ShadePlotType_e,
    Shade_Invalid = BadEnumValue
  } ShadePlotType_e;

/*
 * NOTE: LightingEffect_None is Deprecated.  It must remain
 *       in the list to allow macro processing of older 
 *       (i.e. early v9) macros.
 */
typedef enum
  {
    LightingEffect_Paneled,
    LightingEffect_Gouraud,
    LightingEffect_None,
    END_LightingEffect_e,
    LightingEffect_Invalid = BadEnumValue
  } LightingEffect_e;

typedef enum
  {
    Lines_I,
    Lines_J,
    Lines_K,
    END_IJKLines_e,
    Lines_Invalid = BadEnumValue
  } IJKLines_e;

typedef enum
  {
    IJKCellType_Planes,
    IJKCellType_FacePlanes,
    IJKCellType_Volume,
    END_IJKCellType_e,
    IJKCellType_Invalid = BadEnumValue
  } IJKCellType_e;


/*
 *  Ver 6 used PlaneSet.  Ver 7 uses CellType and Planes variables.
 *
 *   "PlaneSet" in version 6    vs.  IJKPlanes in v7:
 *
 *   'A' = AllPlanes                 CellType = IJKCellType_Volume
 *   'd','e','f','C' = ComboPlanes   CellType = IJKCellType_Planes, IJKPlanes = depends on defC
 *   'F' = Faces Planes Only         CellType = IJKCellType_FacePlanes
 *   'I' = I-Planes                  CellType = IJKCellType_Planes, IJKPlanes = Planes_I
 *   'J' = J-Planes                  CellType = IJKCellType_Planes, IJKPlanes = Planes_J
 *   'K' = K-Planes                  CellType = IJKCellType_Planes, IJKPlanes = Planes_K
 *
 *
 * NOTE: IJKPlanes_e is still used internally in tecplot (and in the TecUtil layer).
 *       it has been relagated to communicating which planes of an IJK zone are in
 *       use.  
 *
 */

typedef enum
  {
    Planes_I,
    Planes_J,
    Planes_K,
    Planes_IJ,     /* deprecated */
    Planes_JK,     /* deprecated */
    Planes_IK,     /* deprecated */
    Planes_IJK,    /* deprecated */
    Planes_Face,   /* not used anymore - retain to maintain enum list spacing */
    Planes_Volume, /* only necessary for probing/blanking */
    Planes_Unused,
    END_IJKPlanes_e,
    Planes_Invalid = BadEnumValue
  } IJKPlanes_e;



typedef enum
  {
    SurfacesToPlot_BoundaryFaces,
    SurfacesToPlot_ExposedCellFaces,
    SurfacesToPlot_IPlanes,
    SurfacesToPlot_JPlanes,
    SurfacesToPlot_KPlanes,
    SurfacesToPlot_IJPlanes,
    SurfacesToPlot_JKPlanes,
    SurfacesToPlot_IKPlanes,
    SurfacesToPlot_IJKPlanes,
    SurfacesToPlot_All,
    END_SurfacesToPlot_e,
    SurfacesToPlot_Invalid = BadEnumValue
  } SurfacesToPlot_e;

typedef enum
  {
    PointsToPlot_SurfacesOnly,
    PointsToPlot_All,
    END_PointsToPlot_e,
    PointsToPlot_Invalid = BadEnumValue
  } PointsToPlot_e;


typedef enum
{
  SliceSurface_XPlanes,
  SliceSurface_YPlanes,
  SliceSurface_ZPlanes,
  SliceSurface_IPlanes,
  SliceSurface_JPlanes,
  SliceSurface_KPlanes,
  END_SliceSurface_e,
  SliceSurface_Invalid = BadEnumValue
} SliceSurface_e;


typedef enum
  {
    Skip_ByIndex,
    Skip_ByFrameUnits,
    END_SkipMode_e,
    Skip_Invalid = BadEnumValue
  } SkipMode_e;




typedef enum
  {
    Boundary_None,
    Boundary_Min,
    Boundary_Max,
    Boundary_Both,
    END_BoundPlotType_e,
    Boundary_Invalid = BadEnumValue
  } BoundPlotType_e;




typedef enum
  {
    ColorMap_SmRainbow,
    ColorMap_LgRainbow,
    ColorMap_Modern,
    ColorMap_GrayScale,
    ColorMap_Wild,
    ColorMap_UserDef,
    ColorMap_TwoColor,
    ColorMap_RawUserDef,
    END_ContourColorMap_e,
    ColorMap_Invalid = BadEnumValue
  } ContourColorMap_e;



typedef enum
  {
    ErrorBar_Up,
    ErrorBar_Down,
    ErrorBar_Left,
    ErrorBar_Right,
    ErrorBar_Horz,
    ErrorBar_Vert,
    ErrorBar_Cross,
    END_ErrorBar_e,
    ErrorBar_Invalid = BadEnumValue
  } ErrorBar_e;



typedef enum
  {
    ContourLineMode_UseZoneLineType,
    ContourLineMode_SkipToSolid,
    ContourLineMode_DashNegative,
    END_ContourLineMode_e,
    ContourLineMode_Invalid = BadEnumValue
  } ContourLineMode_e;



typedef enum
  {
    Panel_Bad,
    Panel_Cell,                     /* FieldZone */
    Panel_Vector,                   /* FieldZone */
    Panel_Scatter,                  /* FieldZone */
    Panel_BoundaryLine,             /* FieldZone */
    Panel_FEBoundaryLine,           /* FieldZone */
    Panel_FEBoundaryCell,           /* FieldZone */
    Panel_NodeLabel,                /* FieldZone */
    Panel_StreamtraceCell,          /* Streamtrace COB          */
    Panel_StreamtraceMarker,        /* StreamtraceMarker COB (Scatter Symbol) */
    Panel_StreamtraceArrowhead,     /* StreamtraceArrowhead COB (Vector) */
    Panel_IsoSurfaceCell,           /* IsoSurface COB */
    Panel_IsoSurfaceFEBoundaryLine, /* IsoSurface COB */
    Panel_SliceCell,                /* Slice COB */
    Panel_SliceVector,              /* Slice COB */
    Panel_SliceBoundaryLine,        /* Slice COB */
    Panel_SliceFEBoundaryLine,      /* Slice COB */
    Panel_Geom,                     /* Misc */
    END_Panel_e,
    Panel_Invalid = BadEnumValue
  } Panel_e;

typedef enum
  {
    MessageBox_Error,
    MessageBox_Warning,
    MessageBox_Information,
    MessageBox_Question,   /* Ok, Cancel buttons */
    MessageBox_YesNo,
    MessageBox_YesNoCancel,
    END_MessageBoxType_e,
    MessageBox_Invalid = BadEnumValue
  } MessageBoxType_e;


typedef enum
  {
    MessageBoxReply_Yes,
    MessageBoxReply_No,
    MessageBoxReply_Cancel,
    MessageBoxReply_Ok,
    END_MessageBoxReply_e,
    MessageBoxReply_Invalid = BadEnumValue
  } MessageBoxReply_e;

typedef enum
  {
    NumberFormat_Integer,
    NumberFormat_FixedFloat,
    NumberFormat_Exponential,
    NumberFormat_BestFloat,
    NumberFormat_SuperScript,
    NumberFormat_CustomLabel,
    NumberFormat_LogSuperScript,
    END_NumberFormat_e,
    NumberFormat_Invalid = BadEnumValue
  } ValueFormat_e;


typedef enum
  {
    BackingStoreMode_QuickAndDirty,
    BackingStoreMode_RealTimeUpdate,
    BackingStoreMode_PeriodicUpdate,
    END_BackingStoreMode_e,
    BackingStoreMode_Invalid = BadEnumValue
  } BackingStoreMode_e;


typedef enum
  {
    TickDirection_In,
    TickDirection_Out,
    TickDirection_Centered,
    END_TickDirection_e,
    TickDirection_Invalid = BadEnumValue
  } TickDirection_e;

typedef enum
  {
    AxisTitlePosition_Left,
    AxisTitlePosition_Center,
    AxisTitlePosition_Right,
    END_AxisTitlePosition_e,
    AxisTitlePosition_Invalid = BadEnumValue
  } AxisTitlePosition_e;

typedef enum
  {
    AxisTitleMode_NoTitle,
    AxisTitleMode_UseVarName,
    AxisTitleMode_UseText,
    END_AxisTitleMode_e,
    AxisTitleMode_Invalid = BadEnumValue
  } AxisTitleMode_e;

typedef enum
  {
    FunctionDependency_XIndependent,
    FunctionDependency_YIndependent,
    END_FunctionDependency_e,
    FunctionDependency_Invalid = BadEnumValue
  } FunctionDependency_e;

typedef enum
  {
    LaunchDialogMode_ModalSync,
    LaunchDialogMode_Modeless,
    LaunchDialogMode_ModalAsync,
    END_LaunchDialogMode_e,
    LaunchDialogMode_Invalid = BadEnumValue
  } LaunchDialogMode_e;


typedef enum
  {
    SelectFileOption_ReadSingleFile,
    SelectFileOption_ReadMultiFile,
    SelectFileOption_AllowMultiFileRead,
    SelectFileOption_WriteFile,
    END_SelectFileOption_e,
    SelectFileOption_Invalid = BadEnumValue
  } SelectFileOption_e;


/*   CURRENTLY NOT USED .... */
typedef enum
  {
    ViewActionDrawMode_NoDraw,
    ViewActionDrawMode_DrawTrace,
    ViewActionDrawMode_DrawFull,
    END_ViewActionDrawMode_e,
    ViewActionDrawMode_Invalid = BadEnumValue
  } ViewActionDrawMode_e;


typedef enum
  {
    FrameAction_PushTop,
    FrameAction_Pop,
    FrameAction_PopAtPosition,
    FrameAction_DeleteTop,
    FrameAction_FitAllToPaper,
    FrameAction_PushByName,
    FrameAction_PopByName,
    FrameAction_Push,
    END_FrameAction_e,
    FrameAction_Invalid = BadEnumValue
  } FrameAction_e;

typedef enum
  {
    DoubleBufferAction_On,
    DoubleBufferAction_Off,
    DoubleBufferAction_Swap,
    END_DoubleBufferAction_e,
    DoubleBufferAction_Invalid = BadEnumValue
  } DoubleBufferAction_e;

typedef enum
  {
    PickAction_CheckToAdd,
    PickAction_AddAll,
    PickAction_AddAllInRegion,
    PickAction_Edit,
    PickAction_Cut,
    PickAction_Copy,
    PickAction_Clear,
    PickAction_Paste,
    PickAction_PasteAtPosition,
    PickAction_Shift,
    PickAction_Magnify,
    PickAction_Push,
    PickAction_Pop,
    PickAction_SetMouseMode,
    PickAction_DeselectAll,
    PickAction_AddZones,
    PickAction_AddXYMaps,
    END_PickAction_e,
    PickAction_Invalid = BadEnumValue
  } PickAction_e;


typedef enum
  {
    ContourLevelAction_Add,
    ContourLevelAction_New,
    ContourLevelAction_DeleteRange,
    ContourLevelAction_Reset,
    ContourLevelAction_DeleteNearest,
    END_ContourLevelAction_e,
    ContourLevelAction_Invalid = BadEnumValue
  } ContourLevelAction_e;

typedef enum
  {
    ContourLabelAction_Add,
    ContourLabelAction_DeleteAll,
    END_ContourLabelAction_e,
    ContourLabelAction_Invalid = BadEnumValue
  } ContourLabelAction_e;

typedef enum
  {
    StreamtraceAction_Add,
    StreamtraceAction_DeleteAll,
    StreamtraceAction_DeleteRange,
    StreamtraceAction_SetTerminationLine,
    StreamtraceAction_ResetDeltaTime,
    END_StreamtraceAction_e,
    StreamtraceAction_Invalid = BadEnumValue
  } StreamtraceAction_e;

typedef enum
  {
    ColorMapControlAction_RedistributeControlPoints,
    ColorMapControlAction_CopyCannedColorMap,
    ColorMapControlAction_ResetToFactoryDefaults,
    END_ColorMapControlAction_e,
    ColorMapControlAction_Invalid = BadEnumValue
  } ColorMapControlAction_e;

typedef enum
  {
    ColorMapDistribution_Continuous,
    ColorMapDistribution_Banded,
    END_ColorMapDistribution_e,
    ColorMapDistribution_Invalid = BadEnumValue
  } ColorMapDistribution_e;

typedef enum
  {
    TecUtilErr_None,
    TecUtilErr_Undetermined,
    END_TecUtilErr_e,
    TecUtilErr_Invalid = BadEnumValue
  } TecUtilErr_e;

/* CORE SOURCE CODE REMOVED */

typedef enum /* Custom exporter error message */
{
  ExportCustReturnCode_Ok,
  ExportCustReturnCode_Failed,
  ExportCustReturnCode_TecplotLocked,
  ExportCustReturnCode_ExporterNotLoaded,     
  ExportCustReturnCode_ExportCallbackFailed,          
  ExportCustReturnCode_NotAnImageExporter,    
  ExportCustReturnCode_NotAFieldDataExporter, 
  END_ExportCustReturnCode_e,
  ExportCustReturnCode_Invalid = BadEnumValue
} ExportCustReturnCode_e;

/* CORE SOURCE CODE REMOVED */

/**
 * COB/Zone types.
 */
typedef enum
  {
    CZType_FieldDataZone,
    CZType_FEBoundaryCOB,
    CZType_IsoSurfaceCOB,
    CZType_SliceCOB,
    CZType_StreamtraceCOB,
    CZType_StreamtraceMarkerCOB,
    CZType_StreamtraceArrowheadCOB,
    END_CZType_e,
    CZType_Invalid = BadEnumValue
  } CZType_e;


/* CORE SOURCE CODE REMOVED */




/****************************************************************
 *                                                              *
 *                     STRUCTURE TYPEDEFS                       *
 *                                                              *
 ****************************************************************/

typedef struct _StringList_s *StringList_pa;
/* CORE SOURCE CODE REMOVED */

/* CORE SOURCE CODE REMOVED */

typedef struct _Set_a *Set_pa;


typedef struct _XYZ_s
  {
    double X;
    double Y;
    double Z;
  } XYZ_s;


/* CORE SOURCE CODE REMOVED */

/* CORE SOURCE CODE REMOVED */

typedef struct _NodeMap_a *NodeMap_pa;

/* CORE SOURCE CODE REMOVED */

typedef struct _PointState_a   *PointState_pa;
typedef struct _ElementState_a *ElementState_pa;
typedef struct _CoarseBoundary_a *CoarseBoundary_pa;

/* used to indicate that no neighboring element exists */
#define NO_NEIGHBORING_ELEMENT (-1)

typedef struct _FaceNeighbor_a *FaceNeighbor_pa;


/* CORE SOURCE CODE REMOVED */

typedef struct _FieldData_a *FieldData_pa;

/* CORE SOURCE CODE REMOVED */


typedef void (STDCALL *ProbeDestination_pf)(Boolean_t IsNearestPoint);

/*
 *  DynamicMenu Functions are called upon a user selecting
 *  a menu item added via AddMenuSubOption.
 */
typedef void (STDCALL *DynamicMenuCallback_pf)(void);

/*
 *  ExtractDestination Functions are called upon successful
 *  completion of an extract polyline or extract discrete points
 *  operation.
 */
typedef void (STDCALL *ExtractDestination_pf)(LgIndex_t NumPts,
                                              double   *XValues,
                                              double   *YValues);

/*
 *  SelectFileOptionsCallback Functions are called when the
 *  "Options" button is pressed in the modal file selection
 *  dialog.
 */
typedef void (STDCALL *SelectFileOptionsCallback_pf)(void);


typedef Boolean_t (STDCALL *DataSetConverter_pf)(char  *DataFName,
                                                 char  *TempBinFName,
                                                 char **MessageString);

typedef Boolean_t (STDCALL *DataSetLoader_pf)(StringList_pa Instructions);

typedef Boolean_t (STDCALL *DataSetLoaderInstructionOverride_pf)(StringList_pa  Instructions);

/*
 *  Extended Curve Fit Functions
 */
typedef void (STDCALL *GetCurveSettingsCallback_pf) (Set_pa        XYMapSet,
                                                     StringList_pa SelectedXYMapSettings);
typedef void (STDCALL *GetAbbreviatedSettingsStringCallback_pf) (EntIndex_t  XYMapNum,
                                                                 char       *CurveSettings,
                                                                 char      **AbbreviatedSettings);
typedef Boolean_t (STDCALL *GetCurveInfoStringCallback_pf) (FieldData_pa RawIndV,
                                                            FieldData_pa RawDepV,
                                                            CoordScale_e IndVCoordScale,
                                                            CoordScale_e DepVCoordScale,
                                                            LgIndex_t    NumRawPts,
                                                            EntIndex_t   XYMapNum,
                                                            char        *CurveSettings,
                                                            char       **CurveInfoString);
typedef Boolean_t (STDCALL *GetXYDataPointsCallback_pf) (FieldData_pa RawIndV,
                                                         FieldData_pa RawDepV,
                                                         CoordScale_e IndVCoordScale,
                                                         CoordScale_e DepVCoordScale,
                                                         LgIndex_t    NumRawPts,
                                                         LgIndex_t    NumCurvePts,
                                                         EntIndex_t   XYMapNum,
                                                         char        *CurveSettings,
                                                         double      *IndCurveValues,
                                                         double      *DepCurveValues);
typedef Boolean_t (STDCALL *GetProbeValueCallback_pf) (FieldData_pa RawIndV,
                                                       FieldData_pa RawDepV,
                                                       CoordScale_e IndVCoordScale,
                                                       CoordScale_e DepVCoordScale,
                                                       LgIndex_t    NumRawPts,
                                                       LgIndex_t    NumCurvePts,
                                                       EntIndex_t   XYMapNum,
                                                       char        *CurveSettings,
                                                       double       ProbeIndValue,
                                                       double      *ProbeDepValue);



#if defined MSWIN
typedef Boolean_t (STDCALL *PreTranslateMessage_pf)(MSG *pMsg);
#endif

/*********************************************************
 * Add-on Timers
 *********************************************************/
typedef Boolean_t (STDCALL *AddOnTimerCallback_pf) (ArbParam_t ClientData);

/* CORE SOURCE CODE REMOVED */




/* CORE SOURCE CODE REMOVED */

#endif /* _GLOBAL_H */
