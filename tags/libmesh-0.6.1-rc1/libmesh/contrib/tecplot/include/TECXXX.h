/*
 * TECXXX.h: Copyright (C) 1988-98 Amtec Engineering, Inc.
 */

#if !defined TECXXX_H_
#define TECXXX_H_

#if !defined CRAY
#  define TECINI  tecini
#  define TECZNE  teczne
#  define TECDAT  tecdat
#  define TECNOD  tecnod
#  define TECGEO  tecgeo
#  define TECTXT  tectxt
#  define TECLAB  teclab
#  define TECFIL  tecfil
#  define TECEND  tecend
#  define TECUSR  tecusr
#endif


#define INTEGER4  int

#if !defined MSWIN
# if defined (_WINDOWS) || defined (WIN32)
#   define MSWIN
# endif /* _WINDOWS || WIN32 */
#endif /* !MSWIN */

#if !defined (EXTERNC)
# if defined (__cplusplus)
#  define EXTERNC extern "C"
# else
#  define EXTERNC
# endif /* __cplusplus */
#endif /* EXTERN_C */

#if !defined (STDCALL)
# if defined MSWIN
#  define STDCALL __stdcall
# else /* !MSWIN */
#  define STDCALL
# endif /* MSWIN */
#endif /* STDCALL */

#if !defined (DLLEXPORT)
# if defined (MSWIN)
#  define DLLEXPORT _declspec (dllexport)
# else
#  define DLLEXPORT
# endif /* MSWIN */
#endif /* DLLEXPORT */

#if !defined (DLLIMPORT)
# if defined (MSWIN)
#  define DLLIMPORT _declspec (dllimport)
# else
#  define DLLIMPORT
# endif /* MSWIN */
#endif /* DLLIMPORT */


#if defined (TECPLOTKERNEL)
# define LIBCALL
# define LIBFUNCTION
#elif defined (MAKEARCHIVE)
# define LIBCALL STDCALL
# define LIBFUNCTION EXTERNC DLLEXPORT
#else /* !TECPLOTKERNAL && !MAKEARCHIVE */
# define LIBCALL STDCALL
# define LIBFUNCTION EXTERNC DLLIMPORT
#endif


LIBFUNCTION INTEGER4 LIBCALL TECINI(char     *Title,
                                    char     *Variables,
                                    char     *FName,
                                    char     *ScratchDir,
                                    INTEGER4 *Debug,
                                    INTEGER4 *VIsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECZNE(char     *ZoneTitle,
                                    INTEGER4 *IMx,
                                    INTEGER4 *JMx,
                                    INTEGER4 *KMx,
                                    char     *ZFormat,
                                    char     *DupList);

LIBFUNCTION INTEGER4 LIBCALL TECDAT(INTEGER4  *N,
                                    void      *FieldData,
                                    INTEGER4  *IsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECNOD(INTEGER4 *NData);

LIBFUNCTION INTEGER4 LIBCALL TECEND(void);

LIBFUNCTION INTEGER4 LIBCALL TECLAB(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECUSR(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECGEO(double    *XPos,
                                    double    *YPos,
                                    double    *ZPos,
                                    INTEGER4  *PosCoordMode,
                                    INTEGER4  *AttachToZone,
                                    INTEGER4  *Zone,
                                    INTEGER4  *Color,
                                    INTEGER4  *FillColor,
                                    INTEGER4  *IsFilled,
                                    INTEGER4  *GeomType,
                                    INTEGER4  *LinePattern,
                                    double    *PatternLength,
                                    double    *LineThickness,
                                    INTEGER4  *NumEllipsePts,
                                    INTEGER4  *ArrowheadStyle,
                                    INTEGER4  *ArrowheadAttachment,
                                    double    *ArrowheadSize,
                                    double    *ArrowheadAngle,
                                    INTEGER4  *Scope,
                                    INTEGER4  *NumSegments,
                                    INTEGER4  *NumSegPts,
                                    float     *XGeomData,
                                    float     *YGeomData,
                                    float     *ZGeomData,
                                    char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECTXT(double    *XPos,
                                    double    *YPos,
                                    INTEGER4  *PosCoordMode,
                                    INTEGER4  *AttachToZone,
                                    INTEGER4  *Zone,
                                    INTEGER4  *BFont,
                                    INTEGER4  *FontHeightUnits,
                                    double    *FontHeight,
                                    INTEGER4  *BoxType,
                                    double    *BoxMargin,
                                    double    *BoxLineThickness,
                                    INTEGER4  *BoxColor,
                                    INTEGER4  *BoxFillColor,
                                    double    *Angle,
                                    INTEGER4  *Anchor,
                                    double    *LineSpacing,
                                    INTEGER4  *TextColor,
                                    INTEGER4  *Scope,
                                    char      *Text,
                                    char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECFIL(INTEGER4 *F);

#endif /* TECXXX_H_ */
