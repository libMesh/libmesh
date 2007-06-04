/*****************************************************************
 *****************************************************************
 *******                                                  ********
 ****** (C) Copyright 1988-2000 by Amtec Engineering, Inc. *******
 *******                                                  ********
 *****************************************************************
 *****************************************************************/
/* CORE SOURCE CODE REMOVED */

#ifndef _MASTER_H_
#define _MASTER_H_

#ifdef NDEBUG
# ifdef _DEBUG
#   error "Both NDEBUG and _DEBUG defined"
# endif
#else
# ifndef _DEBUG
#   define _DEBUG
# endif
#endif

#if defined OPENGL
#  define USE_3D_HARDWARE
#endif

#ifndef THREED
#  define THREED
#endif
#ifndef USEENUM
#  define USEENUM
#endif
#ifndef DDDV
#  define DDDV
#endif

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#if defined QUICKDEMO
#define DEMO
#endif

#if defined ELM5
#define ELM
#endif

#if !defined USEENUM
#define USEENUM
#endif

#if defined MicrosoftC
#define DOS
#endif

#if defined (NETRESTRICT) || defined (ACTIVATOR)
#define RESTRICT
#endif

#if defined CRAYX
#define CRAY
#endif

#if defined IRISX
#define IRIS
#endif

#if defined HPX
#define HPUX
#define HP
#endif

#if defined IBMRS6000X
#define IBMRS6000
#endif

#if defined COMPAQALPHAX
#define COMPAQALPHA
#define COMPAQX
#define COMPAQ
#endif

#if defined DECALPHAX
#define DECALPHA
#define DECX
#endif

#if defined DECX
#define DEC
#endif

#if defined SUNSOLARISX || defined SUNSOLARIS86X
#define SUNX
#endif

#if defined SUNX
#define SUN
#endif

#if defined IRISX || defined CRAYX || defined HPX || defined SUNX || defined CONVEXX
#define UNIXX
#define SYSV
#endif

#if defined DECX || defined LINUX || defined IBMRS6000X || defined COMPAQX
#define UNIXX
#endif


#ifdef MSWIN
/* CORE SOURCE CODE REMOVED */

#ifndef TECPLOTKERNEL
#define Widget int
#endif

#endif /* MSWIN */

#if defined SYSV || defined MSWIN
#include <string.h>
#else
#include <strings.h>
#endif

#if defined (MicrosoftC)
#include <stdlib.h>
#define EXECOS
#ifndef FAR
#define FAR
#endif
#define VOID       void
#endif

#include <sys/types.h>
#include <stdlib.h>

#if defined UNIXX
#define X11
#define MOTIF
#define FAR
#define NEAR
#include <unistd.h>
#endif

/* CORE SOURCE CODE REMOVED */
/*
 * If not building the tecplot kernel then at least
 * include the X Instrinsics.  This will make most
 * development for addons etc work.
 */
#if defined MOTIF
#  if defined TECPLOTKERNEL
#    include <Xm/Xm.h>
#    if XmVERSION == 1 && XmREVISION == 0
       typedef void*   XtPointer;
#    endif
#  else
#    include <X11/Intrinsic.h>
#  endif
#endif

#if defined MOTIF
#define CREATE_DIALOG_PARAMS Widget W
typedef Widget ComboBoxWidget_t;
typedef Widget DropDownListWidget_t;
typedef Widget FileDialogWidget_t;
typedef Widget LabelWidget_t;
typedef Widget ListWidget_t;
typedef Widget OptionMenuWidget_t;
typedef Widget PullDownMenuWidget_t;
typedef Widget ScaleWidget_t;
typedef Widget TextFieldWidget_t;
typedef Widget ToggleWidget_t;
typedef Widget ButtonWidget_t;
typedef Widget GridWidget_t;
#endif
#if defined MSWIN
#include <windows.h>
#define CREATE_DIALOG_PARAMS     CWnd *, LaunchDialogMode_e
typedef Widget ComboBoxWidget_t;
typedef Widget DropDownListWidget_t;
typedef Widget FileDialogWidget_t;
typedef Widget LabelWidget_t;
typedef Widget ListWidget_t;
typedef Widget OptionMenuWidget_t;
typedef Widget PullDownMenuWidget_t;
typedef Widget ScaleWidget_t;
typedef Widget TextFieldWidget_t;
typedef Widget ToggleWidget_t;
typedef Widget ButtonWidget_t;
typedef Widget GridWidget_t;
#endif


#if !defined (TRACE) /* Assume that if TRACE is not defined, then none of the TRACE macros are */
#  if !defined (NDEBUG) /* DEBUG */
#    if defined (MSWIN)
        /* If the add-on is running in debug mode but does not
         * use MFC, then no TRACE macro is available. Thus, to make tracing available,
         * map TRACE to the win32 OutpuDebugString() function.
         */
#       define TRACE(str) do { OutputDebugString(str); } while (0)
#       define TRACE0(str) TRACE(str)
#       define TRACE1(str,a1) do { char s[5000]; sprintf(s,str,a1); OutputDebugString(s); } while (0)
#       define TRACE2(str,a1,a2) do { char s[5000]; sprintf(s,str,a1,a2); OutputDebugString(s); } while (0)
#       define TRACE3(str,a1,a2,a3) do { char s[5000]; sprintf(s,str,a1,a2,a3); OutputDebugString(s); } while (0)
#    else /* !MSWIN */
#       define TRACE(str) do {printf(str); fflush(stdout);} while (0)
#       define TRACE0(str) TRACE(str)
#       define TRACE1(str,a1) do {printf(str,a1); fflush(stdout);} while (0)
#       define TRACE2(str,a1,a2) do {printf(str,a1,a2); fflush(stdout);} while (0)
#       define TRACE3(str,a1,a2,a3) do {printf(str,a1,a2,a3); fflush(stdout);} while (0)
#    endif
#  else /* !DEBUG */
#    define TRACE(str)           ((void)0)
#    define TRACE0(str)          ((void)0)
#    define TRACE1(str,a1)       ((void)0)
#    define TRACE2(str,a1,a2)    ((void)0)
#    define TRACE3(str,a1,a2,a3) ((void)0)
#  endif
#endif /* !defined (TRACE) */


/*
  Platform independent way for add-ons to know how much space
  to allocate for a filename.
*/

#if !defined (MaxCharsFilePath)
# if defined (MSWIN)
#   define MaxCharsFilePath _MAX_PATH
# else
#   define MaxCharsFilePath 512
# endif /* MSWIN */
#endif /* !MaxCharsFilePath */

/* CORE SOURCE CODE REMOVED */

#endif
