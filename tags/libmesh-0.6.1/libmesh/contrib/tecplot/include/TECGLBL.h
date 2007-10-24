/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** (C) Copyright 1989-1998  by AMTEC ENGINEERING INC. *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/
/* 
 * TECGLBL.h  .... GLOBAL include file for all tecutil?.c files.
 */
/* CORE SOURCE CODE REMOVED */

/* CORE SOURCE CODE REMOVED */

#if defined MOTIF
#  define MANAGESTATE
#else /* MSWIN */
# if defined AFX_MANAGE_STATE
#   ifdef _USRDLL
#     define AFXMANAGESTATE AfxGetStaticModuleState()
#   else
#     define AFXMANAGESTATE AfxGetAppModuleState()
#   endif 
#   define MANAGESTATE AFX_MANAGE_STATE(AFXMANAGESTATE);
# else /* If not using MFC, then AFX_MANAGE_STATE is not available */
#   define MANAGESTATE
# endif /* AFX_MANAGE_STATE */
#endif

/* CORE SOURCE CODE REMOVED */


/* flag for TecUtilSet... functions */
#define TECUTILSETNOTMEMBER  (BAD_SET_VALUE + 1) /* int */

/* CORE SOURCE CODE REMOVED */

/* public ArgList structure */
typedef struct _ArgList_s *ArgList_pa;

/* implementation independent ID's for text and geometry */
#define TECUTILBADID 0 /* long */
typedef ArbParam_t     Text_ID;
typedef ArbParam_t     Geom_ID;

/* CORE SOURCE CODE REMOVED */
