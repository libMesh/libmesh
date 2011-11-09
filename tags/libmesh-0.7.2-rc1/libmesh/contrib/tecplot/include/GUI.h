#if !defined(GUI_H_)
#define GUI_H_

/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** (C) Copyright 1989-1999  by AMTEC ENGINEERING INC.********
*******       All Rights Reserved.                       ********
*******                                                  ********
*****************************************************************
*****************************************************************
*/

/* This file must be included *AFTER* "TECADDON.h" */

/* WINDOWS */
#if defined (MSWIN)
# if defined (STDCALL)
#   undef STDCALL
# endif 
# define STDCALL   __stdcall 
# if defined (__cplusplus)
#   define WINGUI_EXTERN_C extern "C"
# else
#   define WINGUI_EXTERN_C
# endif
# if !defined(LINKTOGUI)
#   define LINKTOGUI WINGUI_EXTERN_C __declspec ( dllimport )
# endif
# if !defined (MANAGESTATE)
#   define MANAGESTATE AFX_MANAGE_STATE(AfxGetStaticModuleState());
# endif
#endif

/* MOTIF */
#if !defined (MSWIN)
# if !defined (STDCALL)
#   define STDCALL
# endif
# if !defined (LINKTOGUI)
#   if defined __cplusplus
#     define LINKTOGUI extern "C"
#   else
#     define LINKTOGUI
#   endif
# endif
# if !defined (MANAGESTATE)
#   define MANAGESTATE
# endif
#endif


typedef int (*GUITextCallback_pf)  (const char *TextString);
typedef void (*GUIIntCallback_pf)   (const int *I);
typedef void (*GUIVoidCallback_pf)  (void);


#define BADDIALOGID  -2
#define MAINDIALOGID -1


/*
 Set/Unset the dialog to always be on top of other windows
*/

LINKTOGUI void STDCALL GUI_DialogSetTopmost(int Dialog,int MakeTopmost);

LINKTOGUI int STDCALL GUI_GetVersion(void);

LINKTOGUI int STDCALL GUI_DialogCreateModeless(int                 DialogParent,
                                               int                 Width,
                                               int                 Height,
                                               const char         *Title,
                                               GUIVoidCallback_pf  HelpButtonCallback,
                                               GUIVoidCallback_pf  CloseButtonCallback,
                                               GUIVoidCallback_pf  InitCallback);

LINKTOGUI int STDCALL GUI_DialogCreateModal(int                 DialogParent,
                                            int                 Width,
                                            int                 Height,
                                            const char         *Title,
                                            GUIVoidCallback_pf  HelpButtonCallback,
                                            GUIVoidCallback_pf  OkButtonCallback,
                                            GUIVoidCallback_pf  CancelButtonCallback,
                                            GUIVoidCallback_pf  InitCallback);


LINKTOGUI int STDCALL GUI_ButtonAdd(int                ParentDialog,
                                    int                X,
                                    int                Y,
                                    int                Width,
                                    int                Height,
                                    const char        *Label,
                                    GUIVoidCallback_pf ButtonCallback);

LINKTOGUI void STDCALL GUI_ButtonSetText(int Button,
                                         const char *NewText);




void GUI_Assert(
    const char *expression,  /* text representation of the assertion */
    const char *explanation, /* text description of the assertion */
    const char *utility);    /* name of the utility containing the assertion */



LINKTOGUI void STDCALL GUI_SetSensitivity(int Control,
                                          int IsSensitive);

LINKTOGUI void STDCALL GUI_SetVisibility(int Control,
                                         int IsVisible);

LINKTOGUI int STDCALL GUI_OptionMenuAdd(int               ParentDialog,
                                        int               X,
                                        int               Y,
                                        int               Width,
                                        int               Height,
                                        const char       *OptionList,
                                        GUIIntCallback_pf ValueChangedCallback);

LINKTOGUI void STDCALL GUI_OptionMenuSet(int OptionMenu,
                                         int Selection);

LINKTOGUI int STDCALL GUI_OptionMenuGet(int OptionMenu);


LINKTOGUI int STDCALL GUI_ListAdd(int               ParentDialog,
                                  int               X,
                                  int               Y,
                                  int               Width,
                                  int               Height,
                                  int               IsMultiSelection,
                                  GUIIntCallback_pf ValueChangedCallback);

LINKTOGUI int STDCALL GUI_ListGetItemCount(int List);

LINKTOGUI void STDCALL GUI_ListAppendItem(int  List,
                                          const char *Item);

LINKTOGUI char * STDCALL GUI_ListGetString(int List,
                                           int Position);

LINKTOGUI void STDCALL GUI_ListReplaceItem(int  List,
                                           const char *Item,
                                           int  Position);

LINKTOGUI void STDCALL GUI_ListDeleteAllItems(int List);

LINKTOGUI void STDCALL GUI_ListDeleteItemAtPos(int List,
                                               int Position);

LINKTOGUI void STDCALL GUI_ListDeselectAllItems(int List);

LINKTOGUI void STDCALL GUI_ListSetSelectedItem(int List,
                                               int Position);

LINKTOGUI void STDCALL GUI_ListGetSelectedItems(int List,
                                                int **SelectedItemList,
                                                int *SelectedItemCount);

LINKTOGUI void STDCALL GUI_ListSetSelectedItems(int List,
                                                int *SelectedItemList,
                                                int SelectedItemCount);
LINKTOGUI void STDCALL GUI_ListDeallocItemList(int **ItemList);

LINKTOGUI int STDCALL GUI_ToggleAdd(int                ParentDialog,
                                    int                X,
                                    int                Y,
                                    int                Width,
                                    int                Height,
                                    const char        *Label,
                                    GUIIntCallback_pf  ValueChangedCallback);

LINKTOGUI void STDCALL GUI_ToggleSet(int Toggle,
                                     int SetOn);



LINKTOGUI int STDCALL GUI_ToggleGet(int Toggle);

LINKTOGUI int STDCALL GUI_RadioBoxAdd(int               ParentDialog,
                                      int               X,
                                      int               Y,
                                      int               Width,
                                      int               Height,
                                      const char       *Label1,
                                      const char       *Label2,
                                      const char       *Label3,
                                      const char       *Label4,
                                      const char       *Label5,
                                      GUIIntCallback_pf ValueChangedCallback);

LINKTOGUI void STDCALL GUI_RadioBoxSetToggle(int RadioBox,
                                             int ToggleNumber);


LINKTOGUI int STDCALL GUI_RadioBoxGetToggle(int RadioBox);

LINKTOGUI int STDCALL GUI_LabelAdd(int      ParentDialog,
                                   int      X,
                                   int      Y,
                                   const char *Label);

LINKTOGUI void STDCALL GUI_LabelSetText(int  Label,
                                        const char *LabelString);

LINKTOGUI int STDCALL GUI_TextFieldAdd(int                ParentDialog,
                                       int                X,
                                       int                Y,
                                       int                Width,
                                       int                Height,
                                       GUITextCallback_pf ValueChangedCallback);

LINKTOGUI int STDCALL GUI_TextAdd(int                ParentDialog,
                                  int                X,
                                  int                Y,
                                  int                Width,
                                  int                Height,
                                  int          IsReadOnly,
                                  GUITextCallback_pf ValueChangedCallback);
LINKTOGUI void STDCALL GUI_TextSetInsertPos(int Text,
                                            int Position);
LINKTOGUI void STDCALL GUI_TextSetMinInsertPos(int Text);
LINKTOGUI void STDCALL GUI_TextSetMaxInsertPos(int Text);
LINKTOGUI void STDCALL GUI_TextSetString(int          Text,
                                         const char  *TextString);
LINKTOGUI char * STDCALL GUI_TextGetString(int Text);
LINKTOGUI void STDCALL GUI_TextInsertString(int          Text,
                                            const char  *TextString);
LINKTOGUI void STDCALL GUI_TextAppendString(int          Text,
                                            const char  *TextString);

LINKTOGUI int STDCALL GUI_ScaleAdd(int               ParentDialog,
                                   int               X,
                                   int               Y,
                                   int               Width,
                                   int               Height,
                                   int               ScaleMin,
                                   int               ScaleMax,
                                   int               DecimalPrecision,
                                   GUIIntCallback_pf ValueChangedCallback,
                                   GUIIntCallback_pf DragValueChangedCallback);

LINKTOGUI void STDCALL GUI_ScaleSetValue(int Scale,
                                         int NewValue);

LINKTOGUI void STDCALL GUI_ScaleSetLimits(int Scale,
                                          int ScaleMin,
                                          int ScaleMax,
                                          int DecimalPrecision);

LINKTOGUI int STDCALL GUI_ScaleGetValue(int Scale);

LINKTOGUI int STDCALL GUI_VertSeparatorAdd(int    ParentDialog,
                                           int    X,
                                           int    Y,
                                           int    Height);

LINKTOGUI int STDCALL GUI_HorzSeparatorAdd(int    ParentDialog,
                                           int    X,
                                           int    Y,
                                           int    Width);

LINKTOGUI int STDCALL GUI_FrameAdd(int     ParentDialog,
                                   int     X,
                                   int     Y,
                                   int     Width,
                                   int     Height,
                                   const char *Label /* may be NULL */);

LINKTOGUI void STDCALL GUI_TextFieldSetString(int   TextField,
                                              const char *TextString);
LINKTOGUI char * STDCALL GUI_TextFieldGetString(int TextField);

LINKTOGUI void STDCALL GUI_DialogLaunch(int DialogID);

LINKTOGUI void STDCALL GUI_DialogDismiss(int DialogID);

LINKTOGUI int STDCALL GUI_DialogIsUp(int DialogID);

LINKTOGUI void STDCALL GUI_DialogSetTitle(int DialogID,
                                          const char *NewTitle);

LINKTOGUI void STDCALL GUI_TextAppendString(int          Text,
                                            const char  *TextString);
                                            
LINKTOGUI int STDCALL GUI_MenuBarAdd(int ParentDialog);

LINKTOGUI int STDCALL GUI_MenuAdd(int         ParentMenu,
                                  const char *Label);

LINKTOGUI int STDCALL GUI_MenuAddItem(int                 ParentMenu,
                                      const char         *Label,
                                      GUIVoidCallback_pf  Callback);

LINKTOGUI int STDCALL GUI_MenuAddToggle(int                ParentMenu,
                                        const char        *Label,
                                        GUIIntCallback_pf  Callback);

LINKTOGUI void STDCALL GUI_MenuAddSeparator(int ParentMenu);

LINKTOGUI void STDCALL GUI_MenuItemSetText(int         MenuItem,
                                           const char *NewText);

LINKTOGUI void STDCALL GUI_MenuSetToggle(int MenuItem,
                                         int SetOn);

LINKTOGUI void STDCALL GUI_MenuDeleteItem(int MenuItem);




/* this line must be last */
#endif /* GUI_H_ */


