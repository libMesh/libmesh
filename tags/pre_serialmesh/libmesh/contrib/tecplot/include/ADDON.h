#ifndef _ADDON_H
#define _ADDON_H

/* CORE SOURCE CODE REMOVED */


typedef void (STDCALL * StateChangeAddOnCallback_pf)(StateChange_e StateChange,
                                                     ArbParam_t    CallData);
typedef Boolean_t (STDCALL * MacroCommandAddOnCallback_pf)(char *CommandString,
                                                           char **ErrMsg);
typedef Boolean_t (STDCALL * MopupQueryAddOnCallback_pf)(void);
typedef void (STDCALL * AddOnInitializeFunction_pf)(void);


/* CORE SOURCE CODE REMOVED */

#endif /* _ADDON_H */
