/*
 * NOTICE and LICENSE for Tecplot Input/Output Library (TecIO) - OpenFOAM
 *
 * Copyright (C) 1988-2009 Tecplot, Inc.  All rights reserved worldwide.
 *
 * Tecplot hereby grants OpenCFD limited authority to distribute without
 * alteration the source code to the Tecplot Input/Output library, known 
 * as TecIO, as part of its distribution of OpenFOAM and the 
 * OpenFOAM_to_Tecplot converter.  Users of this converter are also hereby
 * granted access to the TecIO source code, and may redistribute it for the
 * purpose of maintaining the converter.  However, no authority is granted
 * to alter the TecIO source code in any form or manner.
 *
 * This limited grant of distribution does not supersede Tecplot, Inc.'s 
 * copyright in TecIO.  Contact Tecplot, Inc. for further information.
 * 
 * Tecplot, Inc.
 * 3535 Factoria Blvd, Ste. 550
 * Bellevue, WA 98006, USA
 * Phone: +1 425 653 1200
 * http://www.tecplot.com/
 *
 */
#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2008 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#define STRLISTMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "STRUTIL.h"
#include "ALLOC.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "ARRLIST.h"
#include "STRLIST.h"

/* END HEADER */

using namespace tecplot::strutil;

/*
 * This set of functions provide a wrapper around the array list utilities
 * thereby making it aware of item allocation and deallocation. All strings
 * given to the string list and returned to the client are copies. Therefore
 * it is the client's responsibility to deallocate string results when no
 * longer needed.
 */


/*
 * Destructor for cleaning up string allocations.
 *
 * param ItemRef
 *     Reference to the string item to destroy.
 * param ClientData
 *     Any client data needed for destroying the string.
 *
 * return
 *     TRUE is a requirement
 */
static Boolean_t StringListItemDestructor(void       *ItemRef,
                                          ArbParam_t ClientData)
{
    char **StringRef = (char **)ItemRef;

    REQUIRE(VALID_REF(StringRef));
    REQUIRE(VALID_REF(*StringRef) || *StringRef == NULL);

    if (*StringRef != NULL)
    {
        FREE_ARRAY(*StringRef, "string");
        *StringRef = NULL;
    }

    ENSURE(*StringRef == NULL);
    return TRUE;
}

/*
 * String item duplicator.
 *
 * param TargetItemRef
 *     Reference to the string list item to receive the duplicate.
 * param SourceItemRef
 *     Reference to the string list item to duplicate.
 * param ClientData
 *     Any client data required for duplication.
 *
 * return
 *     TRUE if the duplication was a success
 *     FALSE otherwise. If the duplication failed it
 *     is the client's responsibility to cleanup any
 *     partial duplication
 */
static Boolean_t StringListItemDuplicator(void       *TargetItemRef,
                                          void       *SourceItemRef,
                                          ArbParam_t ClientData)
{
    Boolean_t IsOk = TRUE;
    char **TargetStringRef = (char **)TargetItemRef;
    char **SourceStringRef = (char **)SourceItemRef;

    REQUIRE(VALID_REF(TargetStringRef));
    REQUIRE(VALID_REF(SourceStringRef));
    REQUIRE(VALID_REF(*SourceStringRef) || *SourceStringRef == NULL);

    if (*SourceStringRef != NULL)
        IsOk = ((*TargetStringRef = DupString(dontTranslate(*SourceStringRef))) != NULL);
    else
        *TargetStringRef = NULL;

    ENSURE(VALID_REF(*TargetStringRef) || *TargetStringRef == NULL);
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/*
 * Determine if the string list handle and its members are sane.
 */
Boolean_t StringListValid(StringList_pa stringList)
{
    Boolean_t isValid = ArrayListIsValid((ArrayList_pa)stringList);

    if (isValid)
    {
        LgIndex_t stringCount = ArrayListGetCount((ArrayList_pa)stringList);

#if defined PERFORM_EXPENSIVE_STRLIST_TESTS
        {
            for (LgIndex_t index = 0; index < stringCount; index++)
            {
                char *string = ArrayListGetCharPtr((ArrayList_pa)stringList, index);
                if (string != NULL && !VALID_REF(string))
                {
                    isValid = FALSE;
                    break;
                }
            }
#else
        {
            /* Check first and last only */
            if (stringCount > 0)
            {
                char *string = ArrayListGetCharPtr((ArrayList_pa)stringList, 0);
                if (string != NULL && !VALID_REF(string))
                {
                    isValid = FALSE;
                }
            }
            if (isValid && stringCount > 1)
            {
                char *string = ArrayListGetCharPtr((ArrayList_pa)stringList, stringCount - 1);
                if (string != NULL && !VALID_REF(string))
                {
                    isValid = FALSE;
                }
            }
        }
#endif /* PERFORM_SKIP_EXPENSIVE_STRLIST_TESTS */
        }

        ENSURE(VALID_BOOLEAN(isValid));
        return isValid;
    }


    /*
     * Remove all members of the string list.
     */
    void StringListClear(StringList_pa StringList)
    {
        REQUIRE(StringListValid(StringList));

        ArrayListDeleteAllItems((ArrayList_pa)StringList, StringListItemDestructor, 0);

        ENSURE(StringListValid(StringList) && StringListCount(StringList) == 0);
    }


    /*
     * Remove 'Count' strings from the list beginning at the specified offset.
     * The members following the items removed are shifted down accordingly to
     * fill the vacated space.
     */
    void StringListRemoveStrings(StringList_pa StringList,
                                 LgIndex_t     StringOffset,
                                 LgIndex_t     Count)
    {
        REQUIRE(StringListValid(StringList));
        REQUIRE(0 <= StringOffset && StringOffset <= StringListCount(StringList) - 1);
        REQUIRE(1 <= Count && StringOffset + Count <= StringListCount(StringList));

        ArrayListDeleteItems((ArrayList_pa)StringList, StringOffset, Count,
                             StringListItemDestructor, 0);

        ENSURE(StringListValid(StringList));
    }


    /*
     * Remove the string from the list at the specified offset. The members
     * following the item removed are shifted down accordingly to fill the
     * vacated space.
     */
    void StringListRemoveString(StringList_pa StringList,
                                LgIndex_t     StringOffset)
    {
        REQUIRE(StringListValid(StringList));
        REQUIRE(0 <= StringOffset && StringOffset <= StringListCount(StringList) - 1);

        ArrayListDeleteItems((ArrayList_pa)StringList, StringOffset, 1,
                             StringListItemDestructor, 0);

        ENSURE(StringListValid(StringList));
    }


    /*
     * Deallocate the string list members and handle and set the handle to NULL.
     */
    void LIBCALL StringListDealloc(StringList_pa *StringList)
    {
        REQUIRE(VALID_REF(StringList));
        REQUIRE(*StringList == NULL || StringListValid(*StringList));

        if (*StringList != NULL)
            ArrayListDealloc((ArrayList_pa *)StringList, StringListItemDestructor, 0);

        ENSURE(*StringList == NULL);
    }


    /*
     * Allocate a string list handle. A handle of NULL is
     * returned if sufficient memory is not available.
     */
    StringList_pa StringListAlloc(void)
    {
        StringList_pa Result;

        Result = (StringList_pa)ArrayListAlloc(0, ArrayListType_CharPtr, NULL, 0);

        ENSURE(Result == NULL || StringListValid(Result));
        return Result;
    }


    /*
     * Append a copy of the string to the string list. The string list will be
     * expanded to accommodate the additional item. A return value of TRUE
     * indicates the operation was successful otherwise FALSE is returned
     * indicating that sufficient memory was not available for the additional
     * item.
     */
    Boolean_t StringListAppendString(StringList_pa StringList,
                                     const char   *String)
    {
        Boolean_t IsOk;

        REQUIRE(StringListValid(StringList));
        REQUIRE(String == NULL || VALID_REF(String));

        IsOk = StringListSetString(StringList, StringListCount(StringList), String);

        ENSURE(StringListValid(StringList));
        ENSURE(VALID_BOOLEAN(IsOk));
        return IsOk;
    }


    /*
     * Return the number of strings currently in the string list.
     */
    LgIndex_t LIBCALL StringListCount(StringList_pa StringList)
    {
        LgIndex_t Result;

        REQUIRE(StringListValid(StringList));

        Result = ArrayListGetCount((ArrayList_pa)StringList);

        ENSURE(Result >= 0);
        return Result;
    }


    /*
     * Return a copy of the string at the specified offset in the string list.
     */
    char * LIBCALL StringListGetString(StringList_pa StringList,
                                       LgIndex_t     StringOffset)
    {
        char *Result;
        const char *StringRef;

        REQUIRE(StringListValid(StringList));
        REQUIRE(0 <= StringOffset && StringOffset <= StringListCount(StringList) - 1);

        StringRef = StringListGetStringRef(StringList, StringOffset);
        if (StringRef == NULL)
            Result = NULL;
        else
            Result = DupString(dontTranslate(StringRef));

        ENSURE(Result == NULL || VALID_REF(Result));
        return Result;
    }


#if !defined USE_MACROS_FOR_FUNCTIONS
    /*
     * Returns actual string at the specified offset in the string list.  Do not
     * attempt to free this string.  Changing this string should be done with
     * utmost caution.
     */
    const char *StringListGetStringRef_FUNC(StringList_pa StringList,
                                            LgIndex_t     StringOffset)
    {
        const char *Result;

        REQUIRE(StringListValid(StringList));
        REQUIRE(0 <= StringOffset && StringOffset <= StringListCount(StringList) - 1);

        Result = StringListGetStringRef_MACRO(StringList, StringOffset);

        ENSURE(Result == NULL || VALID_REF(Result));
        return Result;
    }
#endif


    /*
     * Place a copy of the specified string at the specified offset. If the offset
     * is beyond the end of the string list it is sized accordingly and the
     * intervening string references between the last item of the original
     * state and the last item of the new state are assigned NULL. If a string
     * already exists at the specified location its resources are released.
     * A return value of TRUE indicates the operation was successful otherwise
     * FALSE is returned indicating that sufficient memory was not available
     * for the additional item at the specified offset.
     */
    Boolean_t StringListSetString(StringList_pa StringList,
                                  LgIndex_t     StringOffset,
                                  const char   *String)
    {
        Boolean_t       IsOk;
        ArrayListItem_u ItemCopy;

        REQUIRE(StringListValid(StringList));
        REQUIRE(StringOffset >= 0);
        REQUIRE(String == NULL || VALID_REF(String));

        if (String != NULL)
        {
            ItemCopy.CharPtr = DupString(dontTranslate(String));
            IsOk = (ItemCopy.CharPtr != NULL);
        }
        else
        {
            ItemCopy.CharPtr = NULL;
            IsOk = TRUE;
        }

        if (IsOk)
            IsOk = ArrayListSetItem((ArrayList_pa)StringList, StringOffset, ItemCopy,
                                    StringListItemDestructor, 0);

        ENSURE(StringListValid(StringList));
        ENSURE(VALID_BOOLEAN(IsOk));
        return IsOk;
    }


    /*
     * Insert a copy of the string into the string list at the specified offset.
     * The string list will be expanded to accommodate the additional item.
     * A return value of TRUE indicates the operation was successful otherwise
     * FALSE is returned indicating that sufficient memory was not available
     * for the additional item.
     */
    Boolean_t StringListInsertString(StringList_pa  StringList,
                                     LgIndex_t      StringOffset,
                                     const char    *String)
    {
        Boolean_t       IsOk;
        ArrayListItem_u ItemCopy;

        REQUIRE(StringListValid(StringList));
        REQUIRE(StringOffset >= 0);
        REQUIRE(String == NULL || VALID_REF(String));

        if (String != NULL)
        {
            ItemCopy.CharPtr = DupString(dontTranslate(String));
            IsOk = (ItemCopy.CharPtr != NULL);
        }
        else
        {
            ItemCopy.CharPtr = NULL;
            IsOk = TRUE;
        }

        if (IsOk)
            IsOk = ArrayListInsertItem(
                       (ArrayList_pa)StringList, StringOffset, ItemCopy);

        ENSURE(StringListValid(StringList));
        ENSURE(VALID_BOOLEAN(IsOk));
        return IsOk;
    }


    /*
     * Return a handle to a duplicate of the specified string list and its contents.
     * A handle of NULL is returned if sufficient memory is not available.
     */
    StringList_pa StringListCopy(StringList_pa StringList)
    {
        StringList_pa Result;

        REQUIRE(StringListValid(StringList));

        Result = (StringList_pa)ArrayListCopy((ArrayList_pa)StringList,
                                              StringListItemDuplicator, 0);

        ENSURE(Result == NULL ||
               (StringListValid(Result) &&
                StringListCount(Result) == StringListCount(StringList)));
        return Result;
    }



    /*
     * Append a copy of the contents of the source list to the target list.
     * A return value of TRUE indicates the operation was successful otherwise
     * FALSE is returned indicating that sufficient memory was not available
     * for the request.
     */
    Boolean_t StringListAppend(StringList_pa Target,
                               StringList_pa Source)
    {
        Boolean_t     IsOk;
        StringList_pa SourceCopy;

        REQUIRE(StringListValid(Target));
        REQUIRE(StringListValid(Source));

        SourceCopy = StringListCopy(Source);
        IsOk = (SourceCopy != NULL);
        if (IsOk)
        {
            ArrayListAppend((ArrayList_pa)Target, (ArrayList_pa)SourceCopy);
            /* deallocate the list but not the string items since Target now owns them */
            ArrayListDealloc((ArrayList_pa *)(void *)&SourceCopy, NULL, 0);
        }

        ENSURE(StringListValid(Target));
        ENSURE(VALID_BOOLEAN(IsOk));
        return IsOk;
    }



    /*
     * Return a new line, '\n', separated string representation of the string list.
     * Caller is responsible for de-allocating the result.
     */
    char *StringListToNLString(StringList_pa StringList)
    {
        char  *Result;
        int Count;
        size_t Length = 0;

        REQUIRE(StringListValid(StringList));

        /* determine the resulting new line, '\n', separated string length */
        Count = StringListCount(StringList);
        if (Count >= 1)
        {
            int Index;
            for (Index = 0, Length = strlen("\n") * (Count - 1);
                 Index < Count;
                 Index++)
            {
                char *String = ArrayListGetCharPtr((ArrayList_pa)StringList, Index);
                if (String != NULL)
                    Length += strlen(String);
            }
        }

        /* create a new line, '\n', separated string */
        Result = ALLOC_ARRAY(Length + 1, char, "new line separated string");
        if (Result != NULL)
        {
            int Index;
            for (Index = 0, strcpy(Result, "");
                 Index < Count;
                 Index++)
            {
                char *String = ArrayListGetCharPtr(
                                   (ArrayList_pa)StringList, Index);

                if (Index != 0)
                    strcat(Result, "\n");

                if (String != NULL)
                    strcat(Result, String);
            }
        }

        ENSURE(Result == NULL || VALID_REF(Result));
        return Result;
    }


    /*
     * Create a string list from the new line, '\n', separated string. The string
     * is copied and therefore owned and managed by the caller.
     */
    StringList_pa StringListFromNLString(const char *String)
    {
        StringList_pa Result;
        LgIndex_t     StartIndex;
        LgIndex_t     EndIndex;

        REQUIRE(VALID_REF(String));

        /* create the string list and scan the entire string */
        Result = StringListAlloc();
        for (StartIndex = EndIndex = 0; Result != NULL; EndIndex++)
        {
            /* end of sub-string ? */
            if (String[EndIndex] == '\n' || String[EndIndex] == '\0')
            {
                /* extract the sub-string and append it to the string list */
                LgIndex_t Length = EndIndex - StartIndex;
                char     *SubString = ALLOC_ARRAY(Length + 1, char, "sub string");
                if (SubString != NULL)
                {
                    CopySubString(SubString, String, StartIndex, Length);
                    StringListAppendString(Result, SubString);

                    FREE_ARRAY(SubString, "sub string");

                    if (String[EndIndex] != '\0')
                        StartIndex = EndIndex + 1;
                    else
                        break; /* nothing left to scan */
                }
                else
                {
                    /* memory allocation failure: bail out */
                    StringListDealloc(&Result);
                    Result = NULL;
                    break;
                }
            }
        }

        ENSURE(Result == NULL || StringListValid(Result));
        return Result;
    }


    /*
     * Return a 'C' string array representation of the string list.
     * Caller is responsible for de-allocating the result.
     */
    char **StringListToArray(StringList_pa StringList)
    {
        char    **Result;

        REQUIRE(StringListValid(StringList));

        Result = (char **)ArrayListToArray((ArrayList_pa)StringList,
                                           StringListItemDuplicator, 0);

        ENSURE(Result == NULL || VALID_REF(Result));
        return Result;
    }



    /*
     * Create a string list from the 'C' string array. The string array
     * is copied and therefore owned and managed by the caller.
     */
    StringList_pa StringListFromArray(const char **StringArray,
                                      LgIndex_t    Count)
    {
        StringList_pa Result;

        REQUIRE((Count == 0 && StringArray == NULL) ||
                (Count >= 1 && VALID_REF(StringArray)));

        Result = (StringList_pa)ArrayListFromArray((void *)StringArray,
                                                   Count, ArrayListType_CharPtr,
                                                   StringListItemDuplicator, 0);

        ENSURE(Result == NULL || StringListValid(Result));
        return Result;
    }



#define ISJOINCHAR(c) ((c == ';') || (c == '+'))

    static void SkipWhiteSpaceOrComma(const char **CPtr)
    {
        REQUIRE(VALID_REF(CPtr) && VALID_REF(*CPtr));
        while (ISWHITESPACE(**CPtr) || (**CPtr == ','))
            (*CPtr)++;
    }

    /*
     * Obtain the next sub-string.  This can be of the form:
     *
     *  [del]any-character-sequence[del]
     *
     *            or
     *
     *   limited-character-sequence
     *
     *  where a limited-character-sequence cannot contain
     * any of the following characters:  +;,<space>
     *
     */
    static Boolean_t GetNextSubString(const char **OriginalCPtr,
                                      char       **NextSubString)
    {
        Boolean_t   IsOk = TRUE;
        const char *CStart;
        const char *CPtr;
        char        InsideDelimiter = '\0';

        REQUIRE(VALID_REF(OriginalCPtr) && (VALID_REF(*OriginalCPtr)));
        REQUIRE(VALID_REF(NextSubString));

        *NextSubString = NULL;

        CPtr = *OriginalCPtr;
        SkipWhiteSpaceOrComma(&CPtr);

        if (*CPtr == '"' || *CPtr == '\'')
        {
            InsideDelimiter = *CPtr;
            CPtr++;
        }

        CStart = CPtr;

        while (*CPtr &&
               ((InsideDelimiter && (*CPtr != InsideDelimiter)) ||
                (!InsideDelimiter && (*CPtr != ',')       &&
                 !ISJOINCHAR(*CPtr)  &&
                 !ISWHITESPACE(*CPtr))))
        {
            if (InsideDelimiter  &&
                (*CPtr == '\\')  &&
                (*(CPtr + 1) == InsideDelimiter))
                CPtr += 2;
            else
                CPtr++;
        }

        if (InsideDelimiter && (*CPtr != InsideDelimiter))
            IsOk = FALSE;


        if (IsOk && CStart < CPtr)
        {
            size_t StrLen = (size_t)(CPtr - CStart);
            *NextSubString = ALLOC_ARRAY(StrLen + 1, char, "GetNextSubString: NextSubString");
            if (*NextSubString)
            {
                char *NPtr = *NextSubString;
                /*
                 * Don't just copy the string because escaped delimiters need to have
                 * the escape removed...
                 */
                while (CStart < CPtr)
                {
                    if ((*CStart == '\\') && (*(CStart + 1) == InsideDelimiter))
                        CStart++;
                    *NPtr++ = *CStart++;
                }
                *NPtr = '\0';
            }
            else
                IsOk = FALSE;
        }

        if (IsOk)
        {
            if (InsideDelimiter)
                CPtr++;
            SkipWhiteSpaceOrComma(&CPtr);
            *OriginalCPtr = CPtr;
        }

        ENSURE(VALID_BOOLEAN(IsOk));
        return IsOk;
    }




    /*
     * Return a string list representation of a compound string.
     *
     * The compound String parameter has the following form:
     *
     * [del]<character-sequence>[del] [GroupJoinCharacter] [del]<character-sequence>[del] [GroupJoinCharacter] .....
     *                       or
     * <nospace-character-sequence> <nospace-character-sequence> ...
     *
     * where:
     *   [del] is an optional single quote or a double quote.  [del] must be used
     *   if <character-sequence> contains spaces, commas, or the plus symbol.
     *
     *   GroupJoinCharacter can be either a "+" or a ";"
     *
     * The GroupJoinCharacter symbol is used to separate character sequences that
     * are to be grouped together. If the GroupJoinCharacter symbol is omitted then
     * a new group is started.
     *
     * Internally, the original string is converted to a list of strings where
     * each string uses newlines to separate one sub-string from the next.
     */
    StringList_pa StringListFromCompound(const char *String)
    {
        const char   *CPtr;
        StringList_pa Result;
        Boolean_t     IsOk = TRUE;
        char         *CurString = NULL;

        REQUIRE(VALID_REF(String));
        SkipWhiteSpaceOrComma(&String);
        REQUIRE(!ISJOINCHAR(*String));

        /* extract character sequences */
        Result = StringListAlloc();
        CPtr   = String;

        while (IsOk && *CPtr != '\0')
        {
            char     *NextSubString = NULL;
            Boolean_t WantsToJoin   = FALSE;

            if (ISJOINCHAR(*CPtr))
            {
                WantsToJoin = TRUE;
                CPtr++;
                SkipWhiteSpaceOrComma(&CPtr);
            }

            IsOk = GetNextSubString(&CPtr,
                                    &NextSubString);

            if (IsOk)
            {
                /*
                 * Tack on the sub-string to the running string.
                 */
                if (WantsToJoin)
                    TackOnChar(&CurString, '\n');
                if (NextSubString != NULL && strlen(NextSubString) != 0)
                    IsOk = TackOnString(&CurString, NextSubString, FALSE, FALSE);
                else if (CurString == NULL)
                    CurString = DupString(dontTranslate(""));
            }

            if (NextSubString != NULL)
                FREE_ARRAY(NextSubString, "StringListFromCompound: NextSubString");

            /*
             * If this is the end of processing or if the next character is
             * not a join character then add the current string to the stringlist.
             */

            if (IsOk && !ISJOINCHAR(*CPtr))
            {
                StringListAppendString(Result, CurString);
                if (CurString != NULL)
                    FREE_ARRAY(CurString, "current string");
                CurString = NULL;
            }
        }

        if (CurString != NULL)
            FREE_ARRAY(CurString, "current string");

        if (!IsOk)
            StringListDealloc(&Result);

        ENSURE(Result == NULL || StringListValid(Result));
        return Result;
    }


    /*
     * Return a compound string representation of a string list.
     *
     * One common usage in Tecplot:
     *   The $!OpenLayout command in tecplot has the sub-option
     *   ALTDATALOADINSTRUCTIONS that has the form:
     *     '"instr-string1" [GroupJoinCharacter] "instr-string2" [+] ...'
     */
    char *StringListToCompound(StringList_pa  StringList,
                               char           GroupJoinCharacter,
                               const char    *CharsToEscape)
    {
        Boolean_t IsOk = TRUE;
        LgIndex_t Index;
        LgIndex_t Count;
        char      *Result = NULL;

        REQUIRE(StringListValid(StringList));
        REQUIRE(StringListCount(StringList) >= 1);
        REQUIRE(ISJOINCHAR(GroupJoinCharacter));
        REQUIRE(VALID_REF(CharsToEscape));

        for (Index = 0, Count = StringListCount(StringList), IsOk = TRUE;
             Index < Count && IsOk;
             Index++)
        {
            char *String = StringListGetString(StringList, Index);

            if (String != NULL && strlen(String) != 0)
            {
                char       *CStart = NULL;
                char       *CEnd = NULL;
                char       *EscapedString = NULL;
                const char *EscChar = NULL;
                char       *StrChar = NULL;

                /* First scan the string and escape any specified characters.  */
                /* Note that the Escape sequence is a double backslash because */
                /* it the first escape escapes the escape for variable usage.  */
                for (StrChar = String; *StrChar != '\0'; StrChar++)
                {
                    for (EscChar = CharsToEscape; *EscChar != '\0'; EscChar++)
                        if (*StrChar == *EscChar)
                        {
                            IsOk = TackOnChar(&EscapedString, '\\');
                            IsOk = TackOnChar(&EscapedString, '\\');
                            break;
                        }
                    IsOk = TackOnChar(&EscapedString, *StrChar);
                }

                CEnd = EscapedString;
                while (IsOk && *CEnd != '\0')
                {
                    int  Len = 0;
                    char *TString;

                    CStart = CEnd;
                    while (*CEnd != '\0' && *CEnd != '\n')
                    {
                        Len++;
                        if (*CEnd == '"')
                            Len++;
                        CEnd++;
                    }

                    TString = ALLOC_ARRAY(Len + 4, char, "temp compound sub-string");
                    if (TString != NULL)
                    {
                        char *TStr;

                        /* prepend the new string with either   */
                        /* a space character or the plus symbol */
                        if (CStart == EscapedString)
                        {
                            if (Index != 0)
                                IsOk = TackOnChar(&Result, ' ');
                        }
                        else
                        {
                            IsOk = TackOnChar(&Result, GroupJoinCharacter);
                        }

                        /* stuff TString and append the new string */
                        TStr = TString;
                        *TStr++ = '"';
                        while (CStart != CEnd)
                        {
                            if (*CStart == '"')
                                *TStr++ = '\\';
                            *TStr++ = *CStart++;
                        }
                        *TStr++ = '"';
                        *TStr = '\0';

                        TackOnString(&Result, TString, FALSE, FALSE);
                        FREE_ARRAY(TString, "StringListToCompound");
                        TString = NULL;
                        if (*CEnd)
                            CEnd++;
                    }
                    else
                    {
                        IsOk = FALSE;
                    }
                }

                if (EscapedString != NULL)
                    FREE_ARRAY(EscapedString, "escaped string");
            }
            else
            {
                /* a null pointer or length of zero indicates an empty sub-string */
                if (Index == 0)
                    TackOnString(&Result, "\"\"", FALSE, FALSE);
                else
                    TackOnString(&Result, " \"\"", FALSE, FALSE);
            }

            if (String != NULL)
                FREE_ARRAY(String, "string list item");
        }

        if (!IsOk)
        {
            if (Result != NULL)
            {
                FREE_ARRAY(Result, "StringListToCompound");
                Result = NULL;
            }
        }

        ENSURE(Result == NULL || VALID_REF(Result));
        return Result;
    }


    /**
     * Holds the comparator function pointer.
     */
    static StringListStringComparator_pf ComparatorFunction = NULL;


    /**
     * Forwards the comparison test to the 'Comparator' supplied to the
     * 'StringListSort' function.
     *
     * param Item1
     *     Item to compare against Item2.
     * param Item2
     *     Item to compare against Item1.
     * param ClientData
     *     Contextual information that was passed to the 'ArrayListQSort' function.
     *
     * return
     *     -1: if Item1 is less than Item2
     *      0: if Item1 is equal to Item2
     *      1: if Item1 is greater than Item2
     */
    static int STDCALL ComparatorProxy(ArrayListItem_u Item1,
                                       ArrayListItem_u Item2,
                                       ArbParam_t      ClientData)
    {
        /* forward the request */
        return ComparatorFunction(Item1.CharPtr, Item2.CharPtr, ClientData);
    }


    /**
     * Compares two strings from a list string. Note that either string may be
     * NULL as StringLists allow for NULL elements.
     *
     * param String1
     *     String to compare against String2.
     * param String2
     *     String to compare against String1.
     * param ClientData
     *     Contextual information that was passed to the 'StringListSort' function.
     *
     * return
     *      - A value less than zero if String1 is less than String2.
     *      - A value of zero if String1 is equal to String2.
     *      - A value greater than zero if String1 is greater than String2.
     */
    static int STDCALL DefaultStrcmpComparator(const char *String1,
                                               const char *String2,
                                               ArbParam_t  ClientData)
    {
        int Result = 0; /* ...quite compiler */

        REQUIRE(VALID_REF(String1) || String1 == NULL);
        REQUIRE(VALID_REF(String2) || String2 == NULL);

        if (String1 != NULL && String2 != NULL)
        {
            Result = strcmp(String1, String2);
            if (Result < 0)
                Result = -1;
            else if (Result > 0)
                Result = 1;
        }
        else if (String1 == NULL && String2 == NULL)
            Result = 0;
        else if (String1 == NULL)
            Result = -1;
        else if (String2 == NULL)
            Result = 1;
        else
            CHECK(FALSE);

        ENSURE((Result == -1) || (Result == 0) || (Result == 1));
        return Result;
    }

    /**
     * Sorts the string list by repeatedly calling the 'Comparator' function until
     * the list is in order.
     *
     * param StringList
     *     String list to sort.
     * param Comparator
     *     Function called to compare two string list strings or NULL for the
     *     default sort. The default sorting handles NULL elements and uses the
     *     system's strcmp utility for comparing valid strings elements.
     * param ClientData
     *     Contextual information that is passed along to the comparator function.
     */
    void StringListSort(StringList_pa                 StringList,
                        StringListStringComparator_pf Comparator,
                        ArbParam_t                    ClientData)
    {
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        REQUIRE(VALID_REF(StringList));
        REQUIRE(VALID_FN_REF(Comparator) || Comparator == NULL);

        /* set up for comparator proxy */
        if (Comparator != NULL)
            ComparatorFunction = Comparator;
        else
            ComparatorFunction = DefaultStrcmpComparator;

        /* sort the array using the comparator proxy to forward */
        /* the comparison request to the supplied comparator    */
        ArrayListQSort((ArrayList_pa)StringList, ComparatorProxy, ClientData);

        /* cleanup */
        ComparatorFunction = NULL;
    }
