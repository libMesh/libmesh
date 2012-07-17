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
#if defined EXTERN
#undef EXTERN
#endif
#if defined SETMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

#ifndef _SET_H_INCLUDED
#define _SET_H_INCLUDED

/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** Copyright (C) 1988-2008 Tecplot, Inc.              *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#define PadOut(X,Y)           ((int)(((X)-1)/(Y)+1)*(Y))
#define SetBitSize            (8*sizeof(SetData_t))
#define SetLastBit            (((unsigned long)1)<<(SetBitSize-1))

#if defined _DEBUG
#  define USE_FUNCTIONS_FOR_SETS
#endif

/* *
 * * NOTE: "Set_pa" is a pointer to an "abstract type",
 * * hence the "_pa". Pointer here is akin to "handle".
 * * Any routines dealing with the internals of Set_pa
 * * or Set_a must be in the same file as these routines
 * */

/* Set_a is intentionally not defined to further
 * deter usage of this private structure */
struct _Set_a
{
    /* * PRIVATE * */
    SetIndex_t size;
    SetData_pt data;
};

/*
 * Checks set for NULL.
 */
#define IsSetNull(Set) ((Set)==NULL)

/**
 * Indicates how many bytes are required to store the set data.
 */
inline size_t SetDataSizeInBytes(Set_pa Set)
{
    REQUIRE(VALID_REF(Set));
    return Set->size / SetBitSize * sizeof(SetData_t);
}

/*
 * Allocates a new empty set.  Returns NULL if not enough memory.
 */
EXTERN Set_pa AllocSet(Boolean_t show_error_msg);

/*
 * Frees all memory associated with set "*set", and
 * sets "*set" to NULL.
 */
EXTERN void DeallocSet(Set_pa *Set);

/**
 * This function adapts the DeallocSet function to work with the
 * ArrayList's deallocation callback.
 */
EXTERN Boolean_t SetItemDestructor(void       *ItemRef,
                                   ArbParam_t  ClientData);
/*
 * Makes sure set "set" can hold at least "max_val" elements.
 * Returns TRUE if successful, FALSE otherwise.  A successful
 * call to ExpandSet() guarentees that any calls to AddToSet()
 * will be successful as long as the elements added are less
 * than "max_val".
 */
EXTERN Boolean_t ExpandSet(Set_pa     Set,
                           SetIndex_t max_val,
                           Boolean_t  show_error_msg);

/*
 * Copies set "src" to set "dst".  Returns TRUE if successful,
 * FALSE if "src" contains elements it is unable to add to "dst".
 */
EXTERN Boolean_t CopySet(Set_pa    dst,
                         Set_pa    src,
                         Boolean_t show_error_msg);

/*
 * Appends set "src" to set "dst".  Returns TRUE if successful,
 * FALSE if "src" contains elements it is unable to add to "dst".
 */
EXTERN Boolean_t AppendSet(Set_pa dst,
                           Set_pa src,
                           Boolean_t show_error_msg);
/*
 * Empties the set "set".
 */
EXTERN void ClearSet(Set_pa Set);

/*
 * Adds "member" to set "set".  Returns TRUE if successful,
 * FALSE otherwise.  AddToSet() can only return FALSE if
 * "member" is greater than any previous member of "set" and
 * also greater that any "max_val" set with ExpandSet().
 */
#if defined USE_FUNCTIONS_FOR_SETS
EXTERN Boolean_t AddToSet(Set_pa     Set,
                          SetIndex_t member,
                          Boolean_t  show_error_msg);
#else
#  if defined __cplusplus
inline Boolean_t AddToSet(Set_pa     Set,
                          SetIndex_t member,
                          Boolean_t  show_error_msg)
{
    if (Set &&
        (member + 1 <= Set->size ||
         ExpandSet(Set, member + 1, show_error_msg)))
    {
        SetIndex_t word = member / SetBitSize;
        SetData_t  bit = (SetData_t)1 << (member % SetBitSize);
        Set->data[word] |= bit;
        return TRUE;
    }
    else
        return FALSE;
} /* AddToSet() */
#  else
#    define AddToSet(Set,member,show_error_msg) \
       (((Set) && \
         ((member)+1 <= (Set)->size || \
         ExpandSet((Set), (member)+1, (show_error_msg)))) \
           ? (((Set)->data[(member) / SetBitSize] |= (SetData_t)1 << ((member) % SetBitSize)), TRUE) \
           : FALSE)
#  endif
#endif

/*
 * Removes "member" from set "set".
 */
EXTERN void RemoveFromSet(Set_pa     Set,
                          SetIndex_t member);

EXTERN void DeleteSetMember(Set_pa     Set,
                            SetIndex_t Member);
EXTERN Boolean_t InsertSetMember(Set_pa     Set,
                                 SetIndex_t Member,
                                 Boolean_t  ShowErrMsg);
/*
 * Test for membership of "member" in set "set".  This is the only
 * function worth making into a macro or inline function.
 */
#if defined USE_FUNCTIONS_FOR_SETS
EXTERN Boolean_t InSet(Set_pa     Set,
                       SetIndex_t member);
#else
#  if defined __cplusplus
inline Boolean_t InSet(Set_pa     Set,
                       SetIndex_t member)
{
    if (Set && (0 <= member && member < Set->size))
    {
        SetIndex_t word = member / SetBitSize;
        SetData_t  bit = (SetData_t)1 << (member % SetBitSize);
        return (Set->data[word]&bit) != 0;
    }
    else
        return FALSE;
} /* InSet() */
#  else
#    define InSet(Set,member)  ((Set && (0<=(member) && (member)<(Set)->size)) \
                                ? ((Set)->data[(member)/SetBitSize]&((SetData_t)1<<((member)%SetBitSize)))!=0 \
                                : FALSE)
#  endif
#endif

/*
 * Returns TRUE if set "set" is empty.
 */
EXTERN Boolean_t IsEmpty(Set_pa Set);

/*
 * Returns TRUE if Set has voids.
 */
EXTERN Boolean_t HasVoids(Set_pa Set);

/*
 * Returns number of members in Set "Set".
 */
EXTERN SetIndex_t MemberCount(Set_pa Set);

/*
 * Returns the next member in set "set" after member "start_at".
 * Use "start_at" of BAD_ZV_VALUE to find first member.
 */
EXTERN SetIndex_t GetNextMember(Set_pa     Set,
                                SetIndex_t start_at);

/*
 * Returns the previous member in set "set" before member
 * "start_at".  Use "start_at" of BAD_ZV_VALUE to find last member.
 */
EXTERN SetIndex_t GetPrevMember(Set_pa     Set,
                                SetIndex_t start_at);

/*
 * Returns TRUE if sets are equal (have same members).  FALSE otherwise.
 */
EXTERN Boolean_t EqualSets(Set_pa  set1,
                           Set_pa  set2);

/*
 * Returns TRUE if all members of childset are contained in parentset.
 */
EXTERN Boolean_t IsSubSet(Set_pa childset,
                          Set_pa parentset);

EXTERN SetIndex_t MemberOffset(Set_pa Set,
                               SetIndex_t    Member);

EXTERN SetIndex_t OffsetMember(Set_pa Set,
                               SetIndex_t    Offset);


EXTERN Boolean_t CopySetMember(Set_pa     DstSet,
                               SetIndex_t DstOffset,
                               Set_pa     SrcSet,
                               SetIndex_t SrcOffset);

EXTERN void ShiftSet(Set_pa     Set,
                     SetIndex_t ShiftPos1,
                     SetIndex_t ShiftPos2,
                     SetIndex_t ShiftAmount);


/*
 *  Handy macros
 */
#define GetFirstSetMember(Set) (GetNextMember((Set), BAD_SET_VALUE))
#define GetLastSetMember(Set)  (GetPrevMember((Set), BAD_SET_VALUE))

#define ForAllMembersInSet(Member, Set) \
            for (Member = GetFirstSetMember((Set)); \
                 Member != BAD_SET_VALUE; \
                 Member = GetNextMember((Set), (Member)))
#define ForAllMembersInReversedSet(Member, Set) \
            for (Member = GetLastSetMember((Set)); \
                 Member != BAD_SET_VALUE; \
                 Member = GetPrevMember((Set), (Member)))

#endif // _SET_H_INCLUDED
