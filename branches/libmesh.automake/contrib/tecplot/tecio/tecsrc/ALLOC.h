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
/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2009 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
#ifndef ALLOC_H
#define ALLOC_H

#include "TASSERT.h"
#if defined __cplusplus
#include <new>
#endif

#if !defined __cplusplus
#define ALLOC_ARRAY(N,Type,str) (Type *)malloc((N)*sizeof(Type))
#define ALLOC_ITEM(Type,str)    (Type *)malloc(sizeof(Type))
#ifdef _DEBUG
/* NOTE: the pointer is set to 0xFFFF after the free for debug   */
/*       versions in the hopes of catching invalid pointer usage */
#define FREE_ARRAY(X,str)  do { free((void *)(X)); *((void **)&(X)) = (void *)0xFFFF; } while (0)
#define FREE_ITEM(X,str)   do { free((void *)(X)); *((void **)&(X)) = (void *)0xFFFF; } while (0)
#else
#define FREE_ARRAY(X,str)  free((void *)(X))
#define FREE_ITEM(X,str)   free((void *)(X))
#endif
#else
#ifdef TRACK_MEMORY_USAGE
extern void initMemoryUsageTracking(void);
extern void cleanUpMemoryUsageTracking(void);
extern void trackMemoryAlloc(size_t size);
extern void trackMemoryFree(size_t size);
extern void trackMemoryClearHighMark(void);
extern void trackMemorySaveHighMark(void);
extern void getMemoryUsage(size_t* memoryInUse,
                           size_t* memoryCurrentHighMark,
                           size_t* memorySavedHighMark,
                           size_t* memoryTotalHighMark);
#endif
/*
* Create a version of new that returns NULL instead
* of throwing std::bad_alloc.  A lot of code is written using
* ALLOC_ITEM and ALLOC_ARRAY that expect a return value of
* NULL on failure instead of the exception. 2008-05-08 CAM
*/
#if defined MSWIN && defined _DEBUG
template <typename T>
inline T *nonExceptionNew(size_t      numItems,
                          const char* fileName,
                          int         lineNumber)
{
    REQUIRE(numItems > 0);
    REQUIRE(VALID_REF(fileName));
    REQUIRE(lineNumber > 0);
    T* result;
    try
    {
#ifdef DEBUG_NEW
#ifdef new
#undef new
#define USING_DEBUG_NEW
#endif
        result = new(fileName, lineNumber) T[numItems];
#ifdef USING_DEBUG_NEW
#define new DEBUG_NEW
#undef USING_DEBUG_NEW
#endif
#else
        result = new T[numItems];
#endif
    }
    catch (std::bad_alloc&)
    {
        result = NULL;
    }
#ifdef TRACK_MEMORY_USAGE
    if (result != NULL)
    {
#ifdef MSWIN
        trackMemoryAlloc(_msize(result));
#else
        trackMemoryAlloc(malloc_usable_size(result));
#endif
    }
#endif
    ENSURE(VALID_REF_OR_NULL(result));
    return result;
}
#define ALLOC_ARRAY(N,Type,str) nonExceptionNew<Type>((N),__FILE__,__LINE__)
#else
template <typename T>
inline T *nonExceptionNew(size_t numItems)
{
    REQUIRE(numItems > 0);
    T *result;
    try
    {
        result = new T[numItems];
    }
    catch (std::bad_alloc&)
    {
        result = NULL;
    }
#ifdef TRACK_MEMORY_USAGE
    if (result != NULL)
    {
#ifdef MSWIN
        trackMemoryAlloc(_msize(result));
#else
        trackMemoryAlloc(malloc_usable_size(result));
#endif
    }
#endif
    ENSURE(VALID_REF_OR_NULL(result));
    return result;
}
#define ALLOC_ARRAY(N,Type,str) nonExceptionNew<Type>((N))
#endif
#define ALLOC_ITEM(Type,str)    ALLOC_ARRAY(1,Type,str)

/*
* Although delete doesn't throw exceptions, this function matches
* nonExceptionNew, and also reports the size of the block if we
* are tracking memory.
*/
template <typename T>
inline void nonExceptionDelete(T* &ptr)
{
#if defined MSWIN && !defined NO_ASSERTS
    CHECK(!IsBadReadPtr((void*)ptr, 1));
#endif
#if defined TRACK_MEMORY_USAGE
    if (ptr != NULL)
    {
#ifdef MSWIN
        trackMemoryFree(_msize(ptr));
#else
        trackMemoryFree(malloc_usable_size(ptr));
#endif
    }
#endif

    delete [] ptr;
#if !defined NO_ASSERTS
    /*
    * NOTE: the pointer is set to 0xFFFF after the free for asserted
    * builds in the hopes of catching invalid pointer usage
    */
    ptr = (T*)(void*)0xFFFF;
#endif
}
#define FREE_ARRAY(ptr,str)  nonExceptionDelete((ptr))
#define FREE_ITEM(ptr,str)   FREE_ARRAY(ptr,str)
#endif

/**
* The following functor can be used to easily deallocate memory from containers
* that hold pointers to allocated objects. For example:
*
* vector<MyObject*> container;
* for (int ii = 0; ii < 10; ii++
*  container.push_back(new MyObject);
* ... do something with the objects ...
* ... now we need to clean up ...
* for_each(container.begin(),
*          container.end(),
*          DeleteItem());
*/
struct DeleteItem
{
    template<typename T>
    void operator()(T*& object)
    {
        delete object;
        object = NULL;
    }
};

#endif /* ALLOC_H */
