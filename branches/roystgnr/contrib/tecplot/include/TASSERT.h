/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** (C) Copyright 1988-2000  by AMTEC ENGINEERING INC. *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/
/*
 * Provide four levels of assertion control. Assertions provide a mechanism
 * to enforce a contract between a client and service provider. The assertions
 * are listed in order of highest to lowest priority. Assertions can be turned
 * off individually by defining the appropriate name (see preprossessor
 * definitions below), however, lower priority assertions should be turned
 * off prior to higher ones. As confidence in the code increases all assertions
 * can be turned off by defining NO_ASSERTS.
 *
 * The assertions defined below have the following meanings:
 *
 *     INVARIANT - Asserts that a property's state is invariant throughout the
 *                 life of the property's scope. For instance in Tecplot
 *                 the 'CurFrame' variable has global scope and must always
 *                 pass the VALID_REF test both before and after any public
 *                 method call. Stating invariant properties of an application
 *                 provides a deeper understanding of the application's state.
 *                 These statements are usually positioned just ahead of the
 *                 preconditions and just after the postconditions.
 *
 *     REQUIRE   - Asserts that a method's preconditions are within their
 *                 valid domains. Preconditions are conditions placed upon
 *                 any state information relied upon for the call. These
 *                 statements should be as close to the top of the method
 *                 as possible (except for assertions on invariant properties).
 *
 *     ENSURE    - Asserts that a method's postconditions are within their
 *                 valid ranges. Postconditions are conditions placed upon
 *                 any state information modified by the call. These
 *                 statements should be as close to the bottom of the method
 *                 (presumably there is only one exit point) as possible
 *                 (except for assertions on invariant properties).
 *
 *     CHECK     - Any other assertion not covered by the above assertions.
 *                 These are often added within a method body to specify
 *                 something that may not be immediately obvious to the reader
 *                 or to validate your assumptions about a call to a 3rd party
 *                 method that does not use runtime assertions for its
 *                 preconditions or postconditions. Obviously if the 3rd party
 *                 method uses assertions then there is no need for the CHECK.
 *
 * Additionally a convenience macro is available to place in code that is
 * pending implementation.
 *
 *     NOT_IMPLEMENTED - Assertion that always fails during runtime for debug
 *                       builds and always fails at compile time for release
 *                       builds.
 */
#if !defined TASSERT_H
#define TASSERT_H

#if defined (MSWIN)
# include <assert.h>
#endif /* MSWIN */

#if !defined TECPLOTKERNEL && !defined STD_ASSERTS
#define STD_ASSERTS
#endif


#if !defined (MSWIN)
#  if defined LINUX
#    define ASSERT(x)
#  else
#    include <assert.h>
#    define ASSERT assert
#  endif
#endif

/* CORE SOURCE CODE REMOVED */

#ifdef MOTIF
/* CORE SOURCE CODE REMOVED */

/* CORE SOURCE CODE REMOVED */
#  define VALID_REF(p)      ( (p)  != NULL )
#  define VALID_FN_REF(fp)  ( (fp) != NULL )
/* CORE SOURCE CODE REMOVED */

/* CORE SOURCE CODE REMOVED */
#endif /* MOTIF */
#ifdef MSWIN
   /* Under Windows, we can use AfxIsValidAddress to check for */
   /* correct addresses with in the program.  In the future,   */
   /* could use AfxIsValidMemBlock to check for valid heap     */
   /* addresses when we know we are to be looking at the heap. */
/* CORE SOURCE CODE REMOVED */
# if defined(_AFX)
/* VC++ 6.0 puts literal strings (stuff in "") in read-only memory */
/* #define VALID_REF(pointer)         AfxIsValidAddress((void*)(pointer),1,TRUE) */
#define VALID_REF(pointer)         AfxIsValidAddress((void*)(pointer),1,FALSE)
#   define VALID_FN_REF(fn_pointer)   AfxIsValidAddress((fn_pointer),1,FALSE)
# else /* !_AFX */
/* VC++ 6.0 puts literal strings (stuff in "") in read-only memory */
/* #   define VALID_REF(p)     ( p != NULL && !IsBadReadPtr((CONST VOID *)p,1) && !IsBadWritePtr((LPVOID)p,1) ) */
#   define VALID_REF(p)     ( p != NULL && !IsBadReadPtr((CONST VOID *)p,1) )
#   define VALID_FN_REF(pf) ( pf != NULL && !IsBadReadPtr((CONST VOID *)pf,1) )
# endif /* _AFX */
/* CORE SOURCE CODE REMOVED */
#endif /* MSWIN */
/* CORE SOURCE CODE REMOVED */
/* other useful validity checks */
#define VALID_BOOLEAN(b)           ((b) == TRUE || (b) == FALSE)
#define VALID_ENUM(value, type)    (0 <= (value) && \
                                         (value) < END_##type)
/* CORE SOURCE CODE REMOVED */

/* Check for a non-zero length string */
#define VALID_STR(str) (VALID_REF(str) && strlen(str) > 0)

/* Check for valid stdio file handle */
#define VALID_FILE_HANDLE(stream) ((stream) != NULL)

/* To check colors and pen numbers */
/* CORE SOURCE CODE REMOVED */

#if defined NO_ASSERTS
/* CORE SOURCE CODE REMOVED */
#   define INVARIANT(EXPR)
#   define REQUIRE(EXPR)
#   define ENSURE(EXPR)
#   define CHECK(EXPR)
#   define VERIFY(EXPR)  ((void)(EXPR))
    /* 
     * Only define IGNORENOTIMPLEMENTED if building a "test" release build
     * that you are fully aware may contain unimplemented features.
     */
#   if defined IGNORENOTIMPLEMENTED
#     define NOT_IMPLEMENTED()
#   else
#     define NOT_IMPLEMENTED() not implemented /* intentionally break the compile */
#   endif
#elif defined STD_ASSERTS
/* CORE SOURCE CODE REMOVED */
#   define INVARIANT(EXPR)       assert(EXPR)
#   define REQUIRE(EXPR)         assert(EXPR)
#   define ENSURE(EXPR)          assert(EXPR)
#   define CHECK(EXPR)           assert(EXPR) 
#   ifndef VERIFY
#     define VERIFY(EXPR)        assert(EXPR)
#   endif /* VERIFY */
#   define NOT_IMPLEMENTED()     assert(!("Not Implemented"))
#else
/* CORE SOURCE CODE REMOVED */
#endif
/* CORE SOURCE CODE REMOVED */

#endif /* TASSERT_H */
