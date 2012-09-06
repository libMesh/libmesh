#ifndef IOUTIL_h__
#define IOUTIL_h__
/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** (C) Copyright 1989-1998  by AMTEC ENGINEERING INC. *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

/* TECADDON.h *must* be included first */
#ifndef ASSERT
# ifndef assert
#   include <assert.h>
# endif
# define ASSERT(exp) assert(exp)
#endif

/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** (C) Copyright 1989-1998  by AMTEC ENGINEERING INC.********
*******       All Rights Reserved.                       ********
*******                                                  ********
*****************************************************************
*****************************************************************
*/



/********************************************************
*
* IOUTIL.h
*
* Addon utility functions
*
*********************************************************
*/
int IOUtilGetVersion(void);

void IOUtilDoSwapBytes(Boolean_t swap);

/*
* --------------------------------------------------------
* Input:  
*   
*   swap: TRUE if the IOUtil routines should swap bytes
*               when reading numbers from files,
*         FALSE otherwise
*
* If this function is not called, IOUtil routines will
* not swap bytes
*
* Error Conditions: None
*
*
*=========================================================*/


Boolean_t IOUtilIsAsciiFile(const char *fname) ;
/*
*
* This functions attempts to determine if
* a given file is an ASCII file by counting
* the number of printable characters in the first 200
* bytes of the file.
*
* NOTE: NOT guarenteed to be correct 100% of the time
* If the ascii file should have only numbers, more reliable results
* will be obtained by using IOUtilIsAsciiNumFile()
* 
*
* The sample size may be adjusted by #defining IOUTIL_ASCII_SAMPLE_SIZE
* to the sample size (i.e., default = 200) and recompiling IOUTIL.c
*
* Input:
*
*   fname:  filename
*
* Error Conditions: returns FALSE if fname is NULL or invalid
*
* Return: TRUE if there is a high (but not 100%)
*   probabilty that 'fname' is an ASCII file
*
*
*
*=================================================================*/

Boolean_t IOUtilIsAsciiNumFile(const char *fname) ;
/*
*
* --> Same as above but more reliable if the file has only numbers
*
*==============================================================*/


Boolean_t IOUtilIsFortranFile(const char *fname);
/*
*
* This functions determines if a file is a FORTRAN
* unformatted binary file.
*
* Input:
*
*   fname:  filename
*
* Error Conditions: returns FALSE if fname is NULL or invalid
*
* Returns:  TRUE if this file is a FORTRAN unformatted binary file
*

*============================================================*/


Boolean_t IOUtilReadBinaryInt4(FILE *f, 
                               LgInteger_t *val);
/*
* ------------------------------------------------------------
* Reads a 4 byte integer from a file, swapping bytes if necessary
*
* Input:
*
*   f:    file to read from
*   val:  value read
*
* Returns:  TRUE if no error, FALSE otherwise
*
*=============================================================*/
Boolean_t IOUtilReadBinaryInt2(FILE *f, 
                               SmInteger_t *val);

Boolean_t IOUtilReadBinaryFloat(FILE *f, 
                                float *val);

Boolean_t IOUtilReadBinaryDouble(FILE *f, 
                                 double *val);


int IOUtilReadAsciiToken(FILE* f,
                         char buffer[], 
                         int len);
/*
* 
* Reads the next token from the file.
* 'tokens' are strings of printable characters separated by one of
* the following:
*
* ' '   (space)
* '\t'
* '\n'
* '\r'
*
* Input:
*
*   f:      file to read from
*   buffer  where to put the token string
*   len     length of the buffer
*
* Returns:  length of the string, or -1 if error
*           buffer[] will have the token
*
*
*=================================================================*/



LgInteger_t IOUtilReadAsciiDoubleArray(FILE *f, 
                                       LgInteger_t count,
                                       double data[]);
/*
*
* Reads an array of size 'count' numbers from an ascii file
*
*
* Input:
*
*   f:      file to read from (must be an ascii file)
*   count:  number of values to read
*   data[]: array of doubles to put the values in
*           (must have at least 'count' elements
*
*
* Returns:
*   number of values read, or -1 if an error
*
*====================================================================*/




LgInteger_t IOUtilReadFortranFloatArray(FILE *f, 
                                        float data[]);
/*
*
* Reads an array of size 'count' floats from a
* FORTRAN unformatted binary file
*
* The array should be a single line of the binary file.
*
* That is, this function expects to
* read the following:
*
* 1) a starting record length of nn
* 2) nn/sizeof(float) floats
* 3) an ending record length of nn
*
* Input:
*
*   f:    file to read from
*   data[]: array of floats to read values from
*
*   if data[] is NULL, then this function
*   returns the size of the start_record and
*   does not change the file position
*
*
* Returns:
*   number of values read, or -1L if an error
*
*====================================================================*/
LgInteger_t IOUtilReadFortranDoubleArray(FILE  *f, 
                                         double data[]);

LgInteger_t IOUtilReadFortranInt4Array(FILE       *f, 
                                       LgInteger_t data[]);


/*====================================================================*/
LgInteger_t IOUtilReadFortranRealArray(FILE  *f, 
                                       int    size, 
                                       double data[]);
/*
* Reads an array of floats (size == 4) or 
* doubles (size == 8) from a fortran 
* unformatted file and puts them
* into the (double) array data[]
*
*/

LgInteger_t IOUtilReadRealArray(FILE     *f, 
                                int       size, 
                                LgIndex_t count, 
                                double    data[]);
/*
* Same as above except the file to be read from
* is pure binary,
* so the caller must supply a count of how
* many values to read
*/

Boolean_t IOUtilIsBlankLine(char *line);
/*
* returns TRUE if line has no printable characters
*
*/

Boolean_t IOUtilFileExists(char *FName);

/*
* reads a single line of text,
* replacing the '\n' with a NULL.
*/
Boolean_t IOUtilReadLine(char *s, 
                         int   n, 
                         FILE *f);

/* 
* --------------------------------------------------------
*
* Some non-io (but useful) functions
*
* ---------------------------------------------------------
*/
char *IOUtilTrimLeft(const char *S); /* removes whitespace from the left of S */

#endif  /* IOUTIL_h__ */


