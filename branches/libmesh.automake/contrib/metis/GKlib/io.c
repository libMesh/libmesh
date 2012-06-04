/*!
\file  io.c
\brief Various file I/O functions.

This file contains various functions that perform I/O.

\date Started 4/10/95
\author George
\version\verbatim $Id$ \endverbatim
*/

#ifdef HAVE_GETLINE
/* Get getline to be defined. */
#define _GNU_SOURCE
#include <stdio.h>
#undef _GNU_SOURCE
#endif

#include <GKlib.h>

/*************************************************************************
* This function opens a file
**************************************************************************/
FILE *gk_fopen(char *fname, char *mode, const char *msg)
{
  FILE *fp;
  char errmsg[8192];

  fp = fopen(fname, mode);
  if (fp != NULL)
    return fp;

  sprintf(errmsg,"file: %s, mode: %s, [%s]", fname, mode, msg);
  perror(errmsg);
  errexit("Failed on gk_fopen()\n");

  return NULL;
}


/*************************************************************************
* This function closes a file
**************************************************************************/
void gk_fclose(FILE *fp)
{
  fclose(fp);
}


/*************************************************************************/
/*! This function is the GKlib implementation of glibc's getline()
    function.
    \returns -1 if the EOF has been reached, otherwise it returns the 
             number of bytes read.
*/
/*************************************************************************/
gk_idx_t gk_getline(char **lineptr, size_t *n, FILE *stream)
{
#ifdef HAVE_GETLINE
  return getline(lineptr, n, stream);
#else
  size_t i;
  int ch;

  if (feof(stream))
    return -1;  

  /* Initial memory allocation if *lineptr is NULL */
  if (*lineptr == NULL || *n == 0) {
    *n = 1024;
    *lineptr = gk_malloc((*n)*sizeof(char), "gk_getline: lineptr");
  }

  /* get into the main loop */
  i = 0;
  while ((ch = getc(stream)) != EOF) {
    (*lineptr)[i++] = (char)ch;

    /* reallocate memory if reached at the end of the buffer. The +1 is for '\0' */
    if (i+1 == *n) { 
      *n = 2*(*n);
      *lineptr = gk_realloc(*lineptr, (*n)*sizeof(char), "gk_getline: lineptr");
    }
      
    if (ch == '\n')
      break;
  }
  (*lineptr)[i] = '\0';

  return (i == 0 ? -1 : i);
#endif
}


/*************************************************************************/
/*! This function reads the contents of a text file and returns it in the
    form of an array of strings.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
*/
/*************************************************************************/
char **gk_readfile(char *fname, gk_idx_t *r_nlines)
{
  size_t lnlen, nlines;
  char *line=NULL, **lines=NULL;
  FILE *fpin;

  gk_getfilestats(fname, &nlines, NULL, NULL, NULL);
  if (nlines > 0) {
    lines = (char **)gk_malloc(nlines*sizeof(char *), "gk_readfile: lines");

    fpin = gk_fopen(fname, "r", "gk_readfile");
    nlines = 0;
    while (gk_getline(&line, &lnlen, fpin) != -1) {
      gk_strtprune(line, "\n\r");
      lines[nlines++] = gk_strdup(line);
    }
    gk_fclose(fpin);
  }

  gk_free((void **)&line, LTERM);

  if (r_nlines != NULL)
    *r_nlines  = nlines;

  return lines;
}


/*************************************************************************/
/*! This function reads the contents of a file and returns it in the
    form of an array of int32_t.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
*/
/*************************************************************************/
int32_t *gk_i32readfile(char *fname, gk_idx_t *r_nlines)
{
  size_t lnlen, nlines;
  char *line=NULL;
  int32_t *array=NULL;
  FILE *fpin;

  gk_getfilestats(fname, &nlines, NULL, NULL, NULL);
  if (nlines > 0) {
    array = gk_i32malloc(nlines, "gk_i32readfile: array");

    fpin = gk_fopen(fname, "r", "gk_readfile");
    nlines = 0;

    while (gk_getline(&line, &lnlen, fpin) != -1) {
      sscanf(line, "%"SCNd32, &array[nlines++]);
    }

    gk_fclose(fpin);
  }

  gk_free((void **)&line, LTERM);

  if (r_nlines != NULL)
    *r_nlines  = nlines;

  return array;
}


/*************************************************************************/
/*! This function reads the contents of a file and returns it in the
    form of an array of int64_t.
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
*/
/*************************************************************************/
int64_t *gk_i64readfile(char *fname, gk_idx_t *r_nlines)
{
  size_t lnlen, nlines;
  char *line=NULL;
  int64_t *array=NULL;
  FILE *fpin;

  gk_getfilestats(fname, &nlines, NULL, NULL, NULL);
  if (nlines > 0) {
    array = gk_i64malloc(nlines, "gk_i64readfile: array");

    fpin = gk_fopen(fname, "r", "gk_readfile");
    nlines = 0;

    while (gk_getline(&line, &lnlen, fpin) != -1) {
      sscanf(line, "%"SCNd64, &array[nlines++]);
    }

    gk_fclose(fpin);
  }

  gk_free((void **)&line, LTERM);

  if (r_nlines != NULL)
    *r_nlines  = nlines;

  return array;
}

