/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * stdheaders.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: stdheaders.h,v 1.3 2005-05-06 17:43:43 roystgnr Exp $
 */


#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#elif defined(__APPLE__)
#include <stdlib.h>
#else
#include <malloc.h>
#endif
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

