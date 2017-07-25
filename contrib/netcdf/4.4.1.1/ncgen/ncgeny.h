/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NC_UNLIMITED_K = 258,
     CHAR_K = 259,
     BYTE_K = 260,
     SHORT_K = 261,
     INT_K = 262,
     FLOAT_K = 263,
     DOUBLE_K = 264,
     UBYTE_K = 265,
     USHORT_K = 266,
     UINT_K = 267,
     INT64_K = 268,
     UINT64_K = 269,
     IDENT = 270,
     TERMSTRING = 271,
     CHAR_CONST = 272,
     BYTE_CONST = 273,
     SHORT_CONST = 274,
     INT_CONST = 275,
     INT64_CONST = 276,
     UBYTE_CONST = 277,
     USHORT_CONST = 278,
     UINT_CONST = 279,
     UINT64_CONST = 280,
     FLOAT_CONST = 281,
     DOUBLE_CONST = 282,
     DIMENSIONS = 283,
     VARIABLES = 284,
     NETCDF = 285,
     DATA = 286,
     TYPES = 287,
     COMPOUND = 288,
     ENUM = 289,
     OPAQUE = 290,
     OPAQUESTRING = 291,
     GROUP = 292,
     PATH = 293,
     FILLMARKER = 294,
     NIL = 295,
     _FILLVALUE = 296,
     _FORMAT = 297,
     _STORAGE = 298,
     _CHUNKSIZES = 299,
     _DEFLATELEVEL = 300,
     _SHUFFLE = 301,
     _ENDIANNESS = 302,
     _NOFILL = 303,
     _FLETCHER32 = 304,
     _NCPROPS = 305,
     _ISNETCDF4 = 306,
     _SUPERBLOCK = 307,
     DATASETID = 308
   };
#endif
/* Tokens.  */
#define NC_UNLIMITED_K 258
#define CHAR_K 259
#define BYTE_K 260
#define SHORT_K 261
#define INT_K 262
#define FLOAT_K 263
#define DOUBLE_K 264
#define UBYTE_K 265
#define USHORT_K 266
#define UINT_K 267
#define INT64_K 268
#define UINT64_K 269
#define IDENT 270
#define TERMSTRING 271
#define CHAR_CONST 272
#define BYTE_CONST 273
#define SHORT_CONST 274
#define INT_CONST 275
#define INT64_CONST 276
#define UBYTE_CONST 277
#define USHORT_CONST 278
#define UINT_CONST 279
#define UINT64_CONST 280
#define FLOAT_CONST 281
#define DOUBLE_CONST 282
#define DIMENSIONS 283
#define VARIABLES 284
#define NETCDF 285
#define DATA 286
#define TYPES 287
#define COMPOUND 288
#define ENUM 289
#define OPAQUE 290
#define OPAQUESTRING 291
#define GROUP 292
#define PATH 293
#define FILLMARKER 294
#define NIL 295
#define _FILLVALUE 296
#define _FORMAT 297
#define _STORAGE 298
#define _CHUNKSIZES 299
#define _DEFLATELEVEL 300
#define _SHUFFLE 301
#define _ENDIANNESS 302
#define _NOFILL 303
#define _FLETCHER32 304
#define _NCPROPS 305
#define _ISNETCDF4 306
#define _SUPERBLOCK 307
#define DATASETID 308




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 136 "ncgen.y"
{
Symbol* sym;
unsigned long  size; /* allow for zero size to indicate e.g. UNLIMITED*/
long           mark; /* track indices into the sequence*/
int            nctype; /* for tracking attribute list type*/
Datalist*      datalist;
NCConstant       constant;
}
/* Line 1529 of yacc.c.  */
#line 164 "ncgeny.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE ncglval;

