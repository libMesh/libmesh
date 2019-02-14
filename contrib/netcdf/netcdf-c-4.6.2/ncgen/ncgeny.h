/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

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

#ifndef YY_NCG_NCGEN_TAB_H_INCLUDED
# define YY_NCG_NCGEN_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int ncgdebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
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
    STRING_K = 270,
    IDENT = 271,
    TERMSTRING = 272,
    CHAR_CONST = 273,
    BYTE_CONST = 274,
    SHORT_CONST = 275,
    INT_CONST = 276,
    INT64_CONST = 277,
    UBYTE_CONST = 278,
    USHORT_CONST = 279,
    UINT_CONST = 280,
    UINT64_CONST = 281,
    FLOAT_CONST = 282,
    DOUBLE_CONST = 283,
    DIMENSIONS = 284,
    VARIABLES = 285,
    NETCDF = 286,
    DATA = 287,
    TYPES = 288,
    COMPOUND = 289,
    ENUM = 290,
    OPAQUE_ = 291,
    OPAQUESTRING = 292,
    GROUP = 293,
    PATH = 294,
    FILLMARKER = 295,
    NIL = 296,
    _FILLVALUE = 297,
    _FORMAT = 298,
    _STORAGE = 299,
    _CHUNKSIZES = 300,
    _DEFLATELEVEL = 301,
    _SHUFFLE = 302,
    _ENDIANNESS = 303,
    _NOFILL = 304,
    _FLETCHER32 = 305,
    _NCPROPS = 306,
    _ISNETCDF4 = 307,
    _SUPERBLOCK = 308,
    _FILTER = 309,
    DATASETID = 310
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 149 "ncgen.y" /* yacc.c:1909  */

Symbol* sym;
unsigned long  size; /* allow for zero size to indicate e.g. UNLIMITED*/
long           mark; /* track indices into the sequence*/
int            nctype; /* for tracking attribute list type*/
Datalist*      datalist;
NCConstant*    constant;

#line 119 "ncgeny.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE ncglval;

int ncgparse (void);

#endif /* !YY_NCG_NCGEN_TAB_H_INCLUDED  */
